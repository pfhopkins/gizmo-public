/* gpu_gravtree.cc
 *
 * GPU gravity walk (mode=0 primary-tree path).  Core walk with
 * PMGRID, ADAPTIVE_GRAVSOFT_FORALL, SYMMETRIZE, EVALPOTENTIAL, plus
 * RT cluster payloads (RT_USE_GRAVTREE, GALSF_FB_FIRE_RT_LONGRANGE,
 * CHIMES_STELLAR_FLUXES, RT_USE_TREECOL_FOR_NH).
 *
 * rt_get_source_luminosity() is not GPU-callable; it is pre-computed on CPU
 * for all local particles into a SharedSpace array (d_src_lum) before kernel
 * launch.  Node stellar luminosities come from the SoA extension.
 * rt_kappa() is KOKKOS_INLINE_FUNCTION and runs on device for RT_LEBRON
 * fac_stellum initialisation.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
/* This TU alone exposes the device-callable gravtree source-payload helpers
 * (rt_get_source_luminosity / sink_lum_bol_core / cr_get_source_injection_rate and
 * their stellar-evolution/cosmology leaves), so the source payload can be evaluated
 * on-device at each local particle-open. The compile-time capability predicate
 * GRAVTREE_SOURCE_LAZY_SUPPORTED gates body visibility; the runtime eager/lazy choice
 * is a separate active-count threshold. Must precede the source-helper header includes. */
#ifdef GRAVTREE_SOURCE_LAZY_SUPPORTED
#define GRAVTREE_SOURCE_DEVICE_TU
#endif
#include "../system/gpu_particles_arena.h"
#include "../core/timestep_functions.h"   /* Hermite source eligibility + prediction, shared verbatim with the host walk */
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"   /* the backend's lane count, which decides how lanes spread over targets */
#include "gpu_gravity_tree.h"
#include "../mesh/gpu_neighbor_list.h"            /* the particle touched set's epoch lifecycle (host side) */
#include "gpu_gravtree.h"
#include "forcetree.h"
#include "gravity_box_distance.h"   /* shared CPU/GPU gravity box-distance SSOT */
#include "gravtree_opening.h"       /* shared CPU/GPU primary-walk acceptance-geometry predicate (SSOT) */
#include "let_data.h"             /* LET_LEAF_TAG_* + grav_classify_node (import topology vocabulary) */

#include "../mesh/kernel.h"
#include "gravtree_force_kernel.h"  /* shared CPU/GPU accepted-source contribution physics (SSOT) */
#include "gravtree_ewald.h"         /* shared CPU/GPU Ewald image-correction trilinear interp (SSOT) */
#include "gravtree_moment_kernel.h" /* the node-motion arithmetic the drift runs, for predicting a node on read */
#include "pm_highres_region.h"      /* pmforce_is_particle_high_res SSOT (device-callable) */
/* gravtree_moment_sources.h (the SSOT source-input fill helper) is included further
 * below, AFTER the device-callable source cores, so that in this DEVICE_TU the helper
 * binds its RT/sink/CR calls to the inline device bodies rather than the proto.h host
 * decls. See the source-core include block after the walk-data struct definitions. */


/* Single gate for the Ewald periodic-image POTENTIAL correction added in the
 * primary walk (item #11): pure-tree periodic gravity with potentials requested.
 * Mirrors the CPU gate at forcetree.cc:2299.  Defined once so the four-flag
 * condition lives in one place and is referenced (not re-spelled) at the
 * table-acquire site, the host helper, and the walk body. */
#if defined(EVALPOTENTIAL) && defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
#define GIZMO_GPU_EWALD_POT_CORRECTION
#endif

/* Globals that live at file-scope in gravtree.cc without a header declaration. */
extern int Ewald_iter;
extern double Costtotal;

#ifdef PMGRID
/* Short-range tables live as file-scope (non-static) globals in forcetree.cc.
 * Table length owned by gravtree_force_kernel.h (shared with the CPU walk). */
#define GIZMO_GPU_GRAVTREE_NTAB GRAVTREE_SHORTRANGE_NTAB
extern float shortrange_table[GIZMO_GPU_GRAVTREE_NTAB];
extern float shortrange_table_potential[GIZMO_GPU_GRAVTREE_NTAB];
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
extern float shortrange_table_tidal[GIZMO_GPU_GRAVTREE_NTAB];
#endif
#endif

/* GPU-callable accessor for the cached force-softening kernel radius. Always available
 * (Pp[p].ForceSoftening is populated for every build by compute_all_force_softening()), so the
 * walk can load a leaf's secondary softening unconditionally -- mirrors the CPU
 * ForceSoftening_KernelRadius(). The actual computation lives in
 * compute_force_softening_kernel_radius(p) in forcetree.cc; new softening physics goes there
 * and is picked up here with no GPU-side change. Pp[p].ForceSoftening is the single source of truth. */
static KOKKOS_INLINE_FUNCTION
double gpu_force_softening_kernel_radius(const struct particle_data *Pp, int p)
{
    return Pp[p].ForceSoftening;
}

/* AGS_zeta field is gated (particle_data.h:331) on
 * ADAPTIVE_GRAVSOFT_FORGAS || AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
 * (the latter auto-defines under FORALL/CBE_INTEGRATOR/DM_FUZZY/SIDM).
 * The accessor must match the field gate — NOT include GALSF_MERGER_STARCLUSTER_PARTICLES
 * alone, which doesn't enable AGS_zeta. */
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
static KOKKOS_INLINE_FUNCTION
double gpu_get_ags_zeta(const struct particle_data *Pp, int p)
{
    return Pp[p].AGS_zeta;
}
#endif

/* Permanent invariant guard (always enforced).  A tree-node multipole must never stand in
 * for a particle leaf where any enabled leaf interaction (AGS softening/zeta, ...) distinguishes
 * them.  In the one-shot LET that risk is a foreign TERMINAL node the shared opening predicate wanted
 * to OPEN but cannot descend (nextnode<0): if it is a tagged real single-particle leaf it is routed
 * through particle-leaf secondary semantics below (legal); otherwise it is an unopenable aggregate
 * that would be silently downgraded to a multipole, which is ILLEGAL until an owner-continuation path
 * exists.  The host hard-surfaces (controlled stop) when g_inv_fterm_aggregate>0 after the walk. */
/* Device-written (Kokkos::atomic_add inside the walk) + host-read (report/reset/controlled-stop).
 * On a GPU compiler the storage MUST be device-addressable, so use the same `__managed__` idiom as
 * gpu_device_error_sentinel.h; on the Mac OpenMP build (no GPU compiler) a plain host static suffices
 * (the Kokkos lambda runs on host threads).  Plain `static long long` here built on Mac but failed the
 * Vista nvcc compile (#20096-D: address of a host variable in device code). */
#if defined(GIZMO_GPU_COMPILER)
static __managed__ long long g_inv_fterm_aggregate = 0;   /* predicate-OPEN foreign terminal, NOT a leaf -> illegal */
static __managed__ long long g_unship_aggregate = 0;      /* ... of which the sender could never have shipped */
#else
static long long g_inv_fterm_aggregate = 0;
static long long g_unship_aggregate = 0;
#endif

/* The device replicas of weight_function_for_weighted_motion_smoothing and
 * ags_gravity_kernel_shared_BITFLAG were collapsed into gravtree_force_kernel.h
 * (grav_weight_function_for_weighted_motion_smoothing / gravtree_ags_kernel_shared_bitflag),
 * shared verbatim with the CPU walk. */

/* RT payload data passed to the GPU walk kernel.
 * src_lum[p * N_RT_FREQ_BINS + kf] = per-particle luminosity precomputed on
 * CPU via rt_get_source_luminosity().  Only populated when RT_USE_GRAVTREE is
 * active.  Sized for [NumPart * N_RT_FREQ_BINS] in SharedSpace. */
#ifdef RT_USE_GRAVTREE
#include "../radiation/rt_functions.h"
struct gpu_rt_walk_data_t {
    MyFloat *src_lum;             /* [NumPart * N_RT_FREQ_BINS] */
#ifdef CHIMES_STELLAR_FLUXES
    double  *src_lum_G0;          /* [NumPart * CHIMES_LOCAL_UV_NBINS] */
    double  *src_lum_ion;         /* [NumPart * CHIMES_LOCAL_UV_NBINS] */
#endif
};
#endif

/* Sink radiation payload.  Pre-computed on CPU before the kernel
 * launches because sink_lum_bol() (and, under SINGLE_STAR_SINK_DYNAMICS,
 * calculate_individual_stellar_luminosity()) are not GPU-callable.
 * bh_lum[p]   = sink_lum_bol(P[p].Sink_Mdot, P[p].Sink_Mass, p) when P[p]
 *               is a valid type-5 sink with Mdot>0, else 0.
 * bh_angle[p] = P[p].Sink_Specific_AngMom (if SINK_FOLLOW_ACCRETED_ANGMOM)
 *               or P[p].GradRho otherwise.  Used for angle-weighted
 *               luminosity at particle leafs.  Node-level sink_lum /
 *               sink_lum_grad already live in the SoA. */
/* COSMIC_RAY_SUBGRID_LEBRON payload.  cr_get_source_injection_rate
 * is not GPU-callable so per-particle injection is precomputed on host.
 * t_max_cr = DMIN(1., evaluate_time_since_t_initial_in_Gyr(All.TimeBegin))/
 * UNIT_TIME_IN_GYR, passed as scalar (independent of target). */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
struct gpu_cr_walk_data_t {
    MyFloat *cr_inject;       /* [NumPart] */
    double   t_max_cr;        /* scalar, in code time units */
};
#endif

#ifdef SINK_PHOTONMOMENTUM
struct gpu_sink_walk_data_t {
    MyFloat       *bh_lum;    /* [NumPart] */
    Vec3<MyFloat> *bh_angle;  /* [NumPart] */
};

/* The device replica of sink_fb_angleweight was collapsed into gravtree_force_kernel.h
 * (grav_sink_fb_angleweight, component args), shared verbatim with the host function. */
#endif /* SINK_PHOTONMOMENTUM */

/* ---- Device-callable source cores for the lazy per-open source evaluation ----
 * This TU opens GRAVTREE_SOURCE_DEVICE_TU, so gravtree_moment_sources.h's fill helper is
 * KOKKOS_INLINE here and its RT/sink/CR calls must see the inline device bodies. Pull the
 * enabled families' *_functions.h cores in FIRST (rt_functions.h is already included above
 * under RT_USE_GRAVTREE and transitively carries the stellar-evolution/cosmology + sink
 * leaves; sink/CR are added here for the non-RT source configs). Host/default TUs compile
 * the helper as static-inline against the proto.h host wrappers and skip this block. */
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
#ifdef SINK_PHOTONMOMENTUM
#include "../sinks/sink_functions.h"                       /* sink_lum_bol_core */
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
#include "../eos/cosmic_ray_fluid/cosmic_ray_functions.h"  /* cr_get_source_injection_rate */
#endif
#endif
#include "gravtree_moment_sources.h" /* SSOT per-particle source-input fill helper (static-inline host / KOKKOS_INLINE device) */

/* Ewald periodic-image POTENTIAL correction for the primary walk.  Under
 * pure-tree periodic gravity with EVALPOTENTIAL the CPU walk (forcetree.cc)
 * adds mass*ewald_pot_corr(dr) to the potential of every accepted interaction;
 * the GPU primary walk previously added only the short-range potential and left
 * the periodic-image term out (the seeded g_d_potcorr table was never read).
 * This POD carries the device mirror of that table + its interpolation scale
 * into the primary walk.  It is passed UNCONDITIONALLY (one struct, optional
 * fields gated once here) rather than as a stacked-#ifdef parameter; the walk
 * reads it only inside the matching compile gate.  In a healthy run 'active' is
 * always 1: an acquire failure hard-stops the primary walk (the build requires
 * the correction).  'active' is 0 only as a NULL-guard for the graceful drain
 * that follows that endrun. */
struct gpu_ewald_pot_data_t {
    const MyFloat *potcorr;   /* flat [(EN+1)^3] Ewald potential-correction table, or NULL */
    double         fac_intp;  /* table interpolation scale (= g_ewald_fac_intp) */
    int            active;    /* 1 iff potcorr is a valid acquired table */
};

/* Hermite pass state for the device walk. HermiteOnlyFlag and TimeBinActive[] are host
 * globals, so the walk cannot read them; the host snapshots them into a few words and
 * passes them by value (the prediction reads the walk's own drift/kick table view).
 * No device allocation and no per-particle array -- this is the few-body regime, where a
 * pass over the particles would cost more than the walk it serves. */
#ifdef HERMITE_INTEGRATION
struct gpu_hermite_walk_data_t {
    struct HermiteWalkState   state;
};
#endif

/* A source the walk found standing AHEAD of the walk time.  No drift produces that state, so it is
 * not a case to route around: the walk reads the source as stored, counts it and keeps the first
 * offender, and the host stops the run once the kernel has finished. */
struct gpu_grav_time_fault_t {
    int n_particles, n_nodes;   /* sources found ahead, counted at every read */
    int claimed;                /* the first reader to raise this writes the offender below */
    int first_no;
    integertime first_ti;
};

static KOKKOS_INLINE_FUNCTION void
gpu_grav_note_source_ahead(struct gpu_grav_time_fault_t *fault, int is_node, int no, integertime ti_source)
{
    Kokkos::atomic_add(is_node ? &fault->n_nodes : &fault->n_particles, 1);
    if(Kokkos::atomic_fetch_add(&fault->claimed, 1) == 0) {fault->first_no = no; fault->first_ti = ti_source;}
}

/* Host, once the walk's kernels have finished: a source found ahead of the walk time means this pass
 * read a state no drift produces, so its forces are not valid.  Report the first offender and stop. */
static void gpu_grav_report_time_fault(const char *walk, const struct gpu_grav_time_fault_t *fault, integertime ti, int errorcode)
{
    if(fault->n_particles == 0 && fault->n_nodes == 0) {return;}
    printf("%s: task %d: %d particle and %d node reads found a source AHEAD of the walk time %lld (first: index %d at %lld); no drift produces that state\n",
           walk, ThisTask, fault->n_particles, fault->n_nodes, (long long) ti, fault->first_no, (long long) fault->first_ti);
    fflush(stdout);
    endrun(errorcode);
}

/* This translation unit's SharedSpace mirror of the drift/gravkick tables, allocated on first
 * use and reused.  Every device walk here reads it, to predict a source to the walk time and for
 * the Hermite source prediction.  One storage, because a device symbol cannot be shared across
 * translation units without relocatable device code and a second copy of the same 16 KB table
 * would have to be kept in step with this one.  Stays NULL on a non-cosmological run: the refresh returns an
 * elapsed-time view that reads no table at all. */
static double *tu_drift_kick_table_dev = NULL;

extern "C" void gpu_gravtree_tables_release(void)
{
    if(tu_drift_kick_table_dev) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(tu_drift_kick_table_dev);
        tu_drift_kick_table_dev = NULL;
    }
}

/* Host: acquire the Ewald tables (idempotent) and fill the potential POD.
 * Returns 0 on success (out->active=1), nonzero if the tables are not ready
 * (out->active=0); the caller treats a nonzero return as a hard stop. */
#ifdef GIZMO_GPU_EWALD_POT_CORRECTION
static int gpu_ewald_acquire_pot_data(struct gpu_ewald_pot_data_t *out);
#endif

/* -------------------------------------------------------------------------
 * Compile-time payload gates.
 * Everything not yet ported remains #error'd so wrong-physics on those
 * configs is caught at compile time rather than producing silent incorrect
 * results.
 * ---------------------------------------------------------------------- */
/* ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION: ported.
 * - Softening lookup: handled by P[i].ForceSoftening cache (single source of truth in
 *   compute_force_softening_kernel_radius()). Sphere-box opening criterion (mirrors
 *   forcetree.cc:2122-2130) handles NEIGHBORS_MUST_BE_COMPUTED auto-defined by this flag.
 * - Tree-node previous-step tidal tensor: SoA tidal_tensorps field (gpu_pseudo_update +
 *   gpu_moment_refresh + let_pack already populate it).
 * - Walk accumulators: tidal_zeta (scalar) + per-pair acc_corr_zeta correction (folded
 *   into acc).  Mirrors forcetree.cc:2481-2526. */
/* SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM: ported via the P[i].ForceSoftening cache.
 * The type-4 mass-based softening override lives in compute_force_softening_kernel_radius()
 * (forcetree.cc:67-69) and is read on GPU as Pp[p].ForceSoftening — single source of truth. */
/* GALSF_MERGER_STARCLUSTER_PARTICLES (type-4 star cluster softening via
 * StarParticleEffectiveSize) and ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT (type-0 softening
 * cap): both live in compute_force_softening_kernel_radius(), and the walk reads the
 * result from the same P[i].ForceSoftening cache. */
/* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE + COMPUTE_JERK_IN_GRAVTREE: ported (ATFU).
 * Tidal tensor accumulation + jerk both mirror forcetree.cc:2081-2292.
 * All sub-cases are ported; see the entries below for details. */
/* COMPUTE_TIDAL_TENSOR + PMGRID: ported.  shortrange_table_tidal is mirrored
 * to SharedSpace once per run by gpu_shortrange_tables_acquire() and consumed
 * in the tidal accumulation block (mirrors forcetree.cc:2538-2549). */
/* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE + ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING: ported.
 * The averaging branch in the inside-softening force-kernel section now sets fac_tidal
 * and averages fac2_tidal alongside fac_accel/fac_pot (mirrors forcetree.cc:2393-2394). */
/* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE + ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION: ported via
 * the same walk-side block that handles ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION above. */
/* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE + GRAVITY_SPHERICAL_SYMMETRY: ported.
 * The shell-theorem fac2_tidal override at the top of the tidal accumulation
 * block handles this case alongside the standard non-spherical formula. */
/* SINK_PHOTONMOMENTUM, SINK_COMPTON_HEATING, SINK_DYNFRICTION_FROMTREE:
 * ported. */
/* SPECIAL_POINT_MOTION + SPECIAL_POINT_WEIGHTED_MOTION: ported.
 * Walk-side accumulation of nearest-special-particle vel/acc lives in the
 * SINK_CALC_DISTANCES branches (leaf-particle and node paths). The Acc_Total_PrevStep
 * field is a member of particle_data and is automatically mirrored in P_dev.
 * Tree-node sink_acc lives in the GravitySoA + populated by gpu_pseudo_update +
 * gpu_moment_refresh.  The weighted variant uses the shared
 * grav_weight_function_for_weighted_motion_smoothing() (gravtree_force_kernel.h). */
/* SINK_CALC_DISTANCES, SINGLE_STAR_SINK_DYNAMICS, SINGLE_STAR_TIMESTEPPING,
 * SINGLE_STAR_FIND_BINARIES, SINGLE_STAR_FB_TIMESTEPLIMIT, SINGLE_STAR_STARFORGE_DEFAULTS:
 * ported. */
/* COSMIC_RAY_SUBGRID_LEBRON: ported. Reads SoA cr_injection at
 * node accepts + per-particle precomputed d_cr_inject at leaf nodes (since
 * cr_get_source_injection_rate is not GPU-callable). */
/* COUNT_MASS_IN_GRAVTREE: ported. tree_mass accumulator declared at function entry,
 * accumulated once per ACCEPTED interaction in the force kernel (r2>0, mass>0;
 * mirrors forcetree.cc -- excludes the target's own leaf), written to
 * P_dev[target].TreeMass at end of walk, scattered back to P[i].TreeMass.  The
 * post-loop +=P[i].Mass in gravtree.cc adds the target's own mass to finalize. */
/* DM_SCALARFIELD_SCREENING: ported. SoA tree-node fields mass_dm + s_dm are populated by
 * gpu_pseudo_update + gpu_moment_refresh + let_pack (already wired). The walk sets per-
 * interaction d_dm and mass_dm_local in both leaf and node branches, then accumulates the
 * Yukawa-screened scalar-field force on non-gas targets after the main force kernel. */
/* GRAVITY_SPHERICAL_SYMMETRY: ported. Box-center sph_center + r_target are computed
 * once at function entry; r_source is set per-interaction (leaf and node branches)
 * to the source distance from the box center. The shell-theorem force law overrides
 * fac_accel + dr right before the acc accumulation, and fac2_tidal at the start of
 * the tidal block (mirrors forcetree.cc:2446-2449 and :2534-2536). */
/* HERMITE_INTEGRATION + NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE: ported.
 * NEIGHBORS_MUST_BE_COMPUTED activates the sphere-box intersection opening criterion
 * in the walk loop above (mirrors forcetree.cc:2122-2130). HERMITE_INTEGRATION
 * additionally requires COMPUTE_JERK_IN_GRAVTREE which is auto-defined and already
 * ported: GravJerk written at line ~1155 and scattered back at line ~1478,
 * read by the Hermite predictor in core/kicks.cc:212 and gravtree.cc:571. */
/* ADAPTIVE_TREEFORCE_UPDATE: pre-walk skip-flag filtering in the
 * dispatcher + jerk accumulation in the walk kernel.  Skip-flag particles
 * (needs_new_treeforce()==0) bypass the GPU walk and use the CPU post-loop
 * jerk extrapolation path at gravtree.cc:512-520 (GravAccel += GravJerk*dt). */
/* Periodic boundary handling:
 *   BOX_PERIODIC + PMGRID   → TreePM.  Long-range forces come from PM; the tree
 *                             walk is short-range-only (rcut-truncated via the
 *                             shortrange-force tables already wired in).  The
 *                             CPU never calls force_treeevaluate_ewald_correction
 *                             in this case (see gravtree.cc:734 gate).  The GPU
 *                             primary walk is therefore already correct.
 *   BOX_PERIODIC + !PMGRID  → pure-tree periodic.  Requires the second Ewald-
 *                             correction walk; the GPU port of that walk lives
 *                             in gpu_ewald_walk_primary (dispatched from
 *                             gravtree.cc after the primary walk).
 *   GRAVITY_NOT_PERIODIC    → non-periodic box, Ewald not relevant. */
/* Pure-tree periodic gravity (BOX_PERIODIC && !GRAVITY_NOT_PERIODIC && !PMGRID)
 * is supported via the GPU Ewald walk (gpu_ewald_walk_primary), dispatched
 * from gravtree.cc on the Ewald_iter==1 pass. Implementation later in this
 * file. */
/* SELFGRAVITY_OFF: gravtree.cc wraps the entire tree dispatch (including GPU
 * dispatch) in #ifndef SELFGRAVITY_OFF, so this TU is never compiled into the
 * walk when gravity is disabled.  No guard needed here. */


/* -------------------------------------------------------------------------
 * The device walk, in units.
 *
 * A walk has a per-call context (the tree and the tables, read-only), a member
 * (one target: its fixed inputs to the opening decision and the pair evaluation,
 * and everything it accumulates), and per-element steps: a member-independent
 * node prelude, a per-member decision, and a per-member evaluation of an accepted
 * element through the SoA / P_dev adapter and the shared pair seam
 * (grav_pair_evaluate_core, gravtree_force_kernel.h -- the CPU walk uses the
 * same seam). gpu_gravtree_walk_one composes them for one target, evaluating each
 * accepted element at encounter; a packet walk composes the same units for several
 * targets sharing one traversal, recording elements and evaluating later. Every
 * unit has one call site per walk so the compiler keeps a member in registers.
 * ---------------------------------------------------------------------- */

/* Per-call inputs, captured by value into the kernel. */
struct gpu_grav_walk_ctx_t {
    int treeBase, treeParticleSlots, maxNodes, maxForeignNodes;   /* pseudos start at treeBase+maxNodes+maxForeignNodes */
    struct particle_data *P_dev;
    struct gas_cell_data *CellP_dev;
    integertime ti;   /* the walk time: every source is read as it stands at this time, predicted there if it is behind */
    struct gpu_gravity_tree_soa_t tree_soa;   /* the SoA handle set, by value (a struct of pointers into SharedSpace) */
    struct extNODE *Extnodes;   /* the node arrays (shifted like Nodes), read only for the pending momenta of a kicked node */
    struct DriftKickTableView tables;   /* the drift/kick factors a prediction reads */
    struct gpu_grav_time_fault_t *time_fault;   /* sources found ahead of the walk time (see its definition) */
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    int is_first_step;   /* hybrid opening: relative criterion applies only after step 0 */
#endif
    grav_pm_shortrange_t pm;   /* PM short-range config (empty when !PMGRID); a member overrides its copy under PM_PLACEHIGHRESREGION */
#ifdef RT_USE_GRAVTREE
    struct gpu_rt_walk_data_t rt_data;
#endif
#ifdef SINK_PHOTONMOMENTUM
    struct gpu_sink_walk_data_t sink_data;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    struct gpu_cr_walk_data_t cr_data;
#endif
    bool use_lazy_source;   /* evaluate the source payload on-device at each local leaf (no dense eager arrays) */
#ifdef HERMITE_INTEGRATION
    struct gpu_hermite_walk_data_t hermite;   /* Hermite pass state, by value (a few words) */
#endif
    struct gpu_ewald_pot_data_t ewald_pot;   /* periodic-image potential correction (read only under the pure-tree-periodic EVALPOTENTIAL gate) */
};

/* One target of a walk: its inputs to the opening decision and to the pair
 * evaluation, fixed for the walk, and everything the walk accumulates for it.
 * Mirrors the CPU walk's member (forcetree.cc). */
/* What the opening decision reads about a target at every node. A packet walk keeps a
 * copy of every member's in team scratch, where every walker reads them. */
struct gpu_grav_open_inputs_t {
    Vec3<double> pos; int ptype; double soft, aold;
#ifdef PMGRID
    double rcut, rcut2;
#endif
    int alive;   /* 0 for a massless target, which takes part in nothing */
};

/* A member is split in two.  The INPUTS are fixed for the walk: what the opening decision and the
 * pair evaluation read about the target.  The SUMS are everything the walk accumulates for it.  A
 * packet team holds one copy of a member's inputs and gives each of the lanes working on that member
 * its own sums, folded together when the packet commits; every field of the sums is listed, with the
 * rule that combines it, in gpu_grav_member_sums_combine below. */
struct gpu_grav_member_inputs_t {
    int target;
    struct gpu_grav_open_inputs_t open;
    /* pair evaluation inputs */
    double pmass, zeta;
    grav_pair_tgt_t tgt;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    Vec3<double> vel;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double sph_center[3];
#endif
#ifdef RT_USE_GRAVTREE
    int valid_gas_particle_for_rt;   /* read into a volatile local at each use: nvc++ constant-propagates a plain gate inside the walk loop */
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    double fac_stellum[N_RT_FREQ_BINS];
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    int cr_active_gate;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double r_for_total_menclosed;
#endif
};

struct gpu_grav_member_sums_t {
    grav_pair_acc_t out;
#ifdef COUNT_MASS_IN_GRAVTREE
    double tree_mass;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double m_enc_in_rcrit;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    double treecol_angular_bins[RT_USE_TREECOL_FOR_NH];
#endif
#ifdef SINK_COMPTON_HEATING
    double incident_flux_agn;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double SubGrid_CosmicRayEnergyDensity;
#endif
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_flux_G0[CHIMES_LOCAL_UV_NBINS], chimes_flux_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#ifdef RT_OTVET
    SymmetricTensor2<double> RT_ET[N_RT_FREQ_BINS];
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    double incident_flux_uv, incident_flux_euv;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    double Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    Vec3<double> Rad_Flux[N_RT_FREQ_BINS];
#endif
#ifdef SINK_CALC_DISTANCES
    grav_sink_prox_accum_t sink_prox;   /* nearest-sink + single-star timestep/binary accumulators (shared with the CPU walk) */
#endif
};

/* Every field of the sums, with the rule that combines two partial sums of the same member into one:
 * `a` takes in `b`.  THE single list for this member: a field added to gpu_grav_member_sums_t is added
 * here in the same edit, or the packet engine's fold silently drops it (the pair and sink-proximity
 * accumulators carry their own rules beside their definitions).  Every rule is a sum, a minimum, or a
 * minimum carrying the fields that describe what attained it; on equal minima `a` keeps its own.
 * target_ptype is the member's type, which the sink-proximity rules need. */
KOKKOS_INLINE_FUNCTION void gpu_grav_member_sums_combine(gpu_grav_member_sums_t &a, const gpu_grav_member_sums_t &b, int target_ptype)
{
    grav_pair_acc_combine(a.out, b.out);
#ifdef COUNT_MASS_IN_GRAVTREE
    a.tree_mass += b.tree_mass;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    a.m_enc_in_rcrit += b.m_enc_in_rcrit;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    for(int k = 0; k < RT_USE_TREECOL_FOR_NH; k++) {a.treecol_angular_bins[k] += b.treecol_angular_bins[k];}
#endif
#ifdef SINK_COMPTON_HEATING
    a.incident_flux_agn += b.incident_flux_agn;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    a.SubGrid_CosmicRayEnergyDensity += b.SubGrid_CosmicRayEnergyDensity;
#endif
#ifdef CHIMES_STELLAR_FLUXES
    for(int k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) {a.chimes_flux_G0[k] += b.chimes_flux_G0[k]; a.chimes_flux_ion[k] += b.chimes_flux_ion[k];}
#endif
#ifdef RT_OTVET
    for(int f = 0; f < N_RT_FREQ_BINS; f++) {for(int k = 0; k < 6; k++) {a.RT_ET[f].data[k] += b.RT_ET[f].data[k];}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    a.incident_flux_uv += b.incident_flux_uv; a.incident_flux_euv += b.incident_flux_euv;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    for(int f = 0; f < N_RT_FREQ_BINS; f++) {a.Rad_E_gamma[f] += b.Rad_E_gamma[f];}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    for(int f = 0; f < N_RT_FREQ_BINS; f++) {for(int k = 0; k < 3; k++) {a.Rad_Flux[f][k] += b.Rad_Flux[f][k];}}
#endif
#ifdef SINK_CALC_DISTANCES
    grav_sink_prox_accum_combine(a.sink_prox, b.sink_prox, target_ptype);
#endif
    (void) target_ptype;
}

/* The single-target walk holds both halves itself. */
struct gpu_grav_member_t {
    gpu_grav_member_inputs_t in;
    gpu_grav_member_sums_t sums;
};

/* What an accepted element carries from its load to the shared evaluation, beyond
 * the pair inputs in grav_pair_src_t: the payload values the walker-local blocks
 * consume. Set on every path that reaches the evaluation. */
struct gpu_grav_src_payload_t {
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    double gasmass;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<double> d_dm;   /* displacement to the DM mass centre (= the total CoM only for a pure-DM leaf) */
    double mass_dm_local;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
    double mass_stellarlum[N_RT_FREQ_BINS];
    Vec3<double> d_stellarlum;
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_mass_stellarlum_G0[CHIMES_LOCAL_UV_NBINS];
    double chimes_mass_stellarlum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#ifdef SINK_PHOTONMOMENTUM
    double mass_sinklumwt_forradfb;
#endif
#endif
};

/* Everything one target contributes to an opening decision, and nothing it contributes to a
 * pair evaluation.  `zeta` and the per-target PM override come out with the inputs because the
 * member prologue below needs them in the same breath.
 * Returns 0 for a massless target, which takes part in nothing. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_open_inputs_init(const gpu_grav_walk_ctx_t &ctx, int target,
                          gpu_grav_open_inputs_t &open, double &pmass, double &zeta,
                          grav_pm_shortrange_t &pm)
{
    struct particle_data *P_dev = ctx.P_dev;
    open.pos = P_dev[target].Pos;
    open.ptype = P_dev[target].Type;
    open.alive = 0;
    pmass = P_dev[target].Mass;
    zeta = 0.0;
    pm = ctx.pm;
    if(pmass <= 0) {return 0;}
    open.alive = 1;
    const int ptype = open.ptype;

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(GALSF_MERGER_STARCLUSTER_PARTICLES)
    double soft = gpu_force_softening_kernel_radius(P_dev, target);
#else
    double soft = All.ForceSoftening[ptype];
#endif
    /* zeta is set unconditionally (matches CPU walk); passed to the shared pair kernel,
     * consumed there only under #if AGS */
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    grav_target_select_soft_and_zeta(ptype, gpu_get_ags_zeta(P_dev, target), soft, zeta);
#endif
    open.soft = soft;
    open.aold = All.ErrTolForceAcc * P_dev[target].OldAcc;

#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
    /* high-res zoom particles use the finer short-range PM cutoff (mirrors forcetree.cc target
     * prologue). The dispatcher passes the coarse-mesh rcut/asmthfac; override per target here. */
    if(pmforce_is_particle_high_res(ptype, open.pos)) {
        pm.rcut = All.Rcut[1]; pm.rcut2 = pm.rcut * pm.rcut; pm.asmthfac = grav_pm_asmthfac(All.Asmth[1]);
    }
#endif
#ifdef PMGRID
    open.rcut = pm.rcut; open.rcut2 = pm.rcut2;
#endif
    return 1;
}

/* Zero one member's sums: every field starts at the identity of the rule that folds it
 * (gpu_grav_member_sums_combine), so a partial that never evaluates anything folds in as nothing. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_member_sums_init(gpu_grav_member_sums_t &sums)
{
    grav_pair_acc_init(sums.out);
#ifdef COUNT_MASS_IN_GRAVTREE
    /* Diagnostic: total mass seen by this target during the walk, summed only
     * over accepted interactions (mirrors forcetree.cc). The walk excludes the
     * target's own leaf (r2==0); the post-loop +=P[i].Mass in gravtree.cc
     * finalizes the sum. */
    sums.tree_mass = 0.0;
#endif
#ifdef SINK_COMPTON_HEATING
    sums.incident_flux_agn = 0.0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    sums.SubGrid_CosmicRayEnergyDensity = 0.0;
#endif
#ifdef SINK_CALC_DISTANCES
    grav_sink_prox_accum_init(sums.sink_prox);
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    {int kb; for(kb=0; kb<RT_USE_TREECOL_FOR_NH; kb++) {sums.treecol_angular_bins[kb]=0.0;}}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    sums.m_enc_in_rcrit = 0.0;
#endif
#ifdef CHIMES_STELLAR_FLUXES
    {int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {sums.chimes_flux_G0[kc]=0; sums.chimes_flux_ion[kc]=0;}}
#endif
#ifdef RT_OTVET
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {sums.RT_ET[kf] = {};}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    sums.incident_flux_uv = 0.0; sums.incident_flux_euv = 0.0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {sums.Rad_E_gamma[kf]=0.0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {sums.Rad_Flux[kf]={};}}
#endif
}

/* Set up one target's inputs for the walk (the CPU walk's target prologue).
 * Returns 0 for a massless target, which takes part in nothing and writes zeros. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_member_inputs_init(const gpu_grav_walk_ctx_t &ctx, int target, gpu_grav_member_inputs_t &in)
{
    struct particle_data *P_dev = ctx.P_dev;
    in.target = target;
    double zeta = 0.0; grav_pm_shortrange_t pm;
    if(!gpu_grav_open_inputs_init(ctx, target, in.open, in.pmass, zeta, pm)) {return 0;}
    in.zeta = zeta;
    const int ptype = in.open.ptype; const double pmass = in.pmass;
    const double soft = in.open.soft;

    /* fed unconditionally to the shared pair kernel (consumed there only under the
     * symmetrize-by-averaging #if); matches the CPU walk's unconditional precompute. */
    const int ags_bitflag_primary = gravtree_ags_kernel_shared_bitflag(ptype);

#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    in.vel = P_dev[target].Vel;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    /* Shell-theorem gravity: forces from any source at r_source > r_target
     * vanish; forces from r_source < r_target use a 1/r^3 enclosed-mass formula
     * pointed toward the box center. Mirrors forcetree.cc. */
    in.sph_center[0] = 0.0; in.sph_center[1] = 0.0; in.sph_center[2] = 0.0;
#ifdef BOX_PERIODIC
    in.sph_center[0] = 0.5 * boxSize_X;
    in.sph_center[1] = 0.5 * boxSize_Y;
    in.sph_center[2] = 0.5 * boxSize_Z;
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    /* per-target CR gate (the host precompute leaves t_max_cr=0 unless All.Time>All.TimeBegin,
     * mirroring the CPU walk's gate) */
    in.cr_active_gate = (ctx.cr_data.t_max_cr > 0) ? 1 : 0;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    in.r_for_total_menclosed = grav_target_menc_radius(soft); /* baseline Rcrit_min applied in the helper */
#endif
#ifdef RT_USE_GRAVTREE
    /* valid-gas RT gate via the shared helper */
    in.valid_gas_particle_for_rt = grav_target_valid_gas_for_rt(ptype, soft, pmass);
#endif /* RT_USE_GRAVTREE */

    /* RT_LEBRON radiation-pressure coupling factor (forcetree.cc). Once per target, before
     * the walk, because it only depends on the target's properties. Skipped when save-flux
     * mode is active (flux is stored and converted to RP after the walk by the caller). */
#if defined(RT_USE_GRAVTREE) && defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {in.fac_stellum[kf]=0.0;}}
    {
        volatile int valid_gas_particle_for_rt = in.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt) {
            double kappa_eff[N_RT_FREQ_BINS]; int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {kappa_eff[kf] = rt_kappa(-1, kf, P_dev, ctx.CellP_dev);}
            grav_target_rt_fac_stellum(soft, pmass, kappa_eff, in.fac_stellum);
        }
    }
#endif

    /* the target's inputs to the shared pair evaluation, fixed for this walk */
    grav_pair_tgt_t tgt{}; tgt.ptype = ptype; tgt.pmass = pmass; tgt.h = soft; tgt.zeta = zeta; tgt.ags_bitflag = ags_bitflag_primary; tgt.pm = pm;
#ifdef SINK_DYNFRICTION_FROMTREE
    tgt.sink_mass = (ptype == 5) ? P_dev[target].Sink_Mass : 0.0;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    {
        SymmetricTensor2<MyFloat> tmp = P_dev[target].tidal_tensorps_prevstep;
        for(int kk = 0; kk < 6; kk++) tgt.i_zeta_tidal_tensorps_prevstep.data[kk] = (double) tmp.data[kk];
    }
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    tgt.pos = in.open.pos; tgt.center[0] = in.sph_center[0]; tgt.center[1] = in.sph_center[1]; tgt.center[2] = in.sph_center[2];
#endif
    in.tgt = tgt;
    return 1;
}

/* Set up one target's member for the single-target walk: its inputs and its zeroed sums. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_member_init(const gpu_grav_walk_ctx_t &ctx, int target, gpu_grav_member_t &mem)
{
    gpu_grav_member_sums_init(mem.sums);
    return gpu_grav_member_inputs_init(ctx, target, mem.in);
}

/* A particle source as a walk at `ti` reads it, with nothing written back: position, velocity and
 * mass predicted to the walk time by the step bodies the drift itself runs (predict_particle_motion,
 * core/timestep_functions.h); a source already at the walk time is read as stored. A source AHEAD of
 * the walk time, which no drift produces, is read as stored and reported. */
struct gpu_grav_particle_now_t { Vec3<double> pos, vel; double mass; };   /* what a walk keeps of a particle source */
static KOKKOS_INLINE_FUNCTION gpu_grav_particle_now_t
gpu_grav_particle_source_at(struct particle_data *P_dev, struct gas_cell_data *CellP_dev, int no, integertime ti,
                            const struct DriftKickTableView &tables, struct gpu_grav_time_fault_t *time_fault)
{
    gpu_grav_particle_now_t out;
    const integertime ti_source = P_dev[no].Ti_current;
    if(ti_source >= ti) {   /* read as stored, without building the predictor (which also loads the cell under MFV) */
        if(ti_source > ti) {gpu_grav_note_source_ahead(time_fault, 0, no, ti_source);}
        out.pos = P_dev[no].Pos; out.vel = P_dev[no].Vel; out.mass = P_dev[no].Mass;
        return out;
    }
    struct particle_motion_prediction motion(P_dev, CellP_dev, no);
    predict_particle_motion(motion, ti, &tables);
    out.pos = motion.pos_; out.vel = motion.vel_; out.mass = motion.mass_;
    return out;
}

/* A node's motion as a walk at `ti` reads it, with nothing written back: the state the drift would
 * leave it in. A node at the walk time is read as stored, its pending kick left pending (a drift folds
 * a kick only when it moves the node forward). A node behind it has its kick folded and its centres
 * and length advanced by the node-motion arithmetic the drift runs (gravtree_moment_kernel.h), on the
 * same two clocks and with the dilation taken at the centre before the advance, as the drift takes
 * them. Everything is read from the mirror except the pending momenta, which live only in the node
 * arrays and are read there only when the mirror says the node was kicked. Returns 1 for a node AHEAD
 * of the walk time, which no drift produces: it is left as stored for the caller to report. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_node_motion_at(const struct gpu_gravity_tree_soa_t &soa, const struct extNODE *ext, int idx, int no,
                        integertime ti, const struct DriftKickTableView &tables, node_motion_prediction &m)
{
    for(int j = 0; j < 3; j++) {m.s_[j] = (MyFloat) soa.s[idx][j]; m.vs_[j] = (MyFloat) soa.node_vs[idx][j]; m.dp_[j] = 0;}
    m.len_ = soa.len[idx]; m.mass_ = (double) soa.mass[idx]; m.vmax_ = (double) soa.vmax[idx];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    for(int j = 0; j < 3; j++) {m.rt_s_[j] = (MyFloat) soa.rt_source_lum_s[idx][j]; m.rt_vs_[j] = (MyFloat) soa.rt_source_lum_vs[idx][j]; m.rt_dp_[j] = 0;}
    m.lum_tot_ = 0;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    for(int j = 0; j < 3; j++) {m.s_dm_[j] = (MyFloat) soa.s_dm[idx][j]; m.vs_dm_[j] = (MyFloat) soa.vs_dm[idx][j]; m.dp_dm_[j] = 0;}
    m.mass_dm_ = (double) soa.mass_dm[idx];
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    for(int j = 0; j < 3; j++) {m.sink_pos_[j] = (MyFloat) soa.sink_pos[idx][j]; m.sink_vel_[j] = (MyFloat) soa.sink_vel[idx][j]; m.sink_dp_[j] = 0;}
    m.sink_mass_ = (double) soa.sink_mass[idx];
#endif
    const integertime node_ti = soa.node_ti[idx];
    if(node_ti > ti) {return 1;}
    if(node_ti == ti) {return 0;}

    if(soa.bitflags[idx] & (1u << BITFLAG_NODEHASBEENKICKED))
    {
        for(int j = 0; j < 3; j++) {m.dp_[j] = ext[no].dp[j];}
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        for(int j = 0; j < 3; j++) {m.rt_dp_[j] = ext[no].rt_source_lum_dp[j];}
        {double l = 0; for(int b = 0; b < N_RT_FREQ_BINS; b++) {l += (double) soa.stellar_lum[(long) idx * N_RT_FREQ_BINS + b];} m.lum_tot_ = l;}
#endif
#ifdef DM_SCALARFIELD_SCREENING
        for(int j = 0; j < 3; j++) {m.dp_dm_[j] = ext[no].dp_dm[j];}
#endif
#ifdef SINK_NODE_MOTION_TRACKED
        for(int j = 0; j < 3; j++) {m.sink_dp_[j] = ext[no].sink_dp[j];}
#endif
        node_motion_fold_kick(m);
    }
    double dilation = 1.0;
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
    dilation = node_timestep_dilation_factor_at(Vec3<double>{(double) m.s_[0], (double) m.s_[1], (double) m.s_[2]});
#endif
    double dt_drift, dt_widen;
    node_motion_intervals(node_ti, ti, dilation, &tables, dt_drift, dt_widen);
    node_motion_advance(m, dt_drift, dt_widen);
    return 0;
}

/* The member-independent view of a tree node: its moments as stored, its geometry at the
 * walk time (gpu_grav_node_motion_at), its import classification, and what the walk does
 * with it before any member decides. */
enum gpu_grav_node_step_t { GPU_GRAV_NODE_SKIP_TO_SIBLING, GPU_GRAV_NODE_DESCEND, GPU_GRAV_NODE_DECIDE };
struct gpu_grav_node_prelude_t {
    int idx;                        /* SoA index no - treeBase */
    Vec3<MyFloat> s_node;           /* centre of mass at the walk time (before any star-target subtraction) */
    MyFloat len_node, msoft_node, mass_node;
    Vec3<MyFloat> center_node;
    node_motion_prediction motion;  /* every moving field at the walk time; s_node and len_node are copied out of it */
    int in_foreign_n;
    int fl_tag, fl_type;            /* foreign-leaf identity (sidecar), or 0 / -1 */
    double fl_zeta, fl_soft;
    grav_node_kind_t node_kind;
    int sibling, nextnode;
};

/* Load a node and decide what is member-independent about it: an empty node is
 * skipped; a single-particle local or openable-foreign node is descended for an
 * exact force; anything else goes to the per-member decision. Mirrors forcetree.cc. */
static KOKKOS_INLINE_FUNCTION gpu_grav_node_step_t
gpu_grav_node_prelude(const gpu_grav_walk_ctx_t &ctx, int no, gpu_grav_node_prelude_t &nd)
{
    const struct gpu_gravity_tree_soa_t *tree_soa = &ctx.tree_soa;
    const int idx = no - ctx.treeBase;
    nd.idx = idx;
    nd.msoft_node = tree_soa->maxsoft[idx];
    nd.mass_node = tree_soa->mass[idx];
    nd.center_node = tree_soa->center[idx];
    nd.sibling = tree_soa->sibling[idx];
    nd.nextnode = tree_soa->nextnode[idx];

    /* LET guard (mirrors forcetree.cc): a foreign node with nextnode < 0 (unreplaced -1
     * sentinel from unpack) must not be opened -- that would exit the walk and skip its force. */
    nd.in_foreign_n = (no >= ctx.treeBase + ctx.maxNodes);

    /* Foreign-leaf identity lookup.  foreign_slot = idx-maxNodes -- the foreign-only sidecar
     * index, EXPLICIT and bounds-checked so it can never be confused with the SoA index idx. */
    nd.fl_tag = 0; nd.fl_type = -1; nd.fl_zeta = 0.0; nd.fl_soft = 0.0;
    if(nd.in_foreign_n && tree_soa->foreign_leaf_tag) {
        int fs = idx - ctx.maxNodes;
        if(fs >= 0 && fs < tree_soa->foreign_leaf_cap) {
            nd.fl_tag  = tree_soa->foreign_leaf_tag[fs];
            nd.fl_type = tree_soa->foreign_leaf_type[fs];
            nd.fl_zeta = (double) tree_soa->foreign_leaf_zeta[fs];
            nd.fl_soft = (double) tree_soa->foreign_leaf_soft[fs];
        }
    }

    /* empty-node skip: a zero-mass node contributes no force -> the sibling. Never
     * converted to a forced multipole; also avoids descending empty nodes. */
    if(nd.mass_node <= 0) {return GPU_GRAV_NODE_SKIP_TO_SIBLING;}

    /* Classify BEFORE anything that could descend, so the wire tag is the one authority on
     * whether this node's children were shipped (mirrors forcetree.cc). */
    nd.node_kind = grav_classify_node(nd.in_foreign_n, nd.fl_tag, nd.nextnode);

    /* single-particle node -> open to its leaf for an exact force.  Only a local node, or a
     * foreign node shipped WITH its children, has anything below it: a foreign node reaching
     * here with the multi-particle bit clear was shipped multipole-only, so its nextnode is
     * the continuation past the subtree and following it would skip the node's mass. */
    if(!(tree_soa->bitflags[idx] & (1 << BITFLAG_MULTIPLEPARTICLES))) {
        if(nd.node_kind == GRAV_NODE_LOCAL || nd.node_kind == GRAV_NODE_FOREIGN_OPENABLE) {return GPU_GRAV_NODE_DESCEND;}
    }

    /* Only a node some member judges needs its geometry, so only such a node is brought to the walk time. */
    if(gpu_grav_node_motion_at(*tree_soa, ctx.Extnodes, idx, no, ctx.ti, ctx.tables, nd.motion)) {
        gpu_grav_note_source_ahead(ctx.time_fault, 1, no, tree_soa->node_ti[idx]);
    }
    nd.s_node = Vec3<MyFloat>{nd.motion.s_[0], nd.motion.s_[1], nd.motion.s_[2]};
    nd.len_node = nd.motion.len_;
    return GPU_GRAV_NODE_DECIDE;
}

/* A node's centre of mass, mass and wrapped separation as ONE member judges and
 * evaluates them. A star target takes no star mass from the tree (star-star pairs are
 * summed exactly in star_direct_gravity_compute(); taking them here too would double
 * them): the tree carries monopoles only, so removing the sinks is exact -- drop their
 * mass and shift the centre of mass to what is left; both terms are on the same clock,
 * since SINK_NODE_MOTION_TRACKED drifts sink_pos with sink_vel exactly as s is drifted
 * with vs. Returns 0 for a pure-star node seen by a star member: nothing is left of it.
 * The decision and the evaluation both derive the geometry through this one function. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_node_member_geometry(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_node_prelude_t &nd, const gpu_grav_open_inputs_t &open,
                              Vec3<MyFloat> &s_node, MyFloat &mass_node, Vec3<double> &dr, double &r2)
{
    s_node = nd.s_node; mass_node = nd.mass_node;
#ifdef SINGLE_STAR_DIRECT_GRAVITY
    if((open.ptype == 5) && (ctx.tree_soa.sink_mass[nd.idx] > 0))
    {
        MyFloat sm = (MyFloat) ctx.tree_soa.sink_mass[nd.idx];
        if(mass_node - sm <= 0) {return 0;} /* pure-star node */
        Vec3<MyFloat> sp = Vec3<MyFloat>{nd.motion.sink_pos_[0], nd.motion.sink_pos_[1], nd.motion.sink_pos_[2]};   /* on the walk's clock (SINGLE_STAR_DIRECT_GRAVITY tracks sink motion) */
        MyFloat mass_nosink = mass_node - sm;
        for(int k = 0; k < 3; k++) {s_node[k] = (s_node[k]*mass_node - sp[k]*sm) / mass_nosink;}
        mass_node = mass_nosink;
    }
#else
    (void) ctx;
#endif
    dr[0] = s_node[0] - open.pos[0];
    dr[1] = s_node[1] - open.pos[1];
    dr[2] = s_node[2] - open.pos[2];
    gravity_box_nearest_image(dr[0], dr[1], dr[2], -1);
    r2 = dr.norm_sq();
    return 1;
}

/* One member's opening decision at a node the prelude handed over, given the geometry
 * it derived: SKIP, OPEN (descend), or ACCEPT (evaluate as a multipole, or with leaf
 * semantics for a tagged foreign leaf). The predicate (gravtree_opening.h) is the single
 * home for the acceptance geometry: PM short-range cull, neighbour sphere-box /
 * softening-open, the angular and relative criteria, the sink-direct gate. The caller
 * owns the wrapped dr/r2 and the foreign LET policy applied here: a terminal foreign
 * node (a tagged single-particle leaf, or an aggregate shipped multipole-only) has no
 * children here, so a predicate OPEN on it means "accept" -- for a truncated aggregate
 * that is an import that no longer covers what this walk asks of it: reported through
 * `note` (0 none, 1 truncated, 2 unshippable) so the caller can count it at once or hold
 * it until its packet succeeds; the host reports the counts once after the walk. */
#define GPU_GRAV_NOTE_NONE 0
#define GPU_GRAV_NOTE_TRUNCATED 1
#define GPU_GRAV_NOTE_UNSHIPPABLE 2
static KOKKOS_INLINE_FUNCTION gravtree_open_t
gpu_grav_node_member_decide(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_node_prelude_t &nd, const gpu_grav_open_inputs_t &open,
                            MyFloat mass_node, double r2, int &note)
{
    note = GPU_GRAV_NOTE_NONE;
    double cen0 = (double)nd.center_node[0] - open.pos[0];
    double cen1 = (double)nd.center_node[1] - open.pos[1];
    double cen2 = (double)nd.center_node[2] - open.pos[2];
#ifdef PMGRID
    double pred_rcut = open.rcut, pred_rcut2 = open.rcut2;
#else
    double pred_rcut = 0.0, pred_rcut2 = 0.0;
#endif
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    int pred_is_first_step = ctx.is_first_step;
#else
    int pred_is_first_step = 0;
#endif
#if (defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)) && defined(SINGLE_STAR_DIRECT_GRAVITY_RADIUS)
    int pred_n_sink = (int)ctx.tree_soa.N_SINK[nd.idx];
#else
    int pred_n_sink = 0;
#endif
    gravtree_open_t pred = gravtree_open_decision_from_distances(
        r2, cen0, cen1, cen2, open.soft, open.soft, open.aold, open.ptype,
        (double)nd.len_node, (double)mass_node, (double)nd.msoft_node,
        pred_rcut, pred_rcut2, pred_n_sink, pred_is_first_step);
    if(pred == GRAV_SKIP_NODE) {return GRAV_SKIP_NODE;}
    if(pred == GRAV_OPEN_NODE && !grav_node_is_terminal(nd.node_kind)) {return GRAV_OPEN_NODE;}
    if(pred == GRAV_OPEN_NODE && nd.node_kind == GRAV_NODE_FOREIGN_TRUNCATED)   {note = GPU_GRAV_NOTE_TRUNCATED;}
    if(pred == GRAV_OPEN_NODE && nd.node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE) {note = GPU_GRAV_NOTE_UNSHIPPABLE;}
    return GRAV_ACCEPT_MULTIPOLE;
}

/* Count an import-completeness note in the shared per-walk ledger (the host reports it
 * once after the walk). A single-target walk counts at once; a packet holds its notes
 * and counts them only when it succeeds. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_note_commit(int n_truncated_or_unshippable, int n_unshippable)
{
    if(n_truncated_or_unshippable) {Kokkos::atomic_add(&g_inv_fterm_aggregate, (long long) n_truncated_or_unshippable);}
    if(n_unshippable) {Kokkos::atomic_add(&g_unship_aggregate, (long long) n_unshippable);}
}

/* Whether a member takes a particle leaf at all: star-star pairs are summed exactly in
 * star_direct_gravity_compute(), so a star member passes a star leaf. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_leaf_member_accepts(const gpu_grav_walk_ctx_t &ctx, int no, const gpu_grav_open_inputs_t &open)
{
#ifdef SINGLE_STAR_DIRECT_GRAVITY
    if((open.ptype == 5) && (ctx.P_dev[no].Type == 5)) {return 0;}
#else
    (void) ctx; (void) no; (void) open;
#endif
    return 1;
}

/* Evaluate one loaded element for a member: the shared pair seam and then the
 * walker-local payload blocks. The element's pair inputs are in src (dr wrapped,
 * r2, mass, secondary softening/type/zeta and the gated per-pair terms); its payload
 * values are in pl. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_evaluate_pair(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_member_inputs_t &in, gpu_grav_member_sums_t &sums, grav_pair_src_t &src, const gpu_grav_src_payload_t &pl)
{
    if(!((src.r2 > 0.0) && (src.mass > 0.0))) {return;}
    /* pair-wise gravity, PM truncation, and the accumulations inside the PM short-range
     * gate, via the shared evaluation (gravtree_force_kernel.h), the single home for the
     * pair physics on both walks */
    grav_pair_result_t res = grav_pair_evaluate_core(in.tgt, src, sums.out);
    const double r = res.r, fac_accel = res.fac_accel;
    (void) r; (void) fac_accel;
#ifdef GIZMO_GPU_EWALD_POT_CORRECTION
    /* Ewald periodic-image potential correction (mirrors forcetree.cc), from dr as the
     * evaluation left it. Pure-tree periodic only; under PMGRID the long-range potential
     * comes from the PM solver.  active is 1 in a healthy run (acquire failure hard-stops
     * the caller); the guard only covers the post-endrun drain. */
    if(ctx.ewald_pot.active) {
        grav_ewald_interp_weights ew = grav_ewald_interp_setup(src.dr[0], src.dr[1], src.dr[2], ctx.ewald_pot.fac_intp);
        sums.out.pot += src.mass * grav_ewald_interp_apply(ctx.ewald_pot.potcorr, ew);
    }
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    /* counted only for accepted interactions (r2>0, mass>0), mirroring forcetree.cc -- the
     * walk excludes the target's own (r2==0) leaf; gravtree.cc adds it back exactly once. */
    sums.tree_mass += src.mass;
#endif

    /* RT cluster payloads.  Structure mirrors forcetree.cc: OUTSIDE the PM short-range
     * gate, so for an out-of-range source fac_accel is the raw un-truncated value here
     * (used by the TREECOL column estimate, which has no PM-side completion);
     * RT_USE_GRAVTREE computes its own fac_rt from d_stellarlum independently. */
#ifdef RT_USE_TREECOL_FOR_NH
    {
        const double angular_bin_size = 4.0 * M_PI / RT_USE_TREECOL_FOR_NH;
        grav_treecol_accumulate(src.dr, r, fac_accel, pl.gasmass, src.mass, angular_bin_size, sums.treecol_angular_bins);
    }
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    /* per-interaction mass accumulation: each visited node contributes its multipole mass when within Rcrit */
    if(r < in.r_for_total_menclosed) {sums.m_enc_in_rcrit += src.mass;}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    grav_cr_lebron_accumulate(in.open.ptype, r, in.open.soft, pl.cr_injection, in.cr_active_gate, ctx.cr_data.t_max_cr, in.tgt.pm, sums.SubGrid_CosmicRayEnergyDensity);
#endif
#ifdef RT_USE_GRAVTREE
    {
        volatile int valid_gas_particle_for_rt = in.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt)
        {
            /* payload formulas in the shared helper; fac_rt computed there from d_stellarlum
             * (may differ from dr when RT_SEPARATELY_TRACK_LUMPOS; otherwise d_stellarlum == dr) */
            grav_rt_src_t rt_src = {}; rt_src.d_stellarlum = pl.d_stellarlum; rt_src.soft = in.open.soft; rt_src.mass_stellarlum = pl.mass_stellarlum;
#ifdef CHIMES_STELLAR_FLUXES
            rt_src.chimes_mass_stellarlum_G0 = pl.chimes_mass_stellarlum_G0; rt_src.chimes_mass_stellarlum_ion = pl.chimes_mass_stellarlum_ion;
#endif
#ifdef SINK_PHOTONMOMENTUM
            rt_src.mass_sinklumwt_forradfb = pl.mass_sinklumwt_forradfb;
#endif
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
            rt_src.fac_stellum = in.fac_stellum;
#endif
            grav_rt_accum_t rt_accum = {};
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            rt_accum.Rad_E_gamma = sums.Rad_E_gamma;
#endif
#ifdef CHIMES_STELLAR_FLUXES
            rt_accum.chimes_flux_G0 = sums.chimes_flux_G0; rt_accum.chimes_flux_ion = sums.chimes_flux_ion;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            rt_accum.incident_flux_uv = &sums.incident_flux_uv; rt_accum.incident_flux_euv = &sums.incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
            rt_accum.incident_flux_agn = &sums.incident_flux_agn;
#endif
#ifdef RT_OTVET
            rt_accum.RT_ET = sums.RT_ET;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
            rt_accum.Rad_Flux = sums.Rad_Flux;
#endif
            grav_rt_payload_accumulate(rt_src, rt_accum, sums.out.acc);
        }
    }
#endif /* RT_USE_GRAVTREE */
#ifdef DM_SCALARFIELD_SCREENING
    /* Yukawa-screened scalar-field force on non-gas targets (shared helper;
     * own table gate keyed on the dm-center distance, outside the main PM gate) */
    if(in.open.ptype != 0)
    {
        Vec3<double> d_dm = pl.d_dm;   /* the helper takes the displacement by non-const reference */
        grav_dm_scalarfield_accumulate(d_dm, pl.mass_dm_local, in.open.soft, in.tgt.pm, sums.out.acc);
    }
#endif
}

/* Load a particle leaf for a member through the P_dev adapter and evaluate it. The
 * source's position, velocity and mass are its state at the walk time
 * (gpu_grav_particle_source_at), except where the Hermite predictor replaces the position
 * and velocity; everything else is read as stored. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_evaluate_leaf(const gpu_grav_walk_ctx_t &ctx, int no, const gpu_grav_member_inputs_t &in, gpu_grav_member_sums_t &sums)
{
    struct particle_data *P_dev = ctx.P_dev;
    grav_pair_src_t src;
    gpu_grav_src_payload_t pl;
    src.h_p = -1.0; src.ptype_sec = -1; src.zeta_sec = 0.0;   /* unconditional, matching the CPU walk: consumed by the shared pair kernel */
#if defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    src.dv = Vec3<double>{0,0,0};
#endif
#ifdef SINK_DYNFRICTION_FROMTREE
    src.m_j_eff_for_df = 0.0;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = 0.0;
#endif

    const gpu_grav_particle_now_t motion = gpu_grav_particle_source_at(P_dev, ctx.CellP_dev, no, ctx.ti, ctx.tables, ctx.time_fault);
    Vec3<double> src_pos = motion.pos;
    Vec3<double> src_vel = motion.vel;   /* unconditional, mirroring forcetree.cc: the sink-proximity block reads it under SINK_CALC_DISTANCES, which several flags reach without the jerk or dynamical-friction terms */
#ifdef HERMITE_INTEGRATION
    /* On a Hermite pass a source the Hermite integrator owns but is not advancing this
       step is second-order wrong where it stands; evaluate it from its own start-of-step
       state instead. Same helper and same conditions as the host walk, so the two agree
       whichever one a step routes to. Single sources only; nothing is written back. */
    if(hermite_source_needs_prediction(no, P_dev, ctx.hermite.state)) {
        hermite_predict_source_state(no, P_dev, ctx.hermite.state, &ctx.tables, src_pos, src_vel);
    }
#endif
    src.dr = src_pos - in.open.pos;
    gravity_box_nearest_image(src.dr[0], src.dr[1], src.dr[2], -1);
    src.r2 = src.dr.norm_sq();
    src.mass = motion.mass;
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
    /* Lazy source payload for this local leaf: evaluate the SSOT helper on-device
     * (identical gates/formula to the eager prefill) instead of reading the dense
     * arrays. Only local particle leaves reach this branch -- foreign leaves and
     * nodes stay moment-backed. */
    struct gravtree_source_inputs_t lazy_src;
    if(ctx.use_lazy_source) { gravtree_fill_particle_source_inputs(no, P_dev, ctx.CellP_dev, &lazy_src); }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    /* per-interaction DM state for this leaf particle (mirrors forcetree.cc) */
    if(in.open.ptype != 0 && P_dev[no].Type == 1) { pl.d_dm = src.dr; pl.mass_dm_local = src.mass; }
    else { pl.d_dm = Vec3<double>{0,0,0}; pl.mass_dm_local = 0; }
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = grav_spherical_symmetry_r_from_center(src_pos[0],src_pos[1],src_pos[2],in.sph_center[0],in.sph_center[1],in.sph_center[2]);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    /* the secondary's previous-step tidal tensor (mirrors forcetree.cc) */
    {
        SymmetricTensor2<MyFloat> tmp = P_dev[no].tidal_tensorps_prevstep;
        for(int kk = 0; kk < 6; kk++) src.j_zeta_tidal_tensorps_prevstep.data[kk] = (double) tmp.data[kk];
    }
#endif
#if defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    src.dv = src_vel - in.vel;
#endif
#ifdef SINK_DYNFRICTION_FROMTREE
    src.m_j_eff_for_df = src.mass;
#endif
    /* secondary (leaf) softening, loaded unconditionally so a pair whose source softening
     * exceeds the target's gets the symmetrized max(h,h_p) force (mirrors forcetree.cc).
     * ptype_sec/zeta_sec stay gated -- only the adaptive symmetrize-by-averaging path uses them. */
    src.h_p = gpu_force_softening_kernel_radius(P_dev, no);
    src.ptype_sec = P_dev[no].Type;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
    if(src.ptype_sec == 0) {src.zeta_sec = gpu_get_ags_zeta(P_dev, no);}
#elif defined(ADAPTIVE_GRAVSOFT_FORALL)
    src.zeta_sec = gpu_get_ags_zeta(P_dev, no);
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    pl.gasmass = (P_dev[no].Type == 0) ? motion.mass : 0.0;
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
    /* gas at the inner edge of a sink's alpha-disk should not see a hole due to
     * the sink (mirrors forcetree.cc leaf branch + the node-moment kernel). */
    if(P_dev[no].Type == 5) {pl.gasmass = (double) P_dev[no].Sink_Mass_Reservoir;}
#endif
#endif
#ifdef RT_USE_GRAVTREE
    /* Load leaf luminosity only for valid gas targets (mirrors forcetree.cc; the RT
     * accumulation is gated the same way, so non-gas targets never read these and
     * skipping the loads avoids the wasted per-leaf table traffic). */
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {pl.mass_stellarlum[kf]=0.0;}}
    pl.d_stellarlum = {};
#ifdef CHIMES_STELLAR_FLUXES
    {int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {pl.chimes_mass_stellarlum_G0[kc]=0; pl.chimes_mass_stellarlum_ion[kc]=0;}}
#endif
#ifdef SINK_PHOTONMOMENTUM
    pl.mass_sinklumwt_forradfb = 0.0;
#endif
    {
        volatile int valid_gas_particle_for_rt = in.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt)
        {
            pl.d_stellarlum = src.dr;
            int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
                pl.mass_stellarlum[kf] = ctx.use_lazy_source ? (lazy_src.rt_active ? lazy_src.src_lum[kf] : (MyFloat)0)
                                                             : ctx.rt_data.src_lum[(long)no * N_RT_FREQ_BINS + kf];
#else
                pl.mass_stellarlum[kf] = ctx.rt_data.src_lum[(long)no * N_RT_FREQ_BINS + kf];
#endif
            }
#ifdef CHIMES_STELLAR_FLUXES
            for(kf=0; kf<CHIMES_LOCAL_UV_NBINS; kf++) {
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
                pl.chimes_mass_stellarlum_G0[kf]  = ctx.use_lazy_source ? (lazy_src.rt_active ? lazy_src.src_lum_G0[kf]  : 0.0) : ctx.rt_data.src_lum_G0[(long)no * CHIMES_LOCAL_UV_NBINS + kf];
                pl.chimes_mass_stellarlum_ion[kf] = ctx.use_lazy_source ? (lazy_src.rt_active ? lazy_src.src_lum_ion[kf] : 0.0) : ctx.rt_data.src_lum_ion[(long)no * CHIMES_LOCAL_UV_NBINS + kf];
#else
                pl.chimes_mass_stellarlum_G0[kf] = ctx.rt_data.src_lum_G0[(long)no * CHIMES_LOCAL_UV_NBINS + kf];
                pl.chimes_mass_stellarlum_ion[kf] = ctx.rt_data.src_lum_ion[(long)no * CHIMES_LOCAL_UV_NBINS + kf];
#endif
            }
#endif
#ifdef SINK_PHOTONMOMENTUM
            /* per-sink-leaf angle-weighted luminosity (shared formula helper) */
            if(P_dev[no].Type == 5) {
                double bhlum_t, bha0, bha1, bha2;
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
                if(ctx.use_lazy_source) {
                    bhlum_t = lazy_src.bh_active ? (double)lazy_src.bh_lum : 0.0;
                    bha0 = lazy_src.bh_active ? (double)lazy_src.bh_angle[0] : 0.0;
                    bha1 = lazy_src.bh_active ? (double)lazy_src.bh_angle[1] : 0.0;
                    bha2 = lazy_src.bh_active ? (double)lazy_src.bh_angle[2] : 0.0;
                } else
#endif
                {
                    bhlum_t = (double) ctx.sink_data.bh_lum[no];
                    bha0 = (double) ctx.sink_data.bh_angle[no][0]; bha1 = (double) ctx.sink_data.bh_angle[no][1]; bha2 = (double) ctx.sink_data.bh_angle[no][2];
                }
                pl.mass_sinklumwt_forradfb = grav_sink_fb_angleweight(bhlum_t, bha0, bha1, bha2, src.dr[0], src.dr[1], src.dr[2]);
            }
#endif
        }
    }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    /* leaf CR source injection (mirrors forcetree.cc) */
#if defined(GRAVTREE_SOURCE_DEVICE_TU)
    pl.cr_injection = ctx.use_lazy_source ? (double) lazy_src.cr_inject : (double) ctx.cr_data.cr_inject[no];
#else
    pl.cr_injection = (double) ctx.cr_data.cr_inject[no];
#endif
#endif
    /* Sink-distance + single-star timestepping tracking on particle leafs via the
     * shared helper (gravtree_force_kernel.h) -- CPU-walk semantics verbatim. */
#ifdef SINK_CALC_DISTANCES
    if((src.r2 > 0) && (src.mass > 0))
    {
        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = in.open.ptype; prox_target.pmass = in.pmass; prox_target.soft = in.open.soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
        prox_target.vel = in.vel;
#endif
        grav_sink_prox_leaf_src_t prox_src = {}; prox_src.src_type = P_dev[no].Type; prox_src.src_mass = motion.mass; prox_src.motion.vel = src_vel;   /* the state this interaction was evaluated at, mirroring forcetree.cc, so (dr, vel) stays a consistent pair on a Hermite pass */
#if defined(SPECIAL_POINT_MOTION) || defined(SPECIAL_POINT_WEIGHTED_MOTION)
        prox_src.motion.acc = P_dev[no].Acc_Total_PrevStep;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        prox_src.motion.max_feedback_vel = P_dev[no].MaxFeedbackVel;
#endif
        grav_sink_prox_leaf_accumulate(src.r2, src.dr, prox_target, prox_src, sums.sink_prox);
    }
#endif /* SINK_CALC_DISTANCES */

    gpu_grav_evaluate_pair(ctx, in, sums, src, pl);
}

/* Load an accepted node for a member through the SoA adapter and evaluate it, given
 * the node's geometry as this member derived it (gpu_grav_node_member_geometry). A
 * tagged real foreign single-particle leaf is consumed with particle-leaf secondary
 * semantics: the node payload supplies mass/h_p and the synthesized RT/sink/CR/tidal
 * terms (singleton-aggregate == particle value); the two leaf-identity fields the
 * moment cannot carry (Type + AGS_zeta) are restored via the shared seam so
 * grav_force_pair applies AGS symmetrization/zeta exactly as on the source's home rank. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_evaluate_node(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_node_prelude_t &nd, const gpu_grav_member_inputs_t &in, gpu_grav_member_sums_t &sums,
                       const Vec3<MyFloat> &s_node, MyFloat mass_node, const Vec3<double> &dr, double r2)
{
    const struct gpu_gravity_tree_soa_t *tree_soa = &ctx.tree_soa;
    const int idx = nd.idx;
    grav_pair_src_t src;
    gpu_grav_src_payload_t pl;
    src.dr = dr; src.r2 = r2; src.mass = mass_node;
    src.h_p = nd.msoft_node; src.ptype_sec = -1; src.zeta_sec = 0.0;
#if defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    src.dv = Vec3<double>{0,0,0};
#endif
#ifdef SINK_DYNFRICTION_FROMTREE
    src.m_j_eff_for_df = 0.0;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = 0.0;
#endif
    if(nd.fl_tag == 1) {
        grav_apply_foreign_leaf_identity(nd.fl_tag, nd.fl_type, nd.fl_zeta, nd.fl_soft, &src.ptype_sec, &src.zeta_sec, &src.h_p);
    }
#ifdef DM_SCALARFIELD_SCREENING
    /* per-interaction DM state for this accepted node (mirrors forcetree.cc): d_dm uses the
     * DM CoM s_dm, NOT the total CoM (s_node). */
    if(in.open.ptype != 0) {
        pl.d_dm[0] = (double)nd.motion.s_dm_[0] - in.open.pos[0];
        pl.d_dm[1] = (double)nd.motion.s_dm_[1] - in.open.pos[1];
        pl.d_dm[2] = (double)nd.motion.s_dm_[2] - in.open.pos[2];
        pl.mass_dm_local = (double)tree_soa->mass_dm[idx];
    } else { pl.d_dm = Vec3<double>{0,0,0}; pl.mass_dm_local = 0; }
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = grav_spherical_symmetry_r_from_center(s_node[0],s_node[1],s_node[2],in.sph_center[0],in.sph_center[1],in.sph_center[2]);
#else
    (void) s_node;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    /* the node's previous-step tidal tensor from the SoA (mirrors forcetree.cc) */
    for(int kk = 0; kk < 6; kk++) {
        src.j_zeta_tidal_tensorps_prevstep.data[kk] = (double) tree_soa->tidal_tensorps[(long)idx * 6 + kk];
    }
#endif
#if defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    src.dv[0] = (double) nd.motion.vs_[0] - in.vel[0];
    src.dv[1] = (double) nd.motion.vs_[1] - in.vel[1];
    src.dv[2] = (double) nd.motion.vs_[2] - in.vel[2];
#endif
#ifdef SINK_DYNFRICTION_FROMTREE
    {
        long np = tree_soa->N_part[idx];
        src.m_j_eff_for_df = (np > 0) ? (src.mass / (double)np) : 0.0;
    }
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    pl.gasmass = tree_soa->gasmass[idx];
#endif
#ifdef RT_USE_GRAVTREE
    /* Load node stellar luminosity only for valid gas targets (mirrors forcetree.cc; the RT
     * accumulation is gated the same way). */
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {pl.mass_stellarlum[kf]=0.0;}}
    pl.d_stellarlum = {};
#ifdef CHIMES_STELLAR_FLUXES
    {int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {pl.chimes_mass_stellarlum_G0[kc]=0; pl.chimes_mass_stellarlum_ion[kc]=0;}}
#endif
#ifdef SINK_PHOTONMOMENTUM
    pl.mass_sinklumwt_forradfb = 0.0;
#endif
    {
        volatile int valid_gas_particle_for_rt = in.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt)
        {
            int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {
                pl.mass_stellarlum[kf] = tree_soa->stellar_lum[idx * N_RT_FREQ_BINS + kf];
            }
#ifdef CHIMES_STELLAR_FLUXES
            for(kf=0; kf<CHIMES_LOCAL_UV_NBINS; kf++) {
                pl.chimes_mass_stellarlum_G0[kf] = tree_soa->chimes_stellar_lum_G0[(long)idx * CHIMES_LOCAL_UV_NBINS + kf];
                pl.chimes_mass_stellarlum_ion[kf] = tree_soa->chimes_stellar_lum_ion[(long)idx * CHIMES_LOCAL_UV_NBINS + kf];
            }
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            pl.d_stellarlum[0] = nd.motion.rt_s_[0] - in.open.pos[0];
            pl.d_stellarlum[1] = nd.motion.rt_s_[1] - in.open.pos[1];
            pl.d_stellarlum[2] = nd.motion.rt_s_[2] - in.open.pos[2];
            gravity_box_nearest_image(pl.d_stellarlum[0], pl.d_stellarlum[1], pl.d_stellarlum[2], -1);
#else
            pl.d_stellarlum = src.dr;
#endif
#ifdef SINK_PHOTONMOMENTUM
            /* node-aggregated sink angle-weighted luminosity (shared formula helper) */
            pl.mass_sinklumwt_forradfb = grav_sink_fb_angleweight((double) tree_soa->sink_lum[idx],
                                                                  (double) tree_soa->sink_lum_grad[idx][0], (double) tree_soa->sink_lum_grad[idx][1], (double) tree_soa->sink_lum_grad[idx][2],
                                                                  pl.d_stellarlum[0], pl.d_stellarlum[1], pl.d_stellarlum[2]);
#endif
        }
    }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    /* node-aggregated CR injection (mirrors forcetree.cc) */
    pl.cr_injection = (double) tree_soa->cr_injection[idx];
#endif
    /* Node-side sink distance + timestepping accumulators via the shared helper
     * (gravtree_force_kernel.h) -- CPU-walk semantics verbatim. The sink_vel/sink_acc
     * SoA fields are populated from Nodes[] by gpu_pseudo_update + gpu_moment_refresh. */
#ifdef SINK_CALC_DISTANCES
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    {
        Vec3<double> node_vs = Vec3<double>{(double)nd.motion.vs_[0], (double)nd.motion.vs_[1], (double)nd.motion.vs_[2]};
        grav_sink_prox_node_specialweighted(src.r2, node_vs, in.open.ptype, sums.sink_prox);
    }
#endif
    if(tree_soa->sink_mass[idx] > 0)
    {
        Vec3<double> sink_dr;
#ifdef SINK_NODE_MOTION_TRACKED
        for(int k = 0; k < 3; k++) {sink_dr[k] = nd.motion.sink_pos_[k] - in.open.pos[k];}   /* moves with its sinks: at the walk time */
#else
        for(int k = 0; k < 3; k++) {sink_dr[k] = tree_soa->sink_pos[idx][k] - in.open.pos[k];}   /* not moved between builds */
#endif
        gravity_box_nearest_image(sink_dr[0], sink_dr[1], sink_dr[2], -1);
        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = in.open.ptype; prox_target.pmass = in.pmass; prox_target.soft = in.open.soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
        prox_target.vel = in.vel;
#endif
        grav_sink_prox_node_src_t prox_src = {}; prox_src.sink_mass = (double) tree_soa->sink_mass[idx];
#if defined(SINGLE_STAR_FIND_BINARIES)
        prox_src.n_sink = (int) tree_soa->N_SINK[idx];
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SPECIAL_POINT_MOTION)
        prox_src.motion.vel = Vec3<double>{(double)nd.motion.sink_vel_[0], (double)nd.motion.sink_vel_[1], (double)nd.motion.sink_vel_[2]};
#endif
#if defined(SPECIAL_POINT_MOTION)
        prox_src.motion.acc = tree_soa->sink_acc[idx];
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        prox_src.motion.max_feedback_vel = tree_soa->MaxFeedbackVel[idx];
#endif
        grav_sink_prox_node_accumulate(src.r2, sink_dr, prox_src, prox_target, sums.sink_prox);
    }
#endif /* SINK_CALC_DISTANCES */

    gpu_grav_evaluate_pair(ctx, in, sums, src, pl);
}

/* Write a completed member's outputs to P_dev / CellP_dev (the host scatter loop in
 * gpu_gravtree_walk_primary copies them to P[] / CellP[]) and return the three the
 * caller collects directly. Mirrors forcetree.cc (mode=0). */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_member_finish(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_member_inputs_t &in, const gpu_grav_member_sums_t &sums, Vec3<double> &acc_out, int &ninter_out, double &pot_out)
{
    struct particle_data *P_dev = ctx.P_dev; const int target = in.target;
#ifdef RT_USE_GRAVTREE
    struct gas_cell_data *CellP_dev = ctx.CellP_dev;
    volatile int valid_gas_particle_for_rt = in.valid_gas_particle_for_rt;   /* nvc++ miscompiles raw boolean gates in device code */
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    {int k; for(k=0; k<RT_USE_TREECOL_FOR_NH; k++) {P_dev[target].ColumnDensityBins[k] = sums.treecol_angular_bins[k];}}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    P_dev[target].MencInRcrit = sums.m_enc_in_rcrit;
#endif
#ifdef RT_USE_GRAVTREE
#ifdef RT_OTVET
    if(valid_gas_particle_for_rt) {
        int k; for(k=0; k<N_RT_FREQ_BINS; k++) {CellP_dev[target].ET[k] = sums.RT_ET[k];}
    } else if(in.open.ptype == 0) {
        int k; for(k=0; k<N_RT_FREQ_BINS; k++) {CellP_dev[target].ET[k] = {};}
    }
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    if(valid_gas_particle_for_rt) {
        CellP_dev[target].Rad_Flux_UV  = sums.incident_flux_uv;
        CellP_dev[target].Rad_Flux_EUV = sums.incident_flux_euv;
    }
#endif
#ifdef CHIMES_STELLAR_FLUXES
    if(valid_gas_particle_for_rt) {
        int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {
            CellP_dev[target].Chimes_G0[kc]          = sums.chimes_flux_G0[kc];
            CellP_dev[target].Chimes_fluxPhotIon[kc] = sums.chimes_flux_ion[kc];
        }
    }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    if(valid_gas_particle_for_rt) {
        int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP_dev[target].Rad_E_gamma[kf] = sums.Rad_E_gamma[kf];}
    }
#endif
#ifdef SINK_COMPTON_HEATING
    if(valid_gas_particle_for_rt) {
        CellP_dev[target].Rad_Flux_AGN = sums.incident_flux_agn;
    }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    if(valid_gas_particle_for_rt) {
        int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP_dev[target].Rad_Flux[kf] = sums.Rad_Flux[kf];}
    }
#endif
#endif /* RT_USE_GRAVTREE */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    if(in.open.ptype == 0) {ctx.CellP_dev[target].SubGrid_CosmicRayEnergyDensity = sums.SubGrid_CosmicRayEnergyDensity;}
#endif
#ifdef SINK_CALC_DISTANCES
    P_dev[target].Min_Distance_to_Sink = sqrt(sums.sink_prox.Min_Distance_to_Sink2);
    P_dev[target].Min_xyz_to_Sink = sums.sink_prox.Min_xyz_to_Sink;
#ifdef SINGLE_STAR_FIND_BINARIES
    P_dev[target].is_in_a_binary = 0;
    P_dev[target].Min_Sink_OrbitalTime = sums.sink_prox.Min_Sink_OrbitalTime;
    if(sums.sink_prox.Min_Sink_OrbitalTime < MAX_REAL_NUMBER) {
        P_dev[target].is_in_a_binary = 1;
        P_dev[target].comp_Mass = sums.sink_prox.comp_Mass;
        P_dev[target].comp_dx = sums.sink_prox.comp_dx;
        P_dev[target].comp_dv = sums.sink_prox.comp_dv;
    }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    P_dev[target].Min_Sink_Approach_Time = sqrt(sums.sink_prox.Min_Sink_Approach_Time);
    P_dev[target].Min_Sink_Freefall_time = sqrt(sqrt(sums.sink_prox.Min_Sink_Freefall_time) / All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    P_dev[target].Min_Sink_FeedbackTime = sqrt(sums.sink_prox.Min_Sink_FeedbackTime);
#endif
#endif
#endif /* SINK_CALC_DISTANCES */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    P_dev[target].tidal_tensorps = sums.out.tidal_tensorps;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    P_dev[target].GravJerk = sums.out.jerk;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    P_dev[target].TreeMass = sums.tree_mass;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    P_dev[target].tidal_zeta = (MyFloat) sums.out.tidal_zeta;
#endif
#ifdef SPECIAL_POINT_MOTION
    P_dev[target].vel_of_nearest_special = Vec3<MyFloat>{(MyFloat)sums.sink_prox.vel_of_nearest_special[0],
                                                         (MyFloat)sums.sink_prox.vel_of_nearest_special[1],
                                                         (MyFloat)sums.sink_prox.vel_of_nearest_special[2]};
    P_dev[target].acc_of_nearest_special = Vec3<MyFloat>{(MyFloat)sums.sink_prox.acc_of_nearest_special[0],
                                                         (MyFloat)sums.sink_prox.acc_of_nearest_special[1],
                                                         (MyFloat)sums.sink_prox.acc_of_nearest_special[2]};
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    P_dev[target].weight_sum_for_special_point_smoothing = (MyFloat) sums.sink_prox.weight_sum_for_special_point_smoothing;
#endif
#endif
    acc_out = sums.out.acc;
    ninter_out = sums.out.ninter;
    pot_out = sums.out.pot;
}

/* -------------------------------------------------------------------------
 * gpu_gravtree_walk_one -- the walk for a single target: the units above composed
 * with each accepted element evaluated at encounter.
 *
 * Returns 1 on success (outputs written), 0 on failure (pseudo-particle hit; host runs the CPU walk
 * for this target).  Every source is read at the walk time without being drifted (the leaf and node
 * loads above).  Mirrors force_treeevaluate().
 * ---------------------------------------------------------------------- */
static KOKKOS_INLINE_FUNCTION int
gpu_gravtree_walk_one(const gpu_grav_walk_ctx_t &ctx, int target, Vec3<double> &acc_out, int &ninter_out, double &pot_out)
{
    gpu_grav_member_t mem;
    if(!gpu_grav_member_init(ctx, target, mem)) {acc_out = Vec3<double>{0,0,0}; ninter_out = 0; pot_out = 0.0; return 1;}
    const struct gpu_gravity_tree_soa_t *tree_soa = &ctx.tree_soa;
    const int treeBase = ctx.treeBase, treeParticleSlots = ctx.treeParticleSlots;
    const int pseudo_start = treeBase + ctx.maxNodes + ctx.maxForeignNodes;   /* foreign-node range below the pseudos */

    int no = treeBase;   /* root */
    while(no >= 0)
    {
        if(no >= treeParticleSlots && no < treeBase) {return 0;} /* gap: malformed tree -- defer; the CPU walk's guard stops loudly */
        if(no < treeParticleSlots) /* particle leaf */
        {
            if(gpu_grav_leaf_member_accepts(ctx, no, mem.in.open)) {gpu_grav_evaluate_leaf(ctx, no, mem.in, mem.sums);}
            no = tree_soa->nextnode_aux[no];
            continue;
        }
        if(no >= pseudo_start) {return 0;} /* pseudo-particle -- remote: host runs the CPU walk for this target */

        gpu_grav_node_prelude_t nd;
        const gpu_grav_node_step_t step = gpu_grav_node_prelude(ctx, no, nd);
        if(step == GPU_GRAV_NODE_SKIP_TO_SIBLING) {no = nd.sibling; continue;}
        if(step == GPU_GRAV_NODE_DESCEND) {no = nd.nextnode; continue;}
        Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
        if(!gpu_grav_node_member_geometry(ctx, nd, mem.in.open, s_node, mass_node, dr, r2)) {no = nd.sibling; continue;} /* pure-star node, star target */
        int note;
        const gravtree_open_t pred = gpu_grav_node_member_decide(ctx, nd, mem.in.open, mass_node, r2, note);
        if(note != GPU_GRAV_NOTE_NONE) {gpu_grav_note_commit(1, (note == GPU_GRAV_NOTE_UNSHIPPABLE) ? 1 : 0);}
        if(pred == GRAV_SKIP_NODE) {no = nd.sibling; continue;}
        if(pred == GRAV_OPEN_NODE) {no = nd.nextnode; continue;}
        gpu_grav_evaluate_node(ctx, nd, mem.in, mem.sums, s_node, mass_node, dr, r2);
        no = nd.sibling;
    }

    gpu_grav_member_finish(ctx, mem.in, mem.sums, acc_out, ninter_out, pot_out);
    return 1;
}


/* -------------------------------------------------------------------------
 * The packet engine: one team of threads walks the tree once for a packet of
 * up to q_dev adjacent targets and evaluates every member's forces.
 *
 * Member m is thread m (m < q_eff). It holds its own target state in registers
 * exactly as the single-target walk does and publishes only its opening inputs to
 * team scratch, where the walker reads them. The traversal is shared: the walker
 * carries one work item -- the next index, the first index NOT in the item, and
 * the members still taking part -- and at every node each member in the item
 * applies its own opening decision. Members that accept or skip a node the packet
 * then descends are absent from that subtree; they re-join through the
 * continuation (sibling, exit, mask) the walker stores when it descends, on its
 * bounded local stack or, when that is full, on the team frontier. Accepted
 * elements go to a bounded record chunk with the mask of the members that
 * accepted them; whenever it fills, and once more when the traversal ends, every
 * member evaluates its own records in append order through the same load and
 * evaluate units as the single-target walk. With one walker the traversal is that
 * walk's depth-first order, so each member's records are the sequence its own
 * walk would have accumulated, in the same order.
 *
 * A packet that meets a pseudo-particle, or a malformed index, or that has no room
 * left for a continuation, fails as a whole: nothing is written for any member. Its
 * members are then walked one at a time by the single-target device walk, which
 * leaves to the host loop exactly the targets that meet the pseudo-particle
 * themselves, as the flat route does. The packet's import notes are held in
 * scratch and counted only when it succeeds.
 *
 * Everything a team holds is in its level-0 scratch, sized for the launch shape
 * in one place (gpu_grav_packet_scratch_plan) so the legality check sees the whole
 * request; nothing here is sized by the configured packet size or by the tree.
 * ---------------------------------------------------------------------- */

#define GRAV_PACKET_MASK_BITS 64
#define GRAV_PACKET_Q_DEV_MAX 256               /* members per packet on the device; a larger configured packet size is walked as several packets */
#define GRAV_PACKET_LOCAL_STACK 16              /* continuations a walker keeps itself; the oldest moves to the frontier when full */
typedef unsigned long long grav_packet_mask_word_t;

struct gpu_grav_walk_item_t { int no, exit; };   /* a work item's indices; its mask words live beside it */

struct gpu_grav_packet_scratch_plan_t {
    int mask_words;
    int local_stack;   /* continuations per walker */
    size_t open_inputs, member_inputs, member_fold, frontier, frontier_masks, records, record_masks, local, local_masks, walker_masks, counters, bytes;
};

/* The scratch a team needs for one launch shape: q_dev members, team_size threads,
 * a frontier of frontier_cap items, a chunk of chunk_cap records and a per-walker
 * continuation stack local_stack deep.
 *
 * The frontier is the caller's: only several walkers share continuations, so a single-walker
 * team asks for none.  Zero leaves the region empty rather than merely unused -- the scratch
 * request is what team_size_max is asked about, so an unused region is a real cost in
 * occupancy, not just in bytes. */
static struct gpu_grav_packet_scratch_plan_t
gpu_grav_packet_scratch_plan(int q_dev, int team_size, int frontier_cap, int chunk_cap, int local_stack, int team_evaluates)
{
    struct gpu_grav_packet_scratch_plan_t p;
    p.mask_words = (q_dev + GRAV_PACKET_MASK_BITS - 1) / GRAV_PACKET_MASK_BITS;
    p.local_stack = local_stack;
    size_t off = 0;
    auto take = [&off](size_t bytes, size_t align) {off = ((off + align - 1) / align) * align; size_t here = off; off += bytes; return here;};
    p.open_inputs    = take((size_t) q_dev * sizeof(gpu_grav_open_inputs_t), alignof(gpu_grav_open_inputs_t));
    /* a flavour whose whole team evaluates also publishes each member's full inputs, which every lane
       working on that member reads, and one partial sum per thread, through which the lanes of a member
       are folded when the packet commits */
    p.member_inputs  = take(team_evaluates ? (size_t) q_dev * sizeof(gpu_grav_member_inputs_t) : 0, alignof(gpu_grav_member_inputs_t));
    p.member_fold    = take(team_evaluates ? (size_t) team_size * sizeof(gpu_grav_member_sums_t) : 0, alignof(gpu_grav_member_sums_t));
    p.frontier       = take((size_t) frontier_cap * sizeof(gpu_grav_walk_item_t), alignof(gpu_grav_walk_item_t));
    p.frontier_masks = take((size_t) frontier_cap * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.records        = take((size_t) chunk_cap * sizeof(grav_walk_record_t), alignof(grav_walk_record_t));
    p.record_masks   = take((size_t) chunk_cap * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.local          = take((size_t) team_size * local_stack * sizeof(gpu_grav_walk_item_t), alignof(gpu_grav_walk_item_t));
    p.local_masks    = take((size_t) team_size * local_stack * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.walker_masks   = take((size_t) team_size * 3 * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));   /* a walker's item mask + its two scratch masks */
    /* 16 ints, deliberately not sized to GRAV_PACKET_CTR_COUNT: the scratch request feeds
     * team_size_max, so shrinking it could change the launch shape and with it what is being
     * measured.  The static_assert that keeps the slack honest sits with the enum below, which
     * is where a counter would be added. */
    p.counters       = take(16 * sizeof(int), alignof(long long));
    p.bytes = off;
    return p;
}

/* Why a packet gave itself up.  These are NOT interchangeable, which is the whole reason they
 * are counted apart: meeting a pseudo-particle is the design working (the host owns those
 * targets), while running out of room for a continuation is the engine hitting its own budget on
 * a divergent traversal -- correctness is safe either way, but the second is a PERFORMANCE CLIFF
 * in exactly the deep clustered subtrees this engine exists to speed up, and a single total would
 * bury it under the ordinary pseudo replays. A malformed index is a broken tree. */
enum {
    GRAV_PACKET_FAIL_NONE = 0,
    GRAV_PACKET_FAIL_MALFORMED_INDEX,    /* an index in the gap between particle slots and the node base */
    GRAV_PACKET_FAIL_PSEUDO,             /* a pseudo-particle: the host walks every member */
    GRAV_PACKET_FAIL_NO_CONTINUATION,    /* local stack and frontier both full */
    GRAV_PACKET_FAIL_RECORD_UNUSABLE,    /* a record the flush could not reproduce */
    GRAV_PACKET_FAIL_NO_PROGRESS,        /* a round in which no walker advanced and none was waiting on the chunk */
    GRAV_PACKET_FAIL_REASONS
};

/* team-scope counters, one int slot each */
enum {
    GRAV_PACKET_CTR_RECORDS = 0,       /* records in the chunk; reserved by atomic increment, clamped to the chunk once per round */
    /* The frontier's three counters.  Slots [0, FR_PUB) are the published stack and do not move
     * for the length of a round; new pushes grow above it and are read by nobody until the round
     * boundary compacts them in.  That is what makes a failed reservation harmless: an overflowing
     * pusher leaves its index above the cap and simply does not write, where a pusher that tried to
     * hand the index back could race a second overflowing pusher and wind the counter back below a
     * slot already holding a live continuation, which the next push would then overwrite. */
    GRAV_PACKET_CTR_FR_PUB,            /* published items, poppable this round */
    GRAV_PACKET_CTR_FR_TAKEN,          /* published items claimed this round */
    GRAV_PACKET_CTR_FR_NEW,            /* reservations made above the published block this round */
    GRAV_PACKET_CTR_FAILED,            /* the packet failed: pseudo-particle, malformed index, or no room for a continuation */
    GRAV_PACKET_CTR_FAIL_REASON,       /* which of the above; read only when FAILED is set */
    GRAV_PACKET_CTR_DONE,              /* the traversal is finished: a team-level verdict taken at the round boundary, never by one walker */
    GRAV_PACKET_CTR_WORKING,           /* walkers still holding an item or a continuation at the end of the round */
    GRAV_PACKET_CTR_STEPPED,           /* nodes stepped this round, by the whole team */
    GRAV_PACKET_CTR_NOTE_INCOMPLETE,   /* import-completeness notes, held until success */
    GRAV_PACKET_CTR_NOTE_UNSHIPPABLE,
    GRAV_PACKET_CTR_COUNT
};

/* The team's counter block is a fixed 16 ints in the scratch plan above, and it is the LAST
 * region taken -- so a counter added past that bound would not overlap another region, it would
 * run off the end of the team's scratch with nothing to say so. */
static_assert(GRAV_PACKET_CTR_COUNT <= 16,
              "the per-team counter block is 16 ints; widen it in gpu_grav_packet_scratch_plan");

/* What a walker does when it reaches a node or a leaf: the per-node decision, and what
 * becomes of the elements that decision accepts.
 *
 * The traversal itself -- the (no, exit, mask) work item, the local LIFO and the frontier, the
 * give-up accounting -- is the same whatever is being decided, so it stays in the engine below
 * and the flavour is a template parameter of it.  One engine with a flavour parameter rather
 * than one walk per flavour is what stops a second device traversal growing up beside this one
 * with its own copy of the periodic-wrap convention and its own rule for which index classes a
 * walk may follow (mesh/device_tree_walk.h states the same invariant for the neighbour walk).
 *
 * A flavour returns YIELD to stop the walker for this round with its item intact; the engine
 * resumes at the same node afterwards, which is why a decision must be reproducible rather than
 * consumed.  CONTINUE means the walker advances as the traversal says.
 */
enum gpu_grav_packet_step_t {
    GRAV_PACKET_STEP_CONTINUE = 0,
    GRAV_PACKET_STEP_YIELD          /* chunk full, or the packet gave up; the engine returns */
};

/* The masked flavour: every member judges every node for
 * itself, and the elements it accepts are reserved in the record chunk for the member-major
 * flush.  Bitwise per member against the single-target walk at n_walkers = 1. */
struct GravPacketMaskedPolicy {
    /* One thread per member evaluates that member's records; see GravPacketMaskedTeamPolicy for the
       flavour that spreads them over a team wider than the packet. */
    static constexpr bool team_evaluates   = false;
    /* Compile-time so the walker's ring index stays a mask rather than a division. */
    static constexpr int  local_stack      = GRAV_PACKET_LOCAL_STACK;
    /* Members share one traversal and then evaluate its records in parallel, so the packet is
       what makes that sharing worth having: the configured size. */
    static constexpr int  packet_size      = TREE_QUERY_PACKET_SIZE;

    /* Which members take this particle leaf, and the record that the flush evaluates. */
    template <class Engine>
    KOKKOS_INLINE_FUNCTION gpu_grav_packet_step_t
    visit_leaf(const Engine &e, int no, const gpu_grav_open_inputs_t *open,
               const grav_packet_mask_word_t *mask, grav_packet_mask_word_t *accept_mask,
               int *ctr, grav_walk_record_t *records, grav_packet_mask_word_t *rmasks) const
    {
        const int W = e.plan.mask_words;
        Engine::mask_clear(accept_mask, W);
        for(int m = 0; m < e.q_dev; m++) {
            if(!Engine::mask_test(mask, m)) {continue;}
            if(gpu_grav_leaf_member_accepts(e.ctx, no, open[m])) {Engine::mask_set(accept_mask, m);}
        }
        if(!Engine::mask_any(accept_mask, W)) {return GRAV_PACKET_STEP_CONTINUE;}
        const int r = e.record_reserve(ctr);
        if(r < 0) {return GRAV_PACKET_STEP_YIELD;}   /* chunk full: the flush follows; resume here */
        records[r].no = no; records[r].kind = GRAV_NODE_LOCAL; records[r].leaf_tag = LET_LEAF_TAG_NODE;
        Engine::mask_copy(rmasks + (size_t) r * W, accept_mask, W);
        return GRAV_PACKET_STEP_CONTINUE;
    }

    /* Which members accept this node's multipole and which must descend into it, and the record
     * for the accepting ones.  `open_mask` is the engine's: it decides the continuation from it. */
    template <class Engine>
    KOKKOS_INLINE_FUNCTION gpu_grav_packet_step_t
    visit_node(const Engine &e, int no, const gpu_grav_node_prelude_t &nd, const gpu_grav_open_inputs_t *open,
               const grav_packet_mask_word_t *mask, grav_packet_mask_word_t *accept_mask,
               grav_packet_mask_word_t *open_mask, int *ctr, grav_walk_record_t *records,
               grav_packet_mask_word_t *rmasks, int &n_note, int &n_unship) const
    {
        const int W = e.plan.mask_words;
        Engine::mask_clear(accept_mask, W); Engine::mask_clear(open_mask, W);
        n_note = 0; n_unship = 0;
        for(int m = 0; m < e.q_dev; m++) {
            if(!Engine::mask_test(mask, m)) {continue;}
            Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
            if(!gpu_grav_node_member_geometry(e.ctx, nd, open[m], s_node, mass_node, dr, r2)) {continue;}   /* a star member is done with a pure-star node */
            int note;
            const gravtree_open_t pred = gpu_grav_node_member_decide(e.ctx, nd, open[m], mass_node, r2, note);
            if(note != GPU_GRAV_NOTE_NONE) {n_note++; if(note == GPU_GRAV_NOTE_UNSHIPPABLE) {n_unship++;}}
            if(pred == GRAV_SKIP_NODE) {continue;}
            if(pred == GRAV_OPEN_NODE) {Engine::mask_set(open_mask, m); continue;}
            Engine::mask_set(accept_mask, m);
        }
        if(Engine::mask_any(accept_mask, W)) {
            const int r = e.record_reserve(ctr);
            if(r < 0) {return GRAV_PACKET_STEP_YIELD;}   /* chunk full: resume at this node after the flush (its decisions are re-made identically) */
            records[r].no = no; records[r].kind = nd.node_kind; records[r].leaf_tag = nd.fl_tag;
            Engine::mask_copy(rmasks + (size_t) r * W, accept_mask, W);
        }
        return GRAV_PACKET_STEP_CONTINUE;
    }
};

/* The same decisions and records, for a team WIDER than its packet (the cooperative schedule):
 * once the chunk is full every thread evaluates, the team divided into one group of lanes per
 * member, each lane holding a partial sum that is folded into the member's first lane when the
 * packet commits.
 *
 * It is a separate flavour, not a runtime switch, so that each schedule is compiled for the work it
 * actually does: a team of one thread per member needs neither the published inputs nor the fold,
 * and their code changes how the compiler lays out the member state for the WHOLE kernel (one kernel
 * serving both schedules ran a quarter slower per call at large N under the FIRE physics). */
struct GravPacketMaskedTeamPolicy : GravPacketMaskedPolicy {
    static constexpr bool team_evaluates = true;
};

template <class Policy>
struct GpuGravPacketWalk {
    using TeamMember = Kokkos::TeamPolicy<>::member_type;
    gpu_grav_walk_ctx_t ctx;
    const int *d_idx;   /* candidates in ActiveParticleList order */
    int n_cand, q_dev, frontier_cap, chunk_cap;
    /* How many of the team's threads traverse, and how far each one goes before returning to the
     * round boundary.
     *
     * n_walkers = 1 is the single-walker traversal: the frontier publishes on the push and pops
     * from its top, which is depth-first order and bitwise against the single-target walk.
     * n_walkers = T is the cooperative traversal -- every thread walks, and an opened subtree is
     * offered to the others through the frontier.  Members are still threads 0..q_eff-1, so a wide
     * team with a narrow packet is a team where most threads only traverse, which is the shape a
     * step with few targets and a deep tree wants.
     *
     * steps_per_round bounds a walker's run so that what it published becomes visible to an idle
     * walker: too large and the others wait, too small and the round boundary is the cost.  It is
     * a column of the launch table, not a knob. */
    int n_walkers, steps_per_round;
    struct gpu_grav_packet_scratch_plan_t plan;
    Vec3<double> *d_acc; int *d_ninter; double *d_pot; int *d_failed;
    int *d_fail_by_reason;   /* [GRAV_PACKET_FAIL_REASONS] packets given up, by reason */
    /* The flavour, by value; it carries no state, only the per-node decisions. */
    Policy policy;

    /* Give the packet up, and say why.  The first reason recorded is kept: it is the one that
     * actually stopped the traversal, and a later lane writing over it would report a symptom. */
    /* Several threads can reach a give-up at once, so the flag is CLAIMED rather than assigned:
     * whoever claims it writes the reason, and the reason the report carries is the one that
     * actually stopped a traversal first.
     *
     * Always claimed, never a cheaper unsynchronised path for a single walker: this is reached
     * from the record-evaluation loop as well as the traversal, and that loop runs on every
     * member thread whatever n_walkers is. Keying the cheap path on n_walkers would be sound only
     * while the two happen to move together, and the failure it would produce -- FAILED set with
     * the reason still zero, tallied into the unused slot -- is invisible in the output and
     * corrupts the very counter the continuation-budget question is read from. It costs an atomic
     * on a path that by construction runs at most once per packet. */
    KOKKOS_INLINE_FUNCTION void fail(int *ctr, int reason) const
    {
        if(Kokkos::atomic_compare_exchange(&ctr[GRAV_PACKET_CTR_FAILED], 0, 1) == 0) {
            ctr[GRAV_PACKET_CTR_FAIL_REASON] = reason;
        }
    }

    /* Add to a team counter several walkers may be adding to at once.  Plain arithmetic with one
     * walker, so the single-walker traversal keeps the instruction sequence it was gated with. */
    KOKKOS_INLINE_FUNCTION void ctr_add(int *slot, int by) const
    {
        if(n_walkers == 1) {*slot += by;} else {Kokkos::atomic_fetch_add(slot, by);}
    }

    KOKKOS_INLINE_FUNCTION static void mask_clear(grav_packet_mask_word_t *m, int words) {for(int w = 0; w < words; w++) {m[w] = 0ULL;}}
    KOKKOS_INLINE_FUNCTION static void mask_copy(grav_packet_mask_word_t *dst, const grav_packet_mask_word_t *src, int words) {for(int w = 0; w < words; w++) {dst[w] = src[w];}}
    KOKKOS_INLINE_FUNCTION static int  mask_test(const grav_packet_mask_word_t *m, int b) {return (int)((m[b / GRAV_PACKET_MASK_BITS] >> (b % GRAV_PACKET_MASK_BITS)) & 1ULL);}
    KOKKOS_INLINE_FUNCTION static void mask_set(grav_packet_mask_word_t *m, int b) {m[b / GRAV_PACKET_MASK_BITS] |= (1ULL << (b % GRAV_PACKET_MASK_BITS));}
    KOKKOS_INLINE_FUNCTION static int  mask_any(const grav_packet_mask_word_t *m, int words) {for(int w = 0; w < words; w++) {if(m[w]) {return 1;}} return 0;}
    KOKKOS_INLINE_FUNCTION static int  mask_equal(const grav_packet_mask_word_t *a, const grav_packet_mask_word_t *b, int words) {for(int w = 0; w < words; w++) {if(a[w] != b[w]) {return 0;}} return 1;}

    /* Offer a continuation to the other walkers, or keep it if there is no room.  Returns 1 when
     * the item is on the frontier and the caller is rid of it.
     *
     * With one walker the frontier is a plain stack that publishes on the push, which is the
     * depth-first overflow store the single-walker traversal is gated on.  With several, a push
     * reserves above the published block and a reservation past the cap is simply not taken: the
     * counter is left where it is and clamped once, at the round boundary, by the one thread doing
     * the compaction.  Handing a failed index back would be the race -- two overflowing pushers
     * each decrementing can leave the counter below a slot that already holds a live continuation. */
    KOKKOS_INLINE_FUNCTION int frontier_push(int *ctr, gpu_grav_walk_item_t *frontier, grav_packet_mask_word_t *fmasks,
                                             const gpu_grav_walk_item_t &item, const grav_packet_mask_word_t *item_mask) const
    {
        const int W = plan.mask_words;
        int slot;
        if(n_walkers == 1) {
            if(ctr[GRAV_PACKET_CTR_FR_PUB] == frontier_cap) {return 0;}
            slot = ctr[GRAV_PACKET_CTR_FR_PUB]++;
        } else {
            const int n = Kokkos::atomic_fetch_add(&ctr[GRAV_PACKET_CTR_FR_NEW], 1);
            slot = ctr[GRAV_PACKET_CTR_FR_PUB] + n;
            if(slot >= frontier_cap) {return 0;}
        }
        frontier[slot] = item;
        mask_copy(fmasks + (size_t) slot * W, item_mask, W);
        return 1;
    }

    /* Take a published item, if this walker gets one.  A claim past the published block is not
     * handed back either: the compaction clamps it against what was actually published. */
    KOKKOS_INLINE_FUNCTION int frontier_pop(int *ctr, const gpu_grav_walk_item_t *frontier, const grav_packet_mask_word_t *fmasks,
                                            int &no, int &exit, grav_packet_mask_word_t *mask) const
    {
        const int W = plan.mask_words;
        const int pub = ctr[GRAV_PACKET_CTR_FR_PUB];
        if(pub <= 0) {return 0;}
        int slot;
        if(n_walkers == 1) {
            slot = --ctr[GRAV_PACKET_CTR_FR_PUB];
        } else {
            const int k = Kokkos::atomic_fetch_add(&ctr[GRAV_PACKET_CTR_FR_TAKEN], 1);
            if(k >= pub) {return 0;}
            slot = pub - 1 - k;   /* newest first, over a block that does not move during the round */
        }
        no = frontier[slot].no; exit = frontier[slot].exit;
        mask_copy(mask, fmasks + (size_t) slot * W, W);
        return 1;
    }

    /* Reserve a record slot.  Past the chunk nothing is written and nothing is handed back: the
     * walker yields with its item intact, the flush follows, and the node is reached again and
     * decided identically.  The clamp at the round boundary is what keeps the count honest. */
    KOKKOS_INLINE_FUNCTION int record_reserve(int *ctr) const
    {
        const int r = (n_walkers == 1) ? ctr[GRAV_PACKET_CTR_RECORDS]++
                                       : Kokkos::atomic_fetch_add(&ctr[GRAV_PACKET_CTR_RECORDS], 1);
        return (r < chunk_cap) ? r : -1;
    }

    /* One member's evaluation of a recorded element: the leaf through the P_dev adapter,
     * a node re-derived from its index through the same prelude and geometry the walker
     * used (the same statements, so the same values), then the shared evaluation. */
    /* Returns 1 when the element was evaluated, 0 when the record cannot be reproduced.
     *
     * A record's accept mask is built only from members whose prelude reached the decision and
     * whose geometry resolved, so re-deriving here reaches the same two answers -- the node index
     * and the member's inputs are the same, and nothing between the decision and the flush moves
     * them. Both are still tested rather than assumed: the invariant lives seventy lines away in
     * the walker, and the cost of it being wrong is not a skipped contribution but an evaluation
     * on uninitialised geometry. Declining hands the whole packet to the replay chain that
     * already exists for every other reason a packet cannot be completed on the device. */
    KOKKOS_INLINE_FUNCTION int evaluate_record(int no, const gpu_grav_member_inputs_t &in, gpu_grav_member_sums_t &sums) const
    {
        if(no < ctx.treeParticleSlots) {gpu_grav_evaluate_leaf(ctx, no, in, sums); return 1;}
        gpu_grav_node_prelude_t nd;
        const gpu_grav_node_step_t step = gpu_grav_node_prelude(ctx, no, nd);
        if(step != GPU_GRAV_NODE_DECIDE) {return 0;}   /* an accepted node is one the prelude handed to the decision */
        Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
        if(!gpu_grav_node_member_geometry(ctx, nd, in.open, s_node, mass_node, dr, r2)) {return 0;}   /* a member with the bit set is never a star seeing a pure-star node */
        gpu_grav_evaluate_node(ctx, nd, in, sums, s_node, mass_node, dr, r2);
        return 1;
    }

    /* The walker: advance its item until its step budget runs out, the chunk is full, it runs out
     * of work, or the packet fails. All of its state persists across the round boundary and the
     * flush that follows: the item in the caller's variables, the continuations in scratch.
     *
     * Running out of work is not the end of the traversal -- another walker may be about to
     * publish -- so the walker simply comes back idle and the round boundary decides. */
    KOKKOS_INLINE_FUNCTION void walk(int &no, int &exit, grav_packet_mask_word_t *mask, int &item_live,
                                     grav_packet_mask_word_t *accept_mask, grav_packet_mask_word_t *open_mask,
                                     int &local_head, int &local_count,
                                     gpu_grav_walk_item_t *local, grav_packet_mask_word_t *lmasks,
                                     const gpu_grav_open_inputs_t *open,
                                     gpu_grav_walk_item_t *frontier, grav_packet_mask_word_t *fmasks,
                                     grav_walk_record_t *records, grav_packet_mask_word_t *rmasks, int *ctr) const
    {
        const int W = plan.mask_words;
        const int treeBase = ctx.treeBase, treeParticleSlots = ctx.treeParticleSlots;
        const int pseudo_start = treeBase + ctx.maxNodes + ctx.maxForeignNodes;
        int budget = steps_per_round, stepped = 0;

        while(1)
        {
            /* the item is complete: resume the newest continuation, else one the frontier has
               published, else come back idle and let the round boundary decide */
            if(!item_live || no == exit) {
                item_live = 0;
                if(local_count > 0) {
                    const int slot = (local_head + local_count - 1) % Policy::local_stack;
                    no = local[slot].no; exit = local[slot].exit; mask_copy(mask, lmasks + (size_t) slot * W, W);
                    local_count--; item_live = 1; continue;
                }
                if(frontier_pop(ctr, frontier, fmasks, no, exit, mask)) {item_live = 1; continue;}
                break;
            }
            if(budget <= 0) {break;}   /* back to the round boundary, item intact */
            budget--; stepped++;
            if(no < 0) {item_live = 0; continue;}   /* the single-target walk's end-of-walk sentinel: nothing beyond it */
            if(no >= treeParticleSlots && no < treeBase) {fail(ctr, GRAV_PACKET_FAIL_MALFORMED_INDEX); break;}   /* gap: malformed tree; the host walk stops loudly */
            if(no < treeParticleSlots) /* particle leaf: per member, the star-star pass */
            {
                if(policy.visit_leaf(*this, no, open, mask, accept_mask, ctr, records, rmasks) == GRAV_PACKET_STEP_YIELD) {break;}
                no = ctx.tree_soa.nextnode_aux[no];
                continue;
            }
            if(no >= pseudo_start) {fail(ctr, GRAV_PACKET_FAIL_PSEUDO); break;}   /* pseudo-particle: the host walks every member */

            gpu_grav_node_prelude_t nd;
            const gpu_grav_node_step_t step = gpu_grav_node_prelude(ctx, no, nd);
            if(step == GPU_GRAV_NODE_SKIP_TO_SIBLING) {no = nd.sibling; continue;}
            if(step == GPU_GRAV_NODE_DESCEND) {no = nd.nextnode; continue;}

            /* the flavour judges the node and records what it accepts; the engine owns only
               what the resulting open mask means for the traversal */
            int n_note = 0, n_unship = 0;
            if(policy.visit_node(*this, no, nd, open, mask, accept_mask, open_mask, ctr, records, rmasks,
                                 n_note, n_unship) == GRAV_PACKET_STEP_YIELD) {break;}
            if(n_note)   {ctr_add(&ctr[GRAV_PACKET_CTR_NOTE_INCOMPLETE], n_note);}
            if(n_unship) {ctr_add(&ctr[GRAV_PACKET_CTR_NOTE_UNSHIPPABLE], n_unship);}
            if(!mask_any(open_mask, W)) {no = nd.sibling; continue;}
            /* The packet descends for the openers; the others re-join at the node's sibling
               through the continuation (sibling, exit, mask).
             *
             * WHERE that continuation goes is the difference between one walker and many, and
             * it is the whole of the cooperative traversal.
             *
             * One walker keeps it, because the only thing the frontier can do for it is hold
             * an overflow, and keeping the newest locally and spilling the OLDEST is what
             * preserves depth-first order.  It also descends in place when nobody was excluded,
             * since there is nothing to hand anyone.
             *
             * Several walkers OFFER it, on every open and even when nobody was excluded: a
             * subtree nobody can reach is a subtree the other walkers sit idle through, and
             * the sibling continuation is precisely the independent piece of work to give
             * away.  Only when the frontier is full does the walker keep it, and only when its
             * own stack is full too does the packet give up -- counted, and handed to the
             * replay, never a spin.
             *
             * An empty continuation (the node's sibling IS the item's exit) is not worth a
             * slot in either regime. */
            const int share = (n_walkers > 1);
            const int have_continuation = (nd.sibling != exit);
            const int all_descend = mask_equal(open_mask, mask, W);
            if(!share && all_descend) {no = nd.nextnode; continue;}   /* same item, deeper */
            if(have_continuation) {
                gpu_grav_walk_item_t cont; cont.no = nd.sibling; cont.exit = exit;
                int placed = 0;
                if(share) {placed = frontier_push(ctr, frontier, fmasks, cont, mask);}
                if(!placed) {
                    if(local_count == Policy::local_stack) {
                        /* the oldest of this walker's own continuations makes room; with one
                           walker that is the spill that keeps the local stack depth-first */
                        gpu_grav_walk_item_t oldest = local[local_head];
                        if(!frontier_push(ctr, frontier, fmasks, oldest, lmasks + (size_t) local_head * W)) {
                            fail(ctr, GRAV_PACKET_FAIL_NO_CONTINUATION); break;
                        }
                        local_head = (local_head + 1) % Policy::local_stack; local_count--;
                    }
                    const int slot = (local_head + local_count) % Policy::local_stack;
                    local[slot] = cont; mask_copy(lmasks + (size_t) slot * W, mask, W);
                    local_count++;
                }
            }
            no = nd.nextnode; exit = nd.sibling; mask_copy(mask, open_mask, W);
        }
        if(stepped) {ctr_add(&ctr[GRAV_PACKET_CTR_STEPPED], stepped);}
    }

    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember &team) const
    {
        const int W = plan.mask_words;
        const int p = team.league_rank();
        const int first = p * q_dev;
        const int q_eff = (n_cand - first < q_dev) ? (n_cand - first) : q_dev;
        const int t = team.team_rank();
        char *scratch = (char *) team.team_scratch(0).get_shmem(plan.bytes);
        gpu_grav_open_inputs_t   *open      = (gpu_grav_open_inputs_t *)   (scratch + plan.open_inputs);
        gpu_grav_walk_item_t     *frontier  = (gpu_grav_walk_item_t *)     (scratch + plan.frontier);
        grav_packet_mask_word_t  *fmasks    = (grav_packet_mask_word_t *)  (scratch + plan.frontier_masks);
        grav_walk_record_t       *records   = (grav_walk_record_t *)       (scratch + plan.records);
        grav_packet_mask_word_t  *rmasks    = (grav_packet_mask_word_t *)  (scratch + plan.record_masks);
        gpu_grav_walk_item_t     *local     = (gpu_grav_walk_item_t *)     (scratch + plan.local) + (size_t) t * Policy::local_stack;
        grav_packet_mask_word_t  *lmasks    = (grav_packet_mask_word_t *)  (scratch + plan.local_masks) + (size_t) t * Policy::local_stack * W;
        grav_packet_mask_word_t  *wmasks    = (grav_packet_mask_word_t *)  (scratch + plan.walker_masks) + (size_t) t * 3 * W;
        int                      *ctr       = (int *)                      (scratch + plan.counters);

        /* the member this thread owns, if any; every thread publishes an entry so the
           walker's loop over q_dev members reads only initialised inputs */
        gpu_grav_member_inputs_t in;
        gpu_grav_member_sums_t sums;
        const int have_member = (t < q_eff);
        /* The lanes that EVALUATE a member are not the lanes that traverse.  In the flavour whose team
         * evaluates, every thread walks and, once the chunk is full, every thread also evaluates, the
         * team being divided into one group of `lanes` threads per member, so a single target's
         * records are spread over the team instead of queuing on one thread while the rest wait at
         * the barrier.  Otherwise lanes == 1: one thread per member. */
        const int lanes    = (Policy::team_evaluates && q_eff > 0) ? (team.team_size() / q_eff) : 1;
        const int serves   = (q_eff > 0 && t < lanes * q_eff) ? 1 : 0;
        const int member   = serves ? (t / lanes) : 0;   /* which member this thread evaluates for */
        const int sub_lane = serves ? (t % lanes) : 0;
        if(have_member) {
            if constexpr (Policy::team_evaluates) {   /* published, for every lane working on this member */
                gpu_grav_member_inputs_t *member_inputs = (gpu_grav_member_inputs_t *) (scratch + plan.member_inputs);
                (void) gpu_grav_member_inputs_init(ctx, d_idx[first + t], member_inputs[t]);
                open[t] = member_inputs[t].open;
            } else {               /* this thread is the member's only lane */
                (void) gpu_grav_member_inputs_init(ctx, d_idx[first + t], in);
                open[t] = in.open;
            }
        } else if(t < q_dev) {
            open[t].alive = 0;
        }
        if(t == 0) {for(int c = 0; c < GRAV_PACKET_CTR_COUNT; c++) {ctr[c] = 0;}}
        team.team_barrier();
        /* every thread holds a partial sum, at the identity until it evaluates something, and a
           thread that evaluates holds its member's inputs */
        gpu_grav_member_sums_init(sums);
        if constexpr (Policy::team_evaluates) {
            if(serves) {in = ((const gpu_grav_member_inputs_t *) (scratch + plan.member_inputs))[member];}
        }

        /* Every thread below n_walkers traverses; the root item starts with one of them and the
         * rest join as subtrees are offered. Members are threads below q_eff, so a team wider than
         * the packet is a team whose extra threads only ever walk -- which is the point when a
         * rank has few targets and a deep tree. */
        const int is_walker = (t < n_walkers);
        int no = ctx.treeBase, exit = -1;
        grav_packet_mask_word_t *mask = wmasks, *accept_mask = wmasks + W, *open_mask = wmasks + 2 * W; mask_clear(mask, W);
        for(int m = 0; m < q_dev; m++) {if(open[m].alive) {mask_set(mask, m);}}
        int item_live = (t == 0) ? mask_any(mask, W) : 0;   /* a packet of massless targets walks nothing */
        int local_head = 0, local_count = 0;

        while(1)
        {
            if(is_walker && !ctr[GRAV_PACKET_CTR_FAILED]) {
                walk(no, exit, mask, item_live, accept_mask, open_mask, local_head, local_count, local, lmasks, open, frontier, fmasks, records, rmasks, ctr);
                if(item_live || local_count > 0) {ctr_add(&ctr[GRAV_PACKET_CTR_WORKING], 1);}
            }
            team.team_barrier();
            if(ctr[GRAV_PACKET_CTR_FAILED]) {break;}

            /* The round boundary, and the only place the frontier changes shape.  What survives is
             * the published block the pops did not reach, and on top of it this round's new pushes
             * moved down to meet it -- so the stack stays newest-on-top and holds no gap.  When
             * nothing was popped the two blocks are already adjacent and nothing moves at all. */
            /* read before the compaction writes anything, so no thread is reading a counter
               another is settling */
            const int chunk_full = (ctr[GRAV_PACKET_CTR_RECORDS] >= chunk_cap);
            const int fr_pub    = ctr[GRAV_PACKET_CTR_FR_PUB];
            const int fr_taken  = (ctr[GRAV_PACKET_CTR_FR_TAKEN] < fr_pub) ? ctr[GRAV_PACKET_CTR_FR_TAKEN] : fr_pub;
            const int survivors = fr_pub - fr_taken;
            const int room      = frontier_cap - fr_pub;
            const int fr_new    = (ctr[GRAV_PACKET_CTR_FR_NEW] < room) ? ctr[GRAV_PACKET_CTR_FR_NEW] : room;
            if(fr_taken > 0 && fr_new > 0) {
                /* ONE thread, ASCENDING, and both halves of that are load-bearing.
                 *
                 * The source and destination blocks OVERLAP whenever more items were pushed than
                 * popped, which is the ordinary case: the slot this step writes,
                 * survivors + j, is the source of item j - fr_taken, an EARLIER index.  Ascending
                 * order is therefore safe by construction -- that earlier item has already been
                 * copied out before anything writes over it -- while a parallel copy has no order
                 * at all, so one lane can destroy a continuation before the lane that owns it has
                 * read it.  The subtree in the destroyed slot is then lost and the overwriting one
                 * is walked twice, with no failure raised, because both are perfectly well-formed
                 * items.  Measured that way: 79% of evrard's particles wrong at the first snapshot,
                 * and per-step interaction counts moving in BOTH directions.
                 *
                 * It costs one thread a run over at most the frontier once per round, against T
                 * walkers each taking up to steps_per_round node steps; if a round boundary ever
                 * shows up in a measurement, the fix is a ring buffer that never moves anything,
                 * not a parallel copy. */
                if(t == 0) {
                    for(int j = 0; j < fr_new; j++) {
                        frontier[survivors + j] = frontier[fr_pub + j];
                        mask_copy(fmasks + (size_t) (survivors + j) * W, fmasks + (size_t) (fr_pub + j) * W, W);
                    }
                }
            }
            team.team_barrier();

            if(t == 0) {
                ctr[GRAV_PACKET_CTR_FR_PUB] = survivors + fr_new;
                ctr[GRAV_PACKET_CTR_FR_TAKEN] = 0; ctr[GRAV_PACKET_CTR_FR_NEW] = 0;
                if(ctr[GRAV_PACKET_CTR_RECORDS] > chunk_cap) {ctr[GRAV_PACKET_CTR_RECORDS] = chunk_cap;}
                ctr[GRAV_PACKET_CTR_DONE] = (ctr[GRAV_PACKET_CTR_WORKING] == 0 && ctr[GRAV_PACKET_CTR_FR_PUB] == 0);
                /* A round in which nobody advanced and nobody was waiting on the chunk cannot
                 * happen: a walker holding an item steps at least once, and a team holding nothing
                 * is done.  It is tested rather than asserted because the alternative to noticing
                 * is a kernel that never returns, and the packet has a counted way out. */
                if(ctr[GRAV_PACKET_CTR_STEPPED] == 0 && !ctr[GRAV_PACKET_CTR_DONE] && !chunk_full) {
                    fail(ctr, GRAV_PACKET_FAIL_NO_PROGRESS);
                }
                ctr[GRAV_PACKET_CTR_STEPPED] = 0;
            }
            team.team_barrier();
            if(ctr[GRAV_PACKET_CTR_FAILED]) {break;}

            /* the chunk is full, or the traversal has finished: every member evaluates its records */
            if(ctr[GRAV_PACKET_CTR_RECORDS] > 0 && (chunk_full || ctr[GRAV_PACKET_CTR_DONE])) {
                /* the member's lanes take its records in turn: lane k of the member takes
                   records k, k + lanes, k + 2*lanes, ... that carry the member's bit */
                if(serves && in.open.alive) {
                    const int n_rec = ctr[GRAV_PACKET_CTR_RECORDS];
                    for(int r = sub_lane; r < n_rec; r += lanes) {
                        if(!mask_test(rmasks + (size_t) r * W, member)) {continue;}
                        /* A record this member cannot reproduce fails the whole packet, exactly as a
                           pseudo-particle does: nothing this team computed is
                           committed, and the replay walks every member again. */
                        if(!evaluate_record(records[r].no, in, sums)) {fail(ctr, GRAV_PACKET_FAIL_RECORD_UNUSABLE); break;}
                    }
                }
                team.team_barrier();
                if(t == 0) {ctr[GRAV_PACKET_CTR_RECORDS] = 0;}
            }
            team.team_barrier();
            if(ctr[GRAV_PACKET_CTR_FAILED]) {break;}
            if(ctr[GRAV_PACKET_CTR_DONE]) {break;}
            if(t == 0) {ctr[GRAV_PACKET_CTR_WORKING] = 0;}
            team.team_barrier();
        }

        /* commit, or discard everything */
        if(ctr[GRAV_PACKET_CTR_FAILED]) {
            /* Once per packet, not per member: the leader owns the tally, so the count is packets
             * given up rather than targets replayed, which is the quantity the budget question
             * asks.  d_failed is the replay's channel and is rewritten by it, so the reason needs
             * this one of its own. */
            if(t == 0 && d_fail_by_reason) {
                const int r = ctr[GRAV_PACKET_CTR_FAIL_REASON];
                Kokkos::atomic_fetch_add(&d_fail_by_reason[(r > 0 && r < GRAV_PACKET_FAIL_REASONS) ? r : 0], 1);
            }
            if(have_member && d_failed) {d_failed[first + t] = 1;}
            return;
        }
        /* Fold each member's partial sums into its first lane: every thread leaves its whole partial
           in team scratch, and the first lane of each member takes in the others in lane order.  A
           thread working on no member leaves the identity it was initialised to.  `lanes` is the same
           on every thread of the team, so the whole team takes this branch together or not at all. */
        if constexpr (Policy::team_evaluates) {
            if(lanes > 1) {
                gpu_grav_member_sums_t *partials = (gpu_grav_member_sums_t *) (scratch + plan.member_fold);
                partials[t] = sums;
                team.team_barrier();
                if(serves && sub_lane == 0) {for(int j = 1; j < lanes; j++) {gpu_grav_member_sums_combine(sums, partials[t + j], in.open.ptype);}}
            }
        }
        if(t == 0) {gpu_grav_note_commit(ctr[GRAV_PACKET_CTR_NOTE_INCOMPLETE], ctr[GRAV_PACKET_CTR_NOTE_UNSHIPPABLE]);}
        if(serves && sub_lane == 0) {
            Vec3<double> acc = Vec3<double>{0,0,0}; int ninter = 0; double pot = 0.0;
            if(in.open.alive) {gpu_grav_member_finish(ctx, in, sums, acc, ninter, pot);}
            d_acc[first + member] = acc; d_ninter[first + member] = ninter; d_pot[first + member] = pot; d_failed[first + member] = 0;
        }
    }
};

/* The shape the last primary walk on this rank used, in full, so a pricing arm records what it
 * measured rather than what it asked for. */
static struct gpu_grav_packet_shape_t g_packet_shape = {GRAV_PACKET_MODE_NONE, 0, 0, 0, 0, 0, 0, -1, -1, 0};

/* gravity_tree() writes it into the per-call timings record, so a run's readout says what was
 * measured -- including the table row, which is how a backend that could not launch the requested
 * shape is seen instead of silently changing the experiment. */
extern "C" void gpu_gravtree_packet_shape(struct gpu_grav_packet_shape_t *out) {*out = g_packet_shape;}

/* Packets given up on the last primary walk, by reason, alongside the shape: a cliff shows up as
 * a NO_CONTINUATION count that grows while the shape stays put, which a total could not show. */
/* The public slot count and the engine's reason list must not drift apart: the report indexes
 * slots by hand, so a new reason added without widening the header would be counted and never
 * printed. */
static_assert(GRAV_PACKET_FAIL_REASONS == GRAV_PACKET_FAIL_REASON_SLOTS,
              "gpu_gravtree.h's GRAV_PACKET_FAIL_REASON_SLOTS must match the engine's reason count");
static long long g_packet_fail[GRAV_PACKET_FAIL_REASONS] = {0};
extern "C" void gpu_gravtree_packet_failures(long long *out, int n)
{
    for(int r = 0; r < n; r++) {out[r] = (r < GRAV_PACKET_FAIL_REASONS) ? g_packet_fail[r] : 0;}
}
extern "C" int gpu_gravtree_packet_failure_reasons(void) {return GRAV_PACKET_FAIL_REASONS;}

/* Host: choose the launch shape and run the packet engine over the candidates.
 * Returns 0 on success (outputs and d_failed filled per candidate), 1 if no legal shape exists
 * for this build, in which case the caller uses the single-target walk.
 *
 * Only the shape decides the scratch: the frontier exists only for several walkers, and the
 * team-evaluation regions only for the team flavour -- the scratch request is what the legality
 * bound reads, so a region left in "because it is unused anyway" would still cost occupancy.
 *
 * The team is NOT the packet width: it is the row's, and the members are the configured packet
 * size when that is narrower, a configured size above the team being walked as several packets of
 * team-many. That ceiling is deliberate -- every member owns a thread, so no thread carries
 * several complete target states across rounds -- and it does mean a row chosen for occupancy
 * caps how many targets can share one traversal, which is why both numbers are reported. */
/* HOW THE DEVICE SPREADS ITS LANES OVER THIS CALL'S TARGETS.
 *
 * ⛔ This decides NOTHING about whether the device is used. That decision is made upstream from
 * the rank's active count, exactly as it always has been, and nothing here is consulted for it or
 * may change it. Every path below ends on the device; none of them can route work to the host.
 *
 * The question here is only which device schedule fits the call:
 *
 *   A WHOLE TEAM PER TARGET IS AFFORDABLE -- the cooperative schedule. Giving each of a handful of
 *     targets one lane leaves almost the whole device idle while those few lanes each grind down a
 *     long serial traversal, which is the defect the cooperative walk exists to remove. So each
 *     target gets a TEAM, and that team's lanes share its descent.
 *
 *   OTHERWISE -- the ordinary schedule. Each lane takes a target and walks it alone, covering
 *     several targets in turn when there are more targets than lanes.
 *
 * The key is structural -- lanes the machine has, against targets THIS RANK must walk -- not a
 * tuned target count, so the same code chooses on a laptop thread pool and on a GPU three orders
 * of magnitude wider:
 *
 *     cooperative  iff  lane_count / n_targets >= GRAV_COOP_MAX_LANES_PER_TARGET
 *     team              = GRAV_COOP_MAX_LANES_PER_TARGET, one target per team
 *
 * The ordinary schedule is one independent pointer-chasing walk per lane, so a wavefront's lanes
 * diverge across unrelated paths, and with fewer targets than lanes a call lasts as long as its
 * slowest target's walk. Where that is the case the masked packet schedule below, several targets
 * sharing ONE traversal, is faster; once a rank's targets reach about half the device's lanes the
 * independent walks keep the device busy and the ordinary schedule is the faster one, so dense
 * calls take it from there (gpu_grav_dense_walks_flat).
 *
 * ⛔ Capacity exhaustion inside a cooperative team is a RARE, COUNTED safety valve and must never
 * be how an ordinary call gets handled: a packet that gives up is re-walked by the device
 * single-target walk, which is the ordinary schedule, so correctness is safe and the cost is the
 * traversal it threw away. The production schedules are therefore required to show zero
 * give-ups; a deliberately undersized capacity is a gate-only positive control. */
enum gpu_grav_sched_mode_t {
    GRAV_SCHED_FLAT = 0,      /* one lane per target, the target's own serial walk: the fallback, and the replay of a packet that gives up */
    GRAV_SCHED_PACKET,        /* one walker per team, members sharing its traversal: the dense schedule */
    GRAV_SCHED_COOPERATIVE    /* one target per team, the team's lanes sharing that target's traversal */
};

/* A schedule is the WHOLE shape, because its pieces are not independent: a wide team is what gives
 * a target with few companions more than one traversing lane, but it also multiplies the
 * per-walker continuation storage, and that scratch is paid for in teams resident per compute
 * unit. The capacities below are a starting bracket, not a measurement -- the B3 pricing campaign
 * sets them, and reads `steps_per_round` (how long a walker runs before what it published becomes
 * visible to an idle one) over {16, 32, 64, 128}. */
struct gpu_grav_sched_row_t {
    enum gpu_grav_sched_mode_t mode;
    int team;              /* threads in the team */
    int n_walkers;         /* of those, how many traverse */
    int frontier_mul;      /* frontier items per walker */
    int chunk;             /* records held between flushes */
    int steps_per_round;
};

/* The cooperative schedules, widest first.
 *
 * ⛔ ONLY THE FIRST ROW IS SELECTED IN PRODUCTION, and deliberately so: a partially wide team is
 * the shape that was MEASURED TO LOSE. On a 128-rank clustered zoom at the stock routing
 * threshold, ranks taking the device route hold ~1e4 targets against ~2e5 lanes -- about 20 lanes
 * per target, so a team of ~16 -- and on the six steps where that happened the imbalance row cost
 * 101.8 s where the same steps without it cost 2-4 s apiece. Cooperation there buys a little
 * traversal sharing on a device that was already fed, pays the barrier and frontier for it, and
 * gives up the packet sharing (q_dev drops to 1) into the bargain; and because the schedule is
 * chosen per rank, only the handful of ranks that qualify slow down, so the whole of it lands in
 * the imbalance row while the other hundred ranks wait.
 *
 * The narrower rows therefore stay here as the shape a MEASURED dispatch table would fill in, and
 * as the only way to exercise cooperation on a host backend whose lane count is its thread count.
 * They are not reachable from the production criterion, and a team the backend cannot launch falls
 * to GRAV_SCHED_FLAT -- still on the device -- rather than quietly narrowing, which would answer a
 * pricing arm with a shape nobody chose. */
static const struct gpu_grav_sched_row_t g_grav_coop_rows[] = {
    {GRAV_SCHED_COOPERATIVE, 64, 64, 4, 256, 32},
    {GRAV_SCHED_COOPERATIVE, 32, 32, 4, 256, 64},
    {GRAV_SCHED_COOPERATIVE, 16, 16, 4, 256, 64},
    {GRAV_SCHED_COOPERATIVE,  8,  8, 4, 256, 64},
    {GRAV_SCHED_COOPERATIVE,  4,  4, 4, 256, 64},
    {GRAV_SCHED_COOPERATIVE,  2,  2, 4, 256, 64},
};
static const int g_grav_coop_n_rows = (int) (sizeof(g_grav_coop_rows) / sizeof(g_grav_coop_rows[0]));

/* The most lanes worth putting on one target's traversal. Beyond some width the walkers spend
 * more of the round waiting at the barrier and contending for the frontier than they save on the
 * descent, and the shared frontier has to hold a continuation for every one of them. Internal,
 * and the pricing campaign's to set. */
#define GRAV_COOP_MAX_LANES_PER_TARGET 64

/* A DENSE call -- one where the device cannot give every target a whole team -- walks as MASKED
 * PACKETS while the rank's targets are fewer than about half the device's lanes: TREE_QUERY_PACKET_SIZE
 * neighbouring targets share one traversal, each still judging every node for itself.  On a 128-rank
 * zoom, with and without the FIRE physics, that is 27-45% faster than one independent traversal per
 * lane at 1e2-1e4 targets a rank.
 *
 * At or above half the lanes the call takes the single-target walk instead.  Independent per-lane
 * walks then already keep the device busy, so sharing a traversal stops paying for its bookkeeping:
 * on a 128-rank FIRE zoom (MI250X) the single-target walk was 25-30% faster on fully active calls and
 * level with packets just below the threshold, where a call's time is set by its slowest target's
 * walk rather than by the number of targets.  The single-target walk is also the schedule when a
 * packet row cannot be launched, and the replay for a packet that gives up. */
static int gpu_grav_dense_walks_flat(int n_cand)
{
    return (2LL * (long long) n_cand >= (long long) gizmo_gpu_lane_count()) ? 1 : 0;
}

/* The cooperative row this call takes, or -1 for the ordinary schedule.
 *
 * The test is FULL-TEAM ELIGIBILITY, not a fitted crossover: cooperate only where every target can
 * be given a whole team. It is stated that way because it is the honest reading of what is known.
 * A rank with ~20 lanes per target was measured to lose badly, and above the cap the team size is
 * the cap whatever the ratio is -- so raising the bar past the cap would change no team allocation
 * anywhere, only exclude a middle band that no per-rank measurement has yet adjudicated. When one
 * does, this becomes a table with narrower rows in it.
 *
 * ⛔ The ratio is this RANK's, from its own target count, because the schedule is this rank's. A
 * global active count divided by the rank count is not this quantity and must not be substituted
 * for it: the ranks are precisely what is inhomogeneous here. */
/* ⚠ The test below asks for a whole team per TARGET, while a cooperative team now serves a whole
 * PACKET of them -- so it demands more lanes than the shape strictly needs, by the packet size.
 * That is deliberate for now: it engages cooperation only where cooperation was already measured
 * to win, and adds the member sharing on top. Relaxing it toward what the shape actually costs
 * (a full team per PACKET) is a MEASUREMENT, not an adjustment, and the last time a criterion here
 * was chosen by argument it cost 42 s of imbalance. */
static int gpu_grav_coop_first_row(int n_cand)
{
    if(n_cand <= 0) {return -1;}
    if(gizmo_gpu_lane_count() / n_cand < GRAV_COOP_MAX_LANES_PER_TARGET) {return -1;}
    return 0;
}

/* Run one schedule row, or report that the backend cannot launch it.
 *
 * Returns 0 on success, 1 when this row is not launchable here -- which is a statement about the
 * row and the backend only. What the caller does about it (the next named row, then the ordinary
 * device schedule) is the caller's, and in no case is it to send work to the host. */
template <class Policy>
static int gpu_grav_packet_launch_row(GpuGravPacketWalk<Policy> &f, const struct gpu_grav_sched_row_t &r,
                                      int row_requested, int row_effective, const char *kernel_name,
                                      struct gpu_grav_packet_shape_t *shape_out)
{
    {
        int team = r.team;
        if(team > GRAV_PACKET_Q_DEV_MAX) {team = GRAV_PACKET_Q_DEV_MAX;}
        /* Members per packet: the configured packet size, or the team when that is narrower
         * (78b.1), for EVERY evaluating schedule including the cooperative one.
         *
         * The cooperative schedule is a team of walkers on a packet, not on a single target. Both
         * levels are wanted and they compose: the members share ONE traversal of the tree and then
         * evaluate its records in parallel, while the team's threads share the descent that
         * traversal makes. Forcing one member here -- which this did -- discards the first level
         * entirely, so the tree was walked once per target rather than once per packet and the
         * flush ran on a single lane of the team. */
        f.q_dev = (Policy::packet_size < team) ? Policy::packet_size : team;
        f.n_walkers = (r.n_walkers < team) ? r.n_walkers : team;
        f.steps_per_round = r.steps_per_round;
        f.frontier_cap = (f.n_walkers > 1)
                             ? ((r.frontier_mul * team > 16) ? r.frontier_mul * team : 16) : 0;
        f.chunk_cap    = r.chunk;
        f.plan = gpu_grav_packet_scratch_plan(f.q_dev, team, f.frontier_cap, f.chunk_cap, Policy::local_stack, Policy::team_evaluates);
        /* the legality bound is asked of a probe policy carrying the same scratch request: a
           policy constructed at an illegal team size throws before it can be asked anything */
        Kokkos::TeamPolicy<> probe(1, 1, 1);
        probe.set_scratch_size(0, Kokkos::PerTeam(f.plan.bytes));
        const int hw = probe.team_size_max(f, Kokkos::ParallelForTag());
        if(hw < team) {return 1;}   /* this row does not fit this backend */
        const int league = (f.n_cand + f.q_dev - 1) / f.q_dev;
        Kokkos::TeamPolicy<> policy(league, team, 1);
        policy.set_scratch_size(0, Kokkos::PerTeam(f.plan.bytes));
        Kokkos::parallel_for(kernel_name, policy, f);
        Kokkos::fence();
        gizmo_gpu_check_last_error(kernel_name, league);
        if(shape_out) {
            shape_out->mode = (int) r.mode;
            shape_out->team = team; shape_out->q_dev = f.q_dev;
            shape_out->n_walkers = f.n_walkers; shape_out->frontier = f.frontier_cap;
            shape_out->chunk = f.chunk_cap; shape_out->steps_per_round = f.steps_per_round;
            shape_out->row_requested = row_requested; shape_out->row_effective = row_effective;
            shape_out->scratch_bytes = (long long) f.plan.bytes;
        }
        return 0;
    }
}

/* Choose and run a device schedule for this call's targets.
 *
 * Returns 0 when the engine ran and filled the outputs, 1 when the call should take the ORDINARY
 * device schedule -- the single-target walk the caller already has. ⛔ A 1 here never means "use
 * the host": the host/device decision was taken before this function was reached and is not
 * revisited by it. */
template <class Policy>
static int gpu_grav_packet_launch(GpuGravPacketWalk<Policy> &f, const char *kernel_name,
                                  struct gpu_grav_packet_shape_t *shape_out)
{
    if constexpr (Policy::team_evaluates) {
        /* The flavour for a team wider than its packet is taken only where the criterion gives
         * every target a whole team, so this is the cooperative row.  ONE attempt, at the row the
         * criterion chose. There is deliberately no step-down to a narrower cooperative team: a
         * partial team is the shape that was measured to lose, so a backend that cannot launch
         * the full one takes the ordinary device schedule instead. */
        const int first_row = gpu_grav_coop_first_row(f.n_cand);
        if(first_row < 0) {return 1;}
        return gpu_grav_packet_launch_row(f, g_grav_coop_rows[first_row], first_row, first_row,
                                          kernel_name, shape_out);
    } else {
        /* Not enough lanes for a whole team per target, but fewer targets than half the lanes:
         * MASKED PACKETS, where TREE_QUERY_PACKET_SIZE targets adjacent in the active list, and
         * therefore adjacent in space, share ONE traversal while each still judges every node for
         * itself.  One independent traversal per lane instead has every lane chasing its own
         * pointer chain, so a wavefront's lanes diverge across unrelated paths and the node loads
         * they have in common are never shared.  From half the lanes up the call returns to the
         * ordinary schedule (gpu_grav_dense_walks_flat). */
        if(Policy::packet_size > 1 && !gpu_grav_dense_walks_flat(f.n_cand)) {
            const int t = (Policy::packet_size < GRAV_PACKET_Q_DEV_MAX) ? Policy::packet_size
                                                                        : GRAV_PACKET_Q_DEV_MAX;
            const struct gpu_grav_sched_row_t dense = {GRAV_SCHED_PACKET, t, 1, 0, 256, 64};
            return gpu_grav_packet_launch_row(f, dense, -1, -1, kernel_name, shape_out);
        }
        return 1;
    }
}

template <class Policy>
static int gpu_gravtree_walk_packets_as(const gpu_grav_walk_ctx_t &ctx, const int *d_idx, int n_cand,
                                        Vec3<double> *d_acc, int *d_ninter, double *d_pot, int *d_failed,
                                        int *d_fail_by_reason, const char *kernel_name)
{
    GpuGravPacketWalk<Policy> f;
    f.ctx = ctx; f.d_idx = d_idx; f.n_cand = n_cand;
    f.d_acc = d_acc; f.d_ninter = d_ninter; f.d_pot = d_pot; f.d_failed = d_failed;
    f.d_fail_by_reason = d_fail_by_reason;
    return gpu_grav_packet_launch(f, kernel_name, &g_packet_shape);
}

/* Where every target can be given a whole team, the team evaluates each member's records over
 * several lanes; elsewhere each member has one thread.  The two are separate kernels so that each
 * is compiled for the evaluation it does (GravPacketMaskedTeamPolicy says why). */
static int gpu_gravtree_walk_packets(const gpu_grav_walk_ctx_t &ctx, const int *d_idx, int n_cand,
                                     Vec3<double> *d_acc, int *d_ninter, double *d_pot, int *d_failed,
                                     int *d_fail_by_reason)
{
    if(gpu_grav_coop_first_row(n_cand) >= 0) {
        return gpu_gravtree_walk_packets_as<GravPacketMaskedTeamPolicy>(ctx, d_idx, n_cand, d_acc, d_ninter, d_pot, d_failed,
                                                                        d_fail_by_reason, "gravtree_walk_packets_team");
    }
    return gpu_gravtree_walk_packets_as<GravPacketMaskedPolicy>(ctx, d_idx, n_cand, d_acc, d_ninter, d_pot, d_failed,
                                                                d_fail_by_reason, "gravtree_walk_packets");
}

#ifdef PMGRID
/* SharedSpace mirrors of the PM short-range lookup tables. Seeded once and kept
 * for the run: force_treeallocate() rebuilds the host tables from
 * u = 3.0/NTAB*(i+0.5) via erfc/exp, a pure function of the table index, so
 * every rebuild writes identical values and the device copy never goes stale.
 * Same acquire-once shape as gpu_ewald_tables_acquire() below. */
static float *g_d_shortrange_tab       = NULL;
static float *g_d_shortrange_pot_tab   = NULL;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
static float *g_d_shortrange_tidal_tab = NULL;
#endif
static int    g_shortrange_tables_ready = 0;

static int gpu_shortrange_tables_acquire(void)
{
    if(g_shortrange_tables_ready) return 0;
    const size_t sz = GIZMO_GPU_GRAVTREE_NTAB * sizeof(float);
    g_d_shortrange_tab = (float *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    g_d_shortrange_pot_tab = (float *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    if(!g_d_shortrange_tab || !g_d_shortrange_pot_tab) {
        printf("gpu_shortrange_tables_acquire: shortrange table alloc failed\n");
        endrun(913205);
        return 1;   /* soft bad-stop: caller drains without walking */
    }
    memcpy(g_d_shortrange_tab,     shortrange_table,           sz);
    memcpy(g_d_shortrange_pot_tab, shortrange_table_potential, sz);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    g_d_shortrange_tidal_tab = (float *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    if(!g_d_shortrange_tidal_tab) {
        printf("gpu_shortrange_tables_acquire: shortrange tidal table alloc failed\n");
        endrun(913206);
        return 1;
    }
    memcpy(g_d_shortrange_tidal_tab, shortrange_table_tidal, sz);
#endif
    g_shortrange_tables_ready = 1;
    return 0;
}
#endif

/* Fused per-call scratch for the gravity walks.
 *
 * Both walks need several per-target arrays for the duration of one call.
 * Allocated separately, each array costs its own device allocation, free and
 * page advice on every call; carved from one block they cost a single set.
 * The lifetime is unchanged: the block is acquired at the top of the walk and
 * released before it returns.
 *
 * This plan is the only place the sizes and offsets are computed, so the
 * allocation size and the carved pointers cannot drift apart. Each offset is
 * aligned explicitly for its element type rather than relying on the field
 * order. `with_potential_and_interactions` selects the primary walk's full
 * set; the Ewald correction walk needs only the target index, the failure
 * flag and the acceleration. */
/* Offset value for a member this plan does not carry. */
#define GRAV_WALK_SCRATCH_ABSENT ((size_t) -1)

struct grav_walk_scratch_plan
{
    size_t bytes;
    size_t acc, pot, idx, failed, ninter;
    size_t fail_by_reason;   /* [GRAV_PACKET_FAIL_REASONS], written by the engine's team leaders */
    size_t time_fault;       /* [1] gpu_grav_time_fault_t: sources found ahead of the walk time */
};

static inline size_t grav_walk_scratch_align(size_t offset, size_t alignment)
{
    return ((offset + alignment - 1) / alignment) * alignment;
}

static struct grav_walk_scratch_plan
grav_walk_scratch_plan_for(int num_targets, int with_potential_and_interactions)
{
    struct grav_walk_scratch_plan plan;
    const size_t n = (size_t) num_targets;
    size_t offset = 0;

    offset = grav_walk_scratch_align(offset, alignof(Vec3<double>));
    plan.acc = offset;  offset += n * sizeof(Vec3<double>);

    /* Poison, not 0: zero is the valid offset of `acc`, so an absent member
     * carved by mistake would alias the acceleration array silently. */
    plan.pot = plan.ninter = GRAV_WALK_SCRATCH_ABSENT;
    if(with_potential_and_interactions) {
        offset = grav_walk_scratch_align(offset, alignof(double));
        plan.pot = offset;  offset += n * sizeof(double);
    }
    offset = grav_walk_scratch_align(offset, alignof(int));
    plan.idx = offset;  offset += n * sizeof(int);
    offset = grav_walk_scratch_align(offset, alignof(int));
    plan.failed = offset;  offset += n * sizeof(int);
    if(with_potential_and_interactions) {
        offset = grav_walk_scratch_align(offset, alignof(int));
        plan.ninter = offset;  offset += n * sizeof(int);
    }
    /* Not per target: one tally for the whole launch, carved from the same block so it needs no
     * allocation of its own and is host-readable the moment the kernel has been waited on. */
    offset = grav_walk_scratch_align(offset, alignof(int));
    plan.fail_by_reason = offset;  offset += (size_t) GRAV_PACKET_FAIL_REASONS * sizeof(int);
    offset = grav_walk_scratch_align(offset, alignof(struct gpu_grav_time_fault_t));
    plan.time_fault = offset;  offset += sizeof(struct gpu_grav_time_fault_t);
    plan.bytes = offset;
    return plan;
}

/* Every device walk predicts nodes from their MIRROR, so the mirror must describe each node's state
 * at that node's own time.  A host lazy drift advances a node and leaves its mirror behind, and
 * claims it; the claims are answered here, before every device walk, whatever else has certified
 * the tree -- O(claims), nothing when there are none.  When the claim record cannot say which nodes
 * they are (it overflowed, or its fail-safe fired) the full mirror refresh answers instead: the only
 * sweep that also rewrites nodes already at this time.  Never expected; reported.  Returns nonzero
 * when not even that succeeded, after requesting the stop: the caller must not walk. */
static int gpu_grav_answer_node_claims(integertime ti, const char *walk)
{
    if(gpu_node_dirty_bring_gravity_current(ti) == 0) {return 0;}
    static int reported_claim_fallback = 0;
    if(!reported_claim_fallback) {
        reported_claim_fallback = 1;
        printf("%s: task %d: the node claim record could not name the host-drifted nodes; refreshing every node mirror instead (reported once)\n", walk, ThisTask);
        fflush(stdout);
    }
    if(gpu_force_drift_nodes_ex(ti, /*refresh_mirrors_already_current=*/1) == 0) {return 0;}
    endrun(929702);
    return 1;
}

extern "C" int gpu_gravtree_walk_primary(int *host_candidates_left)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();
    int num_active_total = (int) ActiveParticleList.size();
    /* Cleared before any early return, so the shape reported for THIS call is this call's: a
     * host-routed or empty call reports 0/0 rather than whatever the previous device call ran.
     * The give-up tally is cleared with it and for the same reason -- it was cleared only where
     * the engine runs, so a host-routed step re-reported, and re-reduced, the previous device
     * call's counts against gpu_gravtree.h's promise that they read zero on such a call. */
    {struct gpu_grav_packet_shape_t cleared = {GRAV_PACKET_MODE_NONE, 0, 0, 0, 0, 0, 0, -1, -1, 0}; g_packet_shape = cleared;}
    for(int r = 0; r < GRAV_PACKET_FAIL_REASONS; r++) {g_packet_fail[r] = 0;}
    /* How many candidates this walk leaves to the host loop: every active until this walk has
     * selected and taken some. The host loop sizes its per-thread packet workspace from it. */
    if(host_candidates_left) {*host_candidates_left = (num_active_total > 0) ? num_active_total : 0;}
    if(Ewald_iter > 0) {return 0;}
    if(num_active_total <= 0) {return 0;}

    /* The CPU walk (forcetree.cc) drifts a particle or a node whose Ti_current is stale at the
     * moment it opens it. This walk drifts nothing: the targets are the active set, drifted here
     * before the routing decision so that the candidacy test and both routes read current
     * targets, and every source it reaches -- particle or node -- is read as it stands at the walk
     * time, predicted there without being written back (gpu_grav_particle_source_at,
     * gpu_grav_node_motion_at). */
    /* Host-side wrapper in the GPU TU must use the out-of-line host accessor
     * `gizmo_host_ti_current()` (defined in core/predict.cc) rather than a
     * bare All.Ti_Current read, so the host-snapshot intent at this call
     * site stays correct even when the device-pass redirect is active. */
    integertime ti_curr_host = gizmo_host_ti_current();
    move_particles(ti_curr_host); /* drifts the ACTIVE set; every other particle is still at its own Ti_current */
    /* The SoA must exist before anything below reads or repairs the node mirror. */
    gpu_gravity_tree_acquire(MaxNodes + 1, Nodes_base, Extnodes_base);

    /* Select the particles this walk will actually cover BEFORE any device work: a step that
     * turns out to walk nothing, or that is small enough to be cheaper on the host, must not pay
     * for it.  The selection reads only particle state (ActiveParticleList, ProcessedFlag and the
     * candidacy predicate), all of which move_particles above has already brought to
     * ti_curr_host. */
    int *idx_host = (int *) mymalloc("gpu_grav_idx", num_active_total * sizeof(int));
    int num_active = 0;
    for(int a = 0; a < num_active_total; a++) {
        int i = ActiveParticleList[a];
        if(ProcessedFlag[i]) {continue;}
        /* SSOT pre-walk candidacy (Mass>0 + Hermite eligibility + needs_new_treeforce):
         * the GPU pre-pass must use the same candidate set as the CPU primary walk +
         * finalization, or it leaves a non-candidate's GravAccel fresh-written but raw.
         * ProcessedFlag is intentionally left unset on a candidacy skip here (matches
         * the CPU primary walk: a cached/extrapolated particle is finalized later). */
        if(!gravity_treewalk_candidate_prewalk(i, a)) {continue;}
        idx_host[num_active++] = i;
    }
    if(host_candidates_left) {*host_candidates_left = num_active;}
    if(num_active <= 0) {myfree(idx_host); return 0;}

    /* Few enough candidates that the host walk is the cheaper route (gravity_walk_route_to_host).
     * Returning with ProcessedFlag untouched leaves every candidate to the host loop in gravtree.cc. */
    if(gravity_walk_route_to_host(num_active)) {myfree(idx_host); return 0;}

    /* Acquire the arena (P_dev + CellP_dev in SharedSpace) */
    gpu_particles_arena_set_site("gpu_gravtree_walk_primary");
    gpu_particles_arena_acquire(NumPart, P, CellP);
    struct particle_data    *P_dev    = gpu_particles_arena_P();
    struct gas_cell_data    *CellP_dev = gpu_particles_arena_CellP();

    int min_nodes = MaxNodes + 1;
    gpu_gravity_tree_acquire(min_nodes, Nodes_base, Extnodes_base);
    /* soa->nextnode_aux aliases UVM Nextnode[] (set by
     * force_treeallocate); no per-walk memcpy needed. */
    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    /* CellP is legitimately NULL on a gas-free (DM-only) problem (TotN_gas==0).
     * Every CellP_dev use in this walk + post-walk scatter is gas-gated
     * (device: valid_gas_particle_for_rt / ptype==0; host scatter: Type==0),
     * so a null CellP_dev is safe there. Only require it when gas exists. */
    const bool need_cellp = (All.TotN_gas > 0);   /* host read; bare All.* is safe in GPU-TU host code (all-mirror) */
    if(!P_dev || !soa || (need_cellp && !CellP_dev)) {
        printf("gpu_gravtree_walk_primary: failed to acquire arena or tree SoA\n");
        endrun(913200);
        myfree(idx_host);   /* LIFO mymalloc cleanup before drain */
        return 1;
    }

    int treeBase = All.TreeNodeIndexBase;
    int treeParticleSlots_snap = All.TreeParticleSlots;
    int maxNodes_snap = MaxNodes;
    int maxForeignNodes_snap = MaxForeignNodes;    /* LET */
    const struct gpu_gravity_tree_soa_t soa_snap = *soa;
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    /* host-evaluate the first-step predicate once; captured by value into the device walk */
    int is_first_step_snap = (All.Ti_Current == 0 && RestartFlag != 1);
#endif

#ifdef PMGRID
    double rcut_snap     = All.Rcut[0];
    double rcut2_snap    = rcut_snap * rcut_snap;
    double asmthfac_snap = 0.5 / All.Asmth[0] * (GIZMO_GPU_GRAVTREE_NTAB / 3.0);
    /* shortrange_table is a host global (forcetree.cc); the SharedSpace mirror the
     * kernel reads is seeded once per run, not per call. */
    if(gpu_shortrange_tables_acquire() != 0) {
        myfree(idx_host);
        return 1;
    }
#endif
    /* read-only PM short-range config captured by value into the device walk (empty when
     * !PMGRID; per-target PLACEHIGHRESREGION override happens inside the walk on its copy). */
    grav_pm_shortrange_t pm_snap{};
#ifdef PMGRID
    pm_snap.rcut = rcut_snap; pm_snap.rcut2 = rcut2_snap; pm_snap.asmthfac = asmthfac_snap;
    pm_snap.shortrange_tab = g_d_shortrange_tab;
#ifdef EVALPOTENTIAL
    pm_snap.shortrange_pot_tab = g_d_shortrange_pot_tab;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    pm_snap.shortrange_tidal_tab = g_d_shortrange_tidal_tab;
#endif
#endif

    /* Scratch arrays for per-target results, carved from one allocation */
    const struct grav_walk_scratch_plan scratch = grav_walk_scratch_plan_for(num_active, 1);
    char *scratch_block = (char *) gizmo_gpu_alloc_shared(scratch.bytes, "gravity_walk");
    if(!scratch_block) {
        printf("gpu_gravtree_walk_primary: kokkos_malloc failed\n");
        endrun(913201);
        myfree(idx_host);   /* LIFO mymalloc cleanup before drain */
        return 1;
    }
    int          *d_idx    = (int *)          (scratch_block + scratch.idx);
    int          *d_failed = (int *)          (scratch_block + scratch.failed);
    Vec3<double> *d_acc    = (Vec3<double> *) (scratch_block + scratch.acc);
    int          *d_ninter = (int *)          (scratch_block + scratch.ninter);
    double       *d_pot    = (double *)       (scratch_block + scratch.pot);
    int          *d_fail_by_reason = (int *) (scratch_block + scratch.fail_by_reason);
    memcpy(d_idx, idx_host, num_active * sizeof(int));
    memset(d_failed, 0, num_active * sizeof(int));
    memset(d_fail_by_reason, 0, GRAV_PACKET_FAIL_REASONS * sizeof(int));
    for(int r = 0; r < GRAV_PACKET_FAIL_REASONS; r++) {g_packet_fail[r] = 0;}

    /* The tables a prediction reads, refreshed once for the call: every source is read at the walk
     * time from wherever it stands, and the Hermite source prediction shares them. */
    struct DriftKickTableView tables_snap;
    if(drift_kick_table_mirror_refresh(&tu_drift_kick_table_dev, &tables_snap) != 0) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
        myfree(idx_host);   /* soft bad-stop already requested: no launch without the tables a prediction reads */
        return 1;
    }
    struct gpu_grav_time_fault_t *time_fault = (struct gpu_grav_time_fault_t *) (scratch_block + scratch.time_fault);
    memset(time_fault, 0, sizeof(struct gpu_grav_time_fault_t));

    gpu_grav_walk_ctx_t ctx{};
    ctx.treeBase = treeBase; ctx.treeParticleSlots = treeParticleSlots_snap;
    ctx.maxNodes = maxNodes_snap; ctx.maxForeignNodes = maxForeignNodes_snap;
    ctx.P_dev = P_dev; ctx.CellP_dev = CellP_dev; ctx.ti = ti_curr_host; ctx.tree_soa = soa_snap;
    ctx.Extnodes = Extnodes; ctx.tables = tables_snap; ctx.time_fault = time_fault;
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    ctx.is_first_step = is_first_step_snap;
#endif
    ctx.pm = pm_snap;

    if(gpu_grav_answer_node_claims(ti_curr_host, "gpu_gravtree_walk_primary") != 0) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
        myfree(idx_host);   /* LIFO mymalloc cleanup before drain */
        return 1;   /* soft bad-stop: never walk a mirror that could not be brought up to date; drains at next poll */
    }


    /* Per-particle gravity source inputs (RT luminosity / sink bolometric luminosity /
     * CR injection), evaluated via the shared SSOT helper gravtree_fill_particle_source_
     * inputs(). Two modes, same physics:
     *   EAGER: one host pass over all NumPart fills dense SharedSpace arrays the kernel
     *          reads by particle id.
     *   LAZY:  no dense arrays and no O(NumPart) pass -- the kernel evaluates the same
     *          helper on-device at each local particle-open, so the cost scales with the
     *          active set rather than NumPart.
     * Selected by the active-count threshold below; the result is identical either way. */
    bool use_lazy_source = false;
#if defined(GRAVTREE_SOURCE_LAZY_SUPPORTED) && (defined(RT_USE_GRAVTREE) || defined(SINK_PHOTONMOMENTUM) || defined(COSMIC_RAY_SUBGRID_LEBRON))
    /* Use lazy when the step touches few particles -- either an absolute count or a small
     * fraction of the local pool. Hard-coded (no env var, no parameter). */
    {
        const long   GRAVTREE_SOURCE_LAZY_CAP  = 256;
        const double GRAVTREE_SOURCE_LAZY_FRAC = 0.01;
        use_lazy_source = ((long)num_active <= GRAVTREE_SOURCE_LAZY_CAP) ||
                          ((double)num_active < GRAVTREE_SOURCE_LAZY_FRAC * (double)NumPart);
    }
#endif

#ifdef RT_USE_GRAVTREE
    MyFloat *d_src_lum = NULL;
#ifdef CHIMES_STELLAR_FLUXES
    double *d_src_lum_G0 = NULL, *d_src_lum_ion = NULL;
#endif
    if(!use_lazy_source) {
        long sz = (long)NumPart * N_RT_FREQ_BINS * sizeof(MyFloat);
        d_src_lum = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
        if(!d_src_lum) {printf("gpu_gravtree_walk_primary: d_src_lum alloc failed\n"); endrun(913202); myfree(idx_host); return 1;}
        memset(d_src_lum, 0, sz);
#ifdef CHIMES_STELLAR_FLUXES
        long szc = (long)NumPart * CHIMES_LOCAL_UV_NBINS * sizeof(double);
        d_src_lum_G0  = (double *) gizmo_gpu_alloc_shared(szc, "gravity_walk");
        d_src_lum_ion = (double *) gizmo_gpu_alloc_shared(szc, "gravity_walk");
        if(!d_src_lum_G0 || !d_src_lum_ion) {printf("gpu_gravtree_walk_primary: CHIMES lum alloc failed\n"); endrun(913203); myfree(idx_host); return 1;}
        memset(d_src_lum_G0,  0, szc);
        memset(d_src_lum_ion, 0, szc);
#endif
    }
#endif /* RT_USE_GRAVTREE */

#ifdef SINK_PHOTONMOMENTUM
    MyFloat       *d_bh_lum   = NULL;
    Vec3<MyFloat> *d_bh_angle = NULL;
    if(!use_lazy_source) {
        long sz_lum  = (long)NumPart * sizeof(MyFloat);
        long sz_ang  = (long)NumPart * sizeof(Vec3<MyFloat>);
        d_bh_lum   = (MyFloat *)       gizmo_gpu_alloc_shared(sz_lum, "gravity_walk");
        d_bh_angle = (Vec3<MyFloat> *) gizmo_gpu_alloc_shared(sz_ang, "gravity_walk");
        if(!d_bh_lum || !d_bh_angle) {printf("gpu_gravtree_walk_primary: bh_lum alloc failed\n"); endrun(913210); myfree(idx_host); return 1;}
        memset(d_bh_lum,   0, sz_lum);
        memset(d_bh_angle, 0, sz_ang);
    }
#endif /* SINK_PHOTONMOMENTUM */

#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyFloat *d_cr_inject = NULL;
    /* per-step CR-age scalar: needed by the walk's CR gate in BOTH modes (cr_active_gate
     * + grav_cr_lebron_accumulate), so it is computed unconditionally, NOT inside the
     * eager-only dense-array block. Lazy skips only the d_cr_inject ARRAY. */
    double   t_max_cr    = 0.0;
    if(All.Time > All.TimeBegin) {
        double t_gyr = evaluate_time_since_t_initial_in_Gyr(All.TimeBegin);
        if(t_gyr > 1.0) {t_gyr = 1.0;}
        t_max_cr = t_gyr / UNIT_TIME_IN_GYR;     /* per-step scalar; computed once, not per-particle */
    }
    if(!use_lazy_source) {
        long sz = (long)NumPart * sizeof(MyFloat);
        d_cr_inject = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
        if(!d_cr_inject) {printf("gpu_gravtree_walk_primary: cr_inject alloc failed\n"); endrun(913211); myfree(idx_host); return 1;}
        memset(d_cr_inject, 0, sz);
    }
#endif /* COSMIC_RAY_SUBGRID_LEBRON */

    /* EAGER only: single host pass over all NumPart, gated physics in the shared SSOT
     * helper, copying ONLY active entries into the bulk-zeroed SharedSpace arrays. In
     * LAZY mode this whole O(NumPart) pass is skipped -- the kernel evaluates the same
     * helper on-device at each local particle-open instead. */
#if defined(RT_USE_GRAVTREE) || defined(SINK_PHOTONMOMENTUM) || defined(COSMIC_RAY_SUBGRID_LEBRON)
    if(!use_lazy_source)
    for(int p = 0; p < NumPart; p++) {
        struct gravtree_source_inputs_t in;
        gravtree_fill_particle_source_inputs(p, P, CellP, &in);
#ifdef RT_USE_GRAVTREE
        if(in.rt_active) {
            int kf;
            for(kf = 0; kf < N_RT_FREQ_BINS; kf++) {d_src_lum[(long)p * N_RT_FREQ_BINS + kf] = in.src_lum[kf];}
#ifdef CHIMES_STELLAR_FLUXES
            for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++) {
                d_src_lum_G0[(long)p * CHIMES_LOCAL_UV_NBINS + kf]  = in.src_lum_G0[kf];
                d_src_lum_ion[(long)p * CHIMES_LOCAL_UV_NBINS + kf] = in.src_lum_ion[kf];
            }
#endif
        }
#endif
#ifdef SINK_PHOTONMOMENTUM
        if(in.bh_active) {
            d_bh_lum[p]   = in.bh_lum;
            d_bh_angle[p] = in.bh_angle;
        }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        if(in.cr_inject != 0) {d_cr_inject[p] = in.cr_inject;}
#endif
    }
#endif

#ifdef RT_USE_GRAVTREE
    struct gpu_rt_walk_data_t rt_data_snap;
    rt_data_snap.src_lum = d_src_lum;
#ifdef CHIMES_STELLAR_FLUXES
    rt_data_snap.src_lum_G0  = d_src_lum_G0;
    rt_data_snap.src_lum_ion = d_src_lum_ion;
#endif
#endif /* RT_USE_GRAVTREE */
#ifdef SINK_PHOTONMOMENTUM
    struct gpu_sink_walk_data_t sink_data_snap;
    sink_data_snap.bh_lum   = d_bh_lum;
    sink_data_snap.bh_angle = d_bh_angle;
#endif /* SINK_PHOTONMOMENTUM */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    struct gpu_cr_walk_data_t cr_data_snap;
    cr_data_snap.cr_inject = d_cr_inject;
    cr_data_snap.t_max_cr  = t_max_cr;
#endif /* COSMIC_RAY_SUBGRID_LEBRON */

    /* Release the optional per-call payload buffers. Defined once so the early
     * returns below and the normal exit path cannot drift apart: endrun() only
     * requests a controlled stop and returns, so a return that skips these
     * leaks device memory on every call until the stop is polled -- precisely
     * when device memory is already short. */
    auto release_payload_buffers = [&]() {
#ifdef RT_USE_GRAVTREE
#ifdef CHIMES_STELLAR_FLUXES
        if(d_src_lum_ion) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_src_lum_ion); d_src_lum_ion = NULL;}
        if(d_src_lum_G0)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_src_lum_G0);  d_src_lum_G0  = NULL;}
#endif
        if(d_src_lum) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_src_lum); d_src_lum = NULL;}
#endif
#ifdef SINK_PHOTONMOMENTUM
        if(d_bh_angle) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_bh_angle); d_bh_angle = NULL;}
        if(d_bh_lum)   {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_bh_lum);   d_bh_lum   = NULL;}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        if(d_cr_inject){Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_cr_inject); d_cr_inject = NULL;}
#endif
    };



#ifdef RT_USE_GRAVTREE
    const struct gpu_rt_walk_data_t rt_data_dev = rt_data_snap;
#endif
#ifdef SINK_PHOTONMOMENTUM
    const struct gpu_sink_walk_data_t sink_data_dev = sink_data_snap;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    const struct gpu_cr_walk_data_t cr_data_dev = cr_data_snap;
#endif

    /* Ewald periodic-image potential correction (pure-tree periodic + EVALPOTENTIAL).
     * Acquire the table once (idempotent). This build requires the correction, so a
     * missing table is a hard stop that aborts the primary walk -- never a silent
     * skip of the term. */
    struct gpu_ewald_pot_data_t ewald_pot_snap;
    ewald_pot_snap.potcorr = NULL; ewald_pot_snap.fac_intp = 0.0; ewald_pot_snap.active = 0;
#ifdef GIZMO_GPU_EWALD_POT_CORRECTION
    if(gpu_ewald_acquire_pot_data(&ewald_pot_snap) != 0) {
        printf("gpu_gravtree_walk_primary: Ewald potential-correction table unavailable; EVALPOTENTIAL periodic build requires it\n");
        endrun(913212);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
        release_payload_buffers();
        myfree(idx_host);   /* LIFO mymalloc cleanup; do not launch the walk with the term disabled */
        return 1;
    }
#endif
    const struct gpu_ewald_pot_data_t ewald_pot_dev = ewald_pot_snap;

#ifdef HERMITE_INTEGRATION
    struct gpu_hermite_walk_data_t hermite_dev;
    hermite_dev.state = HermiteWalk;
#endif

    /* Invariant guard: reset the per-walk counter. */
    g_inv_fterm_aggregate = 0; g_unship_aggregate = 0;

    /* The rest of the context: the source payload arrays and the Hermite and Ewald state.  The
     * pointers it holds are SharedSpace / captured-snapshot addresses valid for the launch. */
#ifdef RT_USE_GRAVTREE
    ctx.rt_data = rt_data_dev;
#endif
#ifdef SINK_PHOTONMOMENTUM
    ctx.sink_data = sink_data_dev;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    ctx.cr_data = cr_data_dev;
#endif
    ctx.use_lazy_source = use_lazy_source;
#ifdef HERMITE_INTEGRATION
    ctx.hermite = hermite_dev;
#endif
    ctx.ewald_pot = ewald_pot_dev;

    /* Choose the device schedule for this call. The engine takes the call only when this rank has
     * so few targets that giving each one a lane would leave most of the device idle; otherwise
     * this returns nonzero and the ordinary single-target walk below takes every candidate, which
     * is exactly what this code did before the cooperative walk existed.
     *
     * ⛔ Neither branch is a routing decision. The host/device question was settled before this
     * point and is not reopened here: both of these run on the device. */
    int walked_as_packets = (gpu_gravtree_walk_packets(ctx, d_idx, num_active, d_acc, d_ninter, d_pot, d_failed,
                                                       d_fail_by_reason) == 0);
    if(walked_as_packets) {
        /* The launcher has waited on the kernel, so the tally is complete and in shared space. */
        for(int r = 0; r < GRAV_PACKET_FAIL_REASONS; r++) {g_packet_fail[r] = (long long) d_fail_by_reason[r];}
    } else {
        /* Say so in the artifact rather than leaving the record empty: "the ordinary schedule ran"
           and "no device walk ran" are different facts and a pricing arm has to tell them apart. */
        g_packet_shape.mode = GRAV_PACKET_MODE_FLAT;
    }
    /* The single-target walk: every candidate when the engine did not run, otherwise only the
     * members of packets that failed. A packet fails as a whole when any member meets a
     * pseudo-particle, so its other members are walked here exactly as the flat route would
     * have walked them, and only a target that meets the pseudo-particle itself is left to the
     * host loop -- the same targets the flat route leaves it. */
    {
        const int only_failed = walked_as_packets;
        Kokkos::parallel_for("gravtree_walk_primary", num_active, KOKKOS_LAMBDA(int a) {
            if(only_failed && !d_failed[a]) {return;}
            int target = d_idx[a];
            Vec3<double> acc;
            int ninter;
            double pot;
            int ok = gpu_gravtree_walk_one(ctx, target, acc, ninter, pot);
            if(ok == 1) {
                d_acc[a] = acc;
                d_ninter[a] = ninter;
                d_pot[a] = pot;
                d_failed[a] = 0;
            } else {
                d_failed[a] = 1;   /* met a pseudo-particle */
            }
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("gravtree_walk_primary", num_active);
    }
    gpu_grav_report_time_fault("gpu_gravtree_walk_primary", time_fault, ti_curr_host, 913213);
    /* Import-completeness record.  A foreign node the sender shipped as a childless multipole, which
     * this walk's predicate now wants to descend, means the import no longer covers what the walk
     * asks of it.  Counted in the kernel and folded into the shared ledger here; gravity_tree()
     * reduces it across ranks and decides collectively whether to rebuild and redo, or stop. */
    gravity_note_incomplete_import_count(g_inv_fterm_aggregate);
    gravity_note_unshippable_import(g_unship_aggregate);

    /* Scatter successes back to host; copy RT CellP fields from device mirror */
    int nsucceeded = 0;
    double costtotal_added = 0;
    for(int a = 0; a < num_active; a++) {
        int i = d_idx[a];
        if(!d_failed[a]) {
            P[i].GravAccel = d_acc[a];
#ifdef EVALPOTENTIAL
            P[i].Potential = d_pot[a];
#endif

            /* RT scatter-back: copy outputs written into the SharedSpace
             * device mirrors back to host P[]/CellP[].  On UVM systems this
             * is effectively a same-pointer copy (no-op performance-wise),
             * but kept explicit for correctness on non-UVM targets. */
#ifdef RT_USE_TREECOL_FOR_NH
            {int k; for(k=0; k<RT_USE_TREECOL_FOR_NH; k++) {P[i].ColumnDensityBins[k] = P_dev[i].ColumnDensityBins[k];}}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
            P[i].MencInRcrit = P_dev[i].MencInRcrit;
#endif
#ifdef RT_USE_GRAVTREE
#ifdef RT_OTVET
            if(P[i].Type == 0) {
                int k; for(k=0; k<N_RT_FREQ_BINS; k++) {CellP[i].ET[k] = CellP_dev[i].ET[k];}
            }
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            if(P[i].Type == 0 && P[i].Mass > 0) {
                CellP[i].Rad_Flux_UV  = CellP_dev[i].Rad_Flux_UV;
                CellP[i].Rad_Flux_EUV = CellP_dev[i].Rad_Flux_EUV;
            }
#endif
#ifdef CHIMES_STELLAR_FLUXES
            if(P[i].Type == 0 && P[i].Mass > 0) {
                int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {
                    CellP[i].Chimes_G0[kc]          = CellP_dev[i].Chimes_G0[kc];
                    CellP[i].Chimes_fluxPhotIon[kc]  = CellP_dev[i].Chimes_fluxPhotIon[kc];
                }
            }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            if(P[i].Type == 0 && P[i].Mass > 0) {
                int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP[i].Rad_E_gamma[kf] = CellP_dev[i].Rad_E_gamma[kf];}
            }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
            if(P[i].Type == 0 && P[i].Mass > 0) {
                int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP[i].Rad_Flux[kf] = CellP_dev[i].Rad_Flux[kf];}
            }
#endif
#ifdef SINK_COMPTON_HEATING
            if(P[i].Type == 0 && P[i].Mass > 0) {
                CellP[i].Rad_Flux_AGN = CellP_dev[i].Rad_Flux_AGN;
            }
#endif
#endif /* RT_USE_GRAVTREE */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            if(P[i].Type == 0 && P[i].Mass > 0) {
                CellP[i].SubGrid_CosmicRayEnergyDensity = CellP_dev[i].SubGrid_CosmicRayEnergyDensity;
            }
#endif

            /* Tidal tensor + GravJerk scatter-back (ATFU/jerk path). */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
            P[i].tidal_tensorps = P_dev[i].tidal_tensorps;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
            P[i].GravJerk = P_dev[i].GravJerk;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
            /* Direct assignment: GPU walk visited the entire tree for this target
             * (no MPI export/import partition).  The post-loop += P[i].Mass at
             * gravtree.cc:605 then adds the target's own mass for the diagnostic. */
            P[i].TreeMass = P_dev[i].TreeMass;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
            P[i].tidal_zeta = P_dev[i].tidal_zeta;
#endif
#ifdef SPECIAL_POINT_MOTION
            P[i].vel_of_nearest_special = P_dev[i].vel_of_nearest_special;
            P[i].acc_of_nearest_special = P_dev[i].acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            P[i].weight_sum_for_special_point_smoothing = P_dev[i].weight_sum_for_special_point_smoothing;
#endif
#endif

            /* Sink-distance / single-star timestepping scatter-back */
#ifdef SINK_CALC_DISTANCES
            P[i].Min_Distance_to_Sink = P_dev[i].Min_Distance_to_Sink;
            P[i].Min_xyz_to_Sink      = P_dev[i].Min_xyz_to_Sink;
#ifdef SINGLE_STAR_FIND_BINARIES
            P[i].is_in_a_binary       = P_dev[i].is_in_a_binary;
            P[i].Min_Sink_OrbitalTime = P_dev[i].Min_Sink_OrbitalTime;
            if(P[i].is_in_a_binary) {
                P[i].comp_Mass = P_dev[i].comp_Mass;
                P[i].comp_dx   = P_dev[i].comp_dx;
                P[i].comp_dv   = P_dev[i].comp_dv;
            }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
            P[i].Min_Sink_Approach_Time = P_dev[i].Min_Sink_Approach_Time;
            P[i].Min_Sink_Freefall_time = P_dev[i].Min_Sink_Freefall_time;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            P[i].Min_Sink_FeedbackTime  = P_dev[i].Min_Sink_FeedbackTime;
#endif
#endif
#endif /* SINK_CALC_DISTANCES */

            ProcessedFlag[i] = 1;
            costtotal_added += d_ninter[a];
            if(TakeLevel >= 0) {P[i].GravCost[TakeLevel] = d_ninter[a];}

            /* No arena mirror-update here: under UVM-canonical
             * P_dev = arena_P aliases host P[], so the
             * struct copy P_dev[i] = P[i] would be self-assignment. */

            nsucceeded++;
        }
    }
    Costtotal += costtotal_added;

    /* mark_clean (not invalidate): the per-active-i
     * P_dev[i]=P[i] mirror in the scatter loop above keeps arena coherent
     * for the touched indices; untouched i's were unchanged from acquire-time
     * (kernel output went to d_* buffers, not arena). Verified safe with
     * GIZMO_GPU_ARENA_DEBUG=1. */
    gpu_particles_arena_mark_clean_after_scatter("gpu_gravtree_walk_primary");
    /* No SoA invalidate: the walk wrote nothing to the tree, and the next force_treebuild fully
     * repopulates the SoA via the build pipeline. */

    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
    release_payload_buffers();
    myfree(idx_host);

    if(host_candidates_left) {*host_candidates_left = num_active - nsucceeded;}
    return nsucceeded;
}



/* ========================================================================= *
 * GPU Ewald-correction walk (pure-tree periodic gravity).                   *
 *                                                                            *
 * Mirrors force_treeevaluate_ewald_correction() mode=0 (forcetree.cc:2631).  *
 * Runs as a second pass over active targets when Ewald_iter==1, after the    *
 * primary walk has completed its Ewald_iter==0 pass. Adds the periodic-image *
 * correction to P[target].GravAccel via trilinear interpolation of the Ewald *
 * lookup tables (fcorrx/y/z), with the same opening criteria and            *
 * nearest-image wrapping as the CPU walk. No optional payloads (monopole     *
 * only), no softening kernel.                                                *
 * ========================================================================= */
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)

/* SharedSpace mirrors of the four Ewald correction tables. Seeded once from
 * the CPU-side static tables in forcetree.cc via gizmo_get_ewald_tables().
 * Flat layout: index [i*(EN+1)^2 + j*(EN+1) + k] with EN = GIZMO_EWALD_EN. */
static MyFloat *g_d_fcorrx   = NULL;
static MyFloat *g_d_fcorry   = NULL;
static MyFloat *g_d_fcorrz   = NULL;
static MyFloat *g_d_potcorr  = NULL;
static double   g_ewald_fac_intp = 0.0;
static int      g_ewald_tables_ready = 0;

static int gpu_ewald_tables_acquire(void)
{
    if(g_ewald_tables_ready) return 0;
    const MyFloat *fx, *fy, *fz, *fp;
    double fi;
    gizmo_get_ewald_tables(&fx, &fy, &fz, &fp, &fi);
    long n  = (long)(GIZMO_EWALD_EN + 1) * (GIZMO_EWALD_EN + 1) * (GIZMO_EWALD_EN + 1);
    long sz = n * sizeof(MyFloat);
    g_d_fcorrx  = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    g_d_fcorry  = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    g_d_fcorrz  = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    g_d_potcorr = (MyFloat *) gizmo_gpu_alloc_shared(sz, "gravity_walk");
    if(!g_d_fcorrx || !g_d_fcorry || !g_d_fcorrz || !g_d_potcorr) {
        printf("gpu_ewald_tables_acquire: kokkos_malloc failed\n");
        endrun(914101);
        return 1;   /* soft bad-stop: caller skips the Ewald walk on not-ready tables */
    }
    memcpy(g_d_fcorrx,  fx, sz);
    memcpy(g_d_fcorry,  fy, sz);
    memcpy(g_d_fcorrz,  fz, sz);
    memcpy(g_d_potcorr, fp, sz);
    g_ewald_fac_intp = fi;
    g_ewald_tables_ready = 1;
    return 0;
}

#ifdef GIZMO_GPU_EWALD_POT_CORRECTION
static int gpu_ewald_acquire_pot_data(struct gpu_ewald_pot_data_t *out)
{
    out->potcorr = NULL; out->fac_intp = 0.0; out->active = 0;
    if(gpu_ewald_tables_acquire() != 0) return 1;   /* table acquire failed; caller hard-stops */
    out->potcorr  = g_d_potcorr;
    out->fac_intp = g_ewald_fac_intp;
    out->active   = 1;
    return 0;
}
#endif

/* Device-side Ewald walk for a single target. Returns 1 on success (acc
 * written), 0 if a pseudo-particle was encountered (defer to CPU). Every source
 * is read at the walk time without being drifted, through the same two readers
 * as the primary walk (position and mass only: the correction needs no more). */
static KOKKOS_INLINE_FUNCTION int
gpu_ewald_walk_one(int target,
                   int treeBase, int treeParticleSlots, int maxNodes, int maxForeignNodes,    /* LET */
                   struct particle_data *P_dev, struct gas_cell_data *CellP_dev, integertime ti,
                   const struct gpu_gravity_tree_soa_t *tree_soa, const struct extNODE *ext_nodes,
                   const struct DriftKickTableView &tables, struct gpu_grav_time_fault_t *time_fault,
#ifdef GRAVITY_HYBRID_OPENING_CRIT
                   int is_first_step,   /* hybrid opening: relative criterion applies only after step 0 */
#endif
                   const MyFloat *fcorrx, const MyFloat *fcorry, const MyFloat *fcorrz,
                   double fac_intp, double boxsize, double boxhalf,
                   double errtoltheta, double errtolforceacc,
                   Vec3<double> &acc_out)
{
    Vec3<double> pos = P_dev[target].Pos;
    double aold = errtolforceacc * P_dev[target].OldAcc;
    Vec3<double> acc = {0.0, 0.0, 0.0};

    int no = treeBase; /* root node */
    while(no >= 0)
    {
        double mass = 0.0, len_node = 0.0;
        Vec3<double> dr = {0.0, 0.0, 0.0};
        int is_leaf = 0;
        int idx = 0;
        /* What the imported tree allows this walk to do with the node; set once per node, above every
         * branch that could descend it, exactly as the primary walk does (see let_data.h).  Particle
         * leaves never reach the branches that read it, so a local node is the right default. */
        grav_node_kind_t node_kind = GRAV_NODE_LOCAL;
        if(no >= treeParticleSlots && no < treeBase) {return 0;} /* gap: malformed tree -- defer; the CPU Ewald walk's guard stops loudly */
        if(no < treeParticleSlots) /* particle leaf */
        {
            const gpu_grav_particle_now_t motion = gpu_grav_particle_source_at(P_dev, CellP_dev, no, ti, tables, time_fault);
            dr[0] = motion.pos[0] - pos[0];
            dr[1] = motion.pos[1] - pos[1];
            dr[2] = motion.pos[2] - pos[2];
            mass  = motion.mass;
            is_leaf = 1;
        }
        else if(no >= treeBase + maxNodes + maxForeignNodes) /* pseudo-particle — defer to CPU (foreign-node range below) */
        {
            return 0;
        }
        else /* internal node */
        {
            idx = no - treeBase;
            /* The wire tag is the authority on whether this node's children were shipped.  Read it
             * before either descent below.  Only the tag is needed: the Ewald correction depends on
             * mass and separation alone, so a foreign leaf carries no identity to restore. */
            {
                int in_foreign_n = (no >= treeBase + maxNodes);
                int fl_tag = 0;
                if(in_foreign_n && tree_soa->foreign_leaf_tag) {
                    int fs = idx - maxNodes;
                    if(fs >= 0 && fs < tree_soa->foreign_leaf_cap) {fl_tag = tree_soa->foreign_leaf_tag[fs];}
                }
                node_kind = grav_classify_node(in_foreign_n, fl_tag, tree_soa->nextnode[idx]);
            }
            /* skip single-particle node (open it to its daughter chain).  A terminal foreign node
             * holds one particle whose multipole is that particle exactly, and its nextnode is the
             * continuation past the subtree, so use it where it stands instead. */
            if(!(tree_soa->bitflags[idx] & (1 << BITFLAG_MULTIPLEPARTICLES))) {
                if(node_kind == GRAV_NODE_LOCAL || node_kind == GRAV_NODE_FOREIGN_OPENABLE) {
                    no = tree_soa->nextnode[idx];
                    continue;
                }
            }
            node_motion_prediction motion;
            if(gpu_grav_node_motion_at(*tree_soa, ext_nodes, idx, no, ti, tables, motion)) {
                gpu_grav_note_source_ahead(time_fault, 1, no, tree_soa->node_ti[idx]);
            }
            mass  = tree_soa->mass[idx];
            dr[0] = motion.s_[0] - pos[0];
            dr[1] = motion.s_[1] - pos[1];
            dr[2] = motion.s_[2] - pos[2];
            len_node = motion.len_;   /* the geometric centre does not move between builds; the length is widened */
        }

        /* nearest-image wrap on the displacement (shared SSOT helper) */
        gravity_box_nearest_image(dr[0], dr[1], dr[2], -1);

        if(is_leaf) {
            no = tree_soa->nextnode_aux[no];
        } else {
            /* Opening check + periodic-boundary skip (mirrors forcetree.cc:2769-2842) */
            double r2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if(r2 <= 0) r2 = 1e-300;
            double len = len_node;
            int openflag = 0;
            if(errtoltheta) {
                if(len * len > r2 * errtoltheta * errtoltheta) openflag = 1;
            }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
            else {
#else
            /* hybrid: relative criterion only after step 0 (mirrors forcetree.cc:3489-3493) */
            if(!is_first_step) {
#endif
                if(mass * len * len > r2 * r2 * aold) {
                    openflag = 1;
                } else {
                    double ad0 = tree_soa->center[idx][0] - pos[0], ad1 = tree_soa->center[idx][1] - pos[1], ad2 = tree_soa->center[idx][2] - pos[2];
                    double adx = gravity_box_long_abs_x(ad0, ad1, ad2, -1);
                    double ady = gravity_box_long_abs_y(ad0, ad1, ad2, -1);
                    double adz = gravity_box_long_abs_z(ad0, ad1, ad2, -1);
                    if(adx < 0.60*len && ady < 0.60*len && adz < 0.60*len) openflag = 1;
                }
            }
            if(openflag) {
                /* The multipole may still stand in for the whole node, but only if the node lies
                 * wholly on one side of every periodic boundary -- otherwise its daughters map to
                 * different nearest images -- and only if it is small compared with the box. */
                int must_refine = 0;
                for(int kdim = 0; kdim < 3 && !must_refine; kdim++) {
                    double uk = tree_soa->center[idx][kdim] - pos[kdim];
                    if(uk >  boxhalf) {uk -= boxsize;} else if(uk < -boxhalf) {uk += boxsize;}
                    if(((uk < 0) ? -uk : uk) > 0.5*(boxsize - len)) { must_refine = 1; }
                }
                if(!must_refine && len > 0.20 * boxsize) { must_refine = 1; }  /* cell too large */
                if(must_refine) {
                    /* Only a node shipped with its children can be resolved below.  A foreign leaf is
                     * a single particle, so refining it would change nothing.  A truncated aggregate
                     * is neither: counted on device, reported once by the host after the walk. */
                    if(!grav_node_is_terminal(node_kind)) { no = tree_soa->nextnode[idx]; continue; }
                    if(node_kind == GRAV_NODE_FOREIGN_TRUNCATED || node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE) {
                        Kokkos::atomic_add(&g_inv_fterm_aggregate, 1LL);
                        if(node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE) {Kokkos::atomic_add(&g_unship_aggregate, 1LL);}
                    }
                }
            }
            no = tree_soa->sibling[idx];
        }

        /* Trilinear interp of the Ewald force octant tables via the shared SSOT helper
         * (gravtree_ewald.h): weights from |dr| once, applied to all three tables; the
         * odd-force per-component signs stay here. */
        double signx = (dr[0] < 0) ? +1.0 : -1.0;
        double signy = (dr[1] < 0) ? +1.0 : -1.0;
        double signz = (dr[2] < 0) ? +1.0 : -1.0;
        grav_ewald_interp_weights ew = grav_ewald_interp_setup(dr[0], dr[1], dr[2], fac_intp);
        acc[0] += mass * signx * grav_ewald_interp_apply(fcorrx, ew);
        acc[1] += mass * signy * grav_ewald_interp_apply(fcorry, ew);
        acc[2] += mass * signz * grav_ewald_interp_apply(fcorrz, ew);
    }

    acc_out = acc;
    return 1;
}

extern "C" int gpu_ewald_walk_primary(void)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();
    if(Ewald_iter == 0) {return 0;}
#ifdef PMGRID
    return 0; /* Ewald walk not needed under TreePM (gravtree.cc:734 gate) */
#else
    int num_active_total = (int) ActiveParticleList.size();
    if(num_active_total <= 0) {return 0;}

    const integertime ti_curr_host = gizmo_host_ti_current();

    /* Re-acquire the tree SoA (cache-hit).  soa->nextnode_aux aliases UVM Nextnode[]; no per-walk
     * memcpy.  Sources are read at the walk time, as in the primary walk, so this pass runs whichever
     * route the primary took. */
    int min_nodes = MaxNodes + 1;
    gpu_gravity_tree_acquire(min_nodes, Nodes_base, Extnodes_base);
    const struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa) {return 0;}

    gpu_particles_arena_set_site("gpu_gravtree_walk_ewald");
    gpu_particles_arena_acquire(NumPart, P, CellP);
    struct particle_data *P_dev = gpu_particles_arena_P();
    struct gas_cell_data *CellP_dev = gpu_particles_arena_CellP();   /* NULL on a gas-free problem; read only for a gas source */
    if(!P_dev || ((All.TotN_gas > 0) && !CellP_dev)) {return 0;}

    if(gpu_ewald_tables_acquire() != 0) {return 1;}   /* soft bad-stop: tables not ready, skip walk (idx_host not yet alloc'd); drains at next poll */

    int *idx_host = (int *) mymalloc("gpu_ewald_idx", num_active_total * sizeof(int));
    int num_active = 0;
    for(int a = 0; a < num_active_total; a++) {
        int i = ActiveParticleList[a];
        if(ProcessedFlag[i]) continue;
        /* SSOT pre-walk candidacy (Mass>0 + Hermite eligibility + needs_new_treeforce):
         * parity with the CPU primary loop + finalization (see primary walk). */
        if(!gravity_treewalk_candidate_prewalk(i, a)) continue;
        idx_host[num_active++] = i;
    }
    if(num_active == 0) { myfree(idx_host); return 0; }

    if(gpu_grav_answer_node_claims(ti_curr_host, "gpu_ewald_walk_primary") != 0) {myfree(idx_host); return 1;}
    struct DriftKickTableView tables_snap;
    if(drift_kick_table_mirror_refresh(&tu_drift_kick_table_dev, &tables_snap) != 0) {myfree(idx_host); return 1;}   /* soft bad-stop already requested */

    const struct grav_walk_scratch_plan scratch = grav_walk_scratch_plan_for(num_active, 0);
    char *scratch_block = (char *) gizmo_gpu_alloc_shared(scratch.bytes, "gravity_walk");
    if(!scratch_block) {printf("gpu_ewald_walk_primary: kokkos_malloc failed\n"); endrun(914102); myfree(idx_host); return 1;}
    struct gpu_grav_time_fault_t *time_fault = (struct gpu_grav_time_fault_t *) (scratch_block + scratch.time_fault);
    memset(time_fault, 0, sizeof(struct gpu_grav_time_fault_t));
    const struct extNODE *Extnodes_snap = Extnodes;
    int          *d_idx    = (int *)          (scratch_block + scratch.idx);
    int          *d_failed = (int *)          (scratch_block + scratch.failed);
    Vec3<double> *d_acc    = (Vec3<double> *) (scratch_block + scratch.acc);
    memcpy(d_idx, idx_host, num_active * sizeof(int));
    memset(d_failed, 0, num_active * sizeof(int));

    /* snapshot scalars */
    const int treeBase            = All.TreeNodeIndexBase;
    const int treeParticleSlots_snap = All.TreeParticleSlots;
    const int maxNodes_snap      = MaxNodes;
    const int maxForeignNodes_sn = MaxForeignNodes;    /* LET */
    const double boxsize     = All.BoxSize;
    const double boxhalf     = 0.5 * All.BoxSize;
    const double fac_intp    = g_ewald_fac_intp;
    const double errtoltheta = All.ErrTolTheta;
    const double errtolforceacc = All.ErrTolForceAcc;
    const MyFloat *fcorrx_dev = g_d_fcorrx;
    const MyFloat *fcorry_dev = g_d_fcorry;
    const MyFloat *fcorrz_dev = g_d_fcorrz;
    const struct gpu_gravity_tree_soa_t soa_snap = *soa;
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    /* host-evaluate the first-step predicate once; captured by value into the device walk */
    const int is_first_step_snap = (All.Ti_Current == 0 && RestartFlag != 1);
#endif

    /* Import-completeness guard: reset the per-walk counter (see the primary walk). */
    g_inv_fterm_aggregate = 0; g_unship_aggregate = 0;

    Kokkos::parallel_for("gpu_ewald_walk_primary", num_active, KOKKOS_LAMBDA(int a) {
        int target = d_idx[a];
        Vec3<double> acc;
        int ok = gpu_ewald_walk_one(target, treeBase, treeParticleSlots_snap, maxNodes_snap, maxForeignNodes_sn,
                                     P_dev, CellP_dev, ti_curr_host, &soa_snap, Extnodes_snap, tables_snap, time_fault,
#ifdef GRAVITY_HYBRID_OPENING_CRIT
                                     is_first_step_snap,
#endif
                                     fcorrx_dev, fcorry_dev, fcorrz_dev,
                                     fac_intp, boxsize, boxhalf,
                                     errtoltheta, errtolforceacc,
                                     acc);
        if(ok == 1) {d_acc[a] = acc; d_failed[a] = 0;}
        else        {d_failed[a] = 1;}   /* met a pseudo-particle */
    });
    Kokkos::fence();
    gizmo_gpu_check_last_error("gpu_ewald_walk_primary", num_active);
    gpu_grav_report_time_fault("gpu_ewald_walk_primary", time_fault, ti_curr_host, 914103);
    /* Import-completeness record, same contract and same ledger as the primary walk: a foreign node
     * the sender shipped as a childless multipole, which this walk has to resolve below, means the
     * import no longer covers what the walk asks of it. */
    gravity_note_incomplete_import_count(g_inv_fterm_aggregate);
    gravity_note_unshippable_import(g_unship_aggregate);

    int nsucceeded = 0;
    for(int a = 0; a < num_active; a++) {
        if(!d_failed[a]) {
            int target = idx_host[a];
            P[target].GravAccel[0] += d_acc[a][0];
            P[target].GravAccel[1] += d_acc[a][1];
            P[target].GravAccel[2] += d_acc[a][2];
            ProcessedFlag[target] = 1;
            nsucceeded++;
        }
    }

    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
    myfree(idx_host);
    return nsucceeded;
#endif /* !PMGRID */
}

#else /* !(BOX_PERIODIC && !GRAVITY_NOT_PERIODIC) */

extern "C" int gpu_ewald_walk_primary(void) {return 0;}

#endif


