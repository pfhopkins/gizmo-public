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
#include "gpu_gravity_tree.h"
#include "gpu_gravtree.h"
#include "forcetree.h"
#include "gravity_box_distance.h"   /* shared CPU/GPU gravity box-distance SSOT */
#include "gravtree_opening.h"       /* shared CPU/GPU primary-walk acceptance-geometry predicate (SSOT) */
#include "let_data.h"             /* LET_LEAF_TAG_* + grav_classify_node (import topology vocabulary) */

#include "../mesh/kernel.h"
#include "gravtree_force_kernel.h"  /* shared CPU/GPU accepted-source contribution physics (SSOT) */
#include "gravtree_ewald.h"         /* shared CPU/GPU Ewald image-correction trilinear interp (SSOT) */
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
 * passes them by value, together with the drift/kick table mirror the prediction needs.
 * No device allocation and no per-particle array -- this is the few-body regime, where a
 * pass over the particles would cost more than the walk it serves. */
#ifdef HERMITE_INTEGRATION
struct gpu_hermite_walk_data_t {
    struct HermiteWalkState   state;
    struct DriftKickTableView tables;
};
#endif

#ifdef HERMITE_INTEGRATION
/* SharedSpace mirror of the drift/kick tables, allocated on first use and reused. Stays NULL
   on a non-cosmological run, which every Hermite few-body problem is: the refresh below
   returns an elapsed-time view that reads no table at all. */
static double *hermite_drift_kick_table_dev = NULL;

extern "C" void gpu_gravtree_hermite_release(void)
{
    if(hermite_drift_kick_table_dev) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(hermite_drift_kick_table_dev);
        hermite_drift_kick_table_dev = NULL;
    }
}
#endif

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
 * StarParticleEffectiveSize) handled inline in gpu_force_softening_kernel_radius
 * below. */
/* ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT (type-0 softening cap) handled inline
 * in gpu_force_softening_kernel_radius below. */
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
    struct gpu_gravity_tree_soa_t tree_soa;   /* the SoA handle set, by value (a struct of pointers into SharedSpace) */
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

struct gpu_grav_member_t {
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
    double r_for_total_menclosed, m_enc_in_rcrit;
#endif
    /* accumulators */
    grav_pair_acc_t out;
#ifdef COUNT_MASS_IN_GRAVTREE
    double tree_mass;
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

/* Set up one target's member for the walk (the CPU walk's target prologue).
 * Returns 0 for a massless target, which takes part in nothing and writes zeros. */
static KOKKOS_INLINE_FUNCTION int
gpu_grav_member_init(const gpu_grav_walk_ctx_t &ctx, int target, gpu_grav_member_t &mem)
{
    struct particle_data *P_dev = ctx.P_dev;
    mem.target = target;
    mem.open.pos = P_dev[target].Pos;
    mem.open.ptype = P_dev[target].Type;
    mem.open.alive = 0;
    mem.pmass = P_dev[target].Mass;
    grav_pair_acc_init(mem.out);
    if(mem.pmass <= 0) {return 0;}
    mem.open.alive = 1;
    const int ptype = mem.open.ptype; const double pmass = mem.pmass;

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(GALSF_MERGER_STARCLUSTER_PARTICLES)
    double soft = gpu_force_softening_kernel_radius(P_dev, target);
#else
    double soft = All.ForceSoftening[ptype];
#endif
    double zeta = 0.0;    /* unconditional (matches CPU walk); passed to the shared pair kernel, consumed there only under #if AGS */
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    grav_target_select_soft_and_zeta(ptype, gpu_get_ags_zeta(P_dev, target), soft, zeta);
#endif
    mem.open.soft = soft; mem.zeta = zeta;
    mem.open.aold = All.ErrTolForceAcc * P_dev[target].OldAcc;

    grav_pm_shortrange_t pm = ctx.pm;
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
    /* high-res zoom particles use the finer short-range PM cutoff (mirrors forcetree.cc target
     * prologue). The dispatcher passes the coarse-mesh rcut/asmthfac; override per target here. */
    if(pmforce_is_particle_high_res(ptype, mem.open.pos)) {
        pm.rcut = All.Rcut[1]; pm.rcut2 = pm.rcut * pm.rcut; pm.asmthfac = grav_pm_asmthfac(All.Asmth[1]);
    }
#endif
#ifdef PMGRID
    mem.open.rcut = pm.rcut; mem.open.rcut2 = pm.rcut2;
#endif

    /* fed unconditionally to the shared pair kernel (consumed there only under the
     * symmetrize-by-averaging #if); matches the CPU walk's unconditional precompute. */
    const int ags_bitflag_primary = gravtree_ags_kernel_shared_bitflag(ptype);

#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    mem.vel = P_dev[target].Vel;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    /* Diagnostic: total mass seen by this target during the walk, summed only
     * over accepted interactions (mirrors forcetree.cc). The walk excludes the
     * target's own leaf (r2==0); the post-loop +=P[i].Mass in gravtree.cc
     * finalizes the sum. */
    mem.tree_mass = 0.0;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    /* Shell-theorem gravity: forces from any source at r_source > r_target
     * vanish; forces from r_source < r_target use a 1/r^3 enclosed-mass formula
     * pointed toward the box center. Mirrors forcetree.cc. */
    mem.sph_center[0] = 0.0; mem.sph_center[1] = 0.0; mem.sph_center[2] = 0.0;
#ifdef BOX_PERIODIC
    mem.sph_center[0] = 0.5 * boxSize_X;
    mem.sph_center[1] = 0.5 * boxSize_Y;
    mem.sph_center[2] = 0.5 * boxSize_Z;
#endif
#endif
#ifdef SINK_COMPTON_HEATING
    mem.incident_flux_agn = 0.0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    mem.SubGrid_CosmicRayEnergyDensity = 0.0;
    /* per-target CR gate (the host precompute leaves t_max_cr=0 unless All.Time>All.TimeBegin,
     * mirroring the CPU walk's gate) */
    mem.cr_active_gate = (ctx.cr_data.t_max_cr > 0) ? 1 : 0;
#endif
#ifdef SINK_CALC_DISTANCES
    grav_sink_prox_accum_init(mem.sink_prox);
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    {int kb; for(kb=0; kb<RT_USE_TREECOL_FOR_NH; kb++) {mem.treecol_angular_bins[kb]=0.0;}}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    mem.m_enc_in_rcrit = 0.0; mem.r_for_total_menclosed = grav_target_menc_radius(soft); /* baseline Rcrit_min applied in the helper */
#endif
#ifdef RT_USE_GRAVTREE
#ifdef CHIMES_STELLAR_FLUXES
    {int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {mem.chimes_flux_G0[kc]=0; mem.chimes_flux_ion[kc]=0;}}
#endif
    /* valid-gas RT gate via the shared helper */
    mem.valid_gas_particle_for_rt = grav_target_valid_gas_for_rt(ptype, soft, pmass);
#ifdef RT_OTVET
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.RT_ET[kf] = {};}}
#endif
#endif /* RT_USE_GRAVTREE */
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    mem.incident_flux_uv = 0.0; mem.incident_flux_euv = 0.0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.Rad_E_gamma[kf]=0.0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.Rad_Flux[kf]={};}}
#endif

    /* RT_LEBRON radiation-pressure coupling factor (forcetree.cc). Once per target, before
     * the walk, because it only depends on the target's properties. Skipped when save-flux
     * mode is active (flux is stored and converted to RP after the walk by the caller). */
#if defined(RT_USE_GRAVTREE) && defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.fac_stellum[kf]=0.0;}}
    {
        volatile int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt) {
            double kappa_eff[N_RT_FREQ_BINS]; int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {kappa_eff[kf] = rt_kappa(-1, kf, P_dev, ctx.CellP_dev);}
            grav_target_rt_fac_stellum(soft, pmass, kappa_eff, mem.fac_stellum);
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
    tgt.pos = mem.open.pos; tgt.center[0] = mem.sph_center[0]; tgt.center[1] = mem.sph_center[1]; tgt.center[2] = mem.sph_center[2];
#endif
    mem.tgt = tgt;
    return 1;
}

/* The member-independent view of a tree node: its moments and geometry as stored,
 * its import classification, and what the walk does with it before any member
 * decides. */
enum gpu_grav_node_step_t { GPU_GRAV_NODE_SKIP_TO_SIBLING, GPU_GRAV_NODE_DESCEND, GPU_GRAV_NODE_DECIDE };
struct gpu_grav_node_prelude_t {
    int idx;                        /* SoA index no - treeBase */
    Vec3<MyFloat> s_node;           /* centre of mass as stored (before any star-target subtraction) */
    MyFloat len_node, msoft_node, mass_node;
    Vec3<MyFloat> center_node;
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
    nd.s_node = Vec3<MyFloat>{(MyFloat)tree_soa->s[idx][0], (MyFloat)tree_soa->s[idx][1], (MyFloat)tree_soa->s[idx][2]};
    nd.len_node = tree_soa->len[idx];
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
        Vec3<MyFloat> sp = Vec3<MyFloat>{(MyFloat)ctx.tree_soa.sink_pos[nd.idx][0],
                                         (MyFloat)ctx.tree_soa.sink_pos[nd.idx][1],
                                         (MyFloat)ctx.tree_soa.sink_pos[nd.idx][2]};
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
gpu_grav_evaluate_pair(const gpu_grav_walk_ctx_t &ctx, gpu_grav_member_t &mem, grav_pair_src_t &src, const gpu_grav_src_payload_t &pl)
{
    if(!((src.r2 > 0.0) && (src.mass > 0.0))) {return;}
    /* pair-wise gravity, PM truncation, and the accumulations inside the PM short-range
     * gate, via the shared evaluation (gravtree_force_kernel.h), the single home for the
     * pair physics on both walks */
    grav_pair_result_t res = grav_pair_evaluate_core(mem.tgt, src, mem.out);
    const double r = res.r, fac_accel = res.fac_accel;
    (void) r; (void) fac_accel;
#ifdef GIZMO_GPU_EWALD_POT_CORRECTION
    /* Ewald periodic-image potential correction (mirrors forcetree.cc), from dr as the
     * evaluation left it. Pure-tree periodic only; under PMGRID the long-range potential
     * comes from the PM solver.  active is 1 in a healthy run (acquire failure hard-stops
     * the caller); the guard only covers the post-endrun drain. */
    if(ctx.ewald_pot.active) {
        grav_ewald_interp_weights ew = grav_ewald_interp_setup(src.dr[0], src.dr[1], src.dr[2], ctx.ewald_pot.fac_intp);
        mem.out.pot += src.mass * grav_ewald_interp_apply(ctx.ewald_pot.potcorr, ew);
    }
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    /* counted only for accepted interactions (r2>0, mass>0), mirroring forcetree.cc -- the
     * walk excludes the target's own (r2==0) leaf; gravtree.cc adds it back exactly once. */
    mem.tree_mass += src.mass;
#endif

    /* RT cluster payloads.  Structure mirrors forcetree.cc: OUTSIDE the PM short-range
     * gate, so for an out-of-range source fac_accel is the raw un-truncated value here
     * (used by the TREECOL column estimate, which has no PM-side completion);
     * RT_USE_GRAVTREE computes its own fac_rt from d_stellarlum independently. */
#ifdef RT_USE_TREECOL_FOR_NH
    {
        const double angular_bin_size = 4.0 * M_PI / RT_USE_TREECOL_FOR_NH;
        grav_treecol_accumulate(src.dr, r, fac_accel, pl.gasmass, src.mass, angular_bin_size, mem.treecol_angular_bins);
    }
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    /* per-interaction mass accumulation: each visited node contributes its multipole mass when within Rcrit */
    if(r < mem.r_for_total_menclosed) {mem.m_enc_in_rcrit += src.mass;}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    grav_cr_lebron_accumulate(mem.open.ptype, r, mem.open.soft, pl.cr_injection, mem.cr_active_gate, ctx.cr_data.t_max_cr, mem.tgt.pm, mem.SubGrid_CosmicRayEnergyDensity);
#endif
#ifdef RT_USE_GRAVTREE
    {
        volatile int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
        if(valid_gas_particle_for_rt)
        {
            /* payload formulas in the shared helper; fac_rt computed there from d_stellarlum
             * (may differ from dr when RT_SEPARATELY_TRACK_LUMPOS; otherwise d_stellarlum == dr) */
            grav_rt_src_t rt_src = {}; rt_src.d_stellarlum = pl.d_stellarlum; rt_src.soft = mem.open.soft; rt_src.mass_stellarlum = pl.mass_stellarlum;
#ifdef CHIMES_STELLAR_FLUXES
            rt_src.chimes_mass_stellarlum_G0 = pl.chimes_mass_stellarlum_G0; rt_src.chimes_mass_stellarlum_ion = pl.chimes_mass_stellarlum_ion;
#endif
#ifdef SINK_PHOTONMOMENTUM
            rt_src.mass_sinklumwt_forradfb = pl.mass_sinklumwt_forradfb;
#endif
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
            rt_src.fac_stellum = mem.fac_stellum;
#endif
            grav_rt_accum_t rt_accum = {};
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
            rt_accum.Rad_E_gamma = mem.Rad_E_gamma;
#endif
#ifdef CHIMES_STELLAR_FLUXES
            rt_accum.chimes_flux_G0 = mem.chimes_flux_G0; rt_accum.chimes_flux_ion = mem.chimes_flux_ion;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            rt_accum.incident_flux_uv = &mem.incident_flux_uv; rt_accum.incident_flux_euv = &mem.incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
            rt_accum.incident_flux_agn = &mem.incident_flux_agn;
#endif
#ifdef RT_OTVET
            rt_accum.RT_ET = mem.RT_ET;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
            rt_accum.Rad_Flux = mem.Rad_Flux;
#endif
            grav_rt_payload_accumulate(rt_src, rt_accum, mem.out.acc);
        }
    }
#endif /* RT_USE_GRAVTREE */
#ifdef DM_SCALARFIELD_SCREENING
    /* Yukawa-screened scalar-field force on non-gas targets (shared helper;
     * own table gate keyed on the dm-center distance, outside the main PM gate) */
    if(mem.open.ptype != 0)
    {
        Vec3<double> d_dm = pl.d_dm;   /* the helper takes the displacement by non-const reference */
        grav_dm_scalarfield_accumulate(d_dm, pl.mass_dm_local, mem.open.soft, mem.tgt.pm, mem.out.acc);
    }
#endif
}

/* Load a particle leaf for a member through the P_dev adapter and evaluate it. The
 * source state is the drifted state, except where the Hermite predictor replaces it. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_evaluate_leaf(const gpu_grav_walk_ctx_t &ctx, int no, gpu_grav_member_t &mem)
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

    Vec3<double> src_pos = P_dev[no].Pos;
    Vec3<double> src_vel = P_dev[no].Vel;   /* unconditional, mirroring forcetree.cc: the sink-proximity block reads it under SINK_CALC_DISTANCES, which several flags reach without the jerk or dynamical-friction terms */
#ifdef HERMITE_INTEGRATION
    /* On a Hermite pass a source the Hermite integrator owns but is not advancing this
       step is second-order wrong where it stands; evaluate it from its own start-of-step
       state instead. Same helper and same conditions as the host walk, so the two agree
       whichever one a step routes to. Single sources only; nothing is written back. */
    if(hermite_source_needs_prediction(no, P_dev, ctx.hermite.state)) {
        hermite_predict_source_state(no, P_dev, ctx.hermite.state, &ctx.hermite.tables, src_pos, src_vel);
    }
#endif
    src.dr = src_pos - mem.open.pos;
    gravity_box_nearest_image(src.dr[0], src.dr[1], src.dr[2], -1);
    src.r2 = src.dr.norm_sq();
    src.mass = P_dev[no].Mass;
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
    if(mem.open.ptype != 0 && P_dev[no].Type == 1) { pl.d_dm = src.dr; pl.mass_dm_local = src.mass; }
    else { pl.d_dm = Vec3<double>{0,0,0}; pl.mass_dm_local = 0; }
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = grav_spherical_symmetry_r_from_center(src_pos[0],src_pos[1],src_pos[2],mem.sph_center[0],mem.sph_center[1],mem.sph_center[2]);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    /* the secondary's previous-step tidal tensor (mirrors forcetree.cc) */
    {
        SymmetricTensor2<MyFloat> tmp = P_dev[no].tidal_tensorps_prevstep;
        for(int kk = 0; kk < 6; kk++) src.j_zeta_tidal_tensorps_prevstep.data[kk] = (double) tmp.data[kk];
    }
#endif
#if defined(SINK_DYNFRICTION_FROMTREE) || defined(COMPUTE_JERK_IN_GRAVTREE)
    src.dv = src_vel - mem.vel;
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
    pl.gasmass = (P_dev[no].Type == 0) ? P_dev[no].Mass : 0.0;
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
        volatile int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
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
        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = mem.open.ptype; prox_target.pmass = mem.pmass; prox_target.soft = mem.open.soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
        prox_target.vel = mem.vel;
#endif
        grav_sink_prox_leaf_src_t prox_src = {}; prox_src.src_type = P_dev[no].Type; prox_src.src_mass = P_dev[no].Mass; prox_src.motion.vel = src_vel;   /* the state this interaction was evaluated at, mirroring forcetree.cc, so (dr, vel) stays a consistent pair on a Hermite pass */
#if defined(SPECIAL_POINT_MOTION) || defined(SPECIAL_POINT_WEIGHTED_MOTION)
        prox_src.motion.acc = P_dev[no].Acc_Total_PrevStep;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        prox_src.motion.max_feedback_vel = P_dev[no].MaxFeedbackVel;
#endif
        grav_sink_prox_leaf_accumulate(src.r2, src.dr, prox_target, prox_src, mem.sink_prox);
    }
#endif /* SINK_CALC_DISTANCES */

    gpu_grav_evaluate_pair(ctx, mem, src, pl);
}

/* Load an accepted node for a member through the SoA adapter and evaluate it, given
 * the node's geometry as this member derived it (gpu_grav_node_member_geometry). A
 * tagged real foreign single-particle leaf is consumed with particle-leaf secondary
 * semantics: the node payload supplies mass/h_p and the synthesized RT/sink/CR/tidal
 * terms (singleton-aggregate == particle value); the two leaf-identity fields the
 * moment cannot carry (Type + AGS_zeta) are restored via the shared seam so
 * grav_force_pair applies AGS symmetrization/zeta exactly as on the source's home rank. */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_evaluate_node(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_node_prelude_t &nd, gpu_grav_member_t &mem,
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
    if(mem.open.ptype != 0) {
        pl.d_dm[0] = (double)tree_soa->s_dm[idx][0] - mem.open.pos[0];
        pl.d_dm[1] = (double)tree_soa->s_dm[idx][1] - mem.open.pos[1];
        pl.d_dm[2] = (double)tree_soa->s_dm[idx][2] - mem.open.pos[2];
        pl.mass_dm_local = (double)tree_soa->mass_dm[idx];
    } else { pl.d_dm = Vec3<double>{0,0,0}; pl.mass_dm_local = 0; }
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    src.r_source = grav_spherical_symmetry_r_from_center(s_node[0],s_node[1],s_node[2],mem.sph_center[0],mem.sph_center[1],mem.sph_center[2]);
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
    src.dv[0] = (double) tree_soa->node_vs[idx][0] - mem.vel[0];
    src.dv[1] = (double) tree_soa->node_vs[idx][1] - mem.vel[1];
    src.dv[2] = (double) tree_soa->node_vs[idx][2] - mem.vel[2];
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
        volatile int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
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
            pl.d_stellarlum[0] = tree_soa->rt_source_lum_s[idx][0] - mem.open.pos[0];
            pl.d_stellarlum[1] = tree_soa->rt_source_lum_s[idx][1] - mem.open.pos[1];
            pl.d_stellarlum[2] = tree_soa->rt_source_lum_s[idx][2] - mem.open.pos[2];
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
        Vec3<double> node_vs = Vec3<double>{(double)tree_soa->node_vs[idx][0], (double)tree_soa->node_vs[idx][1], (double)tree_soa->node_vs[idx][2]};
        grav_sink_prox_node_specialweighted(src.r2, node_vs, mem.open.ptype, mem.sink_prox);
    }
#endif
    if(tree_soa->sink_mass[idx] > 0)
    {
        Vec3<double> sink_dr;
        sink_dr[0] = tree_soa->sink_pos[idx][0] - mem.open.pos[0];
        sink_dr[1] = tree_soa->sink_pos[idx][1] - mem.open.pos[1];
        sink_dr[2] = tree_soa->sink_pos[idx][2] - mem.open.pos[2];
        gravity_box_nearest_image(sink_dr[0], sink_dr[1], sink_dr[2], -1);
        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = mem.open.ptype; prox_target.pmass = mem.pmass; prox_target.soft = mem.open.soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
        prox_target.vel = mem.vel;
#endif
        grav_sink_prox_node_src_t prox_src = {}; prox_src.sink_mass = (double) tree_soa->sink_mass[idx];
#if defined(SINGLE_STAR_FIND_BINARIES)
        prox_src.n_sink = (int) tree_soa->N_SINK[idx];
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SPECIAL_POINT_MOTION)
        prox_src.motion.vel = tree_soa->sink_vel[idx];
#endif
#if defined(SPECIAL_POINT_MOTION)
        prox_src.motion.acc = tree_soa->sink_acc[idx];
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        prox_src.motion.max_feedback_vel = tree_soa->MaxFeedbackVel[idx];
#endif
        grav_sink_prox_node_accumulate(src.r2, sink_dr, prox_src, prox_target, mem.sink_prox);
    }
#endif /* SINK_CALC_DISTANCES */

    gpu_grav_evaluate_pair(ctx, mem, src, pl);
}

/* Write a completed member's outputs to P_dev / CellP_dev (the host scatter loop in
 * gpu_gravtree_walk_primary copies them to P[] / CellP[]) and return the three the
 * caller collects directly. Mirrors forcetree.cc (mode=0). */
static KOKKOS_INLINE_FUNCTION void
gpu_grav_member_finish(const gpu_grav_walk_ctx_t &ctx, const gpu_grav_member_t &mem, Vec3<double> &acc_out, int &ninter_out, double &pot_out)
{
    struct particle_data *P_dev = ctx.P_dev; const int target = mem.target;
#ifdef RT_USE_GRAVTREE
    struct gas_cell_data *CellP_dev = ctx.CellP_dev;
    volatile int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;   /* nvc++ miscompiles raw boolean gates in device code */
#endif
#ifdef RT_USE_TREECOL_FOR_NH
    {int k; for(k=0; k<RT_USE_TREECOL_FOR_NH; k++) {P_dev[target].ColumnDensityBins[k] = mem.treecol_angular_bins[k];}}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    P_dev[target].MencInRcrit = mem.m_enc_in_rcrit;
#endif
#ifdef RT_USE_GRAVTREE
#ifdef RT_OTVET
    if(valid_gas_particle_for_rt) {
        int k; for(k=0; k<N_RT_FREQ_BINS; k++) {CellP_dev[target].ET[k] = mem.RT_ET[k];}
    } else if(mem.open.ptype == 0) {
        int k; for(k=0; k<N_RT_FREQ_BINS; k++) {CellP_dev[target].ET[k] = {};}
    }
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    if(valid_gas_particle_for_rt) {
        CellP_dev[target].Rad_Flux_UV  = mem.incident_flux_uv;
        CellP_dev[target].Rad_Flux_EUV = mem.incident_flux_euv;
    }
#endif
#ifdef CHIMES_STELLAR_FLUXES
    if(valid_gas_particle_for_rt) {
        int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {
            CellP_dev[target].Chimes_G0[kc]          = mem.chimes_flux_G0[kc];
            CellP_dev[target].Chimes_fluxPhotIon[kc] = mem.chimes_flux_ion[kc];
        }
    }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    if(valid_gas_particle_for_rt) {
        int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP_dev[target].Rad_E_gamma[kf] = mem.Rad_E_gamma[kf];}
    }
#endif
#ifdef SINK_COMPTON_HEATING
    if(valid_gas_particle_for_rt) {
        CellP_dev[target].Rad_Flux_AGN = mem.incident_flux_agn;
    }
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    if(valid_gas_particle_for_rt) {
        int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {CellP_dev[target].Rad_Flux[kf] = mem.Rad_Flux[kf];}
    }
#endif
#endif /* RT_USE_GRAVTREE */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    if(mem.open.ptype == 0) {ctx.CellP_dev[target].SubGrid_CosmicRayEnergyDensity = mem.SubGrid_CosmicRayEnergyDensity;}
#endif
#ifdef SINK_CALC_DISTANCES
    P_dev[target].Min_Distance_to_Sink = sqrt(mem.sink_prox.Min_Distance_to_Sink2);
    P_dev[target].Min_xyz_to_Sink = mem.sink_prox.Min_xyz_to_Sink;
#ifdef SINGLE_STAR_FIND_BINARIES
    P_dev[target].is_in_a_binary = 0;
    P_dev[target].Min_Sink_OrbitalTime = mem.sink_prox.Min_Sink_OrbitalTime;
    if(mem.sink_prox.Min_Sink_OrbitalTime < MAX_REAL_NUMBER) {
        P_dev[target].is_in_a_binary = 1;
        P_dev[target].comp_Mass = mem.sink_prox.comp_Mass;
        P_dev[target].comp_dx = mem.sink_prox.comp_dx;
        P_dev[target].comp_dv = mem.sink_prox.comp_dv;
    }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    P_dev[target].Min_Sink_Approach_Time = sqrt(mem.sink_prox.Min_Sink_Approach_Time);
    P_dev[target].Min_Sink_Freefall_time = sqrt(sqrt(mem.sink_prox.Min_Sink_Freefall_time) / All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    P_dev[target].Min_Sink_FeedbackTime = sqrt(mem.sink_prox.Min_Sink_FeedbackTime);
#endif
#endif
#endif /* SINK_CALC_DISTANCES */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    P_dev[target].tidal_tensorps = mem.out.tidal_tensorps;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    P_dev[target].GravJerk = mem.out.jerk;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    P_dev[target].TreeMass = mem.tree_mass;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    P_dev[target].tidal_zeta = (MyFloat) mem.out.tidal_zeta;
#endif
#ifdef SPECIAL_POINT_MOTION
    P_dev[target].vel_of_nearest_special = Vec3<MyFloat>{(MyFloat)mem.sink_prox.vel_of_nearest_special[0],
                                                         (MyFloat)mem.sink_prox.vel_of_nearest_special[1],
                                                         (MyFloat)mem.sink_prox.vel_of_nearest_special[2]};
    P_dev[target].acc_of_nearest_special = Vec3<MyFloat>{(MyFloat)mem.sink_prox.acc_of_nearest_special[0],
                                                         (MyFloat)mem.sink_prox.acc_of_nearest_special[1],
                                                         (MyFloat)mem.sink_prox.acc_of_nearest_special[2]};
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    P_dev[target].weight_sum_for_special_point_smoothing = (MyFloat) mem.sink_prox.weight_sum_for_special_point_smoothing;
#endif
#endif
    acc_out = mem.out.acc;
    ninter_out = mem.out.ninter;
    pot_out = mem.out.pot;
}

/* -------------------------------------------------------------------------
 * gpu_gravtree_walk_one -- the walk for a single target: the units above composed
 * with each accepted element evaluated at encounter.
 *
 * Returns 1 on success (outputs written), 0 on failure (pseudo-particle hit;
 * host runs the CPU walk for this target).  Mirrors force_treeevaluate().
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
            if(gpu_grav_leaf_member_accepts(ctx, no, mem.open)) {gpu_grav_evaluate_leaf(ctx, no, mem);}
            no = tree_soa->nextnode_aux[no];
            continue;
        }
        if(no >= pseudo_start) {return 0;} /* pseudo-particle -- remote: host runs the CPU walk for this target */

        gpu_grav_node_prelude_t nd;
        const gpu_grav_node_step_t step = gpu_grav_node_prelude(ctx, no, nd);
        if(step == GPU_GRAV_NODE_SKIP_TO_SIBLING) {no = nd.sibling; continue;}
        if(step == GPU_GRAV_NODE_DESCEND) {no = nd.nextnode; continue;}
        Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
        if(!gpu_grav_node_member_geometry(ctx, nd, mem.open, s_node, mass_node, dr, r2)) {no = nd.sibling; continue;} /* pure-star node, star target */
        int note;
        const gravtree_open_t pred = gpu_grav_node_member_decide(ctx, nd, mem.open, mass_node, r2, note);
        if(note != GPU_GRAV_NOTE_NONE) {gpu_grav_note_commit(1, (note == GPU_GRAV_NOTE_UNSHIPPABLE) ? 1 : 0);}
        if(pred == GRAV_SKIP_NODE) {no = nd.sibling; continue;}
        if(pred == GRAV_OPEN_NODE) {no = nd.nextnode; continue;}
        gpu_grav_evaluate_node(ctx, nd, mem, s_node, mass_node, dr, r2);
        no = nd.sibling;
    }

    gpu_grav_member_finish(ctx, mem, acc_out, ninter_out, pot_out);
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
#define GRAV_PACKET_MASK_WORDS_MAX 4            /* q_dev <= 256: a launch never asks for more */
#define GRAV_PACKET_LOCAL_STACK 16              /* continuations a walker keeps itself; the oldest moves to the frontier when full */
typedef unsigned long long grav_packet_mask_word_t;

struct gpu_grav_walk_item_t { int no, exit; };   /* a work item's indices; its mask words live beside it */

struct gpu_grav_packet_scratch_plan_t {
    int mask_words;
    size_t open_inputs, frontier, frontier_masks, records, record_masks, local, local_masks, counters, bytes;
};

/* The scratch a team needs for one launch shape: q_dev members, team_size threads,
 * a frontier of frontier_cap items and a chunk of chunk_cap records. */
static struct gpu_grav_packet_scratch_plan_t
gpu_grav_packet_scratch_plan(int q_dev, int team_size, int frontier_cap, int chunk_cap)
{
    struct gpu_grav_packet_scratch_plan_t p;
    p.mask_words = (q_dev + GRAV_PACKET_MASK_BITS - 1) / GRAV_PACKET_MASK_BITS;
    size_t off = 0;
    auto take = [&off](size_t bytes, size_t align) {off = ((off + align - 1) / align) * align; size_t here = off; off += bytes; return here;};
    p.open_inputs    = take((size_t) q_dev * sizeof(gpu_grav_open_inputs_t), alignof(gpu_grav_open_inputs_t));
    p.frontier       = take((size_t) frontier_cap * sizeof(gpu_grav_walk_item_t), alignof(gpu_grav_walk_item_t));
    p.frontier_masks = take((size_t) frontier_cap * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.records        = take((size_t) chunk_cap * sizeof(grav_walk_record_t), alignof(grav_walk_record_t));
    p.record_masks   = take((size_t) chunk_cap * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.local          = take((size_t) team_size * GRAV_PACKET_LOCAL_STACK * sizeof(gpu_grav_walk_item_t), alignof(gpu_grav_walk_item_t));
    p.local_masks    = take((size_t) team_size * GRAV_PACKET_LOCAL_STACK * p.mask_words * sizeof(grav_packet_mask_word_t), alignof(grav_packet_mask_word_t));
    p.counters       = take(16 * sizeof(int), alignof(long long));
    p.bytes = off;
    return p;
}

/* team-scope counters, one int slot each */
enum {
    GRAV_PACKET_CTR_RECORDS = 0,       /* records in the chunk */
    GRAV_PACKET_CTR_FRONTIER,          /* items on the frontier */
    GRAV_PACKET_CTR_FAILED,            /* the packet failed: pseudo-particle, malformed index, or no room for a continuation */
    GRAV_PACKET_CTR_DONE,              /* the traversal is finished */
    GRAV_PACKET_CTR_NOTE_INCOMPLETE,   /* import-completeness notes, held until success */
    GRAV_PACKET_CTR_NOTE_UNSHIPPABLE,
    GRAV_PACKET_CTR_COUNT
};

struct GpuGravPacketWalk {
    using TeamMember = Kokkos::TeamPolicy<>::member_type;
    gpu_grav_walk_ctx_t ctx;
    const int *d_idx;   /* candidates in ActiveParticleList order */
    int n_cand, q_dev, frontier_cap, chunk_cap;
    struct gpu_grav_packet_scratch_plan_t plan;
    Vec3<double> *d_acc; int *d_ninter; double *d_pot; int *d_failed;

    KOKKOS_INLINE_FUNCTION static void mask_clear(grav_packet_mask_word_t *m, int words) {for(int w = 0; w < words; w++) {m[w] = 0ULL;}}
    KOKKOS_INLINE_FUNCTION static void mask_copy(grav_packet_mask_word_t *dst, const grav_packet_mask_word_t *src, int words) {for(int w = 0; w < words; w++) {dst[w] = src[w];}}
    KOKKOS_INLINE_FUNCTION static int  mask_test(const grav_packet_mask_word_t *m, int b) {return (int)((m[b / GRAV_PACKET_MASK_BITS] >> (b % GRAV_PACKET_MASK_BITS)) & 1ULL);}
    KOKKOS_INLINE_FUNCTION static void mask_set(grav_packet_mask_word_t *m, int b) {m[b / GRAV_PACKET_MASK_BITS] |= (1ULL << (b % GRAV_PACKET_MASK_BITS));}
    KOKKOS_INLINE_FUNCTION static int  mask_any(const grav_packet_mask_word_t *m, int words) {for(int w = 0; w < words; w++) {if(m[w]) {return 1;}} return 0;}
    KOKKOS_INLINE_FUNCTION static int  mask_equal(const grav_packet_mask_word_t *a, const grav_packet_mask_word_t *b, int words) {for(int w = 0; w < words; w++) {if(a[w] != b[w]) {return 0;}} return 1;}

    /* One member's evaluation of a recorded element: the leaf through the P_dev adapter,
     * a node re-derived from its index through the same prelude and geometry the walker
     * used (the same statements, so the same values), then the shared evaluation. */
    KOKKOS_INLINE_FUNCTION void evaluate_record(int no, gpu_grav_member_t &mem) const
    {
        if(no < ctx.treeParticleSlots) {gpu_grav_evaluate_leaf(ctx, no, mem); return;}
        gpu_grav_node_prelude_t nd;
        (void) gpu_grav_node_prelude(ctx, no, nd);   /* an accepted node is one the prelude handed to the decision */
        Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
        (void) gpu_grav_node_member_geometry(ctx, nd, mem.open, s_node, mass_node, dr, r2);   /* a member with the bit set is never a star seeing a pure-star node */
        gpu_grav_evaluate_node(ctx, nd, mem, s_node, mass_node, dr, r2);
    }

    /* The walker: advance its item until the chunk is full, the traversal is done, or the
     * packet fails. All of its state persists across the flush that follows: the item in
     * the caller's variables, the continuations in scratch. */
    KOKKOS_INLINE_FUNCTION void walk(int &no, int &exit, grav_packet_mask_word_t *mask, int &item_live,
                                     int &local_head, int &local_count,
                                     gpu_grav_walk_item_t *local, grav_packet_mask_word_t *lmasks,
                                     const gpu_grav_open_inputs_t *open,
                                     gpu_grav_walk_item_t *frontier, grav_packet_mask_word_t *fmasks,
                                     grav_walk_record_t *records, grav_packet_mask_word_t *rmasks, int *ctr) const
    {
        const int W = plan.mask_words;
        const int treeBase = ctx.treeBase, treeParticleSlots = ctx.treeParticleSlots;
        const int pseudo_start = treeBase + ctx.maxNodes + ctx.maxForeignNodes;
        grav_packet_mask_word_t accept_mask[GRAV_PACKET_MASK_WORDS_MAX], open_mask[GRAV_PACKET_MASK_WORDS_MAX];

        while(1)
        {
            /* the item is complete: resume the newest continuation, else one from the frontier, else finish */
            if(!item_live || no == exit) {
                if(local_count > 0) {
                    const int slot = (local_head + local_count - 1) % GRAV_PACKET_LOCAL_STACK;
                    no = local[slot].no; exit = local[slot].exit; mask_copy(mask, lmasks + (size_t) slot * W, W);
                    local_count--; item_live = 1; continue;
                }
                if(ctr[GRAV_PACKET_CTR_FRONTIER] > 0) {
                    const int slot = --ctr[GRAV_PACKET_CTR_FRONTIER];   /* one walker: no contention on the frontier */
                    no = frontier[slot].no; exit = frontier[slot].exit; mask_copy(mask, fmasks + (size_t) slot * W, W);
                    item_live = 1; continue;
                }
                ctr[GRAV_PACKET_CTR_DONE] = 1; return;
            }
            if(no < 0) {item_live = 0; continue;}   /* the single-target walk's end-of-walk sentinel: nothing beyond it */
            if(no >= treeParticleSlots && no < treeBase) {ctr[GRAV_PACKET_CTR_FAILED] = 1; return;}   /* gap: malformed tree; the host walk stops loudly */
            if(no < treeParticleSlots) /* particle leaf: per member, the star-star pass */
            {
                mask_clear(accept_mask, W);
                for(int m = 0; m < q_dev; m++) {
                    if(!mask_test(mask, m)) {continue;}
                    if(gpu_grav_leaf_member_accepts(ctx, no, open[m])) {mask_set(accept_mask, m);}
                }
                if(mask_any(accept_mask, W)) {
                    if(ctr[GRAV_PACKET_CTR_RECORDS] == chunk_cap) {return;}   /* chunk full: the flush follows; resume here */
                    const int r = ctr[GRAV_PACKET_CTR_RECORDS]++;
                    records[r].no = no; records[r].kind = GRAV_NODE_LOCAL; records[r].leaf_tag = LET_LEAF_TAG_NODE;
                    mask_copy(rmasks + (size_t) r * W, accept_mask, W);
                }
                no = ctx.tree_soa.nextnode_aux[no];
                continue;
            }
            if(no >= pseudo_start) {ctr[GRAV_PACKET_CTR_FAILED] = 1; return;}   /* pseudo-particle: the host walks every member */

            gpu_grav_node_prelude_t nd;
            const gpu_grav_node_step_t step = gpu_grav_node_prelude(ctx, no, nd);
            if(step == GPU_GRAV_NODE_SKIP_TO_SIBLING) {no = nd.sibling; continue;}
            if(step == GPU_GRAV_NODE_DESCEND) {no = nd.nextnode; continue;}

            /* every member in the item judges the node for itself */
            mask_clear(accept_mask, W); mask_clear(open_mask, W);
            int n_note = 0, n_unship = 0;
            for(int m = 0; m < q_dev; m++) {
                if(!mask_test(mask, m)) {continue;}
                Vec3<MyFloat> s_node; MyFloat mass_node; Vec3<double> dr; double r2;
                if(!gpu_grav_node_member_geometry(ctx, nd, open[m], s_node, mass_node, dr, r2)) {continue;}   /* a star member is done with a pure-star node */
                int note;
                const gravtree_open_t pred = gpu_grav_node_member_decide(ctx, nd, open[m], mass_node, r2, note);
                if(note != GPU_GRAV_NOTE_NONE) {n_note++; if(note == GPU_GRAV_NOTE_UNSHIPPABLE) {n_unship++;}}
                if(pred == GRAV_SKIP_NODE) {continue;}
                if(pred == GRAV_OPEN_NODE) {mask_set(open_mask, m); continue;}
                mask_set(accept_mask, m);
            }
            if(mask_any(accept_mask, W)) {
                if(ctr[GRAV_PACKET_CTR_RECORDS] == chunk_cap) {return;}   /* chunk full: resume at this node after the flush (its decisions are re-made identically) */
                const int r = ctr[GRAV_PACKET_CTR_RECORDS]++;
                records[r].no = no; records[r].kind = nd.node_kind; records[r].leaf_tag = nd.fl_tag;
                mask_copy(rmasks + (size_t) r * W, accept_mask, W);
            }
            ctr[GRAV_PACKET_CTR_NOTE_INCOMPLETE] += n_note; ctr[GRAV_PACKET_CTR_NOTE_UNSHIPPABLE] += n_unship;
            if(!mask_any(open_mask, W)) {no = nd.sibling; continue;}
            if(mask_equal(open_mask, mask, W)) {no = nd.nextnode; continue;}   /* everyone descends: same item, deeper */
            /* the packet descends for the openers; the others re-join at the node's sibling
               through the continuation (sibling, exit, mask). The walker keeps the newest
               continuations itself; the oldest moves to the frontier when there is no room,
               and the frontier is popped after the local ones, so the depth-first order is kept. */
            if(local_count == GRAV_PACKET_LOCAL_STACK) {
                if(ctr[GRAV_PACKET_CTR_FRONTIER] == frontier_cap) {ctr[GRAV_PACKET_CTR_FAILED] = 1; return;}
                const int fslot = ctr[GRAV_PACKET_CTR_FRONTIER]++;
                frontier[fslot] = local[local_head]; mask_copy(fmasks + (size_t) fslot * W, lmasks + (size_t) local_head * W, W);
                local_head = (local_head + 1) % GRAV_PACKET_LOCAL_STACK; local_count--;
            }
            {
                const int slot = (local_head + local_count) % GRAV_PACKET_LOCAL_STACK;
                local[slot].no = nd.sibling; local[slot].exit = exit; mask_copy(lmasks + (size_t) slot * W, mask, W);
                local_count++;
            }
            no = nd.nextnode; exit = nd.sibling; mask_copy(mask, open_mask, W);
        }
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
        gpu_grav_walk_item_t     *local     = (gpu_grav_walk_item_t *)     (scratch + plan.local) + (size_t) t * GRAV_PACKET_LOCAL_STACK;
        grav_packet_mask_word_t  *lmasks    = (grav_packet_mask_word_t *)  (scratch + plan.local_masks) + (size_t) t * GRAV_PACKET_LOCAL_STACK * W;
        int                      *ctr       = (int *)                      (scratch + plan.counters);

        /* the member this thread owns, if any; every thread publishes an entry so the
           walker's loop over q_dev members reads only initialised inputs */
        gpu_grav_member_t mem;
        const int have_member = (t < q_eff);
        if(have_member) {
            (void) gpu_grav_member_init(ctx, d_idx[first + t], mem);
            open[t] = mem.open;
        } else if(t < q_dev) {
            open[t].alive = 0;
        }
        if(t == 0) {for(int c = 0; c < GRAV_PACKET_CTR_COUNT; c++) {ctr[c] = 0;}}
        team.team_barrier();

        /* stage 1: one walker, thread 0 -- the depth-first order of the single-target walk */
        const int is_walker = (t == 0);
        int no = ctx.treeBase, exit = -1;
        grav_packet_mask_word_t mask[GRAV_PACKET_MASK_WORDS_MAX]; mask_clear(mask, W);
        for(int m = 0; m < q_dev; m++) {if(open[m].alive) {mask_set(mask, m);}}
        int item_live = mask_any(mask, W);   /* a packet of massless targets walks nothing */
        int local_head = 0, local_count = 0;

        while(1)
        {
            if(is_walker && !ctr[GRAV_PACKET_CTR_FAILED]) {
                walk(no, exit, mask, item_live, local_head, local_count, local, lmasks, open, frontier, fmasks, records, rmasks, ctr);
            }
            team.team_barrier();
            if(ctr[GRAV_PACKET_CTR_FAILED]) {break;}
            /* the chunk is full, or the traversal has finished: every member evaluates its records */
            if(have_member && mem.open.alive) {
                const int n_rec = ctr[GRAV_PACKET_CTR_RECORDS];
                for(int r = 0; r < n_rec; r++) {
                    if(!mask_test(rmasks + (size_t) r * W, t)) {continue;}
                    evaluate_record(records[r].no, mem);
                }
            }
            team.team_barrier();
            if(ctr[GRAV_PACKET_CTR_DONE]) {break;}
            if(t == 0) {ctr[GRAV_PACKET_CTR_RECORDS] = 0;}
            team.team_barrier();
        }

        /* commit, or discard everything */
        if(ctr[GRAV_PACKET_CTR_FAILED]) {
            if(have_member) {d_failed[first + t] = 1;}
            return;
        }
        if(t == 0) {gpu_grav_note_commit(ctr[GRAV_PACKET_CTR_NOTE_INCOMPLETE], ctr[GRAV_PACKET_CTR_NOTE_UNSHIPPABLE]);}
        if(have_member) {
            Vec3<double> acc = Vec3<double>{0,0,0}; int ninter = 0; double pot = 0.0;
            if(mem.open.alive) {gpu_grav_member_finish(ctx, mem, acc, ninter, pot);}
            d_acc[first + t] = acc; d_ninter[first + t] = ninter; d_pot[first + t] = pot; d_failed[first + t] = 0;
        }
    }
};

/* Host: choose the launch shape and run the packet engine over the candidates.
 * Returns 0 on success (outputs and d_failed filled per candidate), 1 if no legal shape
 * exists for this build (the caller uses the single-target walk). */
static int gpu_gravtree_walk_packets(const gpu_grav_walk_ctx_t &ctx, const int *d_idx, int n_cand,
                                     Vec3<double> *d_acc, int *d_ninter, double *d_pot, int *d_failed)
{
    GpuGravPacketWalk f;
    f.ctx = ctx; f.d_idx = d_idx; f.n_cand = n_cand;
    f.d_acc = d_acc; f.d_ninter = d_ninter; f.d_pot = d_pot; f.d_failed = d_failed;
    /* the packet is as wide as the team: every member has its own thread. Start from the
       configured packet size and halve until the backend can launch the team with the
       scratch it asks for (a legality bound only, gpu_dispatch_templates.h). */
    int want = TREE_QUERY_PACKET_SIZE; if(want > 256) {want = 256;}
    int team = 1; while((team << 1) <= want) {team <<= 1;}
    for(; team >= 1; team >>= 1) {
        f.q_dev = team;
        f.frontier_cap = (2 * team > 16) ? 2 * team : 16;
        f.chunk_cap = (8 * team > 256) ? 8 * team : 256;
        f.plan = gpu_grav_packet_scratch_plan(f.q_dev, team, f.frontier_cap, f.chunk_cap);
        const int league = (n_cand + f.q_dev - 1) / f.q_dev;
        /* the legality bound is asked of a probe policy carrying the same scratch request: a
           policy constructed at an illegal team size throws before it can be asked anything */
        Kokkos::TeamPolicy<> probe(1, 1, 1);
        probe.set_scratch_size(0, Kokkos::PerTeam(f.plan.bytes));
        const int hw = probe.team_size_max(f, Kokkos::ParallelForTag());
        if(hw < team) {continue;}
        Kokkos::TeamPolicy<> policy(league, team, 1);
        policy.set_scratch_size(0, Kokkos::PerTeam(f.plan.bytes));
        Kokkos::parallel_for("gravtree_walk_packets", policy, f);
        Kokkos::fence();
        gizmo_gpu_check_last_error("gravtree_walk_packets", league);
        return 0;
    }
    return 1;
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
    plan.bytes = offset;
    return plan;
}

extern "C" int gpu_gravtree_walk_primary(int *host_candidates_left)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();
    int num_active_total = (int) ActiveParticleList.size();
    /* How many candidates this walk leaves to the host loop: every active until this walk has
     * selected and taken some. The host loop sizes its per-thread packet workspace from it. */
    if(host_candidates_left) {*host_candidates_left = (num_active_total > 0) ? num_active_total : 0;}
    if(Ewald_iter > 0) {return 0;}
    if(num_active_total <= 0) {return 0;}

    /* The CPU walk (forcetree.cc) JIT-drifts particles and nodes whose
     * Ti_current is stale, at the point of encounter in the walk. Without
     * this, inactive particles (outside the currently active timebin) hold
     * stale positions since GIZMO only drifts active bins per sync-point
     * (run.cc:629) and only rebuilds the tree occasionally. The GPU walk
     * cannot call these host-only helpers from inside the Kokkos kernel,
     * so we apply the drift once up-front here: drift all particles, then
     * drift all nodes whose Ti_current lags All.Ti_Current.  The node drift
     * loop is a single GPU kernel that mutates UVM Nodes/Extnodes AND the
     * SoA mirror in one pass — no host loop, no AoS->SoA reseed afterwards.
     * Cost is O(active drifted nodes) with GPU parallelism over
     * Numnodestree (early-out when Ti_current matches). */
    /* Host-side wrapper in the GPU TU must use the out-of-line host accessor
     * `gizmo_host_ti_current()` (defined in core/predict.cc) rather than a
     * bare All.Ti_Current read, so the host-snapshot intent at this call
     * site stays correct even when the device-pass redirect is active. */
    integertime ti_curr_host = gizmo_host_ti_current();
    move_particles(ti_curr_host); /* drifts all P[], invalidates arena */
    /* SoA must exist before the drift kernel — it writes mirror fields. */
    gpu_gravity_tree_acquire(MaxNodes + 1, Nodes_base, Extnodes_base);

    /* Select the particles this walk will actually cover BEFORE drifting any node: the
     * node drift below exists only to serve this walk, so a step that turns out to walk
     * nothing, or that is small enough to be cheaper on the host, must not pay for it.
     * The selection reads only particle state (ActiveParticleList, ProcessedFlag and the
     * candidacy predicate), all of which move_particles above has already brought to
     * ti_curr_host, so it is independent of the node drift that now follows it. */
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

    /* Few enough candidates that the host walk, which drifts nodes only as it opens
     * them, beats this walk plus the all-node drift it requires. Returning with
     * ProcessedFlag untouched leaves every candidate to the host loop in gravtree.cc. */
    if(gravity_walk_route_to_host(num_active)) {myfree(idx_host); return 0;}

    /* Already-current geometry needs no sweep, and asking for one when a host
     * lazy drift armed the latch earlier in the step would fail rather than
     * no-op.  A tree built since that drift is current and its mirror was
     * rewritten with it. */
    if(!gpu_gravity_tree_nodes_current_at(ti_curr_host)
            && gpu_force_drift_nodes(ti_curr_host) != 0) {
        myfree(idx_host);   /* LIFO mymalloc cleanup before drain */
        endrun(929702);
        return 1;   /* soft bad-stop: skip walk on un-drifted nodes; drains at next poll */
    }

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

    /* Scratch arrays for per-target results, carved from one allocation */
    const struct grav_walk_scratch_plan scratch = grav_walk_scratch_plan_for(num_active, 1);
    char *scratch_block = (char *) gizmo_gpu_alloc_shared(scratch.bytes, "gravity_walk");
    if(!scratch_block) {
        printf("gpu_gravtree_walk_primary: kokkos_malloc failed\n");
        endrun(913201);
        release_payload_buffers();
        myfree(idx_host);   /* LIFO mymalloc cleanup before drain */
        return 1;
    }
    int          *d_idx    = (int *)          (scratch_block + scratch.idx);
    int          *d_failed = (int *)          (scratch_block + scratch.failed);
    Vec3<double> *d_acc    = (Vec3<double> *) (scratch_block + scratch.acc);
    int          *d_ninter = (int *)          (scratch_block + scratch.ninter);
    double       *d_pot    = (double *)       (scratch_block + scratch.pot);
    memcpy(d_idx, idx_host, num_active * sizeof(int));
    memset(d_failed, 0, num_active * sizeof(int));

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
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
        release_payload_buffers();
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
    hermite_dev.tables = drift_kick_table_view(NULL, NULL, 0.0, 0.0, All.Timebase_interval, 0);
    /* the mirror is only built on a Hermite pass: the predictor returns immediately when
       HermiteOnlyFlag is 0, so an ordinary pass would be copying a table nothing reads */
    if(HermiteOnlyFlag && drift_kick_table_mirror_refresh(&hermite_drift_kick_table_dev, &hermite_dev.tables) != 0) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_block);
        release_payload_buffers();
        myfree(idx_host);   /* soft bad-stop already requested: no launch without the tables the prediction reads */
        return 1;
    }
#endif

    /* Invariant guard: reset the per-walk counter. */
    g_inv_fterm_aggregate = 0; g_unship_aggregate = 0;

    /* the per-call context the walk reads, captured by value; the pointers it holds are
     * SharedSpace / captured-snapshot addresses valid for the launch */
    gpu_grav_walk_ctx_t ctx;
    ctx.treeBase = treeBase; ctx.treeParticleSlots = treeParticleSlots_snap; ctx.maxNodes = maxNodes_snap; ctx.maxForeignNodes = maxForeignNodes_snap;
    ctx.P_dev = P_dev; ctx.CellP_dev = CellP_dev; ctx.tree_soa = soa_snap;
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    ctx.is_first_step = is_first_step_snap;
#endif
    ctx.pm = pm_snap;
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

    int walked_as_packets = 0;
#ifdef GX_B2_FORCE_ENGINE
    /* SPIKE (stage-1 gate arm, torn down when the dispatch between the single-target walk and
       the packet engine is decided from measurement): route every candidate through the engine. */
    walked_as_packets = (gpu_gravtree_walk_packets(ctx, d_idx, num_active, d_acc, d_ninter, d_pot, d_failed) == 0);
#endif
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
            if(ok) {
                d_acc[a] = acc;
                d_ninter[a] = ninter;
                d_pot[a] = pot;
                d_failed[a] = 0;
            } else {
                d_failed[a] = 1;
            }
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("gravtree_walk_primary", num_active);
    }
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
    /* No SoA invalidate: the next force_treebuild fully repopulates the SoA
     * via the build pipeline; the next pre-walk drift mutates only the
     * stale-Ti_current nodes via gpu_force_drift_nodes (UVM AoS + SoA in one
     * kernel).  No host-side reseed scaffolding remains. */

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
 * written), 0 if a pseudo-particle was encountered (defer to CPU). */
static KOKKOS_INLINE_FUNCTION int
gpu_ewald_walk_one(int target,
                   int treeBase, int treeParticleSlots, int maxNodes, int maxForeignNodes,    /* LET */
                   struct particle_data *P_dev,
                   const struct gpu_gravity_tree_soa_t *tree_soa,
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
        double mass = 0.0;
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
            dr[0] = P_dev[no].Pos[0] - pos[0];
            dr[1] = P_dev[no].Pos[1] - pos[1];
            dr[2] = P_dev[no].Pos[2] - pos[2];
            mass  = P_dev[no].Mass;
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
            mass  = tree_soa->mass[idx];
            dr[0] = tree_soa->s[idx][0] - pos[0];
            dr[1] = tree_soa->s[idx][1] - pos[1];
            dr[2] = tree_soa->s[idx][2] - pos[2];
        }

        /* nearest-image wrap on the displacement (shared SSOT helper) */
        gravity_box_nearest_image(dr[0], dr[1], dr[2], -1);

        if(is_leaf) {
            no = tree_soa->nextnode_aux[no];
        } else {
            /* Opening check + periodic-boundary skip (mirrors forcetree.cc:2769-2842) */
            double r2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if(r2 <= 0) r2 = 1e-300;
            double len = tree_soa->len[idx];
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

    /* This walk consumes the node mirror (len, center, the node threading) but never
     * refreshes it -- it relies on the primary walk having swept the nodes earlier in
     * the same evaluation. When the primary walk ran on the host instead, no sweep
     * happened and the mirror still holds geometry from an earlier time, so walking it
     * here would write Ewald corrections from stale node extents and mark the particles
     * processed, leaving the host correction loop nothing to fix. Decline instead: with
     * ProcessedFlag untouched, that loop takes every particle and drifts nodes as it
     * opens them.
     *
     * Both halves of the test are needed. The certification records that a sweep ran at
     * this time; it does not survive a host lazy drift that happens afterwards at the
     * same time, which refreshes a node in the AoS while leaving its mirror behind. The
     * second clause covers exactly that ordering. */
    if(!gpu_gravity_tree_nodes_current_at(gizmo_host_ti_current())) {
        return 0;
    }

    /* Particles have already been drifted by the primary walk (Ewald_iter==0);
     * the tree SoA mirror is still valid. Just re-acquire (cache-hit).
     * soa->nextnode_aux aliases UVM Nextnode[]; no per-walk memcpy. */
    int min_nodes = MaxNodes + 1;
    gpu_gravity_tree_acquire(min_nodes, Nodes_base, Extnodes_base);
    const struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa) {return 0;}

    gpu_particles_arena_set_site("gpu_gravtree_walk_ewald");
    gpu_particles_arena_acquire(NumPart, P, CellP);
    struct particle_data *P_dev = gpu_particles_arena_P();
    if(!P_dev) {return 0;}

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

    const struct grav_walk_scratch_plan scratch = grav_walk_scratch_plan_for(num_active, 0);
    char *scratch_block = (char *) gizmo_gpu_alloc_shared(scratch.bytes, "gravity_walk");
    if(!scratch_block) {printf("gpu_ewald_walk_primary: kokkos_malloc failed\n"); endrun(914102); myfree(idx_host); return 1;}
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
                                     P_dev, &soa_snap,
#ifdef GRAVITY_HYBRID_OPENING_CRIT
                                     is_first_step_snap,
#endif
                                     fcorrx_dev, fcorry_dev, fcorrz_dev,
                                     fac_intp, boxsize, boxhalf,
                                     errtoltheta, errtolforceacc,
                                     acc);
        if(ok) {d_acc[a] = acc; d_failed[a] = 0;}
        else   {d_failed[a] = 1;}
    });
    Kokkos::fence();
    gizmo_gpu_check_last_error("gpu_ewald_walk_primary", num_active);
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


