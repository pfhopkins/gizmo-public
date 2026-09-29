#ifndef GRAVTREE_MOMENT_KERNEL_H
#define GRAVTREE_MOMENT_KERNEL_H

/* Node multipole/payload CONSTRUCTION physics, and node MOTION between builds (the drift arithmetic,
 * node_motion_fold_kick / node_motion_advance near the end of this file) — the single home for the per-node moment
 * accumulation formulas shared by the live GPU construction venues: the local-tree refresh
 * (gpu_moment_refresh.cc, atomic into shared scratch) and the topnode re-sum
 * (gpu_pseudo_update.cc::topnode_resum_node_, plain into a local accumulator). Companion to
 * gravtree_force_kernel.h (the per-interaction CONTRIBUTION) and gravtree_opening.h (the node
 * opening DECISION).
 *
 * Contract (mirrors the contract in gravtree_force_kernel.h): helpers operate on small POD accumulators and
 * source structs. Each venue OWNS loading its sources (from P_dev / scratch / SoA) and storing the
 * finished node (to scratch / SoA / AoS), casting to/from its native storage types; the helpers own
 * ONLY the physics/numerics between. This is what lets one physics body serve both the atomic and
 * the plain venue and stay bit-exact: the only per-venue difference is the WRITE POLICY (plain add /
 * Kokkos atomic add) injected as a template argument, plus the caller-side load casts.
 *
 * Kokkos-header discipline: this header includes NO
 * Kokkos header and names NO Kokkos symbol. The plain write policy (moment_plain_ops) lives here
 * unconditionally; the atomic policy (moment_atomic_ops) is defined in gpu_moment_refresh.cc where
 * Kokkos is already in scope, and passed in as the Ops template argument. forcetree.cc therefore
 * acquires NO Kokkos dependency from including this header (it consumes only the plain per-particle
 * hmax/vmax bound primitives at force_add_element_to_tree).
 *
 * Three source forms (_source_from_leaf / _source_from_node split):
 *   - from_particle           : a leaf particle's contribution (gas/sink/RT/sink-distance gates).
 *   - from_child_raw          : a child whose moments are still UN-normalized (Σ m x), as produced
 *                               by gpu_moment_refresh's scratch before the single end-of-pass
 *                               normalize — added directly.
 *   - from_child_normalized   : a child whose moments are already normalized (COM), as stored in the
 *                               SoA — re-weighted by the child's own mass / luminosity / sink-mass to
 *                               reconstruct the Σ-weighted contribution.
 * All three produce a moment_node_accum, consumed by ONE accumulate body (moment_accum_apply): sum
 * fields add, max fields fmax, COM-vector fields add (pre-weighted by the source builder). Fields a
 * given source does not touch are left at the neutral element (0), which is a no-op for + and for fmax
 * (every max field here is >= 0 by construction), so the blind per-field apply is bit-exact.
 *
 * finalize (Σ m x -> COM divide, BITFLAG_MULTIPLEPARTICLES patch) turns the accumulated sums into the
 * stored normalized payloads. Normalization follows the legacy CPU anchor: COM-style vectors divide by
 * their weight, the tidal tensor keeps its 1/(mass+eps) reciprocal-multiply, and the RT luminosity
 * denominator is the bin-order stellar_lum sum. Those divide-vs-reciprocal and order choices are
 * FP-equivalent normalization detail, not load-bearing physics.
 *
 * Adding a payload — checklist (every step in the SAME commit):
 *   1. add the field to moment_node_accum + moment_node_ref below, under its exact #ifdef.
 *   2. zero it in moment_accum_zero.
 *   3. produce its contribution in from_particle (with the leaf gate) AND from_child_raw AND
 *      from_child_normalized (with the child weight).
 *   4. apply it in moment_accum_apply with the correct op (Ops::add / Ops::fmax).
 *   5. finalize it in moment_finalize (the normalize/divide), if it has a normalized form.
 *   6. wire venue storage: scratch View / SoA / AoS load+store at each caller.
 *
 * A node field that DRIFTS between builds (a centre moved by a velocity, with a pending momentum folded
 * into that velocity) is added to the node-motion accessor(s), node_motion_fold_kick and
 * node_motion_advance, and nowhere else: both drift venues run those.
 *
 * The four LOCKSTEP sites for a payload field are moment_node_accum (values), moment_node_ref
 * (pointer mirror), moment_accum_zero, and moment_accum_apply: each must carry the same fields
 * under the same #ifdef gates. Accum/ref/apply drift is mostly caught at COMPILE time — apply
 * names both a.<field> and r.<field>, so a field missing from one struct (or mistyped between
 * them) fails to build. The residual SILENT risks are forgetting to zero a field or forgetting to
 * apply it; those are caught by the non-vacuous moment-construction comparison gate, which any
 * payload edit must pass. A grouping macro to fold these four lists into a single field manifest
 * was considered and deliberately rejected: it would hide this readable data layout behind macro
 * machinery, and the compile + comparison coverage above already removes most of the alignment
 * hazard. Keep the four lists explicit and field-aligned; on any new payload, run the
 * moment-construction gate with that payload actually populated.
 *
 * Include AFTER declarations/allvars.h (Vec3, particle_data, BITFLAG_*, MIN_REAL_NUMBER, the
 * payload-count macros N_RT_FREQ_BINS / CHIMES_LOCAL_UV_NBINS, and the config flags are consumed from
 * the including TU, mirroring gravtree_opening.h / gravtree_force_kernel.h).
 */

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif


/* ==========================================================================================
 * Write policies. update mode is the only template parameter that distinguishes the two venues.
 *   moment_plain_ops  : sequential, single owner of the slot (topnode re-sum local accumulator,
 *                       gpu_moment_refresh zero kernel where each thread owns its slot k).
 *   moment_atomic_ops : Kokkos atomic add / CAS-max into shared scratch (gpu_moment_refresh's
 *                       particle pass + bottom-up walk). Defined in gpu_moment_refresh.cc.
 * ========================================================================================== */
struct moment_plain_ops {
    template <class T> KOKKOS_INLINE_FUNCTION static void add(T *dst, T v)  { *dst += v; }
    template <class T> KOKKOS_INLINE_FUNCTION static void add_vec3(Vec3<T> *dst, const Vec3<T>& v) { *dst += v; }
    template <class T> KOKKOS_INLINE_FUNCTION static void fmax(T *dst, T v) { if(v > *dst) { *dst = v; } }
    KOKKOS_INLINE_FUNCTION static void add_long(long *dst, long v) { *dst += v; }
    KOKKOS_INLINE_FUNCTION static void add_int (int  *dst, int  v) { *dst += v; }
};


/* ==========================================================================================
 * Per-node accumulator (by value) and a pointer view of one node's storage (write-through). The
 * accum carries the Σ-weighted (un-normalized) payloads; the ref points at whatever storage the
 * venue keeps a node's moments in (a local accum struct for topnode; the scratch View slots for
 * gpu_moment_refresh). bitflags is NOT a summed payload — it is reached only by zero (and, in c2,
 * finalize); it lives in the ref, not the accum.
 * ========================================================================================== */
template <class AccT>
struct moment_node_accum {
    AccT       mass;
    Vec3<AccT> s;        /* Σ m x  (pre-normalize) */
    Vec3<AccT> vs;       /* Σ m v */
    long       Npart;
    AccT       hmax;
    AccT       vmax;
    AccT       divVmax;
    AccT       maxsoft;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    AccT       gasmass;
#endif
#ifdef RT_USE_GRAVTREE
    AccT       stellar_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
    double     chimes_G0 [CHIMES_LOCAL_UV_NBINS];
    double     chimes_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<AccT> rt_s;     /* Σ l x */
    Vec3<AccT> rt_vs;    /* Σ l v */
#endif
#ifdef SINK_PHOTONMOMENTUM
    AccT       sink_lum;
    Vec3<AccT> sink_lum_grad;  /* Σ lum * angle */
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    AccT       cr_inject;
#endif
#ifdef SINK_CALC_DISTANCES
    AccT       sink_mass;
    Vec3<AccT> sink_pos;       /* Σ sm x */
#if defined(SINK_NODE_MOTION_TRACKED)
    int        N_SINK;
    Vec3<AccT> sink_vel;       /* Σ sm v */
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<AccT> sink_acc;       /* Σ sm a */
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    AccT       max_fbvel;
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    AccT       tidal[6];       /* Σ m tidal */
#endif
#ifdef DM_SCALARFIELD_SCREENING
    AccT       mass_dm;
    Vec3<AccT> s_dm;           /* Σ mdm x */
    Vec3<AccT> vs_dm;          /* Σ mdm v */
#endif
};

template <class AccT>
struct moment_node_ref {
    AccT          *mass;
    Vec3<AccT>    *s;
    Vec3<AccT>    *vs;
    long          *Npart;
    AccT          *hmax;
    AccT          *vmax;
    AccT          *divVmax;
    AccT          *maxsoft;
    unsigned int  *bitflags;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    AccT          *gasmass;
#endif
#ifdef RT_USE_GRAVTREE
    AccT          *stellar_lum;   /* base of N_RT_FREQ_BINS */
#ifdef CHIMES_STELLAR_FLUXES
    double        *chimes_G0;      /* base of CHIMES_LOCAL_UV_NBINS */
    double        *chimes_ion;
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<AccT>    *rt_s;
    Vec3<AccT>    *rt_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
    AccT          *sink_lum;
    Vec3<AccT>    *sink_lum_grad;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    AccT          *cr_inject;
#endif
#ifdef SINK_CALC_DISTANCES
    AccT          *sink_mass;
    Vec3<AccT>    *sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
    int           *N_SINK;
    Vec3<AccT>    *sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<AccT>    *sink_acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    AccT          *max_fbvel;
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    AccT          *tidal;          /* base of 6 */
#endif
#ifdef DM_SCALARFIELD_SCREENING
    AccT          *mass_dm;
    Vec3<AccT>    *s_dm;
    Vec3<AccT>    *vs_dm;
#endif
};


/* ==========================================================================================
 * Per-particle bound primitives — the textually-duplicated hmax/vmax formulas. Also consumed by
 * force_add_element_to_tree (forcetree.cc), which incrementally bounds a father node by one inserted
 * leaf. Plain scalars only -> no Kokkos, no storage knowledge.
 * ========================================================================================== */

/* Gas leaf's hmax contribution: its smoothing length capped at the global max kernel radius. */
KOKKOS_INLINE_FUNCTION static double moment_gas_hmax_from_kernelradius(double kernel_radius, double max_kernel_radius)
{
    return (max_kernel_radius < kernel_radius) ? max_kernel_radius : kernel_radius;
}

/* ==========================================================================================
 * Source builders. Each fills a moment_node_accum with one source's Σ-weighted contribution. Fields
 * the source does not touch are left 0 (neutral for + and, since every max field is >= 0, for fmax).
 * ========================================================================================== */

/* Component-wise weighted vector cast: reproduces atomic_add_vec3<AccT, NativeT>(&dst, w * x), i.e.
 * the product w*x is formed in the NATIVE (double) operand precision and only THEN cast to AccT.
 * Doing the multiply in AccT instead would (under GIZMO_MIXED_PRECISION_GRAVITY, AccT=float) round
 * the operands first — a precision downgrade and a bit difference. Mass/positions/velocities are all
 * double in this codebase, so the product is always double here. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static Vec3<AccT> moment_weighted_vec3(double w, double x0, double x1, double x2)
{
    return Vec3<AccT>{(AccT)(w * x0), (AccT)(w * x1), (AccT)(w * x2)};
}

#ifdef RT_SEPARATELY_TRACK_LUMPOS
/* Total stellar luminosity of a node (the RT lum-position normalization weight). The SAME running
 * sum a venue keeps for its rt_s/rt_vs divide MUST use this, so the per-child weight and the
 * normalize denominator stay bit-consistent within the venue. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static double moment_child_total_luminosity(const moment_node_accum<AccT>& c)
{
    double l = 0;
    for(int b = 0; b < N_RT_FREQ_BINS; b++) { l += (double) c.stellar_lum[b]; }
    return l;
}
#endif

/* Leaf particle. The walk loads the particle fields + per-particle precomputes into this POD (native
 * double, matching particle_data); the helper owns the physics gates (Type, luminosity,
 * sink-distance) and the bound/weight formulas, forming weighted products in double and casting the
 * accumulated contribution to AccT. */
template <class AccT>
struct moment_particle_src {
    double        mass;
    double        pos[3];    /* position (already in node/box frame) */
    double        vel[3];    /* the velocity that enters the node's mass-weighted mean */
    double        motion_bound;  /* particle_motion_speed_bound: how fast it can move, per unit undilated interval */
    int           type;
    double        kernel_radius;
    double        max_kernel_radius;
    double        force_softening;   /* gpu_force_softening_kernelradius(pa) */
    double        particle_divvel;
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
    double        sink_mass_reservoir;
#endif
#ifdef RT_USE_GRAVTREE
    double        src_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
    double        src_lum_G0 [CHIMES_LOCAL_UV_NBINS];
    double        src_lum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    double        bh_lum;
    double        bh_angle[3];
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double        cr_inject;
#endif
#if defined(SPECIAL_POINT_MOTION)
    double        acc_prevstep[3];
#endif
#if defined(SINK_CALC_DISTANCES) && defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    double        max_feedback_vel;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    double        tidal_prevstep[6];
#endif
};

template <class AccT>
KOKKOS_INLINE_FUNCTION static moment_node_accum<AccT> moment_source_from_particle(const moment_particle_src<AccT>& p)
{
    moment_node_accum<AccT> a = {};

    a.mass  = (AccT) p.mass;
    a.s     = moment_weighted_vec3<AccT>(p.mass, p.pos[0], p.pos[1], p.pos[2]);
    a.vs    = moment_weighted_vec3<AccT>(p.mass, p.vel[0], p.vel[1], p.vel[2]);
    a.Npart = (long) 1;

    a.vmax = (AccT) p.motion_bound;

    a.maxsoft = (AccT) p.force_softening;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    if(p.type == 5) {
        if((AccT) p.kernel_radius > a.maxsoft) { a.maxsoft = (AccT) p.kernel_radius; }
    }
#endif

    if(p.type == 0) {
        a.hmax    = (AccT) moment_gas_hmax_from_kernelradius(p.kernel_radius, p.max_kernel_radius);
        a.divVmax = (AccT) p.particle_divvel;
    }

#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    if(p.type == 0) { a.gasmass = (AccT) p.mass; }
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
    if(p.type == 5) { a.gasmass = (AccT) p.sink_mass_reservoir; }
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    a.cr_inject = (AccT) p.cr_inject;
#endif
#ifdef RT_USE_GRAVTREE
    {
        double l_sum = 0;
        for(int b = 0; b < N_RT_FREQ_BINS; b++) {
            a.stellar_lum[b] = (AccT) p.src_lum[b];
            l_sum += p.src_lum[b];
        }
#ifdef CHIMES_STELLAR_FLUXES
        for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) {
            a.chimes_G0 [b] = p.src_lum_G0 [b];
            a.chimes_ion[b] = p.src_lum_ion[b];
        }
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        if(l_sum > 0) {
            a.rt_s  = moment_weighted_vec3<AccT>(l_sum, p.pos[0], p.pos[1], p.pos[2]);
            a.rt_vs = moment_weighted_vec3<AccT>(l_sum, p.vel[0], p.vel[1], p.vel[2]);
        }
#endif
    }
#endif
#ifdef SINK_PHOTONMOMENTUM
    if(p.type == 5 && p.mass > 0 && p.bh_lum > 0) {
        a.sink_lum      = (AccT) p.bh_lum;
        a.sink_lum_grad = moment_weighted_vec3<AccT>(p.bh_lum, p.bh_angle[0], p.bh_angle[1], p.bh_angle[2]);
    }
#endif
#ifdef SINK_CALC_DISTANCES
    if(p.type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) {
        a.sink_mass = (AccT) p.mass;
        a.sink_pos  = moment_weighted_vec3<AccT>(p.mass, p.pos[0], p.pos[1], p.pos[2]);
#if defined(SINK_NODE_MOTION_TRACKED)
        a.N_SINK   = 1;
        a.sink_vel = moment_weighted_vec3<AccT>(p.mass, p.vel[0], p.vel[1], p.vel[2]);
#endif
#ifdef SPECIAL_POINT_MOTION
        a.sink_acc = moment_weighted_vec3<AccT>(p.mass, p.acc_prevstep[0], p.acc_prevstep[1], p.acc_prevstep[2]);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        a.max_fbvel = (AccT) p.max_feedback_vel;
#endif
    }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    for(int kk = 0; kk < 6; kk++) { a.tidal[kk] = (AccT)(p.mass * p.tidal_prevstep[kk]); }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    if(p.type != 0) {
        a.mass_dm = (AccT) p.mass;
        a.s_dm    = moment_weighted_vec3<AccT>(p.mass, p.pos[0], p.pos[1], p.pos[2]);
        a.vs_dm   = moment_weighted_vec3<AccT>(p.mass, p.vel[0], p.vel[1], p.vel[2]);
    }
#endif
    return a;
}

/* Child whose moments are still un-normalized (gpu_moment_refresh scratch, before the single
 * end-of-pass normalize): add its accumulated values directly. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static moment_node_accum<AccT> moment_source_from_child_raw(const moment_node_accum<AccT>& c)
{
    moment_node_accum<AccT> a = c;   /* sum/vec fields copy directly; max fields carried for fmax */
#if defined(SINK_CALC_DISTANCES) && defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    /* max_fbvel only participates when the child carries sink mass (mirrors the gated max). */
    if(!(c.sink_mass > 0)) { a.max_fbvel = (AccT) 0; }
#endif
    return a;
}

/* Child whose moments are already normalized (COM, as stored in the SoA): reconstruct its
 * Σ-weighted contribution by re-multiplying each COM-vector field by the child's own weight. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static moment_node_accum<AccT> moment_source_from_child_normalized(const moment_node_accum<AccT>& c)
{
    moment_node_accum<AccT> a = {};

    a.mass    = c.mass;
    a.s       = c.mass * c.s;
    a.vs      = c.mass * c.vs;
    a.Npart   = c.Npart;
    a.hmax    = c.hmax;
    a.vmax    = c.vmax;
    a.divVmax = c.divVmax;
    a.maxsoft = c.maxsoft;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    a.gasmass = c.gasmass;
#endif
#ifdef RT_USE_GRAVTREE
    for(int b = 0; b < N_RT_FREQ_BINS; b++) { a.stellar_lum[b] = c.stellar_lum[b]; }
#ifdef CHIMES_STELLAR_FLUXES
    for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) { a.chimes_G0[b] = c.chimes_G0[b]; a.chimes_ion[b] = c.chimes_ion[b]; }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    {
        double l_child = moment_child_total_luminosity<AccT>(c);
        a.rt_s  = moment_weighted_vec3<AccT>(l_child, (double)c.rt_s[0],  (double)c.rt_s[1],  (double)c.rt_s[2]);
        a.rt_vs = moment_weighted_vec3<AccT>(l_child, (double)c.rt_vs[0], (double)c.rt_vs[1], (double)c.rt_vs[2]);
    }
#endif
#ifdef SINK_PHOTONMOMENTUM
    a.sink_lum      = c.sink_lum;
    a.sink_lum_grad = c.sink_lum * c.sink_lum_grad;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    a.cr_inject = c.cr_inject;
#endif
#ifdef SINK_CALC_DISTANCES
    a.sink_mass = c.sink_mass;
    a.sink_pos  = c.sink_mass * c.sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
    a.N_SINK   = c.N_SINK;
    a.sink_vel = c.sink_mass * c.sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
    a.sink_acc = c.sink_mass * c.sink_acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    if(c.sink_mass > 0) { a.max_fbvel = c.max_fbvel; }
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    for(int kk = 0; kk < 6; kk++) { a.tidal[kk] = c.mass * c.tidal[kk]; }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    a.mass_dm = c.mass_dm;
    a.s_dm    = c.mass_dm * c.s_dm;
    a.vs_dm   = c.mass_dm * c.vs_dm;
#endif
    return a;
}


/* ==========================================================================================
 * The three named construction ops.
 * ========================================================================================== */

/* Zero a node's accumulators. bitflags keeps only the topology bits (matches the CPU saved_bitflags
 * mask in force_refresh_node_moments). Plain stores: the caller's slot is single-owner here. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static void moment_accum_zero(const moment_node_ref<AccT>& r, unsigned int saved_bitflags)
{
    *r.mass    = (AccT) 0;
    *r.s       = Vec3<AccT>{};
    *r.vs      = Vec3<AccT>{};
    *r.Npart   = (long) 0;
    *r.hmax    = (AccT) 0;
    *r.vmax    = (AccT) 0;
    *r.divVmax = (AccT) 0;
    *r.maxsoft = (AccT) 0;
    *r.bitflags = saved_bitflags & ((1u << BITFLAG_TOPLEVEL) |
                                     (1u << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT) |
                                     (1u << BITFLAG_INTERNAL_TOPLEVEL));
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    *r.gasmass = (AccT) 0;
#endif
#ifdef RT_USE_GRAVTREE
    for(int b = 0; b < N_RT_FREQ_BINS; b++) { r.stellar_lum[b] = (AccT) 0; }
#ifdef CHIMES_STELLAR_FLUXES
    for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) { r.chimes_G0[b] = 0; r.chimes_ion[b] = 0; }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    *r.rt_s  = Vec3<AccT>{};
    *r.rt_vs = Vec3<AccT>{};
#endif
#ifdef SINK_PHOTONMOMENTUM
    *r.sink_lum      = (AccT) 0;
    *r.sink_lum_grad = Vec3<AccT>{};
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    *r.cr_inject = (AccT) 0;
#endif
#ifdef SINK_CALC_DISTANCES
    *r.sink_mass = (AccT) 0;
    *r.sink_pos  = Vec3<AccT>{};
#if defined(SINK_NODE_MOTION_TRACKED)
    *r.N_SINK   = 0;
    *r.sink_vel = Vec3<AccT>{};
#endif
#if defined(SPECIAL_POINT_MOTION)
    *r.sink_acc = Vec3<AccT>{};
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    *r.max_fbvel = (AccT) 0;
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    for(int kk = 0; kk < 6; kk++) { r.tidal[kk] = (AccT) 0; }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    *r.mass_dm = (AccT) 0;
    *r.s_dm    = Vec3<AccT>{};
    *r.vs_dm   = Vec3<AccT>{};
#endif
}

/* Apply a Σ-weighted source contribution to a node: sum fields add, max fields fmax, COM-vector
 * fields add (the source builder already applied the weight). One body; the write policy (plain or
 * atomic) is the only per-venue knob. */
template <class Ops, class AccT>
KOKKOS_INLINE_FUNCTION static void moment_accum_apply(const moment_node_ref<AccT>& r, const moment_node_accum<AccT>& a)
{
    Ops::add(r.mass, a.mass);
    Ops::add_vec3(r.s,  a.s);
    Ops::add_vec3(r.vs, a.vs);
    Ops::add_long(r.Npart, a.Npart);
    Ops::fmax(r.vmax,    a.vmax);
    Ops::fmax(r.maxsoft, a.maxsoft);
    Ops::fmax(r.hmax,    a.hmax);
    Ops::fmax(r.divVmax, a.divVmax);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    Ops::add(r.gasmass, a.gasmass);
#endif
#ifdef RT_USE_GRAVTREE
    for(int b = 0; b < N_RT_FREQ_BINS; b++) { Ops::add(&r.stellar_lum[b], a.stellar_lum[b]); }
#ifdef CHIMES_STELLAR_FLUXES
    for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) {
        Ops::add(&r.chimes_G0 [b], a.chimes_G0 [b]);
        Ops::add(&r.chimes_ion[b], a.chimes_ion[b]);
    }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Ops::add_vec3(r.rt_s,  a.rt_s);
    Ops::add_vec3(r.rt_vs, a.rt_vs);
#endif
#ifdef SINK_PHOTONMOMENTUM
    Ops::add(r.sink_lum, a.sink_lum);
    Ops::add_vec3(r.sink_lum_grad, a.sink_lum_grad);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    Ops::add(r.cr_inject, a.cr_inject);
#endif
#ifdef SINK_CALC_DISTANCES
    Ops::add(r.sink_mass, a.sink_mass);
    Ops::add_vec3(r.sink_pos, a.sink_pos);
#if defined(SINK_NODE_MOTION_TRACKED)
    Ops::add_int(r.N_SINK, a.N_SINK);
    Ops::add_vec3(r.sink_vel, a.sink_vel);
#endif
#if defined(SPECIAL_POINT_MOTION)
    Ops::add_vec3(r.sink_acc, a.sink_acc);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    Ops::fmax(r.max_fbvel, a.max_fbvel);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    for(int kk = 0; kk < 6; kk++) { Ops::add(&r.tidal[kk], a.tidal[kk]); }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Ops::add(r.mass_dm, a.mass_dm);
    Ops::add_vec3(r.s_dm,  a.s_dm);
    Ops::add_vec3(r.vs_dm, a.vs_dm);
#endif
}

/* Convenience wrappers: the three named add ops (§38.16). */
template <class Ops, class AccT>
KOKKOS_INLINE_FUNCTION static void moment_accum_add_particle(const moment_node_ref<AccT>& r, const moment_particle_src<AccT>& p)
{
    moment_accum_apply<Ops, AccT>(r, moment_source_from_particle<AccT>(p));
}

template <class Ops, class AccT>
KOKKOS_INLINE_FUNCTION static void moment_accum_add_child_raw(const moment_node_ref<AccT>& r, const moment_node_accum<AccT>& child)
{
    moment_accum_apply<Ops, AccT>(r, moment_source_from_child_raw<AccT>(child));
}

template <class Ops, class AccT>
KOKKOS_INLINE_FUNCTION static void moment_accum_add_child_normalized(const moment_node_ref<AccT>& r, const moment_node_accum<AccT>& child)
{
    moment_accum_apply<Ops, AccT>(r, moment_source_from_child_normalized<AccT>(child));
}


/* Finalize a fully-accumulated node in place: turn the mass-weighted sums into the stored
 * (normalized) payloads and set the multiple-particles flag. Plain by construction — every venue
 * finalizes a node it solely owns, so no write policy is needed. `center` is the node geometric
 * center, used as the COM fallback when the relevant weight is zero. If `r.bitflags` is null the
 * caller patches the multiple-particles bit itself (the topnode venue writes that bit into two
 * separate storages).
 *
 * Normalization form follows the legacy CPU anchor (force_update_node_recursive): COM-style vectors
 * DIVIDE by their weight; the tidal tensor keeps the 1/(mass+eps) reciprocal-multiply. The RT
 * luminosity-position denominator is the bin-order sum of stellar_lum (the CPU/local-tree form);
 * routing every venue through it drops the topnode child-order running sum. The divide-vs-reciprocal
 * and summation-order choices are FP-equivalent normalization detail, not load-bearing physics. */
template <class AccT>
KOKKOS_INLINE_FUNCTION static void moment_finalize(const moment_node_ref<AccT>& r, const Vec3<AccT>& center)
{
    AccT mass = *r.mass;
    if(mass > 0) {
        *r.s  = *r.s  / mass;
        *r.vs = *r.vs / mass;
    } else {
        *r.s  = center;
        *r.vs = Vec3<AccT>{};
    }
    if(r.bitflags) {
        if(*r.Npart > 1) { *r.bitflags |=  (1u << BITFLAG_MULTIPLEPARTICLES); }
        else             { *r.bitflags &= ~(1u << BITFLAG_MULTIPLEPARTICLES); }
    }
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    {
        double l_tot = 0;
        for(int b = 0; b < N_RT_FREQ_BINS; b++) { l_tot += (double) r.stellar_lum[b]; }
        if(l_tot > 0) {
            *r.rt_s  = *r.rt_s  / ((AccT) l_tot);
            *r.rt_vs = *r.rt_vs / ((AccT) l_tot);
        } else {
            *r.rt_s  = center;
            *r.rt_vs = Vec3<AccT>{};
        }
    }
#endif
#ifdef SINK_PHOTONMOMENTUM
    if(*r.sink_lum > 0) { *r.sink_lum_grad = *r.sink_lum_grad / *r.sink_lum; }
    else                { *r.sink_lum_grad = Vec3<AccT>{0, 0, 1}; }
#endif
#ifdef SINK_CALC_DISTANCES
    if(*r.sink_mass > 0) {
        *r.sink_pos = *r.sink_pos / *r.sink_mass;
#if defined(SINK_NODE_MOTION_TRACKED)
        *r.sink_vel = *r.sink_vel / *r.sink_mass;
#endif
#if defined(SPECIAL_POINT_MOTION)
        *r.sink_acc = *r.sink_acc / *r.sink_mass;
#endif
    }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    if(mass > 0) {
        AccT inv_m = (AccT) (1.0 / ((double)mass + MIN_REAL_NUMBER));   /* reciprocal-multiply: matches the legacy CPU anchor */
        for(int kk = 0; kk < 6; kk++) { r.tidal[kk] *= inv_m; }
    }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    if(*r.mass_dm > 0) {
        *r.s_dm  = *r.s_dm  / *r.mass_dm;
        *r.vs_dm = *r.vs_dm / *r.mass_dm;
    } else {
        *r.s_dm  = center;
        *r.vs_dm = Vec3<AccT>{};
    }
#endif
}


/* ==========================================================================================
 * Node MOTION between tree builds: the one home of a node's drift arithmetic.
 *
 * Between builds a node carries its centre of mass (and, under their flags, its luminosity-weighted,
 * dark-matter and special-particle centres) with a velocity, plus the momentum its members were
 * kicked by since it last moved. Drifting it folds that pending momentum into the velocities, then
 * advances every centre, and widens its length by how far any member can have moved (on the
 * undilated clock, which the caller supplies as dt_widen). WHETHER to fold is the caller's decision:
 * the host lazy drift folds only when it moves a node forward (a kicked node that is already current
 * keeps its velocity and pending momentum until its next drift), while the device sweep's full mirror
 * refresh also folds a current node at dt = 0. Both reach the same state at the next real drift; the
 * node velocity read in between differs.
 *
 * The host lazy drift (force_drift_node) and the device sweep (gpu_node_drift_apply) both run
 * exactly this; each keeps its own bookkeeping (drift factors, time stamps, flags, locking). The
 * arithmetic reaches the node through an accessor that reads and writes each field where it lives
 * (node_motion_in_arrays below, for the node arrays), one field at a time, as the moment helpers
 * above reach their storage through moment_node_ref.
 * ========================================================================================== */
struct node_motion_in_arrays {
    struct NODE *nodes; struct extNODE *ext; int no;
    KOKKOS_INLINE_FUNCTION MyFloat  &s(int j)    const {return nodes[no].u.d.s[j];}
    KOKKOS_INLINE_FUNCTION MyFloat  &vs(int j)   const {return ext[no].vs[j];}
    KOKKOS_INLINE_FUNCTION MyDouble &dp(int j)   const {return ext[no].dp[j];}
    KOKKOS_INLINE_FUNCTION MyFloat  &len()       const {return nodes[no].len;}
    KOKKOS_INLINE_FUNCTION double    mass()      const {return (double) nodes[no].u.d.mass;}
    KOKKOS_INLINE_FUNCTION double    vmax()      const {return (double) ext[no].vmax;}
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    KOKKOS_INLINE_FUNCTION MyFloat  &rt_s(int j)  const {return nodes[no].rt_source_lum_s[j];}
    KOKKOS_INLINE_FUNCTION MyFloat  &rt_vs(int j) const {return ext[no].rt_source_lum_vs[j];}
    KOKKOS_INLINE_FUNCTION MyDouble &rt_dp(int j) const {return ext[no].rt_source_lum_dp[j];}
    KOKKOS_INLINE_FUNCTION double    lum_tot()    const {double l = 0; for(int b = 0; b < N_RT_FREQ_BINS; b++) {l += (double) nodes[no].stellar_lum[b];} return l;}
#endif
#ifdef DM_SCALARFIELD_SCREENING
    KOKKOS_INLINE_FUNCTION MyFloat  &s_dm(int j)  const {return nodes[no].s_dm[j];}
    KOKKOS_INLINE_FUNCTION MyFloat  &vs_dm(int j) const {return ext[no].vs_dm[j];}
    KOKKOS_INLINE_FUNCTION MyDouble &dp_dm(int j) const {return ext[no].dp_dm[j];}
    KOKKOS_INLINE_FUNCTION double    mass_dm()    const {return (double) nodes[no].mass_dm;}
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    KOKKOS_INLINE_FUNCTION MyFloat  &sink_pos(int j) const {return nodes[no].sink_pos[j];}
    KOKKOS_INLINE_FUNCTION MyFloat  &sink_vel(int j) const {return nodes[no].sink_vel[j];}
    KOKKOS_INLINE_FUNCTION MyDouble &sink_dp(int j)  const {return ext[no].sink_dp[j];}
    KOKKOS_INLINE_FUNCTION double    sink_mass()     const {return (double) nodes[no].sink_mass;}
#endif
};

/* Turn the pending momentum into velocity: each velocity is a mean over the particles it describes,
 * so its momentum is divided by their mass (luminosity for the luminosity-weighted centre, the
 * special-type mass for the sink centre). The momentum is then spent. Node masses and luminosities
 * are finite and non-negative; a zero weight means no particles of that kind, so there is no velocity
 * to change and the momentum is dropped. */
template <class Node>
KOKKOS_INLINE_FUNCTION static void node_motion_fold_kick(const Node &n)
{
    const double mass = n.mass();
    const double fac  = (mass > 0) ? (1.0 / mass) : 0.0;
    for(int j = 0; j < 3; j++) {n.vs(j) = (MyFloat)((double) n.vs(j) + fac * (double) n.dp(j)); n.dp(j) = 0;}
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    {
        const double l_tot = n.lum_tot(), fac_lum = (l_tot > 0) ? (1.0 / l_tot) : 0.0;
        for(int j = 0; j < 3; j++) {n.rt_vs(j) = (MyFloat)((double) n.rt_vs(j) + fac_lum * (double) n.rt_dp(j)); n.rt_dp(j) = 0;}
    }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    {
        const double mass_dm = n.mass_dm(), fac_dm = (mass_dm > 0) ? (1.0 / mass_dm) : 0.0;
        for(int j = 0; j < 3; j++) {n.vs_dm(j) = (MyFloat)((double) n.vs_dm(j) + fac_dm * (double) n.dp_dm(j)); n.dp_dm(j) = 0;}
    }
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    {
        const double sink_mass = n.sink_mass(), fac_sink = (sink_mass > 0) ? (1.0 / sink_mass) : 0.0;
        for(int j = 0; j < 3; j++) {n.sink_vel(j) = (MyFloat)((double) n.sink_vel(j) + fac_sink * (double) n.sink_dp(j)); n.sink_dp(j) = 0;}
    }
#endif
}

/* Advance every centre with its velocity over dt_drift, and widen the length by how far any member
 * can have moved over dt_widen. Every centre moves on the same clock as s: a sink centre left behind
 * would give the nearest-sink distance, the sink timestep criteria and the direct-gravity monopole
 * subtraction a position from a different time. */
template <class Node>
KOKKOS_INLINE_FUNCTION static void node_motion_advance(const Node &n, double dt_drift, double dt_widen)
{
    for(int j = 0; j < 3; j++) {
        n.s(j) = (MyFloat)((double) n.s(j) + (double) n.vs(j) * dt_drift);
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        n.rt_s(j) = (MyFloat)((double) n.rt_s(j) + (double) n.rt_vs(j) * dt_drift);
#endif
#ifdef DM_SCALARFIELD_SCREENING
        n.s_dm(j) = (MyFloat)((double) n.s_dm(j) + (double) n.vs_dm(j) * dt_drift);
#endif
#ifdef SINK_NODE_MOTION_TRACKED
        n.sink_pos(j) = (MyFloat)((double) n.sink_pos(j) + (double) n.sink_vel(j) * dt_drift);
#endif
    }
    n.len() = (MyFloat)((double) n.len() + TREE_NODE_WIDENING_DELTA(n.vmax(), dt_widen));
}

/* The gas kernel lengths a node bounds follow the flow's divergence over the drift. The scalar hmax
 * decays or grows with it (its legacy semantics); the per-type bands only ever grow here, because
 * they include sources that do not shrink under drift (e.g. a particle's force softening), so decaying
 * them could under-bound a node prune -- force_update_hmax re-grows them per particle each call. */
KOKKOS_INLINE_FUNCTION static void node_hmax_drift(struct extNODE &ext, double dt_drift_hmax)
{
    double exp_arg = (double) ext.divVmax * dt_drift_hmax / (double) NUMDIMS;
    if(exp_arg < -1.0) {exp_arg = -1.0;}
    if(exp_arg >  1.0) {exp_arg =  1.0;}
    const double decay_fac = exp(exp_arg);
    if(ext.hmax > 0) {ext.hmax = (MyFloat)((double) ext.hmax * decay_fac);}
    if(decay_fac > 1.0) {
        for(int t = 0; t < 6; t++) {
            if(ext.hmax_per_type[t] > 0) {ext.hmax_per_type[t] = (MyFloat)((double) ext.hmax_per_type[t] * decay_fac);}
        }
    }
}

#endif /* GRAVTREE_MOMENT_KERNEL_H */
