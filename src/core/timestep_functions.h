/* timestep_functions.h — Canonical KOKKOS_INLINE_FUNCTION implementations of
 * timestep utility functions.  Single source of truth for both CPU and GPU.
 *
 * timestep_dilation_factor returns the dilation factor frozen for this particle when its
 * timestep was assigned (core/timestep.cc, get_timestep), so a step's physical landing time
 * cannot mutate while it is being taken.  Identical on CPU and GPU.  The live evaluations,
 * for timestep assignment and for tree nodes, are host-only and live in core/timestep.cc.
 *
 * Include order: after allvars.h, proto.h. */
#pragma once

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
#include "../gravity/binary_functions.h"   /* binary_relative_speed_bound, for the motion bound below */
#endif

KOKKOS_INLINE_FUNCTION
double timestep_dilation_factor(int i, const struct particle_data *pp)
{
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
    if(i < 0) {return 1;}
    return pp[i].TimestepDilationFactor;
#else
    (void)i; (void)pp; return 1;
#endif
}

KOKKOS_INLINE_FUNCTION
double unit_integertime_in_physical(int i, struct particle_data *pp)
{
    return (All.Timebase_interval / All.cf_hubble_a) * timestep_dilation_factor(i, pp);
}

KOKKOS_INLINE_FUNCTION
double get_physical_timestep_from_timebin(int bin, int i, struct particle_data *pp)
{
    return GET_INTEGERTIME_FROM_TIMEBIN(bin) * unit_integertime_in_physical(i, pp);
}

KOKKOS_INLINE_FUNCTION
double get_particle_timestep_in_physical(int i, struct particle_data *pp)
{
    return pp[i].integertime_step() * unit_integertime_in_physical(i, pp);
}

/* --- live dilation factors -------------------------------------------------
 * The two nuclear-zoom helpers were file-scope statics in core/timestep.cc and read
 * nothing but All; the node factor below reads All and the node array passed to it.
 * They live here because the node dilation factor is needed wherever a drift factor
 * is, on host and device alike. core/timestep.cc provides their externally-visible
 * host symbols through its non-inline re-include of this header. */

#if defined(USE_TIMESTEP_DILATION_FOR_ZOOMS) && defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
/* smallest physical distance from pos to any of the refinement centers */
KOKKOS_INLINE_FUNCTION
double distance_to_nearest_refinement_center(Vec3<double> pos)
{
    double rmin = MAX_REAL_NUMBER;
    for(int j = 0; j < SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM; j++)
    {
        Vec3<double> p0 = All.SpecialParticle_Position_ForRefinement[j];
        Vec3<double> dp = All.cf_atime * (pos - p0);
        double r = dp.norm(); if(r < rmin) {rmin = r;}
    }
    return rmin;
}

/* dilation amplitude a >= 1 at distance r from the refinement center: unity far away, saturating
   at amax on approach */
KOKKOS_INLINE_FUNCTION
double nuclear_zoom_dilation_amplitude(double r)
{
    double fac_amax = 100.;
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 3)
    fac_amax = 1.e6;
#endif
#endif
    double amax = fac_amax;
    double r_amax = fac_amax * All.ForceSoftening[3]; // modify as needed
    double index = 1;
    if(r < 1.e-10 || isnan(r) || isfinite(r)==0) {r = 1.e-10;}
    return 1. + 1. / (1./amax + pow(r / r_amax, index));
}
#endif


/* live dilation factor at the center of mass of tree node 'no', for drifting the node itself. Nodes
   carry no particle type, so the stars-only restriction is particle-only and does not apply here.

   Nodes also carry no sink distance: Min_Distance_to_Sink is a per-particle result of the gravity
   walk, and there is no node-level equivalent to feed the weighted-motion smoothing. So under
   SPECIAL_POINT_WEIGHTED_MOTION (without the nuclear-zoom term, which does work for nodes) a node
   drifts undilated while the particles it summarises drift at the smoothing weight, leaving its
   center of mass inconsistent with them. Giving nodes that term means carrying a sink distance
   through the tree moments. The weighted-motion module is still in development; this needs
   resolving before it is relied on. */
KOKKOS_INLINE_FUNCTION
double return_node_timestep_dilation_factor_P(int no, const struct NODE *nodes)
{
#if !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS) || defined(DILATION_FOR_STELLAR_KINEMATICS_ONLY)
    (void)no; (void)nodes; return 1;
#else

    if(All.Time <= All.TimeBegin) {return 1;}
    if(no < 0) {return 1;}

    double a = 1;

#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    Vec3<double> pos_node; pos_node = nodes[no].u.d.s;
    a = nuclear_zoom_dilation_amplitude(distance_to_nearest_refinement_center(pos_node));
#else
    (void)nodes;
#endif

    return 1. / a;
#endif
}


/* --- drift and gravitational-kick time factors -----------------------------
 * A drift or kick factor is the time integral over the step, times the particle's
 * or node's timestep dilation. For a non-cosmological run the integral is just
 * the elapsed code time; for a cosmological one it is a lookup in the tables
 * built once by init_drift_table(), which is called only when
 * ComovingIntegrationOn is set -- so on a non-cosmological run the tables are
 * never filled and the pointers below are null.
 *
 * The view carries whatever a caller needs to evaluate a factor without reaching
 * a global: on the host it is filled from the tables themselves, in device code
 * from their shared-memory mirror. */
struct DriftKickTableView
{
    const double *drift;        /*!< DriftTable, or null on a non-cosmological run */
    const double *gravkick;     /*!< GravKickTable, or null on a non-cosmological run */
    double logTimeBegin;        /*!< log of the scale factor the tables start at */
    double logTimeMax;          /*!< log of the scale factor they end at */
    double timebase_interval;   /*!< code time per unit of integer time */
    int comoving;               /*!< nonzero if the tables are live and must be used */
};

/* Assembles a view. Both the host entry points and the device-side mirror build
   their view here, so the meaning of each field is fixed in one place. */
KOKKOS_INLINE_FUNCTION
struct DriftKickTableView drift_kick_table_view(const double *drift, const double *gravkick,
                                                double logTimeBegin, double logTimeMax,
                                                double timebase_interval, int comoving)
{
    struct DriftKickTableView view;
    view.drift = comoving ? drift : NULL;          /* never read when not comoving, and never built */
    view.gravkick = comoving ? gravkick : NULL;
    view.logTimeBegin = logTimeBegin;
    view.logTimeMax = logTimeMax;
    view.timebase_interval = timebase_interval;
    view.comoving = comoving;
    return view;
}

/* The host view, built from the live tables. Host code that needs a factor for an
   arbitrary interval goes through here rather than naming the four globals again. */
static inline struct DriftKickTableView drift_kick_table_view_host(void)
{
    return drift_kick_table_view(DriftTable, GravKickTable, DriftTable_logTimeBegin,
                                 DriftTable_logTimeMax, All.Timebase_interval, All.ComovingIntegrationOn);
}

/* The one interpolator. Returns the time integral over [time0, time1] with no
   dilation applied; callers multiply in the factor for the particle or node
   they are drifting. */
KOKKOS_INLINE_FUNCTION
double drift_kick_table_factor(const double *table, integertime time0, integertime time1,
                               const struct DriftKickTableView *view)
{
    if(!view->comoving) {return (time1 - time0) * view->timebase_interval;}

    double logTimeBegin = view->logTimeBegin, logTimeMax = view->logTimeMax;
    double a1 = logTimeBegin + time0 * view->timebase_interval;
    double a2 = logTimeBegin + time1 * view->timebase_interval;
    double u1, u2, df1, df2; int i1, i2;

    if(logTimeMax > logTimeBegin)
        u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
    else
        u1 = 0;
    i1 = (int) u1;
    if(i1 >= DRIFT_TABLE_LENGTH)
        i1 = DRIFT_TABLE_LENGTH - 1;

    if(i1 <= 1)
        df1 = u1 * table[0];
    else
        df1 = table[i1 - 1] + (table[i1] - table[i1 - 1]) * (u1 - i1);

    if(logTimeMax > logTimeBegin)
        u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
    else
        u2 = 0;
    i2 = (int) u2;
    if(i2 >= DRIFT_TABLE_LENGTH)
        i2 = DRIFT_TABLE_LENGTH - 1;

    if(i2 <= 1)
        df2 = u2 * table[0];
    else
        df2 = table[i2 - 1] + (table[i2] - table[i2 - 1]) * (u2 - i2);

    return df2 - df1;
}

/*! Cosmological prefactor for a drift step between time0 and time1, dilated by
 *  'dilation'. The value returned is \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a)} \f]
 */
KOKKOS_INLINE_FUNCTION
double get_drift_factor_impl(integertime time0, integertime time1, double dilation,
                             const struct DriftKickTableView *view)
{
    return drift_kick_table_factor(view->drift, time0, time1, view) * dilation;
}

KOKKOS_INLINE_FUNCTION
double get_gravkick_factor_impl(integertime time0, integertime time1, double dilation,
                                const struct DriftKickTableView *view)
{
    return drift_kick_table_factor(view->gravkick, time0, time1, view) * dilation;
}

/* --- Widening a spatial bound for the motion of what it holds ----------------
 * A node or tile records the box its members occupied when it was written, at
 * time t_ref, and the largest speed any of them can have (vmax, from
 * particle_motion_speed_bound, which already carries each member's own
 * dilation).  By `now` a member can have travelled vmax * dt on the UNDILATED
 * clock, so the box LENGTH grows by twice that -- the same rule
 * force_drift_node applies eagerly to a tree node.  Every reader of such a
 * bound widens through this one function, so a tree node opened on the device,
 * a tile opened by the BVH walk, and a node advanced on the host cannot disagree
 * about how far a box may have moved.
 *
 * Returns the growth of the LENGTH (a halfwidth grows by half of it).  Zero
 * when t_ref is not before `now`, so a fresh bound is read as written. */
KOKKOS_INLINE_FUNCTION
double motion_bound_widening(double vmax, integertime t_ref, integertime ti_now,
                             const struct DriftKickTableView *view)
{
    if(!(t_ref >= 0 && t_ref < ti_now)) {return 0.0;}   /* zero is a valid timestamp: >= 0, not > 0 */
    return TREE_DRIFT_VELOCITY_PREFAC * vmax * get_drift_factor_impl(t_ref, ti_now, 1.0, view);
}

/* A non-finite or absurd widening is a defect in the bound, never a large
 * number: falling back to the unwidened box would under-include silently, so
 * every reader tests the value and reports rather than narrows.  NaN fails
 * every comparison, hence the explicit form rather than a range check alone. */
KOKKOS_INLINE_FUNCTION
int motion_bound_widening_is_valid(double dl)
{
    return (dl >= 0.0) && (dl < 1.0e30);
}

/* The kernel radius a drift predicts follows the local compression: over an interval in
 * which the volume changes by exp(DivVel*dt), the radius changes by the NUMDIMS-th root of
 * that.  The exponent is capped so a prediction cannot move far from the last solved value,
 * which also gives anything that bounds a radius across a drift one known limit.  The
 * particle drift and the tree-node drift take their radius factor from here. */
static constexpr double KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE = 0.3;

/* The capped change in log-volume over the drift; the predicted density moves by its inverse. */
KOKKOS_INLINE_FUNCTION
double kernel_radius_drift_log_change(double divv_times_dt)
{
    if(divv_times_dt > +KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE) {return +KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE;}
    if(divv_times_dt < -KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE) {return -KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE;}
    return divv_times_dt;
}

/* The factor a radius (or a bound on radii) is multiplied by over the drift. */
KOKKOS_INLINE_FUNCTION
double kernel_radius_drift_factor(double divv_times_dt)
{
    return exp(kernel_radius_drift_log_change(divv_times_dt) / ((double)NUMDIMS));
}

/* The most one drift can grow a kernel radius by (before the radius floors are applied). */
KOKKOS_INLINE_FUNCTION
double kernel_radius_drift_max_growth_factor(void)
{
    return exp(KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE / ((double)NUMDIMS));
}


/* --- 4th-order Hermite integration -----------------------------------------
 * Which particles the Hermite integrator advances, and how a source that is not
 * itself being advanced this step is evaluated by a force walk.
 *
 * These live here rather than beside the Hermite kick routines in core/kicks.cc
 * because both gravity walks -- the host walk in gravity/forcetree.cc and the
 * device walk in gravity/gpu_gravtree.cc -- need them per interaction, and they
 * are built on the drift/kick primitives above. A cross-translation-unit call
 * per (target, source) pair costs about 16% of the walk, so the definitions are
 * inline and the walks include this header. */

#ifdef HERMITE_INTEGRATION

/*! Is particle i one the Hermite integrator advances? The bitmask says which
 *  types are eligible in principle; the tests after it drop particles whose
 *  state the scheme cannot use, either because another integrator owns their
 *  dynamics or because their history is too short or too disturbed to
 *  extrapolate from. */
KOKKOS_INLINE_FUNCTION
int eligible_for_hermite(int i, struct particle_data *pp)
{
    if(!(HERMITE_INTEGRATION & (1 << pp[i].Type))) {return 0;} // hermite flag said to not include these types
#if defined(CBE_INTEGRATOR)
    if(CBE_INTEGRATOR_DOES_TYPE(pp[i].Type)) {return 0;} // CBE moment particles: not compatible with Hermite
#elif defined(DM_FUZZY)
    if(pp[i].Type==1) {return 0;} // fuzzy-DM: not compatible with Hermite
#endif
#if defined(GRAIN_FLUID)
    if((1 << pp[i].Type) & (GRAIN_PTYPES)) {return 0;} // not compatible with these flags for these types
#endif
#if defined(SINK_PARTICLES) || defined(GALSF)
    if(pp[i].StellarAge >= DMAX(All.Time - 2*(get_particle_timestep_in_physical(i, pp)*All.cf_hubble_a), 0)) {return 0;} // if we were literally born yesterday then let things settle down a bit with the less-accurate, but more-robust regular integration
    if(pp[i].AccretedThisTimestep) {return 0;}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if(pp[i].SuperTimestepFlag >= 2) {return 0;}
#endif
    return 1;
}

/*! The per-pass state a gravity walk needs to evaluate Hermite sources, snapshotted
 *  once by the caller. HermiteOnlyFlag and TimeBinActive[] are host globals, and the
 *  active bins are carried as a bitmask so the whole thing is a few words that can be
 *  passed by value into a device kernel: no allocation and no per-particle array. */
struct HermiteWalkState
{
    integertime    ti_current;        /* All.Ti_Current at the walk */
    unsigned long long active_bins;   /* bit b set iff TimeBinActive[b] (TIMEBINS is 60) */
    int            hermite_only;      /* HermiteOnlyFlag: 0 on an ordinary leapfrog pass, 1 predictor, 2 corrector */
};

static_assert(TIMEBINS <= 64, "HermiteWalkState carries the active timebins as a 64-bit mask");

/*! The pass state the gravity walks read, refreshed once per gravity pass in gravity_tree().
 *  Host globals: the walks consume them on the host, and the device walk is handed a copy. */
extern struct HermiteWalkState   HermiteWalk;
extern struct DriftKickTableView HermiteWalkTables;   /*!< assembled only while HermiteOnlyFlag is set */

/*! Take the snapshot. Host only: it reads the host globals the walks cannot see from
 *  device code, which is the whole reason the snapshot exists. */
static inline struct HermiteWalkState hermite_walk_state_snapshot(void)
{
    struct HermiteWalkState hw;
    hw.ti_current = All.Ti_Current;
    hw.hermite_only = HermiteOnlyFlag;
    hw.active_bins = 0;
    for(int bin = 0; bin < TIMEBINS; bin++) {if(TimeBinActive[bin]) {hw.active_bins |= (1ULL << bin);}}
    return hw;
}

/*! Does source 'no' have to be re-predicted before it contributes a force?
 *
 *  Only on a Hermite pass, and only for a source the Hermite integrator owns that is
 *  NOT being advanced this step. Such a source sits at its leapfrog-drifted position
 *  and carries a whole-step-kicked velocity, both wrong at second order in the middle
 *  of its step, which caps the accuracy of the fourth-order corrector reading it.
 *
 *  A source that IS being advanced keeps its live state. On the predictor pass that
 *  state is already synchronous, and its Old* fields cannot be used anyway because
 *  find_timesteps has moved Ti_begstep up to the present, so the interval below would
 *  be empty and would return the start of its previous step. On the corrector pass
 *  do_hermite_prediction has already written the predicted state into Pos and Vel.
 *
 *  The type test comes first because it is a mask compare, while the full predicate
 *  below it reads several fields and computes a timestep. */
KOKKOS_INLINE_FUNCTION
int hermite_source_needs_prediction(int no, struct particle_data *pp, struct HermiteWalkState hw)
{
    if(!hw.hermite_only) {return 0;}
    if(!(HERMITE_INTEGRATION & (1 << pp[no].Type))) {return 0;}
    if(hw.active_bins & (1ULL << pp[no].TimeBin)) {return 0;}
    return eligible_for_hermite(no, pp);
}

/*! Position and velocity of source 'no' at the present time, extrapolated from the state
 *  it held at the start of its own step with the same third-order formula
 *  do_hermite_prediction applies to the particles being advanced. Nothing is written
 *  back: the predicted state exists only for the force evaluation asking for it. */
KOKKOS_INLINE_FUNCTION
void hermite_predict_source_state(int no, struct particle_data *pp, struct HermiteWalkState hw,
                                  const struct DriftKickTableView *tables,
                                  Vec3<double> &pos_pred, Vec3<double> &vel_pred)
{
    double dt_grav = get_gravkick_factor_impl(pp[no].Ti_begstep, hw.ti_current,
                                              timestep_dilation_factor(no, pp), tables);
    pos_pred = pp[no].OldPos + (pp[no].OldVel + (pp[no].Hermite_OldAcc + pp[no].OldJerk * (dt_grav/3)) * (dt_grav/2)) * dt_grav;
    vel_pred = pp[no].OldVel + (pp[no].Hermite_OldAcc + pp[no].OldJerk * (dt_grav/2)) * dt_grav;
}

#endif /* HERMITE_INTEGRATION */

/* The velocity a particle's POSITION advances at: a finite-volume gas cell moves with its
   mesh-generating point, everything else with its own velocity, and nothing moves at all when the
   hydro is frozen.  This is the term drift_particle_impl applies, stated once so that the speed
   bound below and anything that predicts where a particle will be cannot disagree about which
   velocity moves it.  It is the BASE term only: the bound adds a super-timestepped sink's orbital
   motion and the dilation factor on top, because those do not enter a straight-line prediction.
   Curvilinear mesh motion (HYDRO_FIX_MESH_MOTION 2/3) turns this vector as the point moves; a
   predictor may treat it as linear, since that motion assumes a radius of curvature far larger
   than either the inter-particle spacing or the distance moved in a step. */
KOKKOS_INLINE_FUNCTION
Vec3<double> particle_drift_velocity(int i, const struct particle_data *pp, const struct gas_cell_data *cell)
{
#if defined(FREEZE_HYDRO)
    (void)i; (void)pp; (void)cell;
    return Vec3<double>{0, 0, 0};
#else
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(pp[i].Type == 0) {return Vec3<double>{(double)cell[i].ParticleVel[0], (double)cell[i].ParticleVel[1], (double)cell[i].ParticleVel[2]};}
#else
    (void)cell;
#endif
    return Vec3<double>{(double)pp[i].Vel[0], (double)pp[i].Vel[1], (double)pp[i].Vel[2]};
#endif
}

/* The fastest a particle can move along any one coordinate axis, per unit of the UNDILATED drift
   interval.  A box measured when the particle was last drifted still contains it after it has
   grown by this speed times the interval since, which is what the gravity tree and the spatial
   index rely on to search among particles that have not been brought current.  Follows the
   position update in drift_particle_impl term by term -- the mesh velocity moves a finite-volume
   cell, a super-timestepped sink adds its orbital motion about the binary's centre of mass, and
   under dilation the drift covers only the dilated fraction of the interval while the nearest
   special particle's motion is added back over the rest -- so a change to how a particle moves
   belongs there and here together.  Reported per unit undilated interval so that one clock
   serves every box whatever dilation its members carry: the node's own drift may run on a
   dilated clock, but the widening that bounds its members must not. */
KOKKOS_INLINE_FUNCTION
double particle_motion_speed_bound(int i, const struct particle_data *pp, const struct gas_cell_data *cell)
{
    const Vec3<double> v_drift = particle_drift_velocity(i, pp, cell);
    double vx = v_drift[0], vy = v_drift[1], vz = v_drift[2];
    double ax = fabs(vx), ay = fabs(vy), az = fabs(vz);
    double bound = ax; if(ay > bound) {bound = ay;} if(az > bound) {bound = az;}
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3))
    /* The curvilinear mesh motion (advect_mesh_point_P) advances the radius by |v_r| dt and turns
       the point through a chord no longer than |v_t| dt; the two are not orthogonal, so the step is
       within (|v_r| + |v_t|) dt <= sqrt(2) |v| dt along every axis. */
    if(pp[i].Type == 0) {bound = 1.4142135623730951 * sqrt(vx*vx + vy*vy + vz*vz);}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0) && !defined(FREEZE_HYDRO)
    if((pp[i].Type == 5) && (pp[i].SuperTimestepFlag >= 2))
    {
        /* The binary's centre of mass drifts at its own velocity, and the sink moves about it by
           its companion's mass share of the relative motion, whose speed can only rise as far as
           the softened two-body problem allows (binary_relative_speed_bound).  The instantaneous
           relative speed would not do: the box may be left ungrown across a close passage. */
        const double Mtot = pp[i].Mass + pp[i].comp_Mass, share = pp[i].comp_Mass / Mtot;
        double cx = fabs(vx + (double)pp[i].comp_dv[0] * share), cy = fabs(vy + (double)pp[i].comp_dv[1] * share), cz = fabs(vz + (double)pp[i].comp_dv[2] * share);
        bound = cx; if(cy > bound) {bound = cy;} if(cz > bound) {bound = cz;}
        bound += share * binary_relative_speed_bound(i, pp);
    }
#endif
    const double dilation = timestep_dilation_factor(i, pp);   /* 1 unless zoom dilation is active */
    bound *= dilation;
#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    if(dilation < 1.)
    {
        double ux = fabs((double)pp[i].vel_of_nearest_special[0]), uy = fabs((double)pp[i].vel_of_nearest_special[1]), uz = fabs((double)pp[i].vel_of_nearest_special[2]);
        double u = ux; if(uy > u) {u = uy;} if(uz > u) {u = uz;}
        bound += (1. - dilation) * u;
    }
#endif
    return bound;
}

/* Where particle i will be once drifted to ti_now, for a search that reads particles which may not
   have been.  Returns the centre of that position and a half-width around it, per axis:
     CURRENT   -- drifted to ti_now already: the centre is its position, the half-width zero.
     BOUNDED   -- behind: the centre is where the drift will put it and the half-width how far it
                  can be from there, including MOTION_ENVELOPE_ROUNDING.
     UNBOUNDED -- a clock that is negative or ahead of ti_now, a result that is not finite, or a
                  centre past a special boundary the drift would move it back across: a search
                  must keep it as a candidate and leave the drift to handle or reject that state,
                  never narrow around it.
   Over a drift a particle moves in a straight line at the level the code treats it: its own velocity
   (the mesh velocity for a finite-volume cell), a super-timestepped sink's binary centre of mass,
   and under dilation the nearest special particle's motion over the undilated remainder.  These are
   the terms of the position update in drift_particle_impl, computed by the same helpers above, so an
   ordinary particle's centre is exactly where the drift puts it and a search testing it needs no
   allowance beyond rounding -- in particular none for bulk motion, which a bound built from the
   particle's speed would charge to every particle of a moving flow.  What is left carries a residual:
   a sink's orbit about its binary's centre of mass, bounded by the softened two-body speed; and a
   curvilinear mesh cell, which advect_mesh_point_P turns along an arc.  That arc ends within
   r_new |exp(i theta) - 1 - i theta| <= 2 r_new theta = 2 v_t dt of the straight line (theta = v_t dt /
   r_new), and the branches that advance straight differ from it by nothing, so 2 |v| dt covers it. */
/* The straight-line parts of how a drift over dt_drift moves particle i, shared by drift_particle_impl
   and particle_motion_envelope so that the position a search predicts and the one the drift produces
   come from the same arithmetic.  The displacement of a particle moving with its own velocity (the
   mesh velocity for a finite-volume cell; a super-timestepped sink's binary uses its centre-of-mass
   velocity below): */
KOKKOS_INLINE_FUNCTION
Vec3<double> drift_straight_displacement(int i, const struct particle_data *pp, const struct gas_cell_data *cell, double dt_drift)
{
    return particle_drift_velocity(i, pp, cell) * dt_drift;
}

#if (SINGLE_STAR_TIMESTEPPING > 0)
/* A super-timestepped sink drifts with its binary's centre of mass; its orbit about it is added apart. */
KOKKOS_INLINE_FUNCTION
Vec3<double> super_timestepped_sink_com_velocity(int i, const struct particle_data *pp)
{
    return pp[i].Vel + pp[i].comp_dv * (pp[i].comp_Mass/(pp[i].Mass+pp[i].comp_Mass));
}
#endif

/* How a drift moves particle i, on the particle's own (possibly dilated) drift clock: the velocity its
   position advances at in a straight line, and the speed of whatever part of the motion is not one -- a
   super-timestepped sink's orbit about its binary's centre of mass, and a curvilinear mesh cell's turn
   (bounded by 2|v|, see particle_motion_envelope).  The drift's position terms are these, so the
   predicted position and the drifted one agree bit for bit. */
struct particle_straight_motion {
    Vec3<double> velocity;
    double residual_speed;
};

KOKKOS_INLINE_FUNCTION
struct particle_straight_motion particle_straight_line_motion(int i, const struct particle_data *pp, const struct gas_cell_data *cell)
{
    struct particle_straight_motion m;
    m.velocity = Vec3<double>{0, 0, 0};
    m.residual_speed = 0.0;
#if !defined(FREEZE_HYDRO)
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if((pp[i].Type == 5) && (pp[i].SuperTimestepFlag >= 2))
    {
        m.velocity = super_timestepped_sink_com_velocity(i, pp);
        m.residual_speed = (pp[i].comp_Mass/(pp[i].Mass+pp[i].comp_Mass)) * binary_relative_speed_bound(i, pp);
        return m;
    }
#endif
    m.velocity = particle_drift_velocity(i, pp, cell);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3))
    if(pp[i].Type == 0) {m.residual_speed = 2.0 * m.velocity.norm();}
#endif
#else
    (void)i; (void)pp; (void)cell;
#endif
    return m;
}

/* What every drift does to a position after moving it: zero the unused dimensions, and under dilation
   add back the nearest special particle's motion over the undilated remainder of the interval (dt_drift
   carries the dilation, so only the motion relative to the surroundings is dilated). */
KOKKOS_INLINE_FUNCTION
void drift_position_finish(int i, const struct particle_data *pp, double dt_drift, Vec3<MyDouble> &pos)
{
#if (NUMDIMS==1)
    pos[1] = pos[2] = 0;
#endif
#if (NUMDIMS==2)
    pos[2] = 0;
#endif
#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    const double dilation = timestep_dilation_factor(i, pp); /* f = 1/a <= 1 */
    if(dilation < 1.) {pos += pp[i].vel_of_nearest_special * (dt_drift * (1./dilation - 1.));}
#else
    (void)i; (void)pp; (void)dt_drift;
#endif
}

/* Discovery tests a position the drift has not yet produced, possibly in different arithmetic from the
   drift and from the exact test that follows it (a device kernel may fuse a multiply and an add into one
   rounding where the host rounds twice).  Each of those few operations can move a coordinate by one
   rounding of the largest magnitude involved, so every half-width carries this many roundings of that
   scale -- outward only: it can add a candidate the exact test then rejects, never lose one. */
static constexpr double MOTION_ENVELOPE_ROUNDINGS = 16.0;

KOKKOS_INLINE_FUNCTION
double motion_envelope_rounding_floor(const double center[3], double scale)
{
    double m = fabs(scale);
    for(int k = 0; k < 3; k++) {if(fabs(center[k]) > m) {m = fabs(center[k]);}}
    return MOTION_ENVELOPE_ROUNDINGS * 2.220446049250313e-16 * m;   /* DBL_EPSILON */
}

/* The allowance a leaf test adds to a half-width for its own comparison: the separation of the query
   from the particle and its comparison with the reach are rounded at the scale of all three, which the
   exact test that decides the pair may do differently (on another device, or in another form). */
KOKKOS_INLINE_FUNCTION
double motion_envelope_test_slack(const double center[3], const double query[3], double reach)
{
    double m = fabs(reach);
    for(int k = 0; k < 3; k++) {if(fabs(query[k]) > m) {m = fabs(query[k]);}}
    return motion_envelope_rounding_floor(center, m);
}

enum particle_motion_bound_state {
    PARTICLE_MOTION_CURRENT   = 0,
    PARTICLE_MOTION_BOUNDED   = 1,
    PARTICLE_MOTION_UNBOUNDED = 2
};

KOKKOS_INLINE_FUNCTION
int particle_motion_envelope(int i, const struct particle_data *pp, const struct gas_cell_data *cell,
                             integertime ti_now, const struct DriftKickTableView *view,
                             double center[3], double *half_width)
{
    center[0] = (double)pp[i].Pos[0]; center[1] = (double)pp[i].Pos[1]; center[2] = (double)pp[i].Pos[2];
    *half_width = 0.0;
    const integertime ti_i = pp[i].Ti_current;
    if(ti_i == ti_now) {return PARTICLE_MOTION_CURRENT;}
    if(ti_i < 0 || ti_i > ti_now) {return PARTICLE_MOTION_UNBOUNDED;}
    const double dt_drift = get_drift_factor_impl(ti_i, ti_now, timestep_dilation_factor(i, pp), view);
    const Vec3<MyDouble> pos0 = pp[i].Pos;
    Vec3<MyDouble> pos = pos0;
    double residual = 0.0;
#if !defined(FREEZE_HYDRO)
    const struct particle_straight_motion motion = particle_straight_line_motion(i, pp, cell);
    pos += motion.velocity * dt_drift;
    residual = motion.residual_speed * fabs(dt_drift);
#endif
    drift_position_finish(i, pp, dt_drift, pos);
    for(int k = 0; k < 3; k++) {
        center[k] = (double)pos[k];
        if(!(center[k] - center[k] == 0.0)) {return PARTICLE_MOTION_UNBOUNDED;}   /* NaN or Inf, fast-math safe */
    }
    if(!(residual >= 0.0 && residual < 1.0e30)) {return PARTICLE_MOTION_UNBOUNDED;}
    double moved = 0.0;
    for(int k = 0; k < 3; k++) {const double d = fabs(center[k] - (double)pos0[k]); if(d > moved) {moved = d;}}
    const double width = residual + motion_envelope_rounding_floor(center, moved);
#if BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE
    {   /* the sides the drift acts on: reflect or outflow, lower (code 0 or -1) and upper (0 or 1) */
        const double box_upper[3] = {boxSize_X, boxSize_Y, boxSize_Z};
        for(int k = 0; k < NUMDIMS; k++) {
            const int rf = special_boundary_condition_xyz_def_reflect[k], of = special_boundary_condition_xyz_def_outflow[k];
            const int lower = (rf == 0 || rf == -1 || of == 0 || of == -1), upper = (rf == 0 || rf == 1 || of == 0 || of == 1);
            if((lower && center[k] <= width) || (upper && center[k] >= box_upper[k] - width)) {return PARTICLE_MOTION_UNBOUNDED;}
        }
    }
#endif
    *half_width = width;
    return PARTICLE_MOTION_BOUNDED;
}

/* How far a drift carries particle i per unit of the UNDILATED drift interval, so that one clock serves
   every member of a spatial index whatever dilation each carries: over an undilated interval D the drift
   moves it by u D in a straight line, to within residual_speed D.  These are particle_straight_line_motion's
   terms put on the undilated clock -- the drift covers only the dilated fraction f of the interval, and
   under DILATION_FOR_STELLAR_KINEMATICS_ONLY adds back the nearest special particle's motion over the rest
   (drift_position_finish) -- so a position the drift produced at one time, moved by u D, is where the drift
   puts the particle a drift interval D later.  u changes only when the particle is kicked or its motion is
   written directly; a drift leaves it alone. */
KOKKOS_INLINE_FUNCTION
void particle_transport_velocity(int i, const struct particle_data *pp, const struct gas_cell_data *cell,
                                 double u[3], double *residual_speed)
{
    const double f = timestep_dilation_factor(i, pp);
    u[0] = u[1] = u[2] = 0.0;
    *residual_speed = 0.0;
#if !defined(FREEZE_HYDRO)
    const struct particle_straight_motion m = particle_straight_line_motion(i, pp, cell);
    for(int k = 0; k < 3; k++) {u[k] = f * m.velocity[k];}
    *residual_speed = f * m.residual_speed;
#else
    (void)cell;
#endif
#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    if(f < 1.) {for(int k = 0; k < 3; k++) {u[k] += (1. - f) * (double)pp[i].vel_of_nearest_special[k];}}
#endif
#if (NUMDIMS == 1)
    u[1] = u[2] = 0.0;
#endif
#if (NUMDIMS == 2)
    u[2] = 0.0;
#endif
}
