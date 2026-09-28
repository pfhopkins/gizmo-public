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
    double vx = 0, vy = 0, vz = 0;
#if !defined(FREEZE_HYDRO)
    vx = (double)pp[i].Vel[0]; vy = (double)pp[i].Vel[1]; vz = (double)pp[i].Vel[2];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(pp[i].Type == 0) {vx = (double)cell[i].ParticleVel[0]; vy = (double)cell[i].ParticleVel[1]; vz = (double)cell[i].ParticleVel[2];}
#else
    (void)cell;
#endif
#else
    (void)cell;
#endif
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
