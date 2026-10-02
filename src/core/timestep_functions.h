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
double node_timestep_dilation_factor_at(const Vec3<double> &pos_node)
{
#if !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS) || defined(DILATION_FOR_STELLAR_KINEMATICS_ONLY)
    (void)pos_node; return 1;
#else

    if(All.Time <= All.TimeBegin) {return 1;}

    double a = 1;

#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    a = nuclear_zoom_dilation_amplitude(distance_to_nearest_refinement_center(pos_node));
#else
    (void)pos_node;
#endif

    return 1. / a;
#endif
}

/* The same, at the node's centre of mass as the tree array holds it.  A walk reading its own
   copy of a node evaluates node_timestep_dilation_factor_at on that copy's centre instead. */
KOKKOS_INLINE_FUNCTION
double return_node_timestep_dilation_factor_P(int no, const struct NODE *nodes)
{
#if !defined(USE_TIMESTEP_DILATION_FOR_ZOOMS) || defined(DILATION_FOR_STELLAR_KINEMATICS_ONLY)
    (void)no; (void)nodes; return 1;
#else
    if(no < 0) {return 1;}
    Vec3<double> pos_node = {};
#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    pos_node = nodes[no].u.d.s;
#else
    (void)nodes;
#endif
    return node_timestep_dilation_factor_at(pos_node);
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

/* The two intervals a tree node moves over between time0 and time1: its centres on its own
   (dilated) clock, and its length widened on the undilated clock, because vmax already carries each
   member's own dilation.  One home for both on the device: the full sweep, the claimed-subset drift
   and the gravity walk's read-only prediction.  The host lazy drift (force_drift_node) computes the
   same two factors itself, through the host tables. */
KOKKOS_INLINE_FUNCTION
void node_motion_intervals(integertime time0, integertime time1, double dilation,
                           const struct DriftKickTableView *view, double &dt_drift, double &dt_widen)
{
    dt_drift = get_drift_factor_impl(time0, time1, dilation, view);
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
    dt_widen = get_drift_factor_impl(time0, time1, 1.0, view);
#else
    dt_widen = dt_drift;
#endif
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

#if (SINGLE_STAR_TIMESTEPPING > 0)
/* A sink whose binary is integrated over its own super-timestep: it drifts with the binary's centre
   of mass, and its orbit about it is advanced separately. */
KOKKOS_INLINE_FUNCTION int is_super_timestepped_sink(int i, const struct particle_data *pp)
{
    return ((pp[i].Type == 5) && (pp[i].SuperTimestepFlag >= 2)) ? 1 : 0;
}
#endif

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
    if(is_super_timestepped_sink(i, pp))
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


/* ==========================================================================================
 * Particle motion over a drift.  The drift (core/drift_particle_functions.h) and the read-only
 * prediction below run the same bodies; particle_motion_speed_bound above bounds this motion.
 * ========================================================================================== */
/* Where a particle's motion is read and written by the step bodies below: the mesh-point advection,
 * the special boundaries and the position step.  One body serves two
 * callers.  The drift's accessor, particle_motion_in_arrays, reaches into the particle's own arrays,
 * and its hooks carry out the step's writes to other fields (the predicted gas velocity, the momentum
 * ledger, the reflected fluxes, the cell mass).  A read-only prediction, particle_motion_prediction,
 * keeps its own copies of position, velocity, mass and mesh velocity, and its hooks do nothing, so it
 * writes nothing.  Every hook runs at the point in the body where the original write stood.  Fields a
 * step reads but never changes come from stored() in both. */
struct particle_motion_in_arrays {
    struct particle_data *pp; struct gas_cell_data *cell; int i;
    static constexpr int fewbody_mode = 1;   /* the binary's internal state is advanced and kept */
    KOKKOS_INLINE_FUNCTION struct particle_data *particles() const {return pp;}
    KOKKOS_INLINE_FUNCTION const struct particle_data &stored() const {return pp[i];}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Pos)  &pos()  const {return pp[i].Pos;}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Vel)  &vel()  const {return pp[i].Vel;}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Mass) &mass() const {return pp[i].Mass;}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    KOKKOS_INLINE_FUNCTION decltype(gas_cell_data::ParticleVel) &mesh_vel() const {return cell[i].ParticleVel;}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    KOKKOS_INLINE_FUNCTION void fewbody_use_com_acceleration() const {pp[i].GravAccel = pp[i].COM_GravAccel;}
#endif
    KOKKOS_INLINE_FUNCTION void velocity_reflected(int j) const {if(pp[i].Type==0) {cell[i].VelPred[j]=pp[i].Vel[j]; cell[i].HydroAccel[j]=0;}}
    KOKKOS_INLINE_FUNCTION void add_reflected_momentum(int j, double mass_for_dp) const {pp[i].dp[j]+=2*pp[i].Vel[j]*mass_for_dp;}
    KOKKOS_INLINE_FUNCTION void zero_momentum() const {pp[i].dp[0]=pp[i].dp[1]=pp[i].dp[2]=0;}
    KOKKOS_INLINE_FUNCTION void outflow_removed_mass() const {if(pp[i].Type==0) {cell[i].Mass=0;}}
    KOKKOS_INLINE_FUNCTION void reflect_fluxes_at_lower_face(int j) const
    {
        (void)j;
#ifdef RT_EVOLVE_FLUX
        if(pp[i].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(cell[i].Rad_Flux[kf][j]<0) {cell[i].Rad_Flux[kf][j]=-cell[i].Rad_Flux[kf][j]; cell[i].Rad_Flux_Pred[kf][j]=cell[i].Rad_Flux[kf][j];}}}
#endif
#ifdef COSMIC_RAY_FLUID
        if(pp[i].Type==0) {int kf; for(kf=0;kf<N_CR_PARTICLE_BINS;kf++) {if(cell[i].CosmicRayFlux[kf][j]<0) {cell[i].CosmicRayFlux[kf][j]=-cell[i].CosmicRayFlux[kf][j]; cell[i].CosmicRayFluxPred[kf][j]=cell[i].CosmicRayFlux[kf][j];}}}
#endif
    }
    KOKKOS_INLINE_FUNCTION void reflect_fluxes_at_upper_face(int j) const
    {
        (void)j;
#ifdef RT_EVOLVE_FLUX
        if(pp[i].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(cell[i].Rad_Flux[kf][j]>0) {cell[i].Rad_Flux[kf][j]=-cell[i].Rad_Flux[kf][j]; cell[i].Rad_Flux_Pred[kf][j]=cell[i].Rad_Flux[kf][j];}}}
#endif
#ifdef COSMIC_RAY_FLUID
        if(pp[i].Type==0) {int kf; for(kf=0;kf<N_CR_PARTICLE_BINS;kf++) {if(cell[i].CosmicRayFlux[kf][j]>0) {cell[i].CosmicRayFlux[kf][j]=-cell[i].CosmicRayFlux[kf][j]; cell[i].CosmicRayFluxPred[kf][j]=cell[i].CosmicRayFlux[kf][j];}}}
#endif
    }
};

struct particle_motion_prediction {
    struct particle_data *pp; struct gas_cell_data *cell; int i;
    decltype(particle_data::Pos)  pos_;
    decltype(particle_data::Vel)  vel_;
    decltype(particle_data::Mass) mass_;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    decltype(gas_cell_data::ParticleVel) mesh_vel_;
#endif
    static constexpr int fewbody_mode = 0;   /* the binary's internal state is read, not advanced */
    /* Starts from the stored state; the cell is read only for a gas particle. */
    KOKKOS_INLINE_FUNCTION particle_motion_prediction(struct particle_data *pp_in, struct gas_cell_data *cell_in, int i_in)
        : pp(pp_in), cell(cell_in), i(i_in), pos_(pp_in[i_in].Pos), vel_(pp_in[i_in].Vel), mass_(pp_in[i_in].Mass)
    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        mesh_vel_ = {}; if(pp_in[i_in].Type == 0) {mesh_vel_ = cell_in[i_in].ParticleVel;}
#endif
    }
    KOKKOS_INLINE_FUNCTION struct particle_data *particles() const {return pp;}
    KOKKOS_INLINE_FUNCTION const struct particle_data &stored() const {return pp[i];}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Pos)  &pos()  {return pos_;}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Vel)  &vel()  {return vel_;}
    KOKKOS_INLINE_FUNCTION decltype(particle_data::Mass) &mass() {return mass_;}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    KOKKOS_INLINE_FUNCTION decltype(gas_cell_data::ParticleVel) &mesh_vel() {return mesh_vel_;}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    KOKKOS_INLINE_FUNCTION void fewbody_use_com_acceleration() const {}
#endif
    KOKKOS_INLINE_FUNCTION void velocity_reflected(int) const {}
    KOKKOS_INLINE_FUNCTION void add_reflected_momentum(int, double) const {}
    KOKKOS_INLINE_FUNCTION void zero_momentum() const {}
    KOKKOS_INLINE_FUNCTION void outflow_removed_mass() const {}
    KOKKOS_INLINE_FUNCTION void reflect_fluxes_at_lower_face(int) const {}
    KOKKOS_INLINE_FUNCTION void reflect_fluxes_at_upper_face(int) const {}
};


#ifdef HYDRO_MESHLESS_FINITE_VOLUME
template <class Motion>
KOKKOS_INLINE_FUNCTION
void advect_mesh_point_body(Motion &a, double dt)
{
#if (HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3) // cylindrical or spherical coordinates
    /* The mesh velocity is read as a rotation about the centre plus a radial (and, for the
       cylindrical case, a vertical) motion, so that a rigidly rotating mesh stays rigid: the
       point's direction from the centre is turned in the plane it and its tangential velocity
       span, through the angle that velocity sweeps in dt, and its velocity is turned with it. */
    Vec3<double> dp = a.pos(); Vec3<double> dp_offset = {}; // location relative to the centre; the centre is the coordinate origin ...
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) // ... unless a special anchor defines it
    dp_offset = -(a.pos() + a.stored().Min_xyz_to_Sink);   /* Min_xyz_to_Sink = x_sink - x, so dp = x - x_sink */
#elif defined(BOX_PERIODIC) // ... or the box is periodic, when it is the box mid-point
#if (NUMDIMS==1)
    dp_offset[0] = -boxHalf_X;
#elif (NUMDIMS==2)
    dp_offset[0] = -boxHalf_X; dp_offset[1] = -boxHalf_Y;
#else
    dp_offset = Vec3<double>{-boxHalf_X, -boxHalf_Y, -boxHalf_Z};
#endif
#endif
    dp += dp_offset;
    Vec3<double> v = a.mesh_vel();
    Vec3<double> radial = dp;                 // the part of the position that rotates
#if (HYDRO_FIX_MESH_MOTION == 2)
    radial[2] = 0;                            // cylindrical: z advances on its own
#endif
    const double r = radial.norm();
    if(r > 0)
    {
        const Vec3<double> e_r = radial / r;
        const double v_r = dot(v, e_r);
        Vec3<double> v_t = v - v_r * e_r;     // tangential velocity, in the plane of rotation
#if (HYDRO_FIX_MESH_MOTION == 2)
        v_t[2] = 0;
#endif
        const double vt = v_t.norm();
        const double r_new = r + v_r * dt;
        Vec3<double> e_r_new = e_r, e_t_new = {};
        if(r_new <= 0) {a.pos() += v * dt; return;}    /* carried through the axis in one step: the turn is undefined, so advance straight */
        if(vt > 0)
        {
            /* The angle is the tangential distance over the radius the point ends up at, so the
               chord it turns through is no longer than v_t dt, beside a radial advance of
               |v_r| dt: a displacement of at most (|v_r| + |v_t|) dt (particle_motion_speed_bound). */
            const Vec3<double> e_t = v_t / vt;
            const double angle = vt * dt / r_new, c = cos(angle), s = sin(angle);
            e_r_new = c * e_r + s * e_t;      // the direction turned through the swept angle
            e_t_new = c * e_t - s * e_r;      // and the tangential direction with it
        }
        dp = r_new * e_r_new;
        v  = v_r * e_r_new + vt * e_t_new;    // the same speed, turned with the point
#if (HYDRO_FIX_MESH_MOTION == 2)
        dp[2] = a.pos()[2] + dp_offset[2] + a.mesh_vel()[2] * dt;
        v[2]  = a.mesh_vel()[2];
#endif
        a.mesh_vel() = v;
    }
    else {dp += v * dt;}                      // a point at the centre has no direction to turn
    a.pos() = dp - dp_offset;                 // back to the simulation frame
    return;
#endif // ok done with cylindrical/spherical coordinates


    // ok anything else ('normal' coordinates), does down here
    a.pos() += a.mesh_vel() * dt; // for standard grid velocities, this is trivial //
    return;
}

KOKKOS_INLINE_FUNCTION
void advect_mesh_point_P(int i, double dt, struct particle_data *pp, struct gas_cell_data *cell)
{
    particle_motion_in_arrays a = {pp, cell, i};
    advect_mesh_point_body(a, dt);
}

/* A finite-volume cell's mass after a drift of dt_entr: it follows the mass flux, but never falls
   below half of its conserved mass. */
KOKKOS_INLINE_FUNCTION
double mfv_drifted_mass(double mass, double dt_mass, double mass_true, double dt_entr)
{
    return DMAX(mass + dt_mass * dt_entr, 0.5 * mass_true);
}
#endif /* HYDRO_MESHLESS_FINITE_VOLUME */

/* mass_for_dp is the caller's: the kick passes the mass it has just assigned, the drift the
   particle's current mass. */
template <class Motion>
KOKKOS_INLINE_FUNCTION
void apply_special_boundary_conditions_body(Motion &a, double mass_for_dp, int mode)
{
#if BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE
    double box_upper[3]; int j;
    box_upper[0]=boxSize_X; box_upper[1]=boxSize_Y; box_upper[2]=boxSize_Z;
    for(j=0; j<3; j++)
    {
        if(a.pos()[j] <= 0)
        {
            if(special_boundary_condition_xyz_def_reflect[j] == 0 || special_boundary_condition_xyz_def_reflect[j] == -1)
            {
                if(a.vel()[j]<0) {a.vel()[j]=-a.vel()[j]; a.velocity_reflected(j); if(mode==1) {a.add_reflected_momentum(j, mass_for_dp);}}
                a.pos()[j]=DMAX((0.+((double)a.stored().ID)*2.e-8)*box_upper[j], 0.1*a.pos()[j]); // old  was 1e-9, safer on some problems, but can artificially lead to 'trapping' in some low-res tests
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                a.pos()[j]+=3.e-3*boxSize_X; a.vel()[j] += 0.1; /* special because of our wierd boundary condition for this problem, sorry to have so many hacks for this! */
#endif
                a.reflect_fluxes_at_lower_face(j);
            }
            if(special_boundary_condition_xyz_def_outflow[j] == 0 || special_boundary_condition_xyz_def_outflow[j] == -1) {a.mass()=0; a.outflow_removed_mass(); if(mode==1) {a.zero_momentum();}}
        }
        else if (a.pos()[j] >= box_upper[j])
        {
            if(special_boundary_condition_xyz_def_reflect[j] == 0 || special_boundary_condition_xyz_def_reflect[j] == 1)
            {
                if(a.vel()[j]>0) {a.vel()[j]=-a.vel()[j]; a.velocity_reflected(j); if(mode==1) {a.add_reflected_momentum(j, mass_for_dp);}}
                a.pos()[j]=box_upper[j]*(1.-((double)a.stored().ID)*2.e-8);
                a.reflect_fluxes_at_upper_face(j);
            }
            if(special_boundary_condition_xyz_def_outflow[j] == 0 || special_boundary_condition_xyz_def_outflow[j] == 1) {a.mass()=0; a.outflow_removed_mass(); if(mode==1) {a.zero_momentum();}}
        }
    }
#else
    (void)a; (void)mass_for_dp; (void)mode;
#endif
    return;
}

KOKKOS_INLINE_FUNCTION
void apply_special_boundary_conditions_P(int i, double mass_for_dp, int mode, struct particle_data *pp, struct gas_cell_data *cell)
{
    particle_motion_in_arrays a = {pp, cell, i};
    apply_special_boundary_conditions_body(a, mass_for_dp, mode);
}

/* The position a drift of dt_drift gives a particle: its binary's motion for a super-timestepped
   sink, its mesh-generating point's for a finite-volume gas cell, its own velocity otherwise; then
   the unused axes are zeroed and, under DILATION_FOR_STELLAR_KINEMATICS_ONLY, the bulk motion over
   the undilated remainder of the interval is added back.  One body for the drift and for a
   read-only prediction (see particle_motion_in_arrays above).  fewbody_kick_dv returns the
   binary's velocity kick, which the drift also applies to the gas velocity it predicts. */
template <class Motion>
KOKKOS_INLINE_FUNCTION
void particle_position_step(Motion &a, double dt_drift, Vec3<double> &fewbody_kick_dv)
{
    (void)fewbody_kick_dv;
#if !defined(FREEZE_HYDRO)
    /* A finite-volume gas cell moves with its mesh-generating point, a super-timestepped sink
       with its binary, everything else with its own velocity.  The two special cases are
       different particle types, so both must be live in a build that has both. */
#if (SINGLE_STAR_TIMESTEPPING > 0)
    volatile int super_timestepped_sink = is_super_timestepped_sink(a.i, a.particles());
    if(super_timestepped_sink)
    {
        /* The orbit integration reads the binary's stored state, so the centre-of-mass velocity does too. */
        Vec3<double> fewbody_drift_dx = {};
        Vec3<double> COM_Vel = a.vel() + a.stored().comp_dv * (a.stored().comp_Mass/(a.stored().Mass+a.stored().comp_Mass)); //center of mass velocity
        a.pos() += COM_Vel * dt_drift; //center of mass drift
        odeint_super_timestep(a.i, dt_drift, fewbody_kick_dv, fewbody_drift_dx, Motion::fewbody_mode, a.particles()); // do_fewbody_drift
        a.fewbody_use_com_acceleration(); //Overwrite the acceleration with center of mass value
        a.pos() += fewbody_drift_dx; //Keplerian evolution
        a.vel() += fewbody_kick_dv; //move on binary.orbit
    }
    else
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(a.stored().Type==0) {advect_mesh_point_body(a, dt_drift);}
    else
#endif
    {a.pos() += a.vel() * dt_drift;}
#endif // FREEZE_HYDRO clause
#if (NUMDIMS==1)
    a.pos()[1]=a.pos()[2]=0; // force zero-ing
#endif
#if (NUMDIMS==2)
    a.pos()[2]=0; // force zero-ing
#endif

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    double dilation = timestep_dilation_factor(a.i, a.particles()); /* f = 1/a <= 1 */
    if(dilation < 1.) {
        /* the drift above advanced the particle over only the fraction f of the raw interval, since
           dt_drift already carries the f. add back the bulk motion over the remaining (1-f) of the
           raw interval, so that only the motion relative to the surroundings is dilated */
        double cfac = dt_drift * (1./dilation - 1.);
        a.pos() += a.stored().vel_of_nearest_special * cfac;
    }
#endif
}

/* A particle's position, velocity and mass at ti_to, predicted from its stored state (m starts
   there) by the same step bodies the drift runs -- position, finite-volume mass, special
   boundaries, in the drift's order -- with nothing written back.  Nothing else a drift changes is
   predicted: softening, kernel lengths, density, energy, pressure and the gas velocity stay as
   stored.  A particle already at or past ti_to is left as stored; the caller rejects a source
   ahead of its walk time.  A Hermite source's position and velocity are the caller's to replace
   afterwards (hermite_predict_source_state); a super-timestepped sink is never a Hermite source.
   drift_particle_impl runs these steps in this order with other physics between them (the mass step
   sits in its gas block, the boundaries at its end); a change to one composition belongs in both. */
KOKKOS_INLINE_FUNCTION
void predict_particle_motion(struct particle_motion_prediction &m, integertime ti_to, const struct DriftKickTableView *tables)
{
    const int i = m.i; struct particle_data *pp = m.pp; struct gas_cell_data *cell = m.cell;
    const integertime ti_from = pp[i].Ti_current;
    if(ti_to <= ti_from) {return;}
    const double dt_drift = get_drift_factor_impl(ti_from, ti_to, timestep_dilation_factor(i, pp), tables);
    Vec3<double> fewbody_kick_dv = {};
    particle_position_step(m, dt_drift, fewbody_kick_dv);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    if((pp[i].Type == 0) && (m.mass() > 0)) {m.mass() = mfv_drifted_mass(m.mass(), cell[i].DtMass, cell[i].MassTrue, (ti_to - ti_from) * unit_integertime_in_physical(i, pp));}
#else
    (void)cell;
#endif
    apply_special_boundary_conditions_body(m, m.mass(), 0);
}
