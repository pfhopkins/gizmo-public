#ifndef GRAVTREE_FORCE_KERNEL_H
#define GRAVTREE_FORCE_KERNEL_H

/* Accepted-source contribution physics for the gravity tree walk — the single home for the
 * pair force/potential/tidal kernel and the per-payload accumulation formulas shared by the
 * CPU walk (force_treeevaluate, forcetree.cc) and the GPU walk (gpu_gravtree_walk_one,
 * gpu_gravtree.cc). Companion to gravtree_opening.h (the node-opening DECISION); this file
 * owns the CONTRIBUTION between an accepted source and the target.
 *
 * Contract: helpers take plain scalars / small POD structs. The WALKS own loading (each from
 * its own storage: P[]/Nodes[]/Extnodes[] vs P_dev/SoA) and storing; helpers own only the
 * physics/numerics between. Helpers are epoch-agnostic — they consume a payload snapshot and
 * do not define when it refreshes.
 *
 * Include AFTER declarations/allvars.h and mesh/kernel.h (kernel_gravity/kernel_main and the
 * All.* / unit macros are consumed from the including TU, mirroring gravtree_opening.h).
 *
 * Adding a payload — checklist (every step in the SAME commit):
 *   1. target prologue: any per-target derived quantity -> a grav_target_* helper here.
 *   2. source load: walk-side, leaf (P[]/P_dev) AND node (Nodes/Extnodes vs SoA); if the load
 *      embeds a formula (not a raw field read), the formula becomes a helper here.
 *   3. accumulate: a grav_<payload>_accumulate helper here, called by BOTH walks at the same
 *      structural point (inside/outside the PM short-range gate exactly as the physics needs).
 *   4. store: walk-side (CPU walk writes P[] directly; GPU walk writes P_dev/CellP_dev + host scatter).
 *   5. wire/SoA: extend the LET wire + tree SoA + gpu_scatter_foreign_to_soa + drift kernel.
 *   6. moment construction: gpu_moment_refresh + force_treeupdate_pseudos (+ LET synth leaf).
 */

#include "gravity_box_distance.h"
#include "let_data.h"            /* LET_LEAF_TAG_* + grav_classify_node (import topology vocabulary) */


#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

/* Length of the TreePM short-range truncation tables (forcetree.cc NTAB). */
#define GRAVTREE_SHORTRANGE_NTAB 1000


/* ------------------------------------------------------------------------------------------
 * Shared source/target formula helpers (single home for bodies that previously existed as
 * host functions + device replicas).
 * ---------------------------------------------------------------------------------------- */

/* Which particle types a particle of type 'particle_type_primary' shares an adaptive-softening
 * kernel with (body moved verbatim from ags_rkern.cc; the host ags_gravity_kernel_shared_BITFLAG
 * delegates here). Returns SUM(2^n) over the types the primary 'sees'. */
KOKKOS_INLINE_FUNCTION int gravtree_ags_kernel_shared_bitflag(int particle_type_primary)
{
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(!((1 << particle_type_primary) & (ADAPTIVE_GRAVSOFT_FORALL))) {return 0;} /* particle is NOT one of the designated 'adaptive' types */
#endif

#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    if(!((1 << particle_type_primary) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION))) {return ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION;} /* particle is NOT one of the designated 'adaptive' types */
#endif

    if(particle_type_primary == 0) {return 1;} /* gas particles see gas particles */

#if (ADAPTIVE_GRAVSOFT_FORALL & 32) && defined(SINK_PARTICLES)
    if(particle_type_primary == 5) {return 1;} /* sink particle particles are AGS-active, but using sink physics, they see only gas */
#endif

#if defined(GALSF) && ( (ADAPTIVE_GRAVSOFT_FORALL & 16) || (ADAPTIVE_GRAVSOFT_FORALL & 8) || (ADAPTIVE_GRAVSOFT_FORALL & 4) )
    if(All.ComovingIntegrationOn) /* stars [4 for cosmo runs, 2+3+4 for non-cosmo runs] are AGS-active and see baryons (any type) */
    {
        if(particle_type_primary == 4) {return 17;} // 2^0+2^4
    } else {
        if((particle_type_primary == 4)||(particle_type_primary == 2)||(particle_type_primary == 3)) {return 29;} // 2^0+2^2+2^3+2^4
    }
#endif

#ifdef DM_SIDM
    if((1 << particle_type_primary) & (DM_SIDM)) {return DM_SIDM;} /* SIDM particles see other SIDM particles, regardless of type/mass */
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    return (1 << particle_type_primary); /* if we haven't been caught by one of the above checks, we simply return whether or not we see 'ourselves' */
#endif

    return 0;
}

/* Angle-weighted long-range sink luminosity (body moved verbatim from sinks/sink.cc; the host
 * sink_fb_angleweight delegates here). Component args so both walks' angle-vector storage types feed it. */
KOKKOS_INLINE_FUNCTION double grav_sink_fb_angleweight(double sink_lum_input, double sa_x, double sa_y, double sa_z, double dx, double dy, double dz)
{
#ifdef SINGLE_STAR_SINK_DYNAMICS
    return sink_lum_input;
#else
    if(sink_lum_input <= 0) return 0;
    Vec3<double> d{dx, dy, dz};
    double r2 = d.norm_sq(); if(r2 <= 0) return 0;
    if(r2*UNIT_LENGTH_IN_PC*UNIT_LENGTH_IN_PC*All.cf_atime*All.cf_atime < 1) return 0; /* no force at < 1pc */
#if defined(SINK_FB_COLLIMATED)
    Vec3<double> sa{sa_x, sa_y, sa_z};
    double cos_theta = fabs(dot(d, sa)/sqrt(r2*sa.norm_sq())); if(!isfinite(cos_theta)) {cos_theta=1;}
    double wt_normalized = 0.0847655*exp(4.5*cos_theta*cos_theta); // ~exp(-x^2/2*hR^2), normalized appropriately to give the correct total flux, for hR~0.3
    return sink_lum_input * wt_normalized;
#else
    return sink_lum_input;
#endif
#endif /* SINGLE_STAR_SINK_DYNAMICS */
}

#ifdef SPECIAL_POINT_WEIGHTED_MOTION
/* Weight function for the special-point motion smoothing (body moved verbatim from
 * core/predict.cc; the host weight_function_for_weighted_motion_smoothing delegates here). */
KOKKOS_INLINE_FUNCTION double grav_weight_function_for_weighted_motion_smoothing(double r, int mode)
{
    double wt = 1, amax = 1.e2, rmax_pc = 100., slope_index = 1, r_pc = r * All.cf_atime * UNIT_LENGTH_IN_PC;
    if(r_pc < rmax_pc) {wt = DMAX(pow(r_pc / rmax_pc , slope_index) , 1./amax);}
    if(mode) {return 1 - sqrt(wt);} else {return wt;}
}
#endif


/* ------------------------------------------------------------------------------------------
 * Target prologue helpers — per-target derived quantities. Caller loads raw fields.
 * ---------------------------------------------------------------------------------------- */

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
/* Adaptive-softening target selection: keep the adaptive kernel radius + zeta correction only
 * when it exceeds the type floor, else fall back to the fixed softening with zeta=0. 'soft'
 * enters as the walk's loaded kernel radius and may be reset here. */
KOKKOS_INLINE_FUNCTION void grav_target_select_soft_and_zeta(int ptype, double ags_zeta, double &soft, double &zeta)
{
    zeta = 0;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
    if(ptype == 0) {if(soft > All.ForceSoftening[ptype]) {zeta = ags_zeta;} else {soft=All.ForceSoftening[ptype]; zeta=0;}}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    if(soft > All.ForceSoftening[ptype]) {zeta = ags_zeta;} else {soft=All.ForceSoftening[ptype]; zeta=0;}
#endif
}
#endif

#ifdef RT_USE_GRAVTREE
/* Valid-gas gate for every RT payload load/accumulate (both walks store the result; the GPU
 * walk keeps it in a volatile int — nvc++ constant-propagates plain int gates in device code). */
KOKKOS_INLINE_FUNCTION int grav_target_valid_gas_for_rt(int ptype, double soft, double pmass)
{
    return (ptype==0 && soft>0 && pmass>0) ? 1 : 0;
}
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
/* RT_LEBRON on-the-spot radiation-pressure coupling factor, once per target. kappa_eff[kf] is
 * rt_kappa(-1,kf,...) evaluated caller-side (rt_kappa needs the walk's particle pointers). */
KOKKOS_INLINE_FUNCTION void grav_target_rt_fac_stellum(double soft, double pmass, const double *kappa_eff, double *fac_stellum)
{
    double h_eff_phys = soft * pow(VOLUME_NORM_COEFF_FOR_NDIMS/All.DesNumNgb,1./NUMDIMS) * All.cf_atime; // convert from softening kernel extent to effective size, assuming 3D here, and convert to physical code units
    double sigma_particle =  pmass / (h_eff_phys*h_eff_phys); // quick estimate of effective surface density of the target, in physical code units
    double fac_stellum_0 = -All.PhotonMomentum_Coupled_Fraction / (4.*M_PI * C_LIGHT_CODE_REDUCED * sigma_particle * All.G); // this will be multiplied by L/r^2 below, giving acceleration, then extra G because code thinks this is gravity, so put extra G here. everything is in -physical- code units here //
    int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fac_stellum[kf] = fac_stellum_0*(1 - exp(-kappa_eff[kf]*sigma_particle));} // kappa is in physical code units, so sigma_eff_abs should be also -- approximate surface-density through particle (for checking if we enter optically-thick limit)
}
#endif
#endif /* RT_USE_GRAVTREE */

#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
/* Enclosed-mass count radius: target softening with a fixed floor (otherwise the statistics are very noisy). */
KOKKOS_INLINE_FUNCTION double grav_target_menc_radius(double soft)
{
    double r_for_total_menclosed = soft; r_for_total_menclosed = DMAX( r_for_total_menclosed , 0.1/(UNIT_LENGTH_IN_KPC*All.cf_atime) );
    return r_for_total_menclosed;
}
#endif

#ifdef PMGRID
/* TreePM short-range table scaling, once per target (rcut/asmth may be per-target under PM_PLACEHIGHRESREGION). */
KOKKOS_INLINE_FUNCTION double grav_pm_asmthfac(double asmth)
{
    return 0.5 / asmth * (GRAVTREE_SHORTRANGE_NTAB / 3.0);
}
#endif


/* ------------------------------------------------------------------------------------------
 * Core pair kernel: Newtonian / softened force selection, BOTH softening symmetrizations,
 * AGS zeta corrections. Returns the scalar prefactors; the walk owns acc/pot accumulation.
 * Body moved verbatim from force_treeevaluate (forcetree.cc).
 * ---------------------------------------------------------------------------------------- */

struct grav_force_pair_t {
    double fac_accel;
#ifdef EVALPOTENTIAL
    double fac_pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    double fac_tidal, fac2_tidal;
#endif
};

/* PM short-range truncation CONFIG (read-only, per-walk-constant). Unconditional POD
 * passed BY VALUE into the walk and the PM-gated force helpers; its fields are gated
 * once here and it is empty when !PMGRID. Holds configuration only -- the per-pair
 * tabindex and the per-pair table lookups stay in the walk / helper bodies. Both the
 * CPU and GPU walks build one (PM_PLACEHIGHRESREGION overrides the per-target copy's
 * rcut/rcut2/asmthfac). */
struct grav_pm_shortrange_t {
#ifdef PMGRID
    double rcut, rcut2, asmthfac;
    const float *shortrange_tab;
#ifdef EVALPOTENTIAL
    const float *shortrange_pot_tab;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    const float *shortrange_tidal_tab;
#endif
#endif
};

/* Foreign-leaf source seam (shared by the CPU and GPU walks -- the single home for foreign-leaf
 * secondary-source semantics, so the two walks cannot drift).  A foreign single particle is shipped as a
 * synthesized terminal multipole; the node payload carries mass/pos/etc. (singleton-exact), but it CANNOT
 * carry the leaf-identity fields the force needs.  When the accepted foreign node is a tagged real leaf,
 * restore them so the unchanged grav_force_pair reproduces the source's local-leaf interaction on its home
 * rank (preserving force reciprocity):
 *   - h_p <- leaf_force_soft, which the producer fills from ForceSoftening_KernelRadius() -- the SAME
 *     local-leaf softening accessor the home rank's walk uses for h_p (for adaptive types this is the
 *     adaptive kernel radius, not the fixed table value).  Restoring it makes the accepted foreign leaf
 *     use the identical h_p it would on its home rank, rather than the node-payload maxsoft.
 *   - ptype_sec/zeta_sec <- Type and AGS_zeta, consumed by the AGS symmetrize-by-averaging + zeta-correction
 *     paths (no-ops when those are compiled out / zeta==0, but required for general configs).
 * Returns 1 if leaf identity was applied, 0 for a true node/multipole. */
KOKKOS_INLINE_FUNCTION
int grav_apply_foreign_leaf_identity(int leaf_tag, int leaf_type, double leaf_zeta, double leaf_force_soft,
                                     int *ptype_sec, double *zeta_sec, double *h_p)
{
    if(leaf_tag == LET_LEAF_TAG_REAL_LEAF) { *ptype_sec = leaf_type; *zeta_sec = leaf_zeta; *h_p = leaf_force_soft; return 1; }
    return 0;
}


/* zeta/zeta_sec (AGS-softening correction) and ags_bitflag_primary (symmetrize-by-
 * averaging) are passed unconditionally; callers supply 0 when the corresponding flag
 * is off, and the body consumes them only under the matching #if (no behavior change),
 * so the signature stays flag-independent and callers cannot drift out of gate order. */
KOKKOS_INLINE_FUNCTION grav_force_pair_t grav_force_pair(
    double r, double r2, double mass, double h, double h_p,
    int ptype, int ptype_sec, double pmass,
    double zeta, double zeta_sec, int ags_bitflag_primary
    )
{
    grav_force_pair_t out;
    double fac_accel, h_inv, h3_inv, h_p_inv, h_p3_inv, u, u_p;
#ifdef EVALPOTENTIAL
    double fac_pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    double fac_tidal, fac2_tidal;
#endif
    (void) ptype; (void) ptype_sec; (void) pmass;

    /* now we compute the actual pair-wise gravity terms */
    if((r >= h) && (r >= h_p)) // can safely do a purely-Newtonian force (can be done by kernel-gravity as well, but no need to enter all the conditional statements below so just do it here for simplicity //
    {
        fac_accel = mass / (r2 * r);
#ifdef EVALPOTENTIAL
        fac_pot = -mass / r;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        fac_tidal = fac_accel; fac2_tidal = 3.0 * mass / (r2 * r2 * r); /* second derivative of potential needs this factor */
#endif
    }
    else
    {
        double h_grav = h;
#if !defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
        if(h_p > h_grav) {h_grav = h_p;} // in this case, symmetrize by taking the maximum here always
#endif
        h_inv=1./h_grav; h3_inv=h_inv*h_inv*h_inv; u=r*h_inv; // set here to ensure this is using the correct values //
        fac_accel = mass * kernel_gravity(u, h_inv, h3_inv, 1);
#ifdef EVALPOTENTIAL
        fac_pot = mass * kernel_gravity(u, h_inv, h3_inv, -1);
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* second derivatives needed -> calculate them from softened potential */
        fac_tidal = fac_accel; fac2_tidal = mass * kernel_gravity(u, h_inv, h3_inv, 2);  /* save original fac_accel without shortrange_table factor or zeta terms (needed for tidal field calculation) */
#endif

#if defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
        if(h_p > 0) // first, appropriately symmetrize the forces between particles. only do this is secondary is a particle, so has a type and softening! //
        {
            int symmetrize_by_averaging = 0; // default here to symmetrize by taking the maximum, but this will vary below //
            // the 'zeta' terms for conservation with adaptive softening assume kernel-scale forces are averaged to symmetrize, to make them continuous
            if(ptype_sec>=0) {if((1 << ptype_sec) & (ags_bitflag_primary)) {symmetrize_by_averaging=1;}} // symmetrize by averaging only for particles which have a shared AGS structure since this is how our correction terms are derived //
#ifdef SINGLE_STAR_SINK_DYNAMICS
            if((ptype!=0) || (ptype_sec!=0)) {symmetrize_by_averaging=0;} // we don't want to do the symmetrization below for sink interactions because it can create very noisy interactions between tiny sink particles and diffuse gas. However we do want it for gas-gas interactions so we keep the below
#endif
            double prefac_corr_p=1., prefac_corr_orig=1.; // this will give a symmetrized pair by linear averaging
            if(symmetrize_by_averaging==0) {prefac_corr_p=2; prefac_corr_orig=0.;} // symmetrize instead with the old method of simply taking the larger of the pair. here only act if the softening of the particle whose force is being summed is greater than the target //
            if((symmetrize_by_averaging==1) || (h_p>h)) // condition to need to evaluate the alternate particle ('p' side)
            {
                h_p_inv=1./h_p; h_p3_inv=h_p_inv*h_p_inv*h_p_inv; u_p=r*h_p_inv;
                fac_accel = 0.5 * (prefac_corr_orig * fac_accel + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1)); // average with neighbor
#ifdef EVALPOTENTIAL
                fac_pot = 0.5 * (prefac_corr_orig * fac_pot + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, -1)); // average with neighbor
#endif
#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
                fac_tidal = fac_accel; fac2_tidal = 0.5 * (prefac_corr_orig * fac2_tidal + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 2)); // average forces -> average in tidal tensor as well. also save updated fac_tidal
#endif
            }
        } // closes if((h_p > 0)) clause
#endif // closes clause to symmetrize by averaging instead of taking the larger softening

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        double fac_corr = 0; int add_ags_zeta_terms_primary=1, add_ags_zeta_terms_secondary=1; u_p=r/h_p; // define the correction factor but also a clause to see if we should apply any 'correction' term at all
        if((r<=0) || (pmass<=0) || (mass<=0) || (ptype_sec<0)) {add_ags_zeta_terms_primary=0; add_ags_zeta_terms_secondary=0;} // define conditions to add these terms at all (don't go forward if any of these conditions are met)
        if((zeta == 0) || (u >= 1) || (h <= 0)) {add_ags_zeta_terms_primary=0;} // other conditions that mean -dont- use the term for the ab side
        if((zeta_sec == 0) || (u_p >= 1) || (h_p <= 0)) {add_ags_zeta_terms_secondary=0;} // other conditions that mean -dont- use the term for the ba side
        // correction only applies to 'shared-kernel' particles: so this needs to check if these are the same particles for which the 'shared' kernel lengths are computed.
        // ptype_sec>=0 guard: node sources carry ptype_sec=-1, for which both terms are already disabled above; the guard makes that explicit and avoids the legacy walk's `1<<(-1)` shift UB (result was discarded). Byte-identical to the legacy CPU/GPU walks.
        if(ptype_sec >= 0 && (ptype != 0 || ptype_sec != 0)) {
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
            if(!((1 << ptype) & (ADAPTIVE_GRAVSOFT_FORALL)) || !((1 << ptype_sec) & (gravtree_ags_kernel_shared_bitflag(ptype)))) {add_ags_zeta_terms_primary=0;} // primary must be a valid ags particle and 'see' secondary for ab side
            if(!((1 << ptype_sec) & (ADAPTIVE_GRAVSOFT_FORALL)) || !((1 << ptype) & (gravtree_ags_kernel_shared_bitflag(ptype_sec)))) {add_ags_zeta_terms_secondary=0;} // secondary must be a valid ags particle and 'see' primary for ba side
#else
            add_ags_zeta_terms_primary=0; add_ags_zeta_terms_secondary=0; // primary and secondary must be gas for ab side or ba side
#endif
        }
        if(add_ags_zeta_terms_primary) // ab side
        {
            double dWdr, wp; kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
            fac_corr += -(zeta/pmass) * dWdr / r; // go ahead and add the term
        }
        if(add_ags_zeta_terms_secondary) // ba side
        {
            double dWdr, wp; h_p_inv=1./h_p; h_p3_inv=h_p_inv*h_p_inv*h_p_inv; kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
            fac_corr += -(zeta_sec/pmass) * dWdr / r; // go ahead and add the term
        }
        if(!isnan(fac_corr)) {fac_accel += fac_corr;}
#endif
    } // closes r < h (else) clause [where we need to deal with inside-the-softening factors]

    out.fac_accel = fac_accel;
#ifdef EVALPOTENTIAL
    out.fac_pot = fac_pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    out.fac_tidal = fac_tidal; out.fac2_tidal = fac2_tidal;
#endif
    return out;
}


/* ------------------------------------------------------------------------------------------
 * TreePM short-range truncation (tabindex + table application). The OOR semantics follow the
 * CPU reference: an out-of-range source contributes nothing to acc/pot/DF/tidal/jerk (the walk
 * skips that gated block) but its RAW fac_accel still feeds the column-integral payloads
 * (TREECOL) below the gate — there is no PM-side completion for those integrals.
 * ---------------------------------------------------------------------------------------- */
#ifdef PMGRID
KOKKOS_INLINE_FUNCTION int grav_pm_shortrange_tabindex(double asmthfac, double r)
{
    return (int) (asmthfac * r);
}
KOKKOS_INLINE_FUNCTION int grav_pm_shortrange_in_range(int tabindex)
{
    return (tabindex < GRAVTREE_SHORTRANGE_NTAB && tabindex >= 0);
}
KOKKOS_INLINE_FUNCTION void grav_force_apply_pm_truncation(grav_pm_shortrange_t pm, int tabindex,
                                                           double &fac_pot, double &fac_accel)
{
    fac_accel *= pm.shortrange_tab[tabindex];
#ifdef EVALPOTENTIAL
    fac_pot *= pm.shortrange_pot_tab[tabindex];   /* fac_pot is a dead 0 when !EVALPOTENTIAL */
#endif
}
#endif /* PMGRID */


/* ------------------------------------------------------------------------------------------
 * GRAVITY_SPHERICAL_SYMMETRY: shell-theorem force law (enclosed mass toward the box center).
 * ---------------------------------------------------------------------------------------- */
#ifdef GRAVITY_SPHERICAL_SYMMETRY
KOKKOS_INLINE_FUNCTION double grav_spherical_symmetry_r_from_center(double x0, double x1, double x2, double c0, double c1, double c2)
{
    return sqrt(pow(x0 - c0,2) + pow(x1 - c1,2) + pow(x2 - c2,2));
}
/* Overrides dr (to point at the box center) and fac_accel for the shell-theorem force; pot is
 * left as computed, matching the walk sequencing. */
template <typename VecT>
KOKKOS_INLINE_FUNCTION void grav_spherical_symmetry_force_override(double r_source, double r_target, double h, double mass,
                                                                   double c0, double c1, double c2, double p0, double p1, double p2,
                                                                   VecT &dr, double &fac_accel)
{
    if(r_source < r_target) {dr[0] = c0 - p0; dr[1] = c1 - p1; dr[2] = c2 - p2; fac_accel = mass/pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),3);} else {fac_accel = 0;}
}
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
KOKKOS_INLINE_FUNCTION double grav_spherical_symmetry_fac2_tidal_override(double r_source, double r_target, double h, double mass)
{
    if(r_source < r_target) {return 3 * mass / pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),5);} else {return 0;}
}
#endif
#endif /* GRAVITY_SPHERICAL_SYMMETRY */


/* ------------------------------------------------------------------------------------------
 * SINK_DYNFRICTION_FROMTREE: dynamical-friction deflection terms for sink targets.
 * ---------------------------------------------------------------------------------------- */
#if defined(SINK_DYNFRICTION_FROMTREE)
template <typename AccVecT>
KOKKOS_INLINE_FUNCTION void grav_sink_dynfriction_accumulate(const Vec3<double> &dr, const Vec3<double> &dv,
                                                             double fac_accel, double mass, double sink_mass, double m_j_eff_for_df,
                                                             int ptype, AccVecT &acc)
{
    if( (fac_accel>MIN_REAL_NUMBER) && (ptype==5) && (mass>MIN_REAL_NUMBER) )
    {
        double dv2=dv.norm_sq();
        if((dv2 > MIN_REAL_NUMBER) && (sink_mass > MIN_REAL_NUMBER))
        {
            double dv0=sqrt(dv2); Vec3<double> dv_h=dv/dv0;
            double rdotvhat=dot(dr,dv_h);
            Vec3<double> b_im=dr-rdotvhat*dv_h; double b_impact=b_im.norm();
            double a_im=(b_impact*All.cf_atime)*(dv2*All.cf_a2inv)/(All.G*sink_mass), fac_df=fac_accel*b_impact*a_im/(1.+a_im*a_im); // need to convert to fully-physical units to ensure this has the correct dimensions
            /* this is where we can insert an ad-hoc renormalization to avoid double-counting if we have a genuinely very massive BH (so DF is well-resolved) */
            {
                double m_j=m_j_eff_for_df; /* estimate mean mass of the particles in the node */
                if(sink_mass > 14.251*m_j) {fac_df *= DMIN(1.,DMAX(0.,(-1.+3./log10(sink_mass/m_j))/1.6));} /* approximate correction factor estimated by linhao */
            }
            if((m_j_eff_for_df <= MIN_REAL_NUMBER) || (b_impact <= MIN_REAL_NUMBER) || (dv2 <= MIN_REAL_NUMBER)) {fac_df = 0;}
            /* parallel deflection component: dv = V[distant particle/node] - V[bh], sign here is set to accelerate towards V[ext], as needed */
            acc += fac_df * dv_h;
            /* perpendicular deflection component b_im = P[distant particle/node] - P[bh], so positive = accel -towards- P[ext], but this is the residual term (after subtracting the homogeneous term), which points in the opposite direction */
            double fac_df_p = -fac_df / (b_impact * a_im + MIN_REAL_NUMBER);
            if(fabs(fac_df_p)<MAX_REAL_NUMBER && isfinite(fac_df_p)) {acc += fac_df_p * b_im;}
        }
    }
}
#endif /* SINK_DYNFRICTION_FROMTREE */


/* ------------------------------------------------------------------------------------------
 * ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION: variable-softening 'correction' terms (analogous to
 * the AGS zeta terms). The secondary-uses gate carries an explicit ptype_sec>=0 guard (node
 * sources have ptype_sec=-1; a negative shift is undefined — on all production targets the
 * legacy expression reduced to 0, which the guard makes explicit).
 * ---------------------------------------------------------------------------------------- */
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
template <typename ZTensorT, typename AccVecT>
KOKKOS_INLINE_FUNCTION void grav_ags_tidal_criterion_accumulate(double r, double r2, const Vec3<double> &dr,
                                                                double mass, double h, double h_p, int ptype, int ptype_sec,
                                                                double fac_tidal, double fac2_tidal,
                                                                const ZTensorT &i_zeta_tidal_tensorps_prevstep,
                                                                const ZTensorT &j_zeta_tidal_tensorps_prevstep,
                                                                double &tidal_zeta, AccVecT &acc)
{
    int primary_uses_tidal_criterion=0, secondary_uses_tidal_criterion=0;
    if(mass > 0 && r2 > 0)
    {
        if((1 << ptype) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) {primary_uses_tidal_criterion=1;} /* check if the primary particle uses the tidal softening */
        if(ptype_sec >= 0 && ((1 << ptype_sec) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION))) {secondary_uses_tidal_criterion=1;} /* check if the secondary particle uses the tidal softening */
        double prefac_tt=0.5, h_touse=h, u_tt=sqrt(r2)/h_touse; // this corresponds to the result of symmetrizing by averaging
#if !defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
        if(h >= h_p) {prefac_tt=1;} else {prefac_tt=1; h_touse=h_p; u_tt=sqrt(r2)/h_touse;} // this corresponds to adopting the MAX criterion for the softening
#endif
        if(u_tt<1 && prefac_tt>0) {tidal_zeta += prefac_tt * mass * kernel_gravity(u_tt,1./h,1./(h*h*h),0);} // simple sum to calculate this contribution, only from particles inside the kernel of the primary -- this is up here instead of below the if below because it needs to include the 'self' contribution here
    }
    if(primary_uses_tidal_criterion || secondary_uses_tidal_criterion) // primary or secondary has associated correction terms here
    { // now this is correct, but always need to carefully ensure correction terms are only applied in the correct 'direction' if we have a mixed-particle-type pair //
        double h_touse = DMAX(h, h_p), f_b = -r*fac2_tidal, f_a = (6./r)*fac_tidal; // these will give the correct factors for the correction terms below, automatically symmetrized appropriately based on the same symmetry rules we use to define the tidal tensor itself in the first place
        if(r < h_touse)
        {
            double dwk,wk,f_a_corr,u; u=r/h_touse; kernel_main(u,1.,1.,&wk,&dwk,0); // gather the remaining kernel terms (note these derivatives come from the laplacian and its derivative, so can be reconstructed from our usual wk and dwk terms //
            f_a_corr = 4.*M_PI*mass * (dwk - (2./u)*wk) / (h_touse*h_touse*h_touse*h_touse); // default to symmetrize by taking the maximum, here
#if defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
            if(h<h_p) {h_touse=h;} else {h_touse=h_p;}
            u=r/h_touse; kernel_main(u,1.,1.,&wk,&dwk,0);
            f_a_corr = 0.5 * (f_a_corr + 4.*M_PI*mass * (dwk - (2./u)*wk) / pow(h_touse,4)); // symmetrize by averaging since thats what we did above
#endif
            f_a += f_a_corr; // add this to the relevant function to use below
        }
        int ki,kj,kk; double acc_corr_zeta[3]={0}; Vec3<double> rh = dr / r;
        for(kk=0;kk<3;kk++)
        {
            for(ki=0;ki<3;ki++)
            {
                for(kj=0;kj<3;kj++)
                {
                    double q0=rh[ki]*rh[kj]*rh[kk], fb_rh_add=0; /* first compute di dj dk [phi_kernel] -- this is generic */
                    if(ki==kj) {fb_rh_add+=rh[kk];} /* these are the delta_ij terms */
                    if(ki==kk) {fb_rh_add+=rh[kj];}
                    if(kj==kk) {fb_rh_add+=rh[ki];}
                    double qfun = f_a * q0 + f_b * (-3.*q0 + fb_rh_add); /* now double-dot this properly to get the sum we need - note only here need the combination of TT and zeta terms */
                    acc_corr_zeta[kk] += primary_uses_tidal_criterion * i_zeta_tidal_tensorps_prevstep[ki][kj] * qfun; // only non-zero here if primary is a tidal-softening-active particle
                    acc_corr_zeta[kk] += secondary_uses_tidal_criterion * j_zeta_tidal_tensorps_prevstep[ki][kj] * qfun; // only non-zero here if secondary is a tidal-softening-active particle
                }
            }
        }
        acc[0]+=acc_corr_zeta[0]; acc[1]+=acc_corr_zeta[1]; acc[2]+=acc_corr_zeta[2]; // final assignment
    }
}
#endif /* ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION */


/* ------------------------------------------------------------------------------------------
 * COMPUTE_TIDAL_TENSOR_IN_GRAVTREE: tidal-tensor accumulation (PM-truncated or bare).
 * st/sti are shortrange_table[tabindex] / shortrange_table_tidal[tabindex] (in range here —
 * this call sits inside the PM short-range gate).
 * ---------------------------------------------------------------------------------------- */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
template <typename TensorT>
KOKKOS_INLINE_FUNCTION void grav_tidal_tensor_accumulate(const Vec3<double> &dr, double fac_tidal, double fac2_tidal,
                                                         grav_pm_shortrange_t pm, int tabindex,
                                                         TensorT &tidal_tensorps)
{
    /* tidal_tensorps[][] = Matrix of second derivatives of grav. potential, symmetric:
     |Txx Txy Txz|   |tidal_tensorps[0][0] tidal_tensorps[0][1] tidal_tensorps[0][2]|
     |Tyx Tyy Tyz| = |tidal_tensorps[1][0] tidal_tensorps[1][1] tidal_tensorps[1][2]|
     |Tzx Tzy Tzz|   |tidal_tensorps[2][0] tidal_tensorps[2][1] tidal_tensorps[2][2]|  */
#ifdef PMGRID
    /* short-range PM truncation factors at this pair's table index (caller guarantees
     * tabindex is in range under PMGRID -- this helper is called inside the PM gate). */
    float st = pm.shortrange_tab[tabindex], sti = pm.shortrange_tidal_tab[tabindex];
    tidal_tensorps[0][0] += ((-fac_tidal + dr[0] * dr[0] * fac2_tidal) * st) +
    dr[0] * dr[0] * fac2_tidal / 3.0 * sti;
    tidal_tensorps[0][1] += ((dr[0] * dr[1] * fac2_tidal) * st) +
    dr[0] * dr[1] * fac2_tidal / 3.0 * sti;
    tidal_tensorps[0][2] += ((dr[0] * dr[2] * fac2_tidal) * st) +
    dr[0] * dr[2] * fac2_tidal / 3.0 * sti;
    tidal_tensorps[1][1] += ((-fac_tidal + dr[1] * dr[1] * fac2_tidal) * st) +
    dr[1] * dr[1] * fac2_tidal / 3.0 * sti;
    tidal_tensorps[1][2] += ((dr[1] * dr[2] * fac2_tidal) * st) +
    dr[1] * dr[2] * fac2_tidal / 3.0 * sti;
    tidal_tensorps[2][2] += ((-fac_tidal + dr[2] * dr[2] * fac2_tidal) * st) +
    dr[2] * dr[2] * fac2_tidal / 3.0 * sti;
#else
    tidal_tensorps[0][0] += (-fac_tidal + dr[0] * dr[0] * fac2_tidal);
    tidal_tensorps[0][1] += (dr[0] * dr[1] * fac2_tidal);
    tidal_tensorps[0][2] += (dr[0] * dr[2] * fac2_tidal);
    tidal_tensorps[1][1] += (-fac_tidal + dr[1] * dr[1] * fac2_tidal);
    tidal_tensorps[1][2] += (dr[1] * dr[2] * fac2_tidal);
    tidal_tensorps[2][2] += (-fac_tidal + dr[2] * dr[2] * fac2_tidal);
#endif
}
#endif /* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE */


/* ------------------------------------------------------------------------------------------
 * COMPUTE_JERK_IN_GRAVTREE: jerk accumulation (inside the PM short-range gate).
 * ---------------------------------------------------------------------------------------- */
#ifdef COMPUTE_JERK_IN_GRAVTREE
template <typename VecT>
KOKKOS_INLINE_FUNCTION void grav_jerk_accumulate(const Vec3<double> &dv, const Vec3<double> &dr,
                                                 double fac_accel, double fac2_tidal, int ptype, VecT &jerk)
{
    (void) ptype;
#ifndef ADAPTIVE_TREEFORCE_UPDATE // we want the jerk if we're doing lazy force updates
    if(ptype > 0)
#endif
    {
        double dv_dot_dr = dot(dv, dr);
        jerk += fac_accel * dv - dv_dot_dr * fac2_tidal * dr;
    }
}
#endif /* COMPUTE_JERK_IN_GRAVTREE */


/* ------------------------------------------------------------------------------------------
 * Accepted-interaction record, and the shared evaluation of one accepted element.
 *
 * A walk has two phases per target: traversal decides which elements the target interacts
 * with, and evaluation accumulates each of them. What survives from one to the other is only
 * the element's identity and the classification the traversal already computed; every value
 * the evaluation needs is re-derived from that identity through the walker's own storage, in
 * the same statements the walker used before the split, so evaluating an element later reads
 * exactly the state it would have read at encounter. That holds because every source drift
 * happens during traversal, before the record is appended, and drifting is idempotent within
 * a step; Hermite prediction is re-derived at evaluation from that drifted state.
 * ---------------------------------------------------------------------------------------- */
struct grav_walk_record_t {
    int no;                  /* particle index (no < TreeParticleSlots) or node index */
    grav_node_kind_t kind;   /* node only: the walk's classification of this node */
    int leaf_tag;            /* node only: the LET leaf tag read for this node (LET_LEAF_TAG_*) */
};

/* Per-target inputs to the pair evaluation, fixed for the whole walk of one target. */
struct grav_pair_tgt_t {
    int ptype;
    double pmass;
    double h;                     /* target softening kernel radius */
    double zeta;                  /* adaptive-softening correction term (0 when unused) */
    int ags_bitflag;              /* types this target shares an adaptive-softening kernel with */
    grav_pm_shortrange_t pm;      /* PM short-range snapshot incl. any high-res override (empty when !PMGRID) */
#ifdef SINK_DYNFRICTION_FROMTREE
    double sink_mass;             /* sink mass of a sink target, else 0 */
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<double> i_zeta_tidal_tensorps_prevstep;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    Vec3<double> pos;
    double center[3];
#endif
};

/* Per-element inputs, loaded by the walker from its own storage for each accepted element.
 * dr is the wrapped separation source - target and is overridden in place by the spherical-
 * symmetry force law, so the payload blocks that follow see the same dr the accumulation used. */
struct grav_pair_src_t {
    Vec3<double> dr;
    double r2, mass, h_p, zeta_sec;
    int ptype_sec;                /* -1 for a node source */
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<double> dv;              /* source velocity - target velocity */
#endif
#ifdef SINK_DYNFRICTION_FROMTREE
    double m_j_eff_for_df;        /* mean particle mass of the source (its mass for a leaf) */
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double r_source;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<double> j_zeta_tidal_tensorps_prevstep;
#endif
};

/* The accumulators the pair evaluation writes. */
struct grav_pair_acc_t {
    Vec3<double> acc;
    double pot;                   /* accumulated only under EVALPOTENTIAL */
    int ninter;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    SymmetricTensor2<double> tidal_tensorps;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    Vec3<double> jerk;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    double tidal_zeta;
#endif
};

KOKKOS_INLINE_FUNCTION void grav_pair_acc_init(grav_pair_acc_t &out)
{
    out.acc = Vec3<double>{0,0,0}; out.pot = 0.0; out.ninter = 0;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    out.tidal_tensorps = {};
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    out.jerk = Vec3<double>{0,0,0};
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    out.tidal_zeta = 0.0;
#endif
}

/* Take another partial accumulation of the same target into `a`: every field here is a sum. */
KOKKOS_INLINE_FUNCTION void grav_pair_acc_combine(grav_pair_acc_t &a, const grav_pair_acc_t &b)
{
    a.acc += b.acc; a.pot += b.pot; a.ninter += b.ninter;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    for(int k = 0; k < 6; k++) {a.tidal_tensorps.data[k] += b.tidal_tensorps.data[k];}
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    a.jerk += b.jerk;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    a.tidal_zeta += b.tidal_zeta;
#endif
}

/* What the walker's payload blocks consume after the evaluation: the separation and the
 * acceleration prefactor as accumulated (PM-truncated when the source was in range of the
 * short-range table, raw otherwise -- the column-integral payloads have no PM completion). */
struct grav_pair_result_t {
    double r;
    double fac_accel;
};

/* Evaluate one accepted element for one target: the pair force, the PM short-range gate and
 * truncation, the spherical-symmetry override, and every accumulation that lives inside the
 * gate, in the order the CPU walk has always used. The periodic-image potential correction is
 * NOT here: each walker adds it beside its own dr, since the two walks read dr at different
 * points relative to the spherical-symmetry override. Under SINK_DYNFRICTION_FROMTREE together
 * with ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION the device walk formerly added dynamical friction
 * after the tidal correction; both are pure additions into acc, and the CPU order is used for
 * both walks now. */
KOKKOS_INLINE_FUNCTION grav_pair_result_t
grav_pair_evaluate_core(const grav_pair_tgt_t &tgt, grav_pair_src_t &src, grav_pair_acc_t &out)
{
    const double r = sqrt(src.r2);
    double fac_accel;
    double fac_pot = 0;   /* a dead 0 when !EVALPOTENTIAL (consumed only under that gate) */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    double fac_tidal = 0.0, fac2_tidal = 0.0;
#endif
    {
        grav_force_pair_t pair_out = grav_force_pair(r, src.r2, src.mass, tgt.h, src.h_p, tgt.ptype, src.ptype_sec, tgt.pmass,
                                                     tgt.zeta, src.zeta_sec, tgt.ags_bitflag);
        fac_accel = pair_out.fac_accel;
#ifdef EVALPOTENTIAL
        fac_pot = pair_out.fac_pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        fac_tidal = pair_out.fac_tidal; fac2_tidal = pair_out.fac2_tidal;
#endif
    }
    int tabindex = 0;   /* computed and consumed only under PMGRID */
#ifdef PMGRID
    tabindex = grav_pm_shortrange_tabindex(tgt.pm.asmthfac, r);
    if(grav_pm_shortrange_in_range(tabindex))
#endif
    {
#ifdef PMGRID
        grav_force_apply_pm_truncation(tgt.pm, tabindex, fac_pot, fac_accel);
#endif
#ifdef EVALPOTENTIAL
        out.pot += fac_pot;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
        double r_target = grav_spherical_symmetry_r_from_center(tgt.pos[0],tgt.pos[1],tgt.pos[2],tgt.center[0],tgt.center[1],tgt.center[2]);
        grav_spherical_symmetry_force_override(src.r_source, r_target, tgt.h, src.mass, tgt.center[0],tgt.center[1],tgt.center[2],
                                               tgt.pos[0],tgt.pos[1],tgt.pos[2], src.dr, fac_accel);
#endif
        out.acc += fac_accel * src.dr;
#ifdef SINK_DYNFRICTION_FROMTREE
        grav_sink_dynfriction_accumulate(src.dr, src.dv, fac_accel, src.mass, tgt.sink_mass, src.m_j_eff_for_df, tgt.ptype, out.acc);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        grav_ags_tidal_criterion_accumulate(r, src.r2, src.dr, src.mass, tgt.h, src.h_p, tgt.ptype, src.ptype_sec, fac_tidal, fac2_tidal,
                                            tgt.i_zeta_tidal_tensorps_prevstep, src.j_zeta_tidal_tensorps_prevstep, out.tidal_zeta, out.acc);
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
#ifdef GRAVITY_SPHERICAL_SYMMETRY
        fac2_tidal = grav_spherical_symmetry_fac2_tidal_override(src.r_source, r_target, tgt.h, src.mass);
#endif
        grav_tidal_tensor_accumulate(src.dr, fac_tidal, fac2_tidal, tgt.pm, tabindex, out.tidal_tensorps);
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        grav_jerk_accumulate(src.dv, src.dr, fac_accel, fac2_tidal, tgt.ptype, out.jerk);
#endif
    }
    out.ninter++;
    grav_pair_result_t res; res.r = r; res.fac_accel = fac_accel;
    return res;
}


/* ------------------------------------------------------------------------------------------
 * RT_USE_TREECOL_FOR_NH: six-bin angular column-density estimate. Sits OUTSIDE the PM gate by
 * design — fac_accel here is the raw un-truncated value (no PM-side completion exists for the
 * column integral).
 * ---------------------------------------------------------------------------------------- */
#ifdef RT_USE_TREECOL_FOR_NH
KOKKOS_INLINE_FUNCTION void grav_treecol_accumulate(const Vec3<double> &dr, double r, double fac_accel,
                                                    double gasmass, double mass, double angular_bin_size,
                                                    double *treecol_angular_bins)
{
    if(gasmass>0)
    {
        int bin; // Here we do a simple six-bin angular binning scheme
        if((fabs(dr[0]) > fabs(dr[1])) && (fabs(dr[0])>fabs(dr[2]))) {if (dr[0] > 0) {bin = 0;} else {bin=1;}
        } else if (fabs(dr[1])>fabs(dr[2])){if (dr[1] > 0) {bin = 2;} else {bin=3;}
        } else {if (dr[2] > 0) {bin = 4;} else {bin = 5;}}
        treecol_angular_bins[bin] += fac_accel*gasmass*r / (angular_bin_size*mass); // in our binning scheme, we stretch the gas mass over a patch  of the sphere located at radius r subtending solid angle equal to the bin size - thus the area is r^2 * angular_bin_size, so sigma = m/(r^2 * angular bin size) = fac_accel/r / angular bin size. Factor of gasmass / mass corrects the gravitational mass to the gas mass
    }
}
#endif /* RT_USE_TREECOL_FOR_NH */


/* ------------------------------------------------------------------------------------------
 * COSMIC_RAY_SUBGRID_LEBRON: sub-grid CR energy-density accumulation. t_max_cr and the active
 * gate are evaluated caller-side once per target (the time helper is host-only); each walk
 * supplies its legacy gate expression.
 * ---------------------------------------------------------------------------------------- */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
KOKKOS_INLINE_FUNCTION void grav_cr_lebron_accumulate(int ptype, double r, double soft, double cr_injection,
                                                      int cr_active_gate, double t_max_cr,
                                                      grav_pm_shortrange_t pm,
                                                      double &SubGrid_CosmicRayEnergyDensity)
{
    if(ptype==0 && r>0 && cr_injection>0 && cr_active_gate)
    {
        double kappa_0 = All.CosmicRay_Subgrid_Kappa_0, vst_0 = All.CosmicRay_Subgrid_Vstream_0; // in code units
        double r_phys = sqrt(r*r + soft*soft/4.) * All.cf_atime, t_max = t_max_cr; // make sure we're working in physical code units, and assign max time to formation at begin time, and include very crude 'softening' term here to prevent divergennce as r->0: for our default parameters can't be too large here or we get unphysically large CR halos compared to reality, b/c of large effective streaming terms
        double r_max = 0.5*t_max*vst_0 * (1. + sqrt(1. + 16.*kappa_0/(vst_0*vst_0*t_max))); // maximum stream distance
#ifdef PMGRID
        r_max = DMIN(r_max , 0.5*pm.rcut*All.cf_atime); // truncate before reach the boundary of the grid to avoid numerical errors there
#endif
        double fac_cr_distance = 1./(4.*M_PI*r_phys*(kappa_0 + vst_0*r_phys)) * exp(-DMIN(r_phys*r_phys/(1.e-6*r_phys*r_phys+r_max*r_max),50.));
        if(fac_cr_distance>0) {SubGrid_CosmicRayEnergyDensity += fac_cr_distance * cr_injection / All.cf_a3inv;} // convert to appropriate code units for an energy density or pressure
    }
}
#endif /* COSMIC_RAY_SUBGRID_LEBRON */


/* ------------------------------------------------------------------------------------------
 * RT_USE_GRAVTREE payload family: optically-thin incident flux / energy / Eddington-tensor /
 * radiation-pressure accumulation from the loaded source luminosities. Caller gates on the
 * valid-gas target check; outputs are the per-family accumulators (each compile-gated exactly
 * as the walk locals).
 * ---------------------------------------------------------------------------------------- */
#ifdef RT_USE_GRAVTREE
/* Per-pair source inputs (read-only, passed by value). Optional fields gated once here. */
struct grav_rt_src_t {
    Vec3<double> d_stellarlum;
    double soft;
    const double *mass_stellarlum;
#ifdef CHIMES_STELLAR_FLUXES
    const double *chimes_mass_stellarlum_G0, *chimes_mass_stellarlum_ion;
#endif
#ifdef SINK_PHOTONMOMENTUM
    double mass_sinklumwt_forradfb;
#endif
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    const double *fac_stellum; // unused when SAVE_RAD_FLUX path stores fluxes instead
#endif
};

/* Per-target output accumulators: a live VIEW (pointers) over the walk's existing
 * per-target locals. The walk still owns zeroing + writeback; this only adds into them. */
struct grav_rt_accum_t {
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    double *Rad_E_gamma;
#endif
#ifdef CHIMES_STELLAR_FLUXES
    double *chimes_flux_G0, *chimes_flux_ion;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    double *incident_flux_uv, *incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
    double *incident_flux_agn;
#endif
#ifdef RT_OTVET
    SymmetricTensor2<double> *RT_ET;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    Vec3<double> *Rad_Flux;
#endif
};

template <typename AccVecT>
KOKKOS_INLINE_FUNCTION void grav_rt_payload_accumulate(grav_rt_src_t src, grav_rt_accum_t rt_accum, AccVecT &acc)
{
    (void) acc;
    double r2 = src.d_stellarlum.norm_sq(); double r = sqrt(r2); double fac_rt;
    if(r >= src.soft) {fac_rt=1./(r2*r);} else {double h_inv_rt=1./src.soft, h3_inv_rt=h_inv_rt*h_inv_rt*h_inv_rt, u_rt=r*h_inv_rt; fac_rt=kernel_gravity(u_rt,h_inv_rt,h3_inv_rt,1);}
    if((src.soft>r)&&(src.soft>0)) fac_rt *= (r2/(src.soft*src.soft)); // don't allow cross-section > r2
    double fac_intensity; fac_intensity = fac_rt * r * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^2), in -physical- units, since L is physical
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {rt_accum.Rad_E_gamma[kf] += fac_intensity * src.mass_stellarlum[kf];}}
#endif
#ifdef CHIMES_STELLAR_FLUXES
    int chimes_k; double chimes_fac = fac_intensity / (UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS);  // 1/(4 * pi * r^2), in cm^-2
    for (chimes_k = 0; chimes_k < CHIMES_LOCAL_UV_NBINS; chimes_k++)
    {
        rt_accum.chimes_flux_G0[chimes_k] += chimes_fac * src.chimes_mass_stellarlum_G0[chimes_k];   // Habing flux units
        rt_accum.chimes_flux_ion[chimes_k] += chimes_fac * src.chimes_mass_stellarlum_ion[chimes_k]; // cm^-2 s^-1
    }
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    (*rt_accum.incident_flux_uv) += fac_intensity * src.mass_stellarlum[RT_FREQ_BIN_FIRE_UV]; // don't multiply by shortrange_table since that is to prevent 2x-counting by PMgrid (which never happens here) //
    if((src.mass_stellarlum[RT_FREQ_BIN_FIRE_IR]<src.mass_stellarlum[RT_FREQ_BIN_FIRE_UV])&&(src.mass_stellarlum[RT_FREQ_BIN_FIRE_IR]>0)) // if this -isn't- satisfied, no chance you are optically thin to EUV //
    {
        // here, use ratio and linear scaling of escape with tau to correct to the escape fraction for the correspondingly higher EUV kappa: factor ~2000 here comes from the ratio of (kappa_euv/kappa_uv)
        (*rt_accum.incident_flux_euv) += fac_intensity * src.mass_stellarlum[RT_FREQ_BIN_FIRE_UV] * (All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV) *
                                                                                     ((src.mass_stellarlum[RT_FREQ_BIN_FIRE_UV] + src.mass_stellarlum[RT_FREQ_BIN_FIRE_IR]) /
                                                                                      (src.mass_stellarlum[RT_FREQ_BIN_FIRE_UV] + 2042.6*src.mass_stellarlum[RT_FREQ_BIN_FIRE_IR])));
    } else {
        // here, just enforce a minimum escape fraction //
        double m_lum_total = 0; int ks_q; for(ks_q=0;ks_q<N_RT_FREQ_BINS;ks_q++) {m_lum_total += src.mass_stellarlum[ks_q];}
        (*rt_accum.incident_flux_euv) += All.PhotonMomentum_fUV * fac_intensity * m_lum_total;
    }
#endif
#ifdef SINK_PHOTONMOMENTUM
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    rt_accum.Rad_E_gamma[RT_FREQ_BIN_FIRE_IR] += fac_intensity * src.mass_sinklumwt_forradfb;
#endif
#ifdef SINK_COMPTON_HEATING
    (*rt_accum.incident_flux_agn) += fac_intensity * src.mass_sinklumwt_forradfb; // L/(4pi*r*r) analog
#endif
#endif

#ifdef RT_OTVET
    /* use the information we have here from the gravity tree (optically thin incident fluxes) to estimate the Eddington tensor */
    if(r>0)
    {
        double fac_otvet_sum=0; int kf_rt;
        for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++)
        {
            fac_otvet_sum = src.mass_stellarlum[kf_rt];
            fac_otvet_sum *= fac_rt / (1.e-37 + r); // units are not important, since ET will be dimensionless, but final ET should scale as ~luminosity/r^2
            rt_accum.RT_ET[kf_rt] += fac_otvet_sum * outer_product(src.d_stellarlum);
        }
    }
#endif

#ifdef RT_LEBRON /* now we couple radiation pressure [single-scattering] terms within this module */
#ifdef GALSF_FB_FIRE_RT_LONGRANGE /* we only allow the momentum to couple over some distance to prevent bad approximations when the distance between points here is enormous */
    if(r*UNIT_LENGTH_IN_KPC*All.cf_atime > 50.) {fac_rt=0;}
#endif
    int kf_rt; double lum_force_fac=0;
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) /* save the fluxes for use below, where we will calculate their RP normally */
    double fac_flux = -fac_rt * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^3), in -physical- units (except for last r, cancelled by dx_stellum), since L is physical
    for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {rt_accum.Rad_Flux[kf_rt] += src.mass_stellarlum[kf_rt]*fac_flux*src.d_stellarlum;}
#else /* simply apply an on-the-spot approximation and do the absorption and RP force now */
    for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {lum_force_fac += src.mass_stellarlum[kf_rt] * src.fac_stellum[kf_rt];} // add directly to forces. appropriate normalization (and sign) in 'fac_stellum'
#endif
#ifdef SINK_PHOTONMOMENTUM /* divide out PhotoMom_coupled_frac here b/c we have our own SINK_Rad_Mom factor, and don't want to double-count */
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    rt_accum.Rad_Flux[RT_FREQ_BIN_FIRE_IR] += src.mass_sinklumwt_forradfb*fac_flux*src.d_stellarlum;
#elif !defined(RT_DISABLE_RAD_PRESSURE)
    lum_force_fac += (All.Sink_Rad_MomentumFactor / (MIN_REAL_NUMBER + All.PhotonMomentum_Coupled_Fraction)) * src.mass_sinklumwt_forradfb * src.fac_stellum[N_RT_FREQ_BINS-1];
#endif
#endif
    if(lum_force_fac>0) {acc += (fac_rt*lum_force_fac) * src.d_stellarlum;}
#endif /* RT_LEBRON */
}
#endif /* RT_USE_GRAVTREE */


/* ------------------------------------------------------------------------------------------
 * DM_SCALARFIELD_SCREENING: Yukawa-screened scalar-field force on non-gas targets. Evaluated
 * OUTSIDE the main PM short-range gate, with its own table gate keyed on the dm-center
 * distance. Caller gates on (ptype != 0); d_dm is wrapped here (mirrors the walk text).
 * ---------------------------------------------------------------------------------------- */
#ifdef DM_SCALARFIELD_SCREENING
template <typename AccVecT>
KOKKOS_INLINE_FUNCTION void grav_dm_scalarfield_accumulate(Vec3<double> &d_dm, double mass_dm, double h,
                                                           grav_pm_shortrange_t pm,
                                                           AccVecT &acc)
{
    gravity_box_nearest_image(d_dm[0],d_dm[1],d_dm[2],-1);
    double r2 = d_dm.norm_sq();
    double r = sqrt(r2); double fac_dmsf, h_inv_dmsf, h3inv_dmsf, u_dmsf;
    if(r >= h) {fac_dmsf = mass_dm / (r2 * r);} else {
        h_inv_dmsf = 1.0 / h; h3inv_dmsf = h_inv_dmsf * h_inv_dmsf * h_inv_dmsf; u_dmsf = r * h_inv_dmsf;
        fac_dmsf = mass_dm * kernel_gravity(u_dmsf, h_inv_dmsf, h3inv_dmsf, 1);
    }
    /* assemble force with strength, screening length, and target charge.  */
    fac_dmsf *= All.ScalarBeta * (1 + r / All.ScalarScreeningLength) * exp(-r / All.ScalarScreeningLength);
#ifdef PMGRID
    int tabindex = (int) (pm.asmthfac * r);
    if(tabindex < GRAVTREE_SHORTRANGE_NTAB && tabindex >= 0)
#endif
    {
#ifdef PMGRID
        fac_dmsf *= pm.shortrange_tab[tabindex];
#endif
        acc += fac_dmsf * d_dm;
    }
}
#endif /* DM_SCALARFIELD_SCREENING */


/* ------------------------------------------------------------------------------------------
 * SINK_CALC_DISTANCES (+ SINGLE_STAR_* timestepping/binary payloads): nearest-sink tracking +
 * approach/freefall/feedback/orbital timescales. One accumulator POD shared by both walks;
 * leaf and node variants co-located (the node aggregate lacks the per-sink softened-energy
 * branch by construction — only N_SINK==1 nodes feed the binary search).
 * ---------------------------------------------------------------------------------------- */
#ifdef SINK_CALC_DISTANCES
struct grav_sink_prox_accum_t {
    double Min_Distance_to_Sink2;
    Vec3<double> Min_xyz_to_Sink;
#ifdef SPECIAL_POINT_MOTION
    Vec3<double> vel_of_nearest_special, acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    double weight_sum_for_special_point_smoothing;
#endif
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    double Min_Sink_OrbitalTime, comp_Mass; Vec3<double> comp_dx, comp_dv;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    double Min_Sink_Approach_Time, Min_Sink_Freefall_time;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    double Min_Sink_FeedbackTime;
#endif
};

KOKKOS_INLINE_FUNCTION void grav_sink_prox_accum_init(grav_sink_prox_accum_t &prox)
{
    prox.Min_Distance_to_Sink2 = MAX_REAL_NUMBER;
    prox.Min_xyz_to_Sink = {MAX_REAL_NUMBER,MAX_REAL_NUMBER,MAX_REAL_NUMBER};
#ifdef SPECIAL_POINT_MOTION
    prox.vel_of_nearest_special = {}; prox.acc_of_nearest_special = {};
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    prox.weight_sum_for_special_point_smoothing = 0;
#endif
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    prox.Min_Sink_OrbitalTime = MAX_REAL_NUMBER; prox.comp_Mass = 0; prox.comp_dx = {}; prox.comp_dv = {};
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    prox.Min_Sink_Approach_Time = MAX_REAL_NUMBER; prox.Min_Sink_Freefall_time = MAX_REAL_NUMBER;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    prox.Min_Sink_FeedbackTime = MAX_REAL_NUMBER;
#endif
}

/* Take another partial accumulation of the same target into `a`, by the rules the accumulate functions
 * below apply: minima, and the fields describing what attained a minimum move with it.  Under
 * SPECIAL_POINT_WEIGHTED_MOTION the special-point motion of a special-point target is a weighted sum,
 * while for any other target it is the nearest sink's motion, so the target's type decides the rule.
 * On equal minima `a` keeps its own. */
KOKKOS_INLINE_FUNCTION void grav_sink_prox_accum_combine(grav_sink_prox_accum_t &a, const grav_sink_prox_accum_t &b, int target_ptype)
{
    (void) target_ptype;
    if(b.Min_Distance_to_Sink2 < a.Min_Distance_to_Sink2) {
        a.Min_Distance_to_Sink2 = b.Min_Distance_to_Sink2;
        a.Min_xyz_to_Sink = b.Min_xyz_to_Sink;
#ifdef SPECIAL_POINT_MOTION
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
        if(target_ptype != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
        {a.vel_of_nearest_special = b.vel_of_nearest_special; a.acc_of_nearest_special = b.acc_of_nearest_special;}
#endif
    }
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    if(target_ptype == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) {
        a.vel_of_nearest_special += b.vel_of_nearest_special; a.acc_of_nearest_special += b.acc_of_nearest_special;
    }
    a.weight_sum_for_special_point_smoothing += b.weight_sum_for_special_point_smoothing;
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    if(b.Min_Sink_OrbitalTime < a.Min_Sink_OrbitalTime) {
        a.Min_Sink_OrbitalTime = b.Min_Sink_OrbitalTime; a.comp_Mass = b.comp_Mass; a.comp_dx = b.comp_dx; a.comp_dv = b.comp_dv;
    }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    if(b.Min_Sink_Approach_Time < a.Min_Sink_Approach_Time) {a.Min_Sink_Approach_Time = b.Min_Sink_Approach_Time;}
    if(b.Min_Sink_Freefall_time < a.Min_Sink_Freefall_time) {a.Min_Sink_Freefall_time = b.Min_Sink_Freefall_time;}
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    if(b.Min_Sink_FeedbackTime < a.Min_Sink_FeedbackTime) {a.Min_Sink_FeedbackTime = b.Min_Sink_FeedbackTime;}
#endif
}

/* Target-side state shared by the leaf + node sink-proximity helpers (by value). */
struct grav_sink_prox_target_t {
    int ptype; double pmass, soft;
#ifdef SINGLE_STAR_TIMESTEPPING
    Vec3<double> vel;
#endif
};

/* Source kinematics shared by the leaf + node sink-proximity helpers. Field set is the union
 * of what either side reads (vel always; the others gated); the unused side leaves it zero. */
struct grav_sink_source_motion_t {
    Vec3<double> vel;
#if defined(SPECIAL_POINT_MOTION) || defined(SPECIAL_POINT_WEIGHTED_MOTION)
    Vec3<double> acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    double max_feedback_vel;
#endif
};

/* Leaf (single-particle) source. */
struct grav_sink_prox_leaf_src_t {
    int src_type; double src_mass;
    grav_sink_source_motion_t motion;
};

/* Node (aggregate) source. */
struct grav_sink_prox_node_src_t {
    double sink_mass;
#if defined(SINGLE_STAR_FIND_BINARIES)
    int n_sink;
#endif
    grav_sink_source_motion_t motion;
};

/* Leaf-source variant. Caller gates on (r2 > 0 && mass > 0); source fields loaded walk-side. */
KOKKOS_INLINE_FUNCTION void grav_sink_prox_leaf_accumulate(double r2, const Vec3<double> &dr,
                                                           grav_sink_prox_target_t target,
                                                           grav_sink_prox_leaf_src_t src,
                                                           grav_sink_prox_accum_t &prox)
{
    (void) target;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    if(target.ptype == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
    {
        double wt_special = grav_weight_function_for_weighted_motion_smoothing(sqrt(r2),1);
        prox.weight_sum_for_special_point_smoothing += wt_special;
        prox.vel_of_nearest_special += wt_special * src.motion.vel;
        prox.acc_of_nearest_special += wt_special * src.motion.acc;
    }
#endif
    if(src.src_type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) /* found a BH particle in grav calc */
    {
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
        if(target.ptype != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
            if(r2 < prox.Min_Distance_to_Sink2)    /* is this the closest BH part I've found yet? */
            {
                prox.Min_Distance_to_Sink2 = r2;   /* if yes: adjust min bh dist */
                prox.Min_xyz_to_Sink = dr; /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
                prox.vel_of_nearest_special = src.motion.vel;
                prox.acc_of_nearest_special = src.motion.acc;
#endif
            }
#ifdef SINGLE_STAR_TIMESTEPPING
        Vec3<double> sink_dv=src.motion.vel-target.vel; double vSqr=sink_dv.norm_sq(), M_total=src.src_mass+target.pmass, r2soft=SinkParticle_GravityKernelRadius;
        r2soft = DMAX(r2soft, target.soft);
        r2soft *= KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;
        r2soft = r2 + r2soft*r2soft;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        if(target.ptype == 0) {
            double tSqr_fb = r2soft /(src.motion.max_feedback_vel * src.motion.max_feedback_vel + MIN_REAL_NUMBER);
            if(tSqr_fb < prox.Min_Sink_FeedbackTime) {prox.Min_Sink_FeedbackTime = tSqr_fb;}
        } // for gas, add the signal velocity of feedback from the star
#endif
        double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);

        if(tSqr < prox.Min_Sink_Approach_Time) {prox.Min_Sink_Approach_Time = tSqr;}
        if(tff4 < prox.Min_Sink_Freefall_time) {prox.Min_Sink_Freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
        if(target.ptype == 5) // only for BH particles and for non center of mass calculation
        {
            double r_p5=sqrt(r2), specific_energy = 0.5*vSqr - All.G*M_total/r_p5;
            if(r2 < SinkParticle_GravityKernelRadius*SinkParticle_GravityKernelRadius)
            {
                double hinv_p5 = 1. / SinkParticle_GravityKernelRadius;
                specific_energy = 0.5*vSqr + All.G*M_total*kernel_gravity(r_p5*hinv_p5, hinv_p5, hinv_p5*hinv_p5*hinv_p5, -1);
            }
            if (specific_energy < 0)
            {
                double semimajor_axis= -All.G*M_total/(2.*specific_energy);
                double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
                if(t_orbital < prox.Min_Sink_OrbitalTime) /* Save parameters of companion */
                {
                    prox.Min_Sink_OrbitalTime=t_orbital; prox.comp_Mass=src.src_mass;
                    prox.comp_dx=dr; prox.comp_dv=sink_dv;
                }
            } /* specific_energy < 0 */
        } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
    }
}

#ifdef SPECIAL_POINT_WEIGHTED_MOTION
/* Node-source weighted-motion accumulation for special-point primaries (runs for every
 * accepted node, before the sink-aggregate block).
 *
 * SPECIAL_POINT_WEIGHTED_MOTION IS UNFINISHED AND HAS NEVER BEEN USED: a stub kept as the place to
 * implement it properly, not a form to design other code around. The intent is a weighted average:
 * every accepted source adds wt * velocity (and a matching wt * acceleration) to the sums that are
 * divided by weight_sum_for_special_point_smoothing afterwards. The node branch below does not do that
 * yet -- it overwrites the velocity with the node's and zeroes the acceleration, so the result depends
 * on visit order. */
template <typename NodeVsT>
KOKKOS_INLINE_FUNCTION void grav_sink_prox_node_specialweighted(double r2, const NodeVsT &node_vs, int ptype,
                                                                grav_sink_prox_accum_t &prox)
{
    if(ptype == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
    {
        double wt_special = grav_weight_function_for_weighted_motion_smoothing(sqrt(r2),1);
        prox.weight_sum_for_special_point_smoothing += wt_special;
        prox.vel_of_nearest_special = node_vs;
        prox.acc_of_nearest_special = {}; /* accel not computed here; left zeroed */
    }
}
#endif

/* Node-source variant. Caller gates on (node sink_mass > 0) and supplies the wrapped
 * displacement to the node's sink center-of-mass. (motion/n_sink fields in node_src exist
 * exactly when the corresponding NODE payload fields exist.) */
KOKKOS_INLINE_FUNCTION void grav_sink_prox_node_accumulate(double r2, const Vec3<double> &sink_dr,
                                                           grav_sink_prox_node_src_t src,
                                                           grav_sink_prox_target_t target,
                                                           grav_sink_prox_accum_t &prox)
{
    (void) r2; (void) src; (void) target;
    double sink_r2 = sink_dr.norm_sq(); // + (nop->len)*(nop->len);
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    if(target.ptype != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
        if(sink_r2 < prox.Min_Distance_to_Sink2)
        {
            prox.Min_Distance_to_Sink2 = sink_r2; prox.Min_xyz_to_Sink = sink_dr; /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
            prox.vel_of_nearest_special = src.motion.vel;
            prox.acc_of_nearest_special = src.motion.acc;
#endif
        }
#ifdef SINGLE_STAR_TIMESTEPPING
    Vec3<double> sink_dv = src.motion.vel - target.vel; double vSqr=sink_dv.norm_sq(), M_total=src.sink_mass+target.pmass, r2soft;
    r2soft = DMAX(SinkParticle_GravityKernelRadius, target.soft) * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER; r2soft = r2 + r2soft*r2soft;
    double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    if(target.ptype == 0) {
        double tSqr_fb = r2soft /(src.motion.max_feedback_vel * src.motion.max_feedback_vel + MIN_REAL_NUMBER);
        if(tSqr_fb < prox.Min_Sink_FeedbackTime) {prox.Min_Sink_FeedbackTime = tSqr_fb;}
    } // for gas, add the signal velocity of feedback from the star
#endif
    if(tSqr < prox.Min_Sink_Approach_Time) {prox.Min_Sink_Approach_Time = tSqr;}
    if(tff4 < prox.Min_Sink_Freefall_time) {prox.Min_Sink_Freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
    if(target.ptype == 5 && src.n_sink == 1) // only do it if we're looking at a single star in the node
    {
        double specific_energy = 0.5*vSqr - All.G*M_total/sqrt(r2);
        if (specific_energy<0)
        {
            double semimajor_axis= -All.G*M_total/(2.*specific_energy);
            double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
            if(t_orbital < prox.Min_Sink_OrbitalTime) /* Save parameters of companion */
            {
                prox.Min_Sink_OrbitalTime=t_orbital; prox.comp_Mass=src.sink_mass;
                prox.comp_dx = sink_dr; prox.comp_dv = sink_dv;
            }
        } /* specific_energy < 0 */
    } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
}
#endif /* SINK_CALC_DISTANCES */

#endif /* GRAVTREE_FORCE_KERNEL_H */
