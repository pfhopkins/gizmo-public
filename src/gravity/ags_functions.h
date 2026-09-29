/* ags_functions.h — GPU-callable AGS helpers (particle size/volume) and
 * AGS_ATOMIC_{ADD,STORE} dispatch macros for j-side writes in the
 * AGSForce GPU port.
 *
 * Get_Particle_Size_AGS / get_particle_volume_ags used to be host-only
 * with implicit reliance on the global `P`. They are now KOKKOS_INLINE_FUNCTION
 * and take P explicitly so both the CPU tree-walk and the GPU parallel_for
 * lambda can call them uniformly. Inside a GPU kernel the caller passes the
 * SharedSpace P mirror; inside CPU code the caller passes the global P.
 *
 * AGS_ATOMIC_{ADD,STORE} selects backend atomics at compile time. On the GPU
 * path (Kokkos+GPU) we use Kokkos::atomic_{add,store}. On the CPU
 * tree-walk path there is no macro that captures `#pragma omp atomic` (pragmas
 * can't live inside macros), so we provide a plain assignment and the CPU
 * caller is expected to wrap the statement with `#pragma omp atomic` directly.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef AGS_FUNCTIONS_H
#define AGS_FUNCTIONS_H

#include "../declarations/allvars.h"

/* KOKKOS_INLINE_FUNCTION falls back to `inline` outside GPU TUs. */
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif


#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE

/* Effective particle size based on AGS_KernelRadius and NumNgb. */
KOKKOS_INLINE_FUNCTION
double Get_Particle_Size_AGS_P(int i, const struct particle_data *P_arr)
{
#if (NUMDIMS == 1)
    return 2.00000 * P_arr[i].AGS_KernelRadius / P_arr[i].NumNgb;
#endif
#if (NUMDIMS == 2)
    return 1.77245 * P_arr[i].AGS_KernelRadius / P_arr[i].NumNgb;
#endif
#if (NUMDIMS == 3)
    return 1.61199 * P_arr[i].AGS_KernelRadius / P_arr[i].NumNgb;
#endif
}


/* Particle volume = L^NUMDIMS where L = Get_Particle_Size_AGS. */
KOKKOS_INLINE_FUNCTION
double get_particle_volume_ags_P(int j, const struct particle_data *P_arr)
{
    double L_j = Get_Particle_Size_AGS_P(j, P_arr);
#if (NUMDIMS==1)
    return L_j;
#elif (NUMDIMS==2)
    return L_j*L_j;
#else
    return L_j*L_j*L_j;
#endif
}

/* Whether ags_density must solve for this particle's KernelRadius.  Reads only. */
KOKKOS_INLINE_FUNCTION
int ags_density_solves_for_P(int i, const struct particle_data *P_arr)
{
    int default_to_return = 0; // default to not being active - needs to be pro-actively 'activated' by some physics
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    default_to_return = 1;
    if(!((1 << P_arr[i].Type) & (ADAPTIVE_GRAVSOFT_FORALL))) /* particle is NOT one of the designated 'adaptive' types */
    {
        default_to_return = 0;
    } else {default_to_return = 1;} /* particle is AGS-active */
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    if(P_arr[i].Type==0)
    {
        default_to_return = 0; // don't actually need to do the loop //
    }
#endif
#ifdef DM_SIDM
    if((1 << P_arr[i].Type) & (DM_SIDM)) {default_to_return = 1;}
#endif
#if defined(CBE_INTEGRATOR)
    if(CBE_INTEGRATOR_DOES_TYPE(P_arr[i].Type)) {default_to_return = 1;}
#elif defined(DM_FUZZY)
    if(P_arr[i].Type == 1) {default_to_return = 1;}
#endif
    if(P_arr[i].TimeBin < 0) {default_to_return = 0;} /* check our 'marker' for particles which have finished iterating to an KernelRadius solution (if they have, dont do them again) */
    return default_to_return;
}

/* ags_density_solves_for_P, which also sets AGS_KernelRadius and AGS_zeta for particles that do not
 * iterate -- a safety default, idempotent, and touching only this particle's own fields. */
KOKKOS_INLINE_FUNCTION
int ags_density_isactive_P(int i, struct particle_data *P_arr)
{
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(!((1 << P_arr[i].Type) & (ADAPTIVE_GRAVSOFT_FORALL))) /* particle is NOT one of the designated 'adaptive' types */
    {
        P_arr[i].AGS_KernelRadius = All.ForceSoftening[P_arr[i].Type];
        P_arr[i].AGS_zeta = 0;
    }
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    if(P_arr[i].Type==0) {P_arr[i].AGS_KernelRadius = P_arr[i].KernelRadius;} // gas sees gas, these are identical
#endif
    return ags_density_solves_for_P(i, P_arr);
}


/* maximum allowed softening */
KOKKOS_INLINE_FUNCTION
double ags_return_maxsoft_P(int i, const struct particle_data *P_arr)
{
    double maxsoft = All.MaxKernelRadius; // user-specified maximum: nothing is allowed to exceed this
#ifdef PMGRID /* Maximum allowed gravitational softening when using the TreePM method. The quantity is given in units of the scale used for the force split (PM_ASMTH) */
    maxsoft = DMIN(maxsoft, 1e3 * 0.5 * All.Asmth[0]); /* no more than 1/2 the size of the largest PM cell, times a 'safety factor' which can be pretty big */
#endif
#if (ADAPTIVE_GRAVSOFT_FORALL & 32) && defined(SINK_PARTICLES) && !defined(SINGLE_STAR_SINK_DYNAMICS)
    /* MaxAccretionRadius is defined in params.txt in PHYSICAL units. Enforce the
     * recommended-policy invariant "nothing exceeds MaxKernelRadius" via DMIN
     * (never > All.MaxKernelRadius, already the line-1 seed): keeps the Mode-B
     * per-type node band (capped at MaxKernelRadius) a valid upper bound on this
     * sink's AGS leaf reach. No effect where SinkMaxAccretionRadius < MaxKernelRadius
     * (the recommended + every-tested config). */
    if(P_arr[i].Type == 5) {maxsoft = DMIN(maxsoft, All.SinkMaxAccretionRadius / All.cf_atime);}
#endif
    return maxsoft;
}


/* minimum allowed softening */
KOKKOS_INLINE_FUNCTION
double ags_return_minsoft_P(int i, const struct particle_data *P_arr)
{
    double minsoft = All.ForceSoftening[P_arr[i].Type]; // this is the user-specified minimum
#if !defined(ADAPTIVE_GRAVSOFT_FORALL)
    minsoft = DMIN(All.MinKernelRadius, minsoft);
#endif
    return minsoft;
}

/* Fuzzy-DM drift/kick. Lives here rather than in sidm/dm_fuzzy_loop.h because that
 * header includes <Kokkos_Core.hpp> unconditionally and sidm/dm_fuzzy.cc is a host TU;
 * pulling Kokkos into a host TU breaks the CUDA build while compiling fine on the
 * OpenMP backend. Everything below touches only this particle's AGS_* fields.
 * DM_FUZZY implies AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE, so the enclosing gate holds. */
#ifdef DM_FUZZY
KOKKOS_INLINE_FUNCTION
void do_dm_fuzzy_drift_kick_P(int i, double dt, int mode, struct particle_data *pp)
{
    if(mode==0)
    {
        // calculate various energies: quantum potential QP0, 'stored' numerical pressure NQ0, kinetic energy KE0
        double dNQ=pp[i].AGS_Dt_Numerical_QuantumPotential*dt, NQ0=pp[i].AGS_Numerical_QuantumPotential, NQ1=NQ0+dNQ, KE0=0.5*pp[i].Mass*pp[i].Vel.norm_sq()*All.cf_a2inv;
        double f00 = 0.5 * All.ScalarField_hbar_over_mass; // this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass
        double d2rho = pp[i].AGS_Gradients2_Density[0][0] + pp[i].AGS_Gradients2_Density[1][1] + pp[i].AGS_Gradients2_Density[2][2]; // laplacian
        double drho2 = pp[i].AGS_Gradients_Density.norm_sq();
        double QP0 = (f00*f00 / pp[i].AGS_Density) * (d2rho - 0.5*drho2/pp[i].AGS_Density); // quantum 'potential'
        NQ1 = DMAX(0,DMAX(NQ1,0.1*NQ0)); NQ1 = DMIN(NQ1,1.1*DMAX(DMAX(KE0+NQ0,fabs(QP0)),KE0+NQ0+QP0)); // limit kick to not produce unphysical energy over-or-under-shoot
        pp[i].AGS_Numerical_QuantumPotential = NQ1;
    }

#if (DM_FUZZY > 0) /* if using direct-wavefunction integration methods */
    double vol_inv = pp[i].AGS_Density / pp[i].Mass;
    if(mode == 0)
    {
        //double psimag_mass_old = (pp[i].AGS_Psi_Re*pp[i].AGS_Psi_Re + pp[i].AGS_Psi_Im*pp[i].AGS_Psi_Im) * vol_inv;
        pp[i].AGS_Psi_Re += pp[i].AGS_Dt_Psi_Re * dt;
        pp[i].AGS_Psi_Im += pp[i].AGS_Dt_Psi_Im * dt;
        double mass_old = pp[i].Mass, dmass = pp[i].AGS_Dt_Psi_Mass * dt, mass_new = mass_old + dmass;
        dmass = DMIN(DMAX(dmass,-0.5*mass_old),0.5*mass_old);
        mass_new = mass_old + dmass;
        double psimag_mass_new = (pp[i].AGS_Psi_Re*pp[i].AGS_Psi_Re + pp[i].AGS_Psi_Im*pp[i].AGS_Psi_Im) * vol_inv;
#if (DM_FUZZY == 2)
        mass_new = psimag_mass_new; /* uses direct [NON-MASS-CONSERVING] integration of psi field */
#endif
        double psi_corr_fac = sqrt(mass_new / (MIN_REAL_NUMBER + psimag_mass_new));
        pp[i].Mass = mass_new; pp[i].AGS_Psi_Re *= psi_corr_fac; pp[i].AGS_Psi_Im *= psi_corr_fac;

        pp[i].AGS_Density = pp[i].Mass * vol_inv;
        pp[i].AGS_Psi_Re_Pred = pp[i].AGS_Psi_Re;
        pp[i].AGS_Psi_Im_Pred = pp[i].AGS_Psi_Im;
    } else {
        /* in drift mode, AGS_Density should automatically be drifted already by the predictor step, but not the other quantities here */
        pp[i].AGS_Psi_Re_Pred += pp[i].AGS_Dt_Psi_Re * dt;
        pp[i].AGS_Psi_Im_Pred += pp[i].AGS_Dt_Psi_Im * dt;
        pp[i].AGS_Density *= 1. + DMIN(DMAX(pp[i].AGS_Dt_Psi_Mass*dt/pp[i].Mass,-0.5),0.5);
    }
#endif
}
#endif /* DM_FUZZY */

#endif /* AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE */


/* Kokkos atomics inside a device lambda. */
#define AGS_ATOMIC_ADD(ptr, val)   Kokkos::atomic_add((ptr), (val))
#define AGS_ATOMIC_STORE(ptr, val) Kokkos::atomic_store((ptr), (val))

#endif /* AGS_FUNCTIONS_H */
