#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../system/eigen_symmetric.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../eos/cosmic_ray_fluid/cosmic_ray_functions.h"
#include "../mesh/kernel.h"
#include "../mesh/neighbor_list.h"
#include "../mesh/ghost_writeback.h"
extern void hydro_evaluate_gpu(struct particle_data *, struct gas_cell_data *,
                               int, int *, int, int64_t *, int *, int64_t, void *);
#include "../mesh/ghost_symlist_lifecycle.h"

/*! \file hydro_toplevel.c
 *  \brief This contains the "primary" hydro loop, where the hydro fluxes are computed.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* some very useful notes on the hydro variables in comoving integrations:

 v_code = a * v_peculiar/physical (canonical momentum)
 r_code = r_physical / a (comoving coordinates)
 m_code = m_physical
 rho_code = rho_physical * a^3 (from length/mass scaling)
 InternalEnergy_code = InternalEnergy_physical
 Pressure_code = InternalEnergy_code * rho_code * (gamma-1) = Pressure_physical * a^3
 B_code = a*a * B_physical (comoving magnetic fields)
 Phi_code = B_code*v_code = a^3 * Phi_physical (damping field for Dedner divergence cleaning)
    (note: spec egy of phi field is: phi*phi/(2*mu0*rho*ch*ch); compare Bfield is B*B/(mu0*rho);
    so [phi]~[B]*[ch], where ch is the signal velocity used in the damping equation);
 -- Time derivatives (rate of change from hydro forces) here are all
        assumed to end up in *physical* units ---
 HydroAccel, dMomentum are assumed to end up in *physical* units
    (note, this is different from GADGET's convention, where
     HydroAccel is in units of (Pcode/rhocode)/rcode)
 DtInternalEnergy and dInternalEnergy are assumed to end up in *physical* units
 DtMass and dMass are assumed to end up in *physical* units

 -----------------------------------------

 // All.cf_atime = a = 1/(1+z), the cosmological scale factor //
 All.cf_atime = All.Time;
 // All.cf_a2inv is just handy //
 All.cf_a2inv = 1 / (All.Time * All.Time);
 // All.cf_a3inv * Density_code = Density_physical //
 All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
 // time units: proper time dt_phys = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
 code time unit in comoving is dlna, so dt_phys = dt_code / All.cf_hubble_a   //
 All.cf_hubble_a = hubble_function(All.Time); // hubble_function(a) = H(a) = H(z) //
 // dt_code * v_code/r_code = All.cf_hubble_a2 * dt_phys * v_phys/r_phys //
 All.cf_hubble_a2 = All.Time * All.Time * hubble_function(All.Time);

 -----------------------------------------
 A REMINDER ABOUT GIZMO/GADGET VELOCITY UNITS:: (from Volker)

 The IC file should contain the *peculiar* velocity divided by sqrt(a),
 not the *physical* velocity. Let "x" denote comoving
 coordinates and "r=a*x" physical coordinates. Then I call

 comoving velocity: dx/dt
 physical velocity: dr/dt = H(a)*r + a*dx/dt
 peculiar velocity: v = a * dx/dt

 The physical velocity is hence the peculiar velocity plus the Hubble flow.

 The internal velocity variable is not given by dx/d(ln a). Rather, it is given by
 the canonical momentum p = a^2 * dx/dt.
 The IC-file and snapshot files of gadget/GIZMO don't
 contain the variable "p" directly because of historical reasons.
 Instead, they contain the velocity variable
 u = v/sqrt(a) = sqrt(a) * dx/dt = p / a^(3/2), which is just what the
 manual says. (The conversion between u and p is done on the fly when
 reading or writing snapshot files.)

 Also note that d(ln a)/dt is equal to the
 Hubble rate, i.e.: d(ln a)/dt = H(a) = H_0 * sqrt(omega_m/a^3 + omega_v
 + (1 - omega_m - omega_v)/a^2).

 Best wishes,
 Volker

 -----------------------------------------
*/


static double fac_mu, fac_vsic_fix;
/* fac_magnetic_pressure is lifted out of #ifdef MAGNETIC so per-pair hydro
   functions can take it as a plain argument; zero when MAGNETIC is off. */
static double fac_magnetic_pressure;


/* --------------------------------------------------------------------------------- */
/* define the kernel structure -- purely for handy purposes to clean up notation */
/* Shared struct definitions for hydro force (used by both CPU and GPU paths) */
#include "hydro_structs.h"
#ifndef HYDRO_SPH
#include "reimann.h"
#endif


/* GPU dispatcher consumes hydro_data_in/out directly from hydro_structs.h.
   The legacy export-buffer scaffolding (code_block_xchange_initialize.h, INPUT_STRUCT_NAME,
   OUTPUT_STRUCT_NAME, particle2in_hydra, hydro_force_evaluate) was retired in Step 5 Phase D1. */

/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_hydra(struct hydro_data_out *out, int i, int mode, int loop_iteration);
static inline void out2particle_hydra(struct hydro_data_out *out, int i, int mode, int loop_iteration)
{
    int k;
#if defined(MAGNETIC) && defined(MHD_BATTERY_MECHANISMS)
    {
        /* All cell-centered battery sources funnel into a single per-cell
           accumulator dBdt_battery_total, then a single shared limiter applies
           to the SUM (not per-source). Energy conservation is automatic via the
           per-cell line further below: DtInternalEnergy -= dot(B_phys, DtB).
           Each source is a CONSERVATIVE source, NOT a dissipative one, so we
           do NOT use the dot(B, bflux) heating hook (which is for Ohmic). */
        Vec3<double> dBdt_battery_total = {};

        /* Common unit conversion: dB/dt[G/s_phys] -> code-B per code-time, with
           cosmological a^2 lift (cf_a2inv = 1 non-cosmological). */
        const double dBdt_phys_to_code =
            UNIT_TIME_IN_CGS / All.UnitMagneticField_in_gauss / DMAX(All.cf_a2inv, MIN_REAL_NUMBER);

#if (MHD_BATTERY_MECHANISMS & 1)
        /* Tier-1 Biermann (Garaldi+2021 / AREPO style; see Soliman, Hopkins &
           Squire 2025). Computes dB/dt directly from slope-limited cell-centered
           gradients of n_e and T_e:
             dB/dt = (c k_B / (e n_e)) * (grad T_e) x (grad n_e)        [G/s]
           Curl reduces analytically to a clean two-gradient cross product
           because E_Bier is a pure gradient (curls of gradients vanish). */
        {
            const double n_e_cgs = CellP[i].n_e();
            const double T_e_cgs = CellP[i].T_e();
            if((n_e_cgs > 0) && (T_e_cgs > 0)) {
                const Vec3<double> g_Te_code = CellP[i].Gradients.ElectronTemperature;   /* K per code-length * cf_atime */
                const Vec3<double> g_ne_code = CellP[i].Gradients.ElectronNumberDensity; /* cm^-3 per code-length * cf_atime */

                /* physical gradients in cgs: code-coord conversion + cosmological
                   cf_atime strip (gradient pass stores cf_atime * physical) */
                const double inv_L_cgs = 1.0 / (UNIT_LENGTH_IN_CGS * All.cf_atime);
                const Vec3<double> g_Te_phys = g_Te_code * inv_L_cgs;   /* K/cm */
                const Vec3<double> g_ne_phys = g_ne_code * inv_L_cgs;   /* cm^-4 */

                const double prefac = C_LIGHT_CGS * BOLTZMANN_CGS
                    / (ELECTRONCHARGE_CGS * DMAX(n_e_cgs, MIN_REAL_NUMBER));
                const Vec3<double> dBdt_Bier_phys_cgs = prefac * cross(g_Te_phys, g_ne_phys);

                dBdt_battery_total += dBdt_Bier_phys_cgs * dBdt_phys_to_code;
            }
        }
#endif

#if (MHD_BATTERY_MECHANISMS & (2|4|8))
        /* Tier-2: radiative-ionization (bit 2) and dust battery (bits 4, 8).
           These EMFs do NOT reduce analytically to a clean cross product of
           primitive gradients (because their structure isn't pure-gradient
           like Biermann's). Instead, each per-cell builder writes its EMF
           into CellP[i].E_battery_T2_cell, the gradient pass takes its
           slope-limited gradient (a Mat3 tensor), and here we take the curl:
             dB/dt|_T2 = -c * curl(E_battery_T2)                       [phys]
           with the same physical->code unit conversion as Biermann.

           E_battery_T2_cell is in physical-cgs statvolt/cm (= Gauss).
           Gradients.E_battery_T2 stores the slope-limited gradient with
           the gradient pass's storage convention: per-code-length, with
           one factor of cf_atime injected (since x_phys = x_code * UnitLength
           * cf_atime). Strip both factors with inv_L_cgs to get a physical
           cgs gradient, then curl, then -c gives dB/dt in G/s_phys. */
        {
            const double inv_L_cgs = 1.0 / (UNIT_LENGTH_IN_CGS * All.cf_atime);
            Mat3<double> gradE_phys_cgs;
            for(int kr=0; kr<3; kr++) {
                for(int kc=0; kc<3; kc++) {
                    gradE_phys_cgs[kr][kc] = ((double)CellP[i].Gradients.E_battery_T2[kr][kc]) * inv_L_cgs;
                }
            }
            const Vec3<double> curl_E_phys_cgs = gradE_phys_cgs.curl();
            const Vec3<double> dBdt_T2_phys_cgs = -C_LIGHT_CGS * curl_E_phys_cgs;   /* G/s_phys */
            dBdt_battery_total += dBdt_T2_phys_cgs * dBdt_phys_to_code;
        }
#endif

        /* Single shared limiter on |dE_mag| per step over the SUM of all
           battery sources. GIZMO uses Heaviside-Lorentz internally for magnetic
           energy / pressure (no 4pi/8pi factors); see line below where
           DtInternalEnergy -= dot(B_phys, DtB). Compare:
             dE_mag_cell ~ |B_phys . dB| * V_code
           against eps * E_internal AND eps * (P_thermal + P_mag) * V_code.
           Slope-limited gradients already bound small-scale noise; this
           limiter protects against extrapolation past linearization regime
           (e.g. when growing B approaches equipartition). */
        const double dt_code   = get_particle_timestep_in_physical(i);
        const double V_code    = P[i].Mass / DMAX(CellP[i].Density, MIN_REAL_NUMBER);
        const Vec3<double> B_phys_codeunits = CellP[i].Bfield() * All.cf_a2inv;
        const Vec3<double> dB_step          = dBdt_battery_total * dt_code;
        const double dEmag_cell = fabs(dot(B_phys_codeunits, dB_step)) * V_code;

        const double E_internal_cell = P[i].Mass * CellP[i].InternalEnergyPred;
        const double P_thermal_phys  = CellP[i].Pressure * All.cf_a3inv;
        const double P_mag_phys      = 0.5 * B_phys_codeunits.norm_sq();
        const double E_pressure_cell = (P_thermal_phys + P_mag_phys) * V_code;

        const double eps = 0.1;
        const double allowed = eps * DMIN(DMAX(E_internal_cell, MIN_REAL_NUMBER),
                                          DMAX(E_pressure_cell, MIN_REAL_NUMBER));
        if((dEmag_cell > allowed) && (dEmag_cell > 0)) {
            const double scale = allowed / dEmag_cell;
            dBdt_battery_total *= scale;
        }

        out->DtB += dBdt_battery_total;
    }
#endif
    /* these are zero-d out at beginning of hydro loop so should always be added */
    CellP[i].HydroAccel += out->Acc;
    CellP[i].DtInternalEnergy += out->DtInternalEnergy;
#if defined(TWO_TEMPERATURE_PLASMA) && (TWO_TEMPERATURE_PLASMA & 4) && defined(CONDUCTION)
    CellP[i].DtInternalEnergy_FromConduction += out->DtInternalEnergy_FromConduction;
#endif
    //CellP[i].dInternalEnergy += out->dInternalEnergy; //manifest-indiv-timestep-debug//

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    CellP[i].DtMass += out->DtMass;
    CellP[i].dMass += out->dMass;
    CellP[i].GravWorkTerm += out->GravWorkTerm;
#endif
    if(CellP[i].MaxSignalVel < out->MaxSignalVel) {CellP[i].MaxSignalVel = out->MaxSignalVel;}
#ifdef OUTPUT_SHOCK_MACH_NUMBER
    if(CellP[i].ShockMachNumber < out->MaxShockMachNumber) {CellP[i].ShockMachNumber = out->MaxShockMachNumber;}
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    if(CellP[i].MaxKineticEnergyNgb < out->MaxKineticEnergyNgb) {CellP[i].MaxKineticEnergyNgb = out->MaxKineticEnergyNgb;}
#endif
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    for(k=0;k<NUM_METAL_SPECIES;k++) {CellP[i].Dyield[k] += out->Dyield[k];}
#endif

#ifdef CHIMES_TURB_DIFF_IONS
    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++)
      CellP[i].ChimesNIons[k] = DMAX(CellP[i].ChimesNIons[k] + out->ChimesIonsYield[k], 0.5 * CellP[i].ChimesNIons[k]);
#endif

#if defined(RT_SOLVER_EXPLICIT)
#if defined(RT_EVOLVE_ENERGY)
    for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_E_gamma[k] += out->Dt_Rad_E_gamma[k];}
#endif
#if defined(RT_EVOLVE_FLUX)
    for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_Flux[k] += out->Dt_Rad_Flux[k];}
#endif
#if defined(RT_INFRARED)
    CellP[i].Dt_Rad_E_gamma_T_weighted_IR += out->Dt_Rad_E_gamma_T_weighted_IR;
#endif
#if defined(RT_EVOLVE_INTENSITIES)
    for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {CellP[i].Dt_Rad_Intensity[k][k_dir] += out->Dt_Rad_Intensity[k][k_dir];}}
#endif
#endif

#if defined(MAGNETIC)
    /* can't just do DtB += out-> DtB, because for some hydro methods, the induction equation is solved in the density loop; need to simply add it here */
    CellP[i].DtB += out->DtB; CellP[i].Face_Area += out->Face_Area;
    CellP[i].divB += out->divB;
#if defined(DIVBCLEANING_DEDNER)
#ifdef HYDRO_MESHLESS_FINITE_VOLUME // mass-based phi-flux
    CellP[i].DtPhi += out->DtPhi;
#endif
    CellP[i].DtB_PhiCorr += out->DtB_PhiCorr;
#endif // Dedner //
#endif // MAGNETIC //

#ifdef COSMIC_RAY_FLUID
    CellP[i].Face_DivVel_ForAdOps += out->Face_DivVel_ForAdOps;
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
    CellP[i].DtCREgyNewInjectionFromShocks += out->DtCREgyNewInjectionFromShocks;
#endif
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        CellP[i].DtCosmicRayEnergy[k] += out->DtCosmicRayEnergy[k];
#if defined(CRFLUID_EVOLVE_SPECTRUM)
        CellP[i].DtCosmicRay_Number_in_Bin[k] += out->DtCosmicRay_Number_in_Bin[k];
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
        int kAlf; for(kAlf=0;kAlf<2;kAlf++) {CellP[i].DtCosmicRayAlfvenEnergy[k][kAlf] += out->DtCosmicRayAlfvenEnergy[k][kAlf];}
#endif
    }
#endif
}


/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/* This will perform final operations and corrections on the output from the
    hydro routines, AFTER the neighbors have all been checked and summed */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
void hydro_final_operations_and_cleanup(void)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++)
    { int i = ActiveParticleList[_apl]; int k = 0; /* k is declared here (not at function scope) so it is thread-private under OpenMP — a shared k would race between threads and corrupt CR bin index */
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            double dt; dt = get_particle_timestep_in_physical(i);

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* signal velocity needs to include rate of gas flow -over- the resolution element, which can be non-zero here */
            double v2_p = CellP[i].MaxSignalVel*CellP[i].MaxSignalVel + (CellP[i].VelPred - CellP[i].ParticleVel).norm_sq();
            CellP[i].MaxSignalVel = sqrt(v2_p);
#endif

#if defined(MAGNETIC)
            /* need to subtract out the source terms proportional to the (non-zero) B-field divergence; to stabilize the scheme */
            {
                Vec3<double> Bi = CellP[i].Bfield() * All.cf_a2inv;
#ifndef HYDRO_SPH
                /* this part of the induction equation has to do with advection of div-B, it is not present in SPH */
                CellP[i].DtB -= CellP[i].divB * (CellP[i].VelPred/All.cf_atime);
#endif
                CellP[i].HydroAccel -= CellP[i].divB * Bi;
                CellP[i].DtInternalEnergy -= CellP[i].divB * dot(CellP[i].VelPred/All.cf_atime, Bi);
            }

            double magnorm_closure = Get_DtB_FaceArea_Limiter(i, P, CellP);

#if defined(DIVBCLEANING_DEDNER) && !defined(HYDRO_SPH)
            // ok now deal with the divB correction forces and damping fields //
            double tolerance_for_correction,db_vsig_h_norm;
            tolerance_for_correction = 10.0;
            db_vsig_h_norm = 0.1; // can be as low as 0.03 //
            double DtB_PhiCorr=0,DtB_UnCorr=0,db_vsig_h=0,PhiCorr_Norm=1.0;
            for(k=0; k<3; k++)
            {
                DtB_UnCorr += CellP[i].DtB[k] * CellP[i].DtB[k]; // physical units //
                db_vsig_h = db_vsig_h_norm * (CellP[i].BPred[k]*All.cf_atime) * (0.5*CellP[i].MaxSignalVel) / (P[i].Get_Particle_Size()*All.cf_atime);
                DtB_UnCorr += db_vsig_h * db_vsig_h;
                DtB_PhiCorr += CellP[i].DtB_PhiCorr[k] * CellP[i].DtB_PhiCorr[k];
            }

            /* take a high power of these: here we'll use 4, so it works like a threshold */
            DtB_UnCorr*=DtB_UnCorr; DtB_PhiCorr*=DtB_PhiCorr; tolerance_for_correction *= tolerance_for_correction;
            /* now re-normalize the correction term if its unacceptably large */
            if((DtB_PhiCorr > 0)&&(!isnan(DtB_PhiCorr))&&(DtB_UnCorr>0)&&(!isnan(DtB_UnCorr))&&(tolerance_for_correction>0)&&(!isnan(tolerance_for_correction)))
            {

                if(DtB_PhiCorr > tolerance_for_correction * DtB_UnCorr) {PhiCorr_Norm *= tolerance_for_correction * DtB_UnCorr / DtB_PhiCorr;}
                CellP[i].DtB += PhiCorr_Norm * CellP[i].DtB_PhiCorr;
                CellP[i].DtInternalEnergy += PhiCorr_Norm * dot(CellP[i].DtB_PhiCorr, CellP[i].Bfield() * All.cf_a2inv);
            }

#ifdef HYDRO_MESHLESS_FINITE_VOLUME // mass-based phi-flux
            CellP[i].DtPhi *= magnorm_closure;
#else
            CellP[i].DtPhi = 0;
#endif
            if((!isnan(CellP[i].divB))&&(P[i].KernelRadius>0)&&(CellP[i].divB!=0)&&(CellP[i].Density>0))
            {
                double tmp_ded = 0.5 * CellP[i].MaxSignalVel / (fac_mu*All.cf_atime); // has units of v_physical now
                /* do a check to make sure divB isn't something wildly divergent (owing to particles being too close) */
                Vec3<double> Bfield_i = CellP[i].BPred * (CellP[i].Density / P[i].Mass);
                double b2_max = Bfield_i.norm_sq();
                b2_max = 100.0 * fabs( sqrt(b2_max) * All.cf_a2inv * P[i].Mass / (CellP[i].Density*All.cf_a3inv) * 1.0 / (P[i].KernelRadius*All.cf_atime) );
                if(fabs(CellP[i].divB) > b2_max) {CellP[i].divB *= b2_max / fabs(CellP[i].divB);}
                /* ok now can apply this to get the growth rate of phi */
                // CellP[i].DtPhi -= tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * CellP[i].divB;
                CellP[i].DtPhi -= tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * CellP[i].divB * CellP[i].Density*All.cf_a3inv; // mass-based phi-flux
            }
#endif
#endif // MAGNETIC


            /* we calculated the flux of conserved variables: these are used in the kick operation. But for
             intermediate drift operations, we need the primive variables, so reduce to those here
             (remembering that v_phys = v_code/All.cf_atime, for the sake of doing the unit conversions to physical) */
            {Vec3<double> vel_phys = CellP[i].VelPred / All.cf_atime;
            CellP[i].DtInternalEnergy -= dot(vel_phys, CellP[i].HydroAccel);
            /* we solved for total energy flux (and remember, HydroAccel is still momentum -- keep units straight here!) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            CellP[i].DtInternalEnergy += 0.5 * vel_phys.norm_sq() * CellP[i].DtMass;
            CellP[i].HydroAccel -= vel_phys * CellP[i].DtMass; /* we solved for momentum flux */
#endif
            CellP[i].HydroAccel /= P[i].Mass;} /* we solved for momentum flux */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            CellP[i].DtInternalEnergy -= CellP[i].InternalEnergyPred * CellP[i].DtMass;
#endif
#ifdef MAGNETIC
#ifndef HYDRO_SPH
            CellP[i].DtInternalEnergy -= dot(CellP[i].Bfield() * All.cf_a2inv, CellP[i].DtB);
#endif
            CellP[i].DtB *= magnorm_closure;
#endif
            CellP[i].DtInternalEnergy /= P[i].Mass;
            /* ok, now: HydroAccel = dv/dt, DtInternalEnergy = du/dt (energy per unit mass) */

            /* zero out hydrodynamic PdV work terms if the particle is at the maximum smoothing, these will be incorrect */
            if(P[i].KernelRadius >= 0.99*All.MaxKernelRadius) {CellP[i].DtInternalEnergy = 0;}

            // need to explicitly include adiabatic correction from the hubble-flow (for drifting) here //
            if(All.ComovingIntegrationOn) {CellP[i].DtInternalEnergy -= 3*(CellP[i].gamma_eos_value()-1) * CellP[i].InternalEnergyPred * All.cf_hubble_a;}
            // = du/dlna -3*(gamma-1)*u ; then dlna/dt = H(z) =  All.cf_hubble_a //


#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_EVOLVE_FLUX) && !defined(RT_RADPRESSURE_IN_HYDRO) /* // -- moved for OTVET+FLD to drift-kick operation to deal with limiters more accurately -- // */
            /* calculate the radiation pressure force */
            Vec3<double> radacc = {}; int kfreq;
            for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++)
            {
                double vol_inv = CellP[i].Density*All.cf_a3inv/P[i].Mass, f_kappa_abs = rt_absorb_frac_albedo(i,kfreq, P, CellP), flux_mag=0, erad_i=0, flux_corr=0, work_band=0, rmag=0;
                Vec3<double> vel_i = {}, vdot_h = {}, flux_i = {}, radacc_thisband = {};
                erad_i = CellP[i].Rad_E_gamma_Pred[kfreq]*vol_inv; // radiation energy density, needed below
                flux_i = CellP[i].Rad_Flux_Pred[kfreq] * vol_inv; vel_i = CellP[i].VelPred * ((C_LIGHT_CODE_REDUCED/C_LIGHT_CODE)/All.cf_atime); flux_mag = flux_i.norm_sq();
                vdot_h = (RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*C_LIGHT_CODE/C_LIGHT_CODE_REDUCED) * erad_i * (vel_i + CellP[i].ET[kfreq].matvec(vel_i)); // P_rad term + eI term //
                double flux_thin = erad_i * C_LIGHT_CODE_REDUCED; if(flux_mag>0) {flux_mag=sqrt(flux_mag);} else {flux_mag=1.e-20*flux_thin;}
                if(flux_mag > 0) {flux_corr = DMIN(1., flux_thin/flux_mag); // restrict flux here (b/c drifted can exceed physical b/c of integration errors
#if defined(RT_ENABLE_R15_GRADIENTFIX)
                    flux_corr = flux_thin/flux_mag; // set to maximum (optically thin limit)
#endif
                }
                double L_particle=P[i].Get_Particle_Size()*All.cf_atime, Sigma_particle=P[i].Mass/(M_PI*L_particle*L_particle), abs_per_kappa_dt=C_LIGHT_CODE_REDUCED*(CellP[i].Density*All.cf_a3inv)*dt; // effective surface density through particle & fractional absorption over timestep
                int checker_int = 0; // normal default: only use the corrections below for bands which dont re-emit to the same band
                checker_int = 1; // actually here and above now changed to use the slabfac corrections for all bands. in the resolved limit this should still be correct because the re-emitted photons should be isotropic: otherwise you run into linear momentum conservation problems. this is only an issue if the source is at the center of the distribution.
                double slabfac_rp=1; if(check_if_absorbed_photons_can_be_reemitted_into_same_band(kfreq)<=checker_int) {slabfac_rp=slab_averaging_function(f_kappa_abs*CellP[i].Rad_Kappa[kfreq]*Sigma_particle) * slab_averaging_function(f_kappa_abs*CellP[i].Rad_Kappa[kfreq]*abs_per_kappa_dt);} // reduction factor for absorption over dt
                radacc_thisband = (slabfac_rp * (CellP[i].Rad_Kappa[kfreq]/C_LIGHT_CODE_REDUCED)) * (flux_corr * flux_i - vdot_h); rmag = radacc_thisband.norm_sq(); // acceleration term before accounting for the 'work' term, which is calculated separately in the absorption/emission loop
                if(check_if_absorbed_photons_can_be_reemitted_into_same_band(kfreq)<=checker_int && f_kappa_abs > MIN_REAL_NUMBER && rmag > MIN_REAL_NUMBER && dt > 0 && P[i].Mass > 0) { // bands that destroy photons upon absorption (e.g. ionization, dust absorption) should limit the imparted momentum to the total photon momentum available - the flux in the solver normally prevents this but this addresses some edge cases with e.g. pathological ICs, rapidly-varying kappa, etc.
                    rmag=sqrt(rmag); double r_from_abs=f_kappa_abs*rmag, abs_dt=rt_absorption_rate(i,kfreq, P, CellP)*dt, dE_abs=erad_i*(1.-exp(-abs_dt)); if(abs_dt<0.01) {dE_abs=erad_i*abs_dt;}
                    double rmag_max_abs=dE_abs/(vol_inv*P[i].Mass*C_LIGHT_CODE_REDUCED*dt); if(rmag_max_abs<r_from_abs) {double cfac=1.+(rmag_max_abs-r_from_abs)/rmag; if(cfac>0 && cfac<1) {radacc_thisband *= cfac;}}
                }
                radacc += radacc_thisband; work_band = dot(radacc_thisband, vel_i) * P[i].Mass; // PdV work done by photons [absorbed ones are fully-destroyed, so their loss of energy and momentum is already accounted for by their deletion in this limit -- note that we have to be careful about the RSOL factors here! //
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES
                f_kappa_abs = 0;
#endif
                CellP[i].Dt_Rad_E_gamma[kfreq] += (2.*f_kappa_abs-1.)*work_band; // loss/gain term for the radiation field itself
#ifdef TRANSPORT_SUBCYCLE
                CellP[i].Dt_Rad_E_gamma_Work[kfreq] += (2.*f_kappa_abs-1.)*work_band; // save work term separately for subcycle restoration
#endif
                double du_radwork = (C_LIGHT_CODE/C_LIGHT_CODE_REDUCED) * 2.*f_kappa_abs*work_band / P[i].Mass;
                CellP[i].DtInternalEnergy -= du_radwork; // correct for rsol factor above which reduced vel_i by rsol; -only- add back this term for gas
                /* HYDRO_RADWORK_DIAG: trace radiation work contribution to DtInternalEnergy */
#ifdef GIZMO_DEBUG_RT_COOLING
                {static int hrd_n=0; if(hrd_n < 20 && (P[i].ID == 1 || P[i].ID == 100 || P[i].ID == 1000)) {
                    printf("[HYDRO_RADWORK] ID=%llu k=%d DtU_before=%.10e du_radwork=%.10e f_kappa=%.6e work_band=%.6e Rad_Kappa=%.6e flux_mag=%.6e erad=%.6e\n",
                        (unsigned long long)P[i].ID, kfreq, CellP[i].DtInternalEnergy+du_radwork, du_radwork, f_kappa_abs, work_band, CellP[i].Rad_Kappa[kfreq], flux_mag, erad_i);
                    hrd_n++;}}
#endif
            }
            /* now actually set the frequency-integrated cell values as needed */
#ifdef RT_RAD_PRESSURE_OUTPUT
            CellP[i].Rad_Accel = radacc; // physical units, as desired
#else
            CellP[i].HydroAccel += radacc; // physical units, as desired
#endif
#endif
            /* HYDRO_DTU_DIAG: print final DtInternalEnergy after all corrections — by particle ID */
#ifdef GIZMO_DEBUG_RT_COOLING
            {static int hdtu_n=0; if(hdtu_n < 15 && (P[i].ID == 1 || P[i].ID == 100 || P[i].ID == 1000)) {
                printf("[HYDRO_DTU] ID=%llu DtU_final=%.10e u=%.10e\n", (unsigned long long)P[i].ID, CellP[i].DtInternalEnergy, CellP[i].InternalEnergyPred);
                hdtu_n++;}}
#endif /* GIZMO_DEBUG_RT_COOLING */
#ifdef RT_RADPRESSURE_IN_HYDRO
            int kfreq; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {
                double fac = (1./3.) * CellP[i].flux_limiter(kfreq) * CellP[i].Rad_E_gamma_Pred[kfreq] * P[i].Particle_DivVel*All.cf_a2inv * (1.-2.*rt_absorb_frac_albedo(i, kfreq, P, CellP));
                CellP[i].Dt_Rad_E_gamma[kfreq] -= (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * fac; CellP[i].DtInternalEnergy += fac / P[i].Mass; /* exact energy conservation; for appropriate RSOL definitions - careful of terms here where beta arises */
            }
#endif


#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME)) /* update the metal masses from exchange */
            for(k=0;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] = DMAX(P[i].Metallicity[k] + CellP[i].Dyield[k] / P[i].Mass , 0.01*P[i].Metallicity[k]);}
#if defined(GALSF_ISMDUSTCHEM_MODEL) /* update the dust masses from exchange: offsets within Dyield match Metallicity layout */
            {int dc0 = ISMDUSTCHEM_SPECIES_OFFSET_IN_METALLICITY; /* base offset for dustchem in Metallicity/Dyield */
            for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {CellP[i].ISMDustChem_Dust_Metal[k] = DMAX(CellP[i].ISMDustChem_Dust_Metal[k] + CellP[i].Dyield[dc0+k] / P[i].Mass , 0.01*CellP[i].ISMDustChem_Dust_Metal[k]);}
            for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {CellP[i].ISMDustChem_Dust_Source[k] = DMAX(CellP[i].ISMDustChem_Dust_Source[k] + CellP[i].Dyield[dc0+NUM_ISMDUSTCHEM_ELEMENTS+k] / P[i].Mass , 0.01*CellP[i].ISMDustChem_Dust_Source[k]);}
            for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {CellP[i].ISMDustChem_Dust_Species[k] = DMAX(CellP[i].ISMDustChem_Dust_Species[k] + CellP[i].Dyield[dc0+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+k] / P[i].Mass , 0.01*CellP[i].ISMDustChem_Dust_Species[k]);}
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
            for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {
                int l;
                double new_bin_mass, new_bin_number, old_bin_mass;
                for(l=0;l<NUM_ISMDUSTCHEM_SIZE_BINS;l++) {
                    new_bin_number = DMAX(CellP[i].ISMDustChem_Dust_NumberInBin[k][l] + CellP[i].Dyield[dc0+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES+(k*NUM_ISMDUSTCHEM_SIZE_BINS+l)], 0.01*CellP[i].ISMDustChem_Dust_NumberInBin[k][l]);
                    old_bin_mass = get_ISMDustChemEvo_bin_mass(i,k,l, CellP);
                    new_bin_mass = DMAX(old_bin_mass + CellP[i].Dyield[dc0+NUM_ISMDUSTCHEM_ELEMENTS+NUM_ISMDUSTCHEM_SOURCES+NUM_ISMDUSTCHEM_SPECIES+(NUM_ISMDUSTCHEM_SPECIES*NUM_ISMDUSTCHEM_SIZE_BINS)+(k*NUM_ISMDUSTCHEM_SIZE_BINS+l)] * UNIT_MASS_IN_CGS, 0.01*old_bin_mass);
                    update_ISMDustChemEvo_bin_number_and_slope(i,k,l,new_bin_number,new_bin_mass, CellP);
                }
            }
#endif
            } /* close dc0 scope */
#endif
#endif
            
            
#if (defined(COSMIC_RAY_FLUID) && !defined(COOLING_OPERATOR_SPLIT)) || defined(COSMIC_RAY_SUBGRID_LEBRON)
            /* with the spectrum model, we account here the adiabatic heating/cooling of the 'fluid', here, which was solved in the hydro solver but doesn't resolve which portion goes to CRs and which to internal energy, with gamma=GAMMA_COSMICRAY */
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            double P_cr_spec = (1./3.)*CellP[i].SubGrid_CosmicRayEnergyDensity/CellP[i].Density, P_tot_spec = P_cr_spec + (2./3.)*CellP[i].InternalEnergyPred + (1./2.)*pow(CellP[i].Alfven_speed(),2); // just include CR+thermal+magnetic here
            CellP[i].DtInternalEnergy *= (1.-P_cr_spec/P_tot_spec); /* approximate correction, valid to level here [more sophisticated correction can cause problems since the PdV energy isn't actually being taken -out- of the CR field, as it would be if followed explicitly] */
#else
            double gamma_minus_eCR_tmp=0; for(k=0;k<N_CR_PARTICLE_BINS;k++) {gamma_minus_eCR_tmp+=(GAMMA_COSMICRAY(k)-1.)*CellP[i].CosmicRayEnergyPred[k];} // routine below only depends on the total CR energy, not bin-by-bin energies, when we do it this way here
            double dCR_div = CR_calculate_adiabatic_gasCR_exchange_term(i, dt, gamma_minus_eCR_tmp, 1, P, CellP); // this will handle the update below - separate subroutine b/c we want to allow it to appear in a couple different places
            double u0=DMAX(CellP[i].InternalEnergyPred, All.MinEgySpec) , uf=DMAX(u0 - dCR_div/P[i].Mass , All.MinEgySpec); // final updated value of internal energy per above
            CellP[i].DtInternalEnergy += (uf - u0) / (dt + MIN_REAL_NUMBER); // update gas quantities to be used in cooling function
#endif
#endif
#if defined(COSMIC_RAY_FLUID)
            /* energy transfer from CRs to gas due to the streaming instability (mediated by high-frequency Alfven waves, but they thermalize quickly
                (note this is important; otherwise build up CR 'traps' where the gas piles up and cools but is entirely supported by CRs in outer disks) */
#if !defined(CRFLUID_EVOLVE_SCATTERINGWAVES) // handled in separate solver if explicitly evolving the relevant wave families
            for(k=0;k<N_CR_PARTICLE_BINS;k++) {
                double streamfac = fabs(CR_get_streaming_loss_rate_coefficient(i,k, P, CellP));
                CellP[i].DtInternalEnergy += CellP[i].CosmicRayEnergyPred[k] * streamfac / P[i].Mass; // make sure to divide by mass here to get the correct units since DtInternalEnergy has been converted to specific energy units (while CR energies are absolute)
#if !defined(CRFLUID_EVOLVE_SPECTRUM)
                CellP[i].DtCosmicRayEnergy[k] -= cosmicrayfluid_rsol_corrfac(k) * CellP[i].CosmicRayEnergyPred[k] * streamfac; // in the multi-bin formalism, save this operation for the CR cooling ops since can involve bin-to-bin transfer of energy
#endif
            }
#endif
#if defined(MAGNETIC) // only makes sense to include parallel correction below if all these terms enabled //
            /* 'residual' term from parallel scattering of CRs being not-necessarily-in-equilibrium with a two-moment form of the equations */
            double vA_eff=Get_Gas_ion_Alfven_speed_i(i, P, CellP), vol_i=CellP[i].Density*All.cf_a3inv/P[i].Mass, Bmag=0; Vec3<double> bhat = CellP[i].BPred; // define some useful variables
            Bmag = bhat.norm_sq(); // get direction vector for B-field needed below
            if(Bmag>0) {Bmag=sqrt(Bmag); bhat /= Bmag;} // make dimensionless
            if(Bmag>0) {for(k=0;k<N_CR_PARTICLE_BINS;k++) {
                int target_for_cr_betagamma = i; // if this = -1, use the gamma factor at the bin-center for evaluating this, if this = i, use the mean gamma of the bin, weighted by the CR energy -- won't give exactly the same result here
                target_for_cr_betagamma = -1; // the correction terms depend on these being evaluated at their bin-centered locations
                double three_chi = return_cosmic_ray_anisotropic_closure_function_threechi(i,k, CellP);
                double grad_P_dot_B=0, F_dot_B=0, e0_cr=CellP[i].CosmicRayEnergyPred[k]*vol_i, p0_cr=(GAMMA_COSMICRAY(k)-1.)*e0_cr, vA_k=vA_eff*return_CRbin_nuplusminus_asymmetry(i, k, CellP), beta_fac=return_CRbin_beta_factor(target_for_cr_betagamma,k,CellP);
                Vec3<double> gradpcr = CellP[i].Gradients.CosmicRayPressure[k] * (All.cf_a3inv/All.cf_atime);
                grad_P_dot_B = dot(bhat, gradpcr); F_dot_B = dot(bhat, CellP[i].CosmicRayFluxPred[k]) * vol_i;
                if(F_dot_B < 0) {vA_k *= -1;} // needs to have appropriately-matched signage below //
                double gamma_0=return_CRbin_gamma_factor(target_for_cr_betagamma,k,CellP), gamma_fac=gamma_0/(gamma_0-1.); // lorentz factor here, needed in next line, because the loss term here scales with -total- energy, not kinetic energy
                if(beta_fac<0.1) {gamma_fac=2./(beta_fac*beta_fac) -0.5 - 0.125*beta_fac*beta_fac;} // avoid accidental nan
                Vec3<double> fcorr = bhat * (grad_P_dot_B + (gamma_fac*(F_dot_B/cosmicrayfluid_rsol_corrfac(k)) - three_chi*vA_k*(gamma_fac*e0_cr + p0_cr))*(beta_fac*beta_fac)/(3.*CellP[i].CosmicRayDiffusionCoeff[k])) / (CellP[i].Density*All.cf_a3inv); // physical units
                fcorr += (1.-three_chi) * (gradpcr - bhat*grad_P_dot_B) / (CellP[i].Density*All.cf_a3inv); // physical units
                CellP[i].HydroAccel += fcorr; // add correction term back into hydro acceleration terms -- need to check that don't end up with nasty terms for badly-initialized/limited scattering rates above
            }}
#endif
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
            if((CellP[i].DtCREgyNewInjectionFromShocks <= 0) || (CellP[i].DtInternalEnergy <= 0) || (P[i].Mass <= 0)) {CellP[i].DtCREgyNewInjectionFromShocks = 0;} // should never be negative, thats an error from above, or a spurious shock detection if dtinternalenergy summed is negative - don't inject here (additional useful check over kernel)
            if(CellP[i].DtCREgyNewInjectionFromShocks > 0) { // do some checks and adjust the internal energy evolution to ensure total energy conservation now
                double dtThermal = CellP[i].DtInternalEnergy * P[i].Mass; // correct back to total energy units
                CellP[i].DtCREgyNewInjectionFromShocks = DMIN(CellP[i].DtCREgyNewInjectionFromShocks, 0.5*dtThermal); // don't allow more than 1/2 of the total internal energy change to go into CRs (should usually be satisfied but since cr injection pairwise, can be violated if shock not well-resolved and flow is chaotic)
                CellP[i].DtInternalEnergy = (dtThermal - CellP[i].DtCREgyNewInjectionFromShocks) / P[i].Mass; // reset the thermal energy budget appropriately, now total energy will be conserved (just locally shifting from one reservoir to another)
            }
#endif
#endif // COSMIC_RAY_FLUID


#ifdef GALSF_SUBGRID_WINDS
            /* if we have winds, we decouple particles briefly if delaytime>0 */
            if(CellP[i].DelayTime > 0)
            {
                CellP[i].HydroAccel = {};
                CellP[i].DtInternalEnergy = 0; //CellP[i].dInternalEnergy = 0;
                double windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
                windspeed *= fac_mu;
                double rkern_c = pow(All.WindFreeTravelDensFac * All.PhysDensThresh / (CellP[i].Density * All.cf_a3inv), (1. / 3.));
                CellP[i].MaxSignalVel = rkern_c * DMAX((2 * windspeed), CellP[i].MaxSignalVel);
            }
#endif


#ifdef BOX_BND_PARTICLES
            /* this flag signals all particles with id=0 are frozen (boundary particles) */
            if(P[i].ID == 0)
            {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                CellP[i].DtMass = 0;
                CellP[i].dMass = 0;
                CellP[i].GravWorkTerm = {};
#endif
                CellP[i].DtInternalEnergy = 0;//CellP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
                CellP[i].HydroAccel = {};
#ifdef MAGNETIC
                CellP[i].DtB = {};
#ifdef DIVBCLEANING_DEDNER
                CellP[i].DtB_PhiCorr = {};
                CellP[i].DtPhi = 0;
#endif
#endif
#ifdef SPH_BND_BFLD
                CellP[i].B = {};
#endif
            }
#endif

        } // closes P[i].Type==0 check and so closes loop over particles i
    } // for (loop over active particles) //

    
#ifdef TURB_DRIVING
#ifdef TURB_DRIVING_UPDATE_FORCE_ON_TURBUPDATE // if this is enabled, we only update as frequently as the driving phases are recomputed, as set by TurbDrive_TimeBetweenTurbUpdates. To avoid large errors, must be set by-hand to be << lambda_min / V where V is the typical turbulent velocity and lambda_min is the smallest driven wavelength.
    if(new_turbforce_needed_this_timestep()){add_turb_accel();}
#else    
    add_turb_accel(); // update turbulent driving fields and TurbAccel fields at same time as update HydroAccel, here
#endif    
#endif

}




/* this function exists to loop over the hydro variables and do any needed 'pre-processing' before they enter the primary hydro force loop */
void hydro_force_initial_operations_preloop(void)
{
    // Set global factors for comoving integration of hydro //
    fac_mu = 1 / (All.cf_atime); // code_vel * fac_mu = sqrt[code_pressure/code_density] = code_soundspeed //
    fac_vsic_fix = All.cf_hubble_a ; // note also that signal_vel in forms below should be in units of code_soundspeed //
#ifdef MAGNETIC
    fac_magnetic_pressure = 1. / All.cf_atime; // code_Bfield*code_Bfield * fac_magnetic_pressure = code_pressure -- use this to get alfven velocities, etc, as well as comoving units for magnetic integration //
#endif

    /* need to zero out all numbers that can be set -EITHER- by an active particle in the domain, or by one of the neighbors we will get sent */
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i)
#endif
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++)
    {   int i = ActiveParticleList[_apl]; int k = 0; /* thread-private k — same OMP race fix as hydro_final_operations_and_cleanup */
        if(P[i].Type==0)
        {
            CellP[i].MaxSignalVel = MIN_REAL_NUMBER;
#ifdef OUTPUT_SHOCK_MACH_NUMBER
            CellP[i].ShockMachNumber = 0;
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            CellP[i].MaxKineticEnergyNgb = MIN_REAL_NUMBER;
#endif
            CellP[i].DtInternalEnergy = 0; //CellP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
#if defined(TWO_TEMPERATURE_PLASMA) && (TWO_TEMPERATURE_PLASMA & 4) && defined(CONDUCTION)
            CellP[i].DtInternalEnergy_FromConduction = 0;
#endif
            CellP[i].HydroAccel = {};
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            CellP[i].DtMass = 0; CellP[i].dMass = 0; CellP[i].GravWorkTerm = {};
#endif
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
            for(k=0;k<NUM_METAL_SPECIES;k++) {CellP[i].Dyield[k] = 0;}
#endif
#if defined(RT_SOLVER_EXPLICIT)
#if defined(RT_EVOLVE_ENERGY)
            for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_E_gamma[k] = 0;}
#ifdef TRANSPORT_SUBCYCLE
            for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_E_gamma_Work[k] = 0;}
#endif
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_Flux[k] = {};}
#endif
#if defined(RT_INFRARED)
            CellP[i].Dt_Rad_E_gamma_T_weighted_IR = 0;
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[i].Dt_Rad_Flux[k] = {};}
#endif
#if defined(RT_EVOLVE_INTENSITIES)
            for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {CellP[i].Dt_Rad_Intensity[k][k_dir] = 0;}}
#endif
#endif
#ifdef MAGNETIC
            CellP[i].divB = 0; CellP[i].Face_Area = {};
#ifdef DIVBCLEANING_DEDNER
            CellP[i].DtB_PhiCorr = {};
#endif
#ifndef HYDRO_SPH
            CellP[i].DtB = {};
#ifdef DIVBCLEANING_DEDNER
            CellP[i].DtPhi = 0;
#endif
#endif
#endif // magnetic //
#ifdef COSMIC_RAY_FLUID
            CellP[i].Face_DivVel_ForAdOps = 0;
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
            CellP[i].DtCREgyNewInjectionFromShocks = 0;
#endif
            for(k=0;k<N_CR_PARTICLE_BINS;k++)
            {
                CellP[i].DtCosmicRayEnergy[k] = 0;
#if defined(CRFLUID_EVOLVE_SPECTRUM)
                CellP[i].DtCosmicRay_Number_in_Bin[k] = 0;
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
                int kAlf; for(kAlf=0;kAlf<2;kAlf++) {CellP[i].DtCosmicRayAlfvenEnergy[k][kAlf] = 0;}
#endif
            }
#endif
            P[i].wakeup = 0;
        }
    } // for (_apl over ActiveParticleList) //
}





/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the driver routine for the calculation of hydrodynamical
 *  force, fluxes, etc. */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
void hydro_force(void)
{
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    double t_preloop_start = my_second();
    hydro_force_initial_operations_preloop(); /* do initial pre-processing operations as needed before main hydro force loop */
    #ifdef GIZMO_DEBUG_RT_COOLING
    /* HYDRO_PRIM_DIAG: print core primitives BEFORE neighbor loop for target IDs */
    {static int hpd_n=0; if(hpd_n < 10) { hpd_n++;
        for(int _ii=0; _ii<NumPart; _ii++) { if(P[_ii].Type==0 && (P[_ii].ID==1 || P[_ii].ID==100 || P[_ii].ID==1000)) {
            printf("[HYDRO_PRIM] ID=%llu rho=%.10e P=%.10e cs=%.10e Mass=%.10e h=%.10e T=%.6e Gamma=%.6e DtU_pre=%.10e\n",
                (unsigned long long)P[_ii].ID, CellP[_ii].Density, CellP[_ii].Pressure, CellP[_ii].effective_soundspeed(),
                P[_ii].Mass, P[_ii].KernelRadius, CellP[_ii].Temperature, CellP[_ii].Gamma, CellP[_ii].DtInternalEnergy);
            printf("[HYDRO_PRIM] ID=%llu v=%.8e/%.8e/%.8e B=%.8e/%.8e/%.8e NV_T=%.6e/%.6e/%.6e\n",
                (unsigned long long)P[_ii].ID, CellP[_ii].VelPred[0], CellP[_ii].VelPred[1], CellP[_ii].VelPred[2],
                CellP[_ii].BPred[0], CellP[_ii].BPred[1], CellP[_ii].BPred[2],
                CellP[_ii].NV_T[0][0], CellP[_ii].NV_T[1][1], CellP[_ii].NV_T[2][2]);
            printf("[HYDRO_PRIM] ID=%llu grad_rho=%.8e/%.8e/%.8e grad_P=%.8e/%.8e/%.8e u=%.10e\n",
                (unsigned long long)P[_ii].ID, CellP[_ii].Gradients.Density[0], CellP[_ii].Gradients.Density[1], CellP[_ii].Gradients.Density[2],
                CellP[_ii].Gradients.Pressure[0], CellP[_ii].Gradients.Pressure[1], CellP[_ii].Gradients.Pressure[2],
                CellP[_ii].InternalEnergyPred);
            printf("[HYDRO_PRIM] ID=%llu Kappa_Cond=%.10e Eta_Shear=%.10e grad_u=%.8e/%.8e/%.8e grad_cs=%.8e/%.8e/%.8e dt=%.10e\n",
                (unsigned long long)P[_ii].ID, CellP[_ii].Kappa_Conduction, CellP[_ii].Eta_ShearViscosity,
                CellP[_ii].Gradients.InternalEnergy[0], CellP[_ii].Gradients.InternalEnergy[1], CellP[_ii].Gradients.InternalEnergy[2],
                CellP[_ii].Gradients.SoundSpeed[0], CellP[_ii].Gradients.SoundSpeed[1], CellP[_ii].Gradients.SoundSpeed[2],
                get_particle_timestep_in_physical(_ii, P));
        }}
        /* Also print global Temperature/Pressure/SoundSpeed sums to detect if ANY particle differs */
        {double sum_T=0, sum_P=0, sum_cs=0, sum_gamma=0; int ng=0;
        for(int _jj=0; _jj<NumPart; _jj++) { if(P[_jj].Type==0 && P[_jj].Mass>0) {
            sum_T += CellP[_jj].Temperature; sum_P += CellP[_jj].Pressure;
            sum_cs += CellP[_jj].effective_soundspeed(); sum_gamma += CellP[_jj].Gamma; ng++;}}
        printf("[HYDRO_GLOB] ngas=%d sum_T=%.10e sum_P=%.10e sum_cs=%.10e sum_Gamma=%.10e\n", ng, sum_T, sum_P, sum_cs, sum_gamma);}
        fflush(stdout);
    }}
    #endif /* GIZMO_DEBUG_RT_COOLING */
    double t_preloop = timediff(t_preloop_start, my_second());
    double t_malloc = 0, t_xchange_all = 0, t_demalloc = 0;
    double timeall = 0, timecomp = 0, timewait = 0, timecomm = 0;
    /* Neighbor-list path: GPU/Kokkos dispatch over symmetric CSR list */
    {
        ghost_writeback_zero_hydro();
        struct hydro_data_out *hydro_out = (struct hydro_data_out *) mymalloc("hydro_out",
            (gizmo_sym_num_active > 0 ? gizmo_sym_num_active : 1) * sizeof(struct hydro_data_out));
        hydro_evaluate_gpu(P, CellP, NumPart,
                           gizmo_sym_active_indices, gizmo_sym_num_active,
                           gizmo_sym_neighbor_list.offsets,
                           gizmo_sym_neighbor_list.neighbors,
                           gizmo_sym_neighbor_list.total_pairs,
                           (void *)hydro_out);
        for(int aa = 0; aa < gizmo_sym_num_active; aa++)
        {
            int ii = gizmo_sym_active_indices[aa];
            out2particle_hydra(&hydro_out[aa], ii, 0, 0);
        }
        myfree(hydro_out);
        ghost_writeback_hydro();
    }

#ifdef GIZMO_DEBUG_RT_COOLING
    /* HYDRO_RAWDTU_DIAG: print raw DtU from neighbor loop BEFORE post-loop corrections */
    {static int hraw_n=0; if(hraw_n < 10) { hraw_n++;
        for(int _ii=0; _ii<NumPart; _ii++) { if(P[_ii].Type==0 && (P[_ii].ID==1 || P[_ii].ID==100 || P[_ii].ID==1000)) {
            printf("[HYDRO_RAWDTU] ID=%llu DtU_raw=%.10e Acc=%.8e/%.8e/%.8e MaxSig=%.8e\n",
                (unsigned long long)P[_ii].ID, CellP[_ii].DtInternalEnergy,
                CellP[_ii].HydroAccel[0], CellP[_ii].HydroAccel[1], CellP[_ii].HydroAccel[2], CellP[_ii].MaxSignalVel);
        }}
        fflush(stdout);
    }}
#endif /* GIZMO_DEBUG_RT_COOLING */
    double t_postloop_start = my_second();
    hydro_final_operations_and_cleanup(); /* do final operations on results */
    double t_postloop = timediff(t_postloop_start, my_second());
    /* collect timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    /* Neighbor-list path: free symlist + remove ghosts (skipped under TRANSPORT_SUBCYCLE;
       those are cleaned up after the subcycle loop in run.cc). No-op on tree-walk build. */
    gizmo_hydro_cleanup_symlist_and_ghosts();
}
