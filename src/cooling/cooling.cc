/* Standard and Kokkos headers must precede global_data_all_struct.h
 * (its macros may conflict with stdlib names). */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"

#ifdef HYDRO_MULTIFLUID_DM_COOLING
/* DMCoolTables OWNER (mirrors the CoolTables pattern below at line ~94 for
   cooling_tables_t). __managed__ definition must precede ANY include path that
   pulls sidm/dm_fluid_functions.h (most notably eos/eos_functions.h:382), so
   nvcc rejects `extern + __managed__ in same TU`. Builder bodies (dm_*_impl
   + InitCool_dm) appear later in this file, after sidm/dm_fluid_functions.h
   is in scope (they call dm_g_integral et al.). Per the cooling_tables.h:50-57
   rule, this owner TU does NOT extern-declare DMCoolTables anywhere; the
   definition itself is sufficient and any other declaration would conflict. */
#include "../sidm/dm_cooling_tables.h"
#if defined(GIZMO_GPU_COMPILER)
__managed__ struct dm_cooling_tables_t DMCoolTables = {-1.0, 9.0, 0, nullptr, nullptr, nullptr, nullptr, nullptr};
#else
struct dm_cooling_tables_t DMCoolTables = {-1.0, 9.0, 0, nullptr, nullptr, nullptr, nullptr, nullptr};
#endif
#endif /* HYDRO_MULTIFLUID_DM_COOLING */

/* On CUDA, nvcc uses the FIRST declaration's execution-space attributes.
   proto.h declares RT functions as host-only (no __device__).  If proto.h
   is seen first, device calls get stubs returning 0.  Fix: include the RT
   function DEFINITIONS (rt_functions.h with KOKKOS_INLINE_FUNCTION bodies)
   BEFORE proto.h.  rt_functions.h and its dependencies (eos_functions.h)
   need a few forward declarations that normally come from proto.h/cooling.h. */
#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"
#include "../system/gpu_particles_arena.h"
GIZMO_GPU_FUNCTION double sigmoid_sqrt(double x); /* forward decl; defined inline in proto.h */
GIZMO_GPU_FUNCTION double ThermalProperties(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, struct particle_data *pp, struct gas_cell_data *cell);
/* Forward decls for symbols called by eos_functions.h / rt_functions.h under
   various FIRE / COOL_MOLECFRAC / SINGLE_STAR flag combinations.  Must be declared
   before the including headers so overload resolution inside their inline bodies
   succeeds.  Bodies live in proto.h / predict_functions.h / this TU's body. */
/* get_equilibrium_dust_temperature_estimate body moved to radiation/rt_functions.h
 * as KOKKOS_INLINE_FUNCTION so eos_functions.h (and any other __host__ __device__
 * caller) can call it from device passes without nvc++ warning #20011-D. */
GIZMO_GPU_FUNCTION double evaluate_NH_from_GradRho(MyFloat gradrho[3], double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp);
GIZMO_GPU_FUNCTION double evaluate_NH_from_GradRho(const Vec3<MyFloat>& gradrho, double rkern, double rho, double numngb_ndim, double include_h, int target, struct particle_data *pp);
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
/* mesh/kernel.h needed BEFORE eos_functions.h: KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER
 * macro is used inside eos_functions.h's set_eos_pressure_impl Truelove
 * pressure branch. Cannot be included from inside eos_functions.h itself —
 * kernel.h has no header guards and would trigger same-TU
 * redefinition/reinclude errors if pulled in by both eos_functions.h and
 * any other consumer in the same TU. */
#include "../mesh/kernel.h"
#endif
/* Outside the conditional above: the shutdown release is needed in every build,
   not only the ones that pull kernel.h. */
#include "../mesh/gpu_neighbor_list.h"   /* gx_touched_set_release at shutdown */
#include "../gravity/gpu_gravity_tree.h" /* gpu_node_dirty_release at shutdown */
#include "../eos/eos_functions.h"
#include "../eos/hydrogen_molecule_functions.h"
#include "../core/timestep_functions.h"
KOKKOS_FUNCTION double gas_dust_heating_coeff(int i, double T, double Tdust, struct particle_data *pp, struct gas_cell_data *cell);
#include "../radiation/rt_functions.h"
#include "../core/proto.h"
#include "../sidm/dm_fluid_functions.h"      /* set_dark_eos_pressure (called from set_eos_pressure_impl above; also pulled via eos_functions.h transitively but explicit is clearer) */
#include "../sidm/dm_cooling_functions.h"    /* HYDRO_MULTIFLUID_DM_COOLING cooling chain + do_dark_cooling_for_particle (DMCoolTables defined above) */
#include "./cooling.h"
/* evaluate_NH_from_GradRho — single source in predict_functions.h */
#include "../core/predict_functions.h"
/* CR utility functions (cosmic_ray_utilities.cc is not GPU_OBJS) */
#include "../eos/cosmic_ray_fluid/cosmic_ray_functions.h"
#if defined(GALSF_ISMDUSTCHEM_MODEL)
/* Lambda_Dust_HighTemperature_Gas_ISM (called from CoolingRate which is
 * KOKKOS_INLINE_FUNCTION) — migrated to header to clear #20011-D.
 * Phase D 2026-05-21. */
#include "../solids/ism_dust_chemistry_functions.h"
#endif
/* Simple steady-state chemistry — single source in simple_chemistry.h */
#include "./simple_chemistry.h"
/* NOTE: set_eos_pressure and finish_cooling_host_deferred_dust_updates are
 * intentionally NOT called inside do_the_cooling_for_particle (host or device).
 * Their bodies call ThermalProperties (which calls convert_u_to_temp ->
 * hydrogen_molecule chain), doubling the device stack depth and causing CUDA OOM
 * on the H200. Both are instead routed into a separate post_cooling_tail
 * device kernel:
 *   - set_eos_pressure_impl runs in the device tail, in every configuration.
 *   - dust + molecfrac + DelayTimeHII run in the device tail when
 *     GALSF_ISMDUSTCHEM_MODEL. Their host wrapper
 *     finish_cooling_host_deferred_dust_updates is no longer invoked from the
 *     GPU-path scatter loop; only the CPU OMP fallback below still calls it.
 * For non-ISMDUSTCHEM builds, the second-half molecfrac update and the
 * DelayTimeHII countdown live inside do_the_cooling_for_particle (lines ~860+,
 * #if !defined(GALSF_ISMDUSTCHEM_MODEL)) so they always run. The
 * GALSF_ISMDUSTCHEM_MODEL gate on the device-tail dust+molecfrac+HII block is
 * coupling parity with that existing branch, not a physics gate. */

/*!
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM,
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 */

/*! A file like this was originally part of the GADGET3 code developed by Volker Springel. But the code has been almost completely rewritten
 *   and massively modified/extended for GIZMO, primarily by Phil Hopkins, Mike Grudic, and Alex Richings.
 */

/* Cooling tables consolidated into a single struct (cooling_tables.h).
   With GPU offload: __managed__ so device kernels can access directly.
   Without: plain static. Pointer members point to separately allocated memory.
   Included unconditionally: builds without cooling still define the empty view
   this file hands out, and that definition needs the complete type. */
#include "cooling_tables.h"

#ifdef COOLING

#if !defined(CHIMES)
#if defined(GIZMO_GPU_COMPILER)
__managed__ struct cooling_tables_t CoolTables = {-1.0, 9.0, 0, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, 0,0,0,0,0,0,0};
#else
struct cooling_tables_t CoolTables = {-1.0, 9.0, 0, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, 0,0,0,0,0,0,0};
#endif

/* CoolTables fields are accessed directly as CoolTables.Tmin, CoolTables.gJH0, etc.
   No convenience aliases — they caused linker symbol ownership bugs on CUDA. */
#endif /* !CHIMES */

/* Hand out the tables owned by this file, as data. Everything outside this
   file reaches them this way; nothing outside declares the instances. Call it
   on the host immediately before a dispatch, so the view carries whatever the
   tables hold at that moment. */
struct PhysicsTablesView gizmo_physics_tables_view(void)
{
    struct PhysicsTablesView view = {};
#if !defined(CHIMES)
    view.cooling = &CoolTables;
#endif
#ifdef HYDRO_MULTIFLUID_DM_COOLING
    view.dm_cooling = &DMCoolTables;
#endif
    return view;
}

#ifdef HYDRO_MULTIFLUID_DM_COOLING
/* DMCoolTables OWNER's builder bodies. The __managed__ instance is defined
   earlier in this file (above the eos_functions.h include); the integral
   chain helpers used here come from sidm/dm_fluid_functions.h via
   eos_functions.h's transitive include. */

/* dm cooling table builder. Mirrors InitCoolMemory + the cooling-table
   construction loop below for the SM cooling case. Single owner TU (cooling.cc)
   for the kokkos_malloc<SharedSpace>; called once at startup from begrun.cc. */
static int dm_InitCoolMemory_impl(void)
{
    const size_t N = (size_t)(NCOOLTAB_DM + 1);
    int all_tables_obtained = 1;
    #define DMCOOLMEM(name) do { \
        DMCoolTables.name = (double *) gizmo_gpu_alloc_shared(N * sizeof(double), #name); \
        if(!DMCoolTables.name) { \
            all_tables_obtained = 0; \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), \
                     "dm cooling tables: could not allocate %s (%.1f MB)", \
                     #name, (double)(N * sizeof(double)) / (1024.0 * 1024.0)); \
            gizmo_request_controlled_stop(7717, _msg, __FILE__, __LINE__, __FUNCTION__); \
        } \
    } while(0)
    DMCOOLMEM(BetaH0);
    DMCOOLMEM(AlphaHp);
    DMCOOLMEM(AlphaHpRate);
    DMCOOLMEM(GammaeH0);
    DMCOOLMEM(Betaff);
    #undef DMCOOLMEM
    return all_tables_obtained;
}

static void dm_MakeCoolingTable_impl(void)
{
    DMCoolTables.Tmin = (All.MinGasTemp > 0.0) ? log10(All.MinGasTemp) : -1.0;
    DMCoolTables.deltaT = (DMCoolTables.Tmax - DMCoolTables.Tmin) / NCOOLTAB_DM;
    /* unit-conversion bookkeeping for recombination rate normalization */
    const double gtoGeV    = 1.0 / 1.79e-24;
    const double TtoGeV    = 1.0 / 1.16e13;
    const double InGeVcm   = 1.9733e-14;
    const double GeVtoergs = 0.001602;
    for(int i = 0; i <= NCOOLTAB_DM; i++) {
        double T = pow(10.0, DMCoolTables.Tmin + DMCoolTables.deltaT * i);
        double y2 = All.ADM_ElectronMass * pow(All.ADM_FineStructure, 2.0) * pow(C_LIGHT_CGS, 2.0)
                  / (2.0 * BOLTZMANN_CGS * T);
        DMCoolTables.BetaH0[i]      = 7.4e-18 * pow(All.ADM_FineStructure/0.01, 2.0) * sqrt(ELECTRON_MASS_SM_CGS/All.ADM_ElectronMass) * sqrt(1.0e5/T) * dm_g_integral(y2);
        DMCoolTables.Betaff[i]      = 3.7e-27 * (pow(All.ADM_FineStructure/0.01, 3.0) / pow(All.ADM_ElectronMass/ELECTRON_MASS_SM_CGS, 1.5)) * sqrt(T);
        DMCoolTables.AlphaHp[i]     = (pow(All.ADM_FineStructure, 5.0) / sqrt(pow(T,3.0) * All.ADM_ElectronMass))
                                      * sqrt(pow(2.0, 11.0) * M_PI / 27.0)
                                      * (C_LIGHT_CGS * pow(InGeVcm, 2.0) / sqrt(gtoGeV * pow(TtoGeV, 3.0)))
                                      * dm_recomb_rate_integral(y2);
        DMCoolTables.AlphaHpRate[i] = (pow(All.ADM_FineStructure, 5.0) / sqrt(T * All.ADM_ElectronMass))
                                      * sqrt(pow(2.0, 11.0) * M_PI / 27.0)
                                      * (C_LIGHT_CGS * pow(InGeVcm, 2.0) * GeVtoergs / sqrt(gtoGeV * TtoGeV))
                                      * dm_recomb_cool_integral(y2);
        DMCoolTables.GammaeH0[i]    = 2.2e-7 * pow(ELECTRON_MASS_SM_CGS/All.ADM_ElectronMass, 1.5) * sqrt(1.0e5/T) * dm_f_integral(y2);
    }
}

void InitCool_dm(void) { if(dm_InitCoolMemory_impl()) { dm_MakeCoolingTable_impl(); } }
#endif /* HYDRO_MULTIFLUID_DM_COOLING */

#if defined(CHIMES)
int ChimesEqmMode, ChimesUVBMode, ChimesInitIonState, N_chimes_full_output_freq, Chimes_incl_full_output = 1;
double chimes_rad_field_norm_factor, shielding_length_factor, cr_rate;
char ChimesDataPath[DEFAULT_PATH_BUFFERSIZE_TOUSE], ChimesEqAbundanceTable[DEFAULT_PATH_BUFFERSIZE_TOUSE], ChimesPhotoIonTable[DEFAULT_PATH_BUFFERSIZE_TOUSE];
struct gasVariables *ChimesGasVars;
struct globalVariables ChimesGlobalVars;
#ifdef CHIMES_METAL_DEPLETION
struct Chimes_depletion_data_structure *ChimesDepletionData;
#endif
#endif



/* All is redirected to a shared UVM allocation via gpu_all_mirror.h (included above). */

/* forward declarations for functions used before their definitions.
   KOKKOS_FUNCTION marks them __device__ __host__ so nvcc recognizes them as device-callable
   when they appear in call expressions within other KOKKOS_FUNCTION bodies in this TU.
   nvcc may issue advisory warning #20040 (consistent redeclaration) but this is harmless. */
KOKKOS_FUNCTION double DoCooling(double u_old, double rho, double dt, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION double CoolingRateFromU(double u, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION double CoolingRate(double logT, double rho, double n_elec_guess, double *n_elec_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
/* Functions moved to cooling_functions.h for cross-TU GPU callability:
   find_abundances_and_rates, convert_u_to_temp, convert_temp_to_u,
   return_local_gammamultiplier, return_uvb_shieldfac, ThermalProperties.
   Include with NON-INLINE KOKKOS_INLINE_FUNCTION to produce externally-visible
   host symbols that non-GPU TUs (eos.cc, density.cc, etc.) can link against.
   On CUDA, __host__ __device__ (without inline) still allows the GPU kernel
   in this TU to call these functions. The compiler can still inline them
   within this TU even without the inline keyword. */
#undef KOKKOS_INLINE_FUNCTION
#ifdef GIZMO_GPU_COMPILER
#define KOKKOS_INLINE_FUNCTION __host__ __device__
#else
#define KOKKOS_INLINE_FUNCTION
#endif
#include "cooling_functions.h"

#if !defined(CHIMES)
/* Host-facing entry points, and the ones this file's own kernels use.
   The bodies live in cooling_functions.h and take the tables as data; this
   file owns the tables, so it is the one place allowed to name them. Every
   caller elsewhere -- on the host through these symbols, on the device by
   passing its own view -- goes through that single set of bodies. */
KOKKOS_INLINE_FUNCTION
struct PhysicsTablesView cooling_owner_tables_view(void)
{
    struct PhysicsTablesView view = {};
    view.cooling = &CoolTables;
#ifdef HYDRO_MULTIFLUID_DM_COOLING
    view.dm_cooling = &DMCoolTables;
#endif
    return view;
}

KOKKOS_INLINE_FUNCTION
double return_local_gammamultiplier(int target, struct gas_cell_data *cell)
{
    struct PhysicsTablesView tables_view = cooling_owner_tables_view();
    return return_local_gammamultiplier_impl(target, cell, &tables_view);
}

KOKKOS_INLINE_FUNCTION
double find_abundances_and_rates(double logT, double rho, int target, double shieldfac, int return_cooling_mode, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess, double *LambdaExc_return, double *LambdaIon_return, double *LambdaRec_return, double *LambdaFF_return, struct particle_data *pp, struct gas_cell_data *cell)
{
    struct PhysicsTablesView tables_view = cooling_owner_tables_view();
    return find_abundances_and_rates_impl(logT, rho, target, shieldfac, return_cooling_mode, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess, LambdaExc_return, LambdaIon_return, LambdaRec_return, LambdaFF_return, pp, cell, &tables_view);
}

#if defined(EOS_SUBSTELLAR_ISM)
KOKKOS_INLINE_FUNCTION
double convert_temp_to_u(double temp, double rho, int target, double *cv, double *ne, double *nH0, double *nHp, double *nHe0, double *nHep, double *nHepp, double *mu, struct particle_data *pp, struct gas_cell_data *cell)
{
    struct PhysicsTablesView tables_view = cooling_owner_tables_view();
    return convert_temp_to_u_impl(temp, rho, target, cv, ne, nH0, nHp, nHe0, nHep, nHepp, mu, pp, cell, &tables_view);
}

KOKKOS_INLINE_FUNCTION
double convert_u_to_temp(double u, double rho, int target, double *ne, double *nH0, double *nHp, double *nHe0, double *nHep, double *nHepp, double *mu, struct particle_data *pp, struct gas_cell_data *cell)
{
    struct PhysicsTablesView tables_view = cooling_owner_tables_view();
    return convert_u_to_temp_impl(u, rho, target, ne, nH0, nHp, nHe0, nHep, nHepp, mu, pp, cell, &tables_view);
}
#else
KOKKOS_INLINE_FUNCTION
double convert_u_to_temp(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess, struct particle_data *pp, struct gas_cell_data *cell)
{
    struct PhysicsTablesView tables_view = cooling_owner_tables_view();
    return convert_u_to_temp_impl(u, rho, target, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess, pp, cell, &tables_view);
}
#endif /* EOS_SUBSTELLAR_ISM */

#endif /* !CHIMES */
#ifdef TWO_TEMPERATURE_PLASMA
#include "two_temperature_functions.h"
#endif
KOKKOS_FUNCTION double return_electron_fraction_from_heavy_ions(int target, double temperature, double density_cgs, double n_elec_HHe, struct particle_data *pp, struct gas_cell_data *cell);
/* get_equilibrium_dust_temperature_estimate: body in radiation/rt_functions.h */
KOKKOS_FUNCTION double gas_dust_heating_coeff(int i, double T, double Tdust, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION MyFloat get_FUV_G0(int target, MyFloat shieldfac, int mode, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION double evaluate_Compton_heating_cooling_rate(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg, struct gas_cell_data *cell);
KOKKOS_FUNCTION double get_background_radiation_temperature_for_emission_corrections(int target, struct gas_cell_data *cell);
KOKKOS_FUNCTION void update_explicit_molecular_fraction(int i, double dtime_cgs, struct particle_data *pp, struct gas_cell_data *cell);
KOKKOS_FUNCTION double molecfrac_rootfind_function(double fH2, double x00, double x01, double x_b_0, double x_c, double y_a, double G_LW_dt_unshielded);
KOKKOS_FUNCTION double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z, int target);
KOKKOS_FUNCTION double GetLambdaSpecies(long kspecies, long NCOOLTAB_LOCAL, long ki, long kip, long kj, double dki, double dkj, double tmin, double tdiff);
double chimes_convert_u_to_temp(double u, double rho, int target, struct particle_data *pp, struct gas_cell_data *cell);

#if defined(GALSF_ISMDUSTCHEM_MODEL)
static inline void finish_cooling_host_deferred_dust_updates(int i, double dtime, struct particle_data *pp, struct gas_cell_data *cell)
{
    update_dust_processes(i, dtime*UNIT_TIME_IN_MYR*0.001, pp, cell);
#ifdef COOL_MOLECFRAC_NONEQM
    update_explicit_molecular_fraction(i, 0.5*dtime*UNIT_TIME_IN_CGS, pp, cell); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif

#if defined(GALSF_FB_FIRE_RT_HIIHEATING) /* count off time which has passed since ionization 'clock' */
    if(cell[i].DelayTimeHII > 0) {
        cell[i].DelayTimeHII -= dtime;
#if ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION))
        if(cell[i].DelayTimeHII <= 0) {cell[i].DelayTimeHII = -DMAX(fabs(cell[i].DelayTimeHII),fabs(dtime));} // in new versions, allow this to run over to a negative number as a flag to reset the value at the beginning of the -next- timestep
#endif
    }
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
    if(cell[i].DelayTimeHII < 0) {cell[i].DelayTimeHII = 0;} // older versions simply dont allow negative values here
#endif
#endif
}
#endif

/* this is the 'parent' loop to do the cell cooling+chemistry. this is now openmp-parallelized, since the semi-implicit iteration can be a non-negligible cost.
   Uses a gather-dispatch-scatter pattern with compact arrays so that the cooling chain operates on contiguous memory
   indexed 0..N_active-1, enabling future GPU offloading. */
void cooling_parent_routine(void)
{
    PRINT_STATUS("Cooling and Chemistry update");
    /* Step 1: Determine indices of active gas particles eligible for cooling. */
    std::vector<int> cool_indices;
    cool_indices.reserve(ActiveParticleList.size());
    for (int i : ActiveParticleList)
    {
        if(P[i].Type != 0 || P[i].Mass <= 0) {continue;}
#ifdef GALSF_EFFECTIVE_EQS
        if((CellP[i].Density*All.cf_a3inv > All.PhysDensThresh) && ((All.ComovingIntegrationOn==0) || (CellP[i].Density>=All.OverDensThresh))) {continue;} /* no cooling for effective-eos star-forming particles */
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        if(CellP[i].DelayTimeCoolingSNe > 0) {continue;} /* no cooling for particles marked in delayed cooling */
#endif
        cool_indices.push_back(i);
    }
    int N_active = (int) cool_indices.size();
    if(N_active == 0) {return;}

    /* Steps 2-4: Gather / Dispatch / Scatter in batches.
     *
     * Why batching:  CUDA allocates stack_size × N_launched_threads bytes of HBM3 at
     * kernel launch.  stack_size is set to 4096 at init (sufficient for the cooling
     * chain).  Without batching a large N_active (e.g. entire startup timestep) can
     * launch millions of threads, exhausting device memory with stack alone.  Batching
     * to GPU_COOL_BATCH_SIZE caps peak stack at 4096 × GPU_COOL_BATCH_SIZE bytes per
     * launch.  At 32768 that is 128 MB — safe even when many MPI ranks share a GPU.
     *
     * The compact arrays are also kept small (GPU_COOL_BATCH_SIZE × struct_size) rather
     * than N_active × struct_size, further reducing UVM pressure.
     */
#if !defined(CHIMES)
  if(N_active >= GPU_MIN_PARTICLES_FOR_OFFLOAD) {
    GIZMO_GPU_ENSURE_ALL_FRESH(); /* device cooling reads All mirrors; tiny-N CPU path reads host All directly */
    static const int GPU_COOL_BATCH_SIZE = 32768;

    int batch_cap = (N_active < GPU_COOL_BATCH_SIZE) ? N_active : GPU_COOL_BATCH_SIZE;
    /* Whole-struct staging is shared with the other batched device loops: plain host
       buffers filled in parallel, one explicit copy each way, and device memory the
       kernel actually reads. */
    struct ParticleStagingBatch batch = {};
    int staged = particle_staging_acquire(&batch, batch_cap);
    /* dtime is captured on the host at gather time because
       get_particle_timestep_in_physical is host-only. It is the one deliberate
       exception to staging through host-then-device memory: at 8 bytes an element
       against 3080 for the structs, a third buffer plus a third copy would cost more
       than the shared-memory read it saves. */
    double *compact_dtime = (double *) gizmo_gpu_alloc_shared((size_t) batch_cap * sizeof(double), NULL);
    if(!staged || !compact_dtime) {
        /* No batch buffers means this rank cools nothing this step. The count is
         * zeroed rather than returning, so the loop below runs zero batches and
         * any collective past it is still reached by every rank. */
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "cooling: could not stage a batch of %d cells (%.1f MB); "
                 "no cooling is applied on this rank",
                 batch_cap,
                 (double)((size_t) batch_cap * (sizeof(struct particle_data)
                                                + sizeof(struct gas_cell_data)
                                                + sizeof(double))) / (1024.0 * 1024.0));
        gizmo_request_controlled_stop(7717, msg, __FILE__, __LINE__, __FUNCTION__);
        N_active = 0;
    }

    for(int batch_start = 0; batch_start < N_active; batch_start += GPU_COOL_BATCH_SIZE)
    {
        int batch_n = N_active - batch_start;
        if(batch_n > GPU_COOL_BATCH_SIZE) {batch_n = GPU_COOL_BATCH_SIZE;}

        /* Gather batch */
        if(!particle_staging_gather(&batch, cool_indices.data() + batch_start, batch_n, P, CellP)) {break;}
        /* Every index here is gas by construction (the eligibility loop above skips
           Type != 0), so the staging leaves the order alone and every slot has a cell. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int j = 0; j < batch.count; j++)
        {
            double dt = get_particle_timestep_in_physical(batch.index[j], P);
#ifdef TRANSPORT_SUBCYCLE_COOLING
            dt *= All.Transport_Subcycle_dt_fraction; /* cooling fires N times in the subcycle loop, each with dt/N */
#endif
            compact_dtime[j] = dt;
        }

        /* Dispatch batch to GPU */
        {
            struct particle_data *kp = batch.dev_P;
            struct gas_cell_data *kc = batch.dev_Cell;
            gizmo_gpu_kernel_launch("cooling_loop", batch.count, KOKKOS_LAMBDA(int j) {
                do_the_cooling_for_particle(j, kp, kc);
            }, batch_start);
        }


        /* Post-cooling tail device kernel. Three blocks, each independently
         * #ifdef-gated:
         *   - set_eos_pressure_impl, unconditionally: the equation of state is
         *     device-callable for every configuration, tabulated ones included.
         *   - GALSF_ISMDUSTCHEM_MODEL -> update_dust_processes + (gated)
         *     update_explicit_molecular_fraction + (gated) DelayTimeHII
         *     countdown. The outer GALSF_ISMDUSTCHEM_MODEL gate is coupling
         *     parity with the existing branch in do_the_cooling_for_particle
         *     (lines ~860+, #if !defined(GALSF_ISMDUSTCHEM_MODEL)), where
         *     non-ISMDUSTCHEM builds already run the second-half molecfrac
         *     update and DelayTimeHII countdown inline with cooling. It is
         *     NOT a physics descope: every Config still runs molecfrac+HII
         *     exactly once per cooling call, either inline (no-dust) or here
         *     (dust). The kernel body mirrors finish_cooling_host_deferred_dust_updates
         *     byte-for-byte so the GPU path and the CPU OMP fallback (below)
         *     compute the same physics for every supported Config (multi-path
         *     principle).
         * The scatter loop below therefore calls neither the host dust-update
         * finaliser nor set_eos_pressure: the kernel above has done both.
         * Separate kernel from the cooling loop above to keep per-launch
         * device stack depth bounded -- set_eos_pressure_impl and
         * update_dust_processes each call ThermalProperties (->
         * convert_u_to_temp -> hydrogen_molecule chain), which the cooling
         * loop already exercises; nesting would double the stack at the
         * H200 OOM limit. */
        {
            struct particle_data *kp = batch.dev_P;
            struct gas_cell_data *kc = batch.dev_Cell;
            double               *kdt = compact_dtime;
            /* Host-built, captured by value -- see eos_tables_view. */
            struct EosTableView eos_tables = eos_tables_view();
            gizmo_gpu_kernel_launch("post_cooling_tail", batch.count, KOKKOS_LAMBDA(int j) {
                set_eos_pressure_impl(j, kp, kc, &eos_tables);
#if defined(GALSF_ISMDUSTCHEM_MODEL)
                if((kdt[j] > 0) && (kc[j].Mass > 0) && (kp[j].Type == 0)) {
                    update_dust_processes(j, kdt[j] * UNIT_TIME_IN_MYR * 0.001, kp, kc);
#ifdef COOL_MOLECFRAC_NONEQM
                    update_explicit_molecular_fraction(j, 0.5 * kdt[j] * UNIT_TIME_IN_CGS, kp, kc);
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
                    if(kc[j].DelayTimeHII > 0) {
                        kc[j].DelayTimeHII -= kdt[j];
                    #if ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION))
                        if(kc[j].DelayTimeHII <= 0) {kc[j].DelayTimeHII = -DMAX(fabs(kc[j].DelayTimeHII), fabs(kdt[j]));}
                    #endif
                    }
                  #if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
                    if(kc[j].DelayTimeHII < 0) {kc[j].DelayTimeHII = 0;}
                  #endif
#endif /* GALSF_FB_FIRE_RT_HIIHEATING */
                }
#endif /* GALSF_ISMDUSTCHEM_MODEL */
            }, batch_start);
        }

        /* Scatter batch back. The equation of state has already run, in the
         * device tail kernel above, for every configuration -- there is no
         * longer a build in which it has to be finished here on the host. The
         * dust, molecular-fraction and DelayTimeHII updates that the cooling
         * solve hands to that tail run there too. The CPU OMP fallback further
         * below still does its own dust tail host-side. */
        particle_staging_scatter(&batch, P, CellP);
    }
    gpu_particles_arena_invalidate(); /* host P/CellP scattered; arena stale */

    particle_staging_release(&batch);
    if(compact_dtime) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(compact_dtime);}

  } else
#endif
  { /* CPU path: OpenMP-parallel dispatch (also used as fallback for small N on GPU builds) */
    /* Cooled in place. do_the_cooling_for_particle reads and writes only element i
     * of the arrays it is handed, so the host has nothing to gain from compacting
     * first: the copy in and the copy back are pure overhead here. The device path
     * above still stages, because there the compaction is what puts the elements
     * somewhere the kernel can read. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int j = 0; j < N_active; j++)
    {
        do_the_cooling_for_particle(cool_indices[j], P, CellP);
    }
    for(int j = 0; j < N_active; j++)
    {
        int i = cool_indices[j];
        set_eos_pressure(i, P, CellP);
#if defined(GALSF_ISMDUSTCHEM_MODEL)
        double dtime = get_particle_timestep_in_physical(i, P);
#ifdef TRANSPORT_SUBCYCLE_COOLING
        dtime *= All.Transport_Subcycle_dt_fraction; /* cooling is called N times in the subcycle loop, each with dt/N */
#endif
        if((dtime>0)&&(CellP[i].Mass>0)&&(P[i].Type==0)) {finish_cooling_host_deferred_dust_updates(i, dtime, P, CellP);}
#endif
    }
  } /* end CPU/GPU path selection */


#ifdef CHIMES /* CHIMES records some extra timing information here owing to large possible imbalances */
  CPU_Step[CPU_COOLINGSFR] += measure_time(); MPI_Barrier(MPI_COMM_WORLD);
  CPU_Step[CPU_COOLSFRIMBAL] += measure_time(); PRINT_STATUS("CHIMES chemistry and cooling finished");
#endif
}




/* ---- DEVICE-COMPILABLE COOLING FUNCTIONS ----
   KOKKOS_FUNCTION expands to __device__ __host__ when compiled by nvcc (CUDA backend) and
   to nothing on CPU-only builds. It works correctly inside #ifndef/#ifdef preprocessor
   conditional blocks — unlike OpenMP 'declare target' on NVC++ 24.7. */

/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
KOKKOS_FUNCTION
void do_the_cooling_for_particle(int i, struct particle_data *pp, struct gas_cell_data *cell)
{
#ifdef HYDRO_MULTIFLUID_DM
    /* Dark-fluid early return: adiabatic by default. Under HYDRO_MULTIFLUID_DM_COOLING
     * the dispatch runs the consolidated ADM cooling chain in
     * sidm/dm_cooling_functions.h. The whole chain is KOKKOS_INLINE_FUNCTION
     * and reads from the managed DMCoolTables struct (defined at the top of
     * this file alongside CoolTables), so this is device-clean. */
    if(pp[i].FluidType == FLUID_DM) {
#ifdef HYDRO_MULTIFLUID_DM_COOLING
        struct PhysicsTablesView dm_tables_view = cooling_owner_tables_view();
        do_dark_cooling_for_particle(i, pp, cell, &dm_tables_view);
#endif
        return;
    }
#endif
    double unew, dtime = get_particle_timestep_in_physical(i, pp), ne_in, ne_out;
#ifdef TRANSPORT_SUBCYCLE_COOLING
    dtime *= All.Transport_Subcycle_dt_fraction; /* cooling is called N times in the subcycle loop, each with dt/N */
#endif

    if((dtime>0)&&(cell[i].Mass>0)&&(pp[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
        double uold = DMAX(All.MinEgySpec, cell[i].InternalEnergy); int k; k=0; ne_in=0; ne_out=0;
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
        double uion=HIIRegion_Temp/(MEAN_MOLECULAR_WEIGHT_IONIZED*(GAMMA_DEFAULT-1.)*U_TO_TEMP_UNITS); if(cell[i].DelayTimeHII>0) {if(uold<uion) {uold=uion;}} /* u_old should be >= ionized temp if used here [unless using newer model] */
#else
        if(cell[i].DelayTimeHII < 0) { // this cell re-combined at the end of the previous timestep and has not been re-ionized yet, so we need to recombine it correctly given our sub-grid model (at fixed T not fixed U)
            cell[i].DelayTimeHII = 0; cell[i].InternalEnergy *= MEAN_MOLECULAR_WEIGHT_IONIZED/MEAN_MOLECULAR_WEIGHT_ATOMIC; cell[i].MeanMolecularWeight = MEAN_MOLECULAR_WEIGHT_ATOMIC; cell[i].Ne = DMIN(cell[i].Ne , 0.01); // assume efficient recombination here, at fixed temperature, and reset conserved quantities
            cell[i].InternalEnergyPred = cell[i].InternalEnergy;
            cell[i].HI = DMAX(0, DMIN(1, 1. - cell[i].Ne / 1.2)); // keep the neutral fraction consistent with the recombined electron fraction. the temperature is unchanged here by construction
#ifdef EOS_ANCHOR_INTERNALENERGY_IN_DRIFTS
            cell[i].u_anchor = cell[i].InternalEnergy; /* the energy and the molecular weight above were rescaled together precisely to hold the temperature fixed, so the anchor has to move with them, or the cached pair would report the change this branch exists to avoid */
#endif
            }
#endif
#endif

#if defined(RADTRANSFER)
        for(k=0;k<N_RT_FREQ_BINS;k++) {cell[i].Lambda_RadiativeCooling_toRHDBins[k]=0;} // zero these out before cooling subroutine
#endif

#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction(i, 0.5*dtime*UNIT_TIME_IN_CGS, pp, cell); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif

#ifndef COOLING_OPERATOR_SPLIT
        double DtInternalEnergyEffCGS = cell[i].DtInternalEnergy;
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        Vec3<MyDouble> grav_acc = All.cf_a2inv * pp[i].GravAccel;
#ifdef PMGRID
        grav_acc += All.cf_a2inv * pp[i].GravPM;
#endif
        DtInternalEnergyEffCGS -= All.cf_atime * dot(cell[i].GravWorkTerm, grav_acc);
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
        if(DtInternalEnergyEffCGS < 0) {
            double qfac = DMIN(0,DMAX(DMAX(-0.9, exp(DtInternalEnergyEffCGS*dtime/cell[i].InternalEnergy)-1.), All.MinEgySpec/cell[i].InternalEnergy-1.)); // equivalent to saying this wouldn't lower internal energy to below 10% in one timestep
            DtInternalEnergyEffCGS = DMAX(DtInternalEnergyEffCGS , qfac*cell[i].InternalEnergy/dtime );
            double u_gamma_minus_1 = (cell[i].gamma_eos_value()-1.) * cell[i].InternalEnergy, rho = cell[i].Density*All.cf_a3inv, pressure_thermalonly = u_gamma_minus_1 * rho;
            double vA = cell[i].Alfven_speed(), pressure_total = 0.5*vA*vA*rho + cell[i].Pressure*All.cf_a3inv;
            if(pressure_thermalonly < 0.05*pressure_total) {
                double DtInternalEnergyPdV = - u_gamma_minus_1 * (pp[i].Particle_DivVel*All.cf_a2inv); /* change from expansion in PdV term */
                DtInternalEnergyEffCGS = DMAX(DtInternalEnergyEffCGS , DMIN(DtInternalEnergyPdV, 0)); /* limit to PdV expansion change in limit where the thermal energy is small compared to the total */
            }
        }
        DtInternalEnergyEffCGS = DMIN(DtInternalEnergyEffCGS ,  1.e4*cell[i].InternalEnergy/dtime ); // equivalent to saying we cant massively enhance internal energy in a single timestep from the hydro work terms: should be big, since just numerical [shocks are real!]
        /* and convert to cgs before use in the cooling sub-routine */
        DtInternalEnergyEffCGS *= (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC);
        if(cell[i].CoolingIsOperatorSplitThisTimestep==0) {cell[i].DtInternalEnergy = DtInternalEnergyEffCGS;} // if unsplit, send this converted variable to cooling below
#endif

#if !defined(CHIMES)
        ne_in = cell[i].Ne; ne_out = ne_in; /* this variable is not defined if chimes is on */
#endif
        unew = DoCooling(uold, cell[i].Density * All.cf_a3inv, dtime, ne_in, &ne_out, i, pp, cell);
#if !defined(CHIMES)
        cell[i].Ne = ne_out; /* update the free electron variable */
#endif

#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING) /* for older model, set internal energy to minimum level if marked as ionized by stars */
        if(cell[i].DelayTimeHII > 0)
        {
            if(unew<uion) {unew=uion; if(cell[i].DtInternalEnergy<0) cell[i].DtInternalEnergy=0;}
#ifndef CHIMES
            cell[i].Ne = 1.0 + 2.0*yhelium(i, pp); /* fully ionized. note that this gives Ne as free electron fraction per H */
#ifndef COOL_GRACKLE
            cell[i].HI = 0; cell[i].MeanMolecularWeight = MEAN_MOLECULAR_WEIGHT_IONIZED; /* fully ionized, as assumed for the ionized-energy floor above */
            cell[i].Temperature = cell[i].gas_temperature_from_u(unew);
#ifdef EOS_ANCHOR_INTERNALENERGY_IN_DRIFTS
            cell[i].u_anchor = unew;
#endif
#endif
#endif
        }
#endif


#if defined(SINK_THERMALFEEDBACK)
        if(cell[i].Injected_Sink_Energy) {unew += cell[i].Injected_Sink_Energy / cell[i].Mass; cell[i].Injected_Sink_Energy = 0;}
#endif


#if defined(COSMIC_RAY_FLUID) && !defined(CRFLUID_ALT_DISABLE_LOSSES)
        CR_cooling_and_losses(i, cell[i].Ne, cell[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, dtime*UNIT_TIME_IN_CGS , pp, cell);
#endif
        

#if defined(RADTRANSFER) /* account for cooling radiation which should, according to our modules, come out in certain bands */
        double nHcgs = cell[i].nHcgs(); /* hydrogen number dens in cgs units */
        double ratefact = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * nHcgs * nHcgs / (cell[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS) * (dtime*UNIT_TIME_IN_CGS) / (UNIT_SPECEGY_IN_CGS) * cell[i].Mass; /* need to account for RSOL factors in emission/absorption rates */
        double de_u = (unew - cell[i].InternalEnergy) * cell[i].Mass; /* change in the total internal energy of the gas cell [integrating over everything] */
        double de_rad_tot_final = 0, de_rad_tot = 0; for(k=0;k<N_RT_FREQ_BINS;k++) {de_rad_tot += cell[i].Lambda_RadiativeCooling_toRHDBins[k] * ratefact;} /* energy gained by gas needs to be subtracted from radiation. positive lambda means gas cooling (gas energy loss, so radiation energy gain, so positive here) */
        double de_u_rad = -de_rad_tot, de_u_work = de_u - de_u_rad; /* variables for below showing total change in gas energy from radiation, and placeholder for the hydro work term */
        de_u_work = (cell[i].DtInternalEnergy*(UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS)*(PROTONMASS_CGS/HYDROGEN_MASSFRAC)) / nHcgs * ratefact; /* account for hydro work going into the system as an energy source */
#ifndef COOLING_OPERATOR_SPLIT
        de_u_work = DtInternalEnergyEffCGS / nHcgs * ratefact; /* use the combined and rate-limited value which is more accurately computed above */
#endif
        double de_u_radabs=0; /* need to collect absorbed photon energy to know how much energy to limit the 'dumped' energy to */
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            int k_donor = rt_get_donation_target_bin(k); /* this is used to indicate whether in the rad drift-kick loop, absorption is immediately re-radiated or not, ie. whether or not we should account for it here */
            double tau = fabs(rt_absorption_rate(i,k, pp, cell) * dtime), f_abs = 1.-exp(-tau); if(tau<0.01) {f_abs=tau*(1.-tau/2.);} /* fraction of energy absorbed in the timestep */
            double absorpted_rad_energy = DMIN(cell[i].Rad_E_gamma[k],cell[i].Rad_E_gamma_Pred[k]) * f_abs; /* estimate energy from the band that is absorbed in this timestep */
#ifdef RT_INFRARED
            if(k==RT_FREQ_BIN_INFRARED) {
                k_donor = -1; /* we use this below to indicate radiation which hasn't been re-radiated, which is handled in a special way for the adaptive bin here [which by default re-emits to itself], so set this here */
                double opacity_fraction_from_gas_absorption = rt_kappa_adaptive_IR_band(i,cell[i].Dust_Temperature,cell[i].Radiation_Temperature,-1,-1, pp, cell) / (rt_kappa_adaptive_IR_band(i,cell[i].Dust_Temperature,cell[i].Radiation_Temperature,0,0, pp, cell) + MIN_REAL_NUMBER); /* want the opacity from gas absorption as a fraction of total, because this is -not- assumed to re-radiate immediately in the drift/kick routine */
                absorpted_rad_energy *= opacity_fraction_from_gas_absorption;
            }
#endif
            if(k_donor >= 0) {continue;} /* re-emitted immediately, ignore */
            de_u_radabs += fabs(absorpted_rad_energy); /* sum up absorbed photon energy */
        }
        de_u_work += de_u_radabs; /* add this to the energy reservoir represented by the work function */
        double de_u_touse = de_u - de_u_work; /* this is the actual difference between the implicit hydro work+absorption term and the total term, i.e. a corrected de_u_rad, which we use below */

        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            if((fabs(cell[i].Lambda_RadiativeCooling_toRHDBins[k]) > MIN_REAL_NUMBER) && (fabs(de_rad_tot) > MIN_REAL_NUMBER))
            {
                double de_rad = cell[i].Lambda_RadiativeCooling_toRHDBins[k] * ratefact; /* energy gained by gas needs to be subtracted from radiation. positive lambda means gas cooling (gas energy loss, so radiation energy gain, so positive here) */
                if(fabs(de_rad) > MIN_REAL_NUMBER)
                {
                    double de_rad_min = DMIN(DMAX(-0.99*cell[i].Rad_E_gamma[k], -de_u_touse), 0); // don't let the radiation loss take all the radiation energy into negative, or more than the energy gained from cooling+heating
                    double de_rad_max = DMAX(DMIN(10.*unew*cell[i].Mass, -de_u_touse), 0); // don't let the radiation gain take more than some large factor times the current energy, or more than the energy lost from cooling+heating
                    de_rad = DMAX(DMIN(de_rad, de_rad_max), de_rad_min); // limit de_rad appropriately
                    if(fabs(de_rad) > MIN_REAL_NUMBER)
                    {
                        de_rad_tot_final += de_rad; // add to our running total                        
#ifdef RT_INFRARED  
                        if(k==RT_FREQ_BIN_INFRARED) {cell[i].Radiation_Temperature = cell[i].Radiation_Temperature_CoolingWeighted;} // need to also update the IR band temperature measure
#endif
                        double Rad_E_gamma_before = cell[i].Rad_E_gamma[k]; // save for immediate use below
                        cell[i].Rad_E_gamma[k] += de_rad; /* energy gained by gas is lost here (or vice versa if dust is acting as a net coolant) */
                        cell[i].Rad_E_gamma_Pred[k] = cell[i].Rad_E_gamma[k]; /* updated drifted */
#if defined(RT_EVOLVE_INTENSITIES)
                        int k_tmp; for(k_tmp=0;k_tmp<N_RT_INTENSITY_BINS;k_tmp++) {cell[i].Rad_Intensity[k][k_tmp] += de_rad/RT_INTENSITY_BINS_DOMEGA; cell[i].Rad_Intensity_Pred[k][k_tmp] += de_rad/RT_INTENSITY_BINS_DOMEGA;}
#endif
                        int kv; // add leading-order relativistic corrections here, accounting for gas motion in the addition/subtraction to the flux
#if defined(RT_EVOLVE_FLUX)
                        double corrfac = 0; if(Rad_E_gamma_before > 0 && cell[i].Rad_E_gamma[k] > 0) {corrfac = cell[i].Rad_E_gamma[k] / (MIN_REAL_NUMBER + Rad_E_gamma_before);}
                        if(corrfac > 0) {cell[i].Rad_Flux[k] *= corrfac; cell[i].Rad_Flux_Pred[k] *= corrfac;} else {for(kv=0;kv<3;kv++) {double fluxfac = RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*cell[i].VelPred[kv]/All.cf_atime * de_rad; cell[i].Rad_Flux[k][kv] += fluxfac; cell[i].Rad_Flux_Pred[k][kv] += fluxfac;}}
#endif
                        double momfac = 1. - de_rad / (cell[i].Mass * C_LIGHT_CODE*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission [note peculiar units here, its b/c of how we fold in the existing value of v and tilde[u] in our derivation - one rsol factor in denominator needed]
                        pp[i].dp += pp[i].Vel * ((momfac - 1.) * cell[i].Mass); pp[i].Vel *= momfac; cell[i].VelPred *= momfac;
                    }
                }
            }
        }
        
#endif // done with RHD-cooling block update
        

        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated.
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        cell[i].InternalEnergy = unew;
        cell[i].InternalEnergyPred = cell[i].InternalEnergy;

#ifdef TWO_TEMPERATURE_PLASMA
        /* 2-T plasma post-cooling update. Two operations, in order:
           (1) Apply the radiative delta to u_e. Of the cooling-step Δu_total =
               unew - uold, the hydro work component went to ions (default
               f_e=0); the radiative component went to electrons. The cooling
               iteration already evaluated rates at T_e (see the lambda above)
               so this is consistent.
           (2) Spitzer analytic e-i equilibration: redistributes (u_e, u_i)
               at fixed total via exp(-alpha*dt) decay of (T_e - T_i).
           n_i (ion number density used in both energy partition and Spitzer
           rate denominator alpha = (1 + n_e/n_i) nu_eq) = n_H + n_He = n_H * (1
           + yhelium), independent of ionization state since each helium nucleus
           contributes one heavy particle to the kinetic energy budget. */
        {
            const double rho_phys     = cell[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
            const double n_e_phys     = cell[i].n_e_cell;
            const double nH_phys      = cell[i].nHcgs();
            const double n_i_phys     = nH_phys * (1.0 + yhelium(i, pp));
            const double dt_phys      = dtime * UNIT_TIME_IN_CGS;
            const double u_old_code   = uold; /* start-of-step total u in code units */
            const double u_e_old_code = cell[i].u_e_cell;
            /* hydro_du in code units: matches the value the rootfind absorbed. */
            const double nHcgs_local  = nH_phys;
            const double ratefact     = (nHcgs_local * nHcgs_local) / rho_phys;
            const double hydro_du_cgs = ratefact * cell[i].DtInternalEnergy * dt_phys / nHcgs_local;
            const double hydro_du_code = hydro_du_cgs / UNIT_SPECEGY_IN_CGS;
            const double du_total_code = cell[i].InternalEnergy - u_old_code;
            const double du_rad_code   = du_total_code - hydro_du_code;
            /* (1) absorb radiative delta into u_e (always to electrons), plus
               fraction f_e of NON-conduction hydro dissipation, plus the
               full conduction contribution (Spitzer-Härm electron heat under
               TWO_TEMPERATURE_PLASMA bit 2). f_e=0 default -> only radiative
               + conduction (if bit 2 active); f_e=1 -> electrons absorb full du. */
            const double f_e_shock     = All.TwoTemp_ShockElectronFraction;
#if (TWO_TEMPERATURE_PLASMA & 4) && defined(CONDUCTION)
            const double cond_du_cgs   = ratefact * cell[i].DtInternalEnergy_FromConduction * dt_phys / nHcgs_local;
            const double cond_du_code  = cond_du_cgs / UNIT_SPECEGY_IN_CGS;
            const double noncond_du_code = hydro_du_code - cond_du_code;
            double u_e_intermediate    = u_e_old_code + du_rad_code + f_e_shock * noncond_du_code + cond_du_code;
#else
            double u_e_intermediate    = u_e_old_code + du_rad_code + f_e_shock * hydro_du_code;
#endif
            if(!(u_e_intermediate > 0)) {u_e_intermediate = 1e-30 * (u_e_old_code > 0 ? u_e_old_code : 1.0);}
            const double T_e_before = two_temp_T_e_from_u_e(u_e_intermediate, rho_phys, n_e_phys);
            /* (2) Spitzer analytic relaxation at fixed u_total */
            const double u_e_new = two_temp_relax_step(cell[i].InternalEnergy, u_e_intermediate,
                                                       n_e_phys, n_i_phys, rho_phys, dt_phys);
            cell[i].u_e_cell = u_e_new;
            cell[i].T_e_cell = two_temp_T_e_from_u_e(u_e_new, rho_phys, n_e_phys);
        }
#endif

#ifndef COOLING_OPERATOR_SPLIT
        if(cell[i].CoolingIsOperatorSplitThisTimestep==0) {cell[i].DtInternalEnergy=0;} // if unsplit, zero the internal energy change here
        /* when TRANSPORT_SUBCYCLE_COOLING, DtInternalEnergy is saved/restored in run.cc around each cooling call */
#endif

#if !defined(GALSF_ISMDUSTCHEM_MODEL)
#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction(i, 0.5*dtime*UNIT_TIME_IN_CGS, pp, cell); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif
        
#if defined(GALSF_FB_FIRE_RT_HIIHEATING) /* count off time which has passed since ionization 'clock' */
        if(cell[i].DelayTimeHII > 0) {
            cell[i].DelayTimeHII -= dtime;
#if ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION))
            if(cell[i].DelayTimeHII <= 0) {cell[i].DelayTimeHII = -DMAX(fabs(cell[i].DelayTimeHII),fabs(dtime));} // in new versions, allow this to run over to a negative number as a flag to reset the value at the beginning of the -next- timestep
#endif
        }
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
        if(cell[i].DelayTimeHII < 0) {cell[i].DelayTimeHII = 0;} // older versions simply dont allow negative values here
#endif
#endif // GALSF_FB_FIRE_RT_HIIHEATING
#endif // !defined(GALSF_ISMDUSTCHEM_MODEL)

    } // closes if((dt>0)&&(cell[i].Mass>0)&&(pp[i].Type==0)) check
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
KOKKOS_FUNCTION
double DoCooling(double u_old, double rho, double dt, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
    double u, du; u=0; du=0;
    if(rho <= 0 || dt == 0) {return DMAX(u,All.MinEgySpec);}

#ifdef COOL_GRACKLE
#ifndef COOLING_OPERATOR_SPLIT
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below to split the step */
    if(cell[target].CoolingIsOperatorSplitThisTimestep==0)
    {
        du = dt * cell[target].DtInternalEnergy / ( (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC));
        u_old += 0.5*du;
        u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
        /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
         term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
         the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
         to which we then add the heating term from hydro and compute the solution over a full timestep */
        double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(r-1)/log(r);} r=DMAX(0,DMIN(r,1));
        du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
    }
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
    return DMAX(u,All.MinEgySpec);
#endif

#ifdef CHIMES
    chimes_update_gas_vars(target, pp, cell);

    /* Call CHIMES to evolve the chemistry and temperature over
     * the hydro timestep. */
    chimes_network(&(ChimesGasVars[target]), &ChimesGlobalVars);

    // Compute updated internal energy
    u = (double) ChimesGasVars[target].temperature * BOLTZMANN_CGS / ((cell[target].gamma_eos_value()-1) * PROTONMASS_CGS * calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars));
    u /= UNIT_SPECEGY_IN_CGS;  // code units

#ifdef CHIMES_TURB_DIFF_IONS 
    chimes_update_turbulent_abundances(target, 1, pp, cell); 
#endif 

    return DMAX(u, All.MinEgySpec);

#else // CHIMES

    int iter=0, iter_upper=0, iter_lower=0, iter_condition = 0, maxiter_uplo=10*MAXITER; double LambdaNet, ratefact, u_upper, u_lower;
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double u_min = All.MinEgySpec * UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;
    u = u_upper = u_lower =  u_old; /* initialize values */
#ifdef TWO_TEMPERATURE_PLASMA
    /* 2-T plasma: the radiative cooling rate is physically a function of
       T_e, not T_total. Inside the rootfind we track T_e implicitly: of the
       trial du, the hydro-source portion goes to ions and the remaining
       radiative portion goes to electrons (default f_e=0).
       u_e_trial = u_e_old + (du - hydro_du), and T_e_trial = (2/3) rho
       u_e_trial / (n_e k_B). Pass T_e_trial in place of T_total to CoolingRate
       so all rate components (line cooling, free-free, recombination, photo-
       heating, Compton) automatically use the electron temperature. The hydro-
       source addition inside CoolingRate (Q += DtInternalEnergy/nHcgs) is
       T-independent and is preserved unchanged. The bisection still converges
       since dT_e/du > 0 (same monotonicity as the single-T case). */
    const double u_e_old_cgs_2T = cell[target].u_e_cell * UNIT_SPECEGY_IN_CGS;
    const double n_e_phys_2T    = cell[target].n_e_cell;
    const double hydro_du_2T    = ratefact * cell[target].DtInternalEnergy * dt / nHcgs;
    /* Fraction of hydro dissipation deposited into electrons. f_e=0 (default,
       collisionless-shock limit, ions absorb all dissipation) recovers pure
       (du - hydro_du) for u_e_trial; f_e=1 makes electrons absorb the same du as
       u_total (single-fluid limit). */
    const double f_e_2T         = All.TwoTemp_ShockElectronFraction;
    /* Under bit 2 of TWO_TEMPERATURE_PLASMA, conduction is electron heat
       (Spitzer-Härm) and routes 100% to u_e regardless of f_e. The non-
       conduction hydro work continues to follow f_e. Conduction's contribution
       is captured separately in the hydro pair loop and stored in
       cell.DtInternalEnergy_FromConduction. */
#if (TWO_TEMPERATURE_PLASMA & 4) && defined(CONDUCTION)
    const double cond_du_2T       = ratefact * cell[target].DtInternalEnergy_FromConduction * dt / nHcgs;
    const double noncond_du_2T    = hydro_du_2T - cond_du_2T;
    const double hydro_du_to_e_2T = f_e_2T * noncond_du_2T + cond_du_2T;
#else
    const double hydro_du_to_e_2T = f_e_2T * hydro_du_2T;
#endif
#endif
    auto cooling_rootfind_function = [&](double du) { /* control the *relative* error on the *change* in u */
#ifdef TWO_TEMPERATURE_PLASMA
        double u_e_trial = u_e_old_cgs_2T + (du - hydro_du_2T) + hydro_du_to_e_2T;
        if(!(u_e_trial > 0)) {u_e_trial = (u_e_old_cgs_2T > 0) ? 1e-30 * u_e_old_cgs_2T : 1e-30;}
        double T_e_trial = (2./3.) * rho * u_e_trial / (n_e_phys_2T * BOLTZMANN_CGS);
        if(!(T_e_trial > 0)) {T_e_trial = 1.0;}
        return du - ratefact * CoolingRate(log10(T_e_trial), rho, ne_guess, ne_eval, target, pp, cell) * dt;
#else
        return du - ratefact * CoolingRateFromU(u_old+du, rho, ne_guess, ne_eval, target, pp, cell) * dt;
#endif
    };
    double du_net = cooling_rootfind_function(u - u_old), du_net_upper = du_net, du_net_lower = du_net;

    /* bracketing */
    double u_step_fac = 1.1;
    int bracket_iter = 0, skip_rootfind = 0;
    while(du_net_upper * du_net_lower > 0 && bracket_iter<MAXITER){
        if((u_lower <= All.MinEgySpec) && (du_net_lower > 0) && (du_net_upper > 0)){skip_rootfind = 1; break;} // will never find the root because bouncing off the lower limit
        u_upper *= u_step_fac; 
        du_net_upper = cooling_rootfind_function(u_upper - u_old);
	if(du_net*du_net_upper < 0){u_lower = u; du_net_lower = du_net; break;} // let u_upper and u_old be the brackets
	u_lower = DMAX(u_lower/u_step_fac,All.MinEgySpec); // bound u_lower because we don't trust cooling function below this
        du_net_lower = cooling_rootfind_function(u_lower - u_old);
	if(du_net*du_net_lower < 0){u_upper = u; du_net_upper = du_net; break;} // let u_lower an u_old be the brackets
        u_step_fac *= 1.1;
        bracket_iter++;
    }
    if(fabs(du_net_upper) < MIN_REAL_NUMBER || !isfinite(du_net_upper)) {u=u_upper;}
    else if(fabs(du_net_lower) < MIN_REAL_NUMBER || !isfinite(du_net_lower)) {u=u_lower;}
    else {
        if(!skip_rootfind){ // assuming we're not bouncing off the min temp
            if((du_net_upper * du_net_lower >= 0) || isnan(du_net_lower) || isnan(du_net_upper)) {PRINT_WARNING("Could not bracket cooling solution. ID=%lld u_min=%g u=%g u_lower=%g u_upper=%g f_lower=%g f_upper=%g\n", (long long)(long long)target /* particle index */, u_min, u, u_lower,u_upper, du_net_lower, du_net_upper); endrun(10);}
            
            /* core iteration to convergence */
            BrentRootfindResult rfr = brent_rootfind(cooling_rootfind_function,
                                                    u_upper-u_old, u_lower-u_old,
                                                    du_net_upper, du_net_lower,
                                                    1e-2, 1e-15 * u_old, MAXITER);
            u = rfr.x + u_old;

            /* crash condition */
            if((rfr.iter >= MAXITER) || isnan(u)) {
                printf("failed to converge in DoCooling(cell): x_new=%g x_a=%g x_b=%g error=%g u_in=%g u_upper=%g u_lower=%g rho_in=%g dt=%g ne_in=%g ne_out=%g target=%d ID=%ld \n",rfr.x, u_upper-u_old, u_lower-u_old, rfr.error, u_old, u_upper, u_lower, rho,dt,ne_guess,*ne_eval,target, (long)(long long)target); endrun(10);
            }
            u = DMAX(u_min,u);
        } else {u = All.MinEgySpec;}
    }

    double specific_energy_codeunits_toreturn = u / UNIT_SPECEGY_IN_CGS;    /* in internal units */
    /* solve the chemical state once at the converged energy, outside the iteration: this is the key
       chemistry update, and it is the saved thermochemistry that the equation of state, the radiation
       modules and the dust chemistry read back instead of re-solving it themselves. u and rho are cgs
       here. the values carried through the iteration are the initial guess. every saved field comes
       from this one evaluation, so they describe the same state as the energy actually returned:
       the electron fraction left by the last rate evaluation is at a trial energy, not this one */
    {
        double mu_final = 1, ne_final = *ne_eval, nHI = cell[target].HI, nHII = DMAX(0, 1. - nHI), nHeI = 1, nHeII = 0, nHeIII = 0;
        double temp_final = convert_u_to_temp(u, rho, target, &ne_final, &nHI, &nHII, &nHeI, &nHeII, &nHeIII, &mu_final, pp, cell);
        *ne_eval = ne_final; /* what the caller saves, and what seeds the next step */
        cell[target].Ne = ne_final; cell[target].HI = nHI;
        cell[target].MeanMolecularWeight = mu_final; /* the composition half of the cache: it stays valid as the energy moves */
        cell[target].Temperature = temp_final;
#ifdef EOS_ANCHOR_INTERNALENERGY_IN_DRIFTS
        cell[target].u_anchor = u / UNIT_SPECEGY_IN_CGS; /* exact re-anchor at the converged state */
#endif
        cell[target].Gamma = cell[target].gamma_eos_value(); /* the other half, so the pair describes this solve rather than the previous equation-of-state call */
#ifdef RT_CHEM_PHOTOION
        cell[target].HII = nHII;
#ifdef RT_CHEM_PHOTOION_HE
        cell[target].HeI = nHeI; cell[target].HeII = nHeII; cell[target].HeIII = nHeIII;
#endif
#endif
    }

    /* safe return */
    return specific_energy_codeunits_toreturn;
#endif // CHIMES
}



#ifndef CHIMES
/* returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
#if defined(COOL_GRACKLE) && !defined(GALSF_EFFECTIVE_EQS)
    double LambdaNet = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(LambdaNet >= 0) LambdaNet = 0.0;
    return LambdaNet / UNIT_TIME_IN_CGS;
#else
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    double LambdaNet = CoolingRateFromU(u_old, rho, ne_guess, ne_eval, target, pp, cell);
    if(LambdaNet >= 0) {return 0;} /* net heating due to UV background */
    return u_old / (-(nHcgs * nHcgs / rho) * LambdaNet) / UNIT_TIME_IN_CGS;
#endif
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
KOKKOS_FUNCTION
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
    if(fac <= 0) {return 0.01*m_old;} /* the hot phase is actually colder than the cold reservoir! */
    double m, dm, m_lower, m_upper, ratefact, LambdaNet;
    int iter = 0;

    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    fac /= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho * fac;
    m = m_old; m_lower = m; m_upper = m;
    LambdaNet = CoolingRateFromU(u, rho, ne_guess, ne_eval, target, pp, cell);

    /* bracketing */
    if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
        m_upper *= sqrt(1.1); m_lower /= sqrt(1.1);
        while(m_upper - m_old - m_upper * m_upper / m_old * ratefact * CoolingRateFromU(u, rho * m_upper / m_old, ne_guess, ne_eval, target, pp, cell) * dt < 0)
        {
            m_upper *= 1.1; m_lower *= 1.1;
        }
    }
    if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
        m_lower /= sqrt(1.1); m_upper *= sqrt(1.1);
        while(m_lower - m_old - m_lower * m_lower / m_old * ratefact * CoolingRateFromU(u, rho * m_lower / m_old, ne_guess, ne_eval, target, pp, cell) * dt > 0)
        {
            m_upper /= 1.1; m_lower /= 1.1;
        }
    }

    do
    {
        m = 0.5 * (m_lower + m_upper);
        LambdaNet = CoolingRateFromU(u, rho * m / m_old, ne_guess, ne_eval, target, pp, cell);
        if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0) {m_upper = m;} else {m_lower = m;}
        dm = m_upper - m_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("->m= %g\n", m);}
    }
    while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);
    if(iter >= MAXITER) {printf("failed to converge in DoInstabilityCooling(cell): m_in=%g u_in=%g rho=%g dt=%g fac=%g ne_in=%g target=%d ID=%ld\n",m_old,u,rho,dt,fac,ne_guess,target,(long)(long long)target /* particle index */); endrun(11);}
    return m;
}

#endif // !(CHIMES)






#ifdef CHIMES
/* This function converts thermal energy to temperature, using the mean molecular weight computed from the non-equilibrium CHIMES abundances. */
double chimes_convert_u_to_temp(double u, double rho, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
  return u * (cell[target].gamma_eos_value()-1) * PROTONMASS_CGS * ((double) calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars)) / BOLTZMANN_CGS;
}
// CHIMES
#elif  defined(EOS_SUBSTELLAR_ISM)
/*
Returns the internal energy of the gas mixture per unit mass in erg/g:

u = Etot/Mtot = SUM_i(N_i E_i) / SUM_i(N_i m_i) over all species i with energy E_i, mass m_i, and relative abundance N_i

Argument cv will return the specific heat capacity at constant volume du/dT
 */
/* convert_temp_to_u, convert_u_to_temp, find_abundances_and_rates:
   definitions moved to cooling_functions.h (included above). */
#else
/* default (non-CHIMES, non-EOS_SUBSTELLAR_ISM): all definitions in cooling_functions.h */
#endif // CHIMES / EOS_SUBSTELLAR_ISM

#ifndef CHIMES
/*  this function first computes the self-consistent temperature and abundance ratios, and then it calculates (heating rate-cooling rate)/n_h^2 in cgs units */
KOKKOS_FUNCTION
double CoolingRateFromU(double u, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
    double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu; nH0_guess = DMAX(0,DMIN(1,1.-ne_guess/1.2));
    double temp = convert_u_to_temp(u, rho, target, &ne_guess, &nH0_guess, &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess, &mu, pp, cell);
    double Lambda = CoolingRate(log10(temp), rho, ne_guess, ne_eval, target, pp, cell);
    return Lambda;
}


#endif // !(CHIMES)



extern FILE *fd;



#ifndef CHIMES
/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units
 */
KOKKOS_FUNCTION
double CoolingRate(double logT,  double rho, double n_elec_guess, double *n_elec_eval, int target, struct particle_data *pp, struct gas_cell_data *cell)
{
    double n_elec=n_elec_guess, nH0, nHe0, nHp, nHep, nHepp, mu; /* ionization states [computed below] */
    double Lambda, Heat, LambdaFF, LambdaCompton, LambdaExc, LambdaExcH0, LambdaExcHep, LambdaIon, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
    double LambdaRec, LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd, T, T_cmb_radeff, shieldfac, LambdaMol, LambdaMetal, LambdaPElec, LambdaDust, Heat_Ion_from_UVB, Heat_Ion_from_RHD;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    Lambda=0; Heat=0; LambdaMol=0; LambdaFF=0; LambdaRec=0; LambdaExc=0; LambdaIon=0; LambdaMetal=0; LambdaCompton=0; LambdaPElec=0; LambdaDust=0; Heat_Ion_from_UVB=0; Heat_Ion_from_RHD=0; /* make sure these are all initialized to zero */
    if(logT <= CoolTables.Tmin) {logT = CoolTables.Tmin + 0.5 * CoolTables.deltaT;}	/* floor at CoolTables.Tmin */
    if(!isfinite(rho)) {return 0;}
    T = pow(10.0, logT);
    T_cmb_radeff = get_background_radiation_temperature_for_emission_corrections(target, cell); /* CMB temperature, used below */

    /* some blocks below to define useful variables before calculation of cooling rates: */

#ifdef COOL_METAL_LINES_BY_SPECIES
    double Z[NUM_METAL_SPECIES] = {0.}; // gas-phase metallicity
    if(target>=0)
    {
        int k; 
#if defined(GALSF_ISMDUSTCHEM_MODEL)
        for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {Z[k] = DMAX(0.,pp[target].Metallicity[k]-cell[target].ISMDustChem_Dust_Metal[k]);}
#else
        for(k=0;k<NUM_METAL_SPECIES;k++) {Z[k] = pp[target].Metallicity[k];}
#endif
    } else { /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
        int k; for(k=0;k<NUM_METAL_SPECIES;k++) {Z[k]=All.SolarAbundances[k];}
    }
#endif
    double local_gammamultiplier = return_local_gammamultiplier(target, cell);
    shieldfac = return_uvb_shieldfac(target, local_gammamultiplier*CoolTables.gJH0/1.0e-12, nHcgs, logT, cell);
    
#if defined(COOL_LOW_TEMPERATURES)
    double Tdust = 30.; /* set variables needed for dust heating/cooling. if dust cooling not calculated, default to 0 */
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate(target, shieldfac, T, pp, cell);
#endif
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(SINGLE_STAR_SINK_DYNAMICS) && !defined(SINGLE_STAR_FB_RT_HEATING)
    Tdust = DMIN(DMAX(10., T_cmb_radeff),300.); // runs looking at colder clouds, use a colder default dust temp [floored at CMB temperature] //
#endif
#endif


#if defined(RT_CHEM_PHOTOION) || defined(RT_PHOTOELECTRIC)
    double cx_to_kappa = HYDROGEN_MASSFRAC / PROTONMASS_CGS * UNIT_MASS_IN_CGS; // pre-factor for converting cross sections into opacities
#endif
    if(logT < CoolTables.Tmax)
    {
        /* get ionization states for H and He with associated ionization, collision, recombination, and free-free heating/cooling */
        Lambda += find_abundances_and_rates(logT, rho, target, shieldfac, 1, &n_elec, &nH0, &nHp, &nHe0, &nHep, &nHepp, &mu, &LambdaExc, &LambdaIon, &LambdaRec, &LambdaFF, pp, cell); /* adds all of these to our running total for cooling */
        *n_elec_eval = n_elec; /* save this value for the output cycle */
        LambdaCompton = evaluate_Compton_heating_cooling_rate(target,T,nHcgs,n_elec,shieldfac, cell); /* note this can have either sign: heating or cooling */
        if(LambdaCompton > 0) {Lambda += LambdaCompton;}

#if defined(RT_CHEM_PHOTOION) && defined(METALS)
        /* Nebular (photoionized forbidden-line) cooling from O+, O++, N+, S+, Ne+: Kim, Gong, Kim &
           Ostriker 2023 (ApJS 264, 10), Eq. 47. The tabulated metal-line cooling below assumes collisional
           ionization, which under-predicts the forbidden-line cooling of PHOTOIONIZED gas; this supplies
           the coolant that sets the ~1e4 K equilibrium of an HII region. The rate is per nH^2 like the
           other channels, since the volumetric rate goes as n_e*n_H+. METALS is in the guard because the
           metallicity arrays exist only under it, and the term is metallicity-scaled in any case. */
        if(T > 2.0e3 && T < 5.0e4 && n_elec > 0 && nHp > 0 && target >= 0) {
            double T4 = T * 1.0e-4, lnT4 = log(T4), ne_100 = nHcgs * n_elec / 100.; /* electron density in units of 100 cm^-3 */
            double log10_fneb = 0.692 + lnT4*(-0.586 + lnT4*(0.816 + lnT4*(-0.505 + lnT4*(0.118 + lnT4*(0.00766 - 0.00508*lnT4)))));
            double LambdaNeb = n_elec * nHp * (pp[target].Metallicity[0] / All.SolarAbundances[0])
                * 3.68e-23 * exp(-DMIN(3.86/T4, 100.)) / sqrt(T4) * pow(10., log10_fneb)
                / (1. + 0.12*pow(ne_100, 0.38 - 0.12*lnT4)); /* collisional de-excitation: multi-line, so a temperature-dependent power rather than a single critical density */
            LambdaNeb *= 1. / (1. + exp(DMIN(10.*(T - 2.75e4)/1.5e4, 60.))); /* taper out over 2e4-3.5e4 K, handing off to the collisional metal-line tables below so this is not double-counted through the transition */
            if(LambdaNeb > 0) {Lambda += LambdaNeb;}
        }
#endif
        
#ifdef COOL_METAL_LINES_BY_SPECIES
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
        if(CoolTables.J_UV != 0)
#else
        if((CoolTables.J_UV != 0)&&(logT > 4.00))
#endif
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies(nHcgs, logT, Z, target); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= n_elec;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted; ni not used b/c it should vary species-to-species */
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
            if(logT<2) {LambdaMetal *= exp(-DMIN((2.-logT)*(2.-logT)/0.1,40.));}
            if(LambdaMetal > 0) {LambdaMetal *= ((T-T_cmb_radeff)/(T+T_cmb_radeff));}
            if(LambdaMetal > 0) {Lambda += LambdaMetal;}
#else
            Lambda += LambdaMetal;
#endif
#if defined(OUTPUT_COOLRATE_DETAIL)
            if(target >= 0) {cell[target].MetalCoolingRate = LambdaMetal;}
#endif
        }
#endif

#if defined(GALSF_ISMDUSTCHEM_HIGHTEMPDUSTCOOLING)
        LambdaDust = Lambda_Dust_HighTemperature_Gas_ISM(target,T,n_elec, pp, cell);
        Lambda += LambdaDust;
#if defined(OUTPUT_COOLRATE_DETAIL)
        if(target >= 0) {cell[target].DustCoolingRate = LambdaDust;}
#endif
#endif

#ifdef COOL_LOW_TEMPERATURES
        if(logT <= 5.3)
        {
            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
            LambdaMol = 2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));
            LambdaMol *= (1-shieldfac) / (1. + nHcgs/700.); // above the critical density, cooling rate suppressed by ~1/n; use critical density of CO[J(1-0)] as a proxy for this
            double Z_sol=1, truncation_factor=1; /* if don't have actual metallicities, we'll assume solar */
            if(logT>4.5) {double dx=(logT-4.5)/0.20; truncation_factor *= exp(-DMIN(dx*dx,40.));} /* continuous cutoff here just to avoid introducing artificial features in temperature-density */
#ifdef COOL_METAL_LINES_BY_SPECIES
            Z_sol = Z[0] / All.SolarAbundances[0]; /* use actual metallicity for this */
#endif
            LambdaMol *= (1+Z_sol)*(0.001 + 0.1*nHcgs/(1.+nHcgs) + 0.09*nHcgs/(1.+0.1*nHcgs) + Z_sol*Z_sol/(1.0+nHcgs)); // gives very crude estimate of metal-dependent terms //
#if defined(COOL_METAL_LINES_BY_SPECIES) && ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION))
            double column = evaluate_NH_from_GradRho(cell[target].Gradients.Density,pp[target].KernelRadius,cell[target].Density,pp[target].NumNgb,1,target,pp) * UNIT_SURFDEN_IN_CGS; // converts to cgs
            double Z_C = DMAX(1.e-6, Z[2]/All.SolarAbundances[2]), sqrt_T=sqrt(T), ncrit_CO=1.9e4*sqrt_T, Sigma_crit_CO=3.0e-5*T/Z_C, T3=T/1.e3, EXPmax=90.; // carbon abundance (relative to solar and 1/2 factor for original assumed 0.5 depletion), critical density and column
#if defined(GALSF_ISMDUSTCHEM_MODEL)
            Z_C = DMAX(1.e-6, Z[2]/(0.5*All.SolarAbundances[2])); // gas-phase carbon abundance (relative to solar/2, usual assumption implicitly)
#endif
            // TODO: can now get detailed C+, C, O, and CO abundances from SIMPLE_STEADYSTATE_CHEMISTRY routines, and compute the explicit cooling terms for slightly more accurate cooling here
// #ifdef SIMPLE_STEADYSTATE_CHEMISTRY            
// #endif
            double f_Cplus_CCO=1./(1.+nHcgs/3.e3); // very crude estimate used to transition between C+ cooling curve and C/CO [nearly-identical] cooling curves above C+ critical density, where C+ rate rapidly declines        
            double photoelec=get_FUV_G0(target, shieldfac, 0, pp, cell); f_Cplus_CCO = (nHcgs/(340.*DMAX(0.1,photoelec))); f_Cplus_CCO=1./(1.+f_Cplus_CCO*f_Cplus_CCO/sqrt_T); // fco/(1-fco) ~ 0.0022 * ((n/50 cm^-3)/G0)^2 * (100K/T)^(1/2) from Tielens
            double Lambda_Cplus = Z_C * (4.7e-28 * (pow(T,0.15) + 1.04e4*n_elec/sqrt_T) * exp(-DMIN(91.211/T,EXPmax)) + 2.08e-29*exp(-DMIN(23.6/T,EXPmax))); // fit from Barinovs et al., ApJ, 620, 537, 2005, and Wilson & Bell MNRAS 337 1027 2002; assuming factor of 0.5 depletion factor in ISM; rate per C+ relative to solar; + plus [CI]-609 µm line cooling from Hocuk⋆ et al. 2016MNRAS.456.2586H
            double Lambda_CCO = Z_C * T*sqrt_T * 2.73e-31 / (1. + (nHcgs/ncrit_CO)*(1.+1.*DMAX(column,0.017)/Sigma_crit_CO)); // fit from Hollenbach & McKee 1979 for CO (+CH/OH/HCN/OH/HCl/H20/etc., but those don't matter), with slight re-calibration of normalization (factor ~1.4 or so) to better fit the results from the full Glover+Clark network. As Glover+Clark show, if you shift gas out of CO into C+ and O, you have almost no effect on the integrated cooling rate, so this is a surprisingly good approximation without knowing anything about the detailed chemical/molecular state of the gas. uncertainties in e.g. ambient radiation are -much- larger. also note this rate is really carbon-dominated as the limiting abundance, so should probably use that.
            
            /* Large-velocity-gradient limiter for the CO cooling rate from Whitworth & Jaffa arXiv:1811.06814; compute Lambda_HI and use this as an upper bound */
            double gradv_norm_kms_pc = cell[target].velocity_gradient_norm() * UNIT_VEL_IN_KMS / UNIT_LENGTH_IN_PC; // velocity gradient in km/s/pc
            double Lambda_CO_HI = 4.42e-28 * gradv_norm_kms_pc * pow(T,4) / (nHcgs * nHcgs);
            double Lambda_CO = DMIN(Lambda_CO_HI, (1.-f_Cplus_CCO) * Lambda_CCO); // Let the LVG value be the limiter

            double Lambda_Metals = f_Cplus_CCO * Lambda_Cplus + Lambda_CO; // interpolate between both regimes //
            /* in the above Lambda_Metals expression, the column density expression attempts to account for the optically-thick correction in a slab. this is largely redundant (not exactly, b/c this is specific for CO-type molecules) with our optically-thick cooling module already included below, so we will not double-count it here [coefficient set to zero]. But it's included so you can easily turn it back on, if desired, instead of using the module below. */
            double Lambda_H2_thick = (6.7e-19*exp(-DMIN(5.86/T3,EXPmax)) + 1.6e-18*exp(-DMIN(11.7/T3,EXPmax)) + 3.e-24*exp(-DMIN(0.51/T3,EXPmax)) + 9.5e-22*pow(T3,3.76)*exp(-DMIN(0.0022/(T3*T3*T3),EXPmax))/(1.+0.12*pow(T3,2.1))) / nHcgs; // super-critical H2-H cooling rate [per H2 molecule]
            double Lambda_HD_thin = ((1.555e-25 + 1.272e-26*pow(T,0.77))*exp(-DMIN(128./T,EXPmax)) + (2.406e-25 + 1.232e-26*pow(T,0.92))*exp(-DMIN(255./T,EXPmax))) * exp(-DMIN(T3*T3/25.,EXPmax)); // optically-thin HD cooling rate [assuming all D locked into HD at temperatures where this is relevant], per molecule
            double f_molec = 0.5 * Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(CoolTables.gJH0/2.29e-10), pp, cell); // [0.5*f_molec for H2/HD cooling b/c cooling rates above are per molecule, not per nucleon]

            double q = logT - 3., Y_Hefrac=DMAX(0.,DMIN(1.,Z[1])), X_Hfrac=DMAX(0.,DMIN(1.,1.-Y_Hefrac-Z[0])); // variable used below
            double Lambda_H2_thin = DMAX(nH0-2.*f_molec,0) * X_Hfrac * pow(10., DMAX(-103. + 97.59*logT - 48.05*logT*logT + 10.8*logT*logT*logT - 0.9032*logT*logT*logT*logT , -50.)); // sub-critical H2 cooling rate from H2-H collisions [per H2 molecule]; this from Galli & Palla 1998
            Lambda_H2_thin += Y_Hefrac * pow(10., DMAX(-23.6892 + 2.18924*q -0.815204*q*q + 0.290363*q*q*q -0.165962*q*q*q*q + 0.191914*q*q*q*q*q, -50.)); // H2-He; often more efficient than H2-H at very low temperatures (<100 K); this and other H2-x terms below from Glover & Abel 2008
            Lambda_H2_thin += f_molec * X_Hfrac * pow(10., DMAX(-23.9621 + 2.09434*q -0.771514*q*q + 0.436934*q*q*q -0.149132*q*q*q*q -0.0336383*q*q*q*q*q, -50.)); // H2-H2; can be more efficient than H2-H when H2 fraction is order-unity
            Lambda_H2_thin += nHp * X_Hfrac * pow(10., DMAX(-21.7167 + 1.38658*q -0.379153*q*q + 0.114537*q*q*q -0.232142*q*q*q*q + 0.0585389*q*q*q*q*q, -50.)); // H2-H+; very efficient if somehow appreciable H+ fraction remains
            double logLambdaH2_e = -34.2862 -48.5372*q -77.1212*q*q -51.3525*q*q*q -15.1692*q*q*q*q -0.981203*q*q*q*q*q; // H2-e [generally sub-dominant to H2-H+; can dominate if free e- largely from other sources (e.g. Mg, etc.), but in those conditions essentially impossible for H2 cooling to dominate
            if(logT>2.30103) {logLambdaH2_e = -22.1903 + 1.5729*q -0.213351*q*q + 0.961498*q*q*q -0.910232*q*q*q*q + 0.137497*q*q*q*q*q;}
            Lambda_H2_thin += n_elec * X_Hfrac * pow(10., DMAX(logLambdaH2_e, -50.));

            double f_HD = DMIN(0.00126*f_molec , 4.0e-5*nH0); // ratio of HD molecules to H2 molecules: in low limit, HD easier to form so saturates at about 0.13% of H2 molecules, following Galli & Palla 1998, but obviously cannot exceed the cosmic ratio of D/H=4e-5
            double nH_over_ncrit = Lambda_H2_thin / Lambda_H2_thick , Lambda_HD = f_HD * Lambda_HD_thin / (1. + (f_HD/(f_molec+MIN_REAL_NUMBER))*nH_over_ncrit), Lambda_H2 = f_molec * Lambda_H2_thin / (1. + nH_over_ncrit); // correct cooling rates for densities above critical
            double Lambda_Metals_Neutral = nH0 * Lambda_Metals; // finally note our metal terms here are all for atomic or molecular, not ionized (handled in tables above)
            if(!isfinite(Lambda_Metals_Neutral) || Lambda_Metals_Neutral < 0) {Lambda_Metals_Neutral=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_H2) || Lambda_H2 < 0) {Lambda_H2=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_HD) || Lambda_HD < 0) {Lambda_HD=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            LambdaMol = Lambda_Metals_Neutral + Lambda_H2 + Lambda_HD; // combine to get total cooling rate
#endif
            LambdaMol *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            LambdaDust = gas_dust_heating_coeff(target,T,Tdust, pp, cell) * (T-Tdust);// Note our sign convention is such that positive lambda = gas cooling
#if !defined(GALSF_ISMDUSTCHEM_MODEL)
            if(T>3.e5) {double dx=(T-3.e5)/2.e5; LambdaDust *= exp(-DMIN(dx*dx,40.));} /* needs to truncate at high temperatures b/c of dust destruction (in some modules we solve for this explicitly - in that case can protect this more explicitly, but here, we will make a simple approximation, otherwise we run into problems. note this is not sublimation generally, but sputtering, that causes the destruction */
#endif
            LambdaDust *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            if(!isfinite(LambdaDust)) {LambdaDust=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(LambdaMol)) {LambdaMol=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            LambdaMol *= ((T-T_cmb_radeff)/(T+T_cmb_radeff)); // account (approximately) for the CMB temperature 'bath' (could more accurately subtract Lambda(T[cmb]), but that's an approximation as well that can give some odd results owing to not treating the solve for molecules indepedently there, so we use this form instead, which is generally good)
            if(LambdaMol > 0) {Lambda += LambdaMol;}
        }
#endif

        Heat = 0;  /* Now, collect heating terms */

#if ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION)) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
        // here we account for the fact that the local spectrum is softer than the UVB which includes AGN and is hardened by absorption within galaxies. we do this by simply lowering the effective heating rate [mean photon energy absorbed per ionization], which captures the leading-order effect //
        if(CoolTables.J_UV != 0) {Heat_Ion_from_UVB = shieldfac / nHcgs * ((nH0 * CoolTables.epsH0 + nHe0 * CoolTables.epsHe0 + nHep * CoolTables.epsHep) + CoolTables.gJH0*(local_gammamultiplier-1.)*(nH0*2.9 + nHe0*0.44 + nHep*4.2e-4)*1.6e-12);} // this assumes an approximately IMF-averaged mean O-star Teff ~ 40000 K -- note the weights here for this correspond to mean energy per H ionization, so for H is just some energy in eV, but for He is weighted by relative ionization rate: softer spectrum translates to steeper dropoff of these terms
#else
        if(CoolTables.J_UV != 0) {Heat_Ion_from_UVB = local_gammamultiplier * (nH0 * CoolTables.epsH0 + nHe0 * CoolTables.epsHe0 + nHep * CoolTables.epsHep) / nHcgs * shieldfac;} // shieldfac allows for self-shielding from background
#endif
        Heat += Heat_Ion_from_UVB;
        
#if defined(RT_DISABLE_UV_BACKGROUND)
        Heat = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if((target >= 0) && (nHcgs > MIN_REAL_NUMBER))
        {
            int k; double c_light_nH = C_LIGHT_CGS / (nHcgs * UNIT_LENGTH_IN_CGS) * UNIT_ENERGY_IN_CGS; // want physical cgs units for quantities below
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if(RT_BAND_IS_IONIZING(k))
                {
                    double n_gamma_tot = cell[target].rt_photon_number_density(k);
#ifdef RT_INFRARED
                    n_gamma_tot += rt_irband_egydensity_in_band(target,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k], cell) / (DMAX(All.rt_nu_eff_eV[RT_FREQ_BIN_H0],cell[target].Radiation_Temperature/2959.81)*ELECTRONVOLT_IN_ERGS/UNIT_ENERGY_IN_CGS);
#endif
                    double c_nH_time_n_photons_vol = c_light_nH * n_gamma_tot; // gives photon flux
                    double cross_section_ion, kappa_ion, dummy;
                    if(All.rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * All.rt_ion_sigma_HI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = All.rt_ion_G_HI[k] * cross_section_ion * c_nH_time_n_photons_vol;// (egy per photon x cross section x photon flux) :: attenuation factors [already in flux/energy update]: * slab_averaging_function(kappa_ion * Sigma_particle); // egy per photon x cross section x photon flux (w attenuation factors) // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                        Heat_Ion_from_RHD += dummy;
                    }
                    if(All.rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * All.rt_ion_sigma_HeI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = All.rt_ion_G_HeI[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion * Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                        Heat_Ion_from_RHD += dummy;
                    }
                    if(All.rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * All.rt_ion_sigma_HeII[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = All.rt_ion_G_HeII[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion*Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                        Heat_Ion_from_RHD += dummy;
                    }
                }
            }
        }
#endif

        Heat += CR_gas_heating(target, n_elec, nH0, nHcgs, pp, cell); // CR hadronic+Coulomb+ionization heating //
#if defined(COOL_LOW_TEMPERATURES)
        if(LambdaMol<0) {Heat -= LambdaMol;} // Molecular line heating (Trad_mol_cooling_batch > Tgas) //
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION)
        if(LambdaMetal<0) {Heat -= LambdaMetal;} // potential net heating from low-temperature gas-phase metal line absorption //
#endif
#endif
        if(LambdaCompton<0) {Heat -= LambdaCompton;} /* Compton heating rather than cooling */

#if defined(GALSF_FB_FIRE_RT_LONGRANGE) || defined(RT_PHOTOELECTRIC) || defined(RT_ISRF_BACKGROUND)
        /* Photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995); now with 'update' from Wolfire 2005 for PAH [fudge factor 0.5 below] */
        if((target >= 0) && (T < 1.0e6))
        {
            double photoelec = get_FUV_G0(target, shieldfac, 0, pp, cell);
            if(photoelec > 0)
            {
                LambdaPElec = -1.3e-24 * photoelec / nHcgs * (pp[target].Metallicity[0]/All.SolarAbundances[0]) * return_dust_to_metals_ratio_vs_solar(target,0, pp, cell); // negative sign for lambda b/c heating
                double x_photoelec = photoelec * sqrt(T) / (0.5 * (1.0e-12+n_elec) * nHcgs);
                LambdaPElec *= 0.049/(1+pow(x_photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+x_photoelec/5000.);
#if defined(OUTPUT_COOLRATE_DETAIL)
                if(target >= 0) {cell[target].PElecHeatingRate = -LambdaPElec;}
#endif
                Heat -= LambdaPElec;
            }
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>CoolTables.Tmax K */
    {
        /* at high T (fully ionized); only free-free and Compton cooling are present.  Assumes no heating. */
        Heat = LambdaExc = LambdaExcH0 = LambdaExcHep = LambdaIon = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRec = LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;
        nH0 = 0; nHp = 1.0; nHe0 = 0; nHep = 0; nHepp = yhelium(target, pp); n_elec = nHp + 2.0 * nHepp; /* very hot: H and He both fully ionized */
        *n_elec_eval = n_elec; /* save this value for the output cycle */

        LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * n_elec * (1. + sqrt(T/0.4e10)); // free-free (with simplified relativistic correction)
        LambdaCompton = evaluate_Compton_heating_cooling_rate(target,T,nHcgs,n_elec,shieldfac, cell); // Compton
        Lambda = LambdaFF + DMAX(LambdaCompton,0);
    }

#if defined(RT_INFRARED)
    if(target >= 0) { /* attempt to account for gas-phase absorption self-opacity to cooling radiation here in limit where optically-thick regions are poorly resolved: only valid in some limits, so not really general here */
        double gas_self_absorption_opacity = rt_kappa_adaptive_IR_band(target,T,T,-1,-1, pp, cell), surface_density_fromcenter = 0.5 * (cell[target].Density*All.cf_a3inv) * (pp[target].Get_Particle_Size()*All.cf_atime);
        double tau_self = gas_self_absorption_opacity * surface_density_fromcenter, fcorr = 1./(1.+tau_self*tau_self);
        Heat*=fcorr; Lambda*=fcorr; LambdaMetal*=fcorr; LambdaExc*=fcorr; LambdaRec*=fcorr; LambdaIon*=fcorr; LambdaPElec*=fcorr; LambdaFF*=fcorr; LambdaMol*=fcorr; LambdaDust*=fcorr; LambdaCompton*=fcorr;
    }
#endif


#if defined(RT_NUV)
    double Lambda_rad_NUV = LambdaMetal; // most of LambdaMetal coming out in the NUV, as we define it
    Lambda_rad_NUV += LambdaExc; // this represents gas kinetic energy lost to collisional excitation. but each is assumed to produce a cascade back to the ground state, which should emit. note collisional ionization is not included here since the thermal energy is lost to the ionization energy, not to radiation, and the recombination luminosity is tallied separately.
    //Lambda_rad_NUV += LambdaIon; // this is the collisional ionization, distinct from the recombination luminosity. per above, this isn't directly into radiation. if you assumed it did lead to recombination, could use instead of LambdaRec, though would then ignore all but collisional ionization equilibrium. we're usually assuming case B recombination (UV emitted photons re-absorbed), so we'll assume a cascade for these into NUV/optical and other bands [otherwise should be added to photo-ionizing band]. still ignore LambdaRec, because otherwise this will double-count the UV background
    Lambda_rad_NUV += LambdaRec * DMAX(1.-shieldfac,0.) * DMIN(1.,DMAX(0.,Heat_Ion_from_RHD/(Heat_Ion_from_UVB+Heat_Ion_from_RHD+MIN_REAL_NUMBER))); // recombination radiation -- needs to be rate-limited to avoid 2x-counting the UVB, which is what the shieldfac and local_gammamultiplier terms attempt to account for here
#if !defined(RT_PHOTOELECTRIC) // if this module is active, these photons are accounted for explicitly in the photoelectric bands
    Lambda_rad_NUV += LambdaPElec; // otherwise lump it in here as well since it overlaps this band (should deplete it appropriately)
#endif
#if !defined(RT_FREEFREE) // if this module is active, these photons are accounted for explicitly in the free-free bands
    if(logT >= 5.) {Lambda_rad_NUV += LambdaFF;} // this can be coming at a wide range of wavelengths of course, but if not explicitly followed, it will go in one bin or another depending on the gas temperature
#endif
    if(target>=0) {
        int which_bin_to_cool_to = RT_FREQ_BIN_NUV; // default to NUV bin
#if defined(RT_INFRARED)
        cell[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_NUV] = cell[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED] = 0; // set to nil before deciding if we will add radiation here, to avoid double-counting with previous-loop information
        if(cell[target].Radiation_Temperature > 1.e4) {which_bin_to_cool_to = RT_FREQ_BIN_INFRARED;} // our more-accurate effective/adaptive IR band is already covering these wavelengths, rather than do the noisy step of cooling to NUV, re-absorbing with less accurate opacities and down-grading to the grey-band, just dump directly to the grey-band
#endif
#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM) /* may want to make broader, but restrict to these flags for now */
        which_bin_to_cool_to = RT_FREQ_BIN_INFRARED; /* since we're using this here to represent much higher energies, it should always track this bin */
#endif
        cell[target].Lambda_RadiativeCooling_toRHDBins[which_bin_to_cool_to] += Lambda_rad_NUV;} // save this to be used later (include all misc terms that will appear in our NUV radiation umbrella)
#endif
    
#if defined(RT_INFRARED)
    double Lambda_rad_IR = LambdaMol; // include molecular line cooling, cold atomic cooling (in LambdaMol) and dust cooling here, since all coming out in the IR
    Lambda_rad_IR += LambdaCompton; // Compton cooling/heating influences things here as well. Need to be a bit careful because Compton cooling can, in very low-density regions, be off the CMB, which we aren't explicitly evolving. so may want to add a check here for this case. but in the limits of interest that term is small, so not so important here.
#if !defined(RT_FREEFREE) // if this module is active, these photons are accounted for explicitly in the free-free bands
    if(logT < 5.) {Lambda_rad_IR += LambdaFF;} // this can be coming at a wide range of wavelengths of course, but if not explicitly followed, it will go in one bin or another depending on the gas temperature
#endif
    if(target>=0) {cell[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED] += Lambda_rad_IR;} // save this to be used later (include all misc terms that will appear in our IR radiation umbrella)
#endif

#if defined(COOL_LOW_TEMPERATURES)
    /* now add the dust cooling/heating terms - for RT+cooling runs this must be performed after we know all other radiation cooling terms, 
     because we are performing a backward-Euler solve and need to determine the consistent IR radiation energy and Tdust at t+dt, solving
     the equations of gas+dust+radiation energy conservation with 0 dust heat capacity - note that Tdust is always in equilibrium in
     this system, and we should overwrite it with whatever we get here. */
#ifdef RT_INFRARED
    if(target >= 0) {LambdaDust = rt_ir_lambdadust(target, T, pp, cell);} /* This updates Dust_Temperature and Lambda_RadiativeCooling_toRHDBins */
#endif
    if(LambdaDust>0) {Lambda += LambdaDust;} /* add the -positive- Lambda-dust associated with cooling */
    if(LambdaDust<0) {Heat -= LambdaDust;} // Dust collisional heating (Tdust > Tgas) //
#endif
#ifdef RT_COOLING_DUST_ONLY
    if(LambdaDust > 0) {Lambda = LambdaDust; Heat = 0;}
    if(LambdaDust <= 0) {Heat = -LambdaDust; Lambda = 0;}
#endif
    double Q = Heat - Lambda;
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0) {cell[target].CoolingRate = Lambda; cell[target].HeatingRate = Heat;}
#endif
#if defined(COOL_LOW_TEMPERATURES) && !defined(COOL_LOWTEMP_THIN_ONLY)
    /* if we are in the optically thick limit, we need to modify the cooling/heating rates according to the appropriate limits;
        this flag does so by using a simple approximation. we consider the element as if it were a slab, with a column density
        calculated from the simulation properties and the Sobolev approximation. we then assume it develops an equilibrium internal
        temperature structure on a radiative diffusion timescale much faster than the dynamical time, and so the surface radiation
        from a photosphere can be simply related to the local density by the optical depth to infinity. the equations here follow
        Rafikov, 2007 (ApJ, 662, 642):
            denergy/dt/dArea = sigma*T^4 / fc(tau)
            fc(tau) = tau^eta + 1/tau (taking chi, phi~1; the second term describes the optically thin limit, which is calculated above
                more accurately anyways - that was just Kirchoff's Law; so we only need to worry about the first term)
            eta = 4*(gamma-1) / [gamma*(1+alpha+beta*(gamma-1)/gamma)], where gamma=real polytropic index, and alpha/beta follow
                an opacity law kappa=kappa_0 * P^alpha * T^beta. for almost all the regimes of interest, however, eta~1, which is also
                what is obtained for a convectively stable slab. so we will use this.
            now, this gives sigma*T^4/tau * Area_eff / nHcgs as the 'effective' cooling rate in our units of Heat or Lambda above.
                the nHcgs just puts it in the same volumetric terms. The Area_eff must be defined as ~m_particle/surface_density
                to have the same meaning for a slab as assumed in Rafikov (and to integrate correctly over all particles in the slab,
                if/when the slab is resolved). We estimate this in our usual fashion with the Sobolev-type column density
            tau = kappa * surface_density; we estimate kappa ~ 5 cm^2/g * (0.001+Z/Z_solar), as the frequency-integrated kappa for warm
                dust radiation (~150K), weighted by the dust-to-gas ratio (with a floor for molecular absorption). we could make this
                temperature-dependent, though, fairly easily - for this particular problem it won't make much difference
        This rate then acts as an upper limit to the net heating/cooling calculated above (restricts absolute value)
     */
    if( (nHcgs > 0.1) && (target >= 0) )  /* don't bother at very low densities, since youre not optically thick, and protect from target=-1 with GALSF_EFFECTIVE_EQS */
    {
        double surface_density = evaluate_NH_from_GradRho(cell[target].Gradients.Density,pp[target].KernelRadius,cell[target].Density,pp[target].NumNgb,1,target,pp);
        surface_density *=  UNIT_SURFDEN_IN_CGS; // converts to cgs
        double effective_area = 2.3 * PROTONMASS_CGS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
	
	    kappa_eff = rt_kappa_adaptive_IR_band(target,T,T,0,1, pp, cell) / UNIT_SURFDEN_IN_CGS; // will return simple opacity law kappa = 0.1cm^2/g (T/10K)^2, capped at 5 cm^2/g, in code units [convert to physical here]
        if(kappa_eff < 0.1) {kappa_eff=0.1;}
        if(T>1500.){
            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * n_elec; /* Thompson scattering (non-relativistic) */
            double k_molecular = 0.1 * pp[target].Metallicity[0]; /* molecular line opacities */

            double k_Hminus = 1.1e-25 * sqrt(pp[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (pp[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            double k_conductive = 2.6e-7 * n_elec * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
            kappa_eff += 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif

#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0) {cell[target].NetHeatingRateQ = Q;}
#endif
#ifdef OUTPUT_MOLECULAR_FRACTION
    if(target>=0) {cell[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(CoolTables.gJH0/2.29e-10), pp, cell);}
#endif

#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) {if(cell[target].CoolingIsOperatorSplitThisTimestep==0) {Q += cell[target].DtInternalEnergy / nHcgs;}}
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target >= 0) {cell[target].HydroHeatingRate = cell[target].DtInternalEnergy / nHcgs;}
#endif
#endif
    return Q;
} // ends CoolingRate
/* ---- END device-compilable cooling functions (do_the_cooling_for_particle through CoolingRate) ---- */



int InitCoolMemory(void)
{
    /* With Kokkos+CUDA (Kokkos+GPU), allocate cooling table arrays in CUDA unified
       memory via Kokkos so they are accessible from device kernels. The pointer variables
       themselves are __managed__ (see top of file), so assigning to them here is fine.
       Without GPU offload, use the normal mymalloc pool. */
    /* A table that could not be had is reported and leaves the pointer null; the
     * caller stops before MakeCoolingTable would fill it. Every rank asks for the
     * same sizes, so they reach the stop together. */
    int all_tables_obtained = 1;
#define COOLMEM(name, type, n) do { \
        const size_t _bytes = (size_t)(n) * sizeof(type); \
        name = (type *) gizmo_gpu_alloc_shared(_bytes, #name); \
        if(!name) { \
            all_tables_obtained = 0; \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), \
                     "InitCoolMemory: could not allocate cooling table %s (%.1f MB)", \
                     #name, (double) _bytes / (1024.0 * 1024.0)); \
            gizmo_request_controlled_stop(7717, _msg, __FILE__, __LINE__, __FUNCTION__); \
        } \
    } while(0)
    COOLMEM(CoolTables.BetaH0,    double, NCOOLTAB + 1);
    COOLMEM(CoolTables.BetaHep,   double, NCOOLTAB + 1);
    COOLMEM(CoolTables.AlphaHp,   double, NCOOLTAB + 1);
    COOLMEM(CoolTables.AlphaHep,  double, NCOOLTAB + 1);
    COOLMEM(CoolTables.Alphad,    double, NCOOLTAB + 1);
    COOLMEM(CoolTables.AlphaHepp, double, NCOOLTAB + 1);
    COOLMEM(CoolTables.GammaeH0,  double, NCOOLTAB + 1);
    COOLMEM(CoolTables.GammaeHe0, double, NCOOLTAB + 1);
    COOLMEM(CoolTables.GammaeHep, double, NCOOLTAB + 1);
    COOLMEM(CoolTables.Betaff,    double, NCOOLTAB + 1);
#ifdef COOL_METAL_LINES_BY_SPECIES
    long i_nH=41; long i_T=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES;
    COOLMEM(CoolTables.SpCoolTable0, float, kspecies*i_nH*i_T);
    if(All.ComovingIntegrationOn) { COOLMEM(CoolTables.SpCoolTable1, float, kspecies*i_nH*i_T); }
#endif
#undef COOLMEM
    return all_tables_obtained;
}



void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i; double T,Tfact;
    if(All.MinGasTemp > 0.0) {CoolTables.Tmin = log10(All.MinGasTemp);} else {CoolTables.Tmin=-1.0;} // set minimum temperature in this table to some very low value if zero, where none of the cooling approximations above make sense
    CoolTables.deltaT = (CoolTables.Tmax - CoolTables.Tmin) / NCOOLTAB;
    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB; i++)
    {
        CoolTables.BetaH0[i] = CoolTables.BetaHep[i] = CoolTables.Betaff[i] = CoolTables.AlphaHp[i] = CoolTables.AlphaHep[i] = CoolTables.AlphaHepp[i] = CoolTables.Alphad[i] = CoolTables.GammaeH0[i] = CoolTables.GammaeHe0[i] = CoolTables.GammaeHep[i] = 0;
        T = pow(10.0, CoolTables.Tmin + CoolTables.deltaT * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
        if(118348. / T < 70.) {CoolTables.BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;}
        if(473638. / T < 70.) {CoolTables.BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;}

        CoolTables.Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
        //CoolTables.AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
        //CoolTables.AlphaHep[i] = 1.5e-10 * pow(T, -0.6353); /* old Cen92 fit */
        //CoolTables.AlphaHepp[i] = 4. * CoolTables.AlphaHp[i];	/* old Cen92 fit */
        CoolTables.AlphaHp[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        CoolTables.AlphaHep[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        CoolTables.AlphaHepp[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*alphaHp[1,T/Z^2] */

        if(470000.0 / T < 70) {CoolTables.Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));}
        if(157809.1 / T < 70) {CoolTables.GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;}
        if(285335.4 / T < 70) {CoolTables.GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;}
        if(631515.0 / T < 70) {CoolTables.GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;}
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables(void)
{
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables(0);
    }
}

void ReadMultiSpeciesTables(int iT)
{
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES; long i,j,k,r;
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;

    fname=GetMultiSpeciesFilename(iT,0);
    if(ThisTask == 0) printf(" ..opening Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read species cooling table in file `%s'\n", fname); endrun(456); return;}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&CoolTables.SpCoolTable0[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached Cooling EOF! \n");
                }
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename(iT,&fname,1);
     if(!(fdcool = fopen(fname, "r"))) {
     printf(" Cannot read species (He) cooling table in file `%s'\n", fname); endrun(456);}
     for(i=0;i<2;i++)
     for(j=0;j<i_nH;j++)
     for(k=0;k<i_Temp;k++)
     for(l=0;l<i_He;l++)
     fread(&SpCoolTable0_He[i][j][k][l],sizeof(float),1,fdcool);
     fclose(fdcool);
     */
    if (All.ComovingIntegrationOn && i<48) {
        fname=GetMultiSpeciesFilename(iT+1,0);
        if(ThisTask == 0) printf(" ..opening (z+) Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read species 1 cooling table in file `%s'\n", fname); endrun(456); return;}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&CoolTables.SpCoolTable1[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached Cooling EOF! \n");
                    }
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename(iT+1,&fname,1);
         if(!(fdcool = fopen(fname, "r"))) {
         printf(" Cannot read species 1 (He) cooling table in file `%s'\n", fname); endrun(456);}
         for(i=0;i<2;i++)
         for(j=0;j<i_nH;j++)
         for(k=0;k<i_Temp;k++)
         for(l=0;l<i_He;l++)
         fread(&SpCoolTable1_He[i][j][k][l],sizeof(float),1,fdcool);
         fclose(fdcool);
         */
    }
}

char *GetMultiSpeciesFilename(int i, int hk)
{
    static char fname[DEFAULT_PATH_BUFFERSIZE_TOUSE];
    if(i<0) i=0; if(i>48) i=48;
    if(hk==0) {
        snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE,"./spcool_tables/spcool_%d",i);
    } else {
        snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE,"./spcool_tables/spcool_He_%d",i);
    }
    return fname;
}

#endif






/* table input (from file TREECOOL) for ionizing parameters */
#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 250		/* Max # of lines in TREECOOL -- needs to be at least one larger than number of non-zero lines */
static float inlogz[TABLESIZE];
static double gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE]; // upgrade from float to double, should read fine
static double eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE]; // upgrade from float to double, should read fine
static int nheattab;		/* length of table */


void ReadIonizeParams(const char *fname)
{
    int i; FILE *fdcool;
    if(!(fdcool = fopen(fname, "r"))) {printf(" Cannot read ionization table in file `%s'. Make sure the correct TREECOOL file is placed in the code run-time directory, and that any leading comments (e.g. lines preceded by ##) are deleted from the file.\n", fname); endrun(456); return;}
    for(i=0; i<TABLESIZE; i++) {inlogz[i]=100; gH0[i]=0; gHe[i]=0; gHep[i]=0; eH0[i]=0; eHe[i]=0; eHep[i]=0;}
    for(i=0; i<TABLESIZE; i++) {if(fscanf(fdcool, "%g %lg %lg %lg %lg %lg %lg", &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF) {break;}}
    fclose(fdcool);
    for(i=0, nheattab=0; i<TABLESIZE; i++) {if(gH0[i] != 0.0) {nheattab++;} else {break;}} /*  nheattab is the number of entries in the table */
    if(ThisTask == 0) printf(" ..read ionization table [TREECOOL] with %d non-zero UVB entries in file `%s'. Make sure to cite the authors from which the UV background was compiled! (See user guide for the correct references).\n", nheattab, fname);
}


void IonizeParams(void)
{
    IonizeParamsTable();
}



void IonizeParamsTable(void)
{
    int i, ilow, ihi;
    double logz, dzlow, dzhi;
    double redshift;

    if(All.ComovingIntegrationOn)
        {redshift = 1 / All.Time - 1;}
    else
    {
        /* in non-cosmological mode, still use, but adopt z=0 background */
#ifdef RT_ISRF_BACKGROUND
	redshift = All.RadiationBackgroundRedshift;
#else
        redshift = 0;
#endif
        /*
         CoolTables.gJHe0 = CoolTables.gJHep = CoolTables.gJH0 = CoolTables.epsHe0 = CoolTables.epsHep = CoolTables.epsH0 = CoolTables.J_UV = 0;
         return;
         */
    }

    logz = log10(redshift + 1.0);
    ilow = 0;
    if(nheattab <= 0) {CoolTables.gJHe0 = CoolTables.gJHep = CoolTables.gJH0 = CoolTables.epsHe0 = CoolTables.epsHep = CoolTables.epsH0 = CoolTables.J_UV = 0; return;}
    for(i=0; i<nheattab; i++) {if(inlogz[i] < logz) {ilow = i;} else {break;}}
    ihi = i + 1;
    if(ilow >= nheattab) {ihi = i;}
    dzlow = logz - inlogz[ilow];
    dzhi = inlogz[ihi] - logz;
    if((logz > inlogz[nheattab - 1]) || (gH0[ilow] == 0) || (gH0[ihi] == 0) || (ilow > nheattab))
    {
        CoolTables.gJHe0 = CoolTables.gJHep = CoolTables.gJH0 = CoolTables.epsHe0 = CoolTables.epsHep = CoolTables.epsH0 = CoolTables.J_UV = 0;
        return;
    }
    else {CoolTables.J_UV = 1.e-21;}		/* irrelevant as long as it's not 0 */

    CoolTables.gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ihi])) / (dzlow + dzhi));
    CoolTables.gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ihi])) / (dzlow + dzhi));
    CoolTables.gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ihi])) / (dzlow + dzhi));
    CoolTables.epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ihi])) / (dzlow + dzhi));
    CoolTables.epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ihi])) / (dzlow + dzhi));
    CoolTables.epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ihi])) / (dzlow + dzhi));

    return;
}


void SetZeroIonization(void)
{
    CoolTables.gJHe0 = CoolTables.gJHep = CoolTables.gJH0 = 0; CoolTables.epsHe0 = CoolTables.epsHep = CoolTables.epsH0 = 0; CoolTables.J_UV = 0;
}


void IonizeParamsFunction(void)
{
    int i, nint;
    double a0, planck, ev, e0_H, e0_He, e0_Hep;
    double gint, eint, t, tinv, fac, eps;
    double at, beta, s;
    double pi;

#define UVALPHA         1.0
    double Jold = -1.0;
    double redshift;

    CoolTables.J_UV = 0.;
    CoolTables.gJHe0 = CoolTables.gJHep = CoolTables.gJH0 = 0.;
    CoolTables.epsHe0 = CoolTables.epsHep = CoolTables.epsH0 = 0.;


    if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
        redshift = 1 / All.Time - 1;

        if(redshift >= 6) {CoolTables.J_UV = 0.;}
        else
        {
            if(redshift >= 3) {CoolTables.J_UV = 4e-22 / (1 + redshift);}
            else
            {
                if(redshift >= 2) {CoolTables.J_UV = 1e-22;}
                else {CoolTables.J_UV = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);}
            }
        }
        if(CoolTables.J_UV == Jold) {return;}
        Jold = CoolTables.J_UV;
        if(CoolTables.J_UV == 0) {return;}


        a0 = 6.30e-18;
        planck = 6.6262e-27;
        ev = 1.6022e-12;
        e0_H = 13.6058 * ev;
        e0_He = 24.59 * ev;
        e0_Hep = 54.4232 * ev;

        gint = 0.0;
        eint = 0.0;
        nint = 5000;
        at = 1. / ((double) nint);

        for(i = 1; i <= nint; i++)
        {
            t = (double) i;
            t = (t - 0.5) * at;
            tinv = 1. / t;
            eps = sqrt(tinv - 1.);
            fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA + 3.);
            gint += fac * at;
            eint += fac * (tinv - 1.) * at;
        }

        CoolTables.gJH0 = a0 * gint / planck;
        CoolTables.epsH0 = a0 * eint * (e0_H / planck);
        CoolTables.gJHep = CoolTables.gJH0 * pow(e0_H / e0_Hep, UVALPHA) / 4.0;
        CoolTables.epsHep = CoolTables.epsH0 * pow((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

        at = 7.83e-18;
        beta = 1.66;
        s = 2.05;

        CoolTables.gJHe0 = (at / planck) * pow((e0_H / e0_He), UVALPHA) *
        (beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
        CoolTables.epsHe0 = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA) *
        (beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

        pi = M_PI;
        CoolTables.gJH0 *= 4. * pi * CoolTables.J_UV;
        CoolTables.gJHep *= 4. * pi * CoolTables.J_UV;
        CoolTables.gJHe0 *= 4. * pi * CoolTables.J_UV;
        CoolTables.epsH0 *= 4. * pi * CoolTables.J_UV;
        CoolTables.epsHep *= 4. * pi * CoolTables.J_UV;
        CoolTables.epsHe0 *= 4. * pi * CoolTables.J_UV;
    }
}
#endif // !(CHIMES)


/* Check for required cooling data files and provide them from defaults if missing */
static bool file_exists(const char *path) {FILE *f = fopen(path, "r"); if(f) {fclose(f); return true;} return false;}

static void rank0_run_or_die(const char *cmd, const char *errmsg) {
    if(system(cmd) != 0) {printf(" *** %s ***\n", errmsg); endrun(456);}
}

static void ensure_cooling_tables_exist(void) {
    if(ThisTask == 0) {
        if(!file_exists("TREECOOL")) {
            std::string treecool_src = std::string(GIZMO_SOURCE_DIR) + "../data/cooling/TREECOOL";
            printf(" *** WARNING: TREECOOL not found in run directory. Copying default from %s — make sure this is the table you want! ***\n", treecool_src.c_str());
            rank0_run_or_die(("cp '" + treecool_src + "' TREECOOL").c_str(), "Could not copy TREECOOL from source tree. Place it in the run directory manually.");
        }
#ifdef COOL_METAL_LINES_BY_SPECIES
        if(!file_exists(GetMultiSpeciesFilename(0, 0))) {
            const char *url = "https://users.flatironinstitute.org/~mgrudic/gizmo_tests/spcool_tables.tgz";
            printf(" *** WARNING: spcool_tables not found in run directory. Downloading from %s — make sure these are the tables you want! ***\n", url);
            rank0_run_or_die(("curl -fsSL '" + std::string(url) + "' | tar xz").c_str(), "Download/extraction of spcool_tables failed. Place spcool_tables/ in the run directory manually.");
        }
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void InitCool(void)
{
    if(ThisTask == 0) printf("Initializing cooling ...\n");

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

#ifdef COOL_GRACKLE
    InitGrackle();
#endif

#ifdef CHIMES
    snprintf(ChimesGlobalVars.MainDataTablePath, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/chimes_main_data.hdf5", ChimesDataPath);
    snprintf(ChimesGlobalVars.EqAbundanceTablePath, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/EqAbundancesTables/%s", ChimesDataPath, ChimesEqAbundanceTable);
    snprintf(ChimesGlobalVars.PhotoIonTablePath[0], DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/%s", ChimesDataPath, ChimesPhotoIonTable);

    // By default, use 1 spectrum, unless
    // stellar fluxes are enabled.
    ChimesGlobalVars.N_spectra = 1;

#ifdef CHIMES_STELLAR_FLUXES
    ChimesGlobalVars.N_spectra += CHIMES_LOCAL_UV_NBINS;

    int spectrum_idx;
    for (spectrum_idx = 1; spectrum_idx < ChimesGlobalVars.N_spectra; spectrum_idx++)
      snprintf(ChimesGlobalVars.PhotoIonTablePath[spectrum_idx], DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/starburstCrossSections/cross_sections_SB%d.hdf5", ChimesDataPath, spectrum_idx);
#endif

    if (ChimesUVBMode > 0)
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = 0;
	if (ChimesUVBMode == 1)
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
	else
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 1;
      }
    else
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = -1;
	ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
      }

    // Hybrid cooling has not yet been
    // implemented in Gizmo. Switch
    // it off for now.
    ChimesGlobalVars.hybrid_cooling_mode = 0;

    // Set the chimes_exit() function
    // to use the Gizmo-specific
    // chimes_gizmo_exit().
    chimes_exit = &chimes_gizmo_exit;

    // Initialise the CHIMES module.
    init_chimes(&ChimesGlobalVars);

#ifdef CHIMES_METAL_DEPLETION
    chimes_init_depletion_data();
#endif

#else // CHIMES
    if(!InitCoolMemory()) { return; }   /* tables absent: the stop is already requested */
    MakeCoolingTable();
    ensure_cooling_tables_exist();
    ReadIonizeParams("TREECOOL");
    IonizeParams();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables();
#endif
    /* Cooling tables are allocated via Kokkos SharedSpace (__managed__ / UVM),
       so they are automatically accessible from device. No explicit
       #pragma omp target enter data needed (and the CoolTables #define aliases
       are incompatible with the OMP map syntax). */
#endif // CHIMES
}



#ifndef CHIMES
#ifdef COOL_METAL_LINES_BY_SPECIES
/* Metal-line table lookup — KOKKOS_FUNCTION so device kernels can call it */
KOKKOS_FUNCTION
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z, int target)
{
    double ne_over_nh_tbl=1, Lambda=0;
    int k, N_species_active = (int)NUM_LIVE_SPECIES_FOR_COOLTABLES;

    /* pre-calculate the indices for density and temperature, then we just need to call the tables by species */
    int ixmax=40, iymax=175;
    int ix0, iy0, ix1, iy1;
    double dx, dy, dz, mdz;
    long i_T=iymax+1, inHT=i_T*(ixmax+1);
    if(All.ComovingIntegrationOn && All.SpeciesTableInUse<48) {dz=log10(1/All.Time)*48; dz=dz-(int)dz; mdz=1-dz;} else {dz=0; mdz=1;}

    /* A density that is not positive and finite, or a non-finite temperature,
     * cannot index these tables, and it means the cell reached cooling with
     * invalid state -- the real error is wherever that state was produced.
     * Report the values, which are the evidence for finding it, and stop.
     * Catching it here is not optional: log10 of a negative density is a NaN,
     * every comparison against a NaN is false, so a NaN would pass BOTH bounds
     * below and convert to a large negative table index. */
    if(!((nHcgs > 0) && isfinite(nHcgs) && isfinite(logT)))
    {
        printf("GetCoolingRateWSpecies: cell=%d reached the metal-line tables with unusable state "
               "nHcgs=%g logT=%g -- invalid upstream, stopping\n", target, nHcgs, logT);
        endrun(90001019);
        return 0;
    }

    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    /* written so that a NaN, if one ever reaches here, lands on 0 rather than
     * passing through: !(x>=0) is true for a NaN where x<0 is false */
    if(!(dx>=0)) {dx=0;} else {if(dx>ixmax) {dx=ixmax;}}
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) {ix1=ixmax;}
    dx=dx-ix0;
    if(!(dy>=0)) {dy=0;} else {if(dy>iymax) {dy=iymax;}}
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) {iy1=iymax;}
    dy=dy-iy0;
    long index_x0y0=iy0+ix0*i_T, index_x0y1=iy1+ix0*i_T, index_x1y0=iy0+ix1*i_T, index_x1y1=iy1+ix1*i_T;

    ne_over_nh_tbl = GetLambdaSpecies(0,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz);
    if(ne_over_nh_tbl > 0)
    {
        double zfac = 0.0127 / All.SolarAbundances[0];
        for (k=1; k<N_species_active; k++)
        {
            long k_index = k * inHT;
            Lambda += GetLambdaSpecies(k_index,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz) * Z[k+1]/(All.SolarAbundances[k+1]*zfac);
        }
        Lambda /= ne_over_nh_tbl;
    }
    return Lambda;
}


KOKKOS_FUNCTION
double GetLambdaSpecies(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz)
{
    long x0y0 = index_x0y0 + k_index;
    long x0y1 = index_x0y1 + k_index;
    long x1y0 = index_x1y0 + k_index;
    long x1y1 = index_x1y1 + k_index;
    double i1, i2, j1, j2, w1, w2, u1;
    i1 = CoolTables.SpCoolTable0[x0y0];
    i2 = CoolTables.SpCoolTable0[x0y1];
    j1 = CoolTables.SpCoolTable0[x1y0];
    j2 = CoolTables.SpCoolTable0[x1y1];
    if(dz > 0)
    {
        i1 = mdz * i1 + dz * CoolTables.SpCoolTable1[x0y0];
        i2 = mdz * i2 + dz * CoolTables.SpCoolTable1[x0y1];
        j1 = mdz * j1 + dz * CoolTables.SpCoolTable1[x1y0];
        j2 = mdz * j2 + dz * CoolTables.SpCoolTable1[x1y1];
    }
    w1 = i1*(1-dy) + i2*dy;
    w2 = j1*(1-dy) + j2*dy;
    u1 = w1*(1-dx) + w2*dx;
    return u1;
}
#endif // COOL_METAL_LINES_BY_SPECIES
#endif // !(CHIMES)


#ifdef GALSF_FB_FIRE_RT_LONGRANGE
void selfshield_local_incident_uv_flux(void)
{   /* include local self-shielding with the following */
    const int n_act = (int)ActiveParticleList.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int idx = 0; idx < n_act; ++idx)
    {
        int i = ActiveParticleList[idx];
        if(P[i].Type==0)
        {
            if((CellP[i].Rad_Flux_UV>0) && (P[i].KernelRadius>0) && (CellP[i].Density>0) && (P[i].Mass>0) && (All.Time>0))
            {
                CellP[i].Rad_Flux_UV *= UNIT_FLUX_IN_CGS * 1276.19; CellP[i].Rad_Flux_EUV *= UNIT_FLUX_IN_CGS * 1276.19; // convert to Habing units (normalize strength to local MW field in this [narrow] band, so not the 'full' Habing flux)
                double surfdensity = evaluate_NH_from_GradRho(P[i].GradRho,P[i].KernelRadius,CellP[i].Density,P[i].NumNgb,1,i); // in CGS
                double tau_nuv = rt_kappa(i,RT_FREQ_BIN_FIRE_UV, P, CellP) * surfdensity; // optical depth: this part is attenuated by dust //
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
                tau_nuv *= (1.0e-3 + (P[i].Metallicity[0]/All.SolarAbundances[0])*return_dust_to_metals_ratio_vs_solar(i,0, P, CellP)); // if using older FIRE defaults, this was manually added instead of rolled into rt_kappa -- annoying but here for completeness //
#endif
                double tau_euv = 3.7e6 * surfdensity * UNIT_SURFDEN_IN_CGS; // optical depth: 912 angstrom kappa_euv: opacity from neutral gas //
                CellP[i].Rad_Flux_UV *= exp(-DMIN(tau_nuv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                //CellP[i].Rad_Flux_UV *= 0.01 + 0.99/(1.0 + 0.8*tau_nuv + 0.85*tau_nuv*tau_nuv); // attenuate (for clumpy medium with hard-minimum 1% scattering)
                //CellP[i].Rad_Flux_EUV *= exp(-DMIN(tau_euv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                CellP[i].Rad_Flux_EUV *= 0.01 + 0.99/(1.0 + 0.8*tau_euv + 0.85*tau_euv*tau_euv); // attenuate (for clumpy medium with hard-minimum 1% scattering) //
            } else {CellP[i].Rad_Flux_UV = CellP[i].Rad_Flux_EUV = 0;}
#if ((GALSF_FB_FIRE_STELLAREVOLUTION > 2) || !defined(GALSF_FB_FIRE_STELLAREVOLUTION)) && defined(GALSF_FB_FIRE_RT_HIIHEATING) && !defined(CHIMES_HII_REGIONS)
            if(CellP[i].DelayTimeHII > 0)
            {   /* assign typical strong HII region flux + enough flux to maintain cell fully-ionized, regardless (x'safety-factor') */
                double n1000 = CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / 1000.; // density in 1000 cm^-3
                double flux_compactHII = DMAX(0.85*pow(n1000,1./3.) , 1) * 2.6e5*n1000; // set to typical value in HII region or minimum needed to maintain f_neutral < 1e-5-ish, whichever is larger
                CellP[i].Rad_Flux_UV += flux_compactHII; CellP[i].Rad_Flux_EUV += flux_compactHII;
            }
#endif
        }
    }
}
#endif




/* subroutine to update the molecular fraction using our implicit solver for a simple --single-species-- network (just H2) */
KOKKOS_FUNCTION
void update_explicit_molecular_fraction(int i, double dtime_cgs, struct particle_data *pp, struct gas_cell_data *cell)
{
    if(dtime_cgs <= 0) {return;}
#ifdef COOL_MOLECFRAC_NONEQM
    // first define a number of environmental variables that are fixed over this update step
    double fH2_initial = cell[i].MolecularMassFraction_perNeutralH; // initial molecular fraction per H atom, entering this subroutine, needed for update below
    double xn_e=1, nh0=0, nHe0=1, nHepp=0, nhp=0, nHep=0, temperature, mu_meanwt=1, rho=cell[i].Density*All.cf_a3inv, u0=cell[i].InternalEnergy;
    /* this runs inside the cooling solver, where the helium ionization state is wanted at full
       accuracy for the molecular rates below, so solve for it here as this always has [assumes a
       fixed value of fH2 at the previous update value] */
    xn_e = cell[i].Ne; nh0 = cell[i].HI;
    temperature = convert_u_to_temp(u0*UNIT_SPECEGY_IN_CGS, rho*UNIT_DENSITY_IN_CGS, i, &xn_e, &nh0, &nhp, &nHe0, &nHep, &nHepp, &mu_meanwt, pp, cell);
    double T=1, f_dustgas_solar=1, urad_G0=1, xH0=1, x_e=0, nH_cgs=rho*UNIT_DENSITY_IN_NHCGS; // initialize definitions of some variables used below to prevent compiler warnings
    f_dustgas_solar=1; urad_G0=1; // initialize metal/dust and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
#ifdef RT_ISRF_BACKGROUND
    urad_G0 = All.InterstellarRadiationFieldStrength;
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
    if(cell[i].DelayTimeHII > 0) {cell[i].MolecularMassFraction_perNeutralH=cell[i].MolecularMassFraction=0; return;} // force gas flagged as in HII regions to have zero molecular fraction
#endif
    if(temperature > 3.e5) {cell[i].MolecularMassFraction_perNeutralH=cell[i].MolecularMassFraction=0; return;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(nh0, 0.),1.); // get neutral fraction [given by call to this program]
    if(xH0 <= MIN_REAL_NUMBER) {cell[i].MolecularMassFraction_perNeutralH=cell[i].MolecularMassFraction=0; return;} // no neutral gas, no molecules!
    x_e = DMIN(DMAX(xn_e, 0.),2.); // get free electron ratio [number per H nucleon]
    double log_T=log10(T), ln_T=log(T), gamma_12=return_local_gammamultiplier(i, cell)*CoolTables.gJH0/1.0e-12, shieldfac=return_uvb_shieldfac(i,gamma_12,nH_cgs,log_T, cell), urad_from_uvb_in_G0=sqrt(shieldfac)*(CoolTables.gJH0/2.29e-10); // estimate UVB contribution if we have partial shielding, to full photo-dissociation rates //
#ifdef METALS
    f_dustgas_solar=(pp[i].Metallicity[0]/All.SolarAbundances[0])*return_dust_to_metals_ratio_vs_solar(i,0, pp, cell); // this is only used for the dust-phase formation rates below, so just the dust term here
#endif
    /* get incident radiation field from whatever module we are using to track it */
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    urad_G0 = DMAX(cell[i].Rad_Flux_UV, 1.e-10); // note this is ALREADY self-shielded by dust, so we need to be careful about 2x-counting the self-shielding approximation below; hence limit this to a rather sizeable value  //
#endif
#if defined(RT_PHOTOELECTRIC) || defined(RT_LYMAN_WERNER)
    int whichbin = RT_FREQ_BIN_LYMAN_WERNER;
#if !defined(RT_LYMAN_WERNER)
    whichbin = RT_FREQ_BIN_PHOTOELECTRIC; // use photo-electric bin as proxy (very close) if don't evolve LW explicitly
#endif
    urad_G0 = cell[i].Rad_E_gamma[whichbin] * (cell[i].Density*All.cf_a3inv/cell[i].Mass) * UNIT_EGY_DENSITY_IN_HABING; // convert to Habing field //
#endif
    urad_G0 += urad_from_uvb_in_G0; // include whatever is contributed from the meta-galactic background, fed into this routine
    urad_G0 = DMIN(DMAX( urad_G0 , 1.e-10 ) , 1.e10 ); // limit values, because otherwise exponential self-shielding approximation easily artificially gives 0 incident field
#ifdef RT_INFRARED
    urad_G0 += rt_irband_egydensity_in_band(i,11.2,500., cell) * UNIT_EGY_DENSITY_IN_HABING; // add contribution from the adaptive band
#endif
    // define a number of variables needed in the shielding module
    double dx_cell = pp[i].Get_Particle_Size() * All.cf_atime; // cell size
    double surface_density_H2_0 = 5.e14 * PROTONMASS_CGS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
    w0 = 0.035; // actual calibration from Drain, Gnedin, Richings, others: 0.2 is more appropriate as a re-calibration for sims doing local eqm without ability to resolve shielding at higher columns
    //double surface_density_local = xH0 * cell[i].Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. note G0 is -already- attenuated in the pre-processing by dust.
    double surface_density_local = xH0 * evaluate_NH_from_GradRho(pp[i].GradRho,pp[i].KernelRadius,cell[i].Density,pp[i].NumNgb,1,i,pp) * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth to infinity with our Sobolev-type approximation. Note G0 is already attenuated by dust, but we need to include H2 self-shielding, for which it is appropriate to integrate to infinity.
    double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
    double gradv = cell[i].velocity_gradient_norm();
    double dv_turb=gradv*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
    double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, fH2_min, Q_max, Q_min, Q_initial; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]
    double b_time_Mach = 0.5 * dv_turb / (v_thermal_rms/sqrt(3.)); // cs_thermal for molecular [=rms v_thermal / sqrt(3)], dv_turb to full inside dx, assume "b" prefactor for compressive-to-solenoidal ratio corresponding to the 'natural mix' = 0.5. could further multiply by 1.58 if really needed to by extended dvturb to 2h = H, and vthermal from molecular to atomic for the generating field, but not as well-justified
#if defined(GALSF_ISMDUSTCHEM_MODEL) && defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
    cell[i].ISMDustChem_MachNumber = dv_turb / (v_thermal_rms/sqrt(3.)); // dust evolution model uses local mach number/clumping factors for numerous calculations
#endif
    double clumping_factor = 1. + b_time_Mach*b_time_Mach; // this is the exact clumping factor for a standard lognormal PDF with S=ln[1+b^2 Mach^2] //
    double clumping_factor_3 = clumping_factor*clumping_factor*clumping_factor; // clumping factor N for <rho^n>/<rho>^n = clumping factor^(N*(N-1)/2) //
    
    /* evolve dot[nH2]/nH0 = d_dt[fH2[neutral]] = (1/nH0) * (a_Z*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] ) */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, EXPmax=90.; int iter=0; // define some variables for below, including neutral H number density, free electron number, etc.
    double x_p = DMIN(DMAX(nhp , x_e/10.), 2.); // get free H+ fraction [cap because irrelevant to below in very low regime //
    /* use interpolation function from Glover & Abel 2008 [GA08], section 2.1.3, for interpolating between ground state (v=0) and LTE assumptions for states for collisional dissociation rates */
    double XH=HYDROGEN_MASSFRAC, xH2_guess=XH*DMAX(DMIN(cell[i].MolecularMassFraction,1.),0.), xH_guess=DMAX(XH-xH2_guess,0), xHe_guess=nHe0+nHep+nHepp;
    double logT4=log_T-4., ncr_H=pow(10.,3.0-0.416*logT4-0.327*logT4*logT4), ncr_H2=pow(10.,4.845-1.3*logT4+1.62*logT4*logT4), ncr_He=pow(10.,5.0792*(1.-1.23e-5*(T-2000.)));
    double ncrit = 1./(xH_guess/ncr_H + xH2_guess/ncr_H2 + xHe_guess/ncr_He), n_ncrit=nH_cgs/ncrit, f_v0_LTE=1./(1.+n_ncrit), f_LTE_v0 = 1.-f_v0_LTE;

    double b_H2Hp = DMAX(0., -3.3232183e-7 + 3.3735382e-7*ln_T -1.4491368e-7*ln_T*ln_T + 3.4172805e-8*ln_T*ln_T*ln_T -4.7813720e-9*ln_T*ln_T*ln_T*ln_T +3.9731542e-10*ln_T*ln_T*ln_T*ln_T*ln_T -1.8171411e-11*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T +3.5311932e-13*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T) * exp(-DMIN(21237.15/T,EXPmax)) * (nhp*nH_cgs) * clumping_factor; // H2-H+ dissociation, GA08-TableA1-7; note their expression (GA08) has an error where they write log[T] but this gives unphysical values. it should be ln[T], as it is correctly written in the original Savin et al. 2004 paper from which this fitting function is taken
    double b_H2e_v0 = 4.49e-9 * pow(T,0.11) * exp(-DMIN(101858./T,EXPmax)), b_H2e_LTE = 1.91e-9 * pow(T,0.136) * exp(-DMIN(53407.1/T,EXPmax)), b_H2e = pow(10., f_v0_LTE*log10(b_H2e_v0) + f_LTE_v0*log10(b_H2e_LTE)) * (x_e*nH_cgs) * clumping_factor; // collisional H2-e- dissociation; GA08-A1-8
    double b_H2HI_v0 = 6.67e-12 * sqrt_T * exp(-DMIN(1.+63593./T,EXPmax)), b_H2HI_LTE = 3.52e-9 * exp(-DMIN(43900./T,EXPmax)), b_H2HI = pow(10., f_v0_LTE*log10(b_H2HI_v0) + f_LTE_v0*log10(b_H2HI_LTE)) * (xH0*nH_cgs) * clumping_factor; // collisional H2-H dissociation; GA08-A1-10
    double b_H2H2_v0 = 5.996e-30 * pow(T,4.1881) * exp(-DMIN(54657.4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881), b_H2H2_LTE = 1.3e-9 * exp(-DMIN(53300./T,EXPmax)), b_H2H2 = pow(10., f_v0_LTE*log10(b_H2H2_v0) + f_LTE_v0*log10(b_H2H2_LTE)) * (xH0*nH_cgs/2.) * clumping_factor; // collisional H2-H2 dissociation; GA08-A1-10
    double b_H2He_v0_log = -27.029 + 3.801*log_T - 29487./T, b_H2He_LTE_log = -2.729 - 1.75*log_T - 23474./T, b_H2He = pow(10., f_v0_LTE*b_H2He_v0_log + f_LTE_v0*b_H2He_LTE_log) * (nHe0*nH_cgs) * clumping_factor; // collisional H2-He dissociation, GA08-A1-11
    double b_H2Hep = (3.7e-14*exp(DMIN(35./T,EXPmax)) + 7.2e-15) * ((nHep+nHepp)*nH_cgs) * clumping_factor; // collisional H2-He+ dissociation, GA08-A1-24+25
    // D questionable - this will really just convert to HD, should exclude here
    double b_H2D; if(T<=2000.) {b_H2D=pow(10.,-56.4737 +5.88886*log_T +7.19692*log_T*log_T +2.25069*log_T*log_T*log_T -2.16903*log_T*log_T*log_T*log_T +0.317887*log_T*log_T*log_T*log_T*log_T);} else {b_H2D=3.17e-10 * exp(-DMIN(5207./T,EXPmax));}
    b_H2D *= (xH0*2.527e-5*nH_cgs) * clumping_factor; // collisional H2-D dissociation, GA08-A1-37, using D abundance from Cooke, Pettini,& Steidel 2018
    double b_H2Dp = 1.e-9 * (DMAX(0.417 + 0.846*log_T - 0.137*log_T*log_T,  0.)) * (nhp*2.527e-5*nH_cgs) * clumping_factor; // collisional H2-D+ dissociation, GA08-A1-39, using D abundance from Cooke, Pettini,& Steidel 2018
    b_H2D*=1.e-10; b_H2Dp*=1.e-10; // see note above on D: include some non-zero here as a 'safety factor', but the overwhelming fraction is going to HD which we account for implicitly in cooling functions so dont explicitly solve here
    double b_H2ext = b_H2Hp + b_H2e + b_H2He + b_H2Hep + b_H2Dp; b_H2HI += b_H2D; b_H2ext*=1./2.; b_H2HI*=1./2.; b_H2H2*=1./4.; // collect dissociation terms where the secondary (e.g. e- does -not- scale with fmol as we define it here, and those where it does to different powers; 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction, becomes 1/4 in H2-H2 equation
    
    double Tdust = 30.; // need to assume something about dust temperature for reaction rates below for dust-phase formation
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate(i, shieldfac, T, pp, cell);
#endif
    double a_Z = 3.e-18*sqrt_T / ((1. +4.e-2*sqrt(T+Tdust) +2.e-3*T +8.e-6*T*T )*(1. +1.e4/exp(DMIN(EXPmax,600./Tdust)))) * f_dustgas_solar * nH0 * clumping_factor; // dust surface formation (assuming dust-to-metals ratio is 0.5*(Z/solar)*dust-to-gas-relative-to-solar in all regions where this is significant), from Glover & Jappsen 2007

    //double a_GP = (1.833e-18 * pow(T,0.88)) / (1. + x_p*1846.*(1.+T/20000.)/sqrt(T)) * xH0 * x_e*nH_cgs * clumping_factor; // gas-phase formation [Glover & Abel 2008, using fitting functions slightly more convenient and assuming H-->H2 much more rapid than other reactions, from Krumholz & McKee 2010; denominator factor accounts for p+H- -> H + H, instead of H2]; note the Nickerson version of this expression omits the ne,-3 term replacing it with ne from Krumholz+McKee which makes it incorrect by a factor of ~1000 in normalization
    double k1,k2,k5,k15,k16,k17,lnTeV=ln_T-9.35915,R51_n=3.62e-17/nH_cgs; // R51_n is contribution from photo-diss assuming ISRF-like since very low-E threshold (=0.755) photons, serves only as a 'floor' here
    if(T<=6000.) {k1=-17.845 + 0.762*log_T + 0.1523*log_T*log_T - 0.03274*log_T*log_T*log_T;} else {k1=-16.420 + 0.1998*log_T*log_T - 5.447e-3*log_T*log_T*log_T*log_T + 4.0415e-5*log_T*log_T*log_T*log_T*log_T*log_T;}
    k1=pow(10.,DMAX(k1,-50.));
    if(T<=300.) {k2=1.5e-9;} else {k2=4.0e-9 * pow(T,-0.17);}
    k5=5.7e-6/sqrt_T + 6.3e-8 - 9.2e-11*sqrt_T + 4.4e-13*T;
    k15=exp(DMAX(-EXPmax, -1.801849334e1 + 2.36085220e0*lnTeV -2.827443e-1*lnTeV*lnTeV +1.62331664e-2*lnTeV*lnTeV*lnTeV
                 -3.36501203e-2*lnTeV*lnTeV*lnTeV*lnTeV +1.17832978e-2*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                 -1.65619470e-3*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV +1.06827520e-4*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                 -2.63128581e-6*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV));
    if(T<1160.45) {k16=1.46629e-16 * pow(T,1.78186);} else {
        k16=exp(DMAX(-EXPmax, -2.0372609e1 + 1.13944933e0*lnTeV
                     -1.4210135e-1*lnTeV*lnTeV +8.4644554e-3*lnTeV*lnTeV*lnTeV
                     -1.4327641e-3*lnTeV*lnTeV*lnTeV*lnTeV +2.0122503e-4*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                     +8.6639632e-5*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV -2.5850097e-5*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                     +2.4555012e-6*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV -8.0683825e-8*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV));}
    if(T<=8000.) {k17=6.9e-9 * pow(T,-0.35);} else {k17=9.6e-7 * pow(T,-0.90);}
    double x_Hminus = k1 * xH0 * x_e / ((k2+k16)*xH0 + (k5+k17)*x_p + k15*x_e + R51_n); // assuming equilibrium H-, using full set of terms from Glover & Jappsen 2007 // use for H- opacity

    double a_GP = k2 * x_Hminus * nH0 * clumping_factor; // actual rate for 2-body gas-phase formation, given x_Hminus calculated above from local equilibrium
    
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * xH0 * clumping_factor_3; // 3-body collisional formation
    double G_LW = 3.3e-11 * urad_G0 * (1./2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0, 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction
    double xi_cr_H2 = (7.525e-16) * (1./2.); // CR dissociation (+ionization), 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction. Using value favored by Indriolo et al for dense GMCs.
//#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON) || defined(RT_ISRF_BACKGROUND) // scale ionization+dissociation rates with local CR energy density
    xi_cr_H2 = Get_CosmicRayIonizationRate_cgs(i, pp, cell) * (1./2.); // scales following Cummings et al. 2016 to 1.6e-17 per eV/cm^3 [this function now defined for everything; differs slightly from default assumed above but thats just being consistent about the baseline normalization here]
//#endif
    // want to solve the implicit equation: f_f = f_0 + g[f_f]*dt, where g[f_f] = df_dt evaluated at f=f_f, so root-find: dt*g[f_f] + f_0-f_f = 0
    // can write this as a quadtratic: x_a*f^2 - x_b_0*f - xb_LW*f + x_c = 0, where xb_LW is a non-linear function of f accounting for the H2 self-shielding terms
    double G_LW_dt_unshielded = G_LW * dtime_cgs; // LW term without shielding, multiplied by timestep for dimensions needed below
    double x_a_00 = (b_3B + b_H2HI - b_H2H2) * dtime_cgs; // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b_00 = (a_GP + a_Z + 2.*b_3B + b_H2HI + b_H2ext + xi_cr_H2) * dtime_cgs; // terms linear in f [note sign, pulling the -sign out here] -- positive-definite, not including '1' for f term
    double x_c_00 = (a_GP + a_Z + b_3B) * dtime_cgs; // terms independent of f -- positive-definite, not including f0 term

    // use the previous-timestep value of fH2 to guess the shielding term and then compute the resulting fH2
    fH2_tmp=fH2_initial; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); // calculate the shielding term
    double x_a=x_a_00,x_c=x_c_00,y_a=0,x_b=0,y_b=0,z_a=0, x_b_0 = x_b_00 + y_ss*G_LW_dt_unshielded;
    double dfH2_linear = fH2_tmp*fH2_tmp*x_a_00 - fH2_tmp*x_b_0 + x_c_00; // linear derivative term for dependence: ffinal = fH2_tmp + dfH2_linear to linear, explicit order
    int change_in_fH2_is_small = 0; // key to know what to do below
    if(fabs(dfH2_linear) < 0.01*fH2_tmp) {change_in_fH2_is_small=1;} // this should be a sufficient criterion below
    
    if(change_in_fH2_is_small==1) // change is sufficiently small, can linearly approximate, essentially with aleading-order expansion of the non-linear solution for small timesteps
    {
        double order_2_corr = 2.*fH2_tmp*x_a_00 - x_b_0;
        if(fabs(order_2_corr) < 1.) {fH2 = fH2_tmp + dfH2_linear * (1. + order_2_corr);} else {fH2 = fH2_tmp + dfH2_linear;}
    } else { // we do a nonlinear solve
        x_b_0=x_b_00 + 1.; x_c=x_c_00 + fH2_initial; // x_c and x_b re-incorporate their constant terms in this limit to make the math easier
        y_a=x_a/(x_c + MIN_REAL_NUMBER); // convenient to convert to dimensionless variable needed for checking definite-ness
        auto molecfrac_root_fn = [&](double x) { // want to find f_mol such that this is 0
            return molecfrac_rootfind_function(x, x00, x01, x_b_0, x_c, y_a, G_LW_dt_unshielded);
        };
        Q_initial = molecfrac_root_fn(fH2_initial);

        x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        double fH2_mid = fH2;
        double Q_mid = molecfrac_root_fn(fH2_mid);

        // OK now let's let the initial and previous-shielding values by the candidate bracket, and if that fails then find another value to bracket the other end
        double root_x_a, root_x_b, root_f_a, root_f_b;
        root_x_a = fH2_mid; root_f_a = Q_mid;
        root_x_b = fH2_initial; root_f_b = Q_initial;
            //if not bracketing we must try other bounds
        if(root_f_b * root_f_a > 0){
            // lower bound
            x_b=x_b_0+G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // recalculate all terms that depend on the shielding
            fH2_min = DMAX(0,DMIN(1,fH2)); // this serves as a lower-limit for fH2
            Q_min = molecfrac_root_fn(fH2_min);
            root_x_b = fH2_min;  root_f_b = Q_min;
            if(Q_min * Q_mid > 0){
            // upper bound
                fH2_tmp=1.; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
                z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
                fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for fH2
                Q_max = molecfrac_root_fn(fH2_max);
                root_x_b = fH2_max;  root_f_b = Q_max;
                if(Q_max * Q_mid > 0){
                    if(Q_min*Q_max > 0){
                        root_x_a = 0; root_x_b = 1; root_f_a = molecfrac_root_fn(0); root_f_b = molecfrac_root_fn(1.);
                    } else {
                        root_x_a = fH2_min; root_x_b = fH2_max; root_f_a = Q_min; root_f_b = Q_max;
                    }
                }
            }
        }

        if(root_f_a * root_f_b < 0){
            // specify desired relative error in fH2 and call the rootfinder
            BrentRootfindResult rfr = brent_rootfind(molecfrac_root_fn, root_x_a, root_x_b, root_f_a, root_f_b, 1e-3, 0., MAXITER);
            fH2 = rfr.x;
            if(rfr.iter > MAXITER) { PRINT_WARNING("WARNING: Particle %lld did not converge to desired H_2 abundance tolerance\n",(long long)(long long)i); }
        } else { // must be at 0 or 1 within machine precision of solution but not bracketing - choose the closer bracketing value of 0 or 1
            if(fabs(root_f_a) < fabs(root_f_b)){fH2 = root_x_a;} else {fH2 = root_x_b;}
        }
    } // end nonlinear solve part

    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    cell[i].MolecularMassFraction_perNeutralH = fH2; // record variable -- this will be used for the next update, meanwhile the total fraction will be used in various routines through the code
    cell[i].MolecularMassFraction = xH0 * cell[i].MolecularMassFraction_perNeutralH; // record variable -- this is largely what is needed below
#endif
}

KOKKOS_FUNCTION double molecfrac_rootfind_function(double fH2, double x00, double x01, double x_b_0, double x_c, double y_a, double G_LW_dt_unshielded){
    const double x_exp_fac=0.00085, w0 = 0.035, EXPmax=90.;
    double x_ss_1=1.+fH2*x01,
	x_ss_sqrt=sqrt(1.+fH2*x00),
	y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)), 
	x_b=x_b_0+y_ss*G_LW_dt_unshielded, 
	y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
    return 1 + y_a*fH2*fH2 - y_b*fH2; // value of the function we are trying to zero, with the updated value of fH2
}



/* get_equilibrium_dust_temperature_estimate body moved to radiation/rt_functions.h
 * as KOKKOS_INLINE_FUNCTION (single source of truth, device-callable from
 * eos_functions.h's return_dust_to_metals_ratio_vs_solar). The non-inline host
 * external symbol is now emitted by radiation/rt_utilities.cc via its existing
 * #undef KOKKOS_INLINE_FUNCTION re-include of rt_functions.h. */



/* Calculates the coefficient for gas-dust collisional heat transfer, such that LambdaDust = gas_dust_heating_coeff * (T-Tdust) in erg cm^3 s^-1. */
KOKKOS_FUNCTION double gas_dust_heating_coeff(int i, double T, double Tdust, struct particle_data *pp, struct gas_cell_data *cell)
{
    double Z_sol=1;
#ifdef METALS
    if(i>=0) {Z_sol = pp[i].Metallicity[0]/All.SolarAbundances[0];}
#endif
    double fdust = return_dust_to_metals_ratio_vs_solar(i,Tdust, pp, cell); // accounting for dust destruction; we avoid calling the function for this because it can create a circular dependency
    return 1.116e-32 * sqrt(T)*(1.-0.8*exp(-75./T)) * Z_sol * fdust;  // Meijerink & Spaans 2005; Hollenbach & McKee 1979,1989. Assumes 10 Angstrom minimum grain size.
}

/* Computes the normalized FUV flux in Habing units G0 

Mode 0: Account only for dust-shielded FUV flux
Mode 1: Use for the C photoionization rate: apply cross- and self-shielding factor from Tielens & Hollenbach 1985
*/
KOKKOS_FUNCTION MyFloat get_FUV_G0(int target, MyFloat shieldfac, int mode, struct particle_data *pp, struct gas_cell_data *cell)
{
    if(target < 0) {return 1;} /* for non-gas elements assume the standard Habing value unless otherwise specified by special module behavior */
    MyFloat G0 = 0.;
#if defined(GALSF_FB_FIRE_RT_LONGRANGE) && !defined(CHIMES)
    G0 += cell[target].Rad_Flux_UV;
    if(CoolTables.gJH0 > 0 && shieldfac > 0) {G0 += sqrt(shieldfac) * (CoolTables.gJH0 / 2.29e-10);} // uvb contribution //
#endif
    double column = evaluate_NH_from_GradRho(pp[target].GradRho, pp[target].KernelRadius, cell[target].Density, pp[target].NumNgb, 1, target, pp) * UNIT_SURFDEN_IN_CGS; // converts to cgs
#ifdef RT_PHOTOELECTRIC
    G0 += cell[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (cell[target].Density * All.cf_a3inv / cell[target].Mass) * UNIT_EGY_DENSITY_IN_HABING; // convert to Habing field //
#endif
#if defined(RT_ISRF_BACKGROUND) && (!defined(RADTRANSFER) || defined(RT_USE_GRAVTREE))                                                                        // latter flag decides whether we do treecol/sobolev here to get the background intensity // add a constant assumed FUV background, for isolated ISM simulations that don't get FUV from local sources self-consistently    
    G0 += All.InterstellarRadiationFieldStrength * 1.7 * exp(-DMAX(pp[target].Metallicity[0] / All.SolarAbundances[0], 1e-4) * column * 500.);                 // RT_ISRF_BACKGROUND rescales the overal ISRF, factor of 1.7 gives Draine 1978 field in Habing units, extinction factor assumes the same FUV band-integrated dust opacity as RT module
#endif
#ifdef SIMPLE_STEADYSTATE_CHEMISTRY
    if(mode == 1){ // Gong 2017 Eq. 9
        MyFloat tau_C = DMIN(column * pp[target].Metallicity[2] / (12 * PROTONMASS_CGS) * 1.6e-17, 100.);
        MyFloat r_H2 = DMIN(2.8e-22 * cell[target].MolecularMassFraction * column * HYDROGEN_MASSFRAC / (2 * PROTONMASS_CGS), 100.);
        G0 *= exp(-tau_C) * exp(-r_H2) / (1 + r_H2);
    }
#endif
    G0 = DMAX(MIN_REAL_NUMBER, DMIN(G0,1e8));
    return G0;
}

/* this function estimates the free electron fraction from heavy ions, assuming a simple mix of cold molecular gas, Mg, and dust, with the ions from singly-ionized Mg, to prevent artificially low free electron fractions */
KOKKOS_FUNCTION double return_electron_fraction_from_heavy_ions(int target, double temperature, double density_cgs, double n_elec_HHe, struct particle_data *pp, struct gas_cell_data *cell)
{
    if(All.ComovingIntegrationOn) {double rhofac=density_cgs/(1000.*COSMIC_BARYON_DENSITY_CGS); if(rhofac<0.2) {return 0;}} // ignore these reactions in the IGM
    double zeta_cr=1.0e-17, f_dustgas=0.01, n_ion_max=4.1533e-5, XH=HYDROGEN_MASSFRAC; // cosmic ray ionization rate (fixed as constant for non-CR runs) and dust-to-gas ratio
    if(target >= 0) {zeta_cr = Get_CosmicRayIonizationRate_cgs(target, pp, cell);} // convert to ionization rate, using models as in Cummings et al. 2016
#ifdef METALS
    if(target>=0) {f_dustgas=0.5*pp[target].Metallicity[0]*return_dust_to_metals_ratio_vs_solar(target,0, pp, cell) + 1.e-15;} // constant dust-to-metals ratio [floor purely numerical here, needed to avoid a spurious divergence below]
#ifdef COOL_METAL_LINES_BY_SPECIES
    if(target>=0) {n_ion_max = (All.SolarAbundances[6]/24.3)/XH;} // limit, to avoid over-ionization at low metallicities
#endif
#endif
    /* Regime I: highly/photo-ionized, any contributions here would be negligible -- no need to continue */
    if(n_elec_HHe > 0.01) {return n_ion_max;} // contribute something negligible, doesn't matter here //
    double a_grain_micron=0.1, m_ion=24.305*PROTONMASS_CGS, mu_eff=2.38, m_neutrals=mu_eff*PROTONMASS_CGS, m_grain=4.189e-12*(2.4)*a_grain_micron*a_grain_micron*a_grain_micron, ngrain_ngas=(m_neutrals/m_grain)*f_dustgas; // effective size of grains that matter at these densities, and ions [here Mg] that dominate
    double k_ei=9.77e-8, y0=sqrt(m_ion/ELECTRONMASS_CGS);
    double y = exp(1.)*y0, ln_oneplusy=log(1.+y), psi_0 = 1.-ln_oneplusy + ln_oneplusy/(1.+ln_oneplusy) * log(ln_oneplusy*(1.+1./y)); // changed convention: using second-order expansion to solve for psi, psi should be < 0, and electron absorption should be suppressed for larger psi because of grain charge leading to coulomb repulsion
    double k_eg_00=0.0195*a_grain_micron*a_grain_micron*sqrt(temperature), k_eg_0=k_eg_00*exp(psi_0);
    double n_crit=k_ei*zeta_cr/(k_eg_0*ngrain_ngas*k_eg_0*ngrain_ngas), n_eff=density_cgs/m_neutrals;

    /* Regime II: CR-ionized with high enough ionization fraction that gas-phase recombinations dominate */
    if(n_eff < 0.01*n_crit) {return DMIN(n_ion_max , sqrt(zeta_cr / (k_ei * n_eff)) / (XH*mu_eff) );} // CR-dominated off gas with recombination via ions -- put into appropriate units here, but basically just what we would expect here. mu_eff factor transforms from n_e/n_neutrals to n_e/n_H because neutrals if molecular not as numerous
    if(n_eff < 100.*n_crit) {return DMIN(n_ion_max , 0.5 * (sqrt(4.*n_crit/n_eff + 1.) - 1.) * (k_eg_0*ngrain_ngas)/(k_ei*XH*mu_eff));} // interpolates between gas-recombination and dust-dominated regimes (just a simple interpolation function)
    
    /* Regime III: recombination dominated by dust, but dust has a 'normal' efficiency of absorbing grains, and most of the charge is still free+ion (n_ions ~ n_free_electrons >> Z_grains*n_grains, even if -rate- of grain absorption of free e- dominates over k_ei) */
    double psi_fac=16.71/(a_grain_micron*temperature), alpha=zeta_cr*psi_fac/(k_eg_00 * ngrain_ngas*ngrain_ngas * n_eff), alpha_min=0.02, alpha_max=10.; /* Z*psi_fac = Z*e^2/(a_grain*kB*T) to dimensionless [psi] units; alpha=prefactor in charge equation: psi = alpha*(exp[-psi] - y/(1-psi)); this alpha factor gives n_e/n_charge_summed_grains = (this)/(exp[psi]*psi) -- when this gets to smaller than ~10, can start to deplete all the charge onto grains, need to move to the next regime [note here psi only puts you more-safely into free e- regime, so dont need solution below if psi very small, hence ignore it here to be conservative */
    if(alpha > alpha_max) {return DMIN(n_ion_max , zeta_cr / (k_eg_0*ngrain_ngas * XH*mu_eff*n_eff ) );} // regime III limit above: here XH*mu_eff factor again just converts to the correct units. here balancing ionization [zeta_cr*n_neutrals] with absorption by dust grains [k_eg_0*ngrain_ngas*n_neutrals*x_e*n_neutrals]
    
    /* Regime IV: recombination dominated by dust and grains dominate the charge, strongly suppressing the free charges (n_ion ~ Z_grains*n_grains >> n_free_electrons) */
    if(alpha < 1.e-4) {return DMIN(n_ion_max , zeta_cr / (k_eg_00*ngrain_ngas * XH*mu_eff*n_eff) );} // alpha very small means that there is a very tiny number of free e- per dust grain to potentially absorb. this means the mean grain charge must be << 1, so psi->small, negligible correction here, just same expression as above with k_eg_0 (no psi correction, or k_eg calculated without any charge correction here. note in this limit, the ion fraction is larger than the free electron fraction by a ratio sqrt[m_ion/m_e] ~ k_eg/k_ig [see Keith+Wardle 2014, section 3.2, which adopts this limit
    if(alpha < alpha_min) {double psi=0.5*(1.-sqrt(1.+4.*y0*alpha)); return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff) );} // small-psi-limit. when alpha is sufficiently small, here in limit where [already above] can safely neglect k_ei, so solution for x_e is same as above/below, but with a different psi, here need to solve psi = exp[-psi]*alpha - y0*alpha/(1-psi). if alpha is small, psi should also be small, can drop the exp[-psi]*alpha term and this becomes a reasonable approximation here.
    double psi_xmin=0.5*(1.-sqrt(1.+4.*y0*alpha_min)), psi=psi_0 + (psi_xmin-psi_0)*2./(1.+alpha/alpha_min); // this interpolates between the asymptotic limmits at low and high alpha, where we can obtain high-accuracy solutions here. this is a completely ad-hoc fitting function for the numerical solutions in this intermediate range, which only relates to a rather modest range of alpha and is ensured to interpolate correctly to the extremes we already pulled out above.
    return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff));
}




/* this function evaluates Compton heating+cooling rates and synchrotron cooling for thermal gas populations, accounting for the
    explicitly-evolved radiation field if it is evolved (otherwise assuming a standard background), and B-fields if they
    are evolved, as well as the proper relativistic or non-relativistic effects and two-temperature plasma effects. */
KOKKOS_FUNCTION double evaluate_Compton_heating_cooling_rate(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg, struct gas_cell_data *cell)
{
    double Lambda = 0;
    double compton_prefac_eV = 2.16e-35 / nHcgs; // multiply field in eV/cm^3 by this and temperature difference to obtain rate

    double e_CMB_eV=0.262*All.cf_a3inv/All.cf_atime, T_cmb = 2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    Lambda += compton_prefac_eV * n_elec * e_CMB_eV * (T-T_cmb);

    double e_UVB_eV = shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    Lambda += compton_prefac_eV * n_elec * e_UVB_eV * (T-2.e4); // assume very crude approx Compton temp ~2e4 for UVB
    
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    if(target >= 0)
    {
        int k; double E_tot_to_evol_eVcgs = (cell[target].Density*All.cf_a3inv/cell[target].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            double e_tmp = cell[target].Rad_E_gamma_Pred[k] * E_tot_to_evol_eVcgs, Teff = 0;
            
#if defined(GALSF_FB_FIRE_RT_LONGRANGE) /* three-band (UV, OPTICAL, IR) approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
            if(k==RT_FREQ_BIN_FIRE_IR) {Teff=30.;}
            if(k==RT_FREQ_BIN_FIRE_OPT) {Teff=4000.;}
            if(k==RT_FREQ_BIN_FIRE_UV) {Teff=15000.;}
#endif
#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
            if(k==RT_FREQ_BIN_INFRARED) {Teff=cell[target].Radiation_Temperature;}
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
            if(k==RT_FREQ_BIN_OPTICAL_NIR) {Teff=2800.;}
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
            if(k==RT_FREQ_BIN_NUV) {Teff=12000.;}
#endif
#if defined(RT_PHOTOELECTRIC) /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_PHOTOELECTRIC) {Teff=24400.;}
#endif
#if defined(RT_LYMAN_WERNER) /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_LYMAN_WERNER) {Teff=28800.;}
#endif
#if defined(RT_CHEM_PHOTOION) /* Hydrogen and Helium ionizing bands: H0 here */
            if(k==RT_FREQ_BIN_H0) {Teff=2340.*All.rt_nu_eff_eV[k];}
#endif
#if defined(RT_PHOTOION_MULTIFREQUENCY) /* Hydrogen and Helium ionizing bands: He bands */
            if(k==RT_FREQ_BIN_He0 || k==RT_FREQ_BIN_He1 || k==RT_FREQ_BIN_He2) {Teff=2340.*All.rt_nu_eff_eV[k];}
#endif
#if defined(RT_SOFT_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_SOFT_XRAY) {Teff=3.6e6;}
#endif
#if defined(RT_HARD_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_HARD_XRAY) {Teff=1.7e7;}
#endif
            if(Teff < 3.e4) {e_tmp *= n_elec;} // low-energy radiation acts inefficiently on neutrals here
            Lambda += compton_prefac_eV * e_tmp * (T - Teff); // add to compton heating/cooling terms
        }
    }
#else // no explicit RHD terms evolved, so assume a MW-like ISRF instead
    double e_IR_eV=0.31, T_IR=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_OUV_eV=0.66, T_OUV=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
    Lambda += compton_prefac_eV * n_elec * (e_IR_eV*(T-T_IR) + e_OUV_eV*(T-T_OUV));
#endif

#ifdef SINK_COMPTON_HEATING /* custom band to represent (non)relativistic X-ray compton cooling from an AGN source without full RHD */
    if(target >= 0)
    {
        double e_agn = (cell[target].Rad_Flux_AGN * UNIT_FLUX_IN_CGS) / (C_LIGHT_CGS * ELECTRONVOLT_IN_ERGS), T_agn=2.e7; /* approximate from Sazonov et al. */
        Lambda += compton_prefac_eV * e_agn * (T-T_agn); // since the heating here is primarily hard X-rays, and the cooling only relevant for very high temps, do not have an n_elec here
    }
#endif

#ifdef MAGNETIC /* include sychrotron losses as well as long as we're here, since these scale more or less identically just using the magnetic instead of radiation energy */
    if(target >= 0)
    {
        double b_muG = cell[target].Bfield_microGauss(), U_mag_ev=0.0248342*b_muG*b_muG, T_rad_background_at_emission = get_background_radiation_temperature_for_emission_corrections(target, cell);
        Lambda += compton_prefac_eV * U_mag_ev * (T-T_rad_background_at_emission); // synchrotron losses proportional to temperature (non-relativistic here), as inverse compton, just here without needing to worry about "T-T_eff", as if T_eff->0
    }
#endif

    double T_eff_for_relativistic_corr = T; /* used below, but can be corrected */
    if(Lambda > 0) /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
    {
        //double Lambda_limiter_var = 1.483e34 * Lambda*Lambda*T; /* = (Lambda/2.6e-22)^2 * (T/1e9): if this >> 1, follow CAFG and cap at cooling rate assuming equilibrium e- temp from Coulomb exchange balancing compton */
        double Lambda_limiter_var = 2.38e33 * Lambda*Lambda*T; /* = (slightly revised with coefficients recalculated to be more general) */
        if(Lambda_limiter_var > 1 && T < 1.e12) { /* at higher temps, protons start to become relativistic */
            Lambda_limiter_var = 1./pow(Lambda_limiter_var, 1./5.); if(T*Lambda_limiter_var > 1.5e9) {Lambda_limiter_var = 1./pow(Lambda_limiter_var, 1./7.);} /* correction factor is different for the relativistic case */
            if(Lambda_limiter_var < 0.01) {Lambda_limiter_var = 0.01;} /* formulae above not valid for too large a ratio here, because then need to include additional proton terms */
            Lambda*=Lambda_limiter_var; T_eff_for_relativistic_corr*=Lambda_limiter_var;}
    }
    if(T_eff_for_relativistic_corr > 3.e7) {double tau=T_eff_for_relativistic_corr/1.5e9; Lambda*=tau; if(tau < 20.) {Lambda/=1.-exp(-tau);}} /* relativistic correction term, becomes important at >1e9 K, enhancing rate */
    
    return Lambda;
}


/* this function defines an effective background radiation temperature for purposes of computing the emission corrections above */
KOKKOS_FUNCTION double get_background_radiation_temperature_for_emission_corrections(int target, struct gas_cell_data *cell)
{
    double T_cmb = 2.73/All.cf_atime;
#ifdef RT_ISRF_BACKGROUND
    if(!All.ComovingIntegrationOn) {T_cmb *= 1.+All.RadiationBackgroundRedshift;}
#endif
#if defined(RT_INFRARED)
    if(target >= 0)
    {
        double e_tmp_CMB = 0.262*All.cf_a3inv/All.cf_atime;
#ifdef RT_ISRF_BACKGROUND
        if(!All.ComovingIntegrationOn) {e_tmp_CMB *= pow(1.+All.RadiationBackgroundRedshift,4);}
#endif
        double e_tot_to_evol_eVcgs = (cell[target].Density*All.cf_a3inv/cell[target].Mass) * UNIT_PRESSURE_IN_EV;
        double e_tmp_IR = cell[target].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED] * e_tot_to_evol_eVcgs, T_tmp_IR = cell[target].Radiation_Temperature;
        return (e_tmp_IR * T_tmp_IR + e_tmp_CMB * T_cmb) / (e_tmp_IR + e_tmp_CMB); // if evolving IR band, use sum of it plus cmb for background
    }
#endif
    return T_cmb; // default to assume total radiation temperature is cmb-dominated
}



/* ThermalProperties, return_local_gammamultiplier, return_uvb_shieldfac:
   definitions moved to cooling_functions.h for cross-TU GPU callability. */


#ifdef CHIMES
/* This routine updates the ChimesGasVars structure for particle target. */
void chimes_update_gas_vars(int target, struct particle_data *pp, struct gas_cell_data *cell)
{
  double dt = get_particle_timestep_in_physical(target, pp);
  double u_old_cgs = DMAX(All.MinEgySpec, cell[target].InternalEnergy) * UNIT_SPECEGY_IN_CGS;
  double rho_cgs = cell[target].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

#ifdef COOL_METAL_LINES_BY_SPECIES
  double H_mass_fraction = 1.0 - (pp[target].Metallicity[0] + pp[target].Metallicity[1]);
#else
  double H_mass_fraction = XH;
#endif

  ChimesGasVars[target].temperature = (ChimesFloat) chimes_convert_u_to_temp(u_old_cgs, rho_cgs, target, pp, cell);
  ChimesGasVars[target].nH_tot = (ChimesFloat) (H_mass_fraction * rho_cgs / PROTONMASS_CGS);
  ChimesGasVars[target].ThermEvolOn = All.ChimesThermEvolOn;

  // If there is an EoS, need to set TempFloor to that instead.
  ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);
#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
    if (cell[target].DelayTimeHII > 0) {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(HIIRegion_Temp, 10.1);}
        else {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);}
#endif

  // Flag to control how the temperature
  // floor is implemented in CHIMES.
  ChimesGasVars[target].temp_floor_mode = 1;

  // Extragalactic UV background
  ChimesGasVars[target].isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0];
  ChimesGasVars[target].isotropic_photon_density[0] *= chimes_rad_field_norm_factor;

  ChimesGasVars[target].G0_parameter[0] = chimes_table_spectra.G0_parameter[0];
  ChimesGasVars[target].H2_dissocJ[0] =  chimes_table_spectra.H2_dissocJ[0];

#ifdef CHIMES_STELLAR_FLUXES
  int kc;
  for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++)
    {
      ChimesGasVars[target].isotropic_photon_density[kc + 1] = (ChimesFloat) (cell[target].Chimes_fluxPhotIon[kc] / C_LIGHT_CGS);

#ifdef CHIMES_HII_REGIONS
    if(cell[target].DelayTimeHII > 0) {
        ChimesGasVars[target].isotropic_photon_density[kc + 1] += (ChimesFloat) (cell[target].Chimes_fluxPhotIon_HII[kc] / C_LIGHT_CGS);
        ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) ((cell[target].Chimes_G0[kc] + cell[target].Chimes_G0_HII[kc]) / DMAX((cell[target].Chimes_fluxPhotIon[kc] + cell[target].Chimes_fluxPhotIon_HII[kc]), 1.0e-300));
	} else {ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (cell[target].Chimes_G0[kc] / DMAX(cell[target].Chimes_fluxPhotIon[kc], 1.0e-300));}
    ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#else
      ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (cell[target].Chimes_G0[kc] / DMAX(cell[target].Chimes_fluxPhotIon[kc], 1.0e-300));
      ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#endif
    }
#endif

  ChimesGasVars[target].cr_rate = (ChimesFloat) cr_rate;  // For now, assume a constant cr_rate.
  ChimesGasVars[target].hydro_timestep = (ChimesFloat) (dt * UNIT_TIME_IN_CGS);

  ChimesGasVars[target].ForceEqOn = ChimesEqmMode;
  ChimesGasVars[target].divVel = (ChimesFloat) pp[target].Particle_DivVel / UNIT_TIME_IN_CGS;
  if (All.ComovingIntegrationOn)
    {
      ChimesGasVars[target].divVel *= (ChimesFloat) All.cf_a2inv;
      ChimesGasVars[target].divVel += (ChimesFloat) (3 * All.cf_hubble_a / UNIT_TIME_IN_CGS);  /* Term due to Hubble expansion */
    }
  ChimesGasVars[target].divVel = (ChimesFloat) fabs(ChimesGasVars[target].divVel);

#ifndef COOLING_OPERATOR_SPLIT
  if(cell[target].CoolingIsOperatorSplitThisTimestep==0) {ChimesGasVars[target].constant_heating_rate = ChimesGasVars[target].nH_tot * ((ChimesFloat) cell[target].DtInternalEnergy);}
#else
  ChimesGasVars[target].constant_heating_rate = 0.0;
#endif

#ifdef CHIMES_SOBOLEV_SHIELDING
  double surface_density;
  surface_density = evaluate_NH_from_GradRho(pp[target].GradRho,pp[target].KernelRadius,cell[target].Density,pp[target].NumNgb,1,target,pp) * UNIT_SURFDEN_IN_CGS; // converts to cgs
  ChimesGasVars[target].cell_size = (ChimesFloat) (shielding_length_factor * surface_density / rho_cgs);
#else
  ChimesGasVars[target].cell_size = 1.0;
#endif

  ChimesGasVars[target].doppler_broad = 7.1;  // km/s. For now, just set this constant. Thermal broadening is also added within CHIMES.

  ChimesGasVars[target].InitIonState = ChimesInitIonState;

#ifdef CHIMES_HII_REGIONS
  if (cell[target].DelayTimeHII > 0.0) {ChimesGasVars[target].cell_size = 1.0;} // switch of shielding in HII regions
#endif

#ifdef CHIMES_TURB_DIFF_IONS
  chimes_update_turbulent_abundances(target, 0, pp, cell);
#endif

#if defined(COOL_METAL_LINES_BY_SPECIES) && !defined(GALSF_FB_NOENRICHMENT)
  chimes_update_element_abundances(target, pp, cell);
#endif

  return;
}

#ifdef COOL_METAL_LINES_BY_SPECIES
/* This routine re-computes the element abundances from
 * the metallicity array and updates the individual ion
 * abundances accordingly. */
void chimes_update_element_abundances(int i, struct particle_data *pp, struct gas_cell_data *cell)
{
  double H_mass_fraction = 1.0 - (pp[i].Metallicity[0] + pp[i].Metallicity[1]);

  /* Update the element abundances in ChimesGasVars. */
  ChimesGasVars[i].element_abundances[0] = (ChimesFloat) (pp[i].Metallicity[1] / (4.0 * H_mass_fraction));   // He
  ChimesGasVars[i].element_abundances[1] = (ChimesFloat) (pp[i].Metallicity[2] / (12.0 * H_mass_fraction));  // C
  ChimesGasVars[i].element_abundances[2] = (ChimesFloat) (pp[i].Metallicity[3] / (14.0 * H_mass_fraction));  // N
  ChimesGasVars[i].element_abundances[3] = (ChimesFloat) (pp[i].Metallicity[4] / (16.0 * H_mass_fraction));  // O
  ChimesGasVars[i].element_abundances[4] = (ChimesFloat) (pp[i].Metallicity[5] / (20.0 * H_mass_fraction));  // Ne
  ChimesGasVars[i].element_abundances[5] = (ChimesFloat) (pp[i].Metallicity[6] / (24.0 * H_mass_fraction));  // Mg
  ChimesGasVars[i].element_abundances[6] = (ChimesFloat) (pp[i].Metallicity[7] / (28.0 * H_mass_fraction));  // Si
  ChimesGasVars[i].element_abundances[7] = (ChimesFloat) (pp[i].Metallicity[8] / (32.0 * H_mass_fraction));  // S
  ChimesGasVars[i].element_abundances[8] = (ChimesFloat) (pp[i].Metallicity[9] / (40.0 * H_mass_fraction));  // Ca
  ChimesGasVars[i].element_abundances[9] = (ChimesFloat) (pp[i].Metallicity[10] / (56.0 * H_mass_fraction)); // Fe

  ChimesGasVars[i].metallicity = (ChimesFloat) (pp[i].Metallicity[0] / 0.0129);  // In Zsol. CHIMES uses Zsol = 0.0129.
  ChimesGasVars[i].dust_ratio = ChimesGasVars[i].metallicity;

#ifdef CHIMES_METAL_DEPLETION
#ifdef _OPENMP
  int ThisThread = omp_get_thread_num();
#else
  int ThisThread = 0;
#endif
  chimes_compute_depletions(ChimesGasVars[i].nH_tot, ChimesGasVars[i].temperature, ThisThread);

  ChimesGasVars[i].element_abundances[1] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[0]; // C
  ChimesGasVars[i].element_abundances[2] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[1]; // N
  ChimesGasVars[i].element_abundances[3] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[2]; // O
  ChimesGasVars[i].element_abundances[5] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[3]; // Mg
  ChimesGasVars[i].element_abundances[6] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[4]; // Si
  ChimesGasVars[i].element_abundances[7] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[5]; // S
  ChimesGasVars[i].element_abundances[9] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[6]; // Fe

  ChimesGasVars[i].dust_ratio *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDustRatio;
#endif // CHIMES_METAL_DEPLETION

  /* The element abundances may have changed, so use the check_constraint_equations()
   * routine to update the abundance arrays accordingly. If metals have been injected
   * into a particle, it is distributed across all of the metal's atomic/ionic/molecular
   * species, preserving the ion and molecule fractions of that element. */

  check_constraint_equations(&(ChimesGasVars[i]), &ChimesGlobalVars);
}
#endif // COOL_METAL_LINES_BY_SPECIES

#ifdef CHIMES_TURB_DIFF_IONS
/* mode == 0: re-compute the CHIMES abundance array from the ChimesNIons and ChimesHtot
 *            arrays that are used to track turbulent diffusion of ions and molecules.
 * mode == 1: update the ChimeSNIons array from the current CHIMES abundance array.
 */
void chimes_update_turbulent_abundances(int i, int mode, struct particle_data *pp, struct gas_cell_data *cell)
{
  int k_species;
  double NHtot = (1.0 - (pp[i].Metallicity[0] + pp[i].Metallicity[1])) * (cell[i].Mass * UNIT_MASS_IN_CGS) / PROTONMASS_CGS;

  if (mode == 0)
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	ChimesGasVars[i].abundances[k_species] = (ChimesFloat) (cell[i].ChimesNIons[k_species] / NHtot);
    }
  else
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	cell[i].ChimesNIons[k_species] = ((double) ChimesGasVars[i].abundances[k_species]) * NHtot;
    }
}
#endif // CHIMES_TURB_DIFF_IONS

#ifdef CHIMES_METAL_DEPLETION
void chimes_init_depletion_data(void)
{
  int i;

#ifdef _OPENMP
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(maxThreads * sizeof(struct Chimes_depletion_data_structure));
#else
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(sizeof(struct Chimes_depletion_data_structure));
#endif

  // Elements in Jenkins (2009) in the order
  // C, N, O, Mg, Si, P, S, Cl, Ti, Cr, Mn, Fe,
  // Ni, Cu, Zn, Ge, Kr
  // Solar abundances are as mass fractions,
  // taken from the Cloudy default values, as
  // used in CHIMES.
  double SolarAbund[DEPL_N_ELEM] = {2.07e-3, 8.36e-4, 5.49e-3, 5.91e-4, 6.83e-4, 7.01e-6, 4.09e-4, 4.72e-6, 3.56e-6, 1.72e-5, 1.12e-5, 1.1e-3, 7.42e-5, 7.32e-7, 1.83e-6, 2.58e-7, 1.36e-7};

  // Fit parameters, using equation 10 of Jenkins (2009).
  // Where possible, we take the updated fit parameters
  // A2 and B2 from De Cia et al. (2016), which we convert
  // to the Ax and Bx of J09 (with zx = 0). Otherwise, we
  // use the original J09 parameters. We list these in the
  // order Ax, Bx, zx.
  double DeplPars[DEPL_N_ELEM][3] = {{-0.101, -0.193, 0.803}, // C
				     {0.0, -0.109, 0.55},     // N
				     {-0.101, -0.172, 0.0},     // O
				     {-0.412, -0.648, 0.0},     // Mg
				     {-0.426, -0.669, 0.0},     // Si
				     {-0.068, -0.091, 0.0},      // P
				     {-0.189, -0.324, 0.0},     // S
				     {-1.242, -0.314, 0.609}, // Cl
				     {-2.048, -1.957, 0.43},  // Ti
				     {-0.892, -1.188, 0.0},      // Cr
				     {-0.642, -0.923, 0.0},      // Mn
				     {-0.851, -1.287, 0.0},     // Fe
				     {-1.490, -1.829, 0.599}, // Ni
				     {-0.710, -1.102, 0.711}, // Cu
				     {-0.182, -0.274, 0.0},       // Zn
				     {-0.615, -0.725, 0.69},  // Ge
				     {-0.166, -0.332, 0.684}}; // Kr

  int n_thread;
#ifdef _OPENMP
  n_thread = maxThreads;
#else
  n_thread = 1;
#endif

  for (i = 0; i < n_thread; i++)
    {
      memcpy(ChimesDepletionData[i].SolarAbund, SolarAbund, DEPL_N_ELEM * sizeof(double));
      memcpy(ChimesDepletionData[i].DeplPars, DeplPars, DEPL_N_ELEM * 3 * sizeof(double));
      // DustToGasSaturated is the dust to gas ratio when F_star = 1.0, i.e. at maximum depletion onto grains.
      ChimesDepletionData[i].DustToGasSaturated = 5.9688e-03;
    }
}

/* Computes the linear fits as in Jenkins (2009)
 * or De Cia et al. (2016). Note that this returns
 * log10(fraction in the gas phase) */
double chimes_depletion_linear_fit(double nH, double T, double Ax, double Bx, double zx)
{
  // First, compute the parameter F_star, using the
  // best-fit relation from Fig. 16 of Jenkins (2009).
  double F_star = 0.772 + (0.461 * log10(nH));
  double depletion;

  // Limit F_star to be no greater than unity
  if (F_star > 1.0)
    F_star = 1.0;

  // All metals are in the gas phase if T > 10^6 K
  if (T > 1.0e6)
    return 0.0;
  else
    {
      depletion = Bx + (Ax * (F_star - zx));

      // Limit depletion to no greater than 0.0 (remember: it is a log)
      if (depletion > 0.0)
	return 0.0;
      else
	return depletion;
    }
}

void chimes_compute_depletions(double nH, double T, int thread_id)
{
  int i;
  double pars[DEPL_N_ELEM][3];
  memcpy(pars, ChimesDepletionData[thread_id].DeplPars, DEPL_N_ELEM * 3 * sizeof(double));

  // ChimesDepletionFactors array is for the metals
  // in the order C, N, O, Mg, Si, S, Fe. The other
  // metals in CHIMES are not depleted.
  ChimesDepletionData[thread_id].ChimesDepletionFactors[0] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[0][0], pars[0][1], pars[0][2])); // C
  ChimesDepletionData[thread_id].ChimesDepletionFactors[1] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[1][0], pars[1][1], pars[1][2])); // N
  ChimesDepletionData[thread_id].ChimesDepletionFactors[2] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[2][0], pars[2][1], pars[2][2])); // O
  ChimesDepletionData[thread_id].ChimesDepletionFactors[3] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[3][0], pars[3][1], pars[3][2])); // Mg
  ChimesDepletionData[thread_id].ChimesDepletionFactors[4] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[4][0], pars[4][1], pars[4][2])); // Si
  ChimesDepletionData[thread_id].ChimesDepletionFactors[5] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[6][0], pars[6][1], pars[6][2])); // S
  ChimesDepletionData[thread_id].ChimesDepletionFactors[6] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[11][0], pars[11][1], pars[11][2])); // Fe

  // The dust abundance as used in CHIMES will be the
  // metallicity in solar units multiplied by ChimesDustRatio
  ChimesDepletionData[thread_id].ChimesDustRatio = 0.0;
  for (i = 0; i < DEPL_N_ELEM; i++)
    ChimesDepletionData[thread_id].ChimesDustRatio += ChimesDepletionData[thread_id].SolarAbund[i] * (1.0 - pow(10.0, chimes_depletion_linear_fit(nH, T, pars[i][0], pars[i][1], pars[i][2])));

  // The above sum gives the dust to gas mass ratio.
  // Now normalise it by the saturated value.
  ChimesDepletionData[thread_id].ChimesDustRatio /= ChimesDepletionData[thread_id].DustToGasSaturated;
}
#endif // CHIMES_METAL_DEPLETION


void chimes_gizmo_exit(void)
{
  endrun(56275362);
}
#endif // CHIMES

#endif /* COOLING */

/* Kokkos lifecycle and sync functions — OUTSIDE #ifdef COOLING because they
   must exist for any GPU build (e.g. NUCLEAR_NETWORK without COOLING). */
void gizmo_kokkos_initialize(int argc, char *argv[]) {
    Kokkos::initialize(argc, argv);
#if defined(__CUDACC__)
    {size_t stack_size = 65536; /* 64KB per thread — cooling kernel is deeply nested:
         DoCooling → CoolingRateFromU → convert_u_to_temp → find_abundances_and_rates
         plus metal cooling table lookups. 32KB may overflow on some GPU architectures. */
     cudaError_t _e = cudaDeviceSetLimit(cudaLimitStackSize, stack_size);
     if(_e != cudaSuccess) {printf("[GPU] WARNING: cudaDeviceSetLimit(stack,%zu) failed: %s\n", stack_size, cudaGetErrorString(_e)); fflush(stdout);}
     else {printf("[GPU] CUDA thread stack size set to %zu bytes\n", stack_size); fflush(stdout);}}
#elif defined(__HIPCC__)
    {size_t stack_size = 65536;
     hipError_t _e = hipDeviceSetLimit(hipLimitStackSize, stack_size);
     if(_e != hipSuccess) {printf("[GPU] WARNING: hipDeviceSetLimit(stack,%zu) failed: %s\n", stack_size, hipGetErrorString(_e)); fflush(stdout);}
     else {printf("[GPU] HIP thread stack size set to %zu bytes\n", stack_size); fflush(stdout);}}
#endif
}
void gizmo_kokkos_finalize(void)
{
    /* Persistent device-side state is handed back before Kokkos is torn down;
       finalizing with allocations still tracked is what this ordering avoids.
       This is the single chokepoint every shutdown path already routes through. */
    gx_touched_set_release();
    gx_motion_target_release();
    gpu_node_dirty_release();
    Kokkos::finalize();
}
/* Best-effort drain of in-flight device work. Used by the reviewed hard-abort
 * path (core/run.cc) to quiesce the GPU before MPI_Abort, reducing (not
 * eliminating) the Vista CG-stuck wedge risk. Guarded so it is safe to call
 * before Kokkos init or after finalize. */
void gizmo_kokkos_fence(void) { if(Kokkos::is_initialized() && !Kokkos::is_finalized()) { Kokkos::fence(); } }

/* Central registry of per-TU AllDeviceMirror addresses.
 * Each GPU TU's gpu_all_mirror.h include instantiates a private
 * `static __managed__ AllDeviceMirror` plus auto-registration via a
 * __attribute__((constructor)) function AND via the dispatch-boundary
 * GIZMO_GPU_ENSURE_ALL_FRESH() macro. Either path lands here. The
 * registry de-duplicates so double registration is harmless.
 * Meyer's-singleton vector keeps static-init-order across TUs safe. */
#if defined(GIZMO_GPU_COMPILER)
static std::vector<struct global_data_all_processes *> &gizmo_all_device_mirror_registry(void)
{
    static std::vector<struct global_data_all_processes *> reg;
    return reg;
}
extern "C" void gizmo_register_all_device_mirror(struct global_data_all_processes *m)
{
    if(!m) return;  /* HIP: ignore a pre-managed-init NULL &AllDeviceMirror capture */
    auto &reg = gizmo_all_device_mirror_registry();
    for(auto *x : reg) { if(x == m) return; }
    reg.push_back(m);
}
#endif

/* Central host->device sync. Fences before and after the per-TU mirror
   copies: leading fence drains any in-flight kernels reading mirrors,
   trailing fence publishes new values to subsequent launches. On
   host-only Kokkos backends there is no mirror; the function reduces
   to the fences. */
void gizmo_gpu_sync_all(void) {
    Kokkos::fence();
#if defined(GIZMO_GPU_COMPILER)
    const struct global_data_all_processes *host_all = gizmo_host_all_ptr();
    auto &reg = gizmo_all_device_mirror_registry();
    for(auto *m : reg) { if(m) *m = *host_all; }  /* skip any residual NULL (HIP pre-init capture) */
#endif
    Kokkos::fence();
}

#ifndef COOLING
/* Builds without the tables still need this symbol, so that code taking a view
   reads the same in every configuration and no call site has to ask which one
   it is in. The view is empty; the physics that would consult it is compiled
   out alongside the tables. A CHIMES build is NOT one of these: it still
   compiles the definition above, whose body is what the CHIMES gate empties. */
struct PhysicsTablesView gizmo_physics_tables_view(void)
{
    struct PhysicsTablesView view = {};
    return view;
}
#endif
