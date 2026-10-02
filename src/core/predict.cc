#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "../declarations/allvars.h"
#include "../declarations/lifecycle_counters.h"
#include "../core/proto.h"
#include "../system/gpu_particles_arena.h"
#include "../mesh/gpu_neighbor_list.h" /* gizmo_mark_kernel_radius_dirty_* */
#include "../mesh/kernel.h"
#include "../gravity/gravtree_force_kernel.h" /* shared weight-function SSOT (grav_weight_function_for_weighted_motion_smoothing) */
#include "../gravity/binary_functions.h" /* odeint_super_timestep (SINGLE_STAR_TIMESTEPPING) */

/*! Routines for the drift/predict step */

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially in detail (although the highest-level algorithm
 * structure remains essentially the same)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO, and many new
 * options and subroutines added for flexibility with different
 * hydro solvers, timestepping schemes, boundary conditions, and
 * different mesh-motion options added.
 */

void reconstruct_timebins(void)
{
    int i, bin;
    long long glob_sum;
    
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        TimeBinCount[bin] = 0;
        TimeBinCountGas[bin] = 0;
        FirstInTimeBin[bin] = -1;
        LastInTimeBin[bin] = -1;
#ifdef GALSF
        TimeBinSfr[bin] = 0;
#endif
#ifdef SINK_PARTICLES
        TimeBin_Sink_mass[bin] = 0;
        TimeBin_Sink_dynamicalmass[bin] = 0;
        TimeBin_Sink_Mdot[bin] = 0;
        TimeBin_Sink_Medd[bin] = 0;
#endif
    }
    
    for(i = 0; i < NumPart; i++)
    {
        bin = P[i].TimeBin;
        
        if(TimeBinCount[bin] > 0)
        {
            PrevInTimeBin[i] = LastInTimeBin[bin];
            NextInTimeBin[i] = -1;
            NextInTimeBin[LastInTimeBin[bin]] = i;
            LastInTimeBin[bin] = i;
        }
        else
        {
            FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
            PrevInTimeBin[i] = NextInTimeBin[i] = -1;
        }
        TimeBinCount[bin]++;
        if(P[i].Type == 0)
            TimeBinCountGas[bin]++;
        
#ifdef GALSF
        if(P[i].Type == 0)
            TimeBinSfr[bin] += CellP[i].Sfr;
#endif
#ifdef SINK_PARTICLES
        if(P[i].Type == 5)
        {
            TimeBin_Sink_mass[bin] += P[i].Sink_Mass;
            TimeBin_Sink_dynamicalmass[bin] += P[i].Mass;
            TimeBin_Sink_Mdot[bin] += P[i].Sink_Mdot;
            TimeBin_Sink_Medd[bin] += P[i].Sink_Mdot / P[i].Sink_Mass;
        }
#endif
    }
    
    make_list_of_active_particles();
    
    NumForceUpdate = 0;
    for (int i : ActiveParticleList)
    {
        NumForceUpdate++;
        if(i >= NumPart)
        {
            printf("inconsistent active list: i=%d >= NumPart=%d (task=%d)\n", i, NumPart, ThisTask); fflush(stdout);
            endrun(90001003);
            break;   /* graceful: stop scanning the bad list; all ranks still reach sumup_large_ints below, then drain at the next phase poll */
        }
    }
    
    sumup_large_ints(1, &NumForceUpdate, &glob_sum);
    GlobNumForceUpdate = glob_sum;
    gizmo_exit_bad_stop_if_requested("predict:reconstruct_timebins_active_list");  /* drain a bad active-list (line ~93) before callers use GlobNumForceUpdate; collective-safe (every call reaches sumup_large_ints above) */
}










/* Drift-cache state.
 *
 * g_last_full_drift_Ti tracks the time1 of the most recent full-N drift.
 * When move_particles or gizmo_full_drift_to is called with time1 <= this
 * value, the loop is short-circuited (already drifted; drift_particle would
 * bail per-particle). Reset to -1 via gizmo_full_drift_invalidate() if
 * external state advances time without a corresponding drift call.
 *
 * Lazy-drift mode (Attack C): move_particles() iterates only
 * ActiveParticleList instead of all NumPart particles, leaving non-active
 * particles' Ti_current at their previous value. drift_particle() fires
 * lazily for neighbors on-demand via gizmo_lazy_drift_for_neighbor_list()
 * inside gpu_ngb_list_build. Domain_decomp boundaries call
 * gizmo_full_drift_to() explicitly so the decomp + subsequent fresh SIDX
 * build see all-particles-at-time1 P[].Pos. */
static integertime g_last_full_drift_Ti = -1;

extern "C" void gizmo_full_drift_invalidate(void) { g_last_full_drift_Ti = -1; }

/* Time of the most recent full-N drift on this rank; -1 if none has run.
 * force_treebuild reads it to establish that the particles a tree is being
 * built from had all reached that time, which is what makes the resulting
 * node geometry current rather than merely freshly written. move_particles
 * deliberately does not advance it (it drifts only the active set), so this
 * is a proof and not a convention. */
extern "C" integertime gizmo_full_drift_ti(void) { return g_last_full_drift_Ti; }

/* Full-N drift to time1. Idempotent re-entry via g_last_full_drift_Ti cache.
 * Used by:
 *  - run.cc explicitly before each domain_Decomposition_* call so the decomp
 *    + subsequent fresh SIDX build see correct P[].Pos for every particle
 *    (lazy move_particles leaves non-active particles undrifted)
 *  - any code path that needs the full-particle set at a uniform time
 *    (output, restart, box-wrapping) */
void gizmo_full_drift_to(integertime time1)
{
    if(time1 <= g_last_full_drift_Ti) return; /* already drifted — no h change */
    /* The stamp is a CLAIM that every local particle stands at time1, and readers act
       on it by skipping work: the neighbour pool sweep treats the pool as uniform, and
       a tree built at this time may be stamped current. So it is published only when
       the drift reports that it completed with nothing pending. On a failure the stamp
       stays where it was, which costs the callers below a repeat of work they would
       otherwise skip and costs nothing else -- the alternative is a false claim that
       suppresses the very scan that would notice it. */
    if(drift_particles_batch(NULL, NumPart, time1) == 0) {g_last_full_drift_Ti = time1;}
    /* drift_particle just multiplied KernelRadius by kernel_radius_drift_factor for every
     * particle. Mark the whole pool h-dirty so the next
     * NGL build / next ghost_exchange refreshes compact_xyzh.h from current P[].
     * Conservative: covers all types at once. (A future refinement could narrow
     * to only Type 0 + AGS-active types if profiling shows this is too eager.) */
    gizmo_mark_kernel_radius_dirty_range(0, NumPart);
    /* One post-drift arena coherence point. Under UVM-canonical builds the
     * arena IS the host P/CellP (alias), so refresh is a diagnostic-counter
     * no-op; kept for forward compat with non-aliased arena builds. */
    gpu_particles_arena_refresh_from_host(NumPart, P, CellP, "move_particles_post_drift");
}

void move_particles(integertime time1)
{
    /* Tiny-N corridor counter: increments on API entry, including the
     * cache-hit early-return below. Mode B paths in run_neighbor_loop must
     * NOT enter this function regardless of whether it would do work. See
     * declarations/lifecycle_counters.h. */
    g_global_drift_counter++;

    if(time1 <= g_last_full_drift_Ti) {
        return;
    }
    /* When every local particle is active, drifting the active set IS the full drift: take the
       full-drift path so its stamp is published under that path's own proof (the batch completed
       with nothing pending) and a later full-drift request at this time finds nothing left to do. */
    if((int) ActiveParticleList.size() == NumPart) {gizmo_full_drift_to(time1); return;}
    /* Drift only the active-particle set. Non-active particles stay at their
     * previous Ti_current; drift_particle fires lazily via the
     * gpu_ngb_list_build hook when a kernel actually reads them.
     *
     * IMPORTANT: do NOT advance g_last_full_drift_Ti here — full drift
     * has not happened, so a subsequent gizmo_full_drift_to(time1) call
     * (e.g. before domain_decomp) must still execute the full loop. */
    int n_active = 0;
    /* Materialize the active list into a vector once so we can both drive the
     * drift loop AND batch-mark h-dirty after. Vector construction is O(N_active)
     * and trivially parallelizable; the drift loop dominates anyway.
     *
     * OMP over the global iterator ActiveParticleList isn't trivially
     * parallelizable (range-based-for over a custom container). Single-thread
     * host iteration is fine — drift_particle's per-particle work bounds the
     * loop, paralleling here gains little. */
    std::vector<int> active_idx;
    active_idx.reserve(ActiveParticleList.size());
    for(int i : ActiveParticleList) {
        active_idx.push_back(i);
        n_active++;
    }
    drift_particles_batch(active_idx.data(), (int) active_idx.size(), time1);
    /* drift_particle just multiplied KernelRadius for each active particle
     * (predict.cc:160,229). Mark h-dirty for both GPU SIDX tracker and host
     * glt cache via the SSOT helper. */
    if(n_active > 0) gizmo_mark_kernel_radius_dirty_indices(active_idx.data(), n_active);
}














/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */




/* evaluate_NH_from_GradRho: definition now in predict_functions.h (single source of truth).
   Include with non-inline linkage to provide externally-visible symbols. */
#undef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION
#include "predict_functions.h"

/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.  The fold itself is
 *  box_wrap_position_to_primary_image (predict_functions.h); here each x fold
 *  also carries the shearing velocity offset into the velocities.
 */
#ifdef BOX_PERIODIC
void do_box_wrapping(void)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < NumPart; i++)
    {
        struct box_wrap_x_folds folds = box_wrap_position_to_primary_image(P[i].Pos);
#ifdef BOX_SHEARING
        for(int f = 0; f < folds.up_from_below; f++)
        {
            P[i].Vel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
            P[i].dp[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset * P[i].Mass;
            if(P[i].Type==0)
            {
                CellP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                CellP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#endif
            }
        }
        for(int f = 0; f < folds.down_from_above; f++)
        {
            P[i].Vel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
            P[i].dp[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset * P[i].Mass;
            if(P[i].Type==0)
            {
                CellP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                CellP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#endif
            }
        }
#else
        (void)folds;
#endif
    }
}
#endif

#undef KOKKOS_INLINE_FUNCTION
#include "drift_particle_functions.h"

void drift_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
    drift_extra_physics_P(i, tstart, tend, dt_entr, P, CellP);
}

void drift_particle(int i, integertime time1)
{
    /* The bounds are the cached ones init_drift_table() built the tables over; recomputing
       them here would put two libm calls on every drifted particle. */
    struct DriftKickTableView tables = drift_kick_table_view(DriftTable, GravKickTable,
            DriftTable_logTimeBegin, DriftTable_logTimeMax, All.Timebase_interval, All.ComovingIntegrationOn);
    struct EosTableView eos_tables = eos_tables_view();
    drift_particle_impl(i, time1, P, CellP, &tables, &eos_tables);
    /* A particle drifted outside the domain extent would get a wrong Peano key from any tree build or
       repartition that reuses this domain; ask for a full decomposition, which re-measures the extent. */
    if(position_outside_domain_extent(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], DomainCorner[0], DomainCorner[1], DomainCorner[2], DomainLen))
    {
#ifdef _OPENMP
#pragma omp atomic write
#endif
        DomainExtentOutgrownLocal = 1;
    }
}







#ifdef MAGNETIC
/* this function is needed to control volume fluxes of the normal components of B and phi in the 
    -bad- situation where the meshless method 'faces' do not properly close (usually means you are 
    using boundary conditions that you should not) */
double Get_DtB_FaceArea_Limiter(int i, struct particle_data *pp, struct gas_cell_data *cell)
{
#ifdef HYDRO_SPH
    return 1;
#else
    /* define some variables */
    double dt_entr = get_particle_timestep_in_physical(i, P);
    /* check the magnitude of the predicted change in B-fields, vs. B-magnitude */
    Vec3<double> dB = cell[i].DtB * (dt_entr / All.cf_atime); /* converts to code units of Vol_code*B_code = Vol_phys*B_phys/a */
    double dBmag = dB.norm(), Bmag = cell[i].BPred.norm();
    /* also make sure to check the actual pressure, since if P>>B, we will need to allow larger changes in B per timestep */
    double P_BV_units = sqrt(2.*cell[i].Pressure*All.cf_a3inv)*pp[i].Mass/cell[i].Density / All.cf_a2inv;
    /* the above should be in CODE Bcode*Vol_code units! */
    double Bmag_max = DMAX(Bmag, DMIN( P_BV_units, 10.*Bmag ));
    /* now check how accurately the cell is 'closed': the face areas are ideally zero */
    double area_sum = fabs(cell[i].Face_Area[0])+fabs(cell[i].Face_Area[1])+fabs(cell[i].Face_Area[2]);
    /* but this needs to be normalized to the 'expected' area given KernelRadius */
    double area_norm = Get_Particle_Expected_Area(pp[i].KernelRadius * All.cf_atime);
    /* ok, with that in hand, define an error tolerance based on this */
    if(area_norm>0)
    {
        double area_norm_min_threshold = 0.001;
        double area_norm_weight = 200.0;
        if(area_sum/area_norm > area_norm_min_threshold)
        {
            double tol = (All.CourantFac/0.2) * DMAX( 0.01, area_norm/(area_norm_weight * area_sum) );
            tol *= Bmag_max; /* give the limiter dimensions */
            if(dBmag > tol) {return tol/dBmag;} /* now actually check if we exceed this */
        }
    }
    return 1;
#endif
}


#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Gas_PhiField(int i_particle_id) { return Get_Gas_PhiField_P(i_particle_id, P, CellP); }

double INLINE_FUNC Get_Gas_PhiField_DampingTimeInv(int i_particle_id) { return Get_Gas_PhiField_DampingTimeInv_P(i_particle_id, P, CellP); }

#endif // dedner
#endif // magnetic





/* -------------------------------------------------------------------------------------------------------------------------------------
 ------------------- the following routines are not setting the velocity, but instead are useful routines for computation of
 -------------------  various quantities needed in the mesh motion for different coordinate systems or assumed mesh shapes
 ------------------------------------------------------------------------------------------------------------------------------------- */

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
/* time-step the positions of the mesh points. this is trivial except if we are evolving the mesh points in non-cartesian coordinates
    (cylindrical or spherical) based on assumed fixed initial velocities (if HYDRO_FIX_MESH_MOTION=2 or 3),
    in which case we have to convert back and forth. */
void advect_mesh_point(int i, double dt) { advect_mesh_point_P(i, dt, P, CellP); }
#endif




#ifdef SPECIAL_POINT_WEIGHTED_MOTION
double weight_function_for_weighted_motion_smoothing(double r, int mode)
{
    return grav_weight_function_for_weighted_motion_smoothing(r, mode); /* body moved verbatim to gravtree_force_kernel.h (shared with the GPU gravity walk) */
}
#endif
