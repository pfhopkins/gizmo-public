#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include "../system/eigen_symmetric.h"

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/gpu_neighbor_list.h" /* gizmo_mark_kernel_radius_dirty_* */
#include "../system/gpu_particles_arena.h"

#include <vector>
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"


/*! This file contains the operations needed for merging/splitting gas particles/cells on-the-fly in the simulations.
    If more complicated routines, etc. are to be added to determine when (and how) splitting/merging occurs, they should also be
    added here. The split routine should also be the template for spawning new gas particles (collisionless particles are spawned
    much more easily; for those, see the star formation routines). */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#ifdef SINK_WIND_SPAWN
#define MASS_THRESHOLD_FOR_WINDPROMO(i) (DMAX(5.*target_mass_for_wind_spawning(i),0.25*All.MaxMassForParticleSplit))
#endif /* define a mass threshold for this model above which a 'hyper-element' has accreted enough to be treated as 'normal' */


/*! Here we can insert any desired criteria for particle mergers: by default, this will occur
    when particles fall below some minimum mass threshold */
int does_particle_need_to_be_merged(int i)
{
    if(P[i].Mass <= 0) {return 0;}
#ifdef PREVENT_PARTICLE_MERGE_SPLIT
    return 0;
#else
    if(P[i].Type==0) {if(CellP[i].recent_refinement_flag==1) return 0;}
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(check_if_sufficient_mergesplit_time_has_passed(i) == 0) return 0;
#endif
#ifdef GRAIN_RDI_TESTPROBLEM
    return 0;
#endif
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    if(P[i].Type==4) {return evaluate_starstar_merger_for_starcluster_eligibility(i);}
#endif
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(P[i].Type>0) {return 0;} // don't allow merging of collisionless particles [only splitting, in these runs]
    if(is_particle_a_special_zoom_target(i)) {return 0;}
#endif
#ifdef SINK_WIND_SPAWN
    if(P[i].ID==All.SpawnedWindCellID && P[i].Type==0)
    {
        if(P[i].Mass >= MASS_THRESHOLD_FOR_WINDPROMO(i)*target_mass_renormalization_factor_for_mergesplit(i,0)) {return 1;}
    }
#endif
#ifdef PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT
    if(P[i].Type==0)
    {
        double lambda_J = CellP[i].fast_MHD_wavespeed() * sqrt(M_PI / (All.G * CellP[i].Density * All.cf_a3inv));
        if((lambda_J > 4. * PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT * P[i].Get_Particle_Size()*All.cf_atime) && (P[i].Mass < All.MaxMassForParticleSplit)) {return 1;} // de-refine
    }
#endif
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT)
    if(P[i].Type>0) {return 0;} // don't allow merging of collisionless particles [only splitting, in these runs]
    if(P[i].Type==0) {
        double rmin_pc = 700.;
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
        rmin_pc = 10.;
#endif
#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES)
        rmin_pc = 0;
#endif
        if(P[i].Type==0) {if(P[i].Get_Particle_Size()*All.cf_atime*UNIT_LENGTH_IN_PC < rmin_pc) {return 0;}} // if too high-res spatially, this equiv to size for m=7000 msun for nH=1e-3, dont let de-refine
    }
#endif
    if((P[i].Type>0) && (P[i].Mass > 0.5*All.MinMassForParticleMerger*target_mass_renormalization_factor_for_mergesplit(i,0))) {return 0;}
    if(P[i].Mass <= (All.MinMassForParticleMerger*target_mass_renormalization_factor_for_mergesplit(i,0))) {return 1;}
    return 0;
#endif
}


/*! Here we can insert any desired criteria for particle splitting: by default, this will occur
    when particles become too massive, but it could also be done when KernelRadius gets very large, densities are high, etc */
int does_particle_need_to_be_split(int i)
{
    if(P[i].Type != 0) {return 0;} // default behavior: only gas particles split //
#ifdef PREVENT_PARTICLE_MERGE_SPLIT
    return 0;
#else
    if(P[i].Type==0) {if(CellP[i].recent_refinement_flag==1) return 0;}
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(check_if_sufficient_mergesplit_time_has_passed(i) == 0) return 0;
#endif
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    if(P[i].Type==4) {return 0;}
#endif
    if(P[i].Mass >= (All.MaxMassForParticleSplit*target_mass_renormalization_factor_for_mergesplit(i,1))) {return 1;}
#ifdef PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT
    if(P[i].Type == 0)
    {
        double lambda_J = CellP[i].fast_MHD_wavespeed() * sqrt(M_PI / (All.G * CellP[i].Density * All.cf_a3inv));
        if((lambda_J < PARTICLE_MERGE_SPLIT_TRUELOVE_REFINEMENT * P[i].Get_Particle_Size()*All.cf_atime) && (P[i].Mass > 2*All.MinMassForParticleMerger)) {return 1;} // refine
    }
#endif
    return 0;
#endif
}

/*! A multiplicative factor that determines the target mass of a particle for the (de)refinement routines; split_key tells you if this is for a split (1) or merge (0) */
double target_mass_renormalization_factor_for_mergesplit(int i, int split_key)
{
    double ref_factor=1.0;
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
    ref_factor = 1; // need to determine appropriate desired refinement criterion, if resolution is not strictly pre-defined //
#endif
    
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    if(P[i].Type==0)
    {
        double dt_to_ramp_refinement = 0.00001;
        dt_to_ramp_refinement = 1.e-6; // testing [for specific time chosen, hard-coded but change as needed]
        double minimum_refinement_mass_in_solar = 0.003; // aims at 0.01 effective
        
        double mcrit_0=1.*(4000.), T_eff = 1.23 * (5./3.-1.) * U_TO_TEMP_UNITS * CellP[i].InternalEnergyPred, nH_cgs = CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, MJ = 9.e6 * pow( 1 + T_eff/1.e4, 1.5) / sqrt(1.e-12 + nH_cgs);
        if(All.ComovingIntegrationOn) {MJ *= pow(1. + (100.*COSMIC_BARYON_DENSITY_CGS) / (CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS), 3);}
        double m_ref_mJ = 0.001 * MJ;
        int j; double r2=0,r2min=MAX_REAL_NUMBER;
        for(j=0;j<SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM;j++)
        {
            r2 = (P[i].Pos - All.SpecialParticle_Position_ForRefinement[j]).norm_sq() * All.cf_atime * All.cf_atime;
            if(r2<r2min) {r2min=r2;}
        }
        r2 = r2min; // want minimum distance to nearest refinement center
        
        double rbh = sqrt(r2) * UNIT_LENGTH_IN_PC/1000.;
        if(rbh > 1.e-10 && isfinite(rbh) && rbh < 1.e10)
        {
            double mc=1.e10, m_r1=DMIN(mcrit_0, 7.e3), m_r2=10.*m_r1, m_r3=10.*m_r2, r1=1., r2=10., r3=20.;
            if(rbh<r1) {mc=m_r1;} else {if(rbh<r2) {mc=m_r1*exp(log(m_r2/m_r1)*log(rbh/r1)/log(r2/r1));} else
            {if(rbh<r3) {mc=m_r2*exp(log(m_r3/m_r2)*log(rbh/r2)/log(r3/r2));} else {mc=m_r3*pow(rbh/r3,3);}}}
            m_ref_mJ = DMIN(m_ref_mJ , mc);
        }
        double r_pc = rbh*1000.,r0, f0=1, target_slope=1.2;
        double slope=0; slope = target_slope * (1. - exp(-(All.Time - All.TimeBegin) / dt_to_ramp_refinement)); // gradually ramp up refinement from snapshot

        double t_00 = All.TimeBegin;
#if(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=1)
        t_00 = 0.1843705; // testing [for specific time chosen, hard-coded but change as needed]
#endif
#if(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=3)
        t_00 = -0.2; // want this to be negative so skip the ramp step for these runs
#endif
        double dtau = (All.Time - t_00) / dt_to_ramp_refinement, dtdelay=0.1, tfinal=1.;
        if(dtau < dtdelay) {slope=0;} else {slope=target_slope * (1. - exp(- ((dtau-dtdelay) / (tfinal - dtdelay)) ));} // alt model

        if(dtau < dtdelay) {slope=0;} else {if(dtau>tfinal) {slope=target_slope;} else {slope = target_slope * (dtau-dtdelay) / (tfinal - dtdelay);}}
        r0=1000.; if(r_pc<r0) {f0 *= pow(r_pc/r0,slope);}
        if(dtau < 2.*dtdelay) {slope*=0;} else if(dtau < 3.*dtdelay) {slope*=(dtau-2.*dtdelay)/dtdelay;}
        r0=100.; if(r_pc<r0) {f0 *= pow(r_pc/r0,slope);}
        if(dtau < 4.*dtdelay) {slope*=0;} else if(dtau < 5.*dtdelay) {slope*=(dtau-4.*dtdelay)/dtdelay;}
        r0=10.; if(r_pc<r0) {f0 *= pow(r_pc/r0,slope);}
        if(dtau < 6.*dtdelay) {slope*=0;} else if(dtau < 7.*dtdelay) {slope*=(dtau-6.*dtdelay)/dtdelay;}
        r0=1.; if(r_pc<r0) {f0 *= pow(r_pc/r0,slope);}
        
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=2) /* simpler refinement for even smaller-scale simulations */
        mcrit_0 = m_ref_mJ = 0.003; minimum_refinement_mass_in_solar = 4.e-8; f0 = 1; // 'baseline' resolution & maximum resolution target
        r0 = 0.01; if(r_pc < r0) {f0 *= pow(r_pc/r0 , 2);} // simple additional refinement criterion vs radius interior to inner radius
        f0 = DMAX(pow(r_pc / r0, 4) , 1.e-4); // optional alt - using for now, can revert to above.
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=3)
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=4)
        mcrit_0 = m_ref_mJ = 0.009; f0 = 1; r0 = 0.3; if(r_pc > r0) {f0 *= pow(r_pc/r0 , 1.5);}
#endif
        f0 = 1;
        r0 = 0.3; if(r_pc > r0) {f0 = pow(r_pc/r0 , 1.5);}
        double r0p=r0; r0 = 200.; if(r_pc > r0) {f0 = pow(r0/r0p , 1.5) * pow(r_pc/r0 , 0.75);}
        r0 = 0.1; if(r_pc < r0) {f0 *= pow(r_pc/r0 , 1);}
        r0 = 0.01; if(r_pc < r0) {f0 *= pow(r_pc/r0 , 2);}
        r0 = 0.0033; if(r_pc < r0) {f0 *= pow(r_pc/r0 , 2);}
        f0 = DMAX(f0, 1.e-4);
        double f0minfac = 30.;
        r0 = 0.001; if(r_pc < r0) {f0minfac *= pow(r_pc/r0 , 2);}
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES < 4)
        f0minfac = DMAX(f0minfac , 0.005); // may need to be further lowered later
        minimum_refinement_mass_in_solar = 1.e-9;
        f0 = DMIN(DMAX(1.,f0), f0minfac*f0);
#else
        f0minfac = DMAX(f0minfac , 0.015); // may need to be further lowered later
        f0 = DMIN(DMAX(DMAX(1.,f0),1./m_ref_mJ), f0minfac*f0);
#endif
#endif
#endif

        double M_target = f0 * DMAX( mcrit_0, m_ref_mJ ) / UNIT_MASS_IN_SOLAR;
        double M_min_absolute = minimum_refinement_mass_in_solar / UNIT_MASS_IN_SOLAR; // arbitrarily set minimum mass for refinement at any level
        double normal_median_mass = All.MaxMassForParticleSplit / 3.;
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
        ref_factor = DMAX(M_min_absolute / normal_median_mass, DMIN( M_target / normal_median_mass , 0.1/(normal_median_mass*UNIT_MASS_IN_SOLAR)));
#else
        ref_factor = DMAX(M_min_absolute / normal_median_mass, DMIN( M_target / normal_median_mass , 1));
#endif
        return ref_factor;
    }
#endif
    
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT)
    if(P[i].Type==0)
    {
        double mcrit_0=1.*(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT), T_eff = 1.23 * (5./3.-1.) * U_TO_TEMP_UNITS * CellP[i].InternalEnergyPred, nH_cgs = CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS;
        double T_min_jeans = 1.e4; // dont use values below this for thermal jeans mass in refinement criterion
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
        T_min_jeans = 0.; // allow as cold as needed
#endif
        double MJ = 9.e6 * pow( (T_eff + T_min_jeans)/1.e4, 1.5) / sqrt(1.e-12 + nH_cgs); // Jeans mass (in solar), but modified with lower limit for temperature so we refine all cool gas equally, lower limit for numerical convenience for density
        if(All.ComovingIntegrationOn) {MJ *= pow(1. + (100.*COSMIC_BARYON_DENSITY_CGS) / (CellP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS), 3);} // ensure that only cells much denser than cosmic mean are eligible for refinement. use 100x so even cells outside Rvir are potentially eligible
        // to check against hot gas in high-density ISM getting worse than a certain resolution level, we want to check that we don't down-grade the spatial resolution too much
        double m_ref_mJ = 0.001 * MJ;
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
        m_ref_mJ = 0.01 * MJ; // since using real thermal jeans, fully-cold, don't need to be as aggressive here
#endif
#if defined(SINK_CALC_DISTANCES)
        double rbh = P[i].Min_Distance_to_Sink * All.cf_atime; // distance to nearest BH
        if(rbh > 1.e-10 && isfinite(rbh) && rbh < 1.e10)
        {
            double r1=1., r2=10., r3=20.; // boundaries for different distance-based refinement thresholds, in physical kpc
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL)
            r1=0.1; r2=1.; r3=10.; // need to be customized for the problem, here refining more so smaller zone size
#endif
            double mc=1.e10, m_r1=DMIN(mcrit_0, 7.e3), m_r2=10.*m_r1, m_r3=10.*m_r2; r1/=UNIT_LENGTH_IN_KPC; r2/=UNIT_LENGTH_IN_KPC; r3/=UNIT_LENGTH_IN_KPC;
            if(rbh<r1) {mc=m_r1;} else {if(rbh<r2) {mc=m_r1*exp(log(m_r2/m_r1)*log(rbh/r1)/log(r2/r1));} else
                {if(rbh<r3) {mc=m_r2*exp(log(m_r3/m_r2)*log(rbh/r2)/log(r3/r2));} else {mc=m_r3*pow(rbh/r3,3);}}}
            
            m_ref_mJ = DMIN(m_ref_mJ , mc);
        }
#endif
        double M_target = DMAX( mcrit_0, m_ref_mJ ) / UNIT_MASS_IN_SOLAR; // enforce minimum refinement to 7000 Msun, and convert to code units, compare to 0.001xJeans mass, which is designed to target desired levels
        double normal_median_mass = All.MaxMassForParticleSplit / 3.; // code median mass from ICs
        ref_factor = DMAX(1.e-30, DMIN( M_target / normal_median_mass , 1)); // this shouldn't get larger than unity since that would exceed the normal maximum mass
        return ref_factor; // return it
    }
#endif
    
/*!
 #if defined(SINK_CALC_DISTANCES) && !defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) && !defined(SINGLE_STAR_SINK_DYNAMICS)
    ref_factor = DMIN(1.,sqrt(P[i].Min_Distance_to_Sink + 0.0001)); // this is an example of the kind of routine you could use to scale resolution with BH distance //
#endif
 */

    return ref_factor;
}


/*! This is the parent routine to actually determine if mergers/splits need to be performed, and if so, to do them
  modified by Takashi Okamoto (t.t.okamoto@gmail.com) on 20/6/2019
 */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */
void merge_and_split_particles(void)
{
    struct flags_merg_split {
        int flag; // 0 for nothing, -1 for clipping, 1 for merging, 2 for splitting, and 3 marked as merged
        int target_index;
    } *Ptmp;

    int target_for_merger,dummy=0,numngb_inbox,startnode,i,j,n; double threshold_val;
    int n_particles_merged,n_particles_split,n_particles_gas_split,MPI_n_particles_merged,MPI_n_particles_split,MPI_n_particles_gas_split;
    (void)dummy; (void)numngb_inbox; (void)startnode; (void)n; /* may go unused on modern path */
    Gas_split=0; n_particles_merged=0; n_particles_split=0; n_particles_gas_split=0; MPI_n_particles_merged=0; MPI_n_particles_split=0; MPI_n_particles_gas_split=0;
    Ptmp = (struct flags_merg_split *) mymalloc("Ptmp", NumPart * sizeof(struct flags_merg_split));

    // TO: need initialization
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < NumPart; i++) {
      Ptmp[i].flag = 0;
      Ptmp[i].target_index = -1;
    }

    /* Modern path: prebuilt CSR neighbor list. Skip ghosts in the inner loop —
     * merge_particles_ij/split_particle_i operate on local P[] indices only.
     * Bitmask is OR of all merge/split-eligible types (gas + GALSF stars);
     * type-equality filter inside the per-source loop keeps semantics identical
     * to legacy ngb_treefind_variable_targeted with bitflag = 1<<P[i].Type. */
    {
        std::vector<int> ms_src_idx;
        std::vector<double> ms_src_radii;
        ms_src_idx.reserve(64);
        ms_src_radii.reserve(64);

        int merge_split_types_mask = (1 << 0) /* gas */
#if defined(GALSF)
                                   | (1 << 4) /* stars (SF active) */
#endif
                                   ;

        for (int ip = 0; ip < NumPart; ip++) {
            if (P[ip].Mass <= 0) continue;
#if defined(GALSF)
            if (!(((P[ip].Type==0)||(P[ip].Type==4)) && TimeBinActive[P[ip].TimeBin])) continue;
#else
            if (!((P[ip].Type==0) && TimeBinActive[P[ip].TimeBin])) continue;
#endif
            if (!(does_particle_need_to_be_merged(ip) || does_particle_need_to_be_split(ip))) continue;
            if (P[ip].KernelRadius <= 0) continue;
            ms_src_idx.push_back(ip);
            ms_src_radii.push_back(P[ip].KernelRadius);
        }

        int num_src = (int)ms_src_idx.size();
        gpu_neighbor_list_t gnl = {};
        std::vector<int> gnl_neighbors_host;
        int local_count = ghost_get_num_local();
        if (local_count <= 0) local_count = NumPart;
        if (num_src > 0) {
            int num_all = local_count + ghost_get_num_ghosts();
            if (num_all <= 0) num_all = NumPart;
            gpu_particles_arena_acquire(num_all, P, CellP);
            struct particle_data *P_gpu = gpu_particles_arena_P();
            gpu_ngb_list_build(P_gpu, num_all,
                               ms_src_idx.data(), num_src,
                               NGB_SEARCH_ONEWAY, merge_split_types_mask,
                               &gnl, NULL, 1.0, ms_src_radii.data(), NULL, "msplit");
            /* gnl.neighbors is DEVICE_SPACE; host loop below indexes it. */
            if (gnl.total_pairs > 0) {
                gnl_neighbors_host.resize((size_t)gnl.total_pairs);
                gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
            }
        }
        const int *gnl_neighbors = gnl_neighbors_host.empty() ? NULL : gnl_neighbors_host.data();

        for (int aa = 0; aa < num_src; aa++) {
            i = ms_src_idx[aa];
            int64_t nl_start = gnl.offsets[aa], nl_end = gnl.offsets[aa+1];

            if (does_particle_need_to_be_merged(i)) {
                target_for_merger = -1;
                threshold_val = MAX_REAL_NUMBER;
                for (int64_t nn = nl_start; nn < nl_end; nn++) {
                    j = gnl_neighbors[nn];
                    if (j >= local_count) continue; /* skip ghosts (local-only merge) */
                    if (P[j].Type != P[i].Type) continue;
                    double m_eff = P[j].Mass; int do_allow_merger = 0;
                    if ((P[j].Mass >= P[i].Mass) && (P[i].Mass+P[j].Mass < All.MaxMassForParticleSplit)) {do_allow_merger = 1;}
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
                    if (P[i].Type==4 && P[j].Type==4) {m_eff=evaluate_starstar_merger_for_starcluster_particle_pair(i,j); if(m_eff<=0) {do_allow_merger=0;} else {do_allow_merger=1;}}
#endif
#ifdef SINK_WIND_SPAWN
                    if (P[i].ID==All.SpawnedWindCellID && P[i].Type==0) {
                        if (P[i].Mass>=MASS_THRESHOLD_FOR_WINDPROMO(i)) {
                            if ((P[j].ID != All.SpawnedWindCellID) || (P[j].Mass >= MASS_THRESHOLD_FOR_WINDPROMO(j))) {do_allow_merger *= 1;} else {do_allow_merger = 0;}
                        } else if (do_allow_merger) {
                            Vec3<MyDouble> dvel_tmp = P[i].Vel - P[j].Vel; double v2_tmp = dvel_tmp.norm_sq(); double vr_tmp = dot(dvel_tmp, P[i].Pos - P[j].Pos);
                            if (vr_tmp > 0) {do_allow_merger = 0;}
                            if (v2_tmp > 0) {v2_tmp=sqrt(v2_tmp*All.cf_a2inv);} else {v2_tmp=0;}
                            if (v2_tmp > DMIN(CellP[i].effective_soundspeed(),CellP[j].effective_soundspeed())) {do_allow_merger = 0;}
#if !defined(SINK_RIAF_SUBEDDINGTON_MODEL) && !defined(SINGLE_STAR_SINK_DYNAMICS)
                            if (P[j].ID == All.SpawnedWindCellID) {do_allow_merger = 0;}
#if !defined(SINGLE_STAR_FB_JETS) && !defined(SINGLE_STAR_FB_WINDS)
                            if ((v2_tmp > 0.25*All.Sink_outflow_velocity) && (v2_tmp > 0.9*CellP[j].effective_soundspeed())) {do_allow_merger=0;}
#endif
#endif
                        }
                    }
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
                    if (P[i].Type==0 && P[j].Type==0) {
                        if ((P[j].Mass >= P[i].Mass) && (P[i].Mass+P[j].Mass < All.MaxMassForParticleSplit)) {
                            double dti = get_particle_timestep_in_physical(i), dtj=get_particle_timestep_in_physical(j);
                            double dt0 = 1.e-7;
                            if (dti<dt0 || dtj<dt0) {do_allow_merger = 1;}
                        }}
#endif
                    if (P[j].ID==All.SpawnedWindCellID && P[j].Type==0) {m_eff *= 1.0e10;}
#endif
                    if ((j<0)||(j==i)||(P[j].Type!=P[i].Type)||(P[j].Mass<=0)||(Ptmp[j].flag!=0)||(m_eff>=threshold_val)) {do_allow_merger=0;}
                    if (do_allow_merger) {threshold_val=m_eff; target_for_merger=j;}
                }
                if (target_for_merger >= 0) {
                    Ptmp[i].flag = 1; Ptmp[target_for_merger].flag = 3; Ptmp[i].target_index = target_for_merger;
                }
            }
            else if (does_particle_need_to_be_split(i) && (Ptmp[i].flag == 0)) {
                target_for_merger = -1;
                threshold_val = MAX_REAL_NUMBER;
                for (int64_t nn = nl_start; nn < nl_end; nn++) {
                    j = gnl_neighbors[nn];
                    if (j >= local_count) continue;
                    if (P[j].Type != P[i].Type) continue;
                    if ((j>=0)&&(j!=i)&&(P[j].Type==P[i].Type) && (P[j].Mass > 0) && (Ptmp[j].flag == 0)) {
                        Vec3<double> dp = P[i].Pos - P[j].Pos;
                        nearest_xyz(dp);
                        double r2 = dp.norm_sq();
                        if (r2<threshold_val) {threshold_val=r2; target_for_merger=j;}
                    }
                }
                if (target_for_merger >= 0) {
                    Ptmp[i].flag = 2;
                    Ptmp[i].target_index = target_for_merger;
                }
            }
        }

        if (num_src > 0) {
            gpu_ngb_list_free(&gnl, NULL);
            gpu_particles_arena_invalidate();
        }
    }

    // actual merge-splitting loop loop. No tree-walk is allowed below here
    int failed_splits = 0; /* record failed splits to output warning message */
    for (i = 0; i < NumPart; i++) {
        if (Ptmp[i].flag == 1) { // merge this particle
            int did_merge = merge_particles_ij(i, Ptmp[i].target_index);
            if(did_merge == 1) {n_particles_merged++;}
        }
        if (Ptmp[i].flag == 2) {
            int did_split = split_particle_i(i, n_particles_split, Ptmp[i].target_index);
            if(did_split == 1) {n_particles_split++; if(P[i].Type==0) {n_particles_gas_split++;}} else {failed_splits++;}
        }
    }
    if(failed_splits) {printf ("On Task=%d with NumPart=%d we tried and failed to split %d elements, after running out of space (REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPart=%d, REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPartGas=%d ).\n We did split %d total (%d gas) elements. Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, failed_splits, (int)(REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPart), (int)(REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPartGas), n_particles_split, n_particles_gas_split); fflush(stdout);}

#ifdef BOX_PERIODIC
    /* map the particles back onto the box (make sure they get wrapped if they go off the edges). this is redundant here,
     because we only do splits in the beginning of a domain decomposition step, where this will be called as soon as
     the particle re-order is completed. but it is still useful to keep here in case this changes (and to note what needs
     to be done for any more complicated splitting operations */
    do_box_wrapping();
#endif
    myfree(Ptmp);
    
    MPI_Allreduce(&n_particles_merged, &MPI_n_particles_merged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_particles_gas_split, &MPI_n_particles_gas_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(MPI_n_particles_merged > 0 || MPI_n_particles_split > 0)
        {
            printf("Particle split/merge check: %d particles merged, %d particles split (%d gas) \n", MPI_n_particles_merged,MPI_n_particles_split,MPI_n_particles_gas_split);
        }
    }
    /* the reduction or increase of n_part by MPI_n_particles_merged will occur in rearrange_particle_sequence, which -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas += (long long)MPI_n_particles_gas_split;
    Gas_split = n_particles_gas_split; // specific to the local processor //
    NumPart += (n_particles_split - n_particles_gas_split); // specific to the local processor; note the gas split number will be added below, this is just non-gas splits //
    gpu_particles_arena_invalidate(); /* NumPart and P[] indices changed; arena stale */
}




/*! This is the routine that does the particle splitting. Note this is a tricky operation if we're not using meshes to divide the volume,
    so care needs to be taken modifying this so that it's done in a way that is (1) conservative, (2) minimizes perturbations to the
    volumetric quantities of the flow, and (3) doesn't crash the tree or lead to particle 'overlap'
    Modified by Takashi Okamoto on 20/6/2019.  */
int split_particle_i(int i, int n_particles_split, int i_nearest)
{
    double mass_of_new_particle;
    if( ((P[i].Type==0) && (NumPart + n_particles_split + 1 >= (int)(REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPartGas))) || (NumPart + n_particles_split + 1 >= (int)(REDUC_FAC_FOR_MEMORY_IN_DOMAIN*All.MaxPart)) )
    {
        //printf ("On Task=%d with NumPart=%d we tried to split a particle, but there is no space left...(All.MaxPart=%d). Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, All.MaxPart); fflush(stdout);
        return 0;
        endrun(8888);
    }
#ifndef SPAWN_PARTICLES_VIA_SPLITTING
    if(P[i].Type != 0) {printf("Splitting Non-Gas Particle: i=%d ID=%llu Type=%d \n",i,(unsigned long long) P[i].ID,P[i].Type);} //fflush(stdout); endrun(8889);
#endif

    /* here is where the details of the split are coded, the rest is bookkeeping */
    mass_of_new_particle = 0.5;

    int k; double phi,cos_theta;
    k=0;
    phi = 2.0*M_PI*get_random_number(i+1+ThisTask); // random from 0 to 2pi //
    cos_theta = 2.0*(get_random_number(i+3+2*ThisTask)-0.5); // random between 1 to -1 //
    Vec3<double> dp = P[i].Pos - P[i_nearest].Pos;
    nearest_xyz(dp);
    double r_near = dp.norm();
    /* set d_r to the Voronoi-optimal daughter offset: d = (V_daughter)^{1/NDIMS}/2 = 2^{-(NDIMS+1)/NDIMS} * dx_eff,
       placing each daughter at its sub-cell centroid. minimizes O(d^2) perturbations to neighbor density estimates
       while maintaining a glass-like configuration. dimension-general: 0.397*dx_eff in 3D, 0.354 in 2D, 0.25 in 1D */
    double dx_eff = P[i].Get_Particle_Size(), dx_h = KERNEL_CORE_SIZE * P[i].KernelRadius;
    dx_eff = DMAX(DMIN(dx_eff, 3.*dx_h), 0.1*dx_h); // clamp to reasonable range relative to kernel radius
    dx_h = r_near; dx_eff = DMAX(DMIN(dx_eff, 3.*dx_h), 0.1*dx_h); // clamp to reasonable range relative to nearest-neighbor distance
    double d_r = pow(0.5, (NUMDIMS + 1.) / (double)NUMDIMS) * dx_eff; // 2^{-(NDIMS+1)/NDIMS} * dx_eff
    d_r = DMIN(d_r, 0.4 * r_near); // keep displacement well within nearest-neighbor distance
    /*
    double r_near = sqrt(r2_nearest);
    double rkern = P[i].Get_Particle_Size();
    if(rkern < r_near) {rkern = r_near;}
    r_near *= 0.35;
    double d_r = 0.25 * rkern; // needs to be epsilon*KernelRadius where epsilon<<1, to maintain stability //
    d_r = DMAX( DMAX(0.1*r_near , 0.005*rkern) , DMIN(d_r , r_near) ); // use a 'buffer' to limit to some multiple of the distance to the nearest particle //
    */ // the change above appears to cause some numerical instability //
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * ForceSoftening_KernelRadius(i));
#endif
#ifdef BOX_BND_PARTICLES
    if(P[i].Type != 0 && P[i].ID == 0) {d_r *= 1.e-3;}
#endif

    /* find the first non-gas particle and move it to the end of the particle list */
    long j = NumPart + n_particles_split;
    /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
    P[j] = P[i];
    //memcpy(P[j],P[i],sizeof(struct particle_data)); // safer copy to make sure we don't just end up with a pointer re-direct

#ifdef CHIMES
    int abunIndex;
    ChimesGasVars[j] = ChimesGasVars[i];
    allocate_gas_abundances_memory(&(ChimesGasVars[j]), &ChimesGlobalVars);
    for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++) {ChimesGasVars[j].abundances[abunIndex] = ChimesGasVars[i].abundances[abunIndex];}
#endif

    // need to assign new particle a unique ID:
    // new method: preserve the original "ID" field, but assign a unique -child- ID: this is unique up to ~32 *GENERATIONS* of repeated splitting!
    P[j].ID_child_number = P[i].ID_child_number + (MyIDType)(1 << ((int)P[i].ID_generation)); // particle 'i' retains its child number; this ensures uniqueness
    P[i].ID_generation = P[i].ID_generation + 1;
    if(P[i].ID_generation > 30) {P[i].ID_generation=0;} // roll over at 32 generations (unlikely to ever reach this)
    P[j].ID_generation = P[i].ID_generation; // ok, all set!

    /* assign masses to both particles (so they sum correctly), but first record some numbers in case we need them below */
    double mass_before_split = 0, density_before_split = 0, volume_before_split = 0; mass_before_split = P[i].Mass; // save for use below
    if(P[i].Type==0) {density_before_split = CellP[i].Density; volume_before_split = mass_before_split/density_before_split;} // save for use below
    P[j].Mass = mass_of_new_particle * P[i].Mass; CellP[j].Mass = P[j].Mass;
    P[i].Mass -= P[j].Mass; CellP[i].Mass = P[i].Mass;
#ifdef BOX_BND_PARTICLES
    if(P[i].ID==0)  {P[j].ID=1; double m0=P[i].Mass+P[j].Mass; P[i].Mass=P[j].Mass=m0;} CellP[i].Mass = P[i].Mass;
#endif

    /* prepare to shift the particle locations according to the random number we drew above */
    double dx, dy, dz;
#if (NUMDIMS == 1)
    dy=dz=0; dx=d_r; // here the split direction is trivial //
#else
    /* in 2D and 3D its not so trivial how to split the directions */
    double sin_theta = sqrt(1 - cos_theta*cos_theta);
    dx = d_r * sin_theta * cos(phi);
    dy = d_r * sin_theta * sin(phi);
    dz = d_r * cos_theta;
#if (NUMDIMS == 2)
    dz=0; dx=d_r*cos(phi); dy=d_r*sin(phi);
#endif
#endif

    if(P[i].Type==0)
    {
        /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
        CellP[j] = CellP[i];
        //memcpy(CellP[j],CellP[i],sizeof(struct gas_cell_data)); // safer copy to make sure we don't just end up with a pointer re-direct
        /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
        CellP[i].ConditionNumber *= 10.0;
        CellP[j].ConditionNumber = CellP[i].ConditionNumber;
#ifdef MAGNETIC
        /* we evolve the -conserved- VB and Vphi, so this must be partitioned */
        {
            auto B_before_split = CellP[i].B * (1.0 / volume_before_split); /* calculate the real value of B pre-split to know what we need to correctly re-initialize to once the volume partition can be recomputed */
            CellP[j].BField_prerefinement = B_before_split; CellP[i].BField_prerefinement = B_before_split; /* record the real value of B pre-split */
            CellP[j].B = CellP[i].B * mass_of_new_particle; CellP[i].B -= CellP[j].B; /* take a reasonable -guess- for the new updated conserved B */
            CellP[j].BPred = CellP[i].BPred * mass_of_new_particle; CellP[i].BPred -= CellP[j].BPred; CellP[j].DtB = CellP[i].DtB * mass_of_new_particle; CellP[i].DtB -= CellP[j].DtB; /* partition these terms as well, doesn't have much effect if zero them, since before re-calc anyways */
            //CellP[j].BPred = CellP[j].B; CellP[i].BPred = CellP[i].B; CellP[j].DtB = CellP[i].DtB = {}; /* set BPred equal to the conserved B, and zero the time derivatives and cleaning terms */
        }
        CellP[j].divB = mass_of_new_particle * CellP[i].divB; CellP[i].divB -= CellP[j].divB; //CellP[i].divB = CellP[j].divB = 0; /* this will be self-consistently recomputed on the next timestep */
#ifdef DIVBCLEANING_DEDNER
        CellP[j].Phi = mass_of_new_particle * CellP[i].Phi; CellP[i].Phi -= CellP[j].Phi; CellP[j].DtPhi = mass_of_new_particle * CellP[i].DtPhi; CellP[i].DtPhi -= CellP[j].DtPhi; CellP[j].PhiPred = mass_of_new_particle * CellP[i].PhiPred; CellP[i].PhiPred -= CellP[j].PhiPred; /* same partition as above */
        CellP[j].DtB_PhiCorr = CellP[i].DtB_PhiCorr * mass_of_new_particle; CellP[i].DtB_PhiCorr -= CellP[j].DtB_PhiCorr; /* same partition as above */
        //CellP[i].Phi = CellP[i].PhiPred = CellP[i].DtPhi = CellP[j].Phi = CellP[j].PhiPred = CellP[j].DtPhi = 0; for(k=0;k<3;k++) {CellP[i].DtB_PhiCorr[k] = CellP[j].DtB_PhiCorr[k] = 0;} /* zero the time derivatives and cleaning terms: these will be re-calculated self-consistently from the new mesh configuration */
#endif
        /* ideally, particle-splits should be accompanied by a re-partition of the density via the density() call for the particles affected, after the tree-reconstruction, with quantities like B used to re-calculate after */
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            int k_dir; k_dir=0;
            CellP[j].Rad_E_gamma[k] = mass_of_new_particle * CellP[i].Rad_E_gamma[k]; CellP[i].Rad_E_gamma[k] -= CellP[j].Rad_E_gamma[k];
#if defined(RT_EVOLVE_ENERGY)
            CellP[j].Rad_E_gamma_Pred[k] = mass_of_new_particle * CellP[i].Rad_E_gamma_Pred[k]; CellP[i].Rad_E_gamma_Pred[k] -= CellP[j].Rad_E_gamma_Pred[k];
            CellP[j].Dt_Rad_E_gamma[k] = mass_of_new_particle * CellP[i].Dt_Rad_E_gamma[k]; CellP[i].Dt_Rad_E_gamma[k] -= CellP[j].Dt_Rad_E_gamma[k];
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k_dir=0;k_dir<3;k_dir++)
            {
                CellP[j].Rad_Flux[k][k_dir] = mass_of_new_particle * CellP[i].Rad_Flux[k][k_dir]; CellP[i].Rad_Flux[k][k_dir] -= CellP[j].Rad_Flux[k][k_dir];
                CellP[j].Rad_Flux_Pred[k][k_dir] = mass_of_new_particle * CellP[i].Rad_Flux_Pred[k][k_dir]; CellP[i].Rad_Flux_Pred[k][k_dir] -= CellP[j].Rad_Flux_Pred[k][k_dir];
                CellP[j].Dt_Rad_Flux[k][k_dir] = mass_of_new_particle * CellP[i].Dt_Rad_Flux[k][k_dir]; CellP[i].Dt_Rad_Flux[k][k_dir] -= CellP[j].Dt_Rad_Flux[k][k_dir];
            }
#endif
#ifdef RT_EVOLVE_INTENSITIES
            for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++)
            {
                CellP[j].Rad_Intensity[k][k_dir] = mass_of_new_particle * CellP[i].Rad_Intensity[k][k_dir]; CellP[i].Rad_Intensity[k][k_dir] -= CellP[j].Rad_Intensity[k][k_dir];
                CellP[j].Rad_Intensity_Pred[k][k_dir] = mass_of_new_particle * CellP[i].Rad_Intensity_Pred[k][k_dir]; CellP[i].Rad_Intensity_Pred[k][k_dir] -= CellP[j].Rad_Intensity_Pred[k][k_dir];
                CellP[j].Dt_Rad_Intensity[k][k_dir] = mass_of_new_particle * CellP[i].Dt_Rad_Intensity[k][k_dir]; CellP[i].Dt_Rad_Intensity[k][k_dir] -= CellP[j].Dt_Rad_Intensity[k][k_dir];
            }
#endif
        }
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        double dmass = mass_of_new_particle * CellP[i].DtMass;
        CellP[j].DtMass = dmass;
        CellP[i].DtMass -= dmass;
        dmass = mass_of_new_particle * CellP[i].dMass;
        CellP[j].dMass = dmass;
        CellP[i].dMass -= dmass;
        CellP[j].GravWorkTerm = {}; CellP[i].GravWorkTerm = {}; //= mass_of_new_particle * ...; //appears more stable with this zero'd
        CellP[j].MassTrue = mass_of_new_particle * CellP[i].MassTrue;
        CellP[i].MassTrue -= CellP[j].MassTrue;
#endif
#ifdef COSMIC_RAY_FLUID
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
        CellP[j].DtCREgyNewInjectionFromShocks = mass_of_new_particle * CellP[i].DtCREgyNewInjectionFromShocks; CellP[i].DtCREgyNewInjectionFromShocks -= CellP[j].DtCREgyNewInjectionFromShocks;
#endif
#if defined(SINK_CR_INJECTION_AT_TERMINATION)
        CellP[j].Sink_CR_Energy_Available_For_Injection = mass_of_new_particle * CellP[i].Sink_CR_Energy_Available_For_Injection; CellP[i].Sink_CR_Energy_Available_For_Injection -= CellP[j].Sink_CR_Energy_Available_For_Injection;
#endif
        int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {
            CellP[j].CosmicRayEnergy[k_CRegy] = mass_of_new_particle * CellP[i].CosmicRayEnergy[k_CRegy]; CellP[i].CosmicRayEnergy[k_CRegy] -= CellP[j].CosmicRayEnergy[k_CRegy];
            CellP[j].CosmicRayEnergyPred[k_CRegy] = mass_of_new_particle * CellP[i].CosmicRayEnergyPred[k_CRegy]; CellP[i].CosmicRayEnergyPred[k_CRegy] -= CellP[j].CosmicRayEnergyPred[k_CRegy];
            CellP[j].DtCosmicRayEnergy[k_CRegy] = mass_of_new_particle * CellP[i].DtCosmicRayEnergy[k_CRegy]; CellP[i].DtCosmicRayEnergy[k_CRegy] -= CellP[j].DtCosmicRayEnergy[k_CRegy];
#if defined(CRFLUID_EVOLVE_SPECTRUM)
            CellP[j].CosmicRay_Number_in_Bin[k_CRegy] = mass_of_new_particle * CellP[i].CosmicRay_Number_in_Bin[k_CRegy]; CellP[i].CosmicRay_Number_in_Bin[k_CRegy] -= CellP[j].CosmicRay_Number_in_Bin[k_CRegy];
            CellP[j].DtCosmicRay_Number_in_Bin[k_CRegy] = mass_of_new_particle * CellP[i].DtCosmicRay_Number_in_Bin[k_CRegy]; CellP[i].DtCosmicRay_Number_in_Bin[k_CRegy] -= CellP[j].DtCosmicRay_Number_in_Bin[k_CRegy];
#endif
            CellP[j].CosmicRayFlux[k_CRegy] = CellP[i].CosmicRayFlux[k_CRegy] * mass_of_new_particle; CellP[i].CosmicRayFlux[k_CRegy] -= CellP[j].CosmicRayFlux[k_CRegy];
            CellP[j].CosmicRayFluxPred[k_CRegy] = CellP[i].CosmicRayFluxPred[k_CRegy] * mass_of_new_particle; CellP[i].CosmicRayFluxPred[k_CRegy] -= CellP[j].CosmicRayFluxPred[k_CRegy];
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(k=0;k<2;k++)
            {
                CellP[j].CosmicRayAlfvenEnergy[k_CRegy][k] = mass_of_new_particle * CellP[i].CosmicRayAlfvenEnergy[k_CRegy][k]; CellP[i].CosmicRayAlfvenEnergy[k_CRegy][k] -= CellP[j].CosmicRayAlfvenEnergy[k_CRegy][k];
                CellP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k] = mass_of_new_particle * CellP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k]; CellP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k] -= CellP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k];
                CellP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k] = mass_of_new_particle * CellP[i].DtCosmicRayAlfvenEnergy[k_CRegy][k]; CellP[i].DtCosmicRayAlfvenEnergy[k_CRegy][k] -= CellP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k];
            }
#endif
        }
#endif

#if (GALSF_ISMDUSTCHEM_MODEL & 2) && defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
        int l;
        double total_bin_num, total_bin_mass;
        for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {
            for(l=0;l<NUM_ISMDUSTCHEM_SIZE_BINS;l++) {
                total_bin_num = mass_of_new_particle * CellP[i].ISMDustChem_Dust_NumberInBin[k][l]; /* dust grain number conserving */
                total_bin_mass = mass_of_new_particle * get_ISMDustChemEvo_bin_mass(i,k,l, CellP); /* dust grain mass conserving */
                update_ISMDustChemEvo_bin_number_and_slope(j,k,l,total_bin_num,total_bin_mass, CellP); /* update new particle */
                update_ISMDustChemEvo_bin_number_and_slope(i,k,l,total_bin_num,total_bin_mass, CellP); /* update old particle */
            }
        }
#endif
        /* use a better particle shift based on the moment of inertia tensor to place new particles in the direction which is less well-sampled */
#if (NUMDIMS > 1)        
        double norm=0, dp[3]; int m; dp[0]=dp[1]=dp[2]=0;

        // get the eigenvector of NV_T that has the largest eigenvalue (= sparsest sampling direction, if assume equal-mass particles, for our definition of NV_T)
        double nvt[NUMDIMS*NUMDIMS]={0}; for(k=0;k<NUMDIMS;k++) {for(m=0;m<NUMDIMS;m++) {nvt[NUMDIMS*k + m]=CellP[i].NV_T[k][m];}} // auxiliary array to store NV_T in for feeding to GSL eigen routine
        double eigvals_loc[3]={0}, eigvecs_loc[9]={0};
        eigen_symmetric(nvt, eigvals_loc, eigvecs_loc, NUMDIMS);
        int min_eigvec_index = 0; double max_eigval = -MAX_REAL_NUMBER;
        for(k=0;k<NUMDIMS;k++) {if(eigvals_loc[k] > max_eigval) {max_eigval = eigvals_loc[k]; min_eigvec_index=k;}}
        for(k=0;k<NUMDIMS;k++) {dp[k] = eigvecs_loc[NUMDIMS*min_eigvec_index + k];}
        for(k=0;k<NUMDIMS;k++) {norm += dp[k] * dp[k];}
        if(norm > 0)
        {
            norm = 1/sqrt(norm);
            double qq = get_random_number(63432*k + 84*i + 99*j + 358453 + 84537*ThisTask);
            if(qq < 0.5) {norm *= -1.;} // randomly decide which direction along principle axis to orient split (since this is arbitrary, this helps prevent accidental collisions)
            for(k=0;k<NUMDIMS;k++) {dp[k] *= norm;}
            /* project dp onto the isodensity plane (perpendicular to grad_rho): placing daughters on the same
               isodensity surface eliminates the leading O(d . grad_rho) density mismatch between the two daughters,
               which is first-order in d and dominates over the second-order neighbor perturbation the NV_T direction minimizes */
            {
                double grad_rho[3]={0}, grad_rho_norm=0;
                for(k=0;k<NUMDIMS;k++) {grad_rho[k]=CellP[i].Gradients.Density[k]; grad_rho_norm+=grad_rho[k]*grad_rho[k];}
                grad_rho_norm = sqrt(grad_rho_norm);
                if(grad_rho_norm > 0 && CellP[i].Density > 0) {
                    double relative_gradient = grad_rho_norm * P[i].Get_Particle_Size() / CellP[i].Density; // dimensionless: |grad_rho|*dx/rho
                    if(relative_gradient > 0.1) { // non-trivial gradient: project dp onto isodensity plane
                        double dot=0; for(k=0;k<NUMDIMS;k++) {dot += dp[k] * grad_rho[k] / grad_rho_norm;}
                        double dp_perp[3]={0}, perp_norm=0;
                        for(k=0;k<NUMDIMS;k++) {dp_perp[k] = dp[k] - dot*(grad_rho[k]/grad_rho_norm); perp_norm += dp_perp[k]*dp_perp[k];}
                        perp_norm = sqrt(perp_norm);
                        if(perp_norm > 0.1) {for(k=0;k<NUMDIMS;k++) {dp[k] = dp_perp[k] / perp_norm;}} // adopt projected direction; else dp near-parallel to grad_rho, fall through to NV_T direction
                    }
                }
            }
            dx=d_r*dp[0]; dy=d_r*dp[1]; dz=d_r*dp[2];
        }
#endif
        P[i].wakeup = -1; P[j].wakeup = -1; NeedToWakeupParticles_local = 1; /* rather conservative. But we want to update Density and KernelRadius after the particle masses were changed */

    } // closes special operations required only of gas particles

    /* this is allowed to push particles over the 'edges' of periodic boxes, because we will call the box-wrapping routine immediately below.
     but it is important that the periodicity of the box be accounted for in relative positions and that we correct for this before allowing
     any other operations on the particles */
    P[i].Pos[0] += dx; P[j].Pos[0] -= dx; P[i].Pos[1] += dy; P[j].Pos[1] -= dy; P[i].Pos[2] += dz; P[j].Pos[2] -= dz;

#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    P[i].Time_Of_Last_MergeSplit = All.Time; P[j].Time_Of_Last_MergeSplit = All.Time;
#endif
    
    /* Note: New tree construction can be avoided because of  `force_add_element_to_tree()' */
#ifdef PARTICLE_MERGE_SPLIT_EVERY_TIMESTEP    
    long bin = P[i].TimeBin;
    if(FirstInTimeBin[bin] < 0) {FirstInTimeBin[bin]=j; LastInTimeBin[bin]=j; NextInTimeBin[j]=-1; PrevInTimeBin[j]=-1;} /* only particle in this time bin on this task */
    else {NextInTimeBin[j]=FirstInTimeBin[bin]; PrevInTimeBin[j]=-1; PrevInTimeBin[FirstInTimeBin[bin]]=j; FirstInTimeBin[bin]=j;} /* there is already at least one particle; add this one "to the front" of the list */
    force_add_element_to_tree(i, j);
#endif    
    /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
    
#if defined(MHD_CONSERVE_B_ON_REFINEMENT)
    /* flag cells as having just undergone refinement/derefinement for other subroutines to be aware */
    CellP[j].recent_refinement_flag = CellP[i].recent_refinement_flag = 1;
#endif

    return 1; // completed routine successfully
}



/*! Routine to merge particle 'i' into particle 'j'
    The volumetric quantities (density, pressure, gradients, kernel lengths)
    can be re-estimated after the merger operation. but we need to make sure
    all conserved quantities are appropriately dealt with. This also requires some care, to be
    done appropriately, but is a little bit less sensitive and more well-defined compared to
    particle splitting */
int merge_particles_ij(int i, int j)
{
    int k;
    if(P[i].Mass <= 0)
    {
        P[i].Mass = 0; CellP[i].Mass = P[i].Mass;
        return 0;
    }
    if(P[j].Mass <= 0)
    {
        P[j].Mass = 0; CellP[j].Mass = P[j].Mass;
        return 0;
    }
    double mtot = P[j].Mass + P[i].Mass;
    double wt_i = P[i].Mass / mtot;
    double wt_j = P[j].Mass / mtot;
    
    if(((P[i].Type==0)||(P[j].Type==0)) && (P[j].Type!=P[i].Type)) {printf("WARNING: code is trying to merge a gas cell with a non-gas particle. I dont know how to do this. Exiting the merge subroutine."); fflush(stdout); return 0;}

    int swap_ids = 0; if(P[i].Mass > P[j].Mass) {swap_ids = 1;} /* retain the IDs of the more massive progenitor */
#ifdef SINK_WIND_SPAWN
    if(P[i].ID == All.SpawnedWindCellID) {swap_ids = 0;} /* don't copy an agn wind id */
    if(P[j].ID == All.SpawnedWindCellID) {P[j].ID = All.SpawnedWindCellID + 1;} /* offset this to avoid checks through code */
#endif
    if(swap_ids) {P[j].ID=P[i].ID; P[j].ID_child_number=P[i].ID_child_number; P[j].ID_generation=P[i].ID_generation;} /* swap the ids so save the desired set */
    
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    if(P[i].Type==4 && P[j].Type==P[i].Type) /* identify a star-star merger, need to update the effective size -before- updating anything else */
    {
        double m_i=P[i].Mass, m_j=P[j].Mass, dr_i=P[i].StarParticleEffectiveSize*All.cf_atime, dr_j=P[j].StarParticleEffectiveSize*All.cf_atime;
        Vec3<double> dp_star = (P[j].Pos - P[i].Pos) * All.cf_atime; double dr2_ij = dp_star.norm_sq(); // ij position separation (physical units)
        Vec3<double> dv_star = (P[j].Vel - P[i].Vel) / All.cf_atime; double dv2_ij = dv_star.norm_sq(); // ij velocity separation (physical units)
        double phi_prefac=2.01887, phi_tot=-All.G*(0.5*phi_prefac*m_i*m_i/dr_i + 0.5*phi_prefac*m_j*m_j/dr_j + m_i*m_j/sqrt(dr_i*dr_i + dr_j*dr_j + dr2_ij)); // estimate the total gravitational energy of the system here; the '0.5' in phi_prefac accounts for assuming each is virialized, so will sum such that only 1/2 of the potential adds to the total; note phi_prefac is the integral over the kernel functions. its pre-computed here for the cubic spline default choice with this module, but can in general be numerically computed from the kernel functions, though its not so important as long as you are consistent
        double ke_sum = 0.5 * (m_i*m_j/(m_i+m_j)) * dv2_ij; // kinetic energy sum
        double etot_new = DMIN(phi_tot + ke_sum , 0.5 * (-All.G*0.5*phi_prefac*(m_i*m_i/dr_i + m_j*m_j/dr_j))); // updated total energy -- don't let it get too close to zero or positive or this will give a bogus result for the final softening
        double dr_new = -0.5 * All.G*phi_prefac*mtot*mtot / etot_new;  /* solve etot_new = -0.5*All.G*phi_prefac*m_sum*m_sum/dr_new */
        P[j].StarParticleEffectiveSize = dr_new; // particle i will be zero'd out
    }
#endif
    
    // block for merging non-gas particles (much simpler, assume collisionless)
    if((P[i].Type>0)&&(P[j].Type>0))
    {
        Vec3<double> dp = P[j].Pos - P[i].Pos;
        nearest_xyz(dp,-1);
        Vec3<double> p_old_i = P[i].Vel * P[i].Mass;
        Vec3<double> p_old_j = P[j].Vel * P[j].Mass;
        P[j].Pos = P[i].Pos + dp * wt_j; // center-of-mass conserving //
        P[j].Vel = P[j].Vel * wt_j + P[i].Vel * wt_i; // momentum-conserving //
        P[j].GravAccel = P[j].GravAccel * wt_j + P[i].GravAccel * wt_i; // force-conserving //
#ifdef PMGRID
        P[j].GravPM = P[j].GravPM * wt_j + P[i].GravPM * wt_i; // force-conserving //
#endif
        P[j].KernelRadius = pow(pow(P[j].KernelRadius,NUMDIMS)+pow(P[i].KernelRadius,NUMDIMS),1.0/NUMDIMS); // volume-conserving to leading order //
        gizmo_mark_kernel_radius_dirty_indices(&j, 1); /* survivor h changed: mark h-dirty across all caches */
#ifdef METALS
        for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = wt_j*P[j].Metallicity[k] + wt_i*P[i].Metallicity[k];} // metal-mass conserving (includes nuclear species) //
#endif
#ifdef GALSF
        if(P[i].Type==4 && P[j].Type==4) // couple extra potentially-important fields to carry for star particles
        {
#ifdef GALSF_SFR_IMF_SAMPLING
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
            P[j].TimeDistribOfStarFormation = wt_j*P[j].TimeDistribOfStarFormation + wt_i*P[i].TimeDistribOfStarFormation; // average formation time //
            if(P[i].IMF_NumMassiveStars + P[j].IMF_NumMassiveStars > 0) {
                P[j].IMF_WeightedMeanStellarFormationTime = (P[j].IMF_NumMassiveStars*P[j].IMF_WeightedMeanStellarFormationTime + P[i].IMF_NumMassiveStars*P[i].IMF_WeightedMeanStellarFormationTime)/(P[j].IMF_NumMassiveStars + P[i].IMF_NumMassiveStars);}
            else {
                P[j].IMF_WeightedMeanStellarFormationTime = wt_j*P[j].IMF_WeightedMeanStellarFormationTime + wt_i*P[i].IMF_WeightedMeanStellarFormationTime;} // weight by number of massive stars if we still have any since these are still active and needed for stellar evolution, otherwise weight by mass
#endif
            P[j].IMF_NumMassiveStars += P[i].IMF_NumMassiveStars; // O-star number conserving //
#endif
            P[j].StellarAge = wt_j*P[j].StellarAge + wt_i*P[i].StellarAge; // average formation time //
        }
#endif
        /* finally zero out the particle mass so it will be deleted */
        P[i].Mass = 0; CellP[i].Mass = P[i].Mass;
        P[j].Mass = mtot; CellP[j].Mass = P[j].Mass;
        /* momentum shift for passing to tree (so we know how to move it) */
        P[i].dp += P[i].Vel * P[i].Mass - p_old_i;
        P[j].dp += P[j].Vel * P[j].Mass - p_old_j;
        return 1;
    } // closes merger of non-gas particles, only gas particles will see the blocks below //


    // now we have to deal with gas particle mergers //
    double mass_before_merger_i,mass_before_merger_j,volume_before_merger_i,volume_before_merger_j; mass_before_merger_i=P[i].Mass; mass_before_merger_j=P[j].Mass; volume_before_merger_i=P[i].Mass/CellP[i].Density; volume_before_merger_j=P[j].Mass/CellP[j].Density; // save some numbers for potential use below
    if(P[i].TimeBin < P[j].TimeBin)
    {
        P[j].wakeup = -1; NeedToWakeupParticles_local = 1;
    }
    double dm_i=0,dm_j=0,de_i=0,de_j=0,dm_ij,de_ij;
    Vec3<double> dp_i{}, dp_j{}, dp_ij{};
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    dm_i += CellP[i].DtMass;
    dm_j += CellP[j].DtMass;
#endif
    de_i = P[i].Mass * CellP[i].DtInternalEnergy + dm_i*CellP[i].InternalEnergy;
    de_j = P[j].Mass * CellP[j].DtInternalEnergy + dm_j*CellP[j].InternalEnergy;
    dp_i = CellP[i].HydroAccel * P[i].Mass + CellP[i].VelPred * (dm_i / All.cf_atime);
    dp_j = CellP[j].HydroAccel * P[j].Mass + CellP[j].VelPred * (dm_j / All.cf_atime);
    de_i += dot(dp_i, CellP[i].VelPred) / All.cf_atime - 0.5 * dm_i * CellP[i].VelPred.norm_sq() * All.cf_a2inv;
    de_j += dot(dp_j, CellP[j].VelPred) / All.cf_atime - 0.5 * dm_j * CellP[j].VelPred.norm_sq() * All.cf_a2inv;
    dp_ij = dp_i + dp_j;
    dm_ij = dm_i+dm_j;
    de_ij = de_i+de_j;

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    CellP[j].MassTrue += CellP[i].MassTrue;
    CellP[i].MassTrue = 0;
    CellP[j].DtMass=dm_ij;
    CellP[i].DtMass=0;
    CellP[j].dMass = CellP[i].dMass + CellP[j].dMass;
    CellP[i].dMass = 0;
#endif

    /* make sure to update the conserved variables correctly: mass and momentum are easy, energy is non-trivial */
    double egy_old = 0;
    egy_old += mtot * (wt_j*CellP[j].InternalEnergy + wt_i*CellP[i].InternalEnergy); // internal energy //
    /* for periodic boxes, we need to (arbitrarily) pick one position as our coordinate center. we pick i. then everything defined in
        position differences relative to i. the final position will be appropriately box-wrapped after these operations are completed */
    Vec3<double> dp = P[j].Pos - P[i].Pos;
    nearest_xyz(dp,-1);
    Vec3<double> pos_new = P[i].Pos + dp * wt_j;
    Vec3<double> dr_j = (P[i].Pos + dp - pos_new) * All.cf_atime; // displacement of j relative to new pos (physical)
    Vec3<double> dr_i = (P[i].Pos - pos_new) * All.cf_atime;       // displacement of i relative to new pos (physical)

    egy_old += mtot*wt_j * 0.5 * P[j].Vel.norm_sq() * All.cf_a2inv; // kinetic energy (j) //
    egy_old += mtot*wt_i * 0.5 * P[i].Vel.norm_sq() * All.cf_a2inv; // kinetic energy (i) //
    // gravitational energy terms need to be added (including work for moving particles 'together') //
    // Egrav = m*g*h = m * (-grav_acc) * (position relative to zero point) //
    egy_old -= mtot*wt_j * dot(dr_j, P[j].GravAccel) * All.cf_a2inv; // work (j) //
    egy_old -= mtot*wt_i * dot(dr_i, P[i].GravAccel) * All.cf_a2inv; // work (i) //
#ifdef PMGRID
    egy_old -= mtot*wt_j * dot(dr_j, P[j].GravPM) * All.cf_a2inv; // work (j) [PMGRID] //
    egy_old -= mtot*wt_i * dot(dr_i, P[i].GravPM) * All.cf_a2inv; // work (i) [PMGRID] //
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    CellP[j].GravWorkTerm = {}; // since we're accounting for the work above and dont want to accidentally double-count //
#endif


    CellP[j].InternalEnergy = wt_j*CellP[j].InternalEnergy + wt_i*CellP[i].InternalEnergy;
    CellP[j].InternalEnergyPred = wt_j*CellP[j].InternalEnergyPred + wt_i*CellP[i].InternalEnergyPred;
    Vec3<double> p_old_i = P[i].Vel * P[i].Mass;
    Vec3<double> p_old_j = P[j].Vel * P[j].Mass;
    P[j].Pos = pos_new; // center-of-mass conserving //
    P[j].Vel = P[j].Vel * wt_j + P[i].Vel * wt_i; // momentum-conserving //
    CellP[j].VelPred = CellP[j].VelPred * wt_j + CellP[i].VelPred * wt_i; // momentum-conserving //
    P[j].GravAccel = P[j].GravAccel * wt_j + P[i].GravAccel * wt_i; // force-conserving //
#ifdef PMGRID
    P[j].GravPM = P[j].GravPM * wt_j + P[i].GravPM * wt_i; // force-conserving //
#endif
#ifdef MAGNETIC
    // we evolve the conservative variables VB and Vpsi, these should simply add in particle-merge operations //
    { Vec3<double> B_before_split = (CellP[i].B + CellP[j].B) * (1.0 / (volume_before_merger_i + volume_before_merger_j)); /* calculate the real value of B pre-split to know what we need to correctly re-initialize to once the volume partition can be recomputed. here since VB is the conserved variable explicitly integrated, that gets added with the pre-summation volumes (so we take a volume-weighted mean here) */
    CellP[j].BField_prerefinement = B_before_split; CellP[i].BField_prerefinement = B_before_split; } /* record the real value of B pre-split to know what we need to correctly re-initialize to once the volume partition can be recomputed */
    CellP[j].B += CellP[i].B; CellP[j].BPred += CellP[i].BPred; CellP[j].DtB += CellP[i].DtB; /* simply add conserved quantities */
    //CellP[j].DtB = CellP[i].DtB = {}; CellP[j].BPred = CellP[j].B; CellP[i].BPred = CellP[i].B; /* zero the time derivatives and cleaning terms: these will be re-calculated self-consistently from the new mesh configuration */
    CellP[j].divB += CellP[i].divB; // CellP[i].divB = CellP[j].divB = 0; /* this is a conserved quantity (volume-integrated) so can simply be added */
#ifdef DIVBCLEANING_DEDNER
    CellP[j].Phi += CellP[i].Phi; CellP[j].PhiPred += CellP[i].PhiPred; CellP[j].DtPhi += CellP[i].DtPhi; CellP[j].DtB_PhiCorr += CellP[i].DtB_PhiCorr; /* on a merge, can simply add all of these conservative quantities */
    //CellP[i].Phi = CellP[i].PhiPred = CellP[i].DtPhi = CellP[j].Phi = CellP[j].PhiPred = CellP[j].DtPhi = 0; CellP[i].DtB_PhiCorr = CellP[j].DtB_PhiCorr = {}; /* zero the time derivatives and cleaning terms: these will be re-calculated self-consistently from the new mesh configuration */
#endif
#endif

    /* correct our 'guess' for the internal energy with the residual from exact energy conservation */
    double egy_new = mtot * CellP[j].InternalEnergy + mtot * 0.5 * P[j].Vel.norm_sq() * All.cf_a2inv;
    egy_new = (egy_old - egy_new) / mtot; /* this residual needs to be put into the thermal energy */
    if(egy_new < -0.5*CellP[j].InternalEnergy) egy_new = -0.5 * CellP[j].InternalEnergy;
    //CellP[j].InternalEnergy += egy_new; CellP[j].InternalEnergyPred += egy_new;//test during splits
    if(CellP[j].InternalEnergyPred<0.5*CellP[j].InternalEnergy) CellP[j].InternalEnergyPred=0.5*CellP[j].InternalEnergy;


    // now use the conserved variables to correct the derivatives to primitive variables //
    de_ij -= dm_ij * CellP[j].InternalEnergyPred;
    CellP[j].HydroAccel = (dp_ij - CellP[j].VelPred * (dm_ij / All.cf_atime)) * (1.0 / mtot);
    de_ij -= mtot * dot(CellP[j].VelPred, CellP[j].HydroAccel) / All.cf_atime + 0.5 * dm_ij * CellP[j].VelPred.norm_sq() * All.cf_a2inv;
    CellP[j].DtInternalEnergy = de_ij;
    // to be conservative adopt the maximum signal velocity and kernel length //
    CellP[j].MaxSignalVel = sqrt(CellP[j].MaxSignalVel*CellP[j].MaxSignalVel + CellP[i].MaxSignalVel*CellP[i].MaxSignalVel); /* need to be conservative */
    P[j].KernelRadius = pow(pow(P[j].KernelRadius,NUMDIMS)+pow(P[i].KernelRadius,NUMDIMS),1.0/NUMDIMS); /* sum the volume of the two particles */
    gizmo_mark_kernel_radius_dirty_indices(&j, 1); /* survivor h changed: mark h-dirty across all caches */
    CellP[j].ConditionNumber = CellP[j].ConditionNumber + CellP[i].ConditionNumber; /* sum to be conservative */
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    CellP[j].MaxKineticEnergyNgb = DMAX(CellP[j].MaxKineticEnergyNgb,CellP[i].MaxKineticEnergyNgb); /* for the entropy/energy switch condition */
#endif

    // below, we need to take care of additional physics //
#if defined(RADTRANSFER)
    int k_dir;
    for(k=0;k<N_RT_FREQ_BINS;k++)
    {
        CellP[j].ET[k] = wt_j * CellP[j].ET[k] + wt_i * CellP[i].ET[k];
        CellP[j].Rad_E_gamma[k] = CellP[j].Rad_E_gamma[k] + CellP[i].Rad_E_gamma[k]; /* this is a photon number, so its conserved (we simply add) */
#if defined(RT_EVOLVE_ENERGY)
        CellP[j].Rad_E_gamma_Pred[k] = CellP[j].Rad_E_gamma_Pred[k] + CellP[i].Rad_E_gamma_Pred[k];
        CellP[j].Dt_Rad_E_gamma[k] = CellP[j].Dt_Rad_E_gamma[k] + CellP[i].Dt_Rad_E_gamma[k];
#endif
#if defined(RT_EVOLVE_FLUX)
        for(k_dir=0;k_dir<3;k_dir++)
        {
            CellP[j].Rad_Flux[k][k_dir] = CellP[j].Rad_Flux[k][k_dir] + CellP[i].Rad_Flux[k][k_dir];
            CellP[j].Rad_Flux_Pred[k][k_dir] = CellP[j].Rad_Flux_Pred[k][k_dir] + CellP[i].Rad_Flux_Pred[k][k_dir];
            CellP[j].Dt_Rad_Flux[k][k_dir] = CellP[j].Dt_Rad_Flux[k][k_dir] + CellP[i].Dt_Rad_Flux[k][k_dir];
        }
#endif
#ifdef RT_EVOLVE_INTENSITIES
        for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++)
        {
            CellP[j].Rad_Intensity[k][k_dir] = CellP[j].Rad_Intensity[k][k_dir] + CellP[i].Rad_Intensity[k][k_dir];
            CellP[j].Rad_Intensity_Pred[k][k_dir] = CellP[j].Rad_Intensity_Pred[k][k_dir] + CellP[i].Rad_Intensity_Pred[k][k_dir];
            CellP[j].Dt_Rad_Intensity[k][k_dir] = CellP[j].Dt_Rad_Intensity[k][k_dir] + CellP[i].Dt_Rad_Intensity[k][k_dir];
        }
#endif
    }
#endif
#ifdef CHIMES
    double wt_h_i, wt_h_j; // Ratio of hydrogen mass fractions.
#ifdef COOL_METAL_LINES_BY_SPECIES
    wt_h_i = 1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1]);
    wt_h_j = 1.0 - (P[j].Metallicity[0] + P[j].Metallicity[1]);
#else
    wt_h_i = 1.0;
    wt_h_j = 1.0;
#endif
    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) {ChimesGasVars[j].abundances[k] = (ChimesFloat) ((ChimesGasVars[j].abundances[k] * wt_j * wt_h_j) + (ChimesGasVars[i].abundances[k] * wt_i * wt_h_i));}
#endif // CHIMES
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k] = wt_j*P[j].Metallicity[k] + wt_i*P[i].Metallicity[k];} /* metal-mass conserving (includes nuclear species) */
#if defined(GALSF_ISMDUSTCHEM_MODEL)
    for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {CellP[j].ISMDustChem_Dust_Metal[k] = wt_j*CellP[j].ISMDustChem_Dust_Metal[k] + wt_i*CellP[i].ISMDustChem_Dust_Metal[k];} /* dust-mass conserving */
    for(k=0;k<NUM_ISMDUSTCHEM_SOURCES;k++) {CellP[j].ISMDustChem_Dust_Source[k] = wt_j*CellP[j].ISMDustChem_Dust_Source[k] + wt_i*CellP[i].ISMDustChem_Dust_Source[k];} /* dust source-mass conserving */
    for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {CellP[j].ISMDustChem_Dust_Species[k] = wt_j*CellP[j].ISMDustChem_Dust_Species[k] + wt_i*CellP[i].ISMDustChem_Dust_Species[k];} /* dust species-mass conserving */
#if defined(GALSF_ISMDUSTCHEM_GRAINSIZEEVO)
    int l;
    double total_bin_num, total_bin_mass;
    for(k=0;k<NUM_ISMDUSTCHEM_SPECIES;k++) {
        for(l=0;l<NUM_ISMDUSTCHEM_SIZE_BINS;l++) {
            total_bin_num = CellP[j].ISMDustChem_Dust_NumberInBin[k][l] + CellP[i].ISMDustChem_Dust_NumberInBin[k][l]; /* dust grain bin number conserving */
            total_bin_mass = get_ISMDustChemEvo_bin_mass(j,k,l, CellP) + get_ISMDustChemEvo_bin_mass(i,k,l, CellP); /* dust grain bin mass conserving */
            update_ISMDustChemEvo_bin_number_and_slope(j,k,l,total_bin_num,total_bin_mass, CellP);
        }
    }
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
    CellP[j].DtCREgyNewInjectionFromShocks += CellP[i].DtCREgyNewInjectionFromShocks;
#endif
#if defined(SINK_CR_INJECTION_AT_TERMINATION)
    CellP[j].Sink_CR_Energy_Available_For_Injection += CellP[i].Sink_CR_Energy_Available_For_Injection;
#endif
    int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        CellP[j].CosmicRayEnergy[k_CRegy] += CellP[i].CosmicRayEnergy[k_CRegy];
        CellP[j].CosmicRayEnergyPred[k_CRegy] += CellP[i].CosmicRayEnergyPred[k_CRegy];
        CellP[j].DtCosmicRayEnergy[k_CRegy] += CellP[i].DtCosmicRayEnergy[k_CRegy];
#if defined(CRFLUID_EVOLVE_SPECTRUM)
        CellP[j].CosmicRay_Number_in_Bin[k_CRegy] += CellP[i].CosmicRay_Number_in_Bin[k_CRegy];
        CellP[j].DtCosmicRay_Number_in_Bin[k_CRegy] += CellP[i].DtCosmicRay_Number_in_Bin[k_CRegy];
#endif
        CellP[j].CosmicRayFlux[k_CRegy] += CellP[i].CosmicRayFlux[k_CRegy];
        CellP[j].CosmicRayFluxPred[k_CRegy] += CellP[i].CosmicRayFluxPred[k_CRegy];
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
        int kAlfDir; for(kAlfDir=0; kAlfDir<2; kAlfDir++) {
            CellP[j].CosmicRayAlfvenEnergy[k_CRegy][kAlfDir] += CellP[i].CosmicRayAlfvenEnergy[k_CRegy][kAlfDir];
            CellP[j].CosmicRayAlfvenEnergyPred[k_CRegy][kAlfDir] += CellP[i].CosmicRayAlfvenEnergyPred[k_CRegy][kAlfDir];
            CellP[j].DtCosmicRayAlfvenEnergy[k_CRegy][kAlfDir] += CellP[i].DtCosmicRayAlfvenEnergy[k_CRegy][kAlfDir];
        }
#endif
    }
#endif

#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    P[i].Time_Of_Last_MergeSplit = All.Time; P[j].Time_Of_Last_MergeSplit = All.Time;
#endif
    
    /* finally zero out the particle mass so it will be deleted */
    P[i].Mass = 0; CellP[i].Mass = P[i].Mass;
    P[j].Mass = mtot; CellP[j].Mass = P[j].Mass;
    /* momentum shift for passing to tree (so we know how to move it) */
    P[i].dp += P[i].Vel * P[i].Mass - p_old_i;
    P[j].dp += P[j].Vel * P[j].Mass - p_old_j;
    /* call the pressure routine to re-calculate pressure (and sound speeds) as needed */
    set_eos_pressure(j, P, CellP);
#if defined(MHD_CONSERVE_B_ON_REFINEMENT)
    /* flag cells as having just undergone refinement/derefinement for other subroutines to be aware */
    CellP[j].recent_refinement_flag = CellP[i].recent_refinement_flag = 1;
#endif
    return 1;
}


/*!
  This routine swaps the location of two pointers/indices (either to a node or to a particle) in the treewalk needed for neighbor searches and gravity.
  This should be run if you are messing around with the indices of things but don't intend to do a whole treebuild after. - MYG
 */

void swap_treewalk_pointers(int i, int j){
    // walk the tree and any time we see a nextnode or sibling set to i, swap it to j and vice versa
    int no, next, pre_sibling_i=-1, pre_sibling_j=-1, previous_node_i, previous_node_j;
    no = All.MaxPart;

    while(no >= 0){ // walk the whole tree, starting from the root node (=All.MaxPart)
        if(no < All.MaxPart) { // we got a particle
            next=Nextnode[no];
            if(no != i && no != j){ // don't mess with Nextnodes if looking at i or j - handle that separately
                if(next == i) {Nextnode[no] = j; previous_node_i = no;}
                else if(next == j) { Nextnode[no] = i; previous_node_j = no;}
            }
            no = next;
        } else if(no < All.MaxPart + MaxNodes + MaxForeignNodes)  { // we have a local or foreign tree node (Phase 9 LET)
            next = Nodes[no].u.d.nextnode;
            if(next == i) { previous_node_i = no; Nodes[no].u.d.nextnode = j;}
            else if(next == j) { previous_node_j = no; Nodes[no].u.d.nextnode = i;}
            if(Nodes[no].u.d.sibling == i) {Nodes[no].u.d.sibling = j; pre_sibling_i = no;}
            else if(Nodes[no].u.d.sibling == j) { Nodes[no].u.d.sibling = i; pre_sibling_j = no;}
            no = next;
        } else { // pseudoparticle (Phase 9: index shifted up by MaxForeignNodes)
            int pseudo_idx = no - MaxNodes - MaxForeignNodes;
            next = Nextnode[pseudo_idx];
            if(next==i) {Nextnode[pseudo_idx] = j;}
            else if(next == j) {Nextnode[pseudo_idx] = i;}
            no = next;
        }
    }

    if(Nextnode[i] == j){ // handle case if i->j
        Nextnode[i] = Nextnode[j];
        Nextnode[j] = i;
    } else if (Nextnode[j] == i) { // if j->i
        Nextnode[j] = Nextnode[i];
        Nextnode[i] = j;
    } else { // neither i->j nor j->i
        no = Nextnode[i];
        Nextnode[i] = Nextnode[j];
        Nextnode[j] = no;
    }
    no = Father[i];
    Father[i] = Father[j];
    Father[j] = no;
}


/*!
  This routine deletes a particle from the linked list for the treewalk, preserving the lists's integrity. This must be run if you are deleting particles but don't want to do a while treebuild after. - MYG
*/
void remove_particle_from_treewalk(int i){
    int no, next;
    no = All.MaxPart;
    while(no >= 0){ // walk the tree to find anything that might point to i and redirect it to i's nextnode
        if(no < All.MaxPart){
            next = Nextnode[no];
            if(Nextnode[no] == i) {Nextnode[no] = Nextnode[i];}
        } else if (no < All.MaxPart+MaxNodes){
            next = Nodes[no].u.d.nextnode;
            if(next == i) {Nodes[no].u.d.nextnode = Nextnode[i];}
            if(Nodes[no].u.d.sibling == i) {Nodes[no].u.d.sibling = Nextnode[i];}
        } else {
            next = Nextnode[no - MaxNodes];
            if(next == i) {Nextnode[no - MaxNodes] = Nextnode[i];}
        }
        no = next;
    }
}

/*! This is an important routine used throughout -- any time particle masses are variable OR particles can
    be created/destroyed: it goes through the particle list, makes sure they are in the appropriate order (gas
    must all come before collisionless particles, though the collisionless particles can be blocked into any order
    we like), swaps particles as needed to restore the correct ordering, eliminates particles from the main list
    if they have zero or negative mass (i.e. this does the actual 'deletion' operation), and then reconstructs the
    list of particles in each timestep (if particles had to be re-ordered) so that the code will not crash looking for
    the 'wrong' particles (or non-existent particles). In general, if you do any operations that involve particle
    creation, this needs to be called before anything is 'done' with those particles. If you do anything involving particle
    'deletion', the standard procedure should be to set the deleted particle mass to zero, and then let this routine
    (when it is called in standard sequence) do its job and 'clean up' the particle
 */
void rearrange_particle_sequence(void)
{
    int i, j, flag = 0, flag_sum, j_next;
    int count_elim, count_gaselim, count_sink_elim, tot_elim, tot_gaselim, tot_sink_elim;
    struct particle_data psave;
    struct gas_cell_data gascellsave;
#ifdef CHIMES
    struct gasVariables gasVarsSave;
#endif

    int do_loop_check = 0;
    if(Gas_split>0)
    {
        N_gas += Gas_split;
        NumPart += Gas_split;
        Gas_split = 0;
        do_loop_check = 1;
    }
#ifdef GALSF
    if(Stars_converted)
    {
        N_gas -= Stars_converted;
        Stars_converted = 0;
        do_loop_check = 1;
    }
#endif
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
    if(Grains_promoted)
    {
        N_gas += Grains_promoted;
        Grains_promoted = 0;
        do_loop_check = 1;
    }
#endif
    if(NumPart <= N_gas) {do_loop_check=0;}
    if(N_gas <= 0) {do_loop_check=0;}

    /* if more gas than stars, need to be sure the block ordering is correct (gas first, then stars) */
    if(do_loop_check)
    {
        j_next = N_gas;
        for(i = 0; i < N_gas; i++) /* loop over the gas block */
            if(P[i].Type != 0) /* and look for a particle converted to non-gas */
            {
                /* ok found a non-gas particle: */
                for(j = j_next; j < NumPart; j++) /* loop from N_gas to Numpart, to find first labeled as gas */
                    if(P[j].Type == 0) break; /* break on that to record the j of interest */
                if(j >= NumPart) endrun(181170); /* if that j is too large, exit with error */
                
                if(j < NumPart-1) {j_next = j + 1;} else {j_next = NumPart - 1;}
                psave = P[i]; /* otherwise, save the old pointer */
                P[i] = P[j]; /* now set the pointer equal to this new P[j] */
                P[j] = psave; /* now set the P[j] equal to the old, saved pointer */
                /* so we've swapped the two P[i] and P[j] */
                gascellsave = CellP[i];
                CellP[i] = CellP[j];
                CellP[j] = gascellsave;  /* have the gas particle take its gas/fluid cell pointer with it */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i,j);
#endif
#ifdef CHIMES /* swap chimes-specific 'gasvars' structure which is separate from the default code gas cell structure */
                gasVarsSave = ChimesGasVars[i]; ChimesGasVars[i] = ChimesGasVars[j]; ChimesGasVars[j] = gasVarsSave;
                /* Old particle (now at position j) is no longer a gas particle, so delete its abundance array. */
                free_gas_abundances_memory(&(ChimesGasVars[j]), &ChimesGlobalVars);
                ChimesGasVars[j].abundances = NULL; ChimesGasVars[j].isotropic_photon_density = NULL; ChimesGasVars[j].G0_parameter = NULL; ChimesGasVars[j].H2_dissocJ = NULL;
#endif /* CHIMES */

                /* ok we've now swapped the ordering so the gas particle is still inside the block */
                flag = 1;
            }
    }

    count_elim = 0;
    count_gaselim = 0;
    count_sink_elim = 0;
    /* loop over entire block looking for things with zero mass, which need to be eliminated */
    for(i = 0; i < NumPart; i++)
        if(P[i].Mass <= 0 || !isfinite(P[i].Mass))
        {
            P[i].Mass = 0; CellP[i].Mass = P[i].Mass;
            TimeBinCount[P[i].TimeBin]--;
            if(TimeBinActive[P[i].TimeBin]) {NumForceUpdate--;}

            if(P[i].Type == 0)
            {
                TimeBinCountGas[P[i].TimeBin]--;

                P[i] = P[N_gas - 1];
                CellP[i] = CellP[N_gas - 1];
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i, N_gas-1);
#endif
                /* swap with properties of last gas particle (i-- below will force a check of this so its ok) */
#ifdef CHIMES
                free_gas_abundances_memory(&(ChimesGasVars[i]), &ChimesGlobalVars);
                ChimesGasVars[i] = ChimesGasVars[N_gas - 1];
                ChimesGasVars[N_gas - 1].abundances = NULL; ChimesGasVars[N_gas - 1].isotropic_photon_density = NULL; ChimesGasVars[N_gas - 1].G0_parameter = NULL; ChimesGasVars[N_gas - 1].H2_dissocJ = NULL;
#endif

                P[N_gas - 1] = P[NumPart - 1]; /* redirect the final gas pointer to go to the final particle (BH) */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(N_gas - 1, NumPart-1);
                remove_particle_from_treewalk(NumPart - 1);
#endif
                N_gas--; /* shorten the total N_gas count */
                count_gaselim++; /* record that a BH was eliminated */
            }
            else
            {
                if(P[i].Type == 5) {count_sink_elim++;} /* record elimination if BH */
                P[i] = P[NumPart - 1]; /* re-directs pointer for this particle to pointer at final particle -- so we
                                        swap the two; note that ordering -does not- matter among the non-fluid/gas cells
                                        so its fine if this mixes up the list ordering of different particle types */
#ifdef MAINTAIN_TREE_IN_REARRANGE
                swap_treewalk_pointers(i, NumPart - 1);
                remove_particle_from_treewalk(NumPart - 1);
#endif
            }

            NumPart--;
            i--;
            count_elim++;
        }

    MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_sink_elim, &tot_sink_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(count_elim) {flag = 1;}

    if(ThisTask == 0) {if(tot_elim > 0) {printf("Rearrange: Eliminated %d/%d gas/star particles and merged away %d sink particles.\n", tot_gaselim, tot_elim - tot_gaselim - tot_sink_elim, tot_sink_elim);}}

    All.TotNumPart -= tot_elim;
    All.TotN_gas -= tot_gaselim;
#ifdef SINK_PARTICLES
    All.TotSinks -= tot_sink_elim;
#endif

    MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(flag_sum) {reconstruct_timebins();}
}



#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
/* check if a star-star particle pair is eligible to be merged within this model context, when representing essentially degenerate kinematic information (sub-softening) */
double evaluate_starstar_merger_for_starcluster_particle_pair(int i, int j)
{
    if(evaluate_starstar_merger_for_starcluster_eligibility(j)) // already evaluated i, make sure j is also eligible
    {
        double eta_position = 1., eta_velocity2 = 2.; // variables for below, defined for convenience here
        
        // consider separation and relative velocities
        int k; Vec3<double> dp = P[j].Pos - P[i].Pos; // calculate separation
        nearest_xyz(dp,-1); // correct for box appropriately
        double r2 = dp.norm_sq(); // squared position difference
        
        double eps_i = ForceSoftening_KernelRadius(i), eps_j = ForceSoftening_KernelRadius(j), eps_ij = DMAX(eps_i,eps_j);
        double threshold_separation = eta_position * eps_ij; // threshold separation to consider
        if(r2 > threshold_separation*threshold_separation) {return -1;} // only allow if sufficiently close
        
        Vec3<double> dv = P[j].Vel - P[i].Vel; // calculate velocity difference
        NGB_SHEARBOX_BOUNDARY_VELCORR_(P[i].Pos,P[j].Pos,dv,-1); // correct for box appropriately
        double v2 = dv.norm_sq(); // squared velocity difference
        v2 *= All.cf_a2inv; // physical
        
        double threshold_velocity2 = eta_velocity2 * All.G*(P[i].Mass+P[j].Mass)/(All.cf_atime*sqrt(r2)); // threshold 'escape' velocity to consider
        if(v2 > threshold_velocity2) {return -1;} // only allow if relative velocity is sufficiently low (bound)
        
        return r2; // we'll use this as our threshold criterion (considering the closest pair which meets all our criteria)
    }
    return -1; // default to no merger allowed
}

/* check if a star is -potentially- eligible to enter this subroutine for looking for a merger pair companion, which can only occur if it is sufficiently old and in a very dense region, otherwise it will not be allowed */
int evaluate_starstar_merger_for_starcluster_eligibility(int i)
{
    if(All.Time <= All.TimeBegin) {return 0;} // don't allow on first timestep
    if(P[i].Type != 4) {return 0;} // only stars
    if(P[i].Mass*UNIT_MASS_IN_SOLAR > 1.e6) {return 0;} // stop merging beyond a certain point, only want to downgrade the resolution so much
    double dt_gyr = evaluate_stellar_age_Gyr(i), dt_threshold_gyr = 0.05;
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
    dt_gyr -= P[i].TimeDistribOfStarFormation*UNIT_TIME_IN_GYR; // want to be sure we're this far past the last possible star that formed
#endif
    if(dt_gyr < dt_threshold_gyr) {return 0;} // sufficiently old (don't want to do this for extremely young stars as messes up feedback and early dynamics)
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION // need to figure out if the new version of this makes sense
    double r_NGB = 1.25 * pow((All.DesNumNgb*All.G*P[i].Mass)/P[i].tidal_tensor_mag_prev , 1./3.); // kernel size enclosing some target neighbor number in a constant-density medium
    if(r_NGB > 0.5*ForceSoftening_KernelRadius(i)) {return 0;} // sufficiently dense region (need to have effective nearest-neighbor spacing approaching the minimum softening, with some arbitrary threshold we set)
#elif defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
    double h_i=ForceSoftening_KernelRadius(i), tidal_mag=0., fac_self=-P[i].Mass*kernel_gravity(0.,1.,1.,1)/(h_i*h_i*h_i); // get what's needed for tidal tensor computation
    tidal_mag -= P[i].tidal_tensorps.trace(); // want the (negative) trace, actually, and in general this -shouldn't- include the self-contribution
    if(tidal_mag > 0) {
        double ngb_dist = 1.25 * pow( (All.DesNumNgb * All.G * P[i].Mass / tidal_mag) , 1./3. ); // distance to the N'th nearest-neighbor
        if(ngb_dist > h_i) {return 0;} // sufficiently dense region (need to have effective nearest-neighbor spacing approaching the minimum softening, with some arbitrary threshold we set)
    }
#else
    if(get_particle_timestep_in_physical(i)*UNIT_TIME_IN_GYR*1000. > 0.01) {return 0;} // if the particle is taking a very long timestep, it's probably in a very low-density region, so don't allow it to merge (this is a crude proxy for the local density, but it's very fast to check and works well in practice) 
#endif
    return 1; // allow this particle to -consider- the possibility of a merger
}
#endif



#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
/* subroutine to check if too little time has passed since the last merge-split, in which case we won't allow it again */
int check_if_sufficient_mergesplit_time_has_passed(int i)
{
    double N_timesteps_fac = 300.; // require > N timesteps before next merge/split, default was 100, but can be more aggressive - something between 10-100 works well in practice [definitely shorter than 10 can cause problems]
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 2) || defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT)
    N_timesteps_fac = 30.;
#endif
    if(P[i].Time_Of_Last_MergeSplit <= All.TimeBegin) {N_timesteps_fac *= 10. * get_random_number(832LL*i + 890345645LL + 83457LL*ThisTask + 12313403LL*P[i].ID);} // spread initial timing out over a broader range so it doesn't all happen at once after the startup
    double dtime_code = All.Time - P[i].Time_Of_Last_MergeSplit; // time [in code units] since last merge/split
    double dt_incodescale = (get_particle_timestep_in_physical(i) * All.cf_hubble_a) * All.cf_atime; // timestep converted appropriately to code units [physical if non-comoving, else scale factor]
    if(dtime_code < N_timesteps_fac*dt_incodescale) {return 0;} // not enough time passed, prohibit
#if !defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(All.ComovingIntegrationOn) {if(dtime_code < 1.e-8) {return 0;}} // also enforce an absolute time limit (don't use for nuclear zooms since absolute timesteps can become arbitrarily short
#endif
    return 1; // otherwise, if no check so far to reject, allow this merge/split
}
#endif
