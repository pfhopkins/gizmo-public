#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

/* this file handles the FIRE short-range radiation-pressure and
    photo-ionization terms. written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

/* radiation_pressure_winds_consolidated lives in galaxy_sf/radfb_rp_loop.cc
 * (alongside RadFBRPSpec) — the toplevel runner-template caller belongs in
 * the same two-file module as the Spec. This file stays HOST-ONLY
 * (STARFORM_OBJS) to keep HII_heating_singledomain + do_the_local_ionization
 * out of GPU_OBJS — the HII physics is serial-greedy singledomain and
 * doesn't belong in a GPU TU. accel.cc:203 calls
 * `radiation_pressure_winds_consolidated()` via core/proto.h's forward
 * decl; the linker resolves it from radfb_rp_loop.o. */

#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
#include <vector>
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../mesh/mode_b_local_walker.h"
#include "../mesh/neighbor_loop_runner.h"
#include "../system/gpu_particles_arena.h"
#endif

/* radiation_pressure_winds_consolidated moved → radfb_rp_loop.cc (see
 * header comment block above). */





#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
/* Routines for simple FIRE local photo-ionization heating feedback model. This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */


/* NOTE: HII_heating_singledomain outer loop cannot be wrapped in Kokkos::parallel_for
 * because the inner greedy search calls ngb_treefind_variable_targeted, which writes
 * to the global Ngblist buffer.  Making Ngblist thread_local in allvars.cc would
 * unblock this; deferred until that change is made. */

/* SINGLEDOMAIN SEMANTIC — DO NOT add a ghost_writeback for this routine, and do
 * NOT instrument it with ghost_write_detector_begin/end (the detector would
 * abort because no writeback is registered).
 *
 * By design (matching the legacy CPU path), HII heating only ionizes gas
 * particles owned by THIS rank — never ghosts. It accepts a small accuracy loss
 * at domain boundaries in exchange for speed: cross-domain ionization is
 * intentionally outside this routine's semantics. The pre-walk filter at
 * `if(j_cand >= local_count) continue;` enforces local-only.
 *
 * THEREFORE THIS ROUTINE MUST NOT IMPORT GHOSTS. The GPU candidate builder
 * (hii_gpu_path) searches the LOCAL gas pool only (num_all = local_count). A
 * prior version imported a broad all-types ghost pool via gizmo_density_prep_ghosts
 * and then discarded every ghost here — pure wasted work that dominated the FIRE
 * "HII" wall (the all_types ghost-discovery elephant). Removed 2026-06-29.
 * Radiation-pressure feedback is a SEPARATE loop (radfb_rp_loop.cc) that
 * legitimately uses ghosts/writeback — do not conflate the two.
 *
 * If you find yourself thinking "shouldn't we add ghost_writeback_hii?" — NO. */

/* ============================================================================
 * Tiny-N surgical refactor (2026-05-17): HII_heating_singledomain now dispatches
 * between two candidate-list providers based on global source count, then funnels
 * every source through one greedy-ionization helper (SSOT for the per-source
 * physics — RHII expansion, per-iter radius filter, do_the_local_ionization).
 *
 * Why this is NOT a Spec port onto the runner template:
 *   - The per-source body is serial-greedy: later sources observe DelayTimeHII /
 *     InternalEnergy set by earlier sources within the same step. The runner
 *     template's parallel-for kernel shape cannot preserve this ordering.
 *   - Singledomain semantic: ghosts are filtered out post-NL build. No ghost
 *     writeback. The runner contract assumes per-pair work into accumulators +
 *     writeback — wrong shape here.
 *
 * Why two providers:
 *   - On tiny-N steps with 1–2 ionizing sources, the GPU-NL path's unconditional
 *     gizmo_density_prep_ghosts + gpu_particles_arena_acquire(num_all,...) +
 *     gpu_ngb_list_build pay O(Ntot) for O(num_active) work — violates the
 *     "never touch globals on tiny-N path" directive.
 *   - hii_local_path uses mode_b_local_neighbor_walk (tree walk on the existing
 *     Nodes[]) per source: zero arena acquire, zero ghost import.
 *   - hii_gpu_path is the existing GPU-NL machinery, unchanged in this commit.
 *     Threshold dispatch picks per-step based on gizmo_nlr_modeb_threshold_sum_for
 *     ("hii", default 256), overridable from the parameterfile via
 *     NeighborLoopModeBThresholdSum.
 *
 * Tree-validity invariant for hii_local_path (verified 2026-05-17): every step
 * driver branch in core/run.cc:185-192 produces a valid Nodes[] for this step
 * (domain_decomp / domain_decomp_treerebuild / force_update_tree). Then
 * core/accel.cc:75-82 runs density()/ags_density()/force_update_hmax() before
 * compute_stellar_feedback at accel.cc:102 → HII_heating_singledomain. No
 * merge_split occurs between tree-build and HII. The local walker reads
 * Nodes[] which is fresh + h-refreshed at the HII call site. ========================================================================== */

/* Per-source state staged in the pre-pass and consumed by both candidate-list
 * providers + the greedy helper. Separates "which sources fire this step" from
 * "how do we collect each source's candidate gas list" from "the greedy
 * ionization physics." */
struct HIISourcePrep {
    int    i;          /* source particle index */
    double dt;         /* per-source physical timestep */
    double stellum;    /* effective ionizing luminosity (chimes-aware) */
    double R_for_NL;   /* per-source max search radius (1.26·RHIIMAX, h-clamped) */
};

/* Output counters. Filled by the pre-pass (census fields) + the greedy helper
 * (ionization fields). Top-level MPI_Reduce reads this. */
struct HIIStats {
    double total_N_ionizing_part;
    double total_Ndot_ionizing;
    double total_m_ionized;
    double total_N_ionized;
    double avg_RHII;
};

/* SSOT for per-source greedy ionization. Walks the passed-in candidate slice,
 * runs the RHII-expansion do/while, filters by radius each iteration, applies
 * do_the_local_ionization to accepted j. Identical physics under both providers
 * — only `candidates[]` differs. */
static void hii_greedy_ionize_source(const HIISourcePrep& src,
                                     const int *candidates, int num_cand,
                                     double uion,
                                     std::vector<int>& ngb_list_touse,
                                     HIIStats& s)
{
    const int MAX_N_ITERATIONS_HIIFB = 5;
    int    i        = src.i;
    double dt       = src.dt;
    double stellum  = src.stellum;
    Vec3<double> pos = P[i].Pos;
    double rho      = P[i].DensityAroundParticle;
    double h_i      = P[i].KernelRadius;

    double RHII = 4.78e-9 * pow(stellum, 0.333) * pow(rho*All.cf_a3inv*UNIT_DENSITY_IN_CGS, -0.66667);
    RHII /= All.cf_atime * UNIT_LENGTH_IN_CGS;
    double RHIIMAX = 2. * 240.0 * pow(stellum, 0.5) / (All.cf_atime*UNIT_LENGTH_IN_CGS);
    if(RHIIMAX < 2.0*h_i)  RHIIMAX = 2.0*h_i;
    if(RHIIMAX > 10.0*h_i) RHIIMAX = 10.0*h_i;
    double mionizable = VOLUME_NORM_COEFF_FOR_NDIMS * rho * RHII*RHII*RHII;
    double M_ionizing_emitted = (3.05e10 * PROTONMASS_CGS) * stellum * (dt * UNIT_TIME_IN_CGS);
    mionizable = DMIN(mionizable, M_ionizing_emitted/UNIT_MASS_IN_CGS);
    if(RHII > RHIIMAX)   RHII = RHIIMAX;
    if(RHII < 0.3*h_i)   RHII = 0.3*h_i;
    double RHII_initial = RHII;
    double prandom = get_random_number(P[i].ID + 7);

    /* Defense-in-depth gate (already enforced in the pre-pass; preserved here
     * to match current per-source semantics exactly). */
    if(prandom >= 5.0*mionizable/P[i].Mass) return;

    double mionized = 0.0, mion_actual = 0.0;
    int    jnearest = -1; double rnearest = MAX_REAL_NUMBER;
    int    NITER_HIIFB = 0;

    /* Per-iter filtered candidate buffer (caller-owned scratch, reused across
     * sources within one provider call to avoid per-source allocation churn on
     * the GPU path with many HII sources). */
    if((int)ngb_list_touse.capacity() < num_cand) ngb_list_touse.reserve(num_cand > 0 ? num_cand : 1);

    int more_iters = 1;
    do {
        double RHII_2 = RHII*RHII;
        jnearest = -1; rnearest = MAX_REAL_NUMBER;
        double R_search = RHII;
        if(h_i > 0.5*R_search) R_search = 0.5*h_i;
        double R_search_2 = R_search * R_search;

        /* Per-iter radius filter over the provider's candidate slice.
         * Provider guarantees: candidates are local real particle indices
         * (j < num_local; ghosts already filtered upstream). Gas type, positive
         * mass, and per-iter radius are RECHECKED here as defense-in-depth —
         * do not delete these checks even if some provider seems to imply them. */
        int numngb = 0;
        ngb_list_touse.resize(num_cand > 0 ? num_cand : 0);
        for(int nn = 0; nn < num_cand; nn++) {
            int j_cand = candidates[nn];
            if(P[j_cand].Type != 0 || P[j_cand].Mass <= 0) continue;
            Vec3<double> dpc = pos - P[j_cand].Pos; nearest_xyz(dpc, 1);
            if(dpc.norm_sq() > R_search_2) continue;
            ngb_list_touse[numngb++] = j_cand;
        }

        if(numngb > 0) {
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
            qsort(ngb_list_touse.data(), numngb, sizeof(int), compare_densities_for_sort);
#endif
            for(int n = 0; n < numngb; n++) {
                if(mionized >= mionizable) break;
                int j = ngb_list_touse[n];
                if(P[j].Mass <= 0 || P[j].Type != 0) continue;
#ifdef HYDRO_MULTIFLUID_DM
                if(P[j].FluidType == FLUID_DM) continue; /* skip dark-fluid neighbors */
#endif
                if(CellP[j].DelayTimeHII > 0) continue;
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) && !defined(CHIMES_HII_REGIONS)
                if(CellP[j].Ne > 0.8) continue;
#endif
                Vec3<double> dr = pos - P[j].Pos; nearest_xyz(dr, 1);
                double r2 = dr.norm_sq();
                if(r2 > RHII_2) continue;
                double r = sqrt(r2), u = 0;
                int already_ionized = 0;
                if(CellP[j].InternalEnergy < CellP[j].InternalEnergyPred) u = CellP[j].InternalEnergy; else u = CellP[j].InternalEnergyPred;
                if(CellP[j].DelayTimeHII > 0) already_ionized = 1;
#if !defined(CHIMES_HII_REGIONS)
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
                if((CellP[j].Ne > 0.8) || (u > 50.*uion)) already_ionized = 1;
#else
                if(u > uion) already_ionized = 1;
#endif
#endif
                if(already_ionized) continue;
                int do_ionize = 0; double prob = 0;
                if((r <= RHII) && (mionized < mionizable)) {
                    double m_effective = P[j].Mass*(CellP[j].Density/rho);
                    double m_available = mionizable - mionized;
                    if(m_effective <= m_available) { do_ionize = 1; prob = 1.001; }
                    else { prob = m_available/m_effective; if(prandom < prob) do_ionize = 1; }
                    if(do_ionize == 1) {
                        already_ionized = do_the_local_ionization(j, dt, i);
                        s.total_N_ionized += 1;
                        mion_actual += P[j].Mass;
                        s.avg_RHII += P[j].Mass*r*All.cf_atime*UNIT_LENGTH_IN_KPC;
                    }
                    mionized += prob*m_effective;
                }
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
                if((CellP[j].Density < rnearest) && (already_ionized == 0)) { rnearest = CellP[j].Density; jnearest = j; }
#else
                if((r < rnearest) && (already_ionized == 0)) { rnearest = r; jnearest = j; }
#endif
            }
        }

        /* Backstop: ionize jnearest if still photons available. */
        if((mionized < mionizable) && (jnearest >= 0)) {
            int j = jnearest;
            double m_effective = P[j].Mass*(CellP[j].Density/rho);
            double m_available = mionizable - mionized;
            double prob = m_available/m_effective;
            int do_ionize = 0;
            if(prandom < prob) do_ionize = 1;
            if(do_ionize == 1) {
                (void)do_the_local_ionization(j, dt, i);
                s.total_N_ionized += 1;
                mion_actual += P[j].Mass;
                Vec3<double> dr = pos - P[j].Pos; nearest_xyz(dr, 1); double r2 = dr.norm_sq();
                s.avg_RHII += P[j].Mass*sqrt(r2)*All.cf_atime*UNIT_LENGTH_IN_KPC;
            }
            mionized += prob*m_effective;
        }

        /* RHII expansion (mirrors legacy logic exactly) */
        double RHIImultiplier = 1.10;
        more_iters = 0;
        if(mionized < 0.95*mionizable) {
            if((RHII >= DMAX(30.0*RHII_initial, RHIIMAX)) || (NITER_HIIFB >= MAX_N_ITERATIONS_HIIFB)) {
                mionized = 1.001*mionizable;
            } else {
                if(mionized <= 0) RHIImultiplier = 2.0;
                else {
                    RHIImultiplier = pow(mionized/mionizable, -0.333);
                    if(RHIImultiplier > 5.0)  RHIImultiplier = 5.0;
                    if(RHIImultiplier < 1.26) RHIImultiplier = 1.26;
                }
                RHII *= RHIImultiplier;
                if(RHII > 1.26*RHIIMAX) RHII = 1.26*RHIIMAX;
                more_iters = 1;
            }
        }
        NITER_HIIFB++;
    } while(more_iters && mionized < mionizable);

    if(mion_actual > 0) s.total_m_ionized += mion_actual;
}

/* Large-N path: existing GPU NL machinery (ghost-import + arena acquire +
 * gpu_ngb_list_build over num_all). Structural change funnels the per-source
 * body through hii_greedy_ionize_source; behavior otherwise unchanged.
 * Redundant-ghost-import drop deferred pending periodic-wrap candidate-counter
 * pre-validation. */
static void hii_gpu_path(const std::vector<HIISourcePrep>& src,
                         double uion,
                         HIIStats& s)
{
    int num_src = (int)src.size();

    /* SINGLEDOMAIN: no ghost import (see file-top NOTE).  The GPU candidate
     * builder searches LOCAL gas only; any ghosts present from an earlier loop
     * lie at [local_count, NumPart) and are excluded by num_all = local_count
     * (and re-checked by the j_cand >= local_count filter below).  No MPI
     * collective here, so no cross-rank deadlock risk.  Particles are already
     * drifted to the current time by density()/force_update before this call
     * (core/accel.cc), so no re-drift is needed (matches hii_local_path). */
    int local_count = ghost_get_num_local();
    int num_all = local_count;   /* singledomain: LOCAL pool only.  ghost_get_num_local()
                                  * already returns NumPart when no ghost pool is present;
                                  * never fall back to a ghost-inclusive count (would
                                  * re-admit ghosts, violating the local-only invariant). */

    /* Build flat index + radius arrays for the GPU NL builder. */
    std::vector<int>    src_idx_flat;     src_idx_flat.reserve(num_src);
    std::vector<double> src_radii_flat;   src_radii_flat.reserve(num_src);
    /* hii_fb sources are non-gas (Type 4/5/...) searching the gas-only cached SIDX
     * (gpu_step_sidx_ptr): pass their current positions explicitly, as the radii are. */
    std::vector<double> src_pos_flat;     src_pos_flat.reserve((size_t)num_src * 3);
    for(int aa = 0; aa < num_src; aa++) {
        src_idx_flat.push_back(src[aa].i);
        src_radii_flat.push_back(src[aa].R_for_NL);
        src_pos_flat.push_back((double)P[src[aa].i].Pos[0]);
        src_pos_flat.push_back((double)P[src[aa].i].Pos[1]);
        src_pos_flat.push_back((double)P[src[aa].i].Pos[2]);
    }

    gpu_neighbor_list_t gnl = {};
    std::vector<int> gnl_neighbors_host;
    if(num_src > 0) {
        gpu_particles_arena_acquire(num_all, P, CellP);
        struct particle_data *P_gpu = gpu_particles_arena_P();
        gpu_ngb_list_build(P_gpu, num_all,
                           src_idx_flat.data(), num_src,
                           NGB_SEARCH_ONEWAY, 1 /* gas only */,
                           &gnl, gpu_step_sidx_ptr(), 1.0, src_radii_flat.data(), src_pos_flat.data(), "hii_fb");
        /* gnl.neighbors lives in DEVICE_SPACE (CudaSpace). The per-source loop
         * below indexes it from host code; deep_copy once to a host buffer.
         * (Host memcpy from CudaSpace segfaults on GH200.) */
        if(gnl.total_pairs > 0) {
            gnl_neighbors_host.resize((size_t)gnl.total_pairs);
            gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
        }
    }
    const int *gnl_neighbors = gnl_neighbors_host.empty() ? NULL : gnl_neighbors_host.data();

    /* Sequential per-source loop — preserves greedy ionization ordering.
     * Two caller-owned scratch buffers (cand_buf for the post-ghost-filter
     * candidate slice, ngb_scratch for the helper's per-iter radius-filtered
     * list); both grow geometrically on first push past capacity and are
     * reused across all sources in this call. */
    std::vector<int> cand_buf;
    std::vector<int> ngb_scratch;
    for(int aa = 0; aa < num_src; aa++) {
        /* offsets are int64_t (CSR row pointers); per-source candidate count
         * is bounded by num_total < 2^31 so the local int nl_n is safe. */
        int64_t nl_start = gnl.offsets[aa], nl_end = gnl.offsets[aa+1];
        int nl_n = (int)(nl_end - nl_start);
        /* Filter ghosts (singledomain). The helper rechecks Type/Mass/radius. */
        cand_buf.clear();
        if((int)cand_buf.capacity() < nl_n) cand_buf.reserve(nl_n > 0 ? nl_n : 1);
        for(int64_t nn = nl_start; nn < nl_end; nn++) {
            int j_cand = gnl_neighbors[nn];
            if(j_cand >= local_count) continue; /* skip ghosts */
            cand_buf.push_back(j_cand);
        }
        hii_greedy_ionize_source(src[aa], cand_buf.data(), (int)cand_buf.size(),
                                 uion, ngb_scratch, s);
    }

    if(num_src > 0) {
        gpu_ngb_list_free(&gnl);
        gpu_particles_arena_invalidate();
    }
    /* No ghost import (singledomain) -> nothing to clean up here. */
}

/* Tiny-N path: per-source tree walk on the existing Nodes[] (valid + h-refreshed
 * at this call site — see invariant note above). Touches NO globals beyond the
 * tree itself: no gizmo_density_prep_ghosts, no gpu_particles_arena_acquire,
 * no gpu_ngb_list_build. The walk's candidates are drifted and reduced to the
 * exact set by mode_b_drift_and_filter_candidates, so the helper reads drifted P[j]. */
static void hii_local_path(const std::vector<HIISourcePrep>& src,
                           double uion,
                           HIIStats& s)
{
    std::vector<int> candidates;
    std::vector<int> ngb_scratch;
    for(size_t a = 0; a < src.size(); a++) {
        int    i = src[a].i;
        double R = src[a].R_for_NL;
        double pos_arr[3] = { P[i].Pos[0], P[i].Pos[1], P[i].Pos[2] };

        candidates.clear();
        mode_b_local_neighbor_walk(pos_arr, R,
                                   1u << 0,                 /* gas-only type mask */
                                   MODE_B_SEARCH_ONEWAY,
                                   MODE_B_RADIUS_DEFAULT,
                                   candidates);
        mode_b_drift_and_filter_candidates(pos_arr, R, 1u << 0, MODE_B_SEARCH_ONEWAY,
                                           MODE_B_RADIUS_DEFAULT, 1.0, candidates);
        hii_greedy_ionize_source(src[a], candidates.data(), (int)candidates.size(),
                                 uion, ngb_scratch, s);
    }
}

void HII_heating_singledomain(void)    /* this version of the HII routine only communicates with particles on the same processor */
{
#ifdef RT_CHEM_PHOTOION
    return; // the work here is done in the actual RT routines if this switch is enabled //
#endif
    if(All.HIIRegion_fLum_Coupled<=0) {return;}
    if(All.Time<=0) {return;}
    PRINT_STATUS("Local HII-Region photo-heating/ionization calculation");
    double uion = HIIRegion_Temp / (MEAN_MOLECULAR_WEIGHT_IONIZED * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS); /* assume fully-ionized gas with gamma=5/3 */

    HIIStats stats = {};
    std::vector<HIISourcePrep> src;
    src.reserve(64);

    /* Pre-pass: identify candidate sources, accumulate census diagnostics for
     * ALL ionizing candidates (preserves total_N_ionizing_part semantics as a
     * source-population census, not a stochastic sample count), apply the
     * deterministic prandom<5*mionizable/Mass filter, and stage per-source
     * (dt, stellum, R_for_NL) for the provider + helper. */
    for(int ip : ActiveParticleList) {
#ifdef SINK_HII_HEATING
        if(!((P[ip].Type == 5) || is_galsf_stellar_candidate_type(P[ip].Type, All.ComovingIntegrationOn))) continue;
#else
        if(!is_galsf_stellar_candidate_type(P[ip].Type, All.ComovingIntegrationOn)) continue;
#endif
        if(P[ip].Mass <= 0 || !isfinite(P[ip].Mass)) continue;
        double dt_i = get_particle_feedback_timestep_in_physical(ip, P);
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
        if(P[ip].Type == 5) dt_i = P[ip].dt_since_last_gas_search;
#endif
        if(dt_i <= 0) continue;
        double stellum_i = All.HIIRegion_fLum_Coupled * particle_ionizing_luminosity_in_cgs(ip);
#ifdef CHIMES_HII_REGIONS
        stellum_i = chimes_ion_luminosity(evaluate_stellar_age_Gyr(ip)*1000., P[ip].Mass*UNIT_MASS_IN_SOLAR) * 4.68e-11;
#endif
        if(stellum_i <= 0) continue;
        if(P[ip].KernelRadius <= 0 || P[ip].DensityAroundParticle <= 0) continue;
        double h_local = P[ip].KernelRadius;

        /* Census diagnostics for ALL ionizing candidates (before stochastic gate). */
        stats.total_N_ionizing_part += 1;
        stats.total_Ndot_ionizing   += stellum_i * (3.05e10/HYDROGEN_MASSFRAC);

        /* Deterministic prandom prefilter — same formula as the per-source helper.
         * Sources that fail this contribute ZERO work in the helper (helper has
         * the same gate as defense-in-depth). Filtering here keeps the NL build
         * tight; get_random_number(ID+7) is deterministic so identical prandom
         * is seen in the helper. */
        {
            double rho_i = P[ip].DensityAroundParticle;
            double RHII_i = 4.78e-9*pow(stellum_i,0.333)*pow(rho_i*All.cf_a3inv*UNIT_DENSITY_IN_CGS,-0.66667) / (All.cf_atime*UNIT_LENGTH_IN_CGS);
            double RHIIMAX_i = 2.*240.0*pow(stellum_i,0.5)/(All.cf_atime*UNIT_LENGTH_IN_CGS);
            if(RHIIMAX_i < 2.0*h_local)  RHIIMAX_i = 2.0*h_local;
            if(RHIIMAX_i > 10.0*h_local) RHIIMAX_i = 10.0*h_local;
            if(RHII_i > RHIIMAX_i) RHII_i = RHIIMAX_i;
            if(RHII_i < 0.3*h_local) RHII_i = 0.3*h_local;
            double mion_i = VOLUME_NORM_COEFF_FOR_NDIMS*rho_i*RHII_i*RHII_i*RHII_i;
            double M_emit = (3.05e10*PROTONMASS_CGS)*stellum_i*(dt_i*UNIT_TIME_IN_CGS)/UNIT_MASS_IN_CGS;
            mion_i = DMIN(mion_i, M_emit);
            double prandom_i = get_random_number(P[ip].ID + 7);
            if(prandom_i >= 5.0*mion_i/P[ip].Mass) continue;
        }

        /* Per-source max search radius (covers full RHII expansion to 1.26·RHIIMAX). */
        double RHIIMAX_l = 2. * 240.0*pow(stellum_i, 0.5) / (All.cf_atime * UNIT_LENGTH_IN_CGS);
        if(RHIIMAX_l < 2.0*h_local)  RHIIMAX_l = 2.0*h_local;
        if(RHIIMAX_l > 10.0*h_local) RHIIMAX_l = 10.0*h_local;
        double R_for_NL = 1.26 * RHIIMAX_l;
        if(R_for_NL < 0.5*h_local) R_for_NL = 0.5*h_local;

        HIISourcePrep pp; pp.i = ip; pp.dt = dt_i; pp.stellum = stellum_i; pp.R_for_NL = R_for_NL;
        src.push_back(pp);
    }

    /* Collective dispatch: all ranks must agree on which path runs, since the
     * GPU path's ghost_exchange is collective. */
    int num_src = (int)src.size();
    int global_num_src = 0;
    MPI_Allreduce(&num_src, &global_num_src, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(global_num_src > 0) {
        int threshold = gizmo_nlr_modeb_threshold_sum_for("hii", /*spec_default=*/256);
        if(global_num_src <= threshold) {
            hii_local_path(src, uion, stats);
        } else {
            hii_gpu_path(src, uion, stats);
        }
    }

    double totMPI_N_ionizing_part=0,totMPI_Ndot_ionizing=0,totMPI_m_ionized=0,totMPI_avg_RHII=0,totMPI_N_ionized=0;
    MPI_Reduce(&stats.total_N_ionizing_part, &totMPI_N_ionizing_part, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.total_Ndot_ionizing,   &totMPI_Ndot_ionizing,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.total_m_ionized,       &totMPI_m_ionized,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.total_N_ionized,       &totMPI_N_ionized,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.avg_RHII,              &totMPI_avg_RHII,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(totMPI_N_ionizing_part>0)
        {
            if(totMPI_m_ionized>0) {totMPI_avg_RHII /= totMPI_m_ionized;}
            PRINT_STATUS(" ..Nsources=%g with dN/dt=%g/s ionized N=%g (M=%g sol) cells in <R_HII>=%g kpc",totMPI_N_ionizing_part,totMPI_Ndot_ionizing,totMPI_N_ionized,totMPI_m_ionized*UNIT_MASS_IN_SOLAR,totMPI_avg_RHII);
            fprintf(FdHIIHeating,"%.16g %g %g %g %g %g \n",All.Time,totMPI_N_ionizing_part,totMPI_Ndot_ionizing,totMPI_N_ionized,totMPI_m_ionized*UNIT_MASS_IN_SOLAR,totMPI_avg_RHII); fflush(FdHIIHeating);
        }
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {fflush(FdHIIHeating);}
    } // ThisTask == 0
    CPU_Step[CPU_HIIHEATING] += measure_time();
} // void HII_heating_singledomain(void)



int do_the_local_ionization(int target, double dt, int source)
{
#if defined(CHIMES_HII_REGIONS) // set a number of chimes-specific quantities here //
    int k,age_bin=0; double stellar_age_myr=1000.*evaluate_stellar_age_Gyr(source), log_age_Myr=log10(stellar_age_myr); // determine stellar age bin
    if(log_age_Myr<CHIMES_LOCAL_UV_AGE_LOW) {age_bin=0;} else if(log_age_Myr < CHIMES_LOCAL_UV_AGE_MID) {age_bin = (int) floor(((log_age_Myr - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW) + 1);}
    else {age_bin = (int) floor((((log_age_Myr - CHIMES_LOCAL_UV_AGE_MID) / CHIMES_LOCAL_UV_DELTA_AGE_HI) + ((CHIMES_LOCAL_UV_AGE_MID - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW)) + 1); if(age_bin > CHIMES_LOCAL_UV_NBINS - 1) {age_bin = CHIMES_LOCAL_UV_NBINS - 1;}}
    // reset all of the HII-region chimes quantities to null
    for(k=0;k<CHIMES_LOCAL_UV_NBINS;k++) {CellP[target].Chimes_fluxPhotIon_HII[k]=0; CellP[target].Chimes_G0_HII[k]=0;}
    // set the quantities desired for this age bin specifically: need a softened radius, for use here //
    double stellar_mass=P[source].Mass*UNIT_MASS_IN_SOLAR; Vec3<double> dp = P[source].Pos - P[target].Pos;
    nearest_xyz(dp); dp *= All.cf_atime*UNIT_LENGTH_IN_CGS; double r2 = dp.norm_sq(); // separation in cgs
    double eps_cgs=KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER*ForceSoftening_KernelRadius(source)*All.cf_atime*UNIT_LENGTH_IN_CGS; // plummer equivalent softening
    r2+=eps_cgs*eps_cgs; // gravitational Softening (cgs units)
    CellP[target].Chimes_fluxPhotIon_HII[age_bin] = (1.0 - All.Chimes_f_esc_ion) * chimes_ion_luminosity(stellar_age_myr, stellar_mass) / r2; // cgs flux of H-ionising photons per second seen by the star particle
    CellP[target].Chimes_G0_HII[age_bin] = (1.0 - All.Chimes_f_esc_G0) * chimes_G0_luminosity(stellar_age_myr, stellar_mass) / r2; // cgs flux in the 6-13.6 eV band

#else

#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) 
    CellP[target].InternalEnergy = DMAX(CellP[target].InternalEnergy , HIIRegion_Temp / (MEAN_MOLECULAR_WEIGHT_IONIZED * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS)); /* assume fully-ionized gas with gamma=5/3 */
    CellP[target].InternalEnergyPred = CellP[target].InternalEnergy; /* full reset of the internal energy */
#else
    double delta_U_of_ionization = (20.-13.6) * ((ELECTRONVOLT_IN_ERGS / PROTONMASS_CGS) / UNIT_SPECEGY_IN_CGS) * (1.-DMAX(0.,DMIN(1.,CellP[target].Ne/1.5))); /* energy injected per unit mass, in code units, by ionization, assuming each atom absorbs, and mean energy of absorbed photons is given by x=18 eV here (-13.6 for energy of ionization) */
    double Theat_star = 1.38 * 3.2, Z_sol = 1.; // typical IMF-averaged temp of ionizing star=32,000 K, with effective ionization temperature parameter psi=1.38 (temp of ionized e's in units of stellar temp). Then metallicity in solar units.
#ifdef METALS
    Z_sol = P[target].Metallicity[0]/All.SolarAbundances[0]; // set metallicity
#if defined(GALSF_ISMDUSTCHEM_MODEL) 
    Z_sol = (P[target].Metallicity[4]-CellP[target].ISMDustChem_Dust_Metal[4])/All.SolarAbundances[4]; // cooling for HII equilibrium depends predominantly on gas-phase O so need to account for any locked up in dust
#endif
#endif
    double t4_eqm = DMIN( 1.5*Theat_star , 3.85/DMAX(log(DMAX(390.*Z_sol/Theat_star,1.001)),0.01) ); // equilibrium H2 region temperature in 1e4 K: 1.5*Theat = eqm temp for pure-H region, while second expression assumes eqm cooling with O+C, etc, but breaks down at low-Z when metals don't dominate cooling.
    double u_eqm = (t4_eqm * HIIRegion_Temp) / (MEAN_MOLECULAR_WEIGHT_IONIZED * (GAMMA_DEFAULT-1.) * U_TO_TEMP_UNITS); // converted to specific internal energy, assuming full ionization
    double u_post_ion_no_cooling = CellP[target].InternalEnergy + delta_U_of_ionization; // energy after ionization, before any cooling
    double u_final = DMIN( u_post_ion_no_cooling , u_eqm ); // don't heat to higher temperature than intial energy of ionization allows
    CellP[target].InternalEnergy = u_final; CellP[target].InternalEnergyPred = u_final; /* add it */
    /* assign typical strong HII region flux + enough flux to maintain cell fully-ionized, regardless (x'safety-factor'): note this is repeated in the 'self-shield' routine but has to be in both because otherwise this wont appear on the first timestep after which a gas cell is flagged as ionized by this subroutine */
    double n1000 = CellP[target].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / 1000.; // density in 1000 cm^-3
    double flux_compactHII = DMAX(0.85*pow(n1000,1./3.) , 1) * 2.6e5*n1000; // set to typical value in HII region or minimum needed to maintain f_neutral < 1e-5-ish, whichever is larger
    CellP[target].Rad_Flux_UV += flux_compactHII; CellP[target].Rad_Flux_EUV += flux_compactHII;
#endif
    CellP[target].Ne = 1.0 + 2.0*yhelium(target, P); /* set the cell to fully ionized */

#endif

    CellP[target].DelayTimeHII = DMIN(dt, 10./UNIT_TIME_IN_MYR); /* tell the code to flag this in the cooling subroutine */
    return 1;
}

#endif
