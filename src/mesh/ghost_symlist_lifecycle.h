/* ghost_symlist_lifecycle.h — ghost-exchange + symmetric-CSR lifecycle helpers.
 *
 * These helpers absorb the ghost_exchange / symlist setup-and-teardown that
 * used to be copy-pasted at every call site of density(), hydro_gradient_calc(),
 * and hydro_force() (in accel.cc, init.cc, run.cc).  The three top-level
 * functions now include the lifecycle inline, so any caller gets the correct
 * sequencing automatically — you cannot forget to prep or clean up.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef GHOST_SYMLIST_LIFECYCLE_H
#define GHOST_SYMLIST_LIFECYCLE_H

#include "neighbor_list.h"
#include "sfc_tiles.h"
#include "ghost_exchange_spec.h"
#include "ghost_writeback.h"

extern void gpu_build_symmetric_neighbor_list(struct particle_data *P, int num_total,
    int *active_indices, int num_active, neighbor_list_t *out,
    double search_radius_factor);

/* Forward declaration of the out-of-line host All accessor (defined in
 * core/predict.cc). Retained as a stylistic intent-tag for host-snapshot
 * reads; bare All.* in host code is equally correct. */
#ifdef __cplusplus
extern "C" {
#endif
struct global_data_all_processes *gizmo_host_all_ptr(void);
integertime gizmo_host_ti_current(void);
#ifdef __cplusplus
}
#endif

/* Shared search/ghost-inflation factor (only grows when TURB_DIFF_DYNAMIC
   widens the dynamic-diffusion kernel). */
static inline double gizmo_ghost_safety_factor(void)
{
    double f = 1.0;
#ifdef TURB_DIFF_DYNAMIC
    f = DMAX(f, gizmo_host_all_ptr()->TurbDynamicDiffFac);
#endif
    return f;
}

/* Prologue for density(): drift all particles to current time and import
   ghost particles from neighbouring ranks.  Unconditional — both underlying
   routines early-out for NTask==1.

   NB: this is the ALL-TYPES variant. Many callers (sinks, AGS gravity,
   stellar feedback, DM dispersion, fuzzy DM, grain physics, etc.) need
   cross-type ghosts and route here. Pure-hydro callers (density iteration,
   gradient prep/refresh) should call gizmo_hydro_prep_ghosts / use
   ghost_exchange_hydro directly to avoid DM/star kernel pollution of tile
   hmax that inflates hydro search radii by 10-100x. */
static inline void gizmo_density_prep_ghosts(double safety)
{
    double t0 = my_second();
    move_particles(gizmo_host_ti_current());
    double t_drift = timediff(t0, my_second());
    double t1 = my_second();
    ghost_exchange(safety);
    double t_ghost = timediff(t1, my_second());
    cpu_charge_child(CPU_DENSMISC, t_drift);
    cpu_charge_child(CPU_GHOSTIMPORT, t_ghost);
}

/* gizmo_request_filtered_ghost_import — drift + caller-list ghost import.
 *
 * Performs a full global drift via move_particles(All.Ti_Current), then runs
 * ghost_exchange_run with the caller's exact active-source list and search
 * semantics. Result: ghost particles for the caller-named loop appended at
 * P[NumPart..NumPart+N_ghost], ready for kernel walks that read imported
 * ghost state.
 *
 * IMPORTANT NAMING: this routine performs full global drift AND ghost import.
 * "request_filtered" means the GHOST POOL is filtered to the caller's active
 * source list — it does NOT mean the routine is request-driven in the
 * Mode B P2P sense. Callers needing ghostless P2P request/reply should NOT
 * call this; they go through run_neighbor_loop<Spec> which selects the
 * Mode B path that does no global drift and no ghost import.
 *
 * Runner-migrated loops (those whose dispatch flows through
 * mesh/neighbor_loop_runner::run_neighbor_loop<Spec>) MUST NOT call this
 * directly. The runner decides whether to call it based on the chosen path
 * (only paths that need imported ghosts use it; Mode B paths skip it
 * entirely, by design, to preserve the lazy-drift / no-global-mutation
 * tiny-N corridor).
 *
 * Legacy non-runner callers may continue to use this API until each is
 * migrated to the runner. New callers should NOT be added — express the
 * loop as a Spec for the runner.
 *
 * Most loops are symmetric; only use NGB_SEARCH_ONEWAY when the actual
 * kernel predicate is strictly r_ij < h_i (density-like calls).
 */
/* radius_policy + j_radius_scale are REQUIRED args (no defaults) so every
 * caller threads the SSOT supply-side reach contract explicitly.  Runner
 * Mode A passes Spec::radius_policy + nlr_spec_symmetric_j_radius_scale<Spec>();
 * legacy non-runner callers explicitly pass MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES
 * + 1.0 to preserve their pre-policy behavior byte-for-byte.  Compiler enforces
 * explicit-thread at every call site (no silent LEGACY fallback default). */
static inline void gizmo_request_filtered_ghost_import(const char *caller_name,
                                                       int search_mode,
                                                       unsigned int supply_type_mask,
                                                       const int *active_indices,
                                                       int num_active,
                                                       const double *active_radii,
                                                       double safety,
                                                       mode_b_radius_policy_t radius_policy,
                                                       double j_radius_scale,
                                                       int supply_band_dominated)
{
    double t0 = my_second();
    move_particles(gizmo_host_ti_current());
    double t_drift = timediff(t0, my_second());

    int alloc_n = (num_active > 0 ? num_active : 1);
    double (*qpos)[3] = (double (*)[3]) malloc((size_t)alloc_n * sizeof(double[3]));
    double *qh = (double *) malloc((size_t)alloc_n * sizeof(double));
    for(int a = 0; a < num_active; a++) {
        int i = active_indices[a];
        qpos[a][0] = P[i].Pos[0];
        qpos[a][1] = P[i].Pos[1];
        qpos[a][2] = P[i].Pos[2];
        qh[a] = active_radii[a];
    }

    struct ghost_exchange_spec_t sp = {
        0u,
        supply_type_mask,
        search_mode,
        safety,
        caller_name,
        num_active,
        (const double (*)[3]) qpos,
        (const double *) qh,
        radius_policy,
        j_radius_scale,
        supply_band_dominated
    };

    double t1 = my_second();
    ghost_exchange_run(&sp);
    double t_ghost = timediff(t1, my_second());
    free(qh);
    free(qpos);
    cpu_charge_child(CPU_DENSMISC, t_drift);
    cpu_charge_child(CPU_GHOSTIMPORT, t_ghost);
}

/* gizmo_request_filtered_ghost_import_fresh — collective-safe variant.
 *
 * Same semantics as gizmo_request_filtered_ghost_import, with a preceding
 * unconditional ghost_exchange_cleanup so every rank enters the import in
 * the same state. A local ghost-count guard is illegal here: exporter-only
 * ranks can have zero received ghosts while importer ranks have nonzero
 * ghosts, causing only some ranks to enter the collective exchange.
 *
 * Same caller restriction as the base routine: loop/Spec/consumer code whose
 * dispatch flows through run_neighbor_loop<Spec> MUST NOT call this directly —
 * the runner is the sanctioned caller and decides, per the chosen path, whether
 * an imported-ghost prep is needed (Mode B paths skip it entirely). The runner's
 * Mode-A prep (neighbor_loop_runner.cc) calling this IS the intended path, not a
 * violation of the restriction.
 */
static inline int gizmo_request_filtered_ghost_import_fresh(const char *caller_name,
                                                           int search_mode,
                                                           unsigned int supply_type_mask,
                                                           const int *active_indices,
                                                           int num_active,
                                                           const double *active_radii,
                                                           double safety,
                                                           mode_b_radius_policy_t radius_policy,
                                                           double j_radius_scale,
                                                           int supply_band_dominated)
{
    ghost_exchange_cleanup();
    gizmo_request_filtered_ghost_import(caller_name, search_mode, supply_type_mask,
                                        active_indices, num_active, active_radii, safety,
                                        radius_policy, j_radius_scale, supply_band_dominated);
    return 1;
}

/* Hydro-density variant of gizmo_density_prep_ghosts: gas-only pool, gas-only
   active gate, ONE-WAY criterion (r_ij < h_i — correct for density; symmetric
   max(h_i, h_j) is wrong here and over-imports by orders of magnitude, per
   direct waste-ratio measurement). Use from hydro density iteration only. */
static inline void gizmo_hydro_density_prep_ghosts(double safety)
{
    double t0 = my_second();
    move_particles(gizmo_host_ti_current());
    double t_drift = timediff(t0, my_second());
    double t1 = my_second();
    ghost_exchange_hydro_oneway(safety);
    double t_ghost = timediff(t1, my_second());
    cpu_charge_child(CPU_DENSMISC, t_drift);
    cpu_charge_child(CPU_GHOSTIMPORT, t_ghost);
}

/* Fresh broad downstream-handoff import after the runner-based density()
   path — IMPORT ONLY, NO DRIFT. The iterative runner's
   Mode A imports an EXACT-QUERY ghost pool sized to current iter actives+
   radii — narrower than the legacy broad hydro-oneway envelope, and
   unsuitable as the downstream pool consumed by cellcorrections /
   gradients / hydro_force. Legacy density() handled this implicitly
   because its top-of-routine prep_ghosts was already a broad pool that
   persisted through finalize. The runner path cleans its exact-query pool at
   runner-return; this helper then rebuilds the broad downstream pool from
   converged radii.
   Cleanup-first is safe even if no pool exists (ghost_exchange_cleanup
   early-returns on NumPart_before_ghost < 0).

   The name avoids "prep_ghosts" because that historically implies the
   drift+import combo (gizmo_hydro_density_prep_ghosts calls
   move_particles(All.Ti_Current) followed by ghost_exchange_hydro_oneway).
   Including the drift here is WRONG for the post-density handoff: by the
   time density() returns, All.Ti_Current is the current step's source time
   for the about-to-run hydro/cellcorrections drift. A second
   move_particles() at the same Ti_Current is nominally a no-op but breaks
   lazy-drift semantics in subtle ways. Concretely it caused SP4
   drift_particle() abort 'no prediction into past allowed' with
   time1 = 0.5 * dt_step (an earlier-than-Ti_current target) when used
   post-density. */
static inline void gizmo_hydro_density_import_ghosts_fresh_no_drift(double safety)
{
    ghost_exchange_cleanup();
    double t0 = my_second();
    ghost_exchange_hydro_oneway(safety);
    cpu_charge_child(CPU_GHOSTIMPORT, timediff(t0, my_second()));
}

/* Corridor/subcycle-INTERNAL shared-topology build (no direct consumer
   callers: gizmo_hydro_corridor_begin owns the per-step call, the corridor
   refresh full path rebuilds through the epilogue below, and
   transport_subcycle_prepare_topology rebuilds after particle mutations):
   allocate the per-step active-index array, refresh ghost CellP with
   converged density/h, and build the symmetric CSR neighbor list used by
   gradients and hydro_force.

   Multi-rank tiny-N gate: when no rank has any active gas particle, the
   ghost refresh + symlist build are pure overhead (~1-2s/step on Vista
   2-rank). MPI_Allreduce the local active-gas count once; if global == 0,
   skip the collective ghost_exchange and the kernel-side symlist build,
   leaving an empty symlist that the gradient/hydro kernels fast-path on.

   is_mode_b: when the corridor runs Mode B (request-driven, tiny-N), the
   gradient/hydro walkers use ONLY the active-index list below; they never
   read imported ghosts (Mode B walks local particles + P2P request/reply)
   and nothing reads the shared CSR (the only CSR consumers are Mode-A
   external_csr and transport_subcycle, and TRANSPORT_SUBCYCLE forces Mode A).
   So on Mode B the ghost exchange + CSR build are pure dead work — skipped
   here; only the active-index list is built. */
static inline void gizmo_gradients_prep_symlist(double safety, double search_fac, bool is_mode_b)
{
    /* active-index allocation (persists across gradients + hydro; used by the
       request-driven Mode B walkers as their active list) */
    gizmo_sym_num_active = 0;
    for(int ii : ActiveParticleList) {if(P[ii].Type == 0 && P[ii].Mass > 0) gizmo_sym_num_active++;}
    gizmo_sym_active_indices = (int *) mymalloc("sym_active",
        (gizmo_sym_num_active > 0 ? gizmo_sym_num_active : 1) * sizeof(int));
    {int aa = 0; for(int ii : ActiveParticleList) {
        if(P[ii].Type == 0 && P[ii].Mass > 0) gizmo_sym_active_indices[aa++] = ii;}}

    /* Global active-gas count: cached for refresh's matching skip decision.
       Single int Allreduce is microseconds; not measurable next to the
       O(0.1-1s) ghost_exchange we may skip. */
    if(NTask > 1) {
        MPI_Allreduce(&gizmo_sym_num_active, &gizmo_sym_num_active_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    } else {
        gizmo_sym_num_active_global = gizmo_sym_num_active;
    }

    if(is_mode_b) {
        /* Mode B imports no ghosts and reads no shared CSR — both would be dead
           work. But it MUST still enter with no live ghost extension in NumPart:
           dematerialize any ghosts left materialized by an earlier phase (safe
           no-op when none). Skip ONLY the fresh ghost import + the CSR build.
           Other per-step code (conservation sums, timestep, next-step setup)
           iterates NumPart, so leftover ghosts here would corrupt it. */
        ghost_exchange_cleanup();
        return;
    }

    /* When global active-gas count is 0, skip the collective ghost refresh
       (otherwise pure overhead — ~1-2s/step on 2-rank tiny-N). The symlist
       build still runs with num_active=0; it allocates the size-1 backstop
       offsets/neighbors that downstream gradient_evaluate_gpu requires (a
       NULL backing buffer there causes a segfault even though pairs=0). */
    if(gizmo_sym_num_active_global > 0 && NTask > 1) {
        double t_refresh = my_second();
        ghost_exchange_cleanup();
        ghost_exchange_hydro(safety);
        if(ThisTask == 0) {PRINT_STATUS("Ghost refresh before gradients (%.4f s)", timediff(t_refresh, my_second()));}
    }

    /* build the symmetric CSR list (max(h_i,h_j) search radius) */
    double t_sym = my_second();
    gpu_build_symmetric_neighbor_list(P, NumPart, gizmo_sym_active_indices, gizmo_sym_num_active, &gizmo_sym_neighbor_list, search_fac);
    if(ThisTask == 0) {PRINT_STATUS("Symmetric neighbor list: %d active, %lld pairs (%.4f s)",
                                    gizmo_sym_num_active, (long long)gizmo_sym_neighbor_list.total_pairs, timediff(t_sym, my_second()));}
}

/* Corridor-INTERNAL full ghost re-import + CSR rebuild (no direct consumer
   callers: gizmo_hydro_corridor_refresh_ghost_values owns the call as its
   fail-closed full path and republishes the view afterwards). Refreshes
   ghosts so downstream pairs see updated owner fields on both sides, and
   rebuilds the CSR (ghost indices may shift after re-exchange). Skipped
   when global active-gas was 0 in prep — no work in either direction. The
   size-1 backstop symlist from prep stays valid for the (no-op) kernels. */
static inline void gizmo_gradients_refresh_symlist(double safety, double search_fac, bool is_mode_b)
{
    if(gizmo_sym_num_active_global == 0) return;
    if(is_mode_b) {
        /* Mode B built no CSR + no ghosts in prep; nothing to refresh. */
        return;
    }
    if(NTask > 1) {
        double t_refresh = my_second();
        free_neighbor_list(&gizmo_sym_neighbor_list);
        ghost_exchange_cleanup();
        ghost_exchange_hydro(safety);
        gpu_build_symmetric_neighbor_list(P, NumPart, gizmo_sym_active_indices, gizmo_sym_num_active, &gizmo_sym_neighbor_list, search_fac);
        if(ThisTask == 0) {PRINT_STATUS("Ghost refresh + CSR rebuild after gradients: %lld pairs (%.4f s)",
                                        (long long)gizmo_sym_neighbor_list.total_pairs, timediff(t_refresh, my_second()));}
    }
}

/* Epilogue for hydro_force(): free the symlist and remove ghost particles.
   Unconditional — including under TRANSPORT_SUBCYCLE. (The old scheme that
   kept this pool alive for the RT subcycle loop was broken by construction:
   sink accretion, wind spawning, rearrange_particle_sequence, and RT source
   injection all run between hydro_force and the subcycle loop, and their own
   ghost imports/cleanups and particle-array mutations corrupted the retained
   pool and the CSR ghost indices into it. The subcycle now builds its own
   topology immediately before its loop — transport_subcycle_prepare_topology
   in core/transport_subcycle.cc — after all of those have run.) */
static inline void gizmo_hydro_cleanup_symlist_and_ghosts(void)
{
    gizmo_sym_neighbor_list_free();
    ghost_exchange_cleanup();
}

/* init.cc only: density() was called for initial h-convergence, no subsequent
   gradients/hydro in this context, so just tear down ghosts. */
static inline void gizmo_density_init_cleanup(void)
{
    if(NTask > 1) {ghost_exchange_cleanup();}
}

#endif /* GHOST_SYMLIST_LIFECYCLE_H */
