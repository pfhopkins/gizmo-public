/* hydro/hydro_corridor.cc — hydro corridor mode-decision + lifecycle.
 *
 * See hydro/hydro_corridor.h for the API contract and the corridor
 * sequencing rationale.
 *
 * Written by Philip F. Hopkins (phopkins@caltech.edu) for GIZMO. */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/neighbor_list.h"            /* gizmo_sym_* globals */
#include "../mesh/neighbor_loop_runner.h"     /* nlr_external_csr */
#include "../mesh/ghost_symlist_lifecycle.h"  /* gizmo_gradients_prep_symlist / refresh_symlist */
#include "hydro_corridor.h"

/* Default adaptive thresholds (global summed / per-rank-max active gas) used
 * when the optional NeighborLoopModeBThreshold{Sum,Max} parameters are unset.
 * The corridor uses the same DECISION FORM as the runner's per-loop policy:
 *   Mode B iff sum > 0 && sum <= threshold_sum && max <= threshold_max.
 * It does NOT use the same NUMBERS. The corridor's default is 1000/1000 on
 * active GAS; a loop that reaches the runner's own adaptive path without a
 * dispatch override falls back to the Spec's threshold, defaulting to 64/64
 * on that loop's own active count (neighbor_loop_runner.cc, nlr_spec_threshold_*).
 * Consequence, by design and measured: between 65 and 1000 active gas the
 * corridor selects Mode B for gradients/hydro_force/cellcorrections (which do
 * take the override) while density -- which does not set dispatch_override --
 * independently selects Mode A and performs its own one-way ghost import.
 * That import feeds density's own neighbour search; it is not handed downstream.
 * Centralized here so a future tuning lives in one place. */
static constexpr int CORRIDOR_DEFAULT_MODEB_THRESHOLD_SUM = 1000;
static constexpr int CORRIDOR_DEFAULT_MODEB_THRESHOLD_MAX = 1000;

/* Step-scoped mode state. Reset to UNSET by gizmo_hydro_corridor_end().
 * File-static so consumers must go through the accessor; no extern reads.
 *
 * Lifetime invariant: UNSET -> {MODE_A, MODE_B} via decide_mode at the top
 * of compute_hydro_densities_and_forces, then back to UNSET via end at the
 * bottom of the same function. Reading UNSET from a consumer is a bug
 * (consumer placed outside the corridor span). */
static GizmoHydroCorridorMode g_corridor_mode = GizmoHydroCorridorMode::UNSET;

/* Global active-gas counts — collective Allreduce so every rank computes
 * the same mode (corridor coherence requires this). Fills the summed and the
 * per-rank-max active-gas counts (matching the runner's sum/max dispatch
 * inputs). Cheap (two int Allreduces). Caller must be inside MPI_COMM_WORLD
 * scope. */
static void
global_active_gas_counts(int *out_sum, int *out_max)
{
    int local_active_gas = 0;
    for(int ii : ActiveParticleList) {
        if(P[ii].Type == 0 && P[ii].Mass > 0) local_active_gas++;
    }
    if(NTask > 1) {
        MPI_Allreduce(&local_active_gas, out_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_active_gas, out_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    } else {
        *out_sum = local_active_gas;
        *out_max = local_active_gas;
    }
}

void gizmo_hydro_corridor_decide_mode(void)
{
#ifdef TRANSPORT_SUBCYCLE
    /* Transport subcycling reuses the shared symmetric gas CSR
     * (gizmo_sym_neighbor_list) across every subcycle iteration
     * (transport_subcycle_exchange_fluxes consumes it), so the corridor MUST
     * build that CSR — i.e. Mode A — on every hydro step. Mode B is
     * request-driven with no CSR and is therefore incompatible with
     * TRANSPORT_SUBCYCLE; it is not offered in this build, and the mode-force
     * env vars are intentionally not consulted. Forcing Mode A here (rather
     * than at each consumer) is the single policy point; it also lets prep/
     * refresh_symlist skip their Mode-B dead-work path without a subcycle
     * special-case. */
    g_corridor_mode = GizmoHydroCorridorMode::MODE_A;
    return;
#endif
    /* Adaptive dispatch, mirroring the runner's per-loop policy:
       Mode B iff sum > 0 && sum <= threshold_sum && max <= threshold_max.
       Thresholds from the optional NeighborLoopModeBThreshold{Sum,Max} params
       (-1 = unset -> corridor default; any set value overrides, <= 0 disables
       Mode B). TRANSPORT_SUBCYCLE has already forced Mode A above. */
    int sum_active_gas = 0, max_active_gas = 0;
    global_active_gas_counts(&sum_active_gas, &max_active_gas);
    const struct global_data_all_processes *host_all = nlr_host_all_ptr();
    int threshold_sum = (host_all->NeighborLoopModeBThresholdSum != -1)
                            ? host_all->NeighborLoopModeBThresholdSum
                            : CORRIDOR_DEFAULT_MODEB_THRESHOLD_SUM;
    int threshold_max = (host_all->NeighborLoopModeBThresholdMax != -1)
                            ? host_all->NeighborLoopModeBThresholdMax
                            : CORRIDOR_DEFAULT_MODEB_THRESHOLD_MAX;
    bool select_mode_b = (sum_active_gas > 0) && (sum_active_gas <= threshold_sum)
                                             && (max_active_gas <= threshold_max);
    GizmoHydroCorridorMode mode = select_mode_b ? GizmoHydroCorridorMode::MODE_B
                                                : GizmoHydroCorridorMode::MODE_A;
    g_corridor_mode = mode;

}

GizmoHydroCorridorMode gizmo_hydro_corridor_get_mode(void)
{
    return g_corridor_mode;
}

/* Corridor-owned view of the gizmo_sym_* shared CSR. Populated by
 * gizmo_hydro_corridor_begin() in Mode A at ANY rank count (at NTask>1
 * publication additionally requires a live ghost pool — the runner's
 * caller-owned-pool contract); republished by
 * gizmo_hydro_corridor_refresh_ghost_values() after each rebuild; cleared
 * by gizmo_hydro_corridor_end() at corridor exit. Does NOT own the
 * underlying buffer — gizmo_sym_neighbor_list is freed by
 * gizmo_hydro_cleanup_symlist_and_ghosts() after hydro_force(). */
static nlr_external_csr g_corridor_csr;
static bool             g_corridor_csr_valid = false;
/* Import epoch of the ghost pool the published CSR was built from (stamped at
 * begin and at every fallback republish). The value-refresh fast path in
 * gizmo_hydro_corridor_refresh_ghost_values runs ONLY when the CURRENT live
 * pool carries this epoch: a live pool from some other import could pass the
 * refresh count checks by coincidence while the CSR still indexes the old
 * slot layout — that case must take the full rebuild path. */
static unsigned long long g_corridor_pool_epoch = 0;

void gizmo_hydro_corridor_begin(void)
{
    if(g_corridor_mode == GizmoHydroCorridorMode::UNSET) return; /* defensive: decide_mode not run */

    /* Mode B: build ONLY the shared active-index list (the request-driven
     * walkers' row source for gradients/hydro_force) and dematerialize any
     * leftover ghosts — no import, no CSR, no arena work (the runner's Mode-B
     * hard-corridor counters enforce this downstream). Consumers no longer
     * build this list themselves. */
    const double gsl_safety = gizmo_ghost_safety_factor();
    if(g_corridor_mode == GizmoHydroCorridorMode::MODE_B) {
        gizmo_gradients_prep_symlist(gsl_safety, gsl_safety, true);
        return;
    }

    /* Mode A: import the shared gas ghost pool (NTask>1) and build the
     * active list + gizmo_sym_* CSR once for the whole corridor span. */
    gizmo_gradients_prep_symlist(gsl_safety, gsl_safety, false);

    /* Publish the corridor's view ONLY when the runner's caller-owned-pool
     * contract can be satisfied: at NTask>1 that requires a LIVE ghost pool
     * (prep skips the import when no rank has active gas — then nothing is
     * published and consumers take their fallback, which handles empty).
     * The gizmo_sym_* globals live until gizmo_hydro_cleanup_symlist_and_
     * ghosts() after hydro_force(), spanning all three consumers. */
    if(NTask > 1 && !ghost_pool_is_live()) return;

    g_corridor_csr.active_indices = gizmo_sym_active_indices;
    g_corridor_csr.num_active     = gizmo_sym_num_active;
    g_corridor_csr.offsets        = gizmo_sym_neighbor_list.offsets;
    g_corridor_csr.neighbors      = gizmo_sym_neighbor_list.neighbors;
    g_corridor_csr.total_pairs    = gizmo_sym_neighbor_list.total_pairs;
    g_corridor_csr.owner_name     = "corridor:hydro";
    g_corridor_csr_valid          = true;
    g_corridor_pool_epoch         = ghost_provenance_epoch();

}

const nlr_external_csr * gizmo_hydro_corridor_external_csr(void)
{
    return g_corridor_csr_valid ? &g_corridor_csr : nullptr;
}

void gizmo_hydro_corridor_refresh_ghost_values(const char *stage)
{
    /* No-op unless the corridor owns a published Mode-A view this step —
     * Mode B / fallback paths keep their own refresh semantics. At NTask==1
     * there are no ghosts and owner values are read directly: nothing to
     * refresh. */
    if(!g_corridor_csr_valid) return;
    if(NTask <= 1) return;

    /* FAST PATH: value-only in-place refresh of the existing ghost pool along
     * the provenance recorded at import (owners re-pack current P/CellP;
     * receivers overwrite the existing slots). Slot identity is preserved by
     * construction, so the pool, the provenance maps, and the published CSR
     * all stay valid — no discovery, no rebuild, no republish. Permitted ONLY
     * when the current live pool is the very import the CSR was built from
     * (epoch match; see g_corridor_pool_epoch). Any other state — pool torn
     * down by an intervening loop's own ghost lifecycle, epoch mismatch, or
     * a refresh-guard failure — falls CLOSED to the full path below. */
    if(ghost_pool_is_live() && ghost_provenance_epoch() == g_corridor_pool_epoch) {
        int refresh_rc = ghost_refresh_values();
        if(refresh_rc == GHOST_REFRESH_OK) {
            return;
        }
        /* fall through: fail-closed */
    }

    /* FULL PATH (also the fail-closed fallback): ghost teardown + re-import +
     * CSR rebuild, then REPUBLISH (the rebuild reallocates offsets/neighbors —
     * consumers re-fetch the view after every refresh; see hydro_corridor.h).
     * Import side-effects (the new ghost provenance the neighbour index keys
     * on) happen inside the re-import itself, so no side-effect list is
     * maintained here. */
    const double gsl_safety = gizmo_ghost_safety_factor();
    gizmo_gradients_refresh_symlist(gsl_safety, gsl_safety, false);

    /* Republish. The active set is frozen across the corridor span, so
     * active_indices/num_active are pointer-stable; only the CSR arrays
     * moved. Re-stamp the pool epoch: subsequent refreshes may fast-path
     * off this fresh import. */
    g_corridor_csr.offsets     = gizmo_sym_neighbor_list.offsets;
    g_corridor_csr.neighbors   = gizmo_sym_neighbor_list.neighbors;
    g_corridor_csr.total_pairs = gizmo_sym_neighbor_list.total_pairs;
    g_corridor_pool_epoch      = ghost_provenance_epoch();

}

void gizmo_hydro_corridor_end(void)
{
    g_corridor_mode      = GizmoHydroCorridorMode::UNSET;
    g_corridor_csr_valid = false;
    /* Note: do NOT free gizmo_sym_* here. The legacy
     * gizmo_hydro_cleanup_symlist_and_ghosts() (called after
     * hydro_force) still owns the free. Double-free would crash. */
}

