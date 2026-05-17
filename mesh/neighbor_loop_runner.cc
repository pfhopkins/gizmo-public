/* mesh/neighbor_loop_runner.cc — generic NeighborLoopSpec runner.
 *
 * Current scope (3c.2): WritePattern::ActiveReduceOnly,
 * RemoteEvalMode::RemoteComputesAccum, non-iterative, NoIdentity, NoScatter
 * — sufficient for SinkEnv1Spec migration.
 *
 * Mode A (GPU NGL pipeline): IMPLEMENTED in 3c.1. Stages caller-supplied
 * radii and per-call CallScalars host-side pre-arena, runs
 * gpu_particles_arena_acquire + gpu_ngb_list_build, then a tiny Kokkos
 * parallel_for that calls Spec::load_active to fill ActiveData[] in UVM
 * (same device epoch as the legacy lambda's q-packing), then the parametric
 * pair-kernel parallel_for calling Spec::load_neighbor + Spec::pair_kernel.
 * Both launches go through gizmo_gpu_kernel_launch (project's parallel_for
 * + fence + check_last_error).
 *
 * Mode B local + Brute oracle: IMPLEMENTED in 3c.2 (single-rank only). Same
 * Spec::load_active / Spec::load_neighbor / Spec::pair_kernel — host-side
 * KOKKOS_INLINE_FUNCTION invocation. Lazy-drift function-boundary invariant
 * structurally encoded as collect_candidates_pre_drift →
 * lazy_drift_candidates → evaluate_pairs_post_drift; Brute oracle uses the
 * SAME drift epoch as Mode B. Multi-rank peer-to-peer is deferred to 3c.3 — runtime
 * abort guards GIZMO_NLR_FORCE_MODE=B on NTask>1.
 *
 * Architecture binding contract:
 *   ~/.claude/memory/reference_neighbor_loop_contract.md
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) and Claude for GIZMO.
 */

#include <mpi.h>
#include <cstdarg>                                   /* va_list, vfprintf for nlr_warn_once_rank0 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <type_traits>
#include <utility>                                   /* std::declval for SFINAE hook detection */
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"          /* per-TU AllDeviceMirror + device-pass `All` redirect; include before allvars.h */
#include "../declarations/allvars.h"
#include "../declarations/lifecycle_counters.h"      /* g_global_drift_counter etc. for Mode B corridor */
#include "../core/proto.h"
#include "../system/gpu_particles_arena.h"
#include "../declarations/gpu_dispatch_templates.h"  /* gizmo_gpu_kernel_launch */

#include "neighbor_loop_runner.h"
#include "gpu_neighbor_list.h"
#include "kernel.h"  /* MUST precede sink_env1_loop.h (kernel_main, NEAREST_XYZ) */
#include "ghost_writeback.h"             /* ghost_get_num_local */
#include "ghost_symlist_lifecycle.h"     /* gizmo_request_filtered_ghost_import_fresh, ghost_exchange_cleanup */
#include "mode_b_local_walker.h"         /* mode_b_local_neighbor_walk, brute, lazy_drift */

#include <vector>
#include <unordered_map>
#include <cmath>
#include <cctype>

#include "mode_b_p2p_transport.h"  /* mode_b_exchange_queries / _replies */

/* Spec instantiations. Each #include declares one Spec type whose explicit
 * template instantiation appears at the bottom of this file. */
#include "../sinks/sink_env1_loop.h"
#include "../sinks/sink_feed_loop.h"
#include "../sinks/sink_swk_loop.h"
#if defined(SINK_PARTICLES) && defined(SINK_GRAVACCRETION) && (SINK_GRAVACCRETION == 0)
#include "../sinks/sink_env2_loop.h"
#endif

#ifdef GIZMO_NLR_ITER_HARNESS_TEST
#include "test_iter_harness_loop.h"
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
#include "../gravity/ags_density_loop.h"
#endif

#include "../hydro/density_loop.h"

#ifdef GALSF_FB_MECHANICAL
#include "../galaxy_sf/mechfb_loop.h"
#endif

#ifdef GALSF_FB_THERMAL
#include "../galaxy_sf/thermal_fb_loop.h"
#endif

#ifdef GALSF_FB_FIRE_RT_LOCALRP
#include "../galaxy_sf/radfb_rp_loop.h"
#endif
#if defined(GALSF_SUBGRID_WINDS) && (GALSF_SUBGRID_WIND_SCALING==2)
#include "../galaxy_sf/dm_dispersion_loop.h"
#endif
#ifdef TURB_DIFF_DYNAMIC
#include "../turb/difffilter_loop.h"
#endif

/* ============================================================================
 * Shared NLR utility helpers (used by env-config and threshold blocks below).
 * File-scope static; TU-local linkage. Defined here so the threshold helpers
 * (which are file-scope `extern` for runner.h API) can call them.
 * ========================================================================== */

/* One-shot rank-0 warning helper.
 *
 * Cached set keyed by string CONTENT (strcmp), not pointer identity, so
 * dynamically-constructed keys (e.g. per-loop env var names like
 * "GIZMO_SINK_ENV1_MODEB_THRESHOLD_SUM" assembled at runtime) dedupe
 * correctly per distinct env var rather than per category. Costs a 64x192
 * BSS array (~12 KB) and an O(N) lookup per call; both fine for a warning
 * surface.
 *
 * Cap is generous (64 distinct keys); if hit, subsequent warnings are
 * silently dropped — they are diagnostics, not correctness gates. */
static void nlr_warn_once_rank0(const char *key, const char *fmt, ...)
{
    if(ThisTask != 0) return;
    static char seen[64][192];
    static int  seen_n = 0;
    for(int i = 0; i < seen_n; i++) {
        if(strcmp(seen[i], key) == 0) return;
    }
    if(seen_n < 64) {
        size_t n = strlen(key);
        if(n >= sizeof(seen[0])) n = sizeof(seen[0]) - 1;
        memcpy(seen[seen_n], key, n);
        seen[seen_n][n] = '\0';
        seen_n++;
    }
    fprintf(stderr, "[NLR env] ");
    va_list ap; va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    fflush(stderr);
}

static bool nlr_env_is_one(const char *name) {
    const char *e = getenv(name);
    return (e && e[0] == '1' && e[1] == '\0');
}

/* ============================================================================
 * Caller-side helpers (declared in runner.h).
 *
 * Out-of-line so changes to global plumbing (NumPart, P, CellP fetch site,
 * ghost-safety-factor source, etc.) are caller-invisible.
 * ========================================================================== */

/* Thin wrapper over the canonical host accessor; retained for callers
 * that already route through this name. With the per-TU AllDeviceMirror
 * scheme's device-pass-gated `#define All AllDeviceMirror`, bare `All.*`
 * in host code reads the host extern unconditionally, so this
 * indirection is no longer load-bearing — kept for stability of
 * existing call sites. */
const struct global_data_all_processes * nlr_host_all_ptr(void)
{
    return gizmo_host_all_ptr();
}

NlrCommonScalars nlr_common_scalars_from_all(void)
{
    const struct global_data_all_processes *h = nlr_host_all_ptr();
    NlrCommonScalars s;
    s.cf_atime                = h->cf_atime;
    s.cf_a2inv                = h->cf_a2inv;
    s.cf_a3inv                = h->cf_a3inv;
    s.cf_hubble_a             = h->cf_hubble_a;
    s.newton_G                = h->G;
    s.hubble                  = h->HubbleParam;
    s.comoving_integration_on = h->ComovingIntegrationOn;
    return s;
}

neighbor_loop_args nlr_default_args(void)
{
    neighbor_loop_args args;
    args.P                   = P;
    args.CellP               = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
    args.num_total           = NumPart;
    args.active_list         = nullptr;       /* caller fills */
    args.num_active          = 0;             /* caller fills */
    args.aux                 = nullptr;       /* caller fills */
    args.ghost_safety_factor = gizmo_ghost_safety_factor();
    return args;
}

void nlr_free_active_list(int *active_list)
{
    if(active_list) myfree(active_list);
}

/* ============================================================================
 * TESTERS' KNOBS — not production policy
 *
 * The four env vars below are dispatch-threshold overrides that exist purely
 * for testing the runner's Mode A / Mode B selection. They are NOT part of
 * GIZMO's production interface, NOT promoted to params.txt or Config.sh,
 * and not intended for end-user tuning.
 *
 * Production dispatch policy is the constexpr Spec::modeb_threshold_sum and
 * Spec::modeb_threshold_max in each NeighborLoopSpec — those are code-level
 * dispatch policy constants for the loop, decided alongside the physics.
 *
 * Recognized env vars (all integer-valued):
 *   GIZMO_NLR_MODEB_THRESHOLD_SUM            global override; sum-of-active
 *   GIZMO_NLR_MODEB_THRESHOLD_MAX            global override; max-rank-active
 *   GIZMO_<UPPER_LOOP_NAME>_MODEB_THRESHOLD_SUM   per-loop override (e.g.
 *   GIZMO_<UPPER_LOOP_NAME>_MODEB_THRESHOLD_MAX   GIZMO_SINK_ENV1_MODEB_THRESHOLD_SUM)
 *
 * Resolution precedence (first wins; per-loop override beats global):
 *   1. per-loop env GIZMO_<LOOP>_MODEB_THRESHOLD_<SUM|MAX>
 *   2. global env  GIZMO_NLR_MODEB_THRESHOLD_<SUM|MAX>
 *   3. Spec::modeb_threshold_<sum|max> constexpr (code default)
 *
 * Invalid env values (non-integer / negative / > 1e9) are silently ignored
 * and the next precedence level is consulted. This preserves the existing
 * behavior; promoting invalid values to a hard endrun is an explicit policy
 * change deferred until threshold semantics stabilize.
 *
 * When an env override is actually consumed (level 1 or 2), a single rank-0
 * one-shot warning is emitted naming the env var, the resolved value, and
 * the spec_default it overrode, with a "tester only" tag.
 *
 * Force-mode env vars (Pass B.i) take precedence over threshold dispatch:
 * GIZMO_NLR_FORCE_MODE=A|B selects unconditionally; thresholds are not
 * consulted on force paths (except when GIZMO_NLR_DIAG>=1, where the
 * runner does an extra Allreduce to populate PHASE0_NLR num_active_global).
 * ========================================================================== */

static int s_global_threshold_sum_cached = -2;
static int s_global_threshold_max_cached = -2;
static int parse_int_env(const char *name)
{
    const char *e = getenv(name);
    if(!e || !e[0]) return -1;
    char *endp = nullptr;
    long v = strtol(e, &endp, 10);
    if(!endp || *endp != '\0' || v < 0 || v > 1000000000) return -1;
    return (int)v;
}
int gizmo_nlr_default_modeb_threshold_sum(void) {
    if(s_global_threshold_sum_cached == -2) {
        s_global_threshold_sum_cached = parse_int_env("GIZMO_NLR_MODEB_THRESHOLD_SUM");
    }
    return s_global_threshold_sum_cached;
}
int gizmo_nlr_default_modeb_threshold_max(void) {
    if(s_global_threshold_max_cached == -2) {
        s_global_threshold_max_cached = parse_int_env("GIZMO_NLR_MODEB_THRESHOLD_MAX");
    }
    return s_global_threshold_max_cached;
}

/* Per-loop env lookup: GIZMO_<UPPER_LOOP_NAME>_MODEB_THRESHOLD_<SUM|MAX>.
 * Returns the parsed int (or -1 on unset / invalid) AND, on success, writes
 * the constructed env-var name into out_name (size out_cap) so the caller
 * can include it in the tester-knob warning text. */
static int per_loop_threshold_lookup(const char *loop_name, const char *suffix,
                                      char *out_name, size_t out_cap)
{
    if(out_cap > 0) out_name[0] = '\0';
    if(!loop_name || !loop_name[0]) return -1;
    char env_name[160];
    int j = 0;
    const char *prefix = "GIZMO_";
    for(int p = 0; prefix[p] && j < 150; p++) env_name[j++] = prefix[p];
    for(int p = 0; loop_name[p] && j < 150; p++) {
        env_name[j++] = (char)toupper((unsigned char)loop_name[p]);
    }
    const char *mid = "_MODEB_THRESHOLD_";
    for(int p = 0; mid[p] && j < 158; p++) env_name[j++] = mid[p];
    for(int p = 0; suffix[p] && j < 159; p++) env_name[j++] = suffix[p];
    env_name[j] = '\0';
    int v = parse_int_env(env_name);
    if(v >= 0 && out_cap > 0) {
        size_t copy_n = (size_t)j;
        if(copy_n >= out_cap) copy_n = out_cap - 1;
        memcpy(out_name, env_name, copy_n);
        out_name[copy_n] = '\0';
    }
    return v;
}

int gizmo_nlr_modeb_threshold_sum_for(const char *loop_name, int spec_default)
{
    char per_loop_env[160];
    int v = per_loop_threshold_lookup(loop_name, "SUM", per_loop_env, sizeof(per_loop_env));
    if(v >= 0) {
        /* Dedup key = the per-loop env-var name itself, so distinct loops
         * each fire their own warning when their per-loop overrides are used. */
        nlr_warn_once_rank0(per_loop_env,
            "Tester knob in use: %s=%d overrides Spec::modeb_threshold_sum=%d for loop '%s' "
            "(not a production interface).",
            per_loop_env, v, spec_default, loop_name ? loop_name : "?");
        return v;
    }
    v = gizmo_nlr_default_modeb_threshold_sum();
    if(v >= 0) {
        nlr_warn_once_rank0("GIZMO_NLR_MODEB_THRESHOLD_SUM",
            "Tester knob in use: GIZMO_NLR_MODEB_THRESHOLD_SUM=%d overrides "
            "Spec::modeb_threshold_sum=%d (loop '%s'; not a production interface).",
            v, spec_default, loop_name ? loop_name : "?");
        return v;
    }
    return spec_default;
}
int gizmo_nlr_modeb_threshold_max_for(const char *loop_name, int spec_default)
{
    char per_loop_env[160];
    int v = per_loop_threshold_lookup(loop_name, "MAX", per_loop_env, sizeof(per_loop_env));
    if(v >= 0) {
        /* Dedup key = the per-loop env-var name itself (see _sum_for). */
        nlr_warn_once_rank0(per_loop_env,
            "Tester knob in use: %s=%d overrides Spec::modeb_threshold_max=%d for loop '%s' "
            "(not a production interface).",
            per_loop_env, v, spec_default, loop_name ? loop_name : "?");
        return v;
    }
    v = gizmo_nlr_default_modeb_threshold_max();
    if(v >= 0) {
        nlr_warn_once_rank0("GIZMO_NLR_MODEB_THRESHOLD_MAX",
            "Tester knob in use: GIZMO_NLR_MODEB_THRESHOLD_MAX=%d overrides "
            "Spec::modeb_threshold_max=%d (loop '%s'; not a production interface).",
            v, spec_default, loop_name ? loop_name : "?");
        return v;
    }
    return spec_default;
}

/* ============================================================================
 * NLR env config (Pass B.i unified API).
 *
 * Single canonical surface for diagnostic, control, and spike (cross-
 * validation) env vars. Old names are accepted as aliases for one cycle
 * with rank-0 deprecation warnings; explicit retire queued for the next
 * cleanup pass after Pass B.
 *
 * Conflict policy and alias precedence: see the comment block on the
 * declarations at the top of mesh/neighbor_loop_runner.h.
 *
 * All accessors are first-use cached and lock-free (single load per call).
 * Initialization may call endrun on a hard conflict; endrun is collective,
 * and TACC env vars are uniform across ranks, so all ranks reach the same
 * decision and abort together.
 * ========================================================================== */

namespace {

/* nlr_warn_once_rank0 and nlr_env_is_one are defined at file scope above
 * (near the includes) so the threshold block can use them too. */

/* Initialize diag level. Reads new var first, then aliases. */
static int nlr_init_diag_level(void)
{
    const char *raw = getenv("GIZMO_NLR_DIAG");
    int new_set_level = -1;
    if(raw && raw[0]) {
        char *endp = nullptr;
        long v = strtol(raw, &endp, 10);
        if(!endp || *endp != '\0' || v < 0 || v > 3) {
            if(ThisTask == 0) {
                fprintf(stderr, "[NLR env] FATAL: GIZMO_NLR_DIAG=\"%s\" must be in {0,1,2,3}.\n",
                        raw);
                fflush(stderr);
            }
            endrun(81100);
        }
        new_set_level = (int)v;
    }

    bool old_phase0   = nlr_env_is_one("GIZMO_PHASE0_DIAG");
    bool old_dispatch = nlr_env_is_one("GIZMO_NLR_DISPATCH_TRACE");

    int level;
    if(new_set_level >= 0) {
        level = new_set_level;
        if(old_phase0) {
            nlr_warn_once_rank0("alias_phase0_overridden",
                "GIZMO_PHASE0_DIAG ignored; GIZMO_NLR_DIAG=%d takes precedence.", level);
        }
        if(old_dispatch) {
            nlr_warn_once_rank0("alias_dispatch_overridden",
                "GIZMO_NLR_DISPATCH_TRACE ignored; GIZMO_NLR_DIAG=%d takes precedence.", level);
        }
    } else {
        level = 0;
        if(old_phase0) {
            nlr_warn_once_rank0("alias_phase0_deprecated",
                "GIZMO_PHASE0_DIAG=1 is deprecated; use GIZMO_NLR_DIAG=1 instead.");
            if(level < 1) level = 1;
        }
        if(old_dispatch) {
            nlr_warn_once_rank0("alias_dispatch_deprecated",
                "GIZMO_NLR_DISPATCH_TRACE=1 is deprecated; use GIZMO_NLR_DIAG=2 instead.");
            if(level < 2) level = 2;
        }
    }

    if(level == 3) {
        nlr_warn_once_rank0("level_3_reserved",
            "GIZMO_NLR_DIAG=3: level 3 currently has no extra diagnostics; "
            "reserved for future scalar-only extensions. Behaving as level 2.");
    }

    return level;
}

/* Initialize force mode. Conflict cases endrun (collective). */
static NlrForceMode nlr_init_force_mode(void)
{
    const char *raw_new = getenv("GIZMO_NLR_FORCE_MODE");
    bool new_set = (raw_new && raw_new[0]);
    bool old_a   = nlr_env_is_one("GIZMO_NLR_FORCE_MODEA");
    bool old_b   = nlr_env_is_one("GIZMO_NLR_FORCE_MODEB");

    if(new_set && (old_a || old_b)) {
        if(ThisTask == 0) {
            fprintf(stderr, "[NLR env] FATAL: GIZMO_NLR_FORCE_MODE='%s' is set AND "
                    "old GIZMO_NLR_FORCE_MODEA=%d / GIZMO_NLR_FORCE_MODEB=%d are set. "
                    "Use only one. Old names are deprecated.\n",
                    raw_new, (int)old_a, (int)old_b);
            fflush(stderr);
        }
        endrun(81101);
    }
    if(old_a && old_b) {
        if(ThisTask == 0) {
            fprintf(stderr, "[NLR env] FATAL: GIZMO_NLR_FORCE_MODEA and GIZMO_NLR_FORCE_MODEB "
                    "are both set; mutually exclusive.\n");
            fflush(stderr);
        }
        endrun(81102);
    }

    if(new_set) {
        if(raw_new[0] == 'A' && raw_new[1] == '\0') return NlrForceMode::A;
        if(raw_new[0] == 'B' && raw_new[1] == '\0') return NlrForceMode::B;
        if(ThisTask == 0) {
            fprintf(stderr, "[NLR env] FATAL: GIZMO_NLR_FORCE_MODE=\"%s\" must be 'A' or 'B'.\n",
                    raw_new);
            fflush(stderr);
        }
        endrun(81103);
    }

    if(old_a) {
        nlr_warn_once_rank0("alias_force_modea_deprecated",
            "GIZMO_NLR_FORCE_MODEA=1 is deprecated; use GIZMO_NLR_FORCE_MODE=A instead.");
        return NlrForceMode::A;
    }
    if(old_b) {
        nlr_warn_once_rank0("alias_force_modeb_deprecated",
            "GIZMO_NLR_FORCE_MODEB=1 is deprecated; use GIZMO_NLR_FORCE_MODE=B instead.");
        return NlrForceMode::B;
    }
    return NlrForceMode::None;
}

static NlrForceMode nlr_parse_force_mode_value(const char *env_name,
                                               const char *raw_value)
{
    if(!raw_value || !raw_value[0]) return NlrForceMode::None;
    if(raw_value[0] == 'A' && raw_value[1] == '\0') return NlrForceMode::A;
    if(raw_value[0] == 'B' && raw_value[1] == '\0') return NlrForceMode::B;
    if(ThisTask == 0) {
        fprintf(stderr, "[NLR env] FATAL: %s=\"%s\" must be 'A' or 'B'.\n",
                env_name ? env_name : "GIZMO_<LOOP>_FORCE_MODE", raw_value);
        fflush(stderr);
    }
    endrun(81104);
    return NlrForceMode::None;
}

static void nlr_loop_env_name(const char *loop_name, const char *suffix,
                              char *out, int out_size)
{
    if(out_size <= 0) return;
    int j = 0;
    const char *prefix = "GIZMO_";
    for(int p = 0; prefix[p] && j < out_size - 1; p++) out[j++] = prefix[p];
    if(loop_name) {
        for(int p = 0; loop_name[p] && j < out_size - 1; p++) {
            out[j++] = (char)toupper((unsigned char)loop_name[p]);
        }
    }
    if(suffix) {
        for(int p = 0; suffix[p] && j < out_size - 1; p++) out[j++] = suffix[p];
    }
    out[j] = '\0';
}

static int nlr_force_modeb_active_cap_for(const char *loop_name)
{
    char env_name[160];
    nlr_loop_env_name(loop_name, "_FORCE_MODEB_MAX_ACTIVE",
                      env_name, (int)sizeof(env_name));
    int cap = parse_int_env(env_name);
    if(cap >= 0) return cap;
    cap = parse_int_env("GIZMO_NLR_FORCE_MODEB_MAX_ACTIVE");
    return (cap >= 0) ? cap : 100000;
}

static void nlr_abort_if_forced_modeb_too_large(const char *loop_name,
                                                int local_active,
                                                int global_active)
{
    const int cap = nlr_force_modeb_active_cap_for(loop_name);
    if(cap >= 0 && global_active > cap) {
        if(ThisTask == 0) {
            fprintf(stderr,
                    "[NLR FORCE_MODE=B] FATAL: caller=%s requested forced "
                    "Mode B with global_active=%d local_active(rank0)=%d, "
                    "exceeding cap=%d. This prevents accidental full-N "
                    "host-walker/oracle runs during startup density setup. "
                    "Use GIZMO_<LOOP>_FORCE_MODE=A for dense loops, a "
                    "per-loop Mode-B override for tiny loops, or raise "
                    "GIZMO_NLR_FORCE_MODEB_MAX_ACTIVE intentionally.\n",
                    loop_name ? loop_name : "?", global_active,
                    local_active, cap);
            fflush(stderr);
        }
        endrun(81105);
    }
}

static bool nlr_init_spike_accum_dump(void)
{
    bool new_set = nlr_env_is_one("GIZMO_NLR_SPIKE_ACCUM_DUMP");
    bool old_set = nlr_env_is_one("GIZMO_MODE_B_XVAL_DUMP");
    if(new_set) {
        if(old_set) {
            nlr_warn_once_rank0("alias_xval_dump_overridden",
                "GIZMO_MODE_B_XVAL_DUMP ignored; GIZMO_NLR_SPIKE_ACCUM_DUMP takes precedence.");
        }
        return true;
    }
    if(old_set) {
        nlr_warn_once_rank0("alias_xval_dump_deprecated",
            "GIZMO_MODE_B_XVAL_DUMP=1 is deprecated; use GIZMO_NLR_SPIKE_ACCUM_DUMP=1 instead.");
        return true;
    }
    return false;
}

static bool nlr_init_spike_nb_dump(void)
{
    bool new_set = nlr_env_is_one("GIZMO_NLR_SPIKE_NB_DUMP");
    bool old_set = nlr_env_is_one("GIZMO_MODE_B_XVAL_NB_DUMP");
    if(new_set) {
        if(old_set) {
            nlr_warn_once_rank0("alias_xval_nb_dump_overridden",
                "GIZMO_MODE_B_XVAL_NB_DUMP ignored; GIZMO_NLR_SPIKE_NB_DUMP takes precedence.");
        }
        return true;
    }
    if(old_set) {
        nlr_warn_once_rank0("alias_xval_nb_dump_deprecated",
            "GIZMO_MODE_B_XVAL_NB_DUMP=1 is deprecated; use GIZMO_NLR_SPIKE_NB_DUMP=1 instead.");
        return true;
    }
    return false;
}

} /* anonymous namespace */

int gizmo_nlr_diag_level(void)
{
    static int cached = -1;
    if(cached < 0) cached = nlr_init_diag_level();
    return cached;
}

NlrForceMode gizmo_nlr_force_mode(void)
{
    static int  cached_int   = -1;          /* -1 = uninit; 0/1/2 = enum */
    if(cached_int < 0) cached_int = (int)nlr_init_force_mode();
    return (NlrForceMode)cached_int;
}

NlrForceMode gizmo_nlr_force_mode_for(const char *loop_name)
{
    if(!loop_name || !loop_name[0]) return gizmo_nlr_force_mode();
    struct entry_t { const char *name; int cached; };
    static entry_t cache[16] = {
        {nullptr,-1},{nullptr,-1},{nullptr,-1},{nullptr,-1},
        {nullptr,-1},{nullptr,-1},{nullptr,-1},{nullptr,-1},
        {nullptr,-1},{nullptr,-1},{nullptr,-1},{nullptr,-1},
        {nullptr,-1},{nullptr,-1},{nullptr,-1},{nullptr,-1}};
    for(int k = 0; k < 16; k++) {
        if(cache[k].name == loop_name) return (NlrForceMode)cache[k].cached;
        if(cache[k].name == nullptr) {
            char env_name[160];
            nlr_loop_env_name(loop_name, "_FORCE_MODE", env_name, (int)sizeof(env_name));
            const char *raw = getenv(env_name);
            NlrForceMode mode = raw && raw[0]
                ? nlr_parse_force_mode_value(env_name, raw)
                : gizmo_nlr_force_mode();
            cache[k].name = loop_name;
            cache[k].cached = (int)mode;
            return mode;
        }
    }
    char env_name[160];
    nlr_loop_env_name(loop_name, "_FORCE_MODE", env_name, (int)sizeof(env_name));
    const char *raw = getenv(env_name);
    return raw && raw[0] ? nlr_parse_force_mode_value(env_name, raw)
                         : gizmo_nlr_force_mode();
}

bool gizmo_nlr_spike_accum_dump_enabled(void)
{
    static int cached = -1;
    if(cached < 0) cached = nlr_init_spike_accum_dump() ? 1 : 0;
    return cached != 0;
}

bool gizmo_nlr_spike_nb_dump_enabled(void)
{
    static int cached = -1;
    if(cached < 0) cached = nlr_init_spike_nb_dump() ? 1 : 0;
    return cached != 0;
}

/* Adapters — preserve existing call-site names. */
bool gizmo_nlr_phase0_diag_enabled(void)    { return gizmo_nlr_diag_level() >= 1; }
bool gizmo_nlr_dispatch_trace_enabled(void) { return gizmo_nlr_diag_level() >= 2; }
bool gizmo_nlr_xval_dump_enabled(void)      { return gizmo_nlr_spike_accum_dump_enabled(); }
bool gizmo_nlr_xval_nb_dump_enabled(void)   { return gizmo_nlr_spike_nb_dump_enabled(); }

/* Subgroup-audit env gate (step 2c.3 step 10). Off by default; harness +
 * test-mode runs set GIZMO_NLR_SUBGROUP_AUDIT=1 to enable hard-aborting
 * runtime checks for subgroups[] collective-symmetry + no-duplicate-actives. */
static bool gizmo_nlr_subgroup_audit_enabled(void) {
    return nlr_env_is_one("GIZMO_NLR_SUBGROUP_AUDIT");
}

/* ============================================================================
 * NeighborLoopPlan path predicates — single source of truth keyed on path.
 *
 * New paths (future Mode C, dual-tree large-N, etc.) extend the switch
 * statements below; never add fields to NeighborLoopPlan.
 * ========================================================================== */
bool nlr_path_uses_imported_ghosts(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
    }
    return false;
}

bool nlr_path_uses_gpu_arena(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
    }
    return false;
}

bool nlr_path_permits_global_numpart_mutation(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
    }
    return false;
}

bool nlr_path_uses_lazy_drift(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return false;
        case NeighborLoopPlan::Path::ModeB_Local:  return true;
        case NeighborLoopPlan::Path::ModeB_Remote: return true;
    }
    return false;
}

const char *nlr_path_label(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return "gpu_ngl";
        case NeighborLoopPlan::Path::ModeB_Local:  return "mode_b_local";
        case NeighborLoopPlan::Path::ModeB_Remote: return "mode_b_remote";
    }
    return "unknown";
}

/* ============================================================================
 * StageTimer — internal RAII helper for PHASE0 timing (3c.4b)
 *
 * `target == nullptr` means phase0 is off: ctor and dtor do nothing except
 * a single predictable nullptr branch. NO MPI_Wtime call when off — codex
 * constraint that "MPI_Wtime is not zero overhead" addressed at the call
 * site, not just at the gating env var.
 *
 * Targets are accumulators: multiple StageTimer scopes can target the same
 * field (e.g. Mode B remote's dt_collect spans BOTH self and peer pre-drift
 * collection — two scopes accumulate). dt_total spans the whole runner call
 * and is set explicitly at top-level, not via this helper.
 * ========================================================================== */
namespace {
struct StageTimer {
    double *target;
    double t0;
    explicit StageTimer(double *tgt) : target(tgt), t0(0.0) {
        if(target) t0 = MPI_Wtime();
    }
    ~StageTimer() {
        if(target) *target += MPI_Wtime() - t0;
    }
    StageTimer(const StageTimer&) = delete;
    StageTimer& operator=(const StageTimer&) = delete;
};
} /* anonymous namespace */

bool gizmo_nlr_force_mode_b_global(void) { return gizmo_nlr_force_mode() == NlrForceMode::B; }
bool gizmo_nlr_force_mode_a_global(void) { return gizmo_nlr_force_mode() == NlrForceMode::A; }
bool gizmo_nlr_oracle_enabled_global(void) {
    static int cached = -1;
    if(cached < 0) {
        const char *e = getenv("GIZMO_NLR_ORACLE");
        cached = (e && e[0] == '1') ? 1 : 0;
    }
    return cached != 0;
}
bool gizmo_nlr_oracle_enabled_for(const char *loop_name) {
    /* Per-loop env: GIZMO_<UPPERCASE_LOOP_NAME>_ORACLE. Lookup is per-call
     * but loop_name is a Spec::loop_name constexpr literal, so the env-name
     * derivation is deterministic; we cache via a small intrusive map of up
     * to 8 entries (cheap linear scan). 8 covers every Spec we expect through
     * 3c-3g; bump if a future stage exceeds. */
    if(!loop_name || !loop_name[0]) return false;
    struct entry_t { const char *name; int cached; };
    static entry_t cache[8] = {
        {nullptr,0},{nullptr,0},{nullptr,0},{nullptr,0},
        {nullptr,0},{nullptr,0},{nullptr,0},{nullptr,0}};
    for(int k = 0; k < 8; k++) {
        if(cache[k].name == loop_name) return cache[k].cached != 0;
        if(cache[k].name == nullptr) {
            char env_name[128];
            int j = 0;
            const char *prefix = "GIZMO_";
            for(int p = 0; prefix[p] && j < 120; p++) env_name[j++] = prefix[p];
            for(int p = 0; loop_name[p] && j < 120; p++) {
                env_name[j++] = (char)toupper((unsigned char)loop_name[p]);
            }
            const char *suffix = "_ORACLE";
            for(int p = 0; suffix[p] && j < 127; p++) env_name[j++] = suffix[p];
            env_name[j] = '\0';
            const char *e = getenv(env_name);
            cache[k].name   = loop_name;
            cache[k].cached = (e && e[0] == '1') ? 1 : 0;
            return cache[k].cached != 0;
        }
    }
    /* Cache full — fall through with a one-shot lookup, no caching. */
    return false;
}

/* ============================================================================
 * SIDX cache kind resolver (private to runner; see neighbor_loop_runner.h
 * SidxCacheKind enum doc).
 * ========================================================================== */

static gpu_spatial_index_t* nlr_resolve_sidx_cache(SidxCacheKind k,
                                                    const char *loop_name)
{
    switch(k) {
        case SidxCacheKind::AllTypes:
            return gpu_step_sidx_alltypes_ptr();
        case SidxCacheKind::GasOnly:
            return gpu_step_sidx_ptr();
        case SidxCacheKind::None:
            return nullptr;
    }
    /* Unreachable today — kept exhaustive for compiler warnings on enum
     * additions. Other kinds land alongside their first caller. */
    fprintf(stderr, "neighbor_loop_runner: SidxCacheKind=%d not implemented "
            "for loop '%s'\n",
            (int)k, loop_name ? loop_name : "?");
    fflush(stderr);
    endrun(81030);
    return nullptr;
}

/* ============================================================================
 * Mode B local self-rank helpers (3c.2)
 *
 * These three helpers STRUCTURALLY ENCODE the lazy-drift invariant from the
 * neighbor-loop binding contract:
 *
 *     collect_candidates_pre_drift<Spec>   — search backend (tree or brute)
 *                                            runs against possibly-stale P[j]
 *     lazy_drift_candidates<Spec>          — drift_particle on every j to
 *                                            All.Ti_Current (idempotent;
 *                                            duplicate j's between Mode B and
 *                                            Brute oracle are dedupe-free
 *                                            via drift_particle's
 *                                            time1==time0 early-return)
 *     evaluate_pairs_post_drift<Spec>      — calls Spec::pair_kernel via the
 *                                            same KOKKOS_INLINE_FUNCTION
 *                                            Spec::load_active /
 *                                            Spec::load_neighbor used by
 *                                            run_mode_a (host invocation here)
 *
 * Brute oracle uses the SAME drift epoch as Mode B (we collect both candidate
 * sets before drifting either — see run_mode_b_local_with_oracle). Annotation
 * on each helper makes the ordering enforcement structural.
 *
 * Self-rank only in 3c.2: candidates are local real P[] indices in
 * [0, num_local). Cross-rank peer-to-peer comes in 3c.3.
 *
 * Walker buffer sized to num_local (worst-case SYMMETRIC, no h-bound
 * pre-pruning).
 * ========================================================================== */

/* WALK-ONLY. Does NOT mutate P[].Pos/Vel — drift_particle must not be
 * called from inside this helper. (Audited 2026-05-08: walker calls only
 * force_drift_node on tree-internal nodes, which is search-side state.)
 *
 * Per-active variant: walks args.P[i].Pos at radii[aa] for each
 * aa in [0,args.num_active). Self-rank queries.
 */
template <typename Spec>
static void collect_candidates_pre_drift(const neighbor_loop_args& args,
                                          const double *radii,
                                          unsigned int neighbor_type_mask,
                                          DispatchPath backend,
                                          std::vector<std::vector<int>>& per_active_cands)
{
    /* neighbor_type_mask is explicit caller parameter (step 2c.3
     * mask-threading refactor 2026-05-10). Non-iter callers pass
     * Spec::neighbor_type_mask (unchanged behavior); iter dispatch
     * passes sg.j_type_bitmask for per-subgroup walks. */
    const int N = args.num_active;
    const int num_local = ghost_get_num_local();
    /* Per-call ownership: each inner vector owns its capacity for the
     * duration of this call. .assign(N, {}) clobbers any stale residue
     * from a previous call. Stage 4 contract: walker appends via
     * push_back; geometric growth handles any-size match set without
     * imposing the previous full-pool .assign(num_local, 0) cost
     * (~24 MB × N_active on fire_m11i). */
    per_active_cands.assign(N, std::vector<int>{});
    if(num_local <= 0) return;
    for(int aa = 0; aa < N; aa++) {
        const int i = args.active_list[aa];
        const double h_q = radii[aa];
        if(h_q <= 0) continue;
        std::vector<int>& cands = per_active_cands[aa];
        cands.clear();
        if(cands.capacity() == 0) cands.reserve(64); /* small initial; grows geometrically */
        double pos_arr[3] = {(double)args.P[i].Pos[0],
                              (double)args.P[i].Pos[1],
                              (double)args.P[i].Pos[2]};
        if(backend == DispatchPath::ModeB_HostWalker) {
            mode_b_local_neighbor_walk(pos_arr, h_q,
                                        neighbor_type_mask,
                                        Spec::search_mode,
                                        Spec::radius_policy,
                                        cands,
                                        nlr_spec_symmetric_j_radius_scale<Spec>());
        } else if(backend == DispatchPath::Brute_Oracle) {
            mode_b_local_brute_walk(pos_arr, h_q,
                                     neighbor_type_mask,
                                     Spec::search_mode,
                                     Spec::radius_policy,
                                     cands,
                                     nlr_spec_symmetric_j_radius_scale<Spec>());
        } else {
            fprintf(stderr, "neighbor_loop_runner: collect_candidates_pre_drift "
                    "called with non-Mode-B/Brute backend (%d) for loop '%s'\n",
                    (int)backend, Spec::loop_name);
            fflush(stderr);
            endrun(81033);
            return;
        }
    }
}

/* WALK-ONLY. Peer-side variant: walks against the LOCAL pool using each
 * remote query's pos/h_search drawn from peer_actives[k].{pos,h_search}.
 * Used for queries received from other ranks via the peer-to-peer transport.
 *
 * Pulls pos/h_search directly from Spec::ActiveData fields. The current
 * convention is that every Spec exposes `pos` and `h_search` as flat
 * ActiveData fields (matches sink_env1_loop.h template). Generalizing
 * this access pattern (e.g. to a Spec::query_pos/query_h trait pair) is
 * tracked for the runner-template-hardening pass.
 */
template <typename Spec>
static void collect_candidates_for_remote_queries(
    const std::vector<typename Spec::ActiveData>& peer_actives,
    unsigned int neighbor_type_mask,
    DispatchPath backend,
    std::vector<std::vector<int>>& per_query_cands)
{
    /* neighbor_type_mask is explicit caller parameter (step 2c.3
     * mask-threading refactor 2026-05-10). Non-iter callers pass
     * Spec::neighbor_type_mask; iter dispatch passes sg.j_type_bitmask. */
    const int K = (int)peer_actives.size();
    const int num_local = ghost_get_num_local();
    per_query_cands.assign(K, std::vector<int>{});
    if(num_local <= 0) return;
    for(int k = 0; k < K; k++) {
        const auto& active = peer_actives[k];
        const double h_q = (double)active.h_search;
        if(h_q <= 0) continue;
        std::vector<int>& cands = per_query_cands[k];
        cands.clear();
        if(cands.capacity() == 0) cands.reserve(64);
        double pos_arr[3] = {(double)active.pos[0], (double)active.pos[1], (double)active.pos[2]};
        if(backend == DispatchPath::ModeB_HostWalker) {
            mode_b_local_neighbor_walk(pos_arr, h_q,
                                        neighbor_type_mask,
                                        Spec::search_mode,
                                        Spec::radius_policy,
                                        cands,
                                        nlr_spec_symmetric_j_radius_scale<Spec>());
        } else if(backend == DispatchPath::Brute_Oracle) {
            mode_b_local_brute_walk(pos_arr, h_q,
                                     neighbor_type_mask,
                                     Spec::search_mode,
                                     Spec::radius_policy,
                                     cands,
                                     nlr_spec_symmetric_j_radius_scale<Spec>());
        } else {
            fprintf(stderr, "neighbor_loop_runner: collect_candidates_for_remote_queries"
                    " bad backend %d for loop '%s'\n", (int)backend, Spec::loop_name);
            fflush(stderr);
            endrun(81033);
            return;
        }
    }
}

/* SAME drift epoch contract. Idempotent: drift_particle's time1==time0
 * early-return makes calling this multiple times (e.g. once each on
 * self_tree, self_brute, peer_tree, peer_brute candidate sets) safe. */
template <typename Spec>
static void lazy_drift_candidates(std::vector<std::vector<int>>& per_active_cands)
{
    for(auto& v : per_active_cands) {
        if(!v.empty()) {
            mode_b_lazy_drift_candidates(v.data(), (int)v.size());
        }
    }
}

/* Evaluates pair_kernel POST-DRIFT for a precomputed ActiveData[] paired with
 * candidate lists. Caller is responsible for having frozen actives[] before
 * any drift, and for having drifted the candidate union before calling.
 *
 * Decoupled from args.active_list: works for self-rank (actives built from
 * args.active_list) and peer-rank (actives built from received envelopes).
 * Same KOKKOS_INLINE_FUNCTION Spec::load_neighbor + Spec::pair_kernel either
 * way.
 */
template <typename Spec, typename DeviceCtx>
static void evaluate_pairs_post_drift(const DeviceCtx& ctx,
                                       const typename Spec::ActiveData *actives,
                                       int N,
                                       const std::vector<std::vector<int>>& per_active_cands,
                                       typename Spec::AccumData *accums)
{
    using NeighborData = typename Spec::NeighborData;
    using ScatterData  = typename Spec::ScatterData;
    for(int aa = 0; aa < N; aa++) {
        Spec::zero_accum(accums[aa]);
        const auto& a = actives[aa];
        ScatterData s{};                              /* NoScatter for ActiveReduceOnly */
        const auto& cands = per_active_cands[aa];
        for(size_t kk = 0; kk < cands.size(); kk++) {
            int j = cands[kk];
            IdentitySidecar id{};                     /* NoIdentity */
            NeighborData nb = Spec::load_neighbor(ctx, j, id, a);
            Spec::pair_kernel(a, nb, accums[aa], s);
        }
    }
}

/* ============================================================================
 * run_mode_b_local<Spec> — single-rank Mode B path.
 *
 * Three-epoch staging contract (matches Mode A timing exactly except step 3
 * runs host-side per query instead of inside a Kokkos kernel):
 *   (1) host pre-drift: Spec::search_radius  → radii[num_active]
 *   (2) host pre-drift: Spec::populate_call_scalars → CallScalars cs
 *   (3) host post-drift: Spec::load_active per query (KOKKOS_INLINE_FUNCTION
 *                         invoked host-side; same UVM/Pos epoch as Mode A's
 *                         device kernel)
 *
 * No Kokkos kernels in this path → no fence (project rule from contract).
 * No GPU NGL build, no SIDX, no arena_acquire — Mode B's whole point is
 * escaping that overhead for tiny-N. `compare_accum`-gated brute oracle is
 * orchestrated by run_mode_b_local_with_oracle below.
 * ========================================================================== */

/* Helper: build the host-frozen ActiveData[] for self-rank actives. Called
 * BEFORE any drift so the active snapshot matches the legacy pack_query
 * epoch. (Codex nuance: this is NOT the same epoch as Mode A's device-staged
 * load_active post-NGL-build; documented at run_neighbor_loop dispatch site.)
 */
template <typename Spec, typename DeviceCtx>
static void build_self_actives_host_pre_drift(
    const neighbor_loop_args& args,
    const DeviceCtx& ctx,
    const double *radii,
    const typename Spec::CallScalars& cs,
    typename Spec::ActiveData *actives_out)
{
    const int N = args.num_active;
    for(int aa = 0; aa < N; aa++) {
        actives_out[aa] = Spec::load_active(ctx, aa, args.active_list[aa],
                                             radii[aa], cs);
    }
}

/* ============================================================================
 * Diagnostic active-dump emit helper (3c.4a)
 *
 * Iterates per-active and calls Spec::diagnostic_dump_active when the env
 * gate is on AND the spec opts in via the SFINAE-detected hook.
 *
 * Timing invariant: post-Spec::apply_active_writeback, pre-caller-side
 * SinkTempInfo scatter. Fires inside the runner before returning to the
 * caller, so the caller's ghost_write_detector_end / scatter loop have
 * not yet run. (Legacy emit fired post-ghost-detector but pre-scatter; the
 * runner-driven emit fires pre-ghost-detector but still pre-scatter.
 * Accumulator values are byte-identical between the two timing points —
 * writeback is the last write into the per-active accumulator buffer
 * before scatter — so the cross-validation dump line content is unchanged.)
 *
 * Mode A path supplies actives_or_null = d_actives (UVM-resident; host-
 * coherent post-fence). Mode B local/remote supply host-frozen actives.
 * Spec hooks read view.active->{pos,h_search,...} fields when needed.
 * ========================================================================== */
template <typename Spec>
static void runner_emit_active_dumps(const neighbor_loop_args& args,
                                      const typename Spec::AccumData *accums,
                                      const typename Spec::ActiveData *actives_or_null,
                                      const char *path)
{
    if(!gizmo_nlr_xval_dump_enabled()) return;
    if constexpr (spec_has_dump_active<Spec>::value) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for(int aa = 0; aa < args.num_active; aa++) {
            ActiveDumpView<Spec> v;
            v.rank        = rank;
            v.origin_rank = rank;
            v.origin_slot = aa;
            v.active_slot = aa;
            v.path        = path;
            v.call_id     = 0;
            v.args        = &args;
            v.active      = actives_or_null ? &actives_or_null[aa] : nullptr;
            v.accum       = &accums[aa];
            Spec::diagnostic_dump_active(v);
        }
        std::fflush(stdout);
    }
}

template <typename Spec>
static void run_mode_b_local(const neighbor_loop_args& args, const double *radii,
                             RunnerStageTimer *tim = nullptr)
{
    using ActiveData = typename Spec::ActiveData;
    using AccumData  = typename Spec::AccumData;
    using DeviceCtx  = typename Spec::DeviceContext;

    const int N = args.num_active;
    if(N <= 0) {
        /* Local-zero is legitimate (caller may enter unconditionally if
         * global_num_active > 0 — multi-rank Mode B requires that). */
        return;
    }

    /* (1) Radii are runner-staged and passed in; pointer is call-lifetime
     * only (do not store). See run_neighbor_loop contract in the header. */

    /* (2) Host pre-drift: per-call scalar globals into POD. */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    DeviceCtx ctx;
    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; the runner captures it by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);

    /* Freeze active snapshots host-side BEFORE drift. */
    std::vector<ActiveData> actives(N);
    build_self_actives_host_pre_drift<Spec>(args, ctx, radii, cs, actives.data());

    /* Helper layout: collect → drift → evaluate. */
    std::vector<std::vector<int>> cand_modeB;
    {
        StageTimer t(tim ? &tim->dt_collect : nullptr);
        collect_candidates_pre_drift<Spec>(args, radii,
                                            (unsigned int)Spec::neighbor_type_mask,
                                            DispatchPath::ModeB_HostWalker, cand_modeB);
    }
    {
        StageTimer t(tim ? &tim->dt_drift : nullptr);
        lazy_drift_candidates<Spec>(cand_modeB);
    }

    std::vector<AccumData> accums(N);
    {
        StageTimer t(tim ? &tim->dt_walk_self : nullptr);
        evaluate_pairs_post_drift<Spec>(ctx, actives.data(), N, cand_modeB, accums.data());
    }

    /* Host writeback — same code path as Mode A's writeback. */
    {
        StageTimer t(tim ? &tim->dt_writeback : nullptr);
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], accums[aa]);
        }
    }
    runner_emit_active_dumps<Spec>(args, accums.data(), actives.data(), "mode_b");
}

/* run_mode_b_local_with_oracle<Spec>: collects BOTH Mode B and Brute candidate
 * sets BEFORE drifting either, then drifts the union (drift_particle is
 * idempotent so two separate lazy_drift calls dedupe naturally).
 *
 * Mode B path is the result; Brute is the comparison. Mismatches print the
 * line shape `[mode_b ORACLE MISMATCH ...]` (parser-stable across the runner
 * extraction).
 *
 * Brute oracle uses ONLY Mode B / Brute machinery — no GPU NGL, no SIDX, no
 * arena. Codex invariant: "the oracle's pair kernel is the SAME pair_kernel
 * Mode B runs; only the search differs." */
/* Mismatch print cap — shared between local and remote oracle paths so we
 * don't double the cap when both fire (e.g. 3c.3 multi-rank-with-oracle
 * runs the local path under NTask==1 dispatch and the remote path under
 * NTask>1). */
static const long long kMismatchPrintCap = 1024;

template <typename Spec>
static void emit_oracle_mismatch_if_any(int rank, int active_slot,
                                          const typename Spec::AccumData& tree,
                                          const typename Spec::AccumData& brute,
                                          const char *origin_tag,
                                          long long *print_count)
{
    double resid = Spec::compare_accum(tree, brute);
    if(!(resid <= Spec::accum_tolerance)) {  /* NaN-safe */
        if(*print_count < kMismatchPrintCap) {
            fprintf(stderr,
                    "[mode_b ORACLE MISMATCH rank=%d caller=%s origin=%s "
                    "active_slot=%d resid=%g (tol=%g)]\n",
                    rank, Spec::loop_name, origin_tag, active_slot, resid,
                    (double)Spec::accum_tolerance);
            fflush(stderr);
            (*print_count)++;
            if(*print_count == kMismatchPrintCap) {
                fprintf(stderr,
                        "[mode_b ORACLE rank=%d caller=%s] mismatch print cap "
                        "reached (%lld); suppressing further mismatches.\n",
                        rank, Spec::loop_name, (long long)kMismatchPrintCap);
                fflush(stderr);
            }
        }
    }
}

/* ============================================================================
 * Iterative-oracle per-iter mismatch emit helpers (step 2c.4).
 *
 * Codex v2 P2: per-iter compare must cover four things, not just
 * AccumData. These four helpers complement emit_oracle_mismatch_if_any
 * (which handles AccumData via Spec::compare_accum) with IterStatus
 * control-flow / radius / active-set-membership compare lines. All share
 * the kMismatchPrintCap counter via the caller-passed long long*.
 * ========================================================================== */
template <typename Spec>
static void emit_oracle_status_mismatch(int rank, int sg, int slot, int iter,
                                          int prod_status, int oracle_status,
                                          long long *print_count)
{
    if (*print_count < kMismatchPrintCap) {
        fprintf(stderr,
                "[mode_b ORACLE STATUS MISMATCH rank=%d caller=%s sg=%d slot=%d "
                "iter=%d prod=%d oracle=%d]\n",
                rank, Spec::loop_name, sg, slot, iter, prod_status, oracle_status);
        fflush(stderr);
        (*print_count)++;
        if (*print_count == kMismatchPrintCap) {
            fprintf(stderr,
                    "[mode_b ORACLE rank=%d caller=%s] mismatch print cap "
                    "reached (%lld); suppressing further mismatches.\n",
                    rank, Spec::loop_name, (long long)kMismatchPrintCap);
            fflush(stderr);
        }
    }
}

template <typename Spec>
static void emit_oracle_radius_mismatch(int rank, int sg, int slot, int iter,
                                          double h_prod, double h_oracle,
                                          long long *print_count)
{
    if (*print_count < kMismatchPrintCap) {
        fprintf(stderr,
                "[mode_b ORACLE RADIUS MISMATCH rank=%d caller=%s sg=%d slot=%d "
                "iter=%d prod=%g oracle=%g (tol=%g)]\n",
                rank, Spec::loop_name, sg, slot, iter, h_prod, h_oracle,
                nlr_spec_radius_tolerance_v<Spec>);
        fflush(stderr);
        (*print_count)++;
    }
}

template <typename Spec>
static void emit_oracle_membership_count_mismatch(int rank, int sg, int iter,
                                                    int n_prod, int n_oracle,
                                                    long long *print_count)
{
    if (*print_count < kMismatchPrintCap) {
        fprintf(stderr,
                "[mode_b ORACLE MEMBERSHIP COUNT MISMATCH rank=%d caller=%s sg=%d "
                "iter=%d prod=%d oracle=%d]\n",
                rank, Spec::loop_name, sg, iter, n_prod, n_oracle);
        fflush(stderr);
        (*print_count)++;
    }
}

template <typename Spec>
static void emit_oracle_membership_slot_mismatch(int rank, int sg, int iter,
                                                   int k,
                                                   int slot_prod, int slot_oracle,
                                                   long long *print_count)
{
    if (*print_count < kMismatchPrintCap) {
        fprintf(stderr,
                "[mode_b ORACLE MEMBERSHIP SLOT MISMATCH rank=%d caller=%s sg=%d "
                "iter=%d k=%d prod=%d oracle=%d]\n",
                rank, Spec::loop_name, sg, iter, k, slot_prod, slot_oracle);
        fflush(stderr);
        (*print_count)++;
    }
}

template <typename Spec>
static void run_mode_b_local_with_oracle(const neighbor_loop_args& args, const double *radii,
                                         RunnerStageTimer *tim = nullptr)
{
    using ActiveData = typename Spec::ActiveData;
    using AccumData  = typename Spec::AccumData;
    using DeviceCtx  = typename Spec::DeviceContext;

    const int N = args.num_active;
    if(N <= 0) { return; }

    /* Radii runner-staged and passed in; pointer is call-lifetime only. */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    DeviceCtx ctx;
    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; the runner captures it by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);

    /* Freeze actives BEFORE drift (same snapshot for tree and brute). */
    std::vector<ActiveData> actives(N);
    build_self_actives_host_pre_drift<Spec>(args, ctx, radii, cs, actives.data());

    /* Collect BOTH BEFORE any drift — preserves identical pre-drift state for
     * each search backend. PHASE0 timing covers ONLY the Mode B path stages,
     * not the Brute oracle (diagnostic-only, not production cost). */
    std::vector<std::vector<int>> cand_modeB, cand_brute;
    {
        StageTimer t(tim ? &tim->dt_collect : nullptr);
        collect_candidates_pre_drift<Spec>(args, radii,
                                            (unsigned int)Spec::neighbor_type_mask,
                                            DispatchPath::ModeB_HostWalker, cand_modeB);
    }
    collect_candidates_pre_drift<Spec>(args, radii,
                                        (unsigned int)Spec::neighbor_type_mask,
                                        DispatchPath::Brute_Oracle, cand_brute);

    /* Drift the union. drift_particle's early-return on time1==time0 means
     * order doesn't matter and duplicates are free. */
    {
        StageTimer t(tim ? &tim->dt_drift : nullptr);
        lazy_drift_candidates<Spec>(cand_modeB);
    }
    lazy_drift_candidates<Spec>(cand_brute);

    std::vector<AccumData> accums_modeB(N);
    std::vector<AccumData> accums_brute(N);
    /* ORDERING: for j-write Specs, the brute oracle pass must run BEFORE the
     * tree pass so that brute reads the same pre-mutation j-state the tree
     * pass starts from. With the previous (tree-first) order, dry-run blocked
     * the brute's writes but did not restore the post-tree state, so brute
     * read e.g. zeroed Mass / changed SwallowID and produced a different
     * accumulator than tree -- a silent false-pass on the oracle whenever no
     * actual j-writes happened in the test config (codex review 2026-05-10).
     * For non-j-write Specs the order is irrelevant; we use brute-first
     * universally for one consistent oracle shape. */
    DeviceCtx ctx_oracle = ctx;
    if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
        Spec::set_oracle_brute_pass(ctx_oracle, true);
    }
    evaluate_pairs_post_drift<Spec>(ctx_oracle, actives.data(), N, cand_brute, accums_brute.data());
    {
        StageTimer t(tim ? &tim->dt_walk_self : nullptr);
        evaluate_pairs_post_drift<Spec>(ctx, actives.data(), N, cand_modeB, accums_modeB.data());
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    static long long s_mismatch_print_count = 0;
    for(int aa = 0; aa < N; aa++) {
        emit_oracle_mismatch_if_any<Spec>(rank, aa, accums_modeB[aa],
                                            accums_brute[aa], "self",
                                            &s_mismatch_print_count);
    }

    /* Mode B path is the result. */
    {
        StageTimer t(tim ? &tim->dt_writeback : nullptr);
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], accums_modeB[aa]);
        }
    }
    runner_emit_active_dumps<Spec>(args, accums_modeB.data(), actives.data(), "mode_b");
}

/* ============================================================================
 * Mode B remote (multi-rank) helpers (3c.3)
 *
 * STRONG INVARIANT (codex-locked, both oracle and non-oracle paths):
 *   collect-all → drift-union → evaluate-all on each rank.
 *
 * On each rank, the local pool is walked TWICE per call: once for this
 * rank's own self-pair queries (its own actives), and once for queries
 * received from peers. Both candidate sets must be collected pre-drift,
 * then the drift covers their UNION (idempotent so duplicates are free),
 * then evaluation runs post-drift. The non-oracle path does NOT simplify
 * back to per-query drift — single shared structure, oracle off only
 * skips the brute candidate sets and the compare step.
 *
 * Sequence:
 *   stage 1 (active rank) build self radii, cs, frozen actives[]
 *   stage 2 (active rank) build envelopes, ALL peers in broadcast pattern
 *   stage 3 (this rank)   collect self_tree[, self_brute] pre-drift
 *   stage 4 (collective)  exchange queries (peer-to-peer) -> recv envelopes
 *   stage 5 (this rank)   flatten envelopes to peer_actives[] + provenance[]
 *   stage 6 (this rank)   collect peer_tree[, peer_brute] pre-drift
 *   stage 7 (this rank)   drift UNION of (self_tree, self_brute, peer_tree,
 *                                          peer_brute)
 *   stage 8 (this rank)   evaluate self_tree[, self_brute]; oracle compare
 *                                  self before any remote merge so prints
 *                                  isolate local-pool walker bugs
 *   stage 9 (this rank)   evaluate peer_tree[, peer_brute]; oracle compare
 *                                  peer-side BEFORE reply transport (peer-
 *                                  rank prints isolate peer-pool bugs)
 *   stage 10 (collective) exchange replies (tree result is what ships)
 *   stage 11 (active rank) merge replies via Spec::merge_accum, ascending
 *                                  rank for FP-reproducible order
 *   stage 12 (active rank) writeback
 *
 * Codex nuance (acknowledged): the host-frozen actives[] match Mode B /
 * legacy pack_query epoch, NOT Mode A's device-staged post-NGL-build
 * epoch. Oracle in this path is Mode B tree vs brute on the same frozen
 * query — it does NOT cross-validate against Mode A's active epoch.
 * ========================================================================== */

/* ============================================================================
 * mode_b_remote_evaluate_into_buffer<Spec, ORACLE> — extracted helper (step 2b.2).
 *
 * Mechanical refactor of the existing run_mode_b_remote_impl body, extracting
 * stages 1-12 (queries / collect / drift / evaluate / merge / replies / merge)
 * + the env-gated active_dumps diagnostic into a reusable helper. The caller
 * provides the output AccumData buffer; the helper does NOT call
 * Spec::apply_active_writeback. Final-only writeback decision stays with the
 * caller (codex 2b.2 constraint 1: "let the caller decide whether to call
 * apply_active_writeback").
 *
 * Two callers:
 *   - run_mode_b_remote_impl (non-iterative wrapper): allocates per-call
 *     AccumData buffer, calls helper (which emits dumps internally on its
 *     locally-built actives[] view), then runs apply_active_writeback.
 *   - nlr_iter_dispatch_subgroup_mode_b_remote (iterative dispatch): calls
 *     helper with driver-owned compacted AccumData buffer; runner does
 *     final-only apply_active_writeback after the iter loop. (Per-iter
 *     active_dumps emits are harmless: env-gated, only fire under spike-test
 *     flags.)
 *
 * Epoch order is preserved EXACTLY (codex constraint 2): build self/peer
 * queries -> collect all candidate sets pre-drift -> exchange -> drift union
 * -> brute-FIRST if oracle -> evaluate -> exchange replies -> merge replies
 * in deterministic peer order. No re-derivation; line-for-line move from
 * the old impl.
 * ========================================================================== */
template <typename Spec, RemoteHelperMode MODE>
static void mode_b_remote_evaluate_into_buffer(
    const neighbor_loop_args& args,
    const double *radii,
    const typename Spec::CallScalars& cs,
    const typename Spec::DeviceContext& ctx,         /* caller-owned (step 2c.1) */
    unsigned int neighbor_type_mask,                  /* explicit caller param (step 2c.3) */
    typename Spec::AccumData *accums_out,             /* size = args.num_active; caller-owned */
    typename Spec::ActiveData *actives_out = nullptr, /* size = args.num_active OR nullptr */
    RunnerStageTimer *tim = nullptr,
    typename Spec::AccumData *accums_oracle_out = nullptr, /* OracleIterative only: brute accum output */
    const typename Spec::DeviceContext *ctx_oracle = nullptr) /* OracleIterative brute context */
{
    using ActiveData    = typename Spec::ActiveData;
    using AccumData     = typename Spec::AccumData;
    using DeviceCtx     = typename Spec::DeviceContext;     /* needed for oracle ctx copies (Stages 8/9) */
    using Envelope      = NlrQueryEnvelope<ActiveData>;
    using ReplyEnvelope = NlrReplyEnvelope<AccumData>;
    using DualReplyEnvelope = NlrDualReplyEnvelope<AccumData>;

    /* Mode predicates (step 2c.4 SSOT extension). Constexpr so the dead
     * branches are elided per instantiation. */
    constexpr bool RUN_TREE         = (MODE == RemoteHelperMode::Production ||
                                        MODE == RemoteHelperMode::OracleCompare ||
                                        MODE == RemoteHelperMode::OracleIterative);
    constexpr bool RUN_BRUTE        = (MODE == RemoteHelperMode::OracleCompare ||
                                        MODE == RemoteHelperMode::OracleBrutePass ||
                                        MODE == RemoteHelperMode::OracleIterative);
    constexpr bool RUN_INLINE_COMPARE = (MODE == RemoteHelperMode::OracleCompare);
    constexpr bool BRUTE_WRITES_OUT = (MODE == RemoteHelperMode::OracleBrutePass);
    /* OracleIterative: collect both candidate sets in one epoch, output tree
     * -> accums_out and brute -> accums_oracle_out separately (no inline compare). */
    constexpr bool DUAL_OUT         = (MODE == RemoteHelperMode::OracleIterative);

    static_assert(std::is_trivially_copyable<Envelope>::value,
        "NlrQueryEnvelope must be trivially-copyable for byte-level MPI transfer");
    static_assert(std::is_trivially_copyable<ReplyEnvelope>::value,
        "NlrReplyEnvelope must be trivially-copyable for byte-level MPI transfer");
    static_assert(std::is_trivially_copyable<DualReplyEnvelope>::value,
        "NlrDualReplyEnvelope must be trivially-copyable for byte-level MPI transfer");

    const int N    = args.num_active;     /* may be 0 on this rank; collective entry */
    const int nt   = NTask;
    const int rank = ThisTask;

    /* CallScalars and DeviceContext passed in by caller (step 2c.1):
     *   - Iterative: driver-owned via NlrIterDriver, populated once at iter-0 entry.
     *   - Non-iterative wrapper: locally-owned, populated + RAII-cleaned by wrapper.
     * Helper does NOT call populate_call_scalars or populate_device_context. */

    /* Freeze actives host-side pre-drift (same snapshot used for self-pair
     * AND for ship-to-peers; codex requirement to prevent self/remote epoch
     * skew on the active rank). */
    std::vector<ActiveData> actives(N);
    if(N > 0) {
        build_self_actives_host_pre_drift<Spec>(args, ctx, radii, cs,
                                                  actives.data());
    }

    /* Stage 2: build envelopes per peer (broadcast pattern: every peer
     * receives ALL of this rank's actives). Self-pair handled locally;
     * envelopes for self entry stay empty. */
    std::vector<std::vector<Envelope>> queries_per_peer(nt);
    if(N > 0) {
        for(int p = 0; p < nt; p++) {
            if(p == rank) continue;
            queries_per_peer[p].reserve(N);
            for(int aa = 0; aa < N; aa++) {
                Envelope env;
                env.origin_slot = aa;
                env.origin_rank = rank;
                env.active      = actives[aa];
                queries_per_peer[p].push_back(env);
            }
        }
    }

    /* Stage 3: collect SELF candidate sets PRE-DRIFT. */
    std::vector<std::vector<int>> cand_self_tree, cand_self_brute;
    if(N > 0) {
        if (RUN_TREE) {
            StageTimer t(tim ? &tim->dt_collect : nullptr);
            collect_candidates_pre_drift<Spec>(args, radii,
                                                neighbor_type_mask,
                                                DispatchPath::ModeB_HostWalker,
                                                cand_self_tree);
        }
        if (RUN_BRUTE) {
            collect_candidates_pre_drift<Spec>(args, radii,
                                                neighbor_type_mask,
                                                DispatchPath::Brute_Oracle,
                                                cand_self_brute);
        }
    }

    /* Stage 4: exchange queries (collective). Every rank participates even
     * if N == 0 (peers may have queries directed at this rank's pool). */
    auto state = [&]{
        StageTimer t(tim ? &tim->dt_exchange_q : nullptr);
        return mode_b_exchange_queries<Envelope>(queries_per_peer);
    }();

    /* Stage 5: flatten received envelopes and build provenance map.
     * provenance[k] carries:
     *   - source_peer / source_qi: where in state.recv_queries[][] this k
     *     came from (used for unflattening replies back into per-peer arrays).
     *   - origin_slot / origin_rank: copied from the received envelope; ride
     *     into the REPLY envelope so the active rank can merge by slot
     *     without relying on transport ordering (codex Option C — symmetric
     *     query and reply envelopes). */
    std::vector<ActiveData> peer_actives;
    struct Provenance {
        int source_peer; int source_qi;
        int origin_slot; int origin_rank;
    };
    std::vector<Provenance> peer_provenance;
    int total_recv = 0;
    for(int p = 0; p < nt; p++) total_recv += state.recv_counts[p];
    peer_actives.reserve(total_recv);
    peer_provenance.reserve(total_recv);
    for(int p = 0; p < nt; p++) {
        if(p == rank) continue;
        for(int qi = 0; qi < state.recv_counts[p]; qi++) {
            const Envelope& env = state.recv_queries[p][qi];
            /* Sanity: envelope's origin_rank should equal sender p. */
            if(env.origin_rank != p) {
                fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                        "envelope origin_rank=%d but received from peer=%d "
                        "(qi=%d). Transport corruption?\n",
                        rank, Spec::loop_name, env.origin_rank, p, qi);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            peer_actives.push_back(env.active);
            peer_provenance.push_back({p, qi, env.origin_slot, env.origin_rank});
        }
    }

    /* Stage 6: collect PEER candidate sets PRE-DRIFT (against MY local pool). */
    std::vector<std::vector<int>> cand_peer_tree, cand_peer_brute;
    if (RUN_TREE) {
        StageTimer t(tim ? &tim->dt_collect : nullptr);
        collect_candidates_for_remote_queries<Spec>(peer_actives,
                                                     neighbor_type_mask,
                                                     DispatchPath::ModeB_HostWalker,
                                                     cand_peer_tree);
    }
    if (RUN_BRUTE) {
        collect_candidates_for_remote_queries<Spec>(peer_actives,
                                                     neighbor_type_mask,
                                                     DispatchPath::Brute_Oracle,
                                                     cand_peer_brute);
    }

    /* Stage 7: drift the UNION of all candidate sets that touch MY pool. */
    {
        StageTimer t(tim ? &tim->dt_drift : nullptr);
        if (RUN_TREE) {
            if(N > 0) lazy_drift_candidates<Spec>(cand_self_tree);
            lazy_drift_candidates<Spec>(cand_peer_tree);
        }
    }
    if (RUN_BRUTE) {
        if(N > 0) lazy_drift_candidates<Spec>(cand_self_brute);
        lazy_drift_candidates<Spec>(cand_peer_brute);
    }

    /* Stage 8: evaluate SELF post-drift.
     *   Production:       tree -> accums_out.
     *   OracleCompare:    brute dry-run (own buffer) -> tree -> accums_out,
     *                     emit per-slot inline compare.
     *   OracleBrutePass:  brute -> accums_out (caller-owned ctx already
     *                     brute-pass-guarded). No tree, no compare.
     *
     * Brute-FIRST ordering (codex 3d.3 / 2026-05-10): dry-run BEFORE tree
     * so brute reads pre-mutation j-state; without it, tree's j-writes
     * leak into brute's read. For OracleBrutePass the caller owns the
     * brute-pass guard at outer scope (NlrOracleBrutePassGuard at iter
     * dispatch helper entry); helper does NOT touch set_oracle_brute_pass.
     * For OracleCompare the helper still owns a local ctx copy + inline
     * toggle for the dry-run pass (legacy shape preserved). */
    std::vector<AccumData> accums_self_brute;
    if(N > 0) {
        if constexpr (RUN_INLINE_COMPARE) {
            accums_self_brute.assign(N, AccumData{});
            DeviceCtx ctx_oracle_self = ctx;
            if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
                Spec::set_oracle_brute_pass(ctx_oracle_self, true);
            }
            evaluate_pairs_post_drift<Spec>(ctx_oracle_self, actives.data(), N,
                                              cand_self_brute, accums_self_brute.data());
        }
        if constexpr (BRUTE_WRITES_OUT) {
            /* Caller's ctx is already brute-pass-guarded via
             * NlrOracleBrutePassGuard at the iter dispatch helper. Walk
             * directly into accums_out. */
            StageTimer t(tim ? &tim->dt_walk_self : nullptr);
            evaluate_pairs_post_drift<Spec>(ctx, actives.data(), N,
                                              cand_self_brute, accums_out);
        }
        if constexpr (DUAL_OUT) {
            /* OracleIterative: brute-first into accums_oracle_out (same epoch
             * as tree; both candidate sets were collected pre-drift above). */
            DeviceCtx ctx_oracle_dual = (ctx_oracle != nullptr) ? *ctx_oracle : ctx;
            if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
                Spec::set_oracle_brute_pass(ctx_oracle_dual, true);
            }
            evaluate_pairs_post_drift<Spec>(ctx_oracle_dual, actives.data(), N,
                                              cand_self_brute, accums_oracle_out);
        }
        if constexpr (RUN_TREE) {
            StageTimer t(tim ? &tim->dt_walk_self : nullptr);
            evaluate_pairs_post_drift<Spec>(ctx, actives.data(), N,
                                              cand_self_tree, accums_out);
        }
        if constexpr (RUN_INLINE_COMPARE) {
            static long long s_self_mismatch_count = 0;
            for(int aa = 0; aa < N; aa++) {
                emit_oracle_mismatch_if_any<Spec>(rank, aa, accums_out[aa],
                                                    accums_self_brute[aa], "self",
                                                    &s_self_mismatch_count);
            }
        }
    }

    /* Stage 9: evaluate PEER queries post-drift.
     *   Production / OracleCompare: tree result -> peer_replies, shipped
     *     back to home rank.
     *   OracleBrutePass:            brute result -> peer_replies, shipped
     *     back to home rank — peer-side brute trajectory for iterative
     *     oracle.
     *   OracleIterative:            brute -> peer_replies_oracle; tree ->
     *     peer_replies; both are exchanged together in one dual payload. */
    const int K = (int)peer_actives.size();
    std::vector<AccumData> peer_replies(K);
    std::vector<AccumData> peer_replies_brute;
    std::vector<AccumData> peer_replies_oracle;
    if(K > 0) {
        if constexpr (RUN_INLINE_COMPARE) {
            peer_replies_brute.assign(K, AccumData{});
            DeviceCtx ctx_oracle_peer = ctx;
            if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
                Spec::set_oracle_brute_pass(ctx_oracle_peer, true);
            }
            evaluate_pairs_post_drift<Spec>(ctx_oracle_peer, peer_actives.data(), K,
                                              cand_peer_brute, peer_replies_brute.data());
        }
        if constexpr (BRUTE_WRITES_OUT) {
            /* Caller's ctx is already brute-pass-guarded; peer-side brute
             * replies fill peer_replies directly (will then exchange and
             * merge into accums_out alongside the self-side brute result). */
            StageTimer t(tim ? &tim->dt_walk_peer : nullptr);
            evaluate_pairs_post_drift<Spec>(ctx, peer_actives.data(), K,
                                              cand_peer_brute, peer_replies.data());
        }
        if constexpr (DUAL_OUT) {
            /* OracleIterative: peer brute -> peer_replies_oracle, peer tree
             * -> peer_replies. Brute first for j-write specs. */
            peer_replies_oracle.assign(K, AccumData{});
            DeviceCtx ctx_oracle_peer_dual = (ctx_oracle != nullptr) ? *ctx_oracle : ctx;
            if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
                Spec::set_oracle_brute_pass(ctx_oracle_peer_dual, true);
            }
            evaluate_pairs_post_drift<Spec>(ctx_oracle_peer_dual, peer_actives.data(), K,
                                              cand_peer_brute, peer_replies_oracle.data());
        }
        if constexpr (RUN_TREE) {
            StageTimer t(tim ? &tim->dt_walk_peer : nullptr);
            evaluate_pairs_post_drift<Spec>(ctx, peer_actives.data(), K,
                                              cand_peer_tree, peer_replies.data());
        }
        if constexpr (RUN_INLINE_COMPARE) {
            static long long s_peer_mismatch_count = 0;
            for(int k = 0; k < K; k++) {
                emit_oracle_mismatch_if_any<Spec>(rank, peer_provenance[k].origin_slot,
                                                    peer_replies[k],
                                                    peer_replies_brute[k],
                                                    "peer", &s_peer_mismatch_count);
            }
        }
    }

    if constexpr (DUAL_OUT) {
        /* Iterative oracle is temporary port-validation scaffolding. Keep the
         * transport simple and correct: production and brute replies travel in
         * one payload, so we do not run two reply exchanges with identical MPI
         * tags. */
        std::vector<std::vector<DualReplyEnvelope>> replies_per_peer(nt);
        for (int p = 0; p < nt; p++) {
            replies_per_peer[p].assign(state.recv_counts[p], DualReplyEnvelope{});
        }
        for (int k = 0; k < K; k++) {
            const Provenance& pv = peer_provenance[k];
            DualReplyEnvelope& re = replies_per_peer[pv.source_peer][pv.source_qi];
            re.origin_slot  = pv.origin_slot;
            re.origin_rank  = pv.origin_rank;
            re.accum_prod   = peer_replies[k];
            re.accum_oracle = peer_replies_oracle[k];
        }
        auto recv_replies = [&]{
            StageTimer t(tim ? &tim->dt_exchange_r : nullptr);
            return mode_b_exchange_replies<Envelope, DualReplyEnvelope>(replies_per_peer, state);
        }();
        {
            StageTimer t(tim ? &tim->dt_reduce : nullptr);
            if (N > 0 && accums_oracle_out != nullptr) {
                for (int p = 0; p < nt; p++) {
                    if (p == rank) continue;
                    const int q_to_p = state.sent_counts[p];
                    for (int qi = 0; qi < q_to_p; qi++) {
                        const DualReplyEnvelope& re = recv_replies[p][qi];
                        if (re.origin_rank != rank) {
                            fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                    "dual reply envelope origin_rank=%d != ThisTask=%d from "
                                    "peer=%d qi=%d. Transport/peer-side corruption?\n",
                                    rank, Spec::loop_name, re.origin_rank, rank, p, qi);
                            fflush(stderr);
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                        const int slot = re.origin_slot;
                        if (slot < 0 || slot >= N) {
                            fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                    "dual reply envelope slot %d out of range [0,%d) from peer %d.\n",
                                    rank, Spec::loop_name, slot, N, p);
                            fflush(stderr);
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                        Spec::merge_accum(accums_out[slot], re.accum_prod);
                        Spec::merge_accum(accums_oracle_out[slot], re.accum_oracle);
                    }
                }
            }
        }
    } else {
        /* Stage 10: build reply envelopes (origin_slot/rank copied from each
         * received query envelope), unflatten into per-peer arrays via the
         * provenance map, then exchange. Reply envelope makes the active-side
         * merge transport-order independent (codex Option C consistency). */
        std::vector<std::vector<ReplyEnvelope>> replies_per_peer(nt);
        for(int p = 0; p < nt; p++) {
            replies_per_peer[p].assign(state.recv_counts[p], ReplyEnvelope{});
        }
        for(int k = 0; k < K; k++) {
            const Provenance& pv = peer_provenance[k];
            ReplyEnvelope& re = replies_per_peer[pv.source_peer][pv.source_qi];
            re.origin_slot = pv.origin_slot;
            re.origin_rank = pv.origin_rank;
            re.accum       = peer_replies[k];
        }

        auto recv_replies = [&]{
            StageTimer t(tim ? &tim->dt_exchange_r : nullptr);
            return mode_b_exchange_replies<Envelope, ReplyEnvelope>(replies_per_peer, state);
        }();

        /* Stage 11: merge replies into accums_self by envelope.origin_slot.
         * Pinned deterministic order: ascending peer rank (self contribution
         * already in accums_self from stage 8). Asserts each reply envelope's
         * origin_rank == ThisTask — a transport-corruption sanity check. */
        {
            StageTimer t(tim ? &tim->dt_reduce : nullptr);
            for (int p = 0; p < nt; p++) {
                if (p == rank) continue;
                const int q_to_p = state.sent_counts[p];
                for (int qi = 0; qi < q_to_p; qi++) {
                    const ReplyEnvelope& re = recv_replies[p][qi];
                    if(re.origin_rank != rank) {
                        fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                "reply envelope origin_rank=%d != ThisTask=%d from "
                                "peer=%d qi=%d. Transport/peer-side corruption?\n",
                                rank, Spec::loop_name, re.origin_rank, rank, p, qi);
                        fflush(stderr);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    const int slot = re.origin_slot;
                    if(slot < 0 || slot >= N) {
                        fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                "reply envelope slot %d out of range [0,%d) from peer %d.\n",
                                rank, Spec::loop_name, slot, N, p);
                        fflush(stderr);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    Spec::merge_accum(accums_out[slot], re.accum);
                }
            }
        }
    }

    /* Optionally export the actives[] snapshot to caller for post-writeback
     * diagnostic emit (codex 2c.1 review fix: non-iter dump ordering must be
     * apply_active_writeback FIRST, then runner_emit_active_dumps — pre-2b.2
     * timing preserved by exporting actives + letting wrapper emit AFTER its
     * writeback loop). Iter dispatch caller passes nullptr — diagnostic
     * plumbing for the iterative path lands in a later slice. */
    if (actives_out != nullptr && N > 0) {
        for (int aa = 0; aa < N; aa++) actives_out[aa] = actives[aa];
    }

    /* End of helper. Caller decides whether to call apply_active_writeback
     * (final-only for iterative; per-call for non-iterative) and whether
     * to emit active_dumps (after writeback for non-iter; deferred for
     * iter until proper diagnostic plumbing lands). */
}

/* ============================================================================
 * run_mode_b_remote_impl<Spec, ORACLE> — non-iterative thin wrapper (step 2b.2).
 *
 * Allocates per-call AccumData buffer, calls helper, runs final
 * apply_active_writeback + active_dumps emit. Behavior is byte-identical to
 * the pre-2b.2 monolithic impl — same epoch order, same oracle handling,
 * same writeback per active.
 * ========================================================================== */
template <typename Spec, bool ORACLE>
static void run_mode_b_remote_impl(const neighbor_loop_args& args, const double *radii,
                                   RunnerStageTimer *tim = nullptr)
{
    using AccumData    = typename Spec::AccumData;
    using ActiveData   = typename Spec::ActiveData;
    using DeviceCtx    = typename Spec::DeviceContext;

    const int N = args.num_active;

    /* Capture CallScalars per call (non-iter wrapper). */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    /* Build per-call DeviceContext + RAII cleanup guard (step 2c.1: helper
     * no longer builds its own ctx; caller does). Behavior byte-equivalent
     * to pre-2c.1 monolith — guard runs Spec::cleanup_device_context at
     * scope exit, after the writeback loop completes. */
    DeviceCtx ctx;
    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);

    /* Caller-owned per-call AccumData buffer; helper writes into it.
     * Explicit Spec::zero_accum per slot (codex 2b.2 review fix): defensive
     * against future evaluate_pairs variants that don't zero internally.
     * Today evaluate_pairs_post_drift calls Spec::zero_accum at line 812
     * before walking candidates, so this outer zero is redundant — but
     * makes the AccumData contract explicit at the caller level (sentinel-
     * bearing AccumData like sink_feed's Sink_PotentialMinimumOfNeighbors
     * needs Spec::zero_accum, NOT default-construction). */
    std::vector<AccumData> accums_self(N);
    for (int aa = 0; aa < N; aa++) {
        Spec::zero_accum(accums_self[aa]);
    }

    /* Caller-owned actives buffer for post-writeback diagnostic emit.
     * Helper fills this when actives_out != nullptr; pre-2b.2 emit ordering
     * (apply_active_writeback FIRST, runner_emit_active_dumps SECOND) is
     * preserved by emitting from this buffer AFTER the writeback loop below. */
    std::vector<ActiveData> actives_for_dumps(N);

    /* Helper runs Stages 1-12; writeback + dumps are the wrapper's
     * responsibility (preserves pre-2b.2 timing per codex 2c.1 review). */
    /* Step 2c.4 SSOT extension: helper now takes RemoteEvalMode enum
     * (Production / OracleCompare / OracleBrutePass). Non-iter wrapper
     * keeps the legacy bool ORACLE template parameter and translates. */
    constexpr RemoteHelperMode MODE = ORACLE ? RemoteHelperMode::OracleCompare
                                            : RemoteHelperMode::Production;
    mode_b_remote_evaluate_into_buffer<Spec, MODE>(args, radii, cs, ctx,
                                                       (unsigned int)Spec::neighbor_type_mask,
                                                       (N > 0) ? accums_self.data() : nullptr,
                                                       (N > 0) ? actives_for_dumps.data() : nullptr,
                                                       tim);

    /* Stage 12 final: writeback per active. */
    {
        StageTimer t(tim ? &tim->dt_writeback : nullptr);
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], accums_self[aa]);
        }
    }

    /* Post-writeback diagnostic emit (env-gated; pre-2b.2 ordering restored). */
    if (N > 0) {
        runner_emit_active_dumps<Spec>(args, accums_self.data(),
                                         actives_for_dumps.data(), "mode_b");
    }
}

template <typename Spec>
static void run_mode_b_remote(const neighbor_loop_args& args, const double *radii,
                              RunnerStageTimer *tim = nullptr) {
    run_mode_b_remote_impl<Spec, /*ORACLE=*/false>(args, radii, tim);
}
template <typename Spec>
static void run_mode_b_remote_with_oracle(const neighbor_loop_args& args, const double *radii,
                                          RunnerStageTimer *tim = nullptr) {
    run_mode_b_remote_impl<Spec, /*ORACLE=*/true>(args, radii, tim);
}

/* ============================================================================
 * run_mode_a<Spec> — generic Mode A path through the GPU NGL pipeline.
 *
 * Three-epoch staging contract (see neighbor_loop_runner.h doc):
 *   (1) host pre-arena: Spec::search_radius  → radii_uvm[num_active]
 *   (2) host pre-arena: Spec::populate_call_scalars → CallScalars cs
 *   (3) device post-NGL-build: Spec::load_active → d_actives[num_active]
 *
 * Both Kokkos launches go through gizmo_gpu_kernel_launch which wraps
 * parallel_for + fence + check_last_error (declarations/gpu_dispatch_templates.h).
 * No additional explicit fence is needed before host-side d_accums readback
 * — the launch helper already fenced.
 * ========================================================================== */

template <typename Spec>
static void run_mode_a(const neighbor_loop_args& args, const double *radii,
                       RunnerStageTimer *tim = nullptr)
{
    using ActiveData   = typename Spec::ActiveData;
    using AccumData    = typename Spec::AccumData;
    using ScatterData  = typename Spec::ScatterData;
    using NeighborData = typename Spec::NeighborData;
    using CallScalars  = typename Spec::CallScalars;
    using DeviceCtx    = typename Spec::DeviceContext;

    const int N = args.num_active;
    if(N <= 0) {
        /* Defensive: caller already early-outs in collective dispatch when
         * global_num_active == 0, but local-zero with global-positive is
         * legitimate (peer rank with no actives still must enter the
         * collective in Mode B; Mode A has no collective work, just no-op). */
        return;
    }

    /* (1) Host, pre-arena: stage caller-supplied radii into UVM for the
     * device-visible NGL build. Source `radii` is runner-staged on host
     * (call-lifetime only); we copy into shared/UVM so gpu_ngb_list_build
     * and the device pair kernel can read the per-active values. */
    double *radii_uvm = (double *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        N * sizeof(double));
    for(int aa = 0; aa < N; aa++) {
        radii_uvm[aa] = radii[aa];
    }

    /* (2) Host, pre-arena: capture per-call scalar globals into a POD. */
    CallScalars cs = Spec::populate_call_scalars(args);

    /* Arena + freshness (matches sinks/sink_environment_gpu.cc:76,86). The
     * caller is responsible for the args.CellP=NULL-when-no-gas decision;
     * runner does not read All.TotN_gas. */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total, args.P, args.CellP);
    struct particle_data *P_gpu    = gpu_particles_arena_P();
    struct gas_cell_data *CellP_gpu = (args.CellP != nullptr)
                                        ? gpu_particles_arena_CellP() : nullptr;

    /* NGL build using pre-arena radii. SIDX cache resolved from spec. */
    gpu_neighbor_list_t gnl;
    gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                       Spec::loop_name);
    {
        StageTimer t(tim ? &tim->dt_collect : nullptr);
        gpu_ngb_list_build(P_gpu, args.num_total,
                           args.active_list, N,
                           Spec::search_mode,
                           (int)Spec::neighbor_type_mask,
                           &gnl, sidx,
                           1.0, radii_uvm, NULL, Spec::loop_name,
                           nlr_spec_symmetric_j_radius_scale<Spec>());
    }

    /* UVM-allocate ActiveData[] and AccumData[] arrays. */
    ActiveData *d_actives = (ActiveData *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(N * sizeof(ActiveData));
    AccumData *d_accums = (AccumData *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(N * sizeof(AccumData));

    /* Build DeviceContext. Specs that extend Spec::DeviceContext beyond
     * NeighborLoopDeviceContextBase get populate_device_context invoked
     * here (Phase 4.A.0); base-only Specs skip the call (trait check). */
    DeviceCtx ctx;
    ctx.P         = P_gpu;
    ctx.CellP     = CellP_gpu;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; the runner captures it by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);

    /* (3) Device, post-NGL-build: stage ActiveData. Same device epoch as
     * the legacy GPU lambda's q-packing. d_active is gnl-resident. */
    int *d_active_idx = gnl.d_active;
    {
        gizmo_gpu_kernel_launch("nlr_stage_active", N, KOKKOS_LAMBDA(int aa) {
            d_actives[aa] = Spec::load_active(ctx, aa, d_active_idx[aa],
                                              radii_uvm[aa], cs);
        });
    }

    /* Pair-kernel launch — generic over Spec. */
    {
        StageTimer t(tim ? &tim->dt_walk_self : nullptr);
        int64_t *offsets = gnl.offsets;
        int *neighbors = gnl.neighbors;
        gizmo_gpu_kernel_launch(Spec::loop_name, N, KOKKOS_LAMBDA(int aa) {
            Spec::zero_accum(d_accums[aa]);
            const ActiveData& a = d_actives[aa];
            ScatterData s{};                     /* NoScatter for ActiveReduceOnly */
            int64_t start = offsets[aa], end = offsets[aa + 1];
            for(int64_t nn = start; nn < end; nn++) {
                int j = neighbors[nn];
                IdentitySidecar id{};            /* NoIdentity */
                NeighborData nb = Spec::load_neighbor(ctx, j, id, a);
                Spec::pair_kernel(a, nb, d_accums[aa], s);
            }
        });
    }
    /* Both launches fenced internally by gizmo_gpu_kernel_launch. UVM is
     * coherent → host can read d_accums directly. */

    /* (4) Host writeback via spec. */
    {
        StageTimer t(tim ? &tim->dt_writeback : nullptr);
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], d_accums[aa]);
        }
    }

    /* (5) Optional diagnostic dumps (3c.4a). Order matches legacy whole-log
     * shape: NB lines first (legacy emitted from sink_environment_gpu.cc
     * BEFORE returning to sink_environment.cc, where the accumulator dump
     * fired), then accumulator lines. NB dump is Mode A only; the GPU CSR
     * host-copy is allocated only when its env gate is on (zero overhead
     * off). Accumulator dump preserves legacy line shape. */
    if constexpr (spec_has_dump_neighbor_list<Spec>::value) {
        if(gizmo_nlr_xval_nb_dump_enabled() && gnl.total_pairs > 0) {
            std::vector<int> nbrs_host((size_t)gnl.total_pairs);
            gpu_ngb_copy_neighbors_to_host(&gnl, nbrs_host.data());
            int64_t *offsets = gnl.offsets;
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            for(int aa = 0; aa < N; aa++) {
                NeighborListDumpView<Spec> v;
                v.rank          = rank;
                v.origin_rank   = rank;
                v.origin_slot   = aa;
                v.active_slot   = aa;
                v.path          = "gpu_ngl";
                v.call_id       = 1;
                v.args          = &args;
                v.active        = &d_actives[aa];
                v.candidate_ids = nbrs_host.data() + offsets[aa];
                /* per-active neighbor count is bounded by num_total < 2^31 */
                v.n_candidates  = (int)(offsets[aa + 1] - offsets[aa]);
                Spec::diagnostic_dump_neighbor_list(v);
            }
            std::fflush(stdout);
        }
    }
    runner_emit_active_dumps<Spec>(args, d_accums, d_actives, "gpu_ngl");

    /* Cleanup. SIDX cache pointer passed so the free leaves cached storage
     * intact for sink_feed/sink_swk reuse (matches existing
     * sink_environment_gpu.cc:261 idiom). */
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_actives);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
    gpu_ngb_list_free(&gnl, sidx);
    gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
}

/* ============================================================================
 * SFINAE: detect optional Spec::modeb_threshold_sum / _max constexpr.
 * If absent, fall back to caller-supplied default.
 * ========================================================================== */

template <typename Spec, typename = void>
struct nlr_has_threshold_sum : std::false_type {};
template <typename Spec>
struct nlr_has_threshold_sum<Spec, decltype((void)Spec::modeb_threshold_sum)>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_threshold_max : std::false_type {};
template <typename Spec>
struct nlr_has_threshold_max<Spec, decltype((void)Spec::modeb_threshold_max)>
    : std::true_type {};

template <typename Spec>
static int nlr_spec_threshold_sum(int fallback) {
    if constexpr (nlr_has_threshold_sum<Spec>::value) {
        return (int)Spec::modeb_threshold_sum;
    } else {
        return fallback;
    }
}
template <typename Spec>
static int nlr_spec_threshold_max(int fallback) {
    if constexpr (nlr_has_threshold_max<Spec>::value) {
        return (int)Spec::modeb_threshold_max;
    } else {
        return fallback;
    }
}

/* ============================================================================
 * SFINAE: detect optional Spec::uses_* lifecycle traits (default false) and
 * dispatch to the corresponding hook methods (default no-op). Two channels:
 *
 *   detector  : audit/debug — catches illegal kernel writes to imported
 *               ghosts. Trait `uses_ghost_write_detector`.
 *   writeback : physics state propagation — pre-kernel ghost-state snapshot
 *               + post-kernel reverse-comm of any j-side writes. Trait
 *               `uses_ghost_writeback`. Spec body is the per-flag #ifdef
 *               union of all enabled physics-flag writebacks.
 *
 * Each hook is gated on (a) the trait being true AND (b) whether the chosen
 * path imports ghosts (per-Spec audit decided imported-ghost-only for these
 * hooks; future Specs may differ).
 * ========================================================================== */

/* Trait detection: Spec::uses_ghost_write_detector / _ghost_writeback.
 * Absent ⇒ false. */
template <typename Spec, typename = void>
struct nlr_has_uses_ghost_write_detector : std::false_type {};
template <typename Spec>
struct nlr_has_uses_ghost_write_detector<Spec, decltype((void)Spec::uses_ghost_write_detector)>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_uses_ghost_writeback : std::false_type {};
template <typename Spec>
struct nlr_has_uses_ghost_writeback<Spec, decltype((void)Spec::uses_ghost_writeback)>
    : std::true_type {};

template <typename Spec> static constexpr bool nlr_uses_ghost_write_detector_v() {
    if constexpr (nlr_has_uses_ghost_write_detector<Spec>::value) {
        return Spec::uses_ghost_write_detector;
    } else { return false; }
}
template <typename Spec> static constexpr bool nlr_uses_ghost_writeback_v() {
    if constexpr (nlr_has_uses_ghost_writeback<Spec>::value) {
        return Spec::uses_ghost_writeback;
    } else { return false; }
}

/* Hook-method detection: Spec::ghost_write_detector_begin etc. Absent ⇒
 * runner skips the call (compile-time short-circuit; no runtime cost). */
template <typename Spec, typename = void>
struct nlr_has_hook_gwd_begin : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwd_begin<Spec,
    decltype(Spec::ghost_write_detector_begin(std::declval<const neighbor_loop_args&>(),
                                              std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwd_end : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwd_end<Spec,
    decltype(Spec::ghost_write_detector_end(std::declval<const neighbor_loop_args&>(),
                                             std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwb_begin : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwb_begin<Spec,
    decltype(Spec::ghost_writeback_begin(std::declval<const neighbor_loop_args&>(),
                                         std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwb_end : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwb_end<Spec,
    decltype(Spec::ghost_writeback_end(std::declval<const neighbor_loop_args&>(),
                                       std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

/* Dispatch wrappers. Each gates on:
 *   (a) `uses_*` trait true (Spec opted in)
 *   (b) `nlr_path_uses_imported_ghosts(plan.path)` — for SinkEnv1Spec these
 *        hooks are imported-ghost-only per Stage 2c audit. The path gate
 *        IS the policy; the hook trait is the Spec opt-in.
 *   (c) hook method exists (SFINAE; absent ⇒ silent no-op).
 * Future Specs that need a different gating may add their own dispatch
 * helpers; this set covers SinkEnv1Spec E1. */

template <typename Spec>
static void nlr_dispatch_ghost_write_detector_begin(const neighbor_loop_args& args,
                                                    const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_write_detector_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwd_begin<Spec>::value) {
                Spec::ghost_write_detector_begin(args, plan);
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_write_detector_end(const neighbor_loop_args& args,
                                                  const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_write_detector_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwd_end<Spec>::value) {
                Spec::ghost_write_detector_end(args, plan);
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_writeback_begin(const neighbor_loop_args& args,
                                               const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_writeback_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwb_begin<Spec>::value) {
                Spec::ghost_writeback_begin(args, plan);
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_writeback_end(const neighbor_loop_args& args,
                                             const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_writeback_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwb_end<Spec>::value) {
                Spec::ghost_writeback_end(args, plan);
            }
        }
    }
}
/* ============================================================================
 * Public entry: run_neighbor_loop<Spec>
 * ========================================================================== */

template <typename Spec>
void run_neighbor_loop(const neighbor_loop_args& args)
{
    /* ---- Compile-time spec consistency ---- */

    /* WritePattern ↔ AccumData/ScatterData consistency. */
    static_assert(Spec::write_pattern != WritePattern::ActiveReduceOnly ||
                  std::is_same<typename Spec::ScatterData, NoScatter>::value,
        "ActiveReduceOnly requires ScatterData == NoScatter");
    static_assert(Spec::write_pattern != WritePattern::NeighborScatter ||
                  std::is_same<typename Spec::AccumData, NoAccum>::value,
        "NeighborScatter requires AccumData == NoAccum");

    /* POD/device-copy contract — captured by value into Kokkos lambdas
     * and/or staged into UVM arrays. Trivially-copyable is the binding
     * requirement. */
    static_assert(std::is_trivially_copyable<typename Spec::CallScalars>::value,
        "Spec::CallScalars must be trivially-copyable (lambda capture by value)");
    static_assert(std::is_trivially_copyable<typename Spec::ActiveData>::value,
        "Spec::ActiveData must be trivially-copyable (UVM-staged)");
    static_assert(std::is_trivially_copyable<typename Spec::NeighborData>::value,
        "Spec::NeighborData must be trivially-copyable (built per-pair on device)");
    static_assert(std::is_trivially_copyable<typename Spec::AccumData>::value,
        "Spec::AccumData must be trivially-copyable (UVM-staged)");

    /* ---- Dispatch ----
     *
     * Selection precedence (highest first):
     *   1. GIZMO_NLR_FORCE_MODE=A → Mode A unconditionally.
     *   2. GIZMO_NLR_FORCE_MODE=B → Mode B (local if NTask==1, remote else).
     *   3. Threshold dispatch: if (sum_active>0 && sum_active<=TS &&
     *      max_active<=TM) → Mode B; else Mode A. Hierarchy of TS/TM:
     *      per-loop env > global env > Spec::modeb_threshold_{sum,max}
     *      constexpr defaults (64/64 today).
     *
     * Force-mode conflict cases (both old and new vars set, both old vars
     * set together, invalid value) are caught and endrun'd centrally inside
     * gizmo_nlr_force_mode(); see the env-config block earlier in this TU.
     *
     * Codex nuance (active-epoch caveat): Mode B host-frozen actives[] are
     * NOT bit-equivalent to Mode A's device-staged post-neighbor-list-build
     * actives. Oracle in this dispatch is Mode B tree vs Mode B brute on the
     * same frozen query — it does NOT cross-validate Mode A. Mode A vs
     * Mode B active-epoch consistency is a separate concern, deferred.
     */
    const NlrForceMode force_mode = gizmo_nlr_force_mode_for(Spec::loop_name);
    const bool force_a   = (force_mode == NlrForceMode::A);
    const bool force_b   = (force_mode == NlrForceMode::B);
    const bool oracle_on = gizmo_nlr_oracle_enabled_global() ||
                            gizmo_nlr_oracle_enabled_for(Spec::loop_name);

    /* PHASE0 timing scaffolding (3c.4b). Cached env-gate; mid-run env
     * changes do not take effect. When on, ALL MPI_Wtime calls inside the
     * runner are gated; off-path overhead is one branch per StageTimer
     * scope, no MPI_Wtime call. PHASE0_NLR measures only RUNNER-OWNED time:
     * caller-side ghost prep / detector / SinkTempInfo scatter are NOT
     * included. */
    const bool phase0_on = gizmo_nlr_phase0_diag_enabled();
    RunnerStageTimer tim = {};
    RunnerStageTimer *tim_ptr = phase0_on ? &tim : nullptr;
    const double t_runner_start = phase0_on ? MPI_Wtime() : 0.0;

    /* Threshold dispatch. Allreduce sum + max of args.num_active.
     * Skipped when force-mode envs are set (cheap path).
     * PHASE0 num_active_global is captured here when the threshold path
     * already did the Allreduce; on force paths an extra Allreduce is done
     * ONLY when phase0_on (codex constraint #4). */
    bool select_mode_b = force_b;
    int phase0_sum_active = -1;       /* -1 = not yet computed */
    if(!force_a && !force_b) {
        int local_act = args.num_active;
        int sum_act = 0, max_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_act, &max_act, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        phase0_sum_active = sum_act;
        /* Spec::modeb_threshold_{sum,max} via SFINAE; default 64/64. */
        const int spec_default_sum = nlr_spec_threshold_sum<Spec>(64);
        const int spec_default_max = nlr_spec_threshold_max<Spec>(64);
        const int TS = gizmo_nlr_modeb_threshold_sum_for(Spec::loop_name, spec_default_sum);
        const int TM = gizmo_nlr_modeb_threshold_max_for(Spec::loop_name, spec_default_max);
        select_mode_b = (sum_act > 0) && (sum_act <= TS) && (max_act <= TM);
    } else if(phase0_on || force_b) {
        /* Force path: dispatch logic skipped the Allreduce. Do it here for
         * phase0 num_active_global and for the forced-Mode-B size guard. */
        int local_act = args.num_active;
        int sum_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        phase0_sum_active = sum_act;
    }
    if(force_b) {
        nlr_abort_if_forced_modeb_too_large(
            Spec::loop_name, args.num_active, phase0_sum_active);
    }

    /* Compute the execution plan from the dispatch decision. Path is the
     * single-source-of-truth; predicates derive from it. */
    NeighborLoopPlan plan;
    if(force_a) {
        plan.path = NeighborLoopPlan::Path::ModeA_GpuNgl;
    } else if(force_b || select_mode_b) {
        plan.path = (NTask > 1) ? NeighborLoopPlan::Path::ModeB_Remote
                                : NeighborLoopPlan::Path::ModeB_Local;
    } else {
        plan.path = NeighborLoopPlan::Path::ModeA_GpuNgl;
    }
    plan.num_active_global = phase0_sum_active;   /* may be -1 when phase0_on=false */

    /* Optional dispatch trace. Rank-0 only to avoid spam. */
    if(gizmo_nlr_dispatch_trace_enabled()) {
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) {
            const char *src =
                force_a       ? "force" :
                force_b       ? "force" :
                                "threshold";
            fprintf(stderr, "[NLR DISPATCH caller=%s path=%s (%s) NTask=%d local_active=%d oracle=%d]\n",
                    Spec::loop_name, nlr_path_label(plan.path), src,
                    NTask, args.num_active, (int)oracle_on);
            fflush(stderr);
        }
    }

    /* ---- Stage radii once ---- */
    /* Computed via Spec::search_radius. Used for any path-conditional
     * prep (Mode A) AND the chosen walker. Pointer is call-lifetime only;
     * see contract in mesh/neighbor_loop_runner.h. */
    std::vector<double> radii(args.num_active);
    for(int aa = 0; aa < args.num_active; ++aa) {
        radii[aa] = Spec::search_radius(args, aa, args.active_list[aa]);
    }

    /* ---- Hard-corridor counter snapshot (always-on, every build) ---- */
    /* Mode B paths must NOT enter move_particles, ghost_exchange_impl, or
     * gpu_particles_arena_acquire, and must NOT mutate NumPart. Counters
     * live in declarations/lifecycle_counters.h, incremented at API entry
     * by the owner TUs. */
    const uint64_t s_drift0 = g_global_drift_counter;
    const uint64_t s_ghost0 = g_ghost_import_counter;
    const uint64_t s_arena0 = g_gpu_arena_acquire_counter;
    const int      s_np0    = NumPart;

    /* ---- Working copy of args; refreshed after path-conditional prep ---- */
    neighbor_loop_args effective_args = args;

    /* ---- Path-conditional prep on imported-ghost paths ---- */
    /* Mode A's substrate today satisfies its freshness + neighbor-pool
     * requirements via gizmo_request_filtered_ghost_import_fresh (full
     * global drift + ghost import). Mode B paths skip this entirely;
     * peer-local pool + lazy candidate drift is sufficient. The dt_prep_import
     * timer is 0 for Mode B paths (genuine 0 — the API isn't called). */
    if(nlr_path_uses_imported_ghosts(plan.path)) {
        StageTimer t_prep(tim_ptr ? &tim_ptr->dt_prep_import : nullptr);
        gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                   Spec::search_mode,
                                                   Spec::neighbor_type_mask,
                                                   args.active_list,
                                                   args.num_active,
                                                   radii.data(),
                                                   args.ghost_safety_factor);
        /* Ghost import grew NumPart and may have realloc'd P/CellP. Refresh
         * the runner's data view; only paths that imported ghosts read this
         * extended view (Mode B paths use the original args via copy). */
        effective_args.num_total = NumPart;
        effective_args.P         = P;
        effective_args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
    }

    /* ---- Spec lifecycle hooks (begin) ---- */
    /* Mode A imported-ghost ordering invariant (two channels: detector +
     * writeback). The Spec's writeback_begin/_end body is the per-flag
     * #ifdef union of all enabled physics-flag writebacks for that loop:
     *
     *   request_filtered_ghost_import_fresh  (above)
     *   ghost_write_detector_begin           (this hook)
     *   ghost_writeback_begin                (this hook; per-flag union)
     *   <run_mode_a kernel>
     *   ghost_writeback_end                  (reverse order, below)
     *   ghost_write_detector_end
     *   ghost_exchange_cleanup               (below)
     *
     * Each dispatch helper checks both the Spec's `uses_*` trait AND the
     * path-imports-ghosts predicate; on Mode B paths the predicate is
     * false and all four hook calls compile to no-ops. */
    nlr_dispatch_ghost_write_detector_begin<Spec>(effective_args, plan);
    nlr_dispatch_ghost_writeback_begin<Spec>(effective_args, plan);

    /* ---- Path dispatch ---- */
    switch(plan.path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl:
            /* Oracle on Mode A is a no-op (oracle compares Mode B vs Brute;
             * no Brute-vs-Mode-A oracle in scope). */
            run_mode_a<Spec>(effective_args, radii.data(), tim_ptr);
            break;
        case NeighborLoopPlan::Path::ModeB_Local:
            if(oracle_on) run_mode_b_local_with_oracle<Spec>(effective_args, radii.data(), tim_ptr);
            else          run_mode_b_local<Spec>(effective_args, radii.data(), tim_ptr);
            break;
        case NeighborLoopPlan::Path::ModeB_Remote:
            if(oracle_on) run_mode_b_remote_with_oracle<Spec>(effective_args, radii.data(), tim_ptr);
            else          run_mode_b_remote<Spec>(effective_args, radii.data(), tim_ptr);
            break;
    }

    /* ---- Spec lifecycle hooks (end, reverse order) ---- */
    nlr_dispatch_ghost_writeback_end<Spec>(effective_args, plan);
    nlr_dispatch_ghost_write_detector_end<Spec>(effective_args, plan);

    /* ---- Imported-ghost cleanup ---- */
    if(nlr_path_uses_imported_ghosts(plan.path) && NTask > 1) {
        ghost_exchange_cleanup();
    }

    /* ---- Hard-corridor enforcement (Mode B paths) ---- */
    /* HARD ABORT on any counter advance or NumPart change across the path
     * body. Always-on; cheap (4 uint64 compares + one int compare). */
    if(plan.path == NeighborLoopPlan::Path::ModeB_Local ||
       plan.path == NeighborLoopPlan::Path::ModeB_Remote) {
        const bool drift_violation = (g_global_drift_counter      != s_drift0);
        const bool ghost_violation = (g_ghost_import_counter      != s_ghost0);
        const bool arena_violation = (g_gpu_arena_acquire_counter != s_arena0);
        const bool np_violation    = (NumPart != s_np0);
        if(drift_violation || ghost_violation || arena_violation || np_violation) {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                    "[NLR CORRIDOR ABORT rank=%d caller=%s path=%s] Mode B path "
                    "violated tiny-N corridor invariant during run_neighbor_loop. "
                    "Counter deltas: drift=%llu ghost=%llu arena=%llu NumPart_pre=%d NumPart_post=%d\n",
                    rank, Spec::loop_name, nlr_path_label(plan.path),
                    (unsigned long long)(g_global_drift_counter - s_drift0),
                    (unsigned long long)(g_ghost_import_counter - s_ghost0),
                    (unsigned long long)(g_gpu_arena_acquire_counter - s_arena0),
                    s_np0, NumPart);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 81036);
        }
    }

    /* ---- PHASE0_NLR emit ---- */
    /* Stable prefix `PHASE0_NLR`. `caller=` is a field, not part of the
     * token, so future Specs keep the parser regex stable. The new
     * dt_prep_import field measures the runner-internal prep wall (Mode A
     * only; 0 on Mode B paths). Other fields per the path-specific
     * documentation in neighbor_loop_runner.h. */
    if(phase0_on) {
        tim.dt_total = MPI_Wtime() - t_runner_start;
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        static long long s_call_id = 0;
        ++s_call_id;
        std::printf("PHASE0_NLR rank=%d caller=%s path=%s call_id=%lld "
                    "num_active_local=%d num_active_global=%d "
                    "dt_prep_import=%.6g "
                    "dt_collect=%.6g dt_drift=%.6g dt_walk_self=%.6g "
                    "dt_walk_peer=%.6g dt_exchange_q=%.6g dt_exchange_r=%.6g "
                    "dt_reduce=%.6g dt_writeback=%.6g dt_total=%.6g\n",
                    rank, Spec::loop_name, nlr_path_label(plan.path), s_call_id,
                    args.num_active, phase0_sum_active,
                    tim.dt_prep_import,
                    tim.dt_collect, tim.dt_drift, tim.dt_walk_self,
                    tim.dt_walk_peer, tim.dt_exchange_q, tim.dt_exchange_r,
                    tim.dt_reduce, tim.dt_writeback, tim.dt_total);
        std::fflush(stdout);
    }
}

/* ============================================================================
 * run_neighbor_loop_iterative<Spec> — Phase 4.B.0 iterative entry point.
 *
 * STEP 2b (2026-05-10): framework + Mode B LOCAL iter dispatch only.
 * Mode A iter, Mode B remote iter, iterative oracle = HARD STUBS pending
 * step 2b.2 (Mode B remote helper extraction) and step 2c (Mode A
 * cached-CSR + iterative oracle). Per Phil + codex 2026-05-10: no slow /
 * regressed Mode A iterative path; harness in step 3 forces Mode B local
 * (np=1) for framework validation; np=2 + Mode A force = expected hard
 * abort until 2b.2 / 2c.
 *
 * Compile-time spec consistency:
 *   - IterControl == Iterative
 *   - IterScratch trivially copyable
 *   - max_iters >= 1; mode_a_csr_buffer_factor > 1.0
 *   - Spec::after_iter declared and callable (clean diagnostic naming the
 *     missing member if forgot)
 *
 * Runtime checks (fire BEFORE any arena/ghost/session touch):
 *   - num_subgroups >= 1 (caller short-circuits if globally empty)
 *   - num_subgroups > 1 requires Spec::SupportsSubgroups::value
 * ========================================================================== */

/* Helper: compile-time member detector for SupportsSubgroups (default false). */
template <typename Spec, typename = void>
struct nlr_supports_subgroups : std::false_type {};
template <typename Spec>
struct nlr_supports_subgroups<Spec, std::void_t<typename Spec::SupportsSubgroups>>
    : Spec::SupportsSubgroups {};

/* ============================================================================
 * NlrIterDriver<Spec> — constructor + destructor (step 2b).
 * ========================================================================== */
template <typename Spec>
NlrIterDriver<Spec>::NlrIterDriver(const neighbor_loop_args_iterative& a,
                                   const typename Spec::CallScalars& s)
    : args(a), cs(s), iter_index(0),
      local_active_total(0), global_active_total(-1),
      local_active_per_sg (a.num_subgroups, 0),
      global_active_per_sg(a.num_subgroups, 0),
      scratch_uvm(a.num_subgroups, nullptr),
      accum_uvm  (a.num_subgroups, nullptr),
      radii_uvm  (a.num_subgroups, nullptr),
      active_set_uvm  (a.num_subgroups, nullptr),
      active_set_size (a.num_subgroups, 0),
      /* Mode A iterative cached-CSR state (step 2c.2): zero-init per subgroup;
       * UVM allocations land lazily on first Mode A iter dispatch. */
      mode_a_cached_gnl       (a.num_subgroups, gpu_neighbor_list_t{}),
      mode_a_csr_offset_lookup(a.num_subgroups, nullptr),
      mode_a_csr_buffered_h   (a.num_subgroups, nullptr),
      mode_a_csr_valid        (a.num_subgroups, false),
      /* effective_args: copy of base slice; refreshed on Mode A ghost import. */
      effective_args(static_cast<const neighbor_loop_args&>(a))
{
    using IterScratch = typename Spec::IterScratch;
    using AccumData   = typename Spec::AccumData;

    /* Build a base neighbor_loop_args view per subgroup so we can call
     * Spec::search_radius for each (subgroup, slot, particle index) tuple.
     * search_radius is host-only and predates any drift / arena work. */
    for (int sg = 0; sg < args.num_subgroups; sg++) {
        const NlrSubgroup& sgr = args.subgroups[sg];
        const int n = sgr.num_active_local;
        active_set_size[sg] = n;
        if (n <= 0) continue;     /* leave pointers as nullptr; pair_kernel is no-op */

        scratch_uvm[sg]    = (IterScratch *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(IterScratch));
        accum_uvm  [sg]    = (AccumData   *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(AccumData));
        radii_uvm  [sg]    = (double      *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(double));
        active_set_uvm[sg] = (int         *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(int));

        /* Byte-zero IterScratch ONCE — persists across iters per v4 §4. */
        std::memset(scratch_uvm[sg], 0, n * sizeof(IterScratch));
        /* AccumData NOT zeroed here; runner zeros via Spec::zero_accum at start of each iter. */

        /* Build a sub-args view for Spec::search_radius. Mirrors the per-
         * subgroup args the iterative dispatch will hand to lower-level
         * helpers (so search_radius signature stays a base neighbor_loop_args&). */
        neighbor_loop_args sub = args;          /* base slice; aux/CellP/etc. carry through */
        sub.active_list = sgr.active_indices;
        sub.num_active  = n;

        for (int slot = 0; slot < n; slot++) {
            radii_uvm[sg][slot]      = Spec::search_radius(sub, slot, sgr.active_indices[slot]);
            active_set_uvm[sg][slot] = slot;     /* {0..n-1} initial */
        }
        local_active_total += n;
    }

    /* ========================================================================
     * Iterative-oracle field init (step 2c.4).
     *
     * `oracle_enabled` is set from the same env-gate the non-iter oracle
     * uses; Mode A iterative oracle is hard-stubbed in the outer body
     * (run_neighbor_loop_iterative step 5.b) BEFORE this constructor runs,
     * so reaching here with oracle_enabled=true implies Mode B production
     * path. Vectors stay empty when oracle_enabled=false (UVM frees in
     * dtor are guarded on pointer state, so empty vectors are safe).
     * ====================================================================== */
    oracle_enabled = gizmo_nlr_oracle_enabled_global() ||
                     gizmo_nlr_oracle_enabled_for(Spec::loop_name);
    if (oracle_enabled) {
        scratch_oracle_uvm    .assign(args.num_subgroups, nullptr);
        accum_oracle_uvm      .assign(args.num_subgroups, nullptr);
        radii_oracle_uvm      .assign(args.num_subgroups, nullptr);
        active_set_oracle_uvm .assign(args.num_subgroups, nullptr);
        active_set_oracle_size.assign(args.num_subgroups, 0);
        status_prod_uvm       .assign(args.num_subgroups, nullptr);
        status_oracle_uvm     .assign(args.num_subgroups, nullptr);
        new_h_prod_uvm        .assign(args.num_subgroups, nullptr);
        new_h_oracle_uvm      .assign(args.num_subgroups, nullptr);

        for (int sg = 0; sg < args.num_subgroups; sg++) {
            const NlrSubgroup& sgr = args.subgroups[sg];
            const int n = sgr.num_active_local;
            active_set_oracle_size[sg] = n;
            if (n <= 0) continue;

            scratch_oracle_uvm   [sg] = (IterScratch *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(IterScratch));
            accum_oracle_uvm     [sg] = (AccumData   *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(AccumData));
            radii_oracle_uvm     [sg] = (double      *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(double));
            active_set_oracle_uvm[sg] = (int         *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(int));
            status_prod_uvm      [sg] = (int         *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(int));
            status_oracle_uvm    [sg] = (int         *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(int));
            new_h_prod_uvm       [sg] = (double      *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(double));
            new_h_oracle_uvm     [sg] = (double      *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(n * sizeof(double));

            std::memset(scratch_oracle_uvm[sg], 0, n * sizeof(IterScratch));

            neighbor_loop_args sub = args;
            sub.active_list = sgr.active_indices;
            sub.num_active  = n;

            for (int slot = 0; slot < n; slot++) {
                radii_oracle_uvm[sg][slot]      = Spec::search_radius(sub, slot, sgr.active_indices[slot]);
                active_set_oracle_uvm[sg][slot] = slot;
                status_prod_uvm[sg][slot]       = static_cast<int>(IterStatus::NeedsMore);
                status_oracle_uvm[sg][slot]     = static_cast<int>(IterStatus::NeedsMore);
                new_h_prod_uvm[sg][slot]        = radii_uvm[sg][slot];
                new_h_oracle_uvm[sg][slot]      = radii_oracle_uvm[sg][slot];
            }
        }
    }
}

template <typename Spec>
NlrIterDriver<Spec>::~NlrIterDriver()
{
    /* DeviceContext cleanup (codex 2c.1 review): only fire if init actually
     * completed. Stubbed/aborted init paths leave ctx_initialized=false →
     * cleanup_device_context is NOT called, preventing free of unallocated
     * resources. Direct Spec call (NOT NlrDeviceContextCleanupGuard) — the
     * guard is RAII-only and wouldn't compose with conditional init.
     *
     * Ordering (codex 2c.1 hardening note): cleanup_device_context fires
     * BEFORE the per-subgroup UVM frees below. Spec::cleanup_device_context
     * MUST only free resources it OWNS (e.g., UVM arrays it allocated in
     * populate_device_context). It MUST NOT assume the driver's per-subgroup
     * UVM arrays (scratch/accum/radii/active_set) remain valid AFTER cleanup
     * returns — those frees happen below. In practice cleanup hooks like
     * sink_feed's only free their own UVM, so this ordering is fine; the
     * note exists to flag the contract for future Specs. */
    if (ctx_initialized) {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(args, ctx);
        }
    }

    /* Iterative-oracle ctx_oracle cleanup (step 2c.4). Independent
     * lifecycle from production ctx: own populate + own cleanup. Guarded
     * on ctx_oracle_initialized — only fire if init actually completed. */
    if (ctx_oracle_initialized) {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(args, ctx_oracle);
        }
    }

    /* Free Mode A cached CSR/lookup state by POINTER STATE (codex 2c.2-fix
     * blocker #5: mode_a_csr_valid=false can mean "allocated but invalid,
     * pending rebuild" if the rebuild trigger fired but rebuild itself
     * hadn't completed yet; check pointers, not the flag).
     *
     * gpu_ngb_list_free passes the SIDX pointer so the step-persistent SIDX
     * cache is preserved across iterative calls (matches sink_env1/feed/swk
     * idiom). */
    {
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
    }

    /* Arena cleanup ONCE per call (matches non-iter run_mode_a line 1718). */
    if (arena_acquired) {
        gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
        arena_acquired = false;
    }

    /* Imported-ghost cleanup ONCE per call (codex 2c.2-fix; matches non-iter
     * line 2094-2097). Only fires for Mode A multi-rank paths that imported.
     * Runner's Mode A imports an exact-query ghost pool sized to the iter's
     * actives+radii; the caller (e.g. density()) is responsible for any
     * downstream handoff pool it needs AFTER runner-return (post-finalize
     * fresh broad import — not "keep this exact-query pool alive"). */
    if (ghost_import_done && NTask > 1) {
        ghost_exchange_cleanup();
        ghost_import_done = false;
    }

    for (int sg = 0; sg < args.num_subgroups; sg++) {
        if (scratch_uvm[sg])    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_uvm[sg]);
        if (accum_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(accum_uvm[sg]);
        if (radii_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm[sg]);
        if (active_set_uvm[sg]) Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(active_set_uvm[sg]);
    }

    /* Iterative-oracle UVM frees (step 2c.4). Vectors are sized to
     * num_subgroups only when oracle_enabled was true at construction; if
     * empty, the per-sg loop iterates zero times. Pointers are nullptr by
     * default — kokkos_free guarded on non-null. */
    for (size_t sg = 0; sg < scratch_oracle_uvm.size(); sg++) {
        if (scratch_oracle_uvm[sg])    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_oracle_uvm[sg]);
        if (accum_oracle_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(accum_oracle_uvm[sg]);
        if (radii_oracle_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_oracle_uvm[sg]);
        if (active_set_oracle_uvm[sg]) Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(active_set_oracle_uvm[sg]);
        if (status_prod_uvm[sg])       Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(status_prod_uvm[sg]);
        if (status_oracle_uvm[sg])     Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(status_oracle_uvm[sg]);
        if (new_h_prod_uvm[sg])        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(new_h_prod_uvm[sg]);
        if (new_h_oracle_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(new_h_oracle_uvm[sg]);
    }
}

/* ============================================================================
 * Path-specific DeviceContext init (step 2c.1).
 *
 * Mode B: bind to caller's args.P/CellP — Mode B reads owner-local per-particle
 *         state directly (lazy-drift contract). Static_asserts mirror run_mode_b_*
 *         (DeviceContext base + trivially copyable).
 * Mode A: bind to arena-resident gpu_particles_arena_P()/CellP() — Mode A's
 *         GPU NGL build + kernel walk operate on the arena copy. Body lands in
 *         step 2c.2 alongside the arena_acquire call site.
 * ========================================================================== */
template <typename Spec>
void NlrIterDriver<Spec>::initialize_device_context_mode_b()
{
    using DeviceCtx = typename Spec::DeviceContext;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");

    /* Double-init guard (codex 2c.1 hardening): once Mode A init lands in
     * step 2c.2, accidental call sequences (e.g., path mis-route calling
     * both inits) would silently double-populate / double-cleanup. Catch
     * loudly here. */
    if (ctx_initialized) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_b] FATAL: "
                "ctx_initialized=true on entry. Double-init would orphan resources.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        endrun(81209);
    }

    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;

    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    ctx_initialized = true;
}

template <typename Spec>
void NlrIterDriver<Spec>::initialize_device_context_mode_a_after_arena()
{
    using DeviceCtx = typename Spec::DeviceContext;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");

    /* Double-init guard (codex 2c.1 hardening). */
    if (ctx_initialized) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_a_after_arena] FATAL: "
                "ctx_initialized=true on entry. Double-init would orphan resources.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        endrun(81210);
    }
    /* Arena must already have been acquired (codex 2c.2 strict lifecycle). */
    if (!arena_acquired) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_a_after_arena] FATAL: "
                "arena not acquired. Call acquire_arena_and_init_ctx_mode_a() instead.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        endrun(81211);
    }

    /* Bind to arena-resident P_gpu / CellP_gpu. Use the driver's
     * effective_args (codex 2c.2-fix): if ghost import ran above, these
     * point at the refreshed POST-IMPORT global P/CellP/NumPart; otherwise
     * they equal the original base args. CellP=NULL is legitimate (gas-free). */
    ctx.P         = gpu_particles_arena_P();
    ctx.CellP     = (effective_args.CellP != nullptr) ? gpu_particles_arena_CellP() : nullptr;
    ctx.num_total = effective_args.num_total;

    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(effective_args, ctx);
    }
    ctx_initialized = true;
}

template <typename Spec>
void NlrIterDriver<Spec>::acquire_arena_and_init_ctx_mode_a()
{
    if (arena_acquired) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::acquire_arena_and_init_ctx_mode_a] FATAL: "
                "arena already acquired. Single-acquire-per-call contract violated.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        endrun(81212);
    }

    /* === (1) Imported-ghost prep ONCE per call (codex 2c.2-fix + 2c.3 step 9) ===
     * Mode A requires neighbors that live on peer ranks to be imported as
     * ghosts. Mirrors non-iter run_neighbor_loop lines 2037-2052.
     *
     * Multi-subgroup union semantics (step 2c.3 step 9): mask = OR of all
     * subgroups' j_type_bitmask; actives = concatenation of all subgroups'
     * full active_indices (initial iter-0 set; later iters' compactions
     * narrow but don't grow); radii = matching concatenation oversized via
     * mode_a_csr_buffer_factor. Single-subgroup case: trivial length-1
     * union = original behavior. Multi-subgroup over-imports in cross-bm
     * cases (deliberate physics/perf tradeoff per codex plan v2 ch. 6).
     *
     * On rebuild via rebuild_mode_a_arena_and_ctx_for_current_active_union,
     * the union covers only CURRENT compacted actives (post-Converged-removal)
     * so the rebuild import doesn't over-cover converged slots. */
    if (NTask > 1) {
        std::vector<int>    union_actives;
        std::vector<double> union_radii_oversized;
        unsigned int        mask_union = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            mask_union |= args.subgroups[sg].j_type_bitmask;
            const NlrSubgroup& sgr = args.subgroups[sg];
            for (int k = 0; k < sgr.num_active_local; k++) {
                union_actives.push_back(sgr.active_indices[k]);
                union_radii_oversized.push_back(radii_uvm[sg][k] * Spec::mode_a_csr_buffer_factor);
            }
        }
        const int union_n = (int)union_actives.size();
        gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                   Spec::search_mode,
                                                   mask_union,
                                                   (union_n > 0) ? union_actives.data() : nullptr,
                                                   union_n,
                                                   (union_n > 0) ? union_radii_oversized.data() : nullptr,
                                                   args.ghost_safety_factor);
        ghost_import_done = true;

        /* Refresh effective_args from globals (ghost import grew NumPart and
         * may have realloc'd P/CellP). Matches non-iter line 2049-2051. */
        effective_args.num_total = NumPart;
        effective_args.P         = P;
        effective_args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
    }

    /* === (2) Arena acquire ONCE per call, with refreshed args === */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(effective_args.num_total,
                                 effective_args.P, effective_args.CellP);
    arena_acquired = true;
    arena_acquire_count++;       /* step 3 harness counter */

    /* === (3) Bind ctx to arena-resident pointers + populate extended DeviceContext === */
    initialize_device_context_mode_a_after_arena();
}

/* ============================================================================
 * rebuild_mode_a_arena_and_ctx_for_current_active_union — codex 2c.3 plan v2
 * tightening 5: self-sufficient rebuild method (no parameters).
 *
 * Builds the union from driver's own subgroup state. Used by the Mode A
 * outer-iter pre-dispatch invalidation sweep (step 8). Invalidates ALL
 * subgroup CSR caches as part of the lifecycle (pitfall 2: cross-subgroup
 * neighbor indices become stale after arena re-acquire).
 *
 * Union semantics (pitfall 4): mask = OR of all globally-active subgroups'
 * j_type_bitmask; active list = concatenation of current compacted actives
 * across globally-active subgroups; radii = matching concatenation of
 * radii_uvm[sg][slot] * mode_a_csr_buffer_factor. Over-imports in
 * multi-bm cases (each active gets ghosts of all union types instead of
 * just its own mask). Cosmological FIRE typical = single bm-group = no
 * over-import. Documented as deliberate physics/perf tradeoff (codex
 * approved 2c.3 plan v2 change 6).
 *
 * Pre-condition: drv.global_active_per_sg has been computed (set by the
 * outer iter loop's per-iter Allreduce). On iter 0 it may all be 0 (not
 * yet Allreduced); caller should NOT invoke this method at iter 0 — the
 * iter-0 ghost import + arena acquire happens in
 * acquire_arena_and_init_ctx_mode_a per the original path.
 * ========================================================================== */
template <typename Spec>
void NlrIterDriver<Spec>::rebuild_mode_a_arena_and_ctx_for_current_active_union()
{
    if (NTask <= 1) {
        /* Single rank: no ghost pool to manage. Just invalidate all CSR
         * caches; per-subgroup dispatch will rebuild local CSRs on first
         * access against the unchanged arena/pool. */
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
        return;
    }

    /* Build union from current active_set across globally-active subgroups. */
    std::vector<int>    union_actives;
    std::vector<double> union_radii_oversized;
    unsigned int        mask_union = 0;
    for (int sg = 0; sg < args.num_subgroups; sg++) {
        if (global_active_per_sg[sg] <= 0) continue;     /* skip globally-converged */
        mask_union |= args.subgroups[sg].j_type_bitmask;
        const NlrSubgroup& sgr = args.subgroups[sg];
        for (int k = 0; k < active_set_size[sg]; k++) {
            int slot = active_set_uvm[sg][k];
            union_actives.push_back(sgr.active_indices[slot]);
            union_radii_oversized.push_back(radii_uvm[sg][slot] * Spec::mode_a_csr_buffer_factor);
        }
    }

    /* === (0) Invalidate ALL subgroup CSR caches (pitfall 2 — cross-subgroup
     * staleness after arena teardown). Free by pointer state. */
    {
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
    }

    /* === (1) cleanup_device_context for extended Specs === */
    if (ctx_initialized) {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(effective_args, ctx);
        }
        ctx_initialized = false;
    }

    /* === (2) END current arena view BEFORE mutating the global ghost pool === */
    if (arena_acquired) {
        gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
        arena_acquired = false;
    }

    /* === (3) Destroy current ghost pool === */
    if (ghost_import_done) {
        ghost_exchange_cleanup();
        ghost_import_done = false;
    }

    /* === (4) Reimport with UNION radii / actives / mask (deliberate over-import
     * documented as 2c.3 plan v2 change 6 acceptable tradeoff) === */
    {
        const int union_n = (int)union_actives.size();
        gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                   Spec::search_mode,
                                                   mask_union,
                                                   (union_n > 0) ? union_actives.data() : nullptr,
                                                   union_n,
                                                   (union_n > 0) ? union_radii_oversized.data() : nullptr,
                                                   args.ghost_safety_factor);
        ghost_import_done = true;
    }

    /* === (5) Refresh effective_args from post-import globals === */
    effective_args.num_total = NumPart;
    effective_args.P         = P;
    effective_args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;

    /* === (6) Fresh arena view of new pool === */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(effective_args.num_total,
                                 effective_args.P, effective_args.CellP);
    arena_acquired = true;
    arena_acquire_count++;       /* step 3 harness counter */
    csr_rebuild_count++;         /* step 3 harness counter — union rebuild fired */

    /* === (7) Rebind ctx + (8) populate extended ctx for NEW pool === */
    initialize_device_context_mode_a_after_arena();

    /* Codex 2c.3 blocker #4 fix: do NOT re-fire reset_per_iter_device_context
     * here. The outer iter loop's step (a-pre) calls this method BEFORE
     * step (a0), which fires the reset hook exactly once per outer iter
     * over the freshly-populated ctx. Re-firing here would cause a
     * double-reset on rebuild iters and violate the "fires once per outer
     * iter" contract used by the synthetic harness validation. */
}

/* ============================================================================
 * Per-iter, per-subgroup dispatch helpers (step 2b — Mode B local only).
 *
 * Each helper takes the driver, a subgroup index, and runs ONE iteration's
 * pair_kernel walk for that subgroup writing into driver-owned accum_uvm.
 * apply_active_writeback is NOT called here — that's final-only after the
 * outer iter loop completes.
 *
 * The Mode B local helper composes existing lower-level helpers
 * (build_self_actives_host_pre_drift / collect_candidates_pre_drift /
 * lazy_drift_candidates / evaluate_pairs_post_drift). SSOT preserved with
 * existing run_mode_b_local — same helper chain, just driver-owned output
 * buffer instead of stack vector.
 *
 * Iterative oracle iter still hard-stubbed (step 2c.4); Mode A iter lands
 * via the sibling nlr_iter_dispatch_subgroup_mode_a helper (step 2c.2).
 *
 * DEVICE CONTEXT LIFETIME — RESOLVED IN 2c.1 + 2c.2:
 * The dispatch helper now uses driver-owned drv.ctx (populated once via
 * NlrIterDriver::initialize_device_context_mode_b at iter-0 entry, cleaned
 * once in destructor via cleanup_device_context if extended). Per-iter
 * reset is handled at the outer iter loop via Spec::reset_per_iter_device_context.
 * No per-iter / per-subgroup ctx rebuild here; the older 2b "rebuilds DeviceCtx
 * per-iter + per-subgroup" hack is gone.
 * ========================================================================== */
/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_b_remote<Spec> (step 2b.2).
 *
 * Per-iter Mode B REMOTE dispatch for one subgroup. Composes the same
 * extracted helper (mode_b_remote_evaluate_into_buffer<Spec, false>) the
 * non-iterative wrapper uses — SSOT preserved. The driver-owned compacted
 * AccumData buffer is the only difference. apply_active_writeback NOT
 * called here — that's final-only after the outer iter loop.
 *
 * Step 2b.2 scope: non-oracle only. ORACLE template parameter is hard-
 * coded false at the helper invocation; iterative oracle (brute-dry-run-
 * FIRST production-SECOND per iter) lands in step 2c.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_remote(NlrIterDriver<Spec>& drv, int sg)
{
    using AccumData = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];
    /* Note: even if n_compacted == 0 on this rank, the helper MUST be
     * entered because Mode B remote uses collectives (alltoallv exchange
     * of queries / replies). Empty rank participates with N=0. */

    /* Compact active particle indices + radii for the still-active slots. */
    std::vector<int>    active_particle_indices(n_compacted);
    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
        radii_compacted[k]         = drv.radii_uvm[sg][slot];
    }

    neighbor_loop_args sub = drv.args;
    sub.active_list = (n_compacted > 0) ? active_particle_indices.data() : nullptr;
    sub.num_active  = n_compacted;

    /* Driver-owned compacted AccumData buffer. Helper writes into this;
     * we scatter back into driver.accum_uvm[sg][slot] for active slots
     * after the helper returns. Explicit Spec::zero_accum per slot
     * (codex 2b.2 review fix): defensive against future evaluate_pairs
     * variants + makes the per-iter AccumData zero contract explicit at
     * the caller level. */
    std::vector<AccumData> accums_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);
    }

    mode_b_remote_evaluate_into_buffer<Spec, RemoteHelperMode::Production>(
        sub,
        radii_compacted.data(),
        drv.cs,                          /* driver-owned CallScalars */
        drv.ctx,                         /* driver-owned DeviceContext (step 2c.1) */
        (unsigned int)sgr.j_type_bitmask, /* per-subgroup mask (step 2c.3 step 4) */
        (n_compacted > 0) ? accums_compacted.data() : nullptr,
        /*actives_out=*/nullptr,         /* iter path: skip dumps (proper
                                             diagnostic plumbing lands later) */
        /*tim=*/nullptr);

    /* Scatter compacted accums back into driver-owned per-slot accum_uvm.
     * Slots NOT in active_set_uvm keep their last-evaluated value (the
     * v4.4 invariant: converged slots' final accum persists for
     * apply_active_writeback). */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        drv.accum_uvm[sg][slot] = accums_compacted[k];
    }
}

/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_a<Spec>(drv, sg) — Phase 4.B.0 step 2c.2.
 *
 * Per-iter Mode A dispatch for one subgroup. Cached oversized CSR + rebuild
 * trigger preserve legacy density_gpu.cc:140-220 session semantics.
 *
 * Lifecycle assumptions (set up by run_neighbor_loop_iterative + driver):
 *   - drv.acquire_arena_and_init_ctx_mode_a() ran ONCE before any subgroup
 *     dispatch. drv.ctx.P/CellP point at arena-resident P_gpu/CellP_gpu.
 *   - drv.arena_acquired = true; drv.ctx_initialized = true.
 *
 * Per-iter sequence:
 *   1. If !drv.mode_a_csr_valid[sg]: rebuild CSR with current compacted-
 *      active radii × Spec::mode_a_csr_buffer_factor. Frees old gnl/lookup
 *      (passing SIDX). Populates csr_offset_lookup + csr_buffered_h keyed
 *      on subgroup-slot-at-build-time.
 *   2. ghost_write_detector_begin + ghost_writeback_begin (gated on Spec
 *      traits + path → fires for Mode A iter + uses_ghost_writeback=true,
 *      e.g. ags_density wakeup; no-op for non-j-write Specs / harness).
 *   3. Stage d_actives via Spec::load_active (compacted indices).
 *   4. Pair-kernel launch over compacted active set, writing into a
 *      compacted UVM accums buffer.
 *   5. ghost_writeback_end + ghost_write_detector_end.
 *   6. Scatter compacted accums into driver-owned accum_uvm[sg][slot].
 * NO apply_active_writeback — final-only at driver level (post-iter-loop).
 *
 * CSR row-key invariant (codex v4.1 §2.A.1 restatement): rows in
 * cached_gnl are in COMPACTED-AT-BUILD-TIME order. csr_offset_lookup[sg]
 * maps subgroup-slot → compacted-build-time row index. Active-set
 * compaction in later iters does NOT re-key rows; converged slots simply
 * aren't walked. Rebuild creates a fresh lookup and discards the old one.
 *
 * Buffer-exceedance rebuild trigger: handled at outer-iter level AFTER
 * after_iter (post-radii mutation on AdjustRadius). Helper just consumes
 * the csr_valid[sg] flag.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_a(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData   = typename Spec::ActiveData;
    using AccumData    = typename Spec::AccumData;
    using ScatterData  = typename Spec::ScatterData;
    using NeighborData = typename Spec::NeighborData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];

    /* Collective-symmetry: ghost_writeback_begin/end MUST fire on every rank
     * for uses_ghost_writeback Specs (codex 2c.2-fix blocker #1) — even
     * empty-actives ranks, otherwise reverse-comm deadlocks. Hoist the
     * dispatch hooks outside the n_compacted<=0 short-circuit. The
     * kernel-side work between begin/end is gated on n_compacted > 0. */
    std::vector<int> active_particle_indices_iter(n_compacted > 0 ? n_compacted : 0);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices_iter[k] = sgr.active_indices[slot];
    }
    NeighborLoopPlan plan;
    plan.path              = NeighborLoopPlan::Path::ModeA_GpuNgl;
    plan.num_active_global = drv.global_active_total;
    neighbor_loop_args sub = drv.effective_args;
    sub.active_list = (n_compacted > 0) ? active_particle_indices_iter.data() : nullptr;
    sub.num_active  = n_compacted;

    /* ===== (1) Build / rebuild CSR if invalid =====
     * Codex 2c.3 blocker #1 fix: the reimport/re-arena lifecycle lives
     * EXCLUSIVELY in the outer pre-dispatch union rebuild (run_neighbor_loop_iterative
     * step a-pre), which uses the union mask across all globally-active
     * subgroups. Inside this helper, when !mode_a_csr_valid[sg], the ghost
     * pool / arena / ctx are already correct (handled by the outer sweep);
     * we just need to free + rebuild the local CSR for this subgroup. */
    gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                         Spec::loop_name);
    if (!drv.mode_a_csr_valid[sg] && n_compacted > 0) {
        /* Free old cached CSR / lookup state by POINTER STATE (codex
         * 2c.2-fix blocker #5: mode_a_csr_valid=false can mean "allocated
         * but invalid, pending rebuild"; check pointers, not flag). */
        if (drv.mode_a_cached_gnl[sg].offsets != nullptr ||
            drv.mode_a_cached_gnl[sg].neighbors != nullptr) {
            gpu_ngb_list_free(&drv.mode_a_cached_gnl[sg], sidx);
            drv.mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
        }
        if (drv.mode_a_csr_offset_lookup[sg]) {
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_offset_lookup[sg]);
            drv.mode_a_csr_offset_lookup[sg] = nullptr;
        }
        if (drv.mode_a_csr_buffered_h[sg]) {
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_buffered_h[sg]);
            drv.mode_a_csr_buffered_h[sg] = nullptr;
        }

        /* Build oversized CSR with CURRENT active set (compacted) and
         * oversized radii × buffer_factor. */
        std::vector<int>    active_particle_indices_host_build(n_compacted);
        std::vector<double> radii_buffered_host_build(n_compacted);
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            active_particle_indices_host_build[k] = sgr.active_indices[slot];
            radii_buffered_host_build[k]          = drv.radii_uvm[sg][slot] * Spec::mode_a_csr_buffer_factor;
        }
        double *radii_oversized_uvm = (double *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
            n_compacted * sizeof(double));
        for (int k = 0; k < n_compacted; k++) radii_oversized_uvm[k] = radii_buffered_host_build[k];

        gpu_ngb_list_build(drv.ctx.P, drv.ctx.num_total,
                           active_particle_indices_host_build.data(), n_compacted,
                           Spec::search_mode,
                           (int)sgr.j_type_bitmask,
                           &drv.mode_a_cached_gnl[sg], sidx,
                           1.0, radii_oversized_uvm, NULL, Spec::loop_name,
                           nlr_spec_symmetric_j_radius_scale<Spec>());
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_oversized_uvm);

        /* Allocate csr_offset_lookup + csr_buffered_h sized to subgroup max.
         * Lookup keyed on subgroup-slot (codex v4.1 §2.A.1 invariant):
         *   csr_offset_lookup[slot] = build-time row index (0..n_compacted-1),
         *                              or -1 if slot wasn't in build-time set
         *                              (shouldn't happen — rebuild covers all
         *                              current actives, and converged slots
         *                              never re-enter active_set). */
        const int n_max = sgr.num_active_local;
        drv.mode_a_csr_offset_lookup[sg] = (int *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
            n_max * sizeof(int));
        drv.mode_a_csr_buffered_h[sg]    = (double *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
            n_max * sizeof(double));
        for (int s = 0; s < n_max; s++) {
            drv.mode_a_csr_offset_lookup[sg][s] = -1;
            drv.mode_a_csr_buffered_h[sg][s]    = 0.0;
        }
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            drv.mode_a_csr_offset_lookup[sg][slot] = k;
            drv.mode_a_csr_buffered_h[sg][slot]    = radii_buffered_host_build[k];
        }
        drv.mode_a_csr_valid[sg] = true;
        drv.csr_local_rebuild_count++;       /* step 3 harness counter */
    }

    /* ===== (2) Ghost write detector + writeback begin =====
     * UNCONDITIONAL (codex 2c.2-fix blocker #1): fires on every rank for
     * uses_ghost_writeback Specs, INCLUDING empty-actives ranks, to keep
     * MPI reverse-comm collectives in lockstep with non-empty ranks.
     * Plan.path == ModeA_GpuNgl gates the dispatch hooks via
     * nlr_path_uses_imported_ghosts; uses effective_args (post-import). */
    nlr_dispatch_ghost_write_detector_begin<Spec>(sub, plan);
    nlr_dispatch_ghost_writeback_begin<Spec>(sub, plan);

    if (n_compacted > 0) {
        /* ===== (3) Stage d_actives + (4) pair-kernel launch =====
         * CSR ROW LOOKUP USAGE (codex 2c.2-fix blocker #2 verbatim):
         *   slot = active_set[k]
         *   row  = csr_offset_lookup[slot]    (build-time row index)
         *   i    = d_active[row]              (particle index at build-time)
         *   h    = radii_uvm[slot]            (current radius — kernel filter)
         * CSR build uses drv.ctx.num_total (= post-import effective num_total).
         * Walk gnl.offsets[row]..offsets[row+1]. */
        ActiveData *d_actives = (ActiveData *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
            n_compacted * sizeof(ActiveData));
        AccumData *d_accums = (AccumData *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
            n_compacted * sizeof(AccumData));
        {
            auto cs_ref = drv.cs;
            const typename Spec::DeviceContext dctx_local = drv.ctx;
            int    *active_set_arr = drv.active_set_uvm[sg];
            int    *csr_lookup     = drv.mode_a_csr_offset_lookup[sg];
            double *radii_arr      = drv.radii_uvm[sg];
            int     *d_active_arr  = drv.mode_a_cached_gnl[sg].d_active;
            int64_t *offsets       = drv.mode_a_cached_gnl[sg].offsets;
            int     *neighbors     = drv.mode_a_cached_gnl[sg].neighbors;

            gizmo_gpu_kernel_launch("nlr_iter_stage_active", n_compacted, KOKKOS_LAMBDA(int k) {
                int slot = active_set_arr[k];
                int row  = csr_lookup[slot];
                int i    = d_active_arr[row];
                double h = radii_arr[slot];
                d_actives[k] = Spec::load_active(dctx_local, slot, i, h, cs_ref);
            });

            gizmo_gpu_kernel_launch(Spec::loop_name, n_compacted, KOKKOS_LAMBDA(int k) {
                int slot = active_set_arr[k];
                int row  = csr_lookup[slot];
                Spec::zero_accum(d_accums[k]);
                const ActiveData& a = d_actives[k];
                ScatterData s{};
                int64_t start = offsets[row], end = offsets[row + 1];
                for (int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    IdentitySidecar id{};
                    NeighborData nb = Spec::load_neighbor(dctx_local, j, id, a);
                    Spec::pair_kernel(a, nb, d_accums[k], s);
                }
            });
        }

        /* ===== (6) Scatter compacted accums into driver accum_uvm ===== */
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            drv.accum_uvm[sg][slot] = d_accums[k];
        }

        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_actives);
    }

    /* ===== (5) Ghost writeback end + detector end (unconditional) ===== */
    nlr_dispatch_ghost_writeback_end<Spec>(sub, plan);
    nlr_dispatch_ghost_write_detector_end<Spec>(sub, plan);
}

template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_local(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData  = typename Spec::ActiveData;
    using AccumData   = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];
    if (n_compacted <= 0) return;     /* All actives in this subgroup converged. */

    /* Build a base sub-args restricted to the still-active slots in this subgroup.
     * This is per-iter — the active_list/num_active reflect the COMPACTED set
     * after Converged compaction from prior iters. */
    std::vector<int> active_particle_indices(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
    }
    neighbor_loop_args sub = drv.args;
    sub.active_list = active_particle_indices.data();
    sub.num_active  = n_compacted;

    /* Per-active radii in compacted order (helpers expect contiguous radii array). */
    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        radii_compacted[k] = drv.radii_uvm[sg][drv.active_set_uvm[sg][k]];
    }

    /* DeviceContext: driver-owned (step 2c.1). Populated once at iter-0 by
     * NlrIterDriver::initialize_device_context_mode_b(). Per-iter reset hook
     * (if Spec declares it) fires at outer iter level, before this dispatch. */

    /* Stage actives + zero accums host-side (compacted). */
    std::vector<ActiveData> actives_compacted(n_compacted);
    build_self_actives_host_pre_drift<Spec>(sub, drv.ctx, radii_compacted.data(),
                                              drv.cs, actives_compacted.data());

    std::vector<AccumData> accums_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);     /* per-iter zero (codex constraint 2) */
    }

    /* Collect → drift → evaluate (same helper chain as run_mode_b_local).
     * Step 2c.3 step 4: per-subgroup mask threaded via sgr.j_type_bitmask
     * for multi-subgroup walker support. */
    std::vector<std::vector<int>> cand_modeB;
    collect_candidates_pre_drift<Spec>(sub, radii_compacted.data(),
                                         (unsigned int)sgr.j_type_bitmask,
                                         DispatchPath::ModeB_HostWalker, cand_modeB);
    lazy_drift_candidates<Spec>(cand_modeB);
    evaluate_pairs_post_drift<Spec>(drv.ctx, actives_compacted.data(), n_compacted,
                                      cand_modeB, accums_compacted.data());

    /* Scatter compacted accums back into driver-owned per-slot accum_uvm.
     * Slots NOT in active_set_uvm keep their stale values (will not be
     * read this iter — after_iter only fires on still-active slots). */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        drv.accum_uvm[sg][slot] = accums_compacted[k];
    }
}

/* ============================================================================
 * nlr_iter_dispatch_subgroup_oracle_b_local<Spec>(drv, sg) — step 2c.4.
 *
 * Independent brute trajectory for the iterative oracle on Mode B local
 * (NTask==1). Mirrors nlr_iter_dispatch_subgroup_mode_b_local exactly,
 * substituting:
 *   - oracle's compacted active_set + radii (independent trajectory).
 *   - DispatchPath::Brute_Oracle backend in collect_candidates_pre_drift.
 *   - drv.ctx_oracle (brute-pass-guarded via NlrOracleBrutePassGuard) for
 *     evaluate_pairs_post_drift.
 *   - drv.accum_oracle_uvm[sg] as the output buffer (no apply_active_writeback).
 *
 * Per-iter brute-FIRST ordering enforced by the OUTER dispatch loop —
 * caller invokes this helper BEFORE production dispatch (G1 / 3d.3 fix).
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_oracle_b_local(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData = typename Spec::ActiveData;
    using AccumData  = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_oracle_size[sg];
    if (n_compacted <= 0) return;

    std::vector<int> active_particle_indices(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_oracle_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
    }
    neighbor_loop_args sub = drv.args;
    sub.active_list = active_particle_indices.data();
    sub.num_active  = n_compacted;

    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        radii_compacted[k] = drv.radii_oracle_uvm[sg][drv.active_set_oracle_uvm[sg][k]];
    }

    /* RAII brute-pass guard: set_oracle_brute_pass(ctx_oracle, true) on
     * construction, false on scope exit (any path, including endrun).
     * SFINAE-no-op for Specs without the hook. */
    NlrOracleBrutePassGuard<Spec> brute_guard(drv.ctx_oracle);

    std::vector<ActiveData> actives_compacted(n_compacted);
    build_self_actives_host_pre_drift<Spec>(sub, drv.ctx_oracle, radii_compacted.data(),
                                              drv.cs, actives_compacted.data());

    std::vector<AccumData> accums_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);
    }

    std::vector<std::vector<int>> cand_brute;
    collect_candidates_pre_drift<Spec>(sub, radii_compacted.data(),
                                         (unsigned int)sgr.j_type_bitmask,
                                         DispatchPath::Brute_Oracle, cand_brute);
    lazy_drift_candidates<Spec>(cand_brute);
    evaluate_pairs_post_drift<Spec>(drv.ctx_oracle, actives_compacted.data(), n_compacted,
                                      cand_brute, accums_compacted.data());

    /* Scatter brute results back into driver's per-slot oracle AccumData
     * (slot keyspace is the same as production — sgr.active_indices). */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_oracle_uvm[sg][k];
        drv.accum_oracle_uvm[sg][slot] = accums_compacted[k];
    }
    /* brute_guard destructor flips set_oracle_brute_pass(false) here. */
}

/* Local Mode-B iterative oracle combined dispatch.
 *
 * The separate brute-first helper is not sufficient for local iterative
 * density: its lazy drift mutates P[] before production snapshots ActiveData,
 * so oracle and production compare different i-side states. Keep the legacy
 * Mode-B epoch contract explicit here: snapshot both i-side active arrays and
 * collect both candidate lists before either path drifts j-side candidates.
 */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_local_with_oracle(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData = typename Spec::ActiveData;
    using AccumData  = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];

    const int n_prod = drv.active_set_size[sg];
    std::vector<int> active_particle_indices_prod(n_prod > 0 ? n_prod : 0);
    std::vector<double> radii_compacted_prod(n_prod > 0 ? n_prod : 0);
    for(int k = 0; k < n_prod; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices_prod[k] = sgr.active_indices[slot];
        radii_compacted_prod[k] = drv.radii_uvm[sg][slot];
    }
    neighbor_loop_args sub_prod = drv.args;
    sub_prod.active_list = (n_prod > 0) ? active_particle_indices_prod.data() : nullptr;
    sub_prod.num_active  = n_prod;

    const int n_oracle = drv.active_set_oracle_size[sg];
    std::vector<int> active_particle_indices_oracle(n_oracle > 0 ? n_oracle : 0);
    std::vector<double> radii_compacted_oracle(n_oracle > 0 ? n_oracle : 0);
    for(int k = 0; k < n_oracle; k++) {
        int slot = drv.active_set_oracle_uvm[sg][k];
        active_particle_indices_oracle[k] = sgr.active_indices[slot];
        radii_compacted_oracle[k] = drv.radii_oracle_uvm[sg][slot];
    }
    neighbor_loop_args sub_oracle = drv.args;
    sub_oracle.active_list = (n_oracle > 0) ? active_particle_indices_oracle.data() : nullptr;
    sub_oracle.num_active  = n_oracle;

    std::vector<ActiveData> actives_prod(n_prod > 0 ? n_prod : 0);
    if(n_prod > 0) {
        build_self_actives_host_pre_drift<Spec>(sub_prod, drv.ctx,
                                                radii_compacted_prod.data(),
                                                drv.cs, actives_prod.data());
    }

    std::vector<ActiveData> actives_oracle(n_oracle > 0 ? n_oracle : 0);
    if(n_oracle > 0) {
        build_self_actives_host_pre_drift<Spec>(sub_oracle, drv.ctx_oracle,
                                                radii_compacted_oracle.data(),
                                                drv.cs, actives_oracle.data());
    }

    std::vector<AccumData> accums_prod(n_prod > 0 ? n_prod : 0);
    for(int k = 0; k < n_prod; k++) Spec::zero_accum(accums_prod[k]);
    std::vector<AccumData> accums_oracle(n_oracle > 0 ? n_oracle : 0);
    for(int k = 0; k < n_oracle; k++) Spec::zero_accum(accums_oracle[k]);

    std::vector<std::vector<int>> cand_prod;
    collect_candidates_pre_drift<Spec>(sub_prod, radii_compacted_prod.data(),
                                       (unsigned int)sgr.j_type_bitmask,
                                       DispatchPath::ModeB_HostWalker, cand_prod);

    std::vector<std::vector<int>> cand_oracle;
    collect_candidates_pre_drift<Spec>(sub_oracle, radii_compacted_oracle.data(),
                                       (unsigned int)sgr.j_type_bitmask,
                                       DispatchPath::Brute_Oracle, cand_oracle);

    lazy_drift_candidates<Spec>(cand_prod);
    lazy_drift_candidates<Spec>(cand_oracle);

    if(n_prod > 0) {
        evaluate_pairs_post_drift<Spec>(drv.ctx, actives_prod.data(), n_prod,
                                        cand_prod, accums_prod.data());
        for(int k = 0; k < n_prod; k++) {
            int slot = drv.active_set_uvm[sg][k];
            drv.accum_uvm[sg][slot] = accums_prod[k];
        }
    }

    if(n_oracle > 0) {
        NlrOracleBrutePassGuard<Spec> brute_guard(drv.ctx_oracle);
        evaluate_pairs_post_drift<Spec>(drv.ctx_oracle, actives_oracle.data(), n_oracle,
                                        cand_oracle, accums_oracle.data());
        for(int k = 0; k < n_oracle; k++) {
            int slot = drv.active_set_oracle_uvm[sg][k];
            drv.accum_oracle_uvm[sg][slot] = accums_oracle[k];
        }
    }

    /* TEMPORARY oracle accum comparison from host vectors — UVM write broken
     * for Vel fields (same root cause as remote path).  Build slot→accum maps
     * for both prod and oracle, then compare over the production active set.
     * Remove with oracle teardown. */
    if(n_prod > 0) {
        std::unordered_map<int,int> oracle_k_by_slot;
        oracle_k_by_slot.reserve(n_oracle);
        for(int k = 0; k < n_oracle; k++) {
            int slot = drv.active_set_oracle_uvm[sg][k];
            oracle_k_by_slot[slot] = k;
        }
        char origin_tag[40];
        std::snprintf(origin_tag, sizeof(origin_tag), "iter_sg%d_it%d", sg, drv.iter_index);
        for(int k = 0; k < n_prod; k++) {
            int slot = drv.active_set_uvm[sg][k];
            auto it = oracle_k_by_slot.find(slot);
            if(it != oracle_k_by_slot.end()) {
                emit_oracle_mismatch_if_any<Spec>(ThisTask, slot,
                    accums_prod[k], accums_oracle[it->second],
                    origin_tag, &drv.oracle_mismatch_count);
            }
        }
        drv.oracle_accum_compared_in_dispatch = true;
    }
}

/* ============================================================================
 * nlr_iter_dispatch_subgroup_oracle_b_remote<Spec>(drv, sg) — step 2c.4.
 *
 * Independent brute trajectory for the iterative oracle on Mode B remote
 * (NTask>1). Reuses the extracted helper mode_b_remote_evaluate_into_buffer
 * with RemoteHelperMode::OracleBrutePass — same epoch order (queries /
 * collect-brute / drift / brute walk / exchange replies / merge), brute
 * backend throughout. Caller-owned brute-pass guard at this scope.
 *
 * SSOT guardrail (codex 2c.4): no duplicated remote body. The helper
 * handles self+peer brute walks and replies-merge via the existing
 * Stage 4 / 10 / 11 machinery.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_oracle_b_remote(NlrIterDriver<Spec>& drv, int sg)
{
    using AccumData = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_oracle_size[sg];
    /* NOTE: helper is COLLECTIVE — must be entered on every rank for the
     * MPI Allgather/Alltoall stages. Empty-rank participation uses
     * n_compacted=0 with nullptr buffers (codex 2b.2 invariant 6 +
     * 2c.2-fix-v2 blocker 23). */

    std::vector<int> active_particle_indices(n_compacted > 0 ? n_compacted : 0);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_oracle_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
    }
    neighbor_loop_args sub = drv.args;
    sub.active_list = (n_compacted > 0) ? active_particle_indices.data() : nullptr;
    sub.num_active  = n_compacted;

    std::vector<double> radii_compacted(n_compacted > 0 ? n_compacted : 0);
    for (int k = 0; k < n_compacted; k++) {
        radii_compacted[k] = drv.radii_oracle_uvm[sg][drv.active_set_oracle_uvm[sg][k]];
    }

    /* Caller-owned brute-pass guard scopes the entire helper call. */
    NlrOracleBrutePassGuard<Spec> brute_guard(drv.ctx_oracle);

    /* Driver-owned compacted AccumData buffer; scatter back to per-slot
     * oracle AccumData after helper returns. */
    std::vector<AccumData> accums_compacted(n_compacted > 0 ? n_compacted : 0);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);
    }

    mode_b_remote_evaluate_into_buffer<Spec, RemoteHelperMode::OracleBrutePass>(
        sub,
        (n_compacted > 0) ? radii_compacted.data() : nullptr,
        drv.cs,
        drv.ctx_oracle,
        (unsigned int)sgr.j_type_bitmask,
        (n_compacted > 0) ? accums_compacted.data() : nullptr,
        /*actives_out=*/nullptr,
        /*tim=*/nullptr);

    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_oracle_uvm[sg][k];
        drv.accum_oracle_uvm[sg][slot] = accums_compacted[k];
    }
    /* brute_guard destructor flips set_oracle_brute_pass(false) here. */
}

/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_b_remote_with_oracle<Spec>(drv, sg)
 *
 * Combined single-epoch iterative oracle for Mode B remote (NTask>1).
 * Replaces the two-call sequence (oracle_b_remote + mode_b_remote) that
 * ran two independent helper round-trips, causing an epoch skew: the brute
 * oracle drifted candidates before the production tree even collected its
 * candidate set.
 *
 * Uses RemoteHelperMode::OracleIterative to collect tree + brute candidate
 * sets in ONE shared epoch (pre-drift both, drift union, evaluate brute-
 * first then tree), receiving tree -> drv.accum_uvm and brute ->
 * drv.accum_oracle_uvm.  Both use the PRODUCTION active set and radii so
 * the oracle comparison is "same radius, same epoch: does tree match brute?"
 *
 * Oracle radii are synced to production radii after each iter so that the
 * oracle after_iter (b.oracle) runs on the same h_search as production.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_remote_with_oracle(
    NlrIterDriver<Spec>& drv, int sg)
{
    using AccumData = typename Spec::AccumData;

    const NlrSubgroup& sgr    = drv.args.subgroups[sg];
    const int n_compacted     = drv.active_set_size[sg];

    std::vector<int>    active_particle_indices(n_compacted);
    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot                   = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
        radii_compacted[k]         = drv.radii_uvm[sg][slot];
    }

    neighbor_loop_args sub = drv.args;
    sub.active_list = (n_compacted > 0) ? active_particle_indices.data() : nullptr;
    sub.num_active  = n_compacted;

    std::vector<AccumData> accums_prod(n_compacted);
    std::vector<AccumData> accums_oracle(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_prod[k]);
        Spec::zero_accum(accums_oracle[k]);
    }

    mode_b_remote_evaluate_into_buffer<Spec, RemoteHelperMode::OracleIterative>(
        sub,
        (n_compacted > 0) ? radii_compacted.data() : nullptr,
        drv.cs,
        drv.ctx,
        (unsigned int)sgr.j_type_bitmask,
        (n_compacted > 0) ? accums_prod.data()   : nullptr,
        /*actives_out=*/nullptr,
        /*tim=*/nullptr,
        (n_compacted > 0) ? accums_oracle.data() : nullptr,
        &drv.ctx_oracle);

    /* TEMPORARY oracle accum comparison — host vectors known-correct (oracle UVM
     * mismatches); oracle UVM write is broken for Vel fields (Probe3 shows Vel=0
     * despite correct host source).  Compare here from host vectors and set flag
     * so b.compare skips the UVM-based emit_oracle_mismatch_if_any.
     * Remove with oracle teardown. */
    {
        char origin_tag[40];
        std::snprintf(origin_tag, sizeof(origin_tag), "iter_sg%d_it%d", sg, drv.iter_index);
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            emit_oracle_mismatch_if_any<Spec>(ThisTask, slot,
                accums_prod[k], accums_oracle[k], origin_tag, &drv.oracle_mismatch_count);
        }
        drv.oracle_accum_compared_in_dispatch = true;
    }
    /* Scatter both into driver-owned per-slot buffers.  Sync oracle radii to
     * production radii so oracle after_iter uses the same h_search. */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        drv.accum_uvm[sg][slot]        = accums_prod[k];
        drv.accum_oracle_uvm[sg][slot] = accums_oracle[k];
        drv.radii_oracle_uvm[sg][slot] = drv.radii_uvm[sg][slot];
    }
    /* Sync oracle active set to production active set for this iteration so
     * (b.oracle) after_iter iterates the same slots. */
    drv.active_set_oracle_size[sg] = n_compacted;
    for (int k = 0; k < n_compacted; k++) {
        drv.active_set_oracle_uvm[sg][k] = drv.active_set_uvm[sg][k];
    }
}

/* ============================================================================
 * Default on_max_iter_exceeded — runner-supplied hard endrun(1155).
 * Matches legacy hydro/density.cc:602 and gravity/ags_rkern.cc:414 policy.
 * Specs override via `static void on_max_iter_exceeded(const NlrIterDriver<Spec>&);`
 * only when port's legacy policy differs and Phil approved (rare).
 * ========================================================================== */
template <typename Spec>
static void nlr_default_on_max_iter_exceeded(const NlrIterDriver<Spec>& drv)
{
    if (ThisTask == 0) {
        fprintf(stderr,
            "[%s] FATAL: failed to converge in %d iterations (max_iters=%d). "
            "global_active_total=%d remained.\n",
            Spec::loop_name, drv.iter_index, Spec::max_iters, drv.global_active_total);
        fflush(stderr);
    }
    endrun(1155);
}

/* ============================================================================
 * run_neighbor_loop_iterative<Spec> — body (step 2b).
 * ========================================================================== */
template <typename Spec>
void run_neighbor_loop_iterative(const neighbor_loop_args_iterative& args)
{
#ifdef GIZMO_NLR_ITER_HARNESS_TEST
    if constexpr (std::is_same_v<Spec, IterHarnessSpec> ||
                  std::is_same_v<Spec, IterHarnessGhostSpec>) {
        g_iter_harness_telemetry = IterHarnessTelemetry{};
    }
#endif

    /* ===== Compile-time spec consistency ===== */
    static_assert(std::is_same_v<typename Spec::IterControl, Iterative>,
                  "run_neighbor_loop_iterative requires Spec::IterControl = Iterative. "
                  "NotIterative Specs must call run_neighbor_loop instead.");
    static_assert(nlr_spec_has_after_iter_v<Spec>,
                  "Iterative Spec is missing `static IterResult after_iter(const AfterIterContext<Spec>&, const AccumData&);`. "
                  "Required when IterControl = Iterative. See OPEN_phase4_b0_iterative_design.md §3.");
    static_assert(std::is_trivially_copyable_v<typename Spec::IterScratch>,
                  "Spec::IterScratch must be std::is_trivially_copyable_v. "
                  "Use `using IterScratch = NoIterScratch;` for the empty form.");
    static_assert(Spec::max_iters >= 1,
                  "Spec::max_iters must be >= 1.");
    static_assert(Spec::mode_a_csr_buffer_factor > 1.0,
                  "Spec::mode_a_csr_buffer_factor must be > 1.0 "
                  "(legacy DENSITY_H_BUFFER_FACTOR = 1.3).");
    /* TRAP-5 carry-forward (codex 2b review 2026-05-10): same trivially-copyable
     * checks as run_neighbor_loop. Don't let iterative Specs bypass TRAP 5. */
    static_assert(std::is_trivially_copyable_v<typename Spec::CallScalars>,
                  "Spec::CallScalars must be trivially copyable (lambda-captured by value).");
    static_assert(std::is_trivially_copyable_v<typename Spec::ActiveData>,
                  "Spec::ActiveData must be trivially copyable (UVM-staged).");
    static_assert(std::is_trivially_copyable_v<typename Spec::NeighborData>,
                  "Spec::NeighborData must be trivially copyable (built per-pair on device).");
    static_assert(std::is_trivially_copyable_v<typename Spec::AccumData>,
                  "Spec::AccumData must be trivially copyable (UVM-staged).");

    /* ===== Runtime checks BEFORE any state touch ===== */
    if (args.num_subgroups < 1) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: num_subgroups=%d < 1. "
                "Caller must short-circuit when global active total is 0.\n",
                Spec::loop_name, args.num_subgroups);
            fflush(stderr);
        }
        endrun(81200);
    }
    if (args.num_subgroups > 1 && !nlr_supports_subgroups<Spec>::value) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: num_subgroups=%d > 1 "
                "but Spec did not declare `using SupportsSubgroups = std::true_type;` (TRAP 9).\n",
                Spec::loop_name, args.num_subgroups);
            fflush(stderr);
        }
        endrun(81201);
    }
    /* Multi-subgroup walker wired in step 2c.3:
     *   - sg.j_type_bitmask threaded through Mode B local + remote collectors
     *     and Mode A NGL build per-subgroup.
     *   - Mode A multi-subgroup ghost import uses union semantics (mask
     *     OR + concatenation of actives/radii across subgroups). Documented
     *     as deliberate over-import for cross-bm cases (codex plan v2 ch. 6).
     *   - Per-subgroup global activity tracking via global_active_per_sg[].
     *   - Pre-dispatch invalidation sweep rebuilds CSR caches once per iter
     *     when any subgroup invalidated (codex plan v2 tightening 1).
     *   - Multi-subgroup actives partition + collective-symmetry validated
     *     under GIZMO_NLR_SUBGROUP_AUDIT=1 (step 10 audit block above).
     * Per-Spec opt-in still REQUIRED via `using SupportsSubgroups = std::true_type;`
     * (TRAP 9; runtime abort 81201 above catches missing trait). */

    /* ===== Oracle env detection (step 2c.4) =====
     * Mode B oracle paths implemented this slice; Mode A + oracle still
     * hard-stubbed (codex 2c.4-v3: Mode A oracle would catch only
     * cached-CSR/lookup/rebuild bugs (covered by the synthetic harness in
     * step 3) and cannot validate ghost-import completeness (walks the
     * same imported pool as production). Its complexity is not justified
     * for temporary scaffolding code). */
    const bool oracle_enabled = gizmo_nlr_oracle_enabled_global() ||
                                gizmo_nlr_oracle_enabled_for(Spec::loop_name);

    /* ===== Path selection at iter 0 (FIXED for whole call per v3 §0) =====
     * Step 2c.2: integrates with the canonical force-mode / threshold dispatch
     * (mirrors run_neighbor_loop logic at line 1936+). Path is held fixed
     * across all iterations of one iterative call.
     *
     * num_active for threshold uses the UNION across all subgroups (base
     * args.num_active per the doc convention). 2c.2 is single-subgroup-only;
     * 2c.3 multi-subgroup walker mask-threading lands separately. */
    const NlrForceMode force_mode = gizmo_nlr_force_mode_for(Spec::loop_name);
    DispatchPath path;
    int forced_modeb_global_active = -1;
    if (force_mode == NlrForceMode::A) {
        path = DispatchPath::ModeA_GPU_NGL;
    } else if (force_mode == NlrForceMode::B) {
        int local_act = args.num_active;
        MPI_Allreduce(&local_act, &forced_modeb_global_active, 1,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        nlr_abort_if_forced_modeb_too_large(
            Spec::loop_name, args.num_active, forced_modeb_global_active);
        path = DispatchPath::ModeB_HostWalker;
    } else {
        /* Threshold dispatch: Allreduce sum + max of base args.num_active
         * (= union across subgroups). */
        int local_act = args.num_active;
        int sum_act = 0, max_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_act, &max_act, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        const int spec_default_sum = nlr_spec_threshold_sum<Spec>(64);
        const int spec_default_max = nlr_spec_threshold_max<Spec>(64);
        const int TS = gizmo_nlr_modeb_threshold_sum_for(Spec::loop_name, spec_default_sum);
        const int TM = gizmo_nlr_modeb_threshold_max_for(Spec::loop_name, spec_default_max);
        bool select_mode_b = (sum_act > 0) && (sum_act <= TS) && (max_act <= TM);
        path = select_mode_b ? DispatchPath::ModeB_HostWalker : DispatchPath::ModeA_GPU_NGL;
    }

    /* ===== Mode A + oracle hard-stub (step 2c.4) =====
     * Mode A iterative production path is fully ported; only the *oracle*
     * backend on Mode A is intentionally not implemented. Oracle is
     * temporary port-validation scaffolding; Mode A oracle over the
     * imported pool catches only CSR/cached-lookup/rebuild bugs (covered
     * by the synthetic harness in step 3) and cannot validate ghost-import
     * completeness (walks the same imported pool as production). For
     * oracle validation force Mode B: GIZMO_NLR_FORCE_MODE=B GIZMO_NLR_ORACLE=1. */
    if (oracle_enabled && path == DispatchPath::ModeA_GPU_NGL) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: Mode A iterative "
                "oracle intentionally not implemented. The Mode A *production* "
                "path is fully ported and validated separately. Mode A oracle "
                "would only catch cached-CSR/lookup/rebuild bugs (covered by "
                "the synthetic harness in step 3) and cannot validate ghost-"
                "import completeness (walks the same imported pool as "
                "production). For oracle validation, force Mode B: "
                "GIZMO_NLR_FORCE_MODE=B GIZMO_NLR_ORACLE=1.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        endrun(81203);
    }

    if(gizmo_nlr_dispatch_trace_enabled()) {
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) {
            const char *src =
                (force_mode == NlrForceMode::A || force_mode == NlrForceMode::B)
                ? "force" : "threshold";
            fprintf(stderr,
                    "[NLR ITER DISPATCH caller=%s path=%s (%s) NTask=%d "
                    "local_active=%d forced_modeb_global_active=%d oracle=%d]\n",
                    Spec::loop_name,
                    path == DispatchPath::ModeA_GPU_NGL ? "gpu_ngl" : "mode_b",
                    src, NTask, args.num_active,
                    forced_modeb_global_active, (int)oracle_enabled);
            fflush(stderr);
        }
    }

    /* ===== CallScalars captured ONCE for whole call (codex constraint 4) ===== */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    /* ===== Step 2c.3 step 10: GIZMO_NLR_SUBGROUP_AUDIT collective-symmetry +
     * no-duplicate-actives runtime checks. Hard-abort on violation. Off by
     * default; harness / test-mode enables.
     *
     * Check 1: subgroups[] length AND ordered bm-key sequence identical
     * across all ranks (codex plan v2 tightening 4 — MIN/MAX hash, NOT BOR).
     *
     * Check 2: each particle index appears in AT MOST ONE subgroup
     * (pitfall 6 — no duplicate active ownership; would double-apply final
     * writeback state). */
    if (gizmo_nlr_subgroup_audit_enabled()) {
        /* Check 1: 64-bit FNV-1a ordered hash of (i << 32) | bm_key. */
        uint64_t local_hash = 0xCBF29CE484222325ULL;
        for (int i = 0; i < args.num_subgroups; i++) {
            local_hash ^= ((uint64_t)(unsigned int)i << 32) |
                          (uint64_t)args.subgroups[i].j_type_bitmask;
            local_hash *= 0x100000001B3ULL;
        }
        int      local_n = args.num_subgroups;
        int      n_min = local_n, n_max = local_n;
        uint64_t h_min = local_hash, h_max = local_hash;
        if (NTask > 1) {
            MPI_Allreduce(&local_n,    &n_min, 1, MPI_INT,      MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&local_n,    &n_max, 1, MPI_INT,      MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&local_hash, &h_min, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&local_hash, &h_max, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
        }
        if (n_min != n_max || h_min != h_max) {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                "[NLR_ITER SUBGROUP_AUDIT ABORT rank=%d caller=%s] subgroups[] "
                "length/order mismatch across ranks: local_num=%d (global "
                "min=%d max=%d), local_hash=%llx (global min=%llx max=%llx). "
                "Caller must fill subgroups[] from global_bm_presence union "
                "with identical ordering on all ranks.\n",
                rank, Spec::loop_name, local_n, n_min, n_max,
                (unsigned long long)local_hash,
                (unsigned long long)h_min, (unsigned long long)h_max);
            fprintf(stderr, "  rank=%d subgroups[] bm_keys: [", rank);
            for (int i = 0; i < args.num_subgroups; i++) {
                fprintf(stderr, "%s%u", (i > 0 ? "," : ""),
                        args.subgroups[i].j_type_bitmask);
            }
            fprintf(stderr, "]\n");
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 81214);
        }

        /* Check 2: no duplicate particle indices across subgroups (local-rank
         * only — cross-rank is governed by caller's partition logic). */
        {
            std::vector<int> seen;
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                const NlrSubgroup& sgr = args.subgroups[sg];
                for (int k = 0; k < sgr.num_active_local; k++) {
                    seen.push_back(sgr.active_indices[k]);
                }
            }
            std::sort(seen.begin(), seen.end());
            for (size_t k = 1; k < seen.size(); k++) {
                if (seen[k] == seen[k-1]) {
                    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    fprintf(stderr,
                        "[NLR_ITER SUBGROUP_AUDIT ABORT rank=%d caller=%s] "
                        "particle index %d appears in 2+ subgroups. Caller "
                        "partition must satisfy no-duplicate-active-ownership "
                        "invariant (pitfall 6) — each particle in AT MOST ONE "
                        "subgroup.\n",
                        rank, Spec::loop_name, seen[k]);
                    fflush(stderr);
                    MPI_Abort(MPI_COMM_WORLD, 81215);
                }
            }
        }
    }

    /* ===== Mode B hard-corridor counter snapshot (codex 2c.2-fix 2026-05-10) =====
     * Mode B paths MUST NOT enter move_particles / ghost_exchange_impl /
     * gpu_particles_arena_acquire, and MUST NOT mutate NumPart. Always-on
     * enforcement around the iterative call; mirrors non-iter run_neighbor_loop
     * lines 2018-2026 + 2099-2122. Snapshot here; check post-iter-loop.
     * Mode A paths legitimately advance these counters (one ghost_import +
     * one arena_acquire per call); enforcement skipped. */
    const uint64_t s_drift0 = g_global_drift_counter;
    const uint64_t s_ghost0 = g_ghost_import_counter;
    const uint64_t s_arena0 = g_gpu_arena_acquire_counter;
    const int      s_np0    = NumPart;

    /* ===== Driver lifetime wrapped in inner scope (codex 2c.2-fix-v2 hardening 1) =====
     * Driver destructor runs at scope exit; Mode B hard-corridor check fires
     * AFTER scope exit so the check covers driver cleanup hooks too (defensive
     * — current cleanup_device_context contractually MUST NOT touch
     * drift/ghost/arena globals, but inner-scope wrapping makes the invariant
     * maximally airtight). */
    {
    NlrIterDriver<Spec> drv(args, cs);

    /* ===== Path-specific DeviceContext init (step 2c.1 + 2c.2) =====
     * Mode B: bind ctx.P=args.P/CellP directly (no arena).
     * Mode A: acquire_arena_and_init_ctx_mode_a — arena_acquire ONCE per call
     *         then bind ctx to arena-resident P_gpu/CellP_gpu. */
    if (path == DispatchPath::ModeB_HostWalker) {
        drv.initialize_device_context_mode_b();
    } else if (path == DispatchPath::ModeA_GPU_NGL) {
        drv.acquire_arena_and_init_ctx_mode_a();
    } else {
        /* Brute_Oracle or future path — not reachable in step 2c.2. */
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: unhandled path %d.\n",
                Spec::loop_name, (int)path);
            fflush(stderr);
        }
        endrun(81208);
    }

    /* ===== Independent ctx_oracle init (step 2c.4) =====
     * Codex 2c.4-v2 P1 fix: ctx_oracle gets its OWN populate_device_context
     * + cleanup_device_context, never aliased to production ctx. Bound to
     * args.P/CellP (Mode B lazy-drift contract — production reads owner-
     * local args.P[j] directly, oracle reads the same args.* with the
     * brute backend). Only reached on Mode B (Mode A + oracle hard-stubbed
     * above). */
    if (drv.oracle_enabled) {
        drv.ctx_oracle.P         = drv.args.P;
        drv.ctx_oracle.CellP     = drv.args.CellP;
        drv.ctx_oracle.num_total = drv.args.num_total;
        if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
            Spec::populate_device_context(drv.args, drv.ctx_oracle);
        }
        drv.ctx_oracle_initialized = true;
    }

    /* ===== Outer iter loop =====
     *
     * INVARIANT (codex v4.4 fix 2026-05-10): `accum_uvm[sg][slot]` always
     * holds the most recent evaluated AccumData for that slot. Converged
     * slots are NEVER re-touched after they leave active_set; their final
     * iter's result persists in accum_uvm until the post-loop apply_active_writeback
     * reads it.
     *
     * Per-iter zero of AccumData is therefore the dispatch helper's
     * responsibility (it zeroes via Spec::zero_accum on ONLY the compacted
     * active subset before pair_kernel evaluation; see
     * nlr_iter_dispatch_subgroup_mode_b_local).
     * The runner does NOT zero converged slots — doing so would destroy
     * exactly the result apply_active_writeback needs. The earlier "defensive
     * cross-iter zero of all slots" in this position was the bug codex
     * caught in 2b review. */
    /* Phase 4.B.0 v4.3 step 0 (partition-by-subgroup contract): record the
     * iter-0 per-subgroup j_type_bitmask values so the per-iter partition
     * assertion below can verify Spec::active_subgroup_key returns the same
     * key for every active in that subgroup on every iter. Only populated
     * when the Spec opts into actives_partition_by_subgroup. Compile-time
     * no-op for Specs that don't. */
    std::vector<unsigned int> partition_expected_bm_key;
    if constexpr (nlr_spec_actives_partition_by_subgroup_v<Spec>) {
        static_assert(nlr_spec_has_active_subgroup_key_v<Spec>,
                      "Spec::actives_partition_by_subgroup=true requires "
                      "Spec::active_subgroup_key(const DeviceContext&, int, "
                      "const CallScalars&) returning int (= bm key). See "
                      "OPEN_3d_agsdensity_design.md §5a.");
        partition_expected_bm_key.resize(args.num_subgroups);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            partition_expected_bm_key[sg] = args.subgroups[sg].j_type_bitmask;
        }
    }

    for (drv.iter_index = 0; drv.iter_index < Spec::max_iters; drv.iter_index++) {

        /* Phase 4.B.0 v4.3 step 0 — mode_a_rebuild_csr_every_iter correctness
         * fallback (OPEN_3d_agsdensity_design.md §10.2). When the Spec sets
         * this trait true, force-invalidate every subgroup's Mode A CSR
         * cache at the start of every iter > 0, bypassing the buffer-
         * exceedance trigger entirely. The static_assert on
         * mode_a_csr_buffer_factor > 1.0 is unchanged (factor stays a valid
         * number even when unused). Iter 0 is exempt — initial CSR build
         * happens in the per-subgroup dispatch on first invalidation pass. */
        if constexpr (nlr_spec_mode_a_rebuild_csr_every_iter_v<Spec>) {
            if (path == DispatchPath::ModeA_GPU_NGL && drv.iter_index > 0) {
                for (int sg = 0; sg < args.num_subgroups; sg++) {
                    drv.mode_a_csr_valid[sg] = false;
                }
            }
        }

        /* Phase 4.B.0 v4.3 step 0 — partition-by-subgroup debug assertion
         * (OPEN_3d_agsdensity_design.md §5a). Gated under DEBUG or the
         * dedicated env-macro so production builds carry zero overhead.
         * Calls the path-correct ctx accessor via Spec::active_subgroup_key
         * for every active in every globally-active subgroup; mismatch with
         * the iter-0 recorded j_type_bitmask = terminate(). Host-side check
         * before per-iter dispatch. */
#if defined(DEBUG) || defined(GIZMO_NLR_ASSERT_PARTITION)
        if constexpr (nlr_spec_actives_partition_by_subgroup_v<Spec>) {
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                /* Skip globally-converged subgroups (iter > 0 only). */
                if (drv.iter_index > 0 && drv.global_active_per_sg[sg] <= 0) continue;
                const int n_compacted = drv.active_set_size[sg];
                const NlrSubgroup& sgr = args.subgroups[sg];
                const unsigned int expect = partition_expected_bm_key[sg];
                for (int k = 0; k < n_compacted; k++) {
                    const int slot = drv.active_set_uvm[sg][k];
                    const int i    = sgr.active_indices[slot];
                    const int key  = Spec::active_subgroup_key(drv.ctx, i, drv.cs);
                    if (static_cast<unsigned int>(key) != expect) {
                        char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
                        snprintf(buf, sizeof(buf),
                            "[run_neighbor_loop_iterative<%s>] FATAL: Spec contract "
                            "violation: actives_partition_by_subgroup key drift. "
                            "sg=%d iter=%d slot=%d i=%d expected_bm=%u got_key=%d. "
                            "active_subgroup_key MUST be a pure function of "
                            "state that does not change across iterations. See "
                            "OPEN_3d_agsdensity_design.md §5a.\n",
                            Spec::loop_name, sg, drv.iter_index, slot, i,
                            expect, key);
                        terminate(buf);
                    }
                }
            }
        }
#endif

        /* (a-pre) Mode A pre-dispatch invalidation sweep (step 2c.3 step 8).
         * If any subgroup's CSR is invalid (set by the buffer-exceedance
         * trigger in last iter's step b.5), call the no-arg union-rebuild
         * method ONCE before any per-subgroup dispatch this iter. Prevents
         * mixing old/new arena pools across subgroup dispatches in the same
         * iter (codex plan v2 tightening 1). Skips iter 0 (initial arena
         * acquire was via acquire_arena_and_init_ctx_mode_a). */
        if (path == DispatchPath::ModeA_GPU_NGL && drv.iter_index > 0) {
            /* Codex 2c.3 blocker #2 + #3 fix: the sweep must be COLLECTIVE.
             * mode_a_csr_valid is rank-local; a rank with no actives for a
             * globally-active subgroup never builds CSR locally, so its
             * mode_a_csr_valid[sg] stays false while another rank's is true.
             * Without an Allreduce, one rank enters ghost_exchange_cleanup +
             * reimport collectives while the other skips them => deadlock.
             * Also skip globally-converged subgroups (blocker #3): a sg with
             * global_active_per_sg[sg]==0 may have invalid CSR from an
             * earlier buffer-exceedance trigger but doesn't need a rebuild. */
            int local_needs_rebuild = 0;
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                if (drv.global_active_per_sg[sg] <= 0) continue;
                /* Codex 2c.3 post-blocker-fix re-review (2026-05-10): a rank
                 * with zero local actives in this globally-active subgroup
                 * never builds CSR locally, so its mode_a_csr_valid[sg]
                 * stays false. Without this guard, such a rank would vote
                 * "rebuild" every iteration forever — Allreduce would
                 * force a global rebuild that doesn't actually fix the
                 * vote, infinite loop of useless global rebuilds. The rank
                 * still participates collectively in the Allreduce below
                 * (voting 0); it just doesn't request a rebuild for a CSR
                 * it doesn't need. */
                if (drv.active_set_size[sg] <= 0) continue;
                if (!drv.mode_a_csr_valid[sg]) {
                    local_needs_rebuild = 1;
                    break;
                }
            }
            int global_needs_rebuild = local_needs_rebuild;
            if (NTask > 1) {
                MPI_Allreduce(&local_needs_rebuild, &global_needs_rebuild, 1,
                              MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (global_needs_rebuild) {
                drv.rebuild_mode_a_arena_and_ctx_for_current_active_union();
            }
        }

        /* (a0) Per-iter Spec hook: reset_per_iter_device_context (step 2c.1).
         * Optional. Runs HOST-side on every rank before subgroup dispatch.
         * Use case: ags_density's per_iter_wakeup_detected counter zero
         * (lands at 3d.4). Synthetic harness in step 3 may consume this
         * for validation. Hook MUST NOT do MPI (TRAP 7 generalized). */
        if constexpr (nlr_spec_has_reset_per_iter_device_context_v<Spec>) {
            Spec::reset_per_iter_device_context(args, drv.ctx, drv.iter_index);
            /* Codex 2c.4 post-review (2026-05-10): ctx_oracle has its own
             * independent lifecycle and carries the same per-iter counter
             * fields (e.g. ags-style wakeup counters, harness diagnostic
             * counters). Without the symmetric reset here, the oracle's
             * brute walk would read stale-from-last-iter state out of
             * ctx_oracle while production walks with freshly-reset state
             * out of ctx — guaranteed per-iter divergence on j-side
             * counter fields. Fires BEFORE oracle dispatch (a.oracle)
             * which is itself brute-first per iter (G1). Mode A is hard-
             * stubbed so no rebuild complication. */
            if (drv.oracle_enabled) {
                Spec::reset_per_iter_device_context(args, drv.ctx_oracle, drv.iter_index);
            }
        }

        /* (a) Per-subgroup dispatch — fixed path for the call. Step 2b.2
         * implements Mode B local (np=1) and Mode B remote (np>1) via the
         * same DispatchPath::ModeB_HostWalker label; the per-subgroup helper
         * picks local vs remote based on NTask. Mode A iter and iterative
         * oracle still hard-stubbed at the env / path-selection block. */
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            /* Step 2c.3 step 6: skip globally-converged subgroups. Saves
             * per-iter no-op collectives. On iter 0 global_active_per_sg
             * is still 0 (not yet Allreduced); use local count as proxy +
             * global presence-of-actives (subgroups[] is filled from
             * global_bm_presence by the caller, so any subgroup present
             * has SOME rank with non-zero actives — entering the dispatch
             * is the correct lockstep collective participation). */
            const bool sg_globally_active = (drv.iter_index == 0)
                ? true   /* iter 0: trust caller's subgroups[] (global union) */
                : (drv.global_active_per_sg[sg] > 0);
            if (!sg_globally_active) continue;

            /* (a.oracle) Step 2c.4 G1: brute-FIRST per iter. Independent
             * oracle trajectory; Mode B only (Mode A + oracle hard-stubbed
             * at outer entry). Helpers are collective on remote — must be
             * entered on every rank regardless of local n_compacted.
             *
             * Both single-rank and multi-rank paths use a COMBINED helper
             * that collects tree + brute candidate sets in one shared epoch
             * (pre-drift both → drift union → evaluate brute-first, then
             * tree).  The old two-call remote sequence was epoch-skewed:
             * nlr_iter_dispatch_subgroup_oracle_b_remote drifted candidates
             * before nlr_iter_dispatch_subgroup_mode_b_remote even collected
             * its set, causing spurious iter-sg0-it0 mismatches at the first
             * dm_dispersion call. */
            if (drv.oracle_enabled) {
                if (NTask == 1) {
                    nlr_iter_dispatch_subgroup_mode_b_local_with_oracle<Spec>(drv, sg);
                } else {
                    nlr_iter_dispatch_subgroup_mode_b_remote_with_oracle<Spec>(drv, sg);
                }
                continue;
            }

            switch (path) {
                case DispatchPath::ModeB_HostWalker:
                    if (NTask == 1) {
                        nlr_iter_dispatch_subgroup_mode_b_local<Spec>(drv, sg);
                    } else {
                        nlr_iter_dispatch_subgroup_mode_b_remote<Spec>(drv, sg);
                    }
                    break;
                case DispatchPath::ModeA_GPU_NGL:
                    nlr_iter_dispatch_subgroup_mode_a<Spec>(drv, sg);
                    break;
                case DispatchPath::Brute_Oracle:
                    /* Brute_Oracle is the oracle backend, NEVER selected as a
                     * production path. Defense-in-depth abort. */
                    if (ThisTask == 0) {
                        fprintf(stderr,
                            "[run_neighbor_loop_iterative<%s>] FATAL: Brute_Oracle "
                            "as production path at sg=%d iter=%d.\n",
                            Spec::loop_name, sg, drv.iter_index);
                        fflush(stderr);
                    }
                    endrun(81202);
                    break;
            }
        }

        /* Snapshot production active sets BEFORE after_iter compaction so
         * (b.compare) can compare only slots evaluated this iteration, not
         * stale slots from prior iterations whose accum_uvm entries were
         * never updated this iter. */
        std::vector<std::vector<int>> prod_active_snapshot;
        if (drv.oracle_enabled) {
            prod_active_snapshot.resize(args.num_subgroups);
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                int n = drv.active_set_size[sg];
                prod_active_snapshot[sg].assign(
                    drv.active_set_uvm[sg], drv.active_set_uvm[sg] + n);
            }
        }

        /* (b) Per-active Spec::after_iter — collect IterStatus, compact
         * active_set per subgroup, mutate radii on AdjustRadius. */
        drv.local_active_total = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            int n_compacted   = drv.active_set_size[sg];
            int write_idx     = 0;
            for (int k = 0; k < n_compacted; k++) {
                int  slot = drv.active_set_uvm[sg][k];
                int  i    = args.subgroups[sg].active_indices[slot];
                AfterIterContext<Spec> ctx{
                    args, sg, slot, i, drv.iter_index,
                    /* h_search_current */ drv.radii_uvm[sg][slot],
                    /* scalars (ref to driver-owned) */ drv.cs,
                    /* scratch (mutable ref) */          drv.scratch_uvm[sg][slot]
                };
                IterResult r = Spec::after_iter(ctx, drv.accum_uvm[sg][slot]);
                /* Step 2c.4 audit instrumentation: stash status + new_h for
                 * per-iter 4-thing compare against oracle trajectory. */
                if (drv.oracle_enabled) {
                    drv.status_prod_uvm[sg][slot] = static_cast<int>(r.status);
                    drv.new_h_prod_uvm[sg][slot]  =
                        (r.status == IterStatus::AdjustRadius)
                            ? r.new_h_search
                            : drv.radii_uvm[sg][slot];
                }
                switch (r.status) {
                    case IterStatus::Converged:
                        /* Drop slot from compacted active_set; its accum_uvm
                         * stays so apply_active_writeback can read it post-loop. */
                        break;
                    case IterStatus::NeedsMore:
                        drv.active_set_uvm[sg][write_idx++] = slot;
                        break;
                    case IterStatus::AdjustRadius:
                        drv.radii_uvm[sg][slot] = r.new_h_search;
                        drv.active_set_uvm[sg][write_idx++] = slot;
                        break;
                    default:
                        /* Bad enum value — silent drop-as-converged would
                         * mask a Spec bug. Hard abort with full context. */
                        if (ThisTask == 0) {
                            fprintf(stderr,
                                "[run_neighbor_loop_iterative<%s>] FATAL: Spec::after_iter "
                                "returned unknown IterStatus=%d at iter=%d sg=%d slot=%d "
                                "(particle index i=%d). Valid values: 0=Converged, "
                                "1=NeedsMore, 2=AdjustRadius.\n",
                                Spec::loop_name, (int)r.status,
                                drv.iter_index, sg, slot, i);
                            fflush(stderr);
                        }
                        endrun(81206);
                }
            }
            drv.active_set_size[sg]     = write_idx;
            drv.local_active_per_sg[sg] = write_idx;       /* step 2c.3 step 5 */
            drv.local_active_total     += write_idx;
        }

        /* (b.oracle) Step 2c.4: oracle after_iter — symmetric pattern over
         * the brute trajectory's state. Mutates oracle's radii / scratch /
         * active_set and stashes status + new_h into status_oracle_uvm /
         * new_h_oracle_uvm for the 4-thing compare below. */
        if (drv.oracle_enabled) {
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                const NlrSubgroup& sgr_o = args.subgroups[sg];
                int n_compacted = drv.active_set_oracle_size[sg];
                int write_idx_o = 0;
                for (int k = 0; k < n_compacted; k++) {
                    int slot = drv.active_set_oracle_uvm[sg][k];
                    int i    = sgr_o.active_indices[slot];
                    AfterIterContext<Spec> ctx_o{
                        args, sg, slot, i, drv.iter_index,
                        /* h_search_current */ drv.radii_oracle_uvm[sg][slot],
                        /* scalars */         drv.cs,
                        /* scratch */         drv.scratch_oracle_uvm[sg][slot]
                    };
                    IterResult r_o = Spec::after_iter(ctx_o, drv.accum_oracle_uvm[sg][slot]);
                    drv.status_oracle_uvm[sg][slot] = static_cast<int>(r_o.status);
                    drv.new_h_oracle_uvm[sg][slot]  =
                        (r_o.status == IterStatus::AdjustRadius)
                            ? r_o.new_h_search
                            : drv.radii_oracle_uvm[sg][slot];
                    switch (r_o.status) {
                        case IterStatus::Converged: break;
                        case IterStatus::NeedsMore:
                            drv.active_set_oracle_uvm[sg][write_idx_o++] = slot;
                            break;
                        case IterStatus::AdjustRadius:
                            drv.radii_oracle_uvm[sg][slot] = r_o.new_h_search;
                            drv.active_set_oracle_uvm[sg][write_idx_o++] = slot;
                            break;
                        default:
                            if (ThisTask == 0) {
                                fprintf(stderr,
                                    "[run_neighbor_loop_iterative<%s>] FATAL: oracle "
                                    "Spec::after_iter returned unknown IterStatus=%d "
                                    "at iter=%d sg=%d slot=%d (i=%d).\n",
                                    Spec::loop_name, (int)r_o.status,
                                    drv.iter_index, sg, slot, i);
                                fflush(stderr);
                            }
                            endrun(81206);
                    }
                }
                drv.active_set_oracle_size[sg] = write_idx_o;
            }

            /* (b.compare) Step 2c.4 § 5.f: per-iter 4-thing compare.
             *   (1) AccumData via Spec::compare_accum (existing helper).
             *   (2) IterStatus (control-flow divergence).
             *   (3) new_h_search on AdjustRadius (radius-update divergence).
             *   (4) Active-set membership count + per-slot ordering. */
            int rank_now = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_now);
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                /* Codex 2c.4 non-blocking polish (2026-05-10): encode sg + iter
                 * into the origin tag so multi-subgroup iterative mismatch
                 * lines are diagnosable. Format keeps the existing "iter_"
                 * prefix the helper's [mode_b ORACLE MISMATCH ... origin=X] line
                 * uses, extended with sg+iter context. */
                char origin_tag[40];
                std::snprintf(origin_tag, sizeof(origin_tag),
                              "iter_sg%d_it%d", sg, drv.iter_index);
                /* Compare only slots that were active (evaluated) in this
                 * iteration.  The old loop over 0..n_max re-emitted stale
                 * accum_uvm entries for slots that converged in prior iters,
                 * producing identical residuals repeating across iterations. */
                for (int slot : prod_active_snapshot[sg]) {
                    /* Workaround for a Vista UVM-scatter regression observed in
                     * the dm_dispersion port (3d.D, 2026-05-16): the oracle
                     * UVM write of multi-field AccumData silently zeroed the
                     * non-leading fields (e.g. DMDispOut::Vx/Vy/Vz/VelDisp
                     * after Ngb) after scatter; pre-scatter values were
                     * correct. The dispatch helpers below now compare oracle
                     * vs production accums from host vectors before scatter
                     * and set this flag; here we skip the broken UVM-based
                     * compare to avoid false mismatches.
                     * REMOVAL CONDITION: oracle is itself temporary scaffolding
                     * (see feedback_oracle_teardown_mandatory.md); when the
                     * oracle path is removed entirely at end of the port
                     * project, this flag + guard go with it. Until then, do
                     * NOT remove without first proving the UVM-scatter
                     * regression is fixed for multi-field AccumData. */
                    if (!drv.oracle_accum_compared_in_dispatch) {
                        emit_oracle_mismatch_if_any<Spec>(rank_now, slot,
                            drv.accum_uvm[sg][slot],
                            drv.accum_oracle_uvm[sg][slot],
                            origin_tag,
                            &drv.oracle_mismatch_count);
                    }

                    if (drv.status_prod_uvm[sg][slot] != drv.status_oracle_uvm[sg][slot]) {
                        emit_oracle_status_mismatch<Spec>(rank_now, sg, slot,
                            drv.iter_index,
                            drv.status_prod_uvm[sg][slot],
                            drv.status_oracle_uvm[sg][slot],
                            &drv.oracle_mismatch_count);
                    }

                    if (drv.status_prod_uvm[sg][slot] ==
                            static_cast<int>(IterStatus::AdjustRadius) &&
                        drv.status_oracle_uvm[sg][slot] ==
                            static_cast<int>(IterStatus::AdjustRadius)) {
                        double hp = drv.new_h_prod_uvm[sg][slot];
                        double ho = drv.new_h_oracle_uvm[sg][slot];
                        double dh = std::fabs(hp - ho);
                        double tol = nlr_spec_radius_tolerance_v<Spec> * std::fabs(hp);
                        if (!(dh <= tol)) {
                            emit_oracle_radius_mismatch<Spec>(rank_now, sg, slot,
                                drv.iter_index, hp, ho,
                                &drv.oracle_mismatch_count);
                        }
                    }
                }

                int n_prod   = drv.active_set_size[sg];
                int n_oracle = drv.active_set_oracle_size[sg];
                if (n_prod != n_oracle) {
                    emit_oracle_membership_count_mismatch<Spec>(rank_now, sg,
                        drv.iter_index, n_prod, n_oracle,
                        &drv.oracle_mismatch_count);
                } else {
                    for (int k = 0; k < n_prod; k++) {
                        if (drv.active_set_uvm[sg][k] !=
                            drv.active_set_oracle_uvm[sg][k]) {
                            emit_oracle_membership_slot_mismatch<Spec>(rank_now, sg,
                                drv.iter_index, k,
                                drv.active_set_uvm[sg][k],
                                drv.active_set_oracle_uvm[sg][k],
                                &drv.oracle_mismatch_count);
                            break;
                        }
                    }
                }
            }
        }

        /* (b.5) Mode A buffer-exceedance rebuild trigger (step 2c.2).
         * After radii mutation on AdjustRadius, check if any still-active
         * slot's radius now exceeds its buffered radius (built oversized
         * at last build). If so, invalidate that subgroup's cached CSR;
         * next iter's dispatch will rebuild. Mirrors legacy density_gpu.cc:
         * 152-188 rebuild trigger. */
        if (path == DispatchPath::ModeA_GPU_NGL) {
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                if (!drv.mode_a_csr_valid[sg]) continue;
                int n_compacted = drv.active_set_size[sg];
                for (int k = 0; k < n_compacted; k++) {
                    int slot = drv.active_set_uvm[sg][k];
                    if (drv.radii_uvm[sg][slot] > drv.mode_a_csr_buffered_h[sg][slot]) {
                        drv.mode_a_csr_valid[sg] = false;
                        break;
                    }
                }
            }
        }

        /* (c) Spec::after_iter_global — host-only, no-MPI (TRAP 7). */
        if constexpr (nlr_spec_has_after_iter_global_v<Spec>) {
            Spec::after_iter_global(args, drv);
        }

        /* (d) Per-iter MPI Allreduce on PER-SUBGROUP local_active counts
         * (step 2c.3 step 5 + 6). Single array Allreduce-SUM gives per-sg
         * global activity; sum is global_active_total for the break check.
         * Globally-converged subgroups (global_active_per_sg[sg]==0) skip
         * the per-iter dispatch + collective on the next iteration. */
        if (NTask > 1) {
            MPI_Allreduce(drv.local_active_per_sg.data(),
                          drv.global_active_per_sg.data(),
                          args.num_subgroups, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        } else {
            drv.global_active_per_sg = drv.local_active_per_sg;
        }
        drv.global_active_total = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            drv.global_active_total += drv.global_active_per_sg[sg];
        }
        if (drv.global_active_total == 0) break;
    }

    /* ===== Max-iter exceeded check (matches legacy density / ags_density endrun(1155)) ===== */
    if (drv.iter_index >= Spec::max_iters && drv.global_active_total > 0) {
        if constexpr (nlr_spec_has_on_max_iter_exceeded_v<Spec>) {
            Spec::on_max_iter_exceeded(drv);
        } else {
            nlr_default_on_max_iter_exceeded<Spec>(drv);
        }
    }

    /* ===== Final apply_active_writeback (final-only per v3 §11) ===== */
    /* Fires once per active across all subgroups, on the final iteration's
     * accum (whether that iteration was Converged for that slot, or
     * max_iters terminated for everyone). The active_list semantics for
     * apply_active_writeback are the SUBGROUP'S full active_indices —
     * every active particle gets its converged result written back. */
    for (int sg = 0; sg < args.num_subgroups; sg++) {
        const NlrSubgroup& sgr = args.subgroups[sg];
        const int n_total = sgr.num_active_local;
        if (n_total <= 0) continue;
        /* Use effective_args (codex 2c.2-fix): Mode A iterative refreshed
         * P/CellP/num_total after ghost import; apply_active_writeback hooks
         * may read these for correctness. Mode B leaves effective_args ==
         * base args. */
        neighbor_loop_args sub = drv.effective_args;
        sub.active_list = sgr.active_indices;
        sub.num_active  = n_total;
        for (int slot = 0; slot < n_total; slot++) {
            int i = sgr.active_indices[slot];
            /* Codex round-10 fix 2026-05-12: if Spec opts into the
             * iterative-variant hook, route the converged radius +
             * IterScratch through it. Production-only path (oracle has
             * its own accum_oracle_uvm / radii_oracle_uvm; oracle does
             * NOT call apply_active_writeback). Specs that don't declare
             * the iterative hook fall through to the original. */
            if constexpr (nlr_spec_has_apply_active_writeback_iterative_v<Spec>) {
                Spec::apply_active_writeback_iterative(
                    sub, slot, i,
                    drv.accum_uvm[sg][slot],
                    drv.radii_uvm[sg][slot],
                    drv.scratch_uvm[sg][slot]);
            } else {
                Spec::apply_active_writeback(sub, slot, i, drv.accum_uvm[sg][slot]);
            }
        }
    }

#ifdef GIZMO_NLR_ITER_HARNESS_TEST
    /* Step 3 harness telemetry (codex step-3 review 2026-05-10): populate
     * the harness's compile-gated snapshot struct BEFORE driver destruction
     * so per-subgroup arrays are still live. Only fires for IterHarnessSpec
     * via `if constexpr` on Spec::loop_name; the strcmp lifts to a
     * compile-time string compare via SFINAE-equivalent constexpr.
     * Production Specs see this block as dead code at zero overhead. */
    if constexpr (std::is_same_v<Spec, IterHarnessSpec> ||
                  std::is_same_v<Spec, IterHarnessGhostSpec>) {
        g_iter_harness_telemetry.last_iter_index          = drv.iter_index;
        g_iter_harness_telemetry.last_global_active_total = drv.global_active_total;
        g_iter_harness_telemetry.csr_rebuild_count        = drv.csr_rebuild_count;
        g_iter_harness_telemetry.arena_acquire_count      = drv.arena_acquire_count;
        g_iter_harness_telemetry.csr_local_rebuild_count  = drv.csr_local_rebuild_count;
        g_iter_harness_telemetry.final_dummy_jflag_prod   = drv.ctx.dummy_jflag;
        g_iter_harness_telemetry.final_dummy_jflag_oracle =
            drv.oracle_enabled ? drv.ctx_oracle.dummy_jflag : 0;
        g_iter_harness_telemetry.oracle_mismatch_count    = drv.oracle_mismatch_count;
        g_iter_harness_telemetry.oracle_enabled           = drv.oracle_enabled ? 1 : 0;
        /* Sum pair counts across all converged slots in all subgroups.
         * accum_uvm holds the final-iter accum for each slot (never zeroed
         * after convergence — codex v4.4 invariant). For non-oracle runs,
         * accum_oracle_uvm pointers are nullptr; skip the sum. */
        long long pairs_prod = 0, pairs_oracle = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            int n_local = args.subgroups[sg].num_active_local;
            if (drv.accum_uvm[sg]) {
                for (int s = 0; s < n_local; s++) {
                    pairs_prod += drv.accum_uvm[sg][s].n_pairs;
                }
            }
            if (drv.oracle_enabled && drv.accum_oracle_uvm[sg]) {
                for (int s = 0; s < n_local; s++) {
                    pairs_oracle += drv.accum_oracle_uvm[sg][s].n_pairs;
                }
            }
        }
        g_iter_harness_telemetry.total_pairs_prod   = pairs_prod;
        g_iter_harness_telemetry.total_pairs_oracle = pairs_oracle;
        /* Per-slot scratch + status — sg=0, up to 4 slots. */
        for (int s = 0; s < 4; s++) {
            g_iter_harness_telemetry.scratch_passes[s]    = 0;
            g_iter_harness_telemetry.final_status_prod[s] = -1;
        }
        if (args.num_subgroups > 0 && args.subgroups[0].num_active_local > 0) {
            const int n_max = std::min(args.subgroups[0].num_active_local, 4);
            for (int s = 0; s < n_max; s++) {
                g_iter_harness_telemetry.scratch_passes[s] =
                    drv.scratch_uvm[0][s].passes_so_far;
                if (drv.oracle_enabled && drv.status_prod_uvm[0] != nullptr) {
                    g_iter_harness_telemetry.final_status_prod[s] =
                        drv.status_prod_uvm[0][s];
                }
            }
        }
    }
#endif

    }  /* end inner scope: driver destructs HERE (codex 2c.2-fix-v2 hardening 1) */

    /* ===== Mode B hard-corridor enforcement (codex 2c.2-fix blocker #3 +
     * v2 hardening 1) =====
     * Check AFTER final apply_active_writeback AND after driver destruction
     * so the invariant covers the FULL Mode B path including writeback hooks
     * and driver cleanup. Mode A paths legitimately advanced counters; skip
     * enforcement. */
    if (path == DispatchPath::ModeB_HostWalker) {
        const bool drift_violation = (g_global_drift_counter      != s_drift0);
        const bool ghost_violation = (g_ghost_import_counter      != s_ghost0);
        const bool arena_violation = (g_gpu_arena_acquire_counter != s_arena0);
        const bool np_violation    = (NumPart != s_np0);
        if (drift_violation || ghost_violation || arena_violation || np_violation) {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                "[NLR_ITER CORRIDOR ABORT rank=%d caller=%s path=mode_b] Mode B "
                "iterative path violated tiny-N corridor invariant during "
                "run_neighbor_loop_iterative. Counter deltas: drift=%llu "
                "ghost=%llu arena=%llu NumPart_pre=%d NumPart_post=%d\n",
                rank, Spec::loop_name,
                (unsigned long long)(g_global_drift_counter      - s_drift0),
                (unsigned long long)(g_ghost_import_counter      - s_ghost0),
                (unsigned long long)(g_gpu_arena_acquire_counter - s_arena0),
                s_np0, NumPart);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 81213);
        }
    }
    /* Driver destructor frees per-subgroup UVM allocations on scope exit. */
}

/* ============================================================================
 * Explicit template instantiations — one per migrated caller.
 *
 * Forgetting an instantiation = clean linker error at the call site. Each
 * Spec must opt in here and is implicitly acknowledging its declared
 * Spec::sidx_cache_kind; cf. nlr_resolve_sidx_cache.
 *
 * No `run_neighbor_loop_iterative<...>` instantiations in step 2a — there
 * are no callers yet (synthetic harness lands in step 3, ags_density in 3d.4).
 * ========================================================================== */

#ifdef SINK_PARTICLES
template void run_neighbor_loop<SinkEnv1Spec>(const neighbor_loop_args&);
template void run_neighbor_loop<SinkFeedSpec>(const neighbor_loop_args&);
template void run_neighbor_loop<SinkSwkSpec>(const neighbor_loop_args&);
#if defined(SINK_GRAVACCRETION) && (SINK_GRAVACCRETION == 0)
template void run_neighbor_loop<SinkEnv2Spec>(const neighbor_loop_args&);
#endif
#endif

#ifdef GIZMO_NLR_ITER_HARNESS_TEST
/* Harness header is included near top of file (above run_neighbor_loop_iterative
 * body) so the telemetry struct + IterHarnessSpec name are visible inside
 * the template body's `if constexpr (std::is_same_v<Spec, IterHarnessSpec>)`
 * gates. Only the explicit instantiation lives here. */
template void run_neighbor_loop_iterative<IterHarnessSpec>(const neighbor_loop_args_iterative&);
template void run_neighbor_loop_iterative<IterHarnessGhostSpec>(const neighbor_loop_args_iterative&);
#endif

/* 3d.4: AgsDensitySpec — first production iterative Spec consumer of the
 * 4.B.0 runner + step-0 partition-by-subgroup + sticky-call-scope wakeup
 * traits. See gravity/ags_density_loop.{h,cc}. */
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
template void run_neighbor_loop_iterative<AgsDensitySpec>(const neighbor_loop_args_iterative&);
#endif

/* Phase 4 Wave-1: DensitySpec — hydro density runner port (Step 6).
 * Single gas-only subgroup, oracle-safe (no after_iter P/CellP writes),
 * uses apply_active_writeback_iterative for oracle-safe radius
 * channeling (codex round-10 fix). See hydro/density_loop.{h,cc}. */
template void run_neighbor_loop_iterative<DensitySpec>(const neighbor_loop_args_iterative&);

/* Phase 4 Wave-3 / 3e.1: MechFBSpec — mechanical-feedback runner port
 * (milestone 3: physics-complete pair kernel + full state-machine Spec
 * contract). 6-mode iterative state machine over loop_iteration
 * {-2,-1,0,1,2,3}; mode_a_rebuild_csr_every_iter=false preserves the legacy
 * 1-CSR-shared-across-6-modes optimization. Mode A multi-rank ghost-side
 * writes are guarded by Kokkos::abort pending milestone 3.5 (custom
 * MechFBGasDelta ghost-writeback callback + lazy d_gas_iter). See
 * galaxy_sf/mechfb_loop.{h,cc} + OPEN_3d_mechfb_design.md. */
#ifdef GALSF_FB_MECHANICAL
template void run_neighbor_loop_iterative<MechFBSpec>(const neighbor_loop_args_iterative&);
#endif

/* Phase 4 Wave-3 / 3e.2: ThermalFBSpec — thermal-feedback runner port.
 * Non-iterative scatter (Type-4 stars → gas neighbors); ActiveReduceOnly +
 * manifest-bundle ghost_writeback; sink_feed pattern. See
 * galaxy_sf/thermal_fb_loop.{h,cc}. The Spec definition is gated on
 * GALSF_FB_THERMAL in thermal_fb_loop.h, so the explicit instantiation must
 * sit inside the same #ifdef (Wave 2 lesson — non-thermal Configs would
 * otherwise hit an undefined type at this template instantiation site). */
#ifdef GALSF_FB_THERMAL
template void run_neighbor_loop<ThermalFBSpec>(const neighbor_loop_args&);
#endif

/* Phase 4 Wave-3 / radfb_local: RadFBRPSpec — local radiation-pressure winds.
 * Iterative 2-pass (iter 0 wt_sum aggregation; iter 1 kick application).
 * Ghost-writeback bundle with PARTICLE_ADD_VEC3 + new GAS_ADD_VEC3 ops.
 * iter-gating via Aux::iter_index (set in reset_per_iter_device_context);
 * inter-iter wt_sum staged through IterScratch by after_iter_global. See
 * galaxy_sf/radfb_rp_loop.{h,cc} + OPEN_3d_radfb_local_design.md. */
#ifdef GALSF_FB_FIRE_RT_LOCALRP
template void run_neighbor_loop_iterative<RadFBRPSpec>(const neighbor_loop_args_iterative&);
#endif

/* Phase 4 Wave-4 / 3d.D: DMDispersionSpec — DM velocity dispersion runner port.
 * Gas actives (Type==0) iterate bisection on KernelRadiusDM to enclose
 * 64±48 DM (Type==1) neighbors; accumulates unweighted Vel sums for dispersion.
 * ActiveReduceOnly + SidxCacheKind::None (DM tbm matches neither GasOnly nor
 * AllTypes cache). apply_active_writeback_iterative + Aux finalize pattern
 * mirrors DensitySpec exactly. See galaxy_sf/dm_dispersion_loop.{h,cc}. */
#if defined(GALSF_SUBGRID_WINDS) && (GALSF_SUBGRID_WIND_SCALING==2)
template void run_neighbor_loop_iterative<DMDispersionSpec>(const neighbor_loop_args_iterative&);
#endif

/* Wave 5 / difffilter: DiffFilterSpec + DynDiffSpec — TURB_DIFF_DYNAMIC
 * velocity-smoothing + dynamic-Smagorinsky loops. Non-iterative,
 * scaled-symmetric gas-gas (symmetric_neighbor_radius_scale = TurbDynamicDiffFac),
 * pure i-side reduce. See turb/difffilter_loop.{h,cc} +
 * OPEN_3d_difffilter_design.md. */
#ifdef TURB_DIFF_DYNAMIC
template void run_neighbor_loop<DiffFilterSpec>(const neighbor_loop_args&);
template void run_neighbor_loop<DynDiffSpec>(const neighbor_loop_args&);
#endif

