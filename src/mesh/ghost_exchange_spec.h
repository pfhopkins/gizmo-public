/* ghost_exchange_spec.h — public spec struct + GHOST_TYPE_* / NGB_SEARCH_*
 * macros for callers that build their own spec literal in-line and call
 * ghost_exchange_run() (declared in core/proto.h).
 *
 * The spec literal at the caller is the single source of truth for that
 * loop's (supply_type_mask, search_mode, query list). Editing those there
 * is the only edit needed to flip the call's physics. Dispatch keys only on
 * spec fields: an explicit query list (n_queries >= 0), search_mode ==
 * NGB_SEARCH_ONEWAY, or a SYMMETRIC spec with supply_band_dominated set and
 * safety_factor <= 1 selects request-driven; every other SYMMETRIC caller
 * uses tile-overlap/broadcast. No per-caller special-casing lives in
 * ghost_exchange.cc.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */

#ifndef GHOST_EXCHANGE_SPEC_H
#define GHOST_EXCHANGE_SPEC_H

#include "neighbor_list.h"   /* NGB_SEARCH_ONEWAY, NGB_SEARCH_SYMMETRIC */
#include "nlr_radius_policy.h"  /* mode_b_radius_policy_t + MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES */

/* Type bitmask values matching P[].Type encoding. Bit k set ⇔ Type k included.
 * Keep these names literal: Type 0 is usually cells/fluids, Type 5 is often
 * sinks, but Types 1-4 are intentionally physics/config dependent. */
#define GHOST_TYPE_0     (1u << 0)
#define GHOST_TYPE_1     (1u << 1)
#define GHOST_TYPE_2     (1u << 2)
#define GHOST_TYPE_3     (1u << 3)
#define GHOST_TYPE_4     (1u << 4)
#define GHOST_TYPE_5     (1u << 5)
#define GHOST_TYPE_ALL   ((1u << 6) - 1u)

struct ghost_exchange_spec_t {
    /* Legacy: filters ActiveParticleList in tile-overlap impl AND in the
     * request-driven impl when n_queries < 0. New explicit-query callers
     * leave 0. */
    unsigned int request_type_mask;

    /* Bitmask of Types eligible to be SUPPLIED as ghosts (the search pool). */
    unsigned int supply_type_mask;

    /* NGB_SEARCH_ONEWAY (r_ij < h_i) or NGB_SEARCH_SYMMETRIC (r_ij < max(h_i,h_j)). */
    int          search_mode;

    /* Multiplier on per-query h before the predicate. */
    double       safety_factor;

    /* Caller label for diagnostic prints only — dispatch never keys on it. */
    const char  *caller_name;

    /* Caller-owned explicit query list. n_queries >= 0 → use these. n_queries
     * < 0 (sentinel) → use legacy ActiveParticleList scan filtered by
     * request_type_mask. query_h is the RAW h, NOT pre-multiplied by
     * safety_factor (the impl applies it). */
    int                  n_queries;
    const double       (*query_pos)[3];
    const double        *query_h;

    /* SSOT supply-side per-particle reach contract. Mirrors what local
     * Mode A's gpu_ngb_list_build uses for the same Spec
     * so the local CSR and the ghost-import candidate set are computed from
     * the IDENTICAL per-particle h_j formula on every rank.  Without this,
     * runner-driven imports under non-default Specs would silently disagree
     * with the local walk → rank-dependent neighbor selection → invalid
     * multi-rank physics.
     *
     * Spec::radius_policy is a compile-time constexpr; identical on every
     * rank.  j_radius_scale comes from nlr_spec_symmetric_j_radius_scale<Spec>()
     * — also a per-Spec/global value, identical across ranks at the call.
     *
     * Legacy non-runner callers MUST pass MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES
     * + 1.0 explicitly; this preserves their pre-policy behavior byte-for-byte
     * (raw P[j].KernelRadius * safety_factor).  The struct provides no
     * defaults so the compiler enforces explicit-thread at every call site
     * — no LEGACY fallback is allowed in runner-driven paths. */
    mode_b_radius_policy_t  radius_policy;
    double                  j_radius_scale;

    /* Supply-band domination: 1 iff this spec's supply-side reach — under its
     * own (radius_policy, j_radius_scale, safety_factor) — is provably bounded
     * by the per-type node band the sender opener walks against. Only then is
     * routed (request-driven) SYMMETRIC discovery COMPLETE; an unproven reach
     * could exceed the band and silently under-import.
     *
     * 0 is the FAIL-CLOSED default in every sense that matters: a spec left at
     * 0 keeps the broadcast path, i.e. exactly today's behavior. Never set 1
     * without citing the bound at the call site. This is a STRUCTURAL property
     * of the spec, never a caller identity — dispatch must not name callers. */
    int                     supply_band_dominated;
};

#endif /* GHOST_EXCHANGE_SPEC_H */
