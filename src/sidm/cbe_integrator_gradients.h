/*
 * sidm/cbe_integrator_gradients.h
 *
 * CBE pre-force gradient module. Uses the persistent-particle-field shape
 * shared with hydro and DMGrad. Persistent storage holds the spatial
 * gradients of the reconstructed flux-frame basis content, per basis, per
 * slot, per direction:
 *   double P[i].Gradients_CBE_basis_moments[NBASIS][NMOMENTS][3]
 * (declarations/particle_data.h, gated on CBE_INTEGRATOR_WITHGRADIENTS).
 *
 * PRIMITIVE-gradient swap:
 * the field name `Gradients_CBE_basis_moments` is now STALE — the slots
 * store gradients of PRIMITIVES (∂ρ, ∂v_k, ∂S_kl), not gradients of the
 * moment row (m, p_k, T_kl). Slot indexing is identical to the moment
 * row (slot 0 = ρ; slots 1..NUMDIMS = v_k; stress slots = S_kl via the
 * existing cbe_T_idx packing). The field-name rename is out of scope
 * (touches scatter / IO / snapshot); the stale name is documented in
 * declarations/particle_data.h. Pass-0 LSQ accumulates primitive deltas,
 * pass-1 limiter is COMPONENT-SEPARATED: ρ gets
 * BJ + positivity (affects only φ_ρ); v_k gets scale-aware BJ only (no
 * realizability role; affects only φ_v_k); the S block gets per-slot BJ
 * + tensor PSD (single φ_S across all stress slots to preserve the PSD
 * cone). The flux body reconstructs primitive face state then converts
 * back to a moment row via cbe_primitives_to_moments_row before the
 * HLLC vacuum solver consumes it.
 *
 * Standard GIZMO ghost import (gizmo_request_filtered_ghost_import_fresh)
 * carries the gradient field naturally as part of P[]. There is no scratch
 * UVM buffer, no custom Alltoallv, and no raw gradient pointer in any
 * DeviceContext.
 *
 * Toplevel `CBEGrad_gradient_calc()` runs pre-force (called from
 * core/accel.cc, paralleling DMGrad_gradient_calc) and refreshes
 * P[i].Gradients_CBE_basis_moments for every currently-active particle —
 * active set is AGSForce_isactive(i), same predicate as the force consumer
 * AGSForce_calc. Inactive particles retain their previous-step gradient
 * (hydro semantics). Subgroup partitioning by ags_gravity_kernel_shared_BITFLAG
 * is done at the toplevel via per-bm `args.neighbor_type_mask_override`
 * calls — same pattern as DMGrad.
 *
 * Two passes orchestrated at the toplevel via Aux::loop_iteration:
 *   pass 0  — raw MFM-LSQ gradient (M += w*dp*dp^T, B += w*dQ*(x_j - x_i));
 *             apply_active_writeback solves M^{-1}.B with Tikhonov-style
 *             ill-conditioning guard and writes
 *             P[i].Gradients_CBE_basis_moments.
 *   pass 1  — component-separated primitive limiter (per-(m,k) phi from
 *             cbe_bj_phi_pair_scaled + per-component realizability —
 *             positivity for ρ, none for v, PSD for the S block; min over
 *             matched neighbors); writeback rescales
 *             P[i].Gradients_CBE_basis_moments in place by phi.
 *
 * "Pairwise BJ-style conservative limiter": a stricter per-pair matched-
 * basis local-extremum form, NOT the full global-stencil Barth-Jespersen
 * (which takes Q_min/Q_max over the whole neighbor stencil). The name is
 * load-bearing — keep the qualifier.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef CBE_INTEGRATOR_GRADIENTS_H
#define CBE_INTEGRATOR_GRADIENTS_H

/* Kokkos_Core MUST come BEFORE declarations/allvars.h — same convention
 * as ags_force_loop.h / dm_fuzzy_loop.h. */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"

#if defined(CBE_INTEGRATOR_WITHGRADIENTS)

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"
#include "../gravity/ags_functions.h"       /* get_particle_volume_ags_P */
#include "cbe_integrator_functions.h"       /* cbe_build_flux_frame_Q_from_stored_moments,
                                               cbe_cost_v_only, cbe_assign_outgoing_greedy */
/* NOTE: caller TUs MUST include "../mesh/kernel.h" before this header
 * (kernel.h lacks include guards; same convention as ags_force_loop.h). */

int AGSForce_isactive(int i);   /* forward decl; defined in gravity/ags_force_loop.cc */

/* Toplevel — defined in sidm/cbe_integrator_gradients.cc; also declared
 * (guarded) in core/proto.h. */
void CBEGrad_gradient_calc(void);

/* ============================================================================
 * Spec types.
 * ========================================================================== */

/* Per-call cosmology + the external gradient-pass index. populate_call_scalars
 * copies loop_iteration out of Aux; pair_kernel + apply_active_writeback
 * branch on it. Same shape as DMGradCallScalars. */
struct CBEGradCallScalars {
    NlrCommonScalars common;
    int              loop_iteration;     /* 0 = LSQ; 1 = pairwise BJ-style limiter */
};

/* Per-active local fill. Carries x_i, h_i, V_i, the stored moment array U_i,
 * Vel_i for Q_i construction, AND a by-value snapshot of the persistent
 * gradient row P[i].Gradients_CBE_basis_moments. The gradient snapshot is
 * Mode-B-safe envelope-transit: pass-1 pair_kernel
 * reads active.local.Gradients_CBE_basis_moments instead of dereferencing
 * peer-rank storage. Populated unconditionally by load_active — pass 0
 * (LSQ) does not consume it, but the cost is fixed-size (NBASIS*NMOMENTS*3
 * doubles, ≤ 960 B worst case) and the unified shape keeps the load_active
 * branch-free. */
struct CBEGradLocalIn {
    double         Mass;
    double         AGS_KernelRadius;
    Vec3<double>   Pos;
    Vec3<double>   Vel;
    double         V_i;
    double         U_i[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    /* By-value snapshot of P[i].Gradients_CBE_basis_moments — written by
     * the previous CBEGrad_gradient_calc pass (or zero-init if first call).
     * Pass 1 limiter reads this; pass 0 LSQ does not consume it. */
    double         Gradients_CBE_basis_moments[CBE_INTEGRATOR_NBASIS]
                                              [CBE_INTEGRATOR_NMOMENTS][3];
    /* Prev-step (stale) gradient snapshot, for the density-continuity matching cost. */
    double         Gradients_CBE_basis_moments_prev[CBE_INTEGRATOR_NBASIS]
                                                   [CBE_INTEGRATOR_NMOMENTS][3];
};

/* Per-active accumulator. Both passes share one POD:
 *   pass 0 (LSQ):  M (3x3, sum-merge) + B[m][k][d] (sum-merge); phi stays 1.0
 *   pass 1 (row-scalar BJ + cone): phi[m][k] (min-merge, NaN-aware);
 *          M, B stay 0
 * merge_accum unconditionally sums M+B AND mins phi — the unused-pass
 * fields are neutral elements (0 for sum, 1 for min) so the per-pass
 * branch isn't needed in merge_accum.
 *
 * Pass 1 is ROW-SCALAR: a single phi per basis row m is produced per neighbor,
 * and the same value is written to every k in row m. The per-(m,k)
 * storage layout is preserved for ABI stability; after min-merge over
 * neighbors all k in a given row therefore carry the same value, and
 * the per-active writeback rescales the whole row coupled (m, p, T)
 * by that scalar. Mixing different phi across k in a row would move
 * the reconstructed state off the line segment the cone predicate
 * validated and could reintroduce non-realizability.
 *
 * MUST be a POD; zero_accum sets M=0, B=0, phi=1.0 explicitly. */
struct CBEGradOut {
    double M[3][3];
    double B[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS][3];
    double phi[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
};

struct CBEGradActiveState {
    Vec3<double>       pos;
    double             h_search;
    CBEGradLocalIn     local;
    short int          TimeBin;
    int                origin_local_idx;
    int                origin_rank;
};

/* ============================================================================
 * Inline pair-kernel helpers (KOKKOS_INLINE_FUNCTION; instantiated from the
 * GPU TU that holds the runner instantiation).
 * ========================================================================== */

/* Pass-0 pair body — MFM-LSQ M and B accumulation. */
template <typename NeighborT>
KOKKOS_INLINE_FUNCTION
static void cbe_grad_lsq_pair_kernel_body(const CBEGradActiveState& active,
                                           const CBEGradCallScalars& scalars,
                                           const NeighborT&          neighbor,
                                           CBEGradOut&               accum)
{
    const int             j   = neighbor.j;
    struct particle_data &Pj  = *neighbor.neighbor_particle;
    const CBEGradLocalIn& L   = active.local;

    if(!(Pj.Mass > 0) || !(Pj.AGS_KernelRadius > 0)) return;

    /* dp = x_i - x_j. The B accumulation uses (-dp) = (x_j - x_i) so
     * apply_active_writeback solves M . grad = sum_j w * dQ * (x_j - x_i),
     * the design-doc convention. M is even in the sign of dp. */
    Vec3<double> dp;
    dp[0] = L.Pos[0] - Pj.Pos[0];
    dp[1] = L.Pos[1] - Pj.Pos[1];
    dp[2] = L.Pos[2] - Pj.Pos[2];
    NEAREST_XYZ(dp[0], dp[1], dp[2], -1);

    const double r2  = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
    const double h_i = L.AGS_KernelRadius;
    if(r2 <= 0) return;
    /* h_i-only acceptance, DESPITE the Spec's symmetric search topology.
     * Matches the raw-gradient semantics of hydro / density / DMGrad:
     * linear consistency of the LSQ estimator relies on the i-kernel
     * weighting + NV_T-style normalisation. The Spec uses symmetric search
     * because pass 1 (pairwise BJ-style limiter) needs the broader face
     * visibility; pass 0 deliberately re-narrows here. Do NOT widen this
     * to a symmetric r <= h_i+h_j (or r <= max(h_i,h_j)) filter. */
    if(r2 >= h_i * h_i) return;

    const double r = sqrt(r2);
    double hinv_i, hinv3_i, hinv4_i;
    kernel_hinv(h_i, &hinv_i, &hinv3_i, &hinv4_i);

    const double u_i = r * hinv_i;
    double       wk_i, dwk_i;
    kernel_main(u_i, hinv3_i, hinv4_i, &wk_i, &dwk_i, 0);

    /* Volume per neighbor — canonical AGS volume helper. Same call the
     * flux body uses (sidm/cbe_integrator_flux_functions.h). */
    const double V_j = get_particle_volume_ags_P(j, neighbor.P_base);

    const double w = wk_i * V_j;

    /* M_i += w * dp dp^T (outer product even in sign of dp). */
    for(int a = 0; a < 3; a++) {
        for(int b = 0; b < 3; b++) {
            accum.M[a][b] += w * dp[a] * dp[b];
        }
    }

    /* Build Q_i, Q_j with the SSOT helper used by the flux body. */
    const double cf_a3inv = scalars.common.cf_a3inv;
    const double cf_atime = scalars.common.cf_atime;

    double Q_i[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    double Q_j[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    {
        double Vel_i_code[3] = { L.Vel[0], L.Vel[1], L.Vel[2] };
        cbe_build_flux_frame_Q_from_stored_moments(
            L.U_i, Vel_i_code, L.V_i, cf_a3inv, cf_atime, Q_i);
    }
    {
        double Vel_j_code[3] = { (double)Pj.CBE_VelPred[0], (double)Pj.CBE_VelPred[1], (double)Pj.CBE_VelPred[2] };
        cbe_build_flux_frame_Q_from_stored_moments(
            Pj.CBE_basis_moments_pred, Vel_j_code, V_j, cf_a3inv, cf_atime, Q_j);
    }

    /* Basis matching via SSOT helper. Matching API
     * stays on Q (moment rows); the cost helpers cbe_cost_v_only and
     * cbe_cost_trace_w2 derive primitive content internally (v = p/ρ;
     * tr S = T/ρ − v²). Therefore the primitive-gradient swap does not
     * require any matching-API change — only the LSQ delta computed
     * below is now primitive. Same SSOT helper / selectors as the
     * limiter pass and the flux body. Q-cell vs Q-face matching split
     * is unchanged. Single-direction; fired-count NULL. */
    int matched_j_for_i[CBE_INTEGRATOR_NBASIS];
    /* Stale-gradient reconstructed face densities for the density-continuity
     * matching cost. Gradient pass: vF=0 (no face velocity in the stencil).
     * NULL when gradients are not compiled in -> base velocity/trace cost. */
    const double *gcrho_a = (const double*)0, *gcrho_b = (const double*)0;
#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
    double gcrho_fi[CBE_INTEGRATOR_NBASIS], gcrho_fj[CBE_INTEGRATOR_NBASIS];
    {
        double h_j_c = (double)Pj.AGS_KernelRadius; double pden = h_i + h_j_c;
        double psi_i = (pden > 0) ? h_j_c/pden : 0.5; double psi_j = 1.0 - psi_i;
        for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
            double gi = L.Gradients_CBE_basis_moments_prev[m][0][0]*dp[0]
                      + L.Gradients_CBE_basis_moments_prev[m][0][1]*dp[1]
                      + L.Gradients_CBE_basis_moments_prev[m][0][2]*dp[2];
            double gj = Pj.Gradients_CBE_basis_moments_prev[m][0][0]*dp[0]
                      + Pj.Gradients_CBE_basis_moments_prev[m][0][1]*dp[1]
                      + Pj.Gradients_CBE_basis_moments_prev[m][0][2]*dp[2];
            gcrho_fi[m] = Q_i[m][0] - psi_i*gi;
            gcrho_fj[m] = Q_j[m][0] + psi_j*gj;
        }
        gcrho_a = gcrho_fi; gcrho_b = gcrho_fj;
    }
#endif
    const double vface_zero[3] = {0.0, 0.0, 0.0};   /* gradient stencil: no single face normal, no face shift */
    cbe_build_pair_matching(Q_i, Q_j, matched_j_for_i,
                            /*alpha_of_beta_for_b=*/NULL,
                            /*free_slot_fired_count_inout=*/NULL,
                            gcrho_a, gcrho_b, vface_zero);

    /* Primitive-gradient swap: convert Q_i, Q_j to primitive rows once per
     * basis. The LSQ delta below is on PRIMITIVE slots — the persistent
     * P[i].Gradients_CBE_basis_moments field now stores primitive gradients
     * (∂ρ, ∂v, ∂S) packed in the slot layout shared with the moment row
     * (slot 0 = density, slots 1..NUMDIMS = v_k, stress slots = S_kl).
     * See cbe_moments_to_primitives_row. */
    double prim_i[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    double prim_j[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
        cbe_moments_to_primitives_row(Q_i[m], prim_i[m]);
        cbe_moments_to_primitives_row(Q_j[m], prim_j[m]);
    }

    /* B_i[m][k][e] += w * (prim_j_matched − prim_i)[m][k] * (x_j − x_i)[e].
     * With dp = x_i − x_j, (x_j − x_i)[e] = −dp[e].
     *
     * Scratch-row skip: if cell i
     * basis m has ρ_i,m ≤ cbe_rho_active_floor(), leave B[m][:][:] at zero
     * for this pair. The pass-0 writeback then naturally produces a zero
     * gradient row (M_inv·0 = 0) — first-order reconstruction for that
     * basis at the cell. The cell's other (active) bases get second-order
     * primitive reconstruction normally. */
    for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
        if(!(Q_i[m][0] > cbe_rho_active_floor())) continue;
        const int n = matched_j_for_i[m];
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
            const double dprim = prim_j[n][k] - prim_i[m][k];
            const double w_dprim = w * dprim;
            accum.B[m][k][0] += w_dprim * (-dp[0]);
            accum.B[m][k][1] += w_dprim * (-dp[1]);
            accum.B[m][k][2] += w_dprim * (-dp[2]);
        }
    }
    (void)j;
}

/* ----------------------------------------------------------------------------
 * Scale-aware pairwise BJ-style phi for the primitive-gradient limiter.
 * Bounds the reconstructed delta
 * `predicted` so that prim_face = prim_i + predicted stays between prim_i
 * and prim_j_matched.
 *
 * Why scale: an earlier helper used tol = |dQ| * 1e-14 +
 * MIN_REAL_NUMBER. For primitive slots where dprim ≡ 0 exactly (uniform
 * v in cold streams, uniform S in cold/warm-small IC) but predicted
 * carries LSQ roundoff (~1e-15..1e-17), the test fell through to
 * ratio = 0 / predicted = 0 and vetoed the row. The scale-aware tol uses
 * the physical magnitude of the primitive quantity (or its R-row source
 * for stress, since S = R - vv has FP-cancellation noise at the R scale):
 *
 *   ρ slot       :  scale = max(|ρ_i|, |ρ_j|, ρ_floor)
 *   v_k slot     :  scale = max(|v_i,k|, |v_j,k|, sqrt(max(R_kk_i, R_kk_j, 0)))
 *   S_kl slot    :  scale = max(|S_kl_i|, |S_kl_j|, |R_kl_i|, |R_kl_j|)
 *                   with R_kl = T_kl / ρ
 *
 * Logic:
 *   - non-finite anywhere               → 0.0
 *   - |predicted| ≤ tol(scale)          → 1.0  (no meaningful predicted; no constraint)
 *   - |dQ| ≤ tol(scale)                 → 0.0  (predicted heads away from equal neighbor)
 *   - sign(predicted) ≠ sign(dQ)        → 0.0
 *   - same sign, |pred| > |dQ|          → dQ/predicted  (∈ (0, 1))
 *   - same sign, |pred| ≤ |dQ|          → 1.0
 *
 * Caller takes min across (j, m); the limiter body is COMPONENT-SEPARATED:
 * one phi per slot per basis, NOT a single row-scalar.
 * ---------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static double cbe_bj_phi_pair_scaled(double dQ, double predicted, double scale,
                                     double oppsign_tol = 0.0)
{
    if(!isfinite(dQ) || !isfinite(predicted) || !isfinite(scale)) return 0.0;
    const double tol = fabs(scale) * 1e-14 + MIN_REAL_NUMBER;
    /* Quiet-slot: predicted is below scale-relative roundoff → no constraint. */
    if(fabs(predicted) <= tol) return 1.0;
    /* Real-signal predicted but neighbor doesn't differ at this scale → overshoot → 0. */
    if(fabs(dQ) <= tol) return 0.0;
    const double ratio = dQ / predicted;
    if(!isfinite(ratio)) return 0.0;
    if(ratio < 0.0) {
        /* Finite opposite-sign tolerance — a sign flip whose reconstructed jump
         * |predicted| is small vs the local scale is safe to ignore (velocity:
         * modest overshoot is physical). oppsign_tol=0 → hard clamp; ρ/S callers
         * pass 0 so are unchanged, the v caller passes CBE_VOPPSIGN. */
        if(oppsign_tol > 0.0 && fabs(predicted) < oppsign_tol * fabs(scale)) return 1.0;
        return 0.0;
    }
    if(ratio > 1.0) return 1.0;
    return ratio;
}


/* ----------------------------------------------------------------------------
 * S-block PSD trial-φ predicate and bisection (3D stress block). The
 * primitive-gradient limiter applies a single φ_S to the entire stress
 * block to preserve PSD as a tensor (mixing per-slot φ on off-diagonals
 * would move S_face off the PSD cone). 1D NMOMENTS=3 (scalar S_xx) gets
 * a closed-form positivity bound in the pair body without these helpers.
 * ---------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static bool cbe_S_psd_at_phi(const double S_i[3][3],
                              const double S_delta[3][3],
                              double phi)
{
    double S_face[3][3];
    for(int a = 0; a < 3; a++)
        for(int b = 0; b < 3; b++)
            S_face[a][b] = S_i[a][b] + phi * S_delta[a][b];
    return cbe_symNxN_active_principal_minors_nonneg(S_face, /*eps_tol=*/ 0.0);
}

KOKKOS_INLINE_FUNCTION
static double cbe_S_psd_phi_bisect(const double S_i[3][3],
                                    const double S_delta[3][3],
                                    double phi_upper)
{
    if(!(phi_upper > 0)) return 0;
    if(!cbe_S_psd_at_phi(S_i, S_delta, 0.0)) return 0;
    if(cbe_S_psd_at_phi(S_i, S_delta, phi_upper)) return phi_upper;
    double lo = 0;
    double hi = phi_upper;
    for(int it = 0; it < 60; it++) {
        const double scale = DMAX(fabs(lo), fabs(hi));
        if((hi - lo) <= 1e-14 + 1e-14 * scale) break;
        const double mid = 0.5 * (lo + hi);
        if(cbe_S_psd_at_phi(S_i, S_delta, mid)) lo = mid;
        else                                    hi = mid;
    }
    return lo;
}


/* Pass-1 pair body — pairwise BJ-style conservative limiter. Reads
 * active.local.Gradients_CBE_basis_moments (by-value snapshot from
 * load_active — Mode-B-safe, NOT a peer-rank dereference). */
template <typename NeighborT>
KOKKOS_INLINE_FUNCTION
static void cbe_grad_bj_pair_kernel_body(const CBEGradActiveState& active,
                                          const CBEGradCallScalars& scalars,
                                          const NeighborT&          neighbor,
                                          CBEGradOut&               accum)
{
    struct particle_data &Pj  = *neighbor.neighbor_particle;
    const CBEGradLocalIn& L   = active.local;
    const int             j   = neighbor.j;

    if(!(Pj.Mass > 0) || !(Pj.AGS_KernelRadius > 0)) return;

    /* dp = x_i - x_j. Mirrors the LSQ pass and the flux body. */
    Vec3<double> dp;
    dp[0] = L.Pos[0] - Pj.Pos[0];
    dp[1] = L.Pos[1] - Pj.Pos[1];
    dp[2] = L.Pos[2] - Pj.Pos[2];
    NEAREST_XYZ(dp[0], dp[1], dp[2], -1);

    const double r2  = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
    const double h_i = L.AGS_KernelRadius;
    if(r2 <= 0) return;
    const double r   = sqrt(r2);
    const double h_j = (double)Pj.AGS_KernelRadius;

    /* Pair filter mirrors AgsForceSpec so the limiter sees every face the flux
     * body will subsequently reconstruct: the ordinary symmetric support. */
    if(r > h_i && r > h_j) return;

    /* MFM face-position weight (matches the flux body's reconstruction:
     * psi_i = h_j/(h_i+h_j); face offset on i's side = -psi_i * dp). */
    const double psi_i_denom = h_i + h_j;
    if(!(psi_i_denom > 0)) return;
    const double psi_i = h_j / psi_i_denom;

    /* Volume per neighbor — same helper as gradient and flux bodies. */
    const double V_j = get_particle_volume_ags_P(j, neighbor.P_base);

    /* Build Q_i and Q_j with the SSOT helper — same matching keeps the
     * limiter consistent with the LSQ pass and the flux body. */
    const double cf_a3inv = scalars.common.cf_a3inv;
    const double cf_atime = scalars.common.cf_atime;

    double Q_i[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    double Q_j[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    {
        double Vel_i_code[3] = { L.Vel[0], L.Vel[1], L.Vel[2] };
        cbe_build_flux_frame_Q_from_stored_moments(
            L.U_i, Vel_i_code, L.V_i, cf_a3inv, cf_atime, Q_i);
    }
    {
        double Vel_j_code[3] = { (double)Pj.CBE_VelPred[0], (double)Pj.CBE_VelPred[1], (double)Pj.CBE_VelPred[2] };
        cbe_build_flux_frame_Q_from_stored_moments(
            Pj.CBE_basis_moments_pred, Vel_j_code, V_j, cf_a3inv, cf_atime, Q_j);
    }

    /* Basis matching via SSOT helper. Same comment as the LSQ pair body:
     * matching API stays on Q (moment rows); cost helpers derive
     * primitive content internally, so the primitive-gradient swap does
     * not require any matching-API change. Single-direction; counter NULL. */
    int matched_j_for_i[CBE_INTEGRATOR_NBASIS];
    /* Stale-gradient reconstructed face densities for the density-continuity
     * matching cost (gradient pass: vF=0). psi_i = h_j/(h_i+h_j) from above.
     * NULL when gradients are not compiled in -> base velocity/trace cost. */
    const double *gcrho_a = (const double*)0, *gcrho_b = (const double*)0;
#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
    double gcrho_fi[CBE_INTEGRATOR_NBASIS], gcrho_fj[CBE_INTEGRATOR_NBASIS];
    {
        for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
            double gi = L.Gradients_CBE_basis_moments_prev[m][0][0]*dp[0]
                      + L.Gradients_CBE_basis_moments_prev[m][0][1]*dp[1]
                      + L.Gradients_CBE_basis_moments_prev[m][0][2]*dp[2];
            double gj = Pj.Gradients_CBE_basis_moments_prev[m][0][0]*dp[0]
                      + Pj.Gradients_CBE_basis_moments_prev[m][0][1]*dp[1]
                      + Pj.Gradients_CBE_basis_moments_prev[m][0][2]*dp[2];
            gcrho_fi[m] = Q_i[m][0] - psi_i*gi;
            gcrho_fj[m] = Q_j[m][0] + (1.0 - psi_i)*gj;
        }
        gcrho_a = gcrho_fi; gcrho_b = gcrho_fj;
    }
#endif
    const double vface_zero[3] = {0.0, 0.0, 0.0};   /* gradient stencil: no single face normal, no face shift */
    cbe_build_pair_matching(Q_i, Q_j, matched_j_for_i,
                            /*alpha_of_beta_for_b=*/NULL,
                            /*free_slot_fired_count_inout=*/NULL,
                            gcrho_a, gcrho_b, vface_zero);

    /* Primitive-gradient swap: COMPONENT-SEPARATED limiter.
     * The realizability constraint factorizes
     * in primitive variables (ρ²·S ≥ 0) so density, velocity, and stress
     * limiters are INDEPENDENT. The previous row-scalar shape was a
     * leftover from the moment-row coupling — applying one phi to all
     * slots let irrelevant roundoff in grad_v / grad_S veto grad_ρ.
     *
     * Per slot:
     *   ρ slot:    BJ + positivity (ρ_face > ρ_floor) — affects accum.phi[m][0] only.
     *   v_k slot:  BJ only, no realizability role — affects accum.phi[m][1..NDIM].
     *   S block:   per-slot BJ + tensor PSD (single φ_S applied uniformly
     *              across all stress slots to keep S_face on the PSD cone)
     *              — affects accum.phi[m][1+NDIM..NMOMENTS-1].
     *
     * Scale-aware BJ tolerance uses physically meaningful scales:
     *   scale_ρ   = max(|ρ_i|, |ρ_j|, ρ_floor)
     *   scale_v_k = max(|v_i,k|, |v_j,k|, sqrt(max(R_kk_i, R_kk_j, 0)))
     *                with R_kk = T_kk / ρ
     *   scale_S_kl = max(|S_kl_i|, |S_kl_j|, |R_kl_i|, |R_kl_j|)
     * So cold/warm-small slots (uniform v, σ²-tiny S) don't fall through
     * the BJ tol from LSQ roundoff noise — preserved as no-constraint
     * (φ = 1) on those slots.
     *
     * Scratch-row skip: if cell i basis m has ρ_i,m ≤ floor (pass-0
     * wrote zero gradient for this row), leave phi at neutral 1.0. */
    const double rho_floor = cbe_rho_active_floor();
    double prim_i[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    double prim_j[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
        cbe_moments_to_primitives_row(Q_i[m], prim_i[m]);
        cbe_moments_to_primitives_row(Q_j[m], prim_j[m]);
    }

    for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
        if(!(Q_i[m][0] > rho_floor)) continue;
        const int n = matched_j_for_i[m];

        /* Per-slot predicted face delta from cell-i primitive gradient. */
        double g_row[CBE_INTEGRATOR_NMOMENTS];
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
            const double gdotdp = L.Gradients_CBE_basis_moments[m][k][0]*dp[0]
                                + L.Gradients_CBE_basis_moments[m][k][1]*dp[1]
                                + L.Gradients_CBE_basis_moments[m][k][2]*dp[2];
            g_row[k] = -psi_i * gdotdp;
        }

        /* ----- ρ slot (k = 0): BJ + positivity ----- */
        {
            const double scale_rho =
                DMAX(DMAX(fabs(prim_i[m][0]), fabs(prim_j[n][0])), rho_floor);
            const double dprim     = prim_j[n][0] - prim_i[m][0];
            double phi_rho = cbe_bj_phi_pair_scaled(dprim, g_row[0], scale_rho);
            /* Positivity: ρ_face(φ) = ρ_i + φ · g_row[0] > ρ_floor.
             * If g_row[0] ≥ 0, ρ_face only grows — no constraint.
             * If g_row[0] < 0, max φ maintaining ρ_face > floor is
             *   (ρ_i − floor) / |g_row[0]|. */
            if(g_row[0] < 0) {
                const double slack = prim_i[m][0] - rho_floor;
                if(slack <= 0) {
                    phi_rho = 0;
                } else {
                    const double phi_pos = slack / (-g_row[0]);
                    if(phi_pos < phi_rho) phi_rho = phi_pos;
                }
            }
            if(phi_rho < 0) phi_rho = 0;
            if(phi_rho > 1) phi_rho = 1;
            if(phi_rho < accum.phi[m][0]) accum.phi[m][0] = phi_rho;
        }

        /* ----- v slots (k = 1 .. NUMDIMS): BJ only ----- */
        for(int k = 1; k <= NUMDIMS && k < CBE_INTEGRATOR_NMOMENTS; k++) {
            /* R_kk_i = T_kk_i / ρ_i for the v-scale. The T_kk slot is at
             * cbe_T_idx(k-1, k-1) when SECONDMOMENT is on; otherwise the
             * R-component falls back to 0 (no T data). */
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
            const int t_kk_idx = cbe_T_idx(k-1, k-1);
            const double R_kk_i =
                (prim_i[m][0] > rho_floor && t_kk_idx < CBE_INTEGRATOR_NMOMENTS)
                    ? (Q_i[m][t_kk_idx] / prim_i[m][0]) : 0;
            const double R_kk_j =
                (prim_j[n][0] > rho_floor && t_kk_idx < CBE_INTEGRATOR_NMOMENTS)
                    ? (Q_j[n][t_kk_idx] / prim_j[n][0]) : 0;
            const double R_scale = sqrt(DMAX(DMAX(R_kk_i, R_kk_j), 0.0));
#else
            const double R_scale = 0.0;
#endif
            const double scale_v_k =
                DMAX(DMAX(fabs(prim_i[m][k]), fabs(prim_j[n][k])), R_scale);
            const double dprim = prim_j[n][k] - prim_i[m][k];
            /* v-slope BJ limiter with finite opposite-sign tolerance CBE_VOPPSIGN. */
            double phi_v = cbe_bj_phi_pair_scaled(dprim, g_row[k], scale_v_k, CBE_VOPPSIGN);
            if(phi_v < 0) phi_v = 0;
            if(phi_v > 1) phi_v = 1;
            if(phi_v < accum.phi[m][k]) accum.phi[m][k] = phi_v;
        }

#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        /* ----- S block (stress slots): per-slot BJ + tensor PSD ----- */
        {
            const int s_start = 1 + NUMDIMS;
            const int s_end   = CBE_INTEGRATOR_NMOMENTS;

            /* Per-slot BJ over stress block — find tightest (block-min) BJ φ. */
            double phi_S_BJ_block = 1.0;
            for(int k = s_start; k < s_end; k++) {
                const double R_kl_i = (prim_i[m][0] > rho_floor) ? (Q_i[m][k] / prim_i[m][0]) : 0;
                const double R_kl_j = (prim_j[n][0] > rho_floor) ? (Q_j[n][k] / prim_j[n][0]) : 0;
                const double scale_S = DMAX(DMAX(fabs(prim_i[m][k]), fabs(prim_j[n][k])),
                                            DMAX(fabs(R_kl_i),         fabs(R_kl_j)));
                const double dprim = prim_j[n][k] - prim_i[m][k];
                const double phi_BJ_kl = cbe_bj_phi_pair_scaled(dprim, g_row[k], scale_S);
                if(phi_BJ_kl < phi_S_BJ_block) phi_S_BJ_block = phi_BJ_kl;
            }

            /* PSD: tighten φ_S so S_face = S_i + φ_S · S_delta stays PSD on
             * the active NUMDIMS x NUMDIMS block. NMOMENTS=3 (1D) admits a
             * closed-form positivity bound on the scalar S_xx; NMOMENTS=10
             * (3D) requires bisection (off-diagonal coupling makes PSD a
             * non-linear condition in φ). */
            double phi_S;
            if(CBE_INTEGRATOR_NMOMENTS == 3) {
                /* 1D: scalar S_xx ≥ 0. Linear in φ → closed form. */
                const double S_xx_i  = prim_i[m][s_start];
                const double S_delta = g_row[s_start];
                phi_S = phi_S_BJ_block;
                if(S_delta < 0) {
                    if(S_xx_i <= 0) {
                        phi_S = 0;
                    } else {
                        const double phi_pos = S_xx_i / (-S_delta);
                        if(phi_pos < phi_S) phi_S = phi_pos;
                    }
                }
            } else {
                /* 3D (or future NMOMENTS) — load active S_i, S_delta as 3×3
                 * symmetric blocks and bisect on PSD. */
                double S_i_mat[3][3], S_delta_mat[3][3];
                cbe_basis_T_load_active(prim_i[m], S_i_mat);
                for(int a = 0; a < 3; a++)
                    for(int b = 0; b < 3; b++)
                        S_delta_mat[a][b] = 0.0;
                for(int a = 0; a < NUMDIMS; a++) {
                    const int kk = cbe_T_idx(a, a);
                    if(kk < CBE_INTEGRATOR_NMOMENTS) S_delta_mat[a][a] = g_row[kk];
                    for(int b = a + 1; b < NUMDIMS; b++) {
                        const int kab = cbe_T_idx(a, b);
                        if(kab < CBE_INTEGRATOR_NMOMENTS) {
                            S_delta_mat[a][b] = g_row[kab];
                            S_delta_mat[b][a] = S_delta_mat[a][b];
                        }
                    }
                }
                phi_S = cbe_S_psd_phi_bisect(S_i_mat, S_delta_mat, phi_S_BJ_block);
            }
            if(phi_S < 0) phi_S = 0;
            if(phi_S > 1) phi_S = 1;
            /* Apply single φ_S to every stress slot in the row. Per-slot
             * φ would move S_face off the PSD cone on off-diagonals. */
            for(int k = s_start; k < s_end; k++) {
                if(phi_S < accum.phi[m][k]) accum.phi[m][k] = phi_S;
            }
        }
#endif  /* CBE_INTEGRATOR_SECONDMOMENT */
    }
    (void)j;
}

/* ============================================================================
 * CBEGradSpec — runner-template Spec contract. Non-iterative; two passes
 * orchestrated by the toplevel via Aux::loop_iteration.
 * ========================================================================== */
struct CBEGradSpec {
    static constexpr const char *loop_name = "cbe_grad";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (M/B lsq + phi min); no j-write, no atomics -> threaded eval bit-identical */

    /* Symmetric search topology because pass 1 (pairwise BJ-style limiter)
     * needs symmetric neighbor visibility to see every face the flux body
     * will subsequently reconstruct — exactly the role hydro GradientsSpec
     * plays for its slope-limiter / min-max collection. The raw LSQ pass
     * (pass 0) intentionally restricts accumulation to the i-kernel inside
     * the pair body (see cbe_grad_lsq_pair_kernel_body), matching the
     * density / NV_T / hydro raw-gradient semantics (linear consistency of
     * the LSQ estimator relies on the i-kernel normalisation). Per-pass
     * filtering is the right shape; do NOT "fix" pass 0 to symmetric. */
    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    /* Nominal default — DM, matching DMGradSpec convention. CBE runs on
     * AGS-eligible types alongside AGSForce; the real per-call mask is set
     * at the toplevel via args.neighbor_type_mask_override = bm
     * (ags_gravity_kernel_shared_BITFLAG). Conservative default avoids the
     * footgun of "missed override → all types become neighbors". */
    static constexpr unsigned int            neighbor_type_mask = (1u << 1);  /* DM */
    /* CBE-grad pair physics: h_j reach is P[j].AGS_KernelRadius for non-gas (DM)
     * — second symmetric pair body at line 399 reads h_j = Pj.AGS_KernelRadius
     * and filters via r > max(h_i, h_j) at line 408.  Gas excluded by
     * neighbor_type_mask. */
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_NONGAS_AGS;

    static constexpr WritePattern   write_pattern             = WritePattern::ActiveReduceOnly;
    /* SidxCacheKind::None: the per-call DM/SIDM mask varies per bm — same
     * reasoning as DMGradSpec. */
    static constexpr SidxCacheKind  sidx_cache_kind           = SidxCacheKind::None;
    static constexpr bool           uses_ghost_writeback      = false;
    static constexpr bool           uses_ghost_write_detector = false;

    /* AGSForce_isactive is the SAME predicate the force consumer uses —
     * gradients refresh for exactly the particles whose force will read
     * them this step. Inactive particles retain their previous-step
     * P[i].Gradients_CBE_basis_moments (hydro semantics). */
    static bool is_active(int i) { return AGSForce_isactive(i) != 0; }

    using IterControl    = NotIterative;
    using CallScalars    = CBEGradCallScalars;
    using ActiveData     = CBEGradActiveState;
    using AccumData      = CBEGradOut;
    using DeviceContext  = NeighborLoopDeviceContextBase;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;

    /* NeighborData carries `j` and the P base pointer because the pair
     * body needs get_particle_volume_ags_P(j, P_base). */
    struct NeighborData {
        struct particle_data *neighbor_particle;
        int                   j;
        struct particle_data *P_base;
    };

    /* Aux carries the external gradient-pass index (mirror DMGradSpec).
     * populate_call_scalars copies it into CallScalars; pair_kernel +
     * apply_active_writeback branch on scalars.loop_iteration. */
    struct Aux {
        int loop_iteration;
    };

    static_assert(uses_ghost_writeback == false,
        "CBEGrad is pure i-side accumulation; no j-side writes.");

    /* Host hooks — bodies in sidm/cbe_integrator_gradients.cc. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    /* merge_accum — additive for M + B; min-merge for phi (NaN -> 0). Unused-pass
     * fields are at their neutral elements (0 for sum, 1 for min) so this single
     * implementation is correct for both passes. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local, const AccumData& peer)
    {
        for(int a = 0; a < 3; a++)
            for(int b = 0; b < 3; b++)
                local.M[a][b] += peer.M[a][b];
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                for(int d = 0; d < 3; d++)
                    local.B[m][k][d] += peer.B[m][k][d];
                const double p = peer.phi[m][k];
                double      &L = local.phi[m][k];
                if(!isfinite(p)) { L = 0.0; continue; }
                if(p < L) L = p;
            }
        }
    }

    /* ---- Device hooks (header-inline; runner instantiates from GPU TUs). ---- */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        for(int a = 0; a < 3; a++)
            for(int b = 0; b < 3; b++)
                accum.M[a][b] = 0.0;
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                for(int d = 0; d < 3; d++) accum.B[m][k][d] = 0.0;
                /* phi initialized to 1.0 for min-tracking; pass 0 never
                 * touches it so it stays neutral for that pass's merge. */
                accum.phi[m][k] = 1.0;
            }
        }
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& dctx,
                                   int /*active_slot*/, int i,
                                   double h_search,
                                   const CallScalars& scalars)
    {
        ActiveData a;
        a.pos[0]           = (double)dctx.P[i].Pos[0];
        a.pos[1]           = (double)dctx.P[i].Pos[1];
        a.pos[2]           = (double)dctx.P[i].Pos[2];
        a.h_search         = h_search;
        a.TimeBin          = dctx.P[i].TimeBin;
        (void)scalars;
        a.origin_local_idx = i;
        a.origin_rank      = -1;

        CBEGradLocalIn& L = a.local;
        L.Mass             = dctx.P[i].Mass;
        L.AGS_KernelRadius = (double)dctx.P[i].AGS_KernelRadius;
        L.Pos              = a.pos;
        /* Build active-particle gradients from the PREDICTED (drifted) state,
         * matching hydro (which constructs active gradients from VelPred) and
         * the predicted flux reconstruction. pred is maintained everywhere by
         * the predictor (init seed + post-kick reset + per-drift advance). */
        L.Vel[0] = (double)dctx.P[i].CBE_VelPred[0];
        L.Vel[1] = (double)dctx.P[i].CBE_VelPred[1];
        L.Vel[2] = (double)dctx.P[i].CBE_VelPred[2];
        L.V_i    = get_particle_volume_ags_P(i, dctx.P);
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++)
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                L.U_i[m][k] = dctx.P[i].CBE_basis_moments_pred[m][k];

        /* Mode-B-safe by-value snapshot of the persistent gradient row.
         * Pass-1 limiter reads this from
         * active.local; pass-0 LSQ does not consume it but populating
         * unconditionally is fixed-cost and keeps the load-active body
         * branch-free. */
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++)
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                for(int d = 0; d < 3; d++) {
                    L.Gradients_CBE_basis_moments[m][k][d] =
                        dctx.P[i].Gradients_CBE_basis_moments[m][k][d];
                    L.Gradients_CBE_basis_moments_prev[m][k][d] =
                        dctx.P[i].Gradients_CBE_basis_moments_prev[m][k][d];
                }
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& dctx,
                                       int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/)
    {
        NeighborData n;
        n.neighbor_particle = &dctx.P[j];
        n.j                 = j;
        n.P_base            = dctx.P;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData&   active,
                             const NeighborData& neighbor,
                             AccumData&          accum,
                             NoScatter& /*scatter*/,
                            const CallScalars& cs)
    {
        if(cs.loop_iteration <= 0) {
            cbe_grad_lsq_pair_kernel_body(active, cs, neighbor, accum);
        } else {
            cbe_grad_bj_pair_kernel_body (active, cs, neighbor, accum);
        }
    }
};

/* Toplevel runner instantiation. Non-iterative form — paralleling DMGrad. */
extern template void run_neighbor_loop<CBEGradSpec>(const neighbor_loop_args&);

#endif /* CBE_INTEGRATOR_WITHGRADIENTS */
#endif /* CBE_INTEGRATOR_GRADIENTS_H */
