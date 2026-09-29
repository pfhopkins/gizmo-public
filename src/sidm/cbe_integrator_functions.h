/* cbe_integrator_functions.h — GPU-callable CBE (Collisionless Boltzmann
 * Equation) per-basis flux + face helpers + drift-kick. Pure math with no
 * global state; called both from the CPU tree-walk path (via
 * cbe_integrator_flux_functions.h) and from the AGSForce GPU kernel as
 * KOKKOS_INLINE_FUNCTION.
 *
 * SSOT: the single per-basis flux implementation lives here as
 * cbe_flux_tophat_vacuum. There is no CPU duplicate in cbe_integrator.cc.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef CBE_INTEGRATOR_FUNCTIONS_H
#define CBE_INTEGRATOR_FUNCTIONS_H

#include "../declarations/allvars.h"
#include "../declarations/gpu_rng.h"

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#ifdef CBE_INTEGRATOR
/* Per-loop FNV-1a salt for the drift-kick basis-resplit RNG. Mixed into the
 * counter so the CBE drift-kick stream is independent of any other loop that
 * happens to share (Ti_Current, ID, tag). See gpu_rng.h for rationale. */
static constexpr uint64_t CBE_DRIFT_KICK_RNG_SALT = gizmo_loop_rng_salt("cbe_drift_kick");

/* SPD repair eigenvalue floor, relative to
 * trace. Applied identically in face-state Q-clamp and the
 * conservative cell-state repair (both via the SSOT helper
 * cbe_basis_row_project_central_stress_to_PSD). 1e-12 is a small relative
 * floor (well below the intended stress scale, but large enough to
 * deterministically exclude strictly-zero / negative eigenvalues from
 * downstream consumers). */
static constexpr double CBE_SPD_RELATIVE_FLOOR = 1.0e-12;

/* Tunings for the conservative cell-state
 * realizability repair. EPS_M_ABS / EPS_M_REL: empty-slot mass threshold
 * (m_thresh = max(EPS_M_ABS, EPS_M_REL * M_cell / NBASIS)). EPS_CONS:
 * FP-scale residual acceptance bound (residual norm <= EPS_CONS * cell_scale
 * after donation -> accept as roundoff; larger -> deterministic conservative
 * fallback in Step 7). */
static constexpr double CBE_REPAIR_EPS_M_ABS = 1.0e-30;
static constexpr double CBE_REPAIR_EPS_M_REL = 1.0e-12;
static constexpr double CBE_REPAIR_EPS_CONS  = 1.0e-12;

/* Aggregate per-basis outflow-budget limiter.
 * FRAC: per-basis outgoing-mass budget as fraction of m_a, capped per drift
 * call. MIN_DONOR_MASS: floor below which a basis is ineligible as a donor
 * AND as a capped-basis recipient of the (1-f)*excess "uncapped remainder"
 * (treated as floor for both ends). 0.9 leaves headroom for the nfac threshold
 * (-0.75) downstream + FP. */
static constexpr double CBE_OUTFLOW_BUDGET_FRAC           = 0.9;
static constexpr double CBE_OUTFLOW_BUDGET_MIN_DONOR_MASS = 1.0e-30;

/* Per-basis outflow-budget limiter status codes (orthogonal to the
 * cell-state repair codes). 200+ for non-error outcomes; no hard-fails / endruns --
 * donor-limited cases surface via the counter, not abort. */
enum {
    CBE_OUTBUDGET_OK              = 200,   /* nothing exceeded budget          */
    CBE_OUTBUDGET_APPLIED         = 201,   /* at least one basis capped, donors clean */
    CBE_OUTBUDGET_DONOR_LIMITED   = 202,   /* cap applied at fraction f<1; (1-f)*excess unredistributed */
    CBE_OUTBUDGET_NO_DONOR        = 203    /* no eligible donor; this basis a uncapped */
};

/* POD diagnostic counters for the outflow-budget limiter. Allocated only
 * when diagnostics are enabled (matching the cell-state repair counts gating).
 * NOT promoted to a cbe_diagnostics.txt column; accessible via the optional
 * CBE_INTEGRATOR_OUTBUDGET_VERBOSE printf path or via debugger. */
struct CbeOutflowBudgetDiag {
    long long n_capped_bases;          /* SUM: bases that hit f_out cap (status APPLIED) */
    long long n_donor_limited;         /* SUM: bisection capped fraction f < 1 (DONOR_LIMITED) */
    long long n_donor_failures;        /* SUM: NO_DONOR + DONOR_LIMITED bookkeeping */
    double    max_pre_out_frac_seen;   /* MAX: pre-cap dt*Udot_out[a][0]/m_a seen */
};

/* Dimension-aware momentum-slot accessors. Stored basis-moment
 * row layout is row[0] = mass, row[1..NUMDIMS] = p_x..p_{NUMDIMS-1}, and
 * (under CBE_INTEGRATOR_SECONDMOMENT) row[1+NUMDIMS..NMOMENTS-1] = stress
 * tensor components. NMOMENTS shrinks in 1D/2D so the y/z momentum slots
 * literally don't exist in the array. These helpers read them as 0 and
 * silently drop the write, matching the hydro convention where Vel[1] /
 * Vel[2] etc. are simply zero in low-D. The 3-vector dot products and
 * vector intermediates downstream (dp[3], Vel[3], Face_Area_Vec[3])
 * remain 3-wide so the same flux/update math compiles for any NUMDIMS.
 * Compile-time NUMDIMS branch -> zero overhead at NUMDIMS=3. */
KOKKOS_INLINE_FUNCTION double cbe_basis_p_r(const double *row, int k) {
    return (k < NUMDIMS) ? row[1 + k] : 0.0;
}
KOKKOS_INLINE_FUNCTION void cbe_basis_p_w(double *row, int k, double val) {
    if(k < NUMDIMS) row[1 + k] = val;
}
KOKKOS_INLINE_FUNCTION void cbe_basis_p_a(double *row, int k, double val) {
    if(k < NUMDIMS) row[1 + k] += val;
}
KOKKOS_INLINE_FUNCTION void cbe_basis_p_load_3(const double *row, double p3[3]) {
    p3[0] = cbe_basis_p_r(row, 0);
    p3[1] = cbe_basis_p_r(row, 1);
    p3[2] = cbe_basis_p_r(row, 2);
}
KOKKOS_INLINE_FUNCTION void cbe_basis_p_store_3(double *row, const double p3[3]) {
    cbe_basis_p_w(row, 0, p3[0]);
    cbe_basis_p_w(row, 1, p3[1]);
    cbe_basis_p_w(row, 2, p3[2]);
}
KOKKOS_INLINE_FUNCTION void cbe_basis_v_load_3(const double *row, double v3[3]) {
    /* v = p/m for the basis. Returns zeros if m <= 0. */
    const double m = row[0];
    const double inv_m = (m > 0) ? 1.0 / m : 0.0;
    v3[0] = cbe_basis_p_r(row, 0) * inv_m;
    v3[1] = cbe_basis_p_r(row, 1) * inv_m;
    v3[2] = cbe_basis_p_r(row, 2) * inv_m;
}

/* Symmetric stress-slot helpers (defined so the SECONDMOMENT path can be
 * unfenced for arbitrary D later without rewriting the hot path). Layout
 * per the existing convention:
 *   1D NMOMENTS=3 : [m, p_x, T_xx]                                   -> NSTRESS=1
 *   2D NMOMENTS=6 : [m, p_x, p_y, T_xx, T_yy, T_xy]                  -> NSTRESS=3
 *   3D NMOMENTS=10: [m, p_x, p_y, p_z, T_xx, T_yy, T_zz, T_xy, T_xz, T_yz] -> NSTRESS=6
 * Diagonals occupy 1+NUMDIMS..1+2*NUMDIMS-1; off-diagonals (upper triangle
 * in lex order (a,b) with a<b) follow. cbe_T_idx returns the slot index
 * for symmetric T_{a,b}; cbe_basis_T_r/_w guard the same way the
 * momentum accessors do. SECONDMOMENT non-3D still fenced today; these
 * helpers exist so the call sites that use them are unfence-safe. */
#define CBE_NSTRESS (NUMDIMS * (NUMDIMS + 1) / 2)
KOKKOS_INLINE_FUNCTION int cbe_T_idx(int a, int b) {
    if(a > b) { int tmp = a; a = b; b = tmp; }
    if(a == b) return 1 + NUMDIMS + a;
    const int off = a * (2 * NUMDIMS - a - 1) / 2 + (b - a - 1);
    return 1 + 2 * NUMDIMS + off;
}
KOKKOS_INLINE_FUNCTION double cbe_basis_T_r(const double *row, int a, int b) {
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    return row[cbe_T_idx(a, b)];
#else
    (void)row; (void)a; (void)b;
    return 0.0;
#endif
}
KOKKOS_INLINE_FUNCTION void cbe_basis_T_w(double *row, int a, int b, double val) {
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    row[cbe_T_idx(a, b)] = val;
#else
    (void)row; (void)a; (void)b; (void)val;
#endif
}

/* Active-dimensional symmetric stress workspace helpers.
 *
 * Stress math operates on a 3x3 workspace tensor for shape uniformity (matches
 * the 3-vector convention for dp[3]/Vel[3]/Face_Area_Vec[3]), but ONLY the
 * upper-left NUMDIMS x NUMDIMS block carries meaningful values. Off-block
 * entries are kept at exactly zero so that downstream consumers (Jacobi PSD
 * projection, principal-minors realizability test, n.S.n contractions) reduce
 * to the correct active-dim physics without per-D special-casing in callers.
 *
 * These helpers carry the dim-aware layout to all stress consumers so the
 * SECONDMOMENT non-3D fence at declarations/precompiler_logic.h can be lifted.
 *
 * Semantics:
 *   _load_active(row, T)  -- zero T[3][3], then fill the active D x D upper
 *                            triangle (symmetrized) from row[T_idx(a,b)].
 *   _load_active_scaled(row, scale, T)  -- same but pre-multiplied by `scale`
 *                            (useful for R = T/m).
 *   _store_active(row, T) -- write T_ab for active (a<=b)<NUMDIMS back into
 *                            row[T_idx(a,b)]. Inactive entries ignored.
 *   _trace_active(row, scale) -- sum_{a<NUMDIMS} row[T_idx(a,a)] * scale.
 *
 * Realizability + PSD projection helpers operate on the same 3x3 workspace
 * but ignore inactive rows/cols, so e.g. in 1D the Jacobi sweep is a no-op
 * and only S[0][0] participates in the eigenvalue floor / clamp. */
KOKKOS_INLINE_FUNCTION void cbe_basis_T_load_active(
    const double *row, double T[3][3])
{
    for(int a=0; a<3; a++) for(int b=0; b<3; b++) T[a][b] = 0.0;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        T[a][a] = row[cbe_T_idx(a, a)];
        for(int b=a+1; b<NUMDIMS; b++) {
            const double Tab = row[cbe_T_idx(a, b)];
            T[a][b] = Tab; T[b][a] = Tab;
        }
    }
#else
    (void)row;
#endif
}
KOKKOS_INLINE_FUNCTION void cbe_basis_T_load_active_scaled(
    const double *row, double scale, double T[3][3])
{
    cbe_basis_T_load_active(row, T);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        T[a][a] *= scale;
        for(int b=a+1; b<NUMDIMS; b++) { T[a][b] *= scale; T[b][a] = T[a][b]; }
    }
#else
    (void)scale;
#endif
}
KOKKOS_INLINE_FUNCTION void cbe_basis_T_store_active(
    double *row, const double T[3][3])
{
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        row[cbe_T_idx(a, a)] = T[a][a];
        for(int b=a+1; b<NUMDIMS; b++) row[cbe_T_idx(a, b)] = T[a][b];
    }
#else
    (void)row; (void)T;
#endif
}
KOKKOS_INLINE_FUNCTION double cbe_basis_T_trace_active(
    const double *row, double scale)
{
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double tr = 0.0;
    for(int a=0; a<NUMDIMS; a++) tr += row[cbe_T_idx(a, a)];
    return tr * scale;
#else
    (void)row; (void)scale;
    return 0.0;
#endif
}

/* ===== Primitive (ρ, v, S) row helpers — primitive-gradient swap support =====
 *
 * The primitive-gradient swap will make the gradient writer and face
 * reconstruction operate on primitive rows (ρ, v_k, S_kl) rather than moment
 * rows (m, p_k, T_kl). These helpers convert between the two packings; the
 * slot layout of the primitive row matches the existing moment-row packing
 * exactly, slot-for-slot:
 *   slot 0           : ρ (= m)
 *   slots 1..NUMDIMS : v_k = p_k / ρ
 *   stress slots     : S_kl = T_kl / ρ - v_k v_l         (SECONDMOMENT only)
 * The stress slots are indexed by cbe_T_idx() exactly as in the moment row.
 *
 * Identity (load-bearing for the swap):
 *     ρ · T_kl − p_k p_l ≡ ρ² · S_kl
 * Moment-cone realizability (m·T − p·p ⪰ 0) is therefore exactly equivalent
 * to primitive realizability (ρ ≥ ρ_floor, S ⪰ 0). The face limiter exploits
 * this.
 *
 * The persistent storage field is still spelled Gradients_CBE_basis_moments
 * but holds primitive-gradient slots after the swap. The field-name rename
 * is out of scope; see the particle_data.h doc.
 *
 * These helpers are defined but not yet used: no gradient writer / face
 * reader / limiter behavior change from their presence. */

/* Single SSOT for "is this row's density active enough to divide by?".
 * Used at every primitive-derivation site (primitive conversion, scratch-row
 * gradient gating, face-clamp consistency check). Matches the existing
 * MIN_REAL_NUMBER threshold used by cbe_clamp_face_Q and the pair-matching
 * cost helpers (cbe_cost_v_only, cbe_cost_trace_w2). Do NOT introduce a
 * second density floor — grep this name to find every site that gates on
 * "active enough to divide". */
KOKKOS_INLINE_FUNCTION
double cbe_rho_active_floor(void) { return MIN_REAL_NUMBER; }

/* Moment row -> primitive row. Both arrays sized [CBE_INTEGRATOR_NMOMENTS]
 * with the slot packing described above. Scratch rows (Q_row[0] below the
 * active-density floor) get prim_row[0] = Q_row[0] (the raw density, preserved
 * for downstream gating) and all v / S slots set to zero — the caller MUST
 * NOT consume those slots without first re-checking the floor. The primitive
 * conversion never divides by a sub-floor density; the explicit zeros prevent
 * NaN/Inf from leaking into downstream consumers.
 *
 * Inactive dimensions (k >= NUMDIMS) are not packed into the row at all
 * (NMOMENTS shrinks in low-D, matching the moment-row convention). The loop
 * bounds use NUMDIMS so the compiler folds away the inactive iterations. */
KOKKOS_INLINE_FUNCTION
void cbe_moments_to_primitives_row(const double Q_row[CBE_INTEGRATOR_NMOMENTS],
                                   double prim_row[CBE_INTEGRATOR_NMOMENTS])
{
    const double rho = Q_row[0];
    prim_row[0] = rho;
    /* Strict `>` matches the existing active-row predicate used by
     * cbe_clamp_face_Q (Q_face[m][0] <= MIN_REAL_NUMBER -> inactive/clamped).
     * SSOT consistency: same edge convention everywhere. */
    const bool active = (rho > cbe_rho_active_floor());
    const double inv_rho = active ? (1.0 / rho) : 0.0;

    double v_tmp[3] = {0.0, 0.0, 0.0};
    for(int k=0; k<NUMDIMS; k++) {
        const double v_k = cbe_basis_p_r(Q_row, k) * inv_rho;
        v_tmp[k] = v_k;
        cbe_basis_p_w(prim_row, k, active ? v_k : 0.0);
    }

#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        const double T_aa = cbe_basis_T_r(Q_row, a, a);
        const double S_aa = active ? (T_aa * inv_rho - v_tmp[a] * v_tmp[a]) : 0.0;
        cbe_basis_T_w(prim_row, a, a, S_aa);
        for(int b=a+1; b<NUMDIMS; b++) {
            const double T_ab = cbe_basis_T_r(Q_row, a, b);
            const double S_ab = active ? (T_ab * inv_rho - v_tmp[a] * v_tmp[b]) : 0.0;
            cbe_basis_T_w(prim_row, a, b, S_ab);
        }
    }
#endif
}

/* Primitive row -> moment row. Deterministic, no iteration. The inverse of
 * cbe_moments_to_primitives_row at every site where the primitive row was
 * built from a row with ρ > cbe_rho_active_floor(). At or below the floor
 * the caller must NOT use this helper — for sub-floor rows the face
 * reconstruction stays first-order moment, bypassing the round-trip.
 *
 * Identity: applying primitives_to_moments after moments_to_primitives is
 * algebraically exact for any row with ρ > floor. FP cancellation can still
 * occur in the forward step (S = T/ρ − vv for cold or warm-small S); the
 * inverse step just multiplies by ρ and adds vv back, so the round-trip
 * recovers the original moment row only up to the cancellation error in S. */
KOKKOS_INLINE_FUNCTION
void cbe_primitives_to_moments_row(const double prim_row[CBE_INTEGRATOR_NMOMENTS],
                                   double Q_row[CBE_INTEGRATOR_NMOMENTS])
{
    const double rho = prim_row[0];
    Q_row[0] = rho;

    double v_tmp[3] = {0.0, 0.0, 0.0};
    for(int k=0; k<NUMDIMS; k++) {
        const double v_k = cbe_basis_p_r(prim_row, k);
        v_tmp[k] = v_k;
        cbe_basis_p_w(Q_row, k, rho * v_k);
    }

#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        const double S_aa = cbe_basis_T_r(prim_row, a, a);
        cbe_basis_T_w(Q_row, a, a, rho * (S_aa + v_tmp[a] * v_tmp[a]));
        for(int b=a+1; b<NUMDIMS; b++) {
            const double S_ab = cbe_basis_T_r(prim_row, a, b);
            cbe_basis_T_w(Q_row, a, b, rho * (S_ab + v_tmp[a] * v_tmp[b]));
        }
    }
#endif
}

/* Active-dim symmetric PSD test via principal minors of the upper-left
 * NUMDIMS x NUMDIMS block. Reduces to scalar nonnegativity in 1D, the
 * 2x2 Sylvester chain in 2D, and the full 3x3 chain in 3D. eps_tol is
 * relative to trace_pos (active diagonals only) at each minor order. */
KOKKOS_INLINE_FUNCTION
bool cbe_symNxN_active_principal_minors_nonneg(const double M[3][3], double eps_tol)
{
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double trace_pos = MIN_REAL_NUMBER;
    for(int a=0; a<NUMDIMS; a++) { trace_pos += M[a][a]; }
    if(trace_pos < MIN_REAL_NUMBER) trace_pos = MIN_REAL_NUMBER;
    const double eps_diag = eps_tol * trace_pos;
    for(int a=0; a<NUMDIMS; a++) { if(M[a][a] < -eps_diag) return false; }
#if (NUMDIMS >= 2)
    const double eps_2x2 = eps_tol * trace_pos * trace_pos;
    for(int a=0; a<NUMDIMS; a++)
        for(int b=a+1; b<NUMDIMS; b++)
            if(M[a][a]*M[b][b] - M[a][b]*M[a][b] < -eps_2x2) return false;
#endif
#if (NUMDIMS >= 3)
    const double eps_det = eps_tol * trace_pos * trace_pos * trace_pos;
    const double det = M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[1][2])
                     - M[0][1] * (M[0][1]*M[2][2] - M[1][2]*M[0][2])
                     + M[0][2] * (M[0][1]*M[1][2] - M[1][1]*M[0][2]);
    if(det < -eps_det) return false;
#endif
    return true;
#else
    (void)M; (void)eps_tol;
    return true;
#endif
}

/* Active-dim symmetric PSD projection in-place on a 3x3 workspace. Outside
 * the active block, entries are assumed zero and are NOT touched (no eigen-
 * value floor on inactive dims -- the trace correction dT_out only counts
 * active diagonals). Algorithm: Jacobi sweep over active (a<b<NUMDIMS) off-
 * diagonals (no-op in 1D), clamp active diagonals to >= eigenvalue_floor,
 * reconstruct via V D V^T on the active block.
 *
 * 3D path: existing 6-sweep over (0,1)/(0,2)/(1,2), matches the pre-commit
 * cbe_project_central_S_to_PSD behavior. 2D: 6 sweeps over (0,1) only --
 * conservative-overkill on a 2x2 (one sweep would suffice; 6 are cheap).
 * 1D: no off-diagonal sweep; floor on S[0][0].
 *
 * Caller provides the eigenvalue_floor (e.g. CBE_SPD_RELATIVE_FLOOR * trace
 * for the face-clamp policy, or 0 for cell-state repair true-PSD). dT_out
 * (nullable) receives the post-minus-pre active-trace change. Non-finite
 * input triggers the same MIN_REAL_NUMBER fallback on active diagonals as
 * the 3D path used. */
KOKKOS_INLINE_FUNCTION
bool cbe_project_central_S_active_to_PSD(double S[3][3], double eigenvalue_floor,
                                         double *dT_out)
{
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    /* Self-protecting contract: explicitly zero inactive rows/cols on entry
     * so callers that hand in stack workspaces with stale junk cannot
     * pollute the active PSD test or the post-reconstruction store. */
    for(int a=0; a<3; a++) for(int b=0; b<3; b++)
        if(a >= NUMDIMS || b >= NUMDIMS) S[a][b] = 0.0;

    /* Pre-trace over active diagonals. */
    double trace_before = 0.0;
    bool all_finite = true;
    for(int a=0; a<NUMDIMS; a++) {
        if(!isfinite(S[a][a])) { all_finite = false; }
        else trace_before += S[a][a];
        for(int b=a+1; b<NUMDIMS; b++)
            if(!isfinite(S[a][b])) { all_finite = false; }
    }
    if(!all_finite || !isfinite(trace_before)) {
        /* Reset active block to the same MIN_REAL_NUMBER isotropic-floor
         * pattern as the old 3D helper. dT_out reports the post-reset
         * trace; we do NOT compute a delta against a non-finite pre-trace
         * (NaN/Inf would silently poison repair_dT_sum diagnostics). */
        double tr_after = 0.0;
        for(int a=0; a<NUMDIMS; a++) {
            S[a][a] = MIN_REAL_NUMBER;
            tr_after += MIN_REAL_NUMBER;
            for(int b=a+1; b<NUMDIMS; b++) { S[a][b] = 0.0; S[b][a] = 0.0; }
        }
        if(dT_out) *dT_out = tr_after;
        return true;
    }

    /* Working copies for Jacobi: A starts as S on the active block (others 0),
     * V starts as identity on the active block (others 0 -- inactive rows/cols
     * are not rotated). */
    double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double V[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int a=0; a<NUMDIMS; a++) {
        for(int b=0; b<NUMDIMS; b++) A[a][b] = S[a][b];
        V[a][a] = 1.0;
    }

#if (NUMDIMS >= 2)
    /* 6 fixed Jacobi sweeps over active off-diagonals (a<b<NUMDIMS). */
    for(int sweep=0; sweep<6; sweep++) {
        for(int p=0; p<NUMDIMS; p++) for(int q=p+1; q<NUMDIMS; q++) {
            double Apq = A[p][q];
            if(fabs(Apq) < 1.0e-300) continue;
            double theta = (A[q][q] - A[p][p]) / (2.0 * Apq);
            double t;
            if(fabs(theta) > 1.0e15) {
                t = 0.5 / theta;
            } else {
                double sgn = (theta >= 0) ? 1.0 : -1.0;
                t = sgn / (fabs(theta) + sqrt(theta*theta + 1.0));
            }
            const double c   = 1.0 / sqrt(1.0 + t*t);
            const double s   = t * c;
            const double tau = s / (1.0 + c);
            const double App = A[p][p], Aqq = A[q][q];
            A[p][p] = App - t * Apq;
            A[q][q] = Aqq + t * Apq;
            A[p][q] = 0.0; A[q][p] = 0.0;
            for(int r=0; r<NUMDIMS; r++) {
                if(r == p || r == q) continue;
                const double Arp = A[r][p], Arq = A[r][q];
                A[r][p] = Arp - s * (Arq + tau * Arp);
                A[r][q] = Arq + s * (Arp - tau * Arq);
                A[p][r] = A[r][p]; A[q][r] = A[r][q];
            }
            for(int r=0; r<NUMDIMS; r++) {
                const double Vrp = V[r][p], Vrq = V[r][q];
                V[r][p] = Vrp - s * (Vrq + tau * Vrp);
                V[r][q] = Vrq + s * (Vrp - tau * Vrq);
            }
        }
    }
#endif

    /* Diagonals of A post-sweep are the eigenvalues. Check whether any active
     * eigenvalue lies strictly below eigenvalue_floor: if not, S was already
     * PSD relative to the requested floor and we MUST leave it untouched
     * (do not let V*diag*V^T reconstruction roundoff drift S silently while
     * we report dT=0 / floored=false). Matches the pre-refactor 3D helper
     * semantics. */
    bool floored = false;
    double lambda[3] = {0,0,0};
    for(int a=0; a<NUMDIMS; a++) {
        if(A[a][a] < eigenvalue_floor) { floored = true; lambda[a] = eigenvalue_floor; }
        else                            { lambda[a] = A[a][a]; }
    }
    if(!floored) {
        if(dT_out) *dT_out = 0.0;
        return false;
    }

    /* Reconstruct S_active = V * diag(lambda) * V^T on the active block;
     * inactive entries stay at zero. */
    double trace_after = 0.0;
    for(int a=0; a<NUMDIMS; a++) {
        for(int b=a; b<NUMDIMS; b++) {
            double sum = 0.0;
            for(int k=0; k<NUMDIMS; k++) sum += V[a][k] * lambda[k] * V[b][k];
            S[a][b] = sum;
            if(b > a) S[b][a] = sum;
        }
        trace_after += S[a][a];
    }
    /* Floor cannot mathematically reduce trace; protect repair_dT_sum
     * diagnostics against tiny negative reconstruction roundoff (matches
     * the pre-refactor 3D helper's invariant accounting). */
    if(dT_out) {
        double dT = trace_after - trace_before;
        if(!isfinite(dT) || dT < 0.0) dT = 0.0;
        *dT_out = dT;
    }
    return true;
#else
    (void)S; (void)eigenvalue_floor;
    if(dT_out) *dT_out = 0.0;
    return false;
#endif
}
#endif


#ifdef CBE_INTEGRATOR

/* --------------------------------------------------------------------------
 * Basis-pair cost and outgoing-side assignment helpers.
 *
 * Conservation principle: the transfer step in cbe_integrator_flux_functions.h
 * applies -F to source basis and +F to target basis for every transfer.
 * Global mass/p/T conservation holds regardless of which pairing rule is
 * used; different rules produce different mixing patterns.
 *
 * Provides velocity-only cost + greedy outgoing assignment, and a trace-W2
 * cost plus adaptive-free-slot row transform, dispatched by compile-time
 * selectors. One-sided-nearest (= greedy) is the production assignment;
 * a Hungarian assignment is a possible comparator option, not implemented.
 *
 * See also the 'theta uses dp not Face_Area_Vec' note in
 * cbe_integrator_flux_functions.h.
 * -------------------------------------------------------------------------- */

/* Velocity-only squared cost between two basis components:
 *     C = sum over k of (v_a_k - v_b_k)^2
 * Basis-mass denominators floored at MIN_REAL_NUMBER to avoid NaN. */
KOKKOS_INLINE_FUNCTION
double cbe_cost_v_only(const double moments_a[CBE_INTEGRATOR_NMOMENTS],
                       const double moments_b[CBE_INTEGRATOR_NMOMENTS])
{
    double inv_a = 1.0 / DMAX(moments_a[0], MIN_REAL_NUMBER);
    double inv_b = 1.0 / DMAX(moments_b[0], MIN_REAL_NUMBER);
    double c = 0;
    /* 3-vector velocity diff; missing components (k>=NUMDIMS) are 0 via
     * cbe_basis_p_r — same dimension-agnostic algebra as hydro. */
    for(int k=0; k<3; k++) {
        double dv = cbe_basis_p_r(moments_a, k) * inv_a - cbe_basis_p_r(moments_b, k) * inv_b;
        c += dv * dv;
    }
    return c;
}


/* Outgoing-side greedy assignment: given an NBASIS x NBASIS cost matrix C,
 * fill beta_of_alpha[m] = argmin_n C[m][n] for each row m. Multiple alphas
 * may collide on the same beta (asymmetric).
 *
 * Call once with C for a-side outgoing, once with C^T for b-side outgoing. */
KOKKOS_INLINE_FUNCTION
void cbe_assign_outgoing_greedy(
    const double C[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NBASIS],
    int beta_of_alpha[CBE_INTEGRATOR_NBASIS])
{
    for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
        double best = MAX_REAL_NUMBER; int best_n = 0;
        for(int n=0; n<CBE_INTEGRATOR_NBASIS; n++) {
            if(C[m][n] < best) { best = C[m][n]; best_n = n; }
        }
        beta_of_alpha[m] = best_n;
    }
}


/* Trace-W^2 squared cost between two basis components. Rotationally
 * invariant 3D generalization of the 1D reference form
 * C = (v_a - v_b)^2 + (sqrt(S_a) - sqrt(S_b))^2.
 *
 * Velocity term: sum_k (v_a_k - v_b_k)^2, same float-op ordering as
 * cbe_cost_v_only. NMOMENTS=4 builds (no stress slots) reduce to the
 * exact velocity-only sum and are byte-identical to cbe_cost_v_only.
 *
 * Stress term [#if NMOMENTS >= 7]: trace-of-covariance form
 * tr(S) = sum_k (T_kk/m - v_k^2), with each diagonal piece floored at 0
 * before summation to absorb pre-SPD-repair numerical noise. The cost
 * contribution is (sqrt(tr(S_a)) - sqrt(tr(S_b)))^2 -- one scalar add,
 * one sqrt per side, no matrix sqrt. Rotationally invariant; collapses
 * to the 1D reference form exactly when only one stress component is
 * meaningful. A diagonal-only sum would be coordinate-frame dependent and
 * is deliberately not used.
 *
 * Full multivariate Gaussian Wasserstein-2 (needs matrix sqrt of
 * Sigma_a^{1/2} Sigma_b Sigma_a^{1/2}) is NOT this function; if a
 * comparator is ever wanted, it lands as a separate cbe_cost_w2_full
 * symbol and selector.
 *
 * The SSOT pair-matching builder routes to this via the CBE_PAIRING_COST
 * compile-time selector. */
KOKKOS_INLINE_FUNCTION
double cbe_cost_trace_w2(const double moments_a[CBE_INTEGRATOR_NMOMENTS],
                         const double moments_b[CBE_INTEGRATOR_NMOMENTS])
{
    double inv_a = 1.0 / DMAX(moments_a[0], MIN_REAL_NUMBER);
    double inv_b = 1.0 / DMAX(moments_b[0], MIN_REAL_NUMBER);
    double v_a[3], v_b[3];
    double c = 0;
    /* 3-vector velocity diff; missing components zero-padded. */
    for(int k=0; k<3; k++) {
        v_a[k] = cbe_basis_p_r(moments_a, k) * inv_a;
        v_b[k] = cbe_basis_p_r(moments_b, k) * inv_b;
        double dv = v_a[k] - v_b[k];
        c += dv * dv;
    }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    /* Active-block trace of central S = sum_{a<NUMDIMS} (T_aa/m - v_a^2),
     * each piece floored at 0 before sqrt. */
    double tr_S_a = 0, tr_S_b = 0;
    for(int a=0; a<NUMDIMS; a++) {
        double S_a_aa = cbe_basis_T_r(moments_a, a, a) * inv_a - v_a[a]*v_a[a];
        double S_b_aa = cbe_basis_T_r(moments_b, a, a) * inv_b - v_b[a]*v_b[a];
        tr_S_a += DMAX(S_a_aa, 0.0);
        tr_S_b += DMAX(S_b_aa, 0.0);
    }
    double dsqrt_S = sqrt(tr_S_a) - sqrt(tr_S_b);
    c += dsqrt_S * dsqrt_S;
#endif
    return c;
}


/* Calibrated weight for the density-continuity term in the pairing cost
 * (C = C_v + CBE_CRHO_LAMBDA * C_rho). PRODUCTION DEFAULT = 0.25. */
static const double CBE_CRHO_LAMBDA = 0.25;

/* v-slope-limiter finite opposite-sign tolerance: ignore a reconstructed sign
 * flip when |predicted| < CBE_VOPPSIGN * scale_v. PRODUCTION DEFAULT = 0.1. */
static const double CBE_VOPPSIGN = 0.1;

/* Pairing cost with phase-space-sheet (density) continuity, used when
 * stale-gradient reconstructed face densities are available:
 *   C_c   = normalized face-centric velocity distance
 *           |du|^2 / (|ua|^2 + |ub|^2 + 2|ua.ub| + trS_a + trS_b)
 *   C_rho = |rho_fa - rho_fb| / (|rho_fa| + |rho_fb| + 1e-3*rho_cell)
 *   C     = C_c + CBE_CRHO_LAMBDA * C_rho
 * rho_fa/rho_fb are the stale-gradient reconstructed face densities of the two
 * bases; vface is the full face-velocity vector (= v_F_normal * A_hat), subtracted
 * per-component so the face-centric velocities are frame-consistent for an
 * arbitrarily-oriented face normal in 3D (reduces to the 1D k=0 subtraction). */
KOKKOS_INLINE_FUNCTION
double cbe_cost_cc_crho(const double Qa[CBE_INTEGRATOR_NMOMENTS],
                        const double Qb[CBE_INTEGRATOR_NMOMENTS],
                        double rho_fa, double rho_fb, const double vface[3], double rho_cell)
{
    /* Both denominators below are regularized by MIN_REAL_NUMBER, which is a
     * pure divide-by-zero guard and carries no physical scale. A fixed additive
     * constant cannot serve here: the two denominators have DIFFERENT dimensions
     * (velocity^2 for C_c, density for C_rho), so one number cannot be small
     * against both, and whichever scale it is tuned to, rescaling the unit
     * changes the cost. Neither denominator needs a physical floor -- each
     * numerator vanishes with its own denominator (du2 <= (|ua|+|ub|)^2;
     * |rho_fa-rho_fb| <= |rho_fa|+|rho_fb|), so the guarded expression is 0/0
     * only in the exactly-degenerate case and evaluates to 0 there. C_rho keeps
     * its genuine relative regularizer 1e-3*rho_cell. */
    double inv_a = 1.0 / DMAX(Qa[0], MIN_REAL_NUMBER);
    double inv_b = 1.0 / DMAX(Qb[0], MIN_REAL_NUMBER);
    double du2 = 0, na = 0, nb = 0, dot = 0, trSa = 0, trSb = 0;
    for(int k=0; k<NUMDIMS; k++) {
        double va = cbe_basis_p_r(Qa, k) * inv_a;
        double vb = cbe_basis_p_r(Qb, k) * inv_b;
        double ua = va - vface[k];   /* face-centric: v relative to face motion along A_hat */
        double ub = vb - vface[k];
        du2 += (ua-ub)*(ua-ub); na += ua*ua; nb += ub*ub; dot += ua*ub;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        double Sa = cbe_basis_T_r(Qa, k, k) * inv_a - va*va;
        double Sb = cbe_basis_T_r(Qb, k, k) * inv_b - vb*vb;
        trSa += DMAX(Sa, 0.0); trSb += DMAX(Sb, 0.0);
#endif
    }
    double Cc   = du2 / (na + nb + 2.0*fabs(dot) + trSa + trSb + MIN_REAL_NUMBER);
    double Crho = fabs(rho_fa - rho_fb) / (fabs(rho_fa) + fabs(rho_fb) + 1.0e-3*rho_cell + MIN_REAL_NUMBER);
    return Cc + CBE_CRHO_LAMBDA * Crho;
}


/* Calibrated constants for the two-cost free-slot gate (cbe_apply_fs_gate).
 * Method constants (not runtime or compile-time configuration options). */
static const double CBE_FSGATE_TAU    = 0.04;    /* free-slot-open threshold on the continuation cost */
static const double CBE_FSGATE_DELTA  = 0.02;    /* smoothstep half-width */
static const double CBE_FSGATE_FFF    = 1.0e-2;  /* reliability mass scale f_FS (q_b = m_b/(m_b+f_FS*m_cell)) */
/* Continuation-velocity floor, as a DIMENSIONLESS fraction of a locally-derived
 * velocity scale: v_floor = CBE_FSGATE_VFLOOR_FRAC * V_local, where V_local is the
 * mass-weighted rms spread of the basis mean velocities about the pooled bulk
 * velocity of the two cells (see cbe_apply_fs_gate). Keying the floor to the
 * stream separation the pair of cells actually carries, rather than to a fixed
 * velocity in code units, makes the gate invariant under a rescaling of the
 * problem's velocity unit -- as every other constant in the pairing cost already
 * is. Same construction as the BJ limiter's scale-aware tolerance
 * (CBE_VOPPSIGN * scale_v_k). */
static const double CBE_FSGATE_VFLOOR_FRAC = 0.02;

/* Two-cost free-slot gate on an NBASIS x NBASIS cost matrix. Replaces the old
 * always-on psi reweight: the free-slot fallback fires ONLY when a source has
 * no reliable existing CONTINUATION. Density-amplitude continuity (C_rho, the
 * C2 route cost) must NOT decide that — it is reconstruction/limiter-sensitive
 * (a same-stream like-match across a steep density edge has a small velocity
 * gap but large C_rho). So the trigger is velocity-based and independent of the
 * route cost:
 *   trigger (continuation): C_cont(a,b) = |dv|^2 / (|dv|^2 + S_src + v_floor^2)
 *     velocity-center gap normalized by SOURCE dispersion S_src (a broad target
 *     cannot absorb a cold source's gap; same-center heating stays closed).
 *     v_floor = CBE_FSGATE_VFLOOR_FRAC * V_local regularizes the cold-source
 *     limit S_src -> 0, where the dispersion alone sets no scale. V_local^2 is
 *     the mass-weighted variance of the basis mean velocities about the pooled
 *     bulk velocity, taken over BOTH cells: dv is an inter-cell quantity, so a
 *     cell that is internally single-stream must still be scaled against the
 *     velocity gap to its neighbour. Both Q_src and Q_tgt reach here in the same
 *     absolute frame (cbe_build_flux_frame_Q_from_stored_moments bakes each
 *     cell's Vel in at every call site), so the pooled moment is well-defined.
 *     Empty rows carry zero mass and so do not set the scale. In the fully
 *     degenerate case -- every populated row in both cells at one velocity --
 *     V_local and every dv vanish together, and C_cont -> 0 closes the gate,
 *     which is the correct answer: a perfect continuation exists.
 *   reliability:  q_b = m_b / (m_b + f_FS * m_cell_tgt)   (no hard empty cutoff)
 *   C_fit(a) = min_b [ C_cont(a,b) + eta*(1-q_b) ];  eta = tau
 *   w_FS(a)  = smoothstep( (C_fit - (tau-Delta)) / (2 Delta) )  in [0,1]
 *   C_eff(a,b) = C_route(a,b) * ( 1 - w_FS(a)*(1 - psi_ab) ),  psi=m_b/(m_a+m_b)
 * C (= C_route) is the caller's cost matrix (C_c + lambda*C_rho); routing keeps
 * C_rho. Reliability/routing masses are the cell-centered row masses Q[m][0]
 * (src_masses/target_masses), NOT reconstructed face densities. Direction is
 * asymmetric: caller passes (Q_src,Q_tgt,src_masses,target_masses) for the build
 * direction. The old always-on psi reweight is the w_FS==1 limit of this gate.
 *
 * fired_count_inout is NULLABLE; when non-null, increments once per source row
 * that opens the gate (w_FS > 0.5). Bounded by 2*NBASIS per face. */
KOKKOS_INLINE_FUNCTION
void cbe_apply_fs_gate(
    double C[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NBASIS],
    const double Q_src[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double Q_tgt[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double src_masses[CBE_INTEGRATOR_NBASIS],
    const double target_masses[CBE_INTEGRATOR_NBASIS],
    int *fired_count_inout)
{
    const int N = CBE_INTEGRATOR_NBASIS;
    const double tau = CBE_FSGATE_TAU, Delta = DMAX(CBE_FSGATE_DELTA, 1.0e-6);
    const double eta = CBE_FSGATE_TAU, f_FS = CBE_FSGATE_FFF;
    double m_cell_tgt = 0.0, sum_rho = 0.0;
    for(int j=0; j<N; j++) { m_cell_tgt += target_masses[j]; sum_rho += src_masses[j] + target_masses[j]; }
    const double eps_rho = 1.0e-8 * sum_rho;

    /* Local velocity scale V_local for the continuation floor: mass-weighted rms
     * spread of the basis mean velocities about the pooled bulk velocity, over
     * both cells. Two passes over the 2N rows; row mass is the same cell-centered
     * Q[m][0] the caller passed as src_masses/target_masses. */
    double vfloor2 = 0.0;
    {
        const double M_pool = DMAX(sum_rho, MIN_REAL_NUMBER);
        double Vbar[3] = {0,0,0};
        for(int r=0; r<2*N; r++) {
            const double *Qr = (r < N) ? Q_src[r] : Q_tgt[r-N];
            const double mr  = DMAX((r < N) ? src_masses[r] : target_masses[r-N], 0.0);
            const double inv = 1.0 / DMAX(Qr[0], MIN_REAL_NUMBER);
            for(int k=0; k<NUMDIMS; k++) Vbar[k] += mr * cbe_basis_p_r(Qr, k) * inv;
        }
        for(int k=0; k<NUMDIMS; k++) Vbar[k] /= M_pool;
        double V2_local = 0.0;
        for(int r=0; r<2*N; r++) {
            const double *Qr = (r < N) ? Q_src[r] : Q_tgt[r-N];
            const double mr  = DMAX((r < N) ? src_masses[r] : target_masses[r-N], 0.0);
            const double inv = 1.0 / DMAX(Qr[0], MIN_REAL_NUMBER);
            double d2 = 0.0;
            for(int k=0; k<NUMDIMS; k++) { double d = cbe_basis_p_r(Qr, k)*inv - Vbar[k]; d2 += d*d; }
            V2_local += mr * d2;
        }
        V2_local /= M_pool;
        vfloor2 = CBE_FSGATE_VFLOOR_FRAC * CBE_FSGATE_VFLOOR_FRAC * DMAX(V2_local, 0.0);
    }
    double q[CBE_INTEGRATOR_NBASIS];
    for(int j=0; j<N; j++) q[j] = target_masses[j] / (target_masses[j] + f_FS*m_cell_tgt + MIN_REAL_NUMBER);
    for(int m=0; m<N; m++) {
        if(src_masses[m] <= eps_rho) continue;   /* floor-mass source: skip (orthogonal) */
        /* source basis velocity + dispersion (absolute frame), for C_cont */
        double inv_a = 1.0 / DMAX(Q_src[m][0], MIN_REAL_NUMBER);
        double va[3] = {0,0,0}, S_src = 0.0;
        for(int k=0; k<NUMDIMS; k++) va[k] = cbe_basis_p_r(Q_src[m], k) * inv_a;
        for(int k=0; k<NUMDIMS; k++) { double s = cbe_basis_T_r(Q_src[m], k, k) * inv_a - va[k]*va[k]; S_src += DMAX(s, 0.0); }
        /* C_fit = min over reliable targets of C_cont + eta*(1-q) */
        double Cfit = 1.0e30;
        for(int b=0; b<N; b++) {
            double inv_b = 1.0 / DMAX(Q_tgt[b][0], MIN_REAL_NUMBER), dv2 = 0.0;
            for(int k=0; k<NUMDIMS; k++) { double vb = cbe_basis_p_r(Q_tgt[b], k) * inv_b; dv2 += (va[k]-vb)*(va[k]-vb); }
            /* MIN_REAL_NUMBER guards the fully-degenerate case only: with the
             * floor now scaled by V_local, all three denominator terms can
             * vanish together (one velocity, no dispersion), and dv2 vanishes
             * with them -- so this evaluates C_cont = 0, gate closed. */
            double Ccont = dv2 / (dv2 + S_src + vfloor2 + MIN_REAL_NUMBER);
            double c = Ccont + eta*(1.0 - q[b]);
            if(c < Cfit) Cfit = c;
        }
        double tt = (Cfit - (tau - Delta)) / (2.0*Delta);
        tt = (tt < 0.0) ? 0.0 : ((tt > 1.0) ? 1.0 : tt);
        double wFS = tt*tt*(3.0 - 2.0*tt);
        for(int p=0; p<N; p++) {
            double psi = target_masses[p] / (src_masses[m] + target_masses[p] + eps_rho);
            C[m][p] *= (1.0 - wFS*(1.0 - psi));
        }
        if(fired_count_inout && wFS > 0.5) { (*fired_count_inout)++; }
    }
}


/* SSOT pair-matching builder. Single
 * canonical helper used by the flux body, gradient LSQ accumulator, and
 * BJ limiter to produce basis-pair matchings. Replaces three open-coded
 * (cost-matrix build + greedy assignment) patterns; all three call sites
 * reduce to one call to this helper.
 *
 * Cost: dispatched on the compile-time selector CBE_PAIRING_COST.
 *
 * Free-slot: applied iff the compile-time selector
 * CBE_PAIRING_USE_FREE_SLOT is 1. Direction is
 * asymmetric: the a->b pass uses source = Q_a masses, target = Q_b
 * masses; the b->a pass (only built when alpha_of_beta_for_b is
 * non-null) builds a fresh cost matrix and uses source = Q_b masses,
 * target = Q_a masses. The b->a matrix is NOT the transpose of the a->b
 * matrix when free-slot is on, because free-slot transforms each row
 * with the directional mass ratio.
 *
 * Assignment: cbe_assign_outgoing_greedy (one-sided-nearest per side).
 *
 * Outputs:
 *   beta_of_alpha_for_a[m] = b-side basis matched as a's outgoing target
 *                            for a-side basis m. ALWAYS written.
 *   alpha_of_beta_for_b[n] = a-side basis matched as b's outgoing target
 *                            for b-side basis n. Optional -- pass NULL if
 *                            only one direction is needed (gradient/BJ
 *                            limiter call sites). Flux passes a real
 *                            pointer (both directions needed).
 *   free_slot_fired_count_inout: nullable counter (see
 *                            cbe_apply_fs_gate). Flux call
 *                            site passes &out.cbe_pairing_free_
 *                            slot_count; gradient and BJ-limiter call
 *                            sites pass NULL since pre-pass matching is
 *                            not a flux-pairing decision.
 *   rho_face_a/rho_face_b:   per-basis stale-gradient reconstructed face
 *                            densities (NBASIS each). When both non-NULL, the
 *                            cost includes the density-continuity term
 *                            (cbe_cost_cc_crho); NULL on either side falls back
 *                            to the velocity/trace cost. vface is the full
 *                            face-velocity vector (= v_F_normal * A_hat) used by
 *                            the face-centric velocity cost. */
KOKKOS_INLINE_FUNCTION
void cbe_build_pair_matching(
    const double Q_a[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double Q_b[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    int beta_of_alpha_for_a[CBE_INTEGRATOR_NBASIS],
    int alpha_of_beta_for_b[CBE_INTEGRATOR_NBASIS],
    int *free_slot_fired_count_inout,
    const double *rho_face_a, const double *rho_face_b, const double vface[3])
{
    const int N = CBE_INTEGRATOR_NBASIS;
    /* Density-continuity cost is used iff the caller supplied stale-gradient
     * reconstructed face densities for both sides (gradient/flux passes with
     * gradients on). Otherwise fall back to the velocity/trace cost. */
    const bool use_crho = (rho_face_a != (const double*)0 && rho_face_b != (const double*)0);
    double rho_cell_a = 0.0, rho_cell_b = 0.0;
    if(use_crho) { for(int m=0;m<N;m++){ rho_cell_a += Q_a[m][0]; rho_cell_b += Q_b[m][0]; } }

    /* a->b cost matrix. */
    double C_ab[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NBASIS];
    for(int m=0; m<N; m++) {
        for(int n=0; n<N; n++) {
            if(use_crho) { C_ab[m][n] = cbe_cost_cc_crho(Q_a[m], Q_b[n], rho_face_a[m], rho_face_b[n], vface, rho_cell_a); }
#if (CBE_PAIRING_COST == CBE_COST_TRACE_W2)
            else C_ab[m][n] = cbe_cost_trace_w2(Q_a[m], Q_b[n]);
#else
            else C_ab[m][n] = cbe_cost_v_only(Q_a[m], Q_b[n]);
#endif
        }
    }
#if CBE_PAIRING_USE_FREE_SLOT
    {
        double src_masses[CBE_INTEGRATOR_NBASIS];
        double tgt_masses[CBE_INTEGRATOR_NBASIS];
        for(int m=0; m<N; m++) {
            src_masses[m] = Q_a[m][0];
            tgt_masses[m] = Q_b[m][0];
        }
        cbe_apply_fs_gate(C_ab, Q_a, Q_b, src_masses, tgt_masses,
                          free_slot_fired_count_inout);
    }
#else
    (void)free_slot_fired_count_inout;
#endif
    cbe_assign_outgoing_greedy(C_ab, beta_of_alpha_for_a);

    if(alpha_of_beta_for_b) {
        /* b->a cost matrix (separately built; not a transpose of C_ab
         * because free-slot is direction-asymmetric). */
        double C_ba[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NBASIS];
        for(int m=0; m<N; m++) {
            for(int n=0; n<N; n++) {
                if(use_crho) { C_ba[m][n] = cbe_cost_cc_crho(Q_b[m], Q_a[n], rho_face_b[m], rho_face_a[n], vface, rho_cell_b); }
#if (CBE_PAIRING_COST == CBE_COST_TRACE_W2)
                else C_ba[m][n] = cbe_cost_trace_w2(Q_b[m], Q_a[n]);
#else
                else C_ba[m][n] = cbe_cost_v_only(Q_b[m], Q_a[n]);
#endif
            }
        }
#if CBE_PAIRING_USE_FREE_SLOT
        {
            double src_masses[CBE_INTEGRATOR_NBASIS];
            double tgt_masses[CBE_INTEGRATOR_NBASIS];
            for(int m=0; m<N; m++) {
                src_masses[m] = Q_b[m][0];
                tgt_masses[m] = Q_a[m][0];
            }
            cbe_apply_fs_gate(C_ba, Q_b, Q_a, src_masses, tgt_masses,
                              free_slot_fired_count_inout);
        }
#endif
        cbe_assign_outgoing_greedy(C_ba, alpha_of_beta_for_b);
    }
}


/* --------------------------------------------------------------------------
 * SSOT scalar c_x = sqrt(gamma_e * n_hat . S . n_hat), gamma_e = 3.
 * Computed from one basis row of stored raw second moments R = <v(x)v>;
 * central stress is S = R - v(x)v. NMOMENTS=4 (no stress stored) returns 0;
 * NMOMENTS>=7 uses diagonal-only (S_kk = R_kk - v_k^2) since off-diagonal
 * raw moments are absent in 7-moment builds; NMOMENTS>=10 uses the full
 * 3D tensor. Negative n_hat.S.n_hat (should not occur post cbe_clamp_face_Q
 * SPD enforcement; defensive) returns 0.
 *
 * Used by:
 *   - cbe_face_K_and_vn_from_Q -- fills per-basis c_x array for downstream
 *     vacuum-flux residual evaluation and bracket-pad sizing.
 *   - cbe_flux_tophat_vacuum -- gets its scalar c_x from here while still
 *     building the per-basis stress contraction S_n[] locally for the
 *     full-tensor flux.
 * -------------------------------------------------------------------------- */
/* Effective adiabatic index for the CBE signal/sound speed: c = sqrt(gamma_eff * n.S.n).
 * gamma_eff = 3 is the 10-moment (Maxwellian) closure value; SINGLE definition point,
 * shared by the per-basis face signal speed (c_x, below) and the fluid HLLC collisional
 * sound speed (cbe_flux_collisional_hllc). Users may retune here. */
static const double CBE_SIGNAL_GAMMA_EFF = 3.0;

KOKKOS_INLINE_FUNCTION
double cbe_face_normal_stress_speed_from_Qrow(
    const double moments[CBE_INTEGRATOR_NMOMENTS],
    const double n_hat[3])
{
    if(!(moments[0] > MIN_REAL_NUMBER)) return 0;
    const double inv_rho = 1.0 / moments[0];
    /* Always-3-vector velocity. Missing components zero-padded by the helper.
     * Stress block runs only under SECONDMOMENT and loops over the active
     * NUMDIMS x NUMDIMS block via cbe_basis_T_load_active_scaled. The
     * layout is dim-agnostic; the SECONDMOMENT D!=3 fence
     * was lifted and the frame-conversion bookkeeping was
     * moved into the drift-kick absolute round-trip via
     * cbe_relative_row_to_absolute / cbe_absolute_to_relative_row. */
    double v[3]; cbe_basis_v_load_3(moments, v);
    double nSn = 0;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    {
        /* R = T / m on the active block; inactive entries stay at zero so the
         * n^T S n contraction below restricts to active dimensions of n_hat. */
        double R[3][3];
        cbe_basis_T_load_active_scaled(moments, inv_rho, R);
        for(int k=0; k<NUMDIMS; k++) {
            double s = 0;
            for(int l=0; l<NUMDIMS; l++) s += (R[k][l] - v[k]*v[l]) * n_hat[l];
            nSn += n_hat[k] * s;
        }
    }
#endif
    return (nSn > 0) ? sqrt(CBE_SIGNAL_GAMMA_EFF * nSn) : 0;
}


/* --------------------------------------------------------------------------
 * SSOT one-sided vacuum mass-flux density per unit face area, per basis. Branches
 * on the source-side outward normal velocity u_out exactly as the full
 * flux solver (cbe_flux_tophat_vacuum) does for its mass slot; extracting it
 * lets the v_F root-find residual and the deposited flux use bit-identical
 * branching, which is the requirement of the strict-root-found policy
 * (basis-summed F_m at v_F == 0 implies cell-summed mass conservation).
 *
 * F_m branching is scheme-dependent (set below by CBE_INTEGRATOR_RP_GAUSSIAN);
 * for the compact-support top-hat #else scheme:
 *   u_out >=  c_x       -> rho * u_out              (cold F0 supersonic)
 *  -c_x <  u_out <  c_x -> rho * (c_x/4)(1+u_out/c_x)^2   (warm fan)
 *   u_out <= -c_x       -> 0                        (vacuum, no outflow)
 *
 * Cold limit c_x -> 0: F0 branch for u_out > 0, F = 0 for u_out < 0
 * (recovers cold-F0 cleanly when no stress is stored or n.S.n is zero).
 * -------------------------------------------------------------------------- */
/* Per-basis face-flux scalars (SSOT): mass flux density per unit area F_m, and
 * the pressure (pstress) and intrinsic-stress (gomega) flux COEFFICIENTS, all
 * for the SELECTED flux scheme. The single compile switch CBE_INTEGRATOR_RP_GAUSSIAN
 * flips compact-support top-hat <-> exact one-sided Gaussian (kinetic) HERE, in
 * this one place; both
 * consumers -- the v_F-normal root-find residual and the full flux solver
 * cbe_flux_tophat_vacuum -- call this, so they can never disagree on the scheme.
 * Inputs are complete in (rho, u_out, c_x): sigma_x = c_x/sqrt(3), with
 * sigma_x^2 = n.S.n the normal central stress; u_out is the source-side outward
 * normal velocity relative to the (moving) face.
 *
 * The full flux is assembled identically for either scheme (cbe_flux_tophat_vacuum):
 *   F_mass   = F_m
 *   F_mom_k  = v_k F_m + pstress * S_n_k
 *   F_T_kl   = R_kl F_m + pstress (v_k S_n_l + v_l S_n_k) + gomega S_n_k S_n_l
 * Top-hat:  compact-support uniform-in-v approximate solver, warm fan
 *           -c_x < u_out < c_x; pstress and gomega both nonzero (see #else).
 * Gaussian: with xi = u_out/sigma_x, Phi = normal CDF, phi = normal PDF,
 *   F_m = rho (sigma_x phi + u_out Phi);  pstress = rho Phi;  gomega = rho phi / sigma_x.
 * F_m is monotone non-decreasing in u_out for both (top-hat warm-fan slope
 * (rho/2)(1+u') >= 0; Gaussian slope rho*Phi in [0,rho]) -> the residual stays
 * monotone, bisection converges.
 *
 * This is a flux-SCHEME / coefficient model: a new scheme (with different
 * pstress/gomega coefficients) is added as another branch here, and both
 * consumers inherit it for free. */
struct CbeFaceFluxScalars { double F_m; double pstress; double gomega; };

KOKKOS_INLINE_FUNCTION
CbeFaceFluxScalars cbe_face_flux_scalars(double rho, double u_out, double c_x)
{
    CbeFaceFluxScalars s;
#if defined(CBE_INTEGRATOR_RP_GAUSSIAN)
    if(c_x <= MIN_REAL_NUMBER) {                       /* cold / zero-stress: F0 free-stream */
        s.F_m = (u_out > 0.0) ? rho * u_out : 0.0; s.pstress = 0.0; s.gomega = 0.0;
        return s;
    }
    const double sigma_x = c_x / 1.7320508075688772;   /* sqrt(3) */
    const double xi      = u_out / sigma_x;
    const double Phi     = 0.5 * erfc(-xi / 1.4142135623730951);    /* normal CDF, /sqrt(2) */
    const double phi     = exp(-0.5*xi*xi) / 2.5066282746310002;    /* normal PDF, /sqrt(2pi) */
    s.F_m     = rho * (sigma_x * phi + u_out * Phi);
    s.pstress = rho * Phi;
    s.gomega  = rho * phi / sigma_x;
#else
    /* Compact-support "top hat" (uniform-in-v) one-sided kinetic flux: the
     * approximate solver derived from a constant-density velocity DF with the
     * same (rho, <v>, S), exact for the 1D collisionless Riemann problem. The
     * edge of support is c_x = sqrt(3) sigma_x, so the warm fan spans
     * -c_x < u_out < c_x (participation boundary -c_x). With u' = u_out/c_x:
     *   w'     = (c_x/4)(1+u')^2                 -> F_m     = rho w'
     *   zeta'  = ((2-u')/4)(1+u')^2              -> pstress = rho zeta'
     *   omega' = (sqrt(3)/8)(1-u'^2)^2           -> gomega  = rho omega'/sigma_x
     *                                                       = (3/8) rho (1-u'^2)^2 / c_x
     * Matches the exact Gaussian in pressure at u_out=0 (zeta=1/2) and to ~8% in
     * mass flux; the conduction term is ~2x smaller at its peak and truncates
     * exactly to 0 for |u_out| >= c_x, in contrast to the Gaussian's slow tail. */
    if(c_x <= MIN_REAL_NUMBER) {                       /* cold / zero-stress: F0 free-stream */
        s.F_m = (u_out > 0.0) ? rho * u_out : 0.0; s.pstress = 0.0; s.gomega = 0.0;
        return s;
    }
    if(u_out >= c_x)      { s.F_m = rho * u_out; s.pstress = rho; s.gomega = 0.0; }   /* cold F0 */
    else if(u_out > -c_x) {                            /* warm fan */
        const double up   = u_out / c_x;
        const double onep = 1.0 + up;
        const double omu2 = 1.0 - up*up;
        s.F_m     = 0.25  * rho * c_x * onep * onep;
        s.pstress = 0.25  * rho * (2.0 - up) * onep * onep;
        s.gomega  = 0.375 * rho * (omu2 * omu2) / c_x;
    }
    else                  { s.F_m = 0.0; s.pstress = 0.0; s.gomega = 0.0; }           /* vacuum */
#endif
    return s;
}


#if defined(CBE_INTEGRATOR_COLLISIONS)
/* --------------------------------------------------------------------------
 * Fluid HLLC collisional flux for the collisional Riemann term (method paper,
 * "How Collisional Can the Method Go?"). This is the strongly-collisional
 * counterpart of the collisionless one-sided vacuum flux: in the collisional
 * limit the operator-split vacuum RP breaks down (mean free path << face size),
 * and the per-face flux should instead approach the fluid HLLC solution. The
 * two are blended by the mean-free-path weight w_c = exp(-Delta_x rho sigma/mu)
 * at the call site (cbe_integrator_flux_compute_pair): the effective flux is
 *   G_eff_alpha = w_c G_vacuum_alpha + (1-w_c)(m_alpha/m_cell) G_HLLC_collisional
 * so w_c -> 1 (collisionless) recovers the pure vacuum flux exactly, and w_c ->
 * 0 (fluid limit) recovers an MFM-like state: mass flux -> 0 (F_HLLC has zero
 * mass flux, consistent with the CBE volume/face partition), with pressure and
 * PdV work exchanged across the face.
 *
 * All inputs are CELL-bulk (summed over bases), face-normal, face-relative:
 *   rho_K   physical cell density (m_cell/V_cell)
 *   vn_K    (<v_K> . n_hat) - v_F_n     face-relative normal cell velocity
 *   P_K     rho_K Tr(S^com_K)/D         isotropic cell pressure
 *   c_K     sqrt(CBE_SIGNAL_GAMMA_EFF * n.S^com_K.n)   signal speed (same convention as c_x)
 * vn_K are FACE-RELATIVE (vn = <v>.n_hat - v_F_normal); v_F_normal is the face
 * velocity, needed to put the pressure-work energy flux back into the lab frame.
 * Output Gc[] is the PER-UNIT-AREA flux (mass, momentum, 2nd-moment slots) in the
 * n_hat-outward frame; the caller multiplies by |Area| and distributes by mass
 * fraction. The intermediate pressure P_star IS the non-vacuum criterion: a valid
 * (non-vacuum) star state has P_star > 0, and P_star -> 0 continuously at the
 * vacuum boundary, so gating on P_star gives a positive, continuous pressure with
 * no clamp. Returns the F1 flux if P_star > 0, else zero.
 * -------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void cbe_flux_collisional_hllc(double rho_L, double vn_L, double P_L, double c_L,
                               double rho_R, double vn_R, double P_R, double c_R,
                               double v_F_normal, const double n_hat[3],
                               double Gc[CBE_INTEGRATOR_NMOMENTS])
{
    for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) Gc[k] = 0.0;
    const double c_star = DMAX(c_L, c_R);
    const double S_L = DMIN(vn_L, vn_R) - c_star;
    const double S_R = DMAX(vn_L, vn_R) + c_star;
    const double eta_L = rho_L * (S_L - vn_L);
    const double eta_R = rho_R * (S_R - vn_R);
    const double deta  = eta_L - eta_R;
    if(!(fabs(deta) > MIN_REAL_NUMBER)) return;
    const double inv_deta = 1.0 / deta;
    const double S_star = ((P_R - P_L) + (eta_L * vn_L - eta_R * vn_R)) * inv_deta;
    const double P_star = (P_L * eta_R - P_R * eta_L + eta_L * eta_R * (vn_R - vn_L)) * (-inv_deta);
    if(!(P_star > 0.0)) return;                              /* vacuum star region -> zero flux */
    /* F1: zero mass flux; momentum = P_star n_hat; 2nd-moment (energy) is the
     * lab-frame pressure work P_star * (S_star + v_F_normal) * 2/D on the diagonal
     * (isotropic; off-diagonals zero). S_star is the face-relative contact speed;
     * adding v_F_normal converts to the lab frame the moments are stored in. */
    for(int k = 0; k < NUMDIMS; k++) cbe_basis_p_w(Gc, k, P_star * n_hat[k]);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    const double e_diag = 2.0 * P_star * (S_star + v_F_normal) / (double)NUMDIMS;
    for(int a = 0; a < NUMDIMS; a++) cbe_basis_T_w(Gc, a, a, e_diag);
#endif
}
#endif /* CBE_INTEGRATOR_COLLISIONS */


/* --------------------------------------------------------------------------
 * Per-face root-found face-normal velocity v_F_n. Returns net
 * per-unit-area face mass flux at trial v_F_n, summing the one-sided vacuum
 * mass-flux contribution (scheme per cbe_face_flux_scalars) from every basis
 * on both sides.
 *
 * Sign convention matches the deposited flux: i-side
 * basis outflow (u_out_i = v_alpha_n_i - v_F_n) carries mass in +A_hat,
 * j-side basis outflow (u_out_j = v_F_n - v_alpha_n_j) carries mass in
 * -A_hat. Net flux in +A_hat is therefore
 *
 *   r(v_F_n) = sum_m  F_m(K_i[m], u_out_i, c_x_i[m])
 *            - sum_m  F_m(K_j[m], u_out_j, c_x_j[m]).
 *
 * Monotone non-increasing in v_F_n (each per-basis F_m is monotone
 * non-decreasing in u_out: F0 slope rho, F1 slope 3 rho / 4, F=0 slope 0;
 * u_out_i decreases as v_F_n rises; u_out_j increases as v_F_n rises);
 * bisection in cbe_face_solve_v_F_normal converges to the unique zero.
 * K==0 (rho-clamped-inactive) rows contribute 0 by construction.
 * -------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
double cbe_face_mass_residual_per_unit_area(
    double v_F_n,
    const double v_alpha_n_i[CBE_INTEGRATOR_NBASIS],
    const double v_alpha_n_j[CBE_INTEGRATOR_NBASIS],
    const double K_i[CBE_INTEGRATOR_NBASIS],
    const double K_j[CBE_INTEGRATOR_NBASIS],
    const double c_x_i[CBE_INTEGRATOR_NBASIS],
    const double c_x_j[CBE_INTEGRATOR_NBASIS])
{
    double r = 0;
    for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
        const double u_out_i = v_alpha_n_i[m] - v_F_n;
        const double u_out_j = v_F_n - v_alpha_n_j[m];
        r += cbe_face_flux_scalars(K_i[m], u_out_i, c_x_i[m]).F_m;
        r -= cbe_face_flux_scalars(K_j[m], u_out_j, c_x_j[m]).F_m;
    }
    return r;
}


/* Bisection root-find for v_F_n with bracket-widen-up-to-4x and explicit
 * fallback flagged via bracket_ok_out=0. Uses bisection (not Brent) for
 * device portability. Scale-aware termination: stop when (hi-lo) drops below
 * tol_abs + tol_rel * max(|lo|,|hi|), so high-velocity (cosmological)
 * brackets reach machine-eps relative precision rather than absolute
 * 1e-12. Iteration cap = 60.
 * NO silent midpoint -- on bracket failure the caller's analytic
 * fallback v_F is used and bracket_fail_count is incremented. */
KOKKOS_INLINE_FUNCTION
double cbe_face_solve_v_F_normal(
    const double v_alpha_n_i[CBE_INTEGRATOR_NBASIS],
    const double v_alpha_n_j[CBE_INTEGRATOR_NBASIS],
    const double K_i[CBE_INTEGRATOR_NBASIS],
    const double K_j[CBE_INTEGRATOR_NBASIS],
    const double c_x_i[CBE_INTEGRATOR_NBASIS],
    const double c_x_j[CBE_INTEGRATOR_NBASIS],
    double v_F_lo, double v_F_hi,
    double fallback_v_F_n,
    int *bracket_ok_out)
{
    double lo = v_F_lo, hi = v_F_hi;
    double f_lo = cbe_face_mass_residual_per_unit_area(lo, v_alpha_n_i, v_alpha_n_j, K_i, K_j, c_x_i, c_x_j);
    double f_hi = cbe_face_mass_residual_per_unit_area(hi, v_alpha_n_i, v_alpha_n_j, K_i, K_j, c_x_i, c_x_j);
    int bracketed = (f_lo * f_hi <= 0) ? 1 : 0;
    for(int widen = 0; widen < 4 && !bracketed; widen++) {
        double mid  = 0.5 * (lo + hi);
        double half = 4.0 * 0.5 * (hi - lo);
        lo = mid - half; hi = mid + half;
        f_lo = cbe_face_mass_residual_per_unit_area(lo, v_alpha_n_i, v_alpha_n_j, K_i, K_j, c_x_i, c_x_j);
        f_hi = cbe_face_mass_residual_per_unit_area(hi, v_alpha_n_i, v_alpha_n_j, K_i, K_j, c_x_i, c_x_j);
        if(f_lo * f_hi <= 0) { bracketed = 1; break; }
    }
    if(!bracketed) {
        *bracket_ok_out = 0;
        return fallback_v_F_n;
    }
    *bracket_ok_out = 1;
    /* tol_rel = 1e-14: the strict-root-found
     * policy means the per-face mass residual at v_F is bounded by
     * |F_m'(v_F)| * (hi - lo), and |F_m'| ~ sum of active basis K's
     * times Face_Area_Norm. For Hernquist v_F bracket scale ~200 the
     * 1e-12 rel tol left residuals at ~1e-10 (above the col-2 <= 1e-11
     * gate); 1e-14 brings them solidly below. Extra iterations are cheap
     * vs the flux loop. 60-iter cap still allows convergence: 2^60 >> any
     * physically plausible width / 1e-14. */
    const double tol_abs = 1e-12;
    const double tol_rel = 1e-14;
    for(int it=0; it<60; it++) {
        double scale = DMAX(fabs(lo), fabs(hi));
        if((hi - lo) <= tol_abs + tol_rel * scale) break;
        double mid = 0.5 * (lo + hi);
        double f_mid = cbe_face_mass_residual_per_unit_area(mid, v_alpha_n_i, v_alpha_n_j, K_i, K_j, c_x_i, c_x_j);
        if(f_lo * f_mid <= 0) { hi = mid; f_hi = f_mid; }
        else                  { lo = mid; f_lo = f_mid; }
    }
    return 0.5 * (lo + hi);
}


/* SSOT relative <-> absolute moment transforms (finite
 * absolute-update round-trip in cbe_drift_kick).
 *
 *   ABSOLUTE-frame moments of basis b:
 *     m_abs    = m_b
 *     p_abs[a] = m_b * <v_a>_b           = p_rel[a] + m_b * V[a]
 *     T_abs[a,b] = m_b * <v_a v_b>_b     = T_rel[a,b] + V[a] p_rel[b] + p_rel[a] V[b] + m_b V[a] V[b]
 *
 *   RELATIVE-frame storage matches the inverse:
 *     p_rel[a]   = p_abs[a]   - m_b * V[a]
 *     T_rel[a,b] = T_abs[a,b] - V[a] p_rel[a,b] - p_rel[a] V[b] - m_b V[a] V[b]   (uses new p_rel)
 *
 * These two helpers are the SSOT used by the cbe_drift_kick absolute
 * round-trip and by any standalone test code (unit-checked algebra). The row
 * layout is the same in both frames; only the [1..NUMDIMS] and [T_idx] slot
 * VALUES differ. Operations are NUMDIMS-aware and gate the stress block on
 * SECONDMOMENT. */
KOKKOS_INLINE_FUNCTION
void cbe_relative_row_to_absolute(
    const double row_rel[CBE_INTEGRATOR_NMOMENTS],
    const double V[3],
    double *m_abs_out,
    double p_abs_out[3],
    double T_abs_out[3][3])
{
    const double m_b = row_rel[0];
    *m_abs_out = m_b;
    /* Read RELATIVE momentum slots; populate the active dims of p_abs. */
    double p_rel_a[3] = {0,0,0};
    for(int a=0;a<NUMDIMS;a++) p_rel_a[a] = cbe_basis_p_r(row_rel, a);
    for(int a=0;a<3;a++) p_abs_out[a] = 0.0;
    for(int a=0;a<NUMDIMS;a++) p_abs_out[a] = p_rel_a[a] + m_b * V[a];
    /* T_abs = T_rel + V*p_rel + p_rel*V + m*V*V on active block; zero outside. */
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) T_abs_out[a][b] = 0.0;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0;a<NUMDIMS;a++) {
        T_abs_out[a][a] = cbe_basis_T_r(row_rel, a, a)
                        + 2.0 * V[a] * p_rel_a[a]
                        + m_b * V[a] * V[a];
        for(int b=a+1;b<NUMDIMS;b++) {
            const double T_ab = cbe_basis_T_r(row_rel, a, b)
                              + V[a] * p_rel_a[b] + V[b] * p_rel_a[a]
                              + m_b * V[a] * V[b];
            T_abs_out[a][b] = T_ab; T_abs_out[b][a] = T_ab;
        }
    }
#endif
}

KOKKOS_INLINE_FUNCTION
void cbe_absolute_to_relative_row(
    double m_abs,
    const double p_abs_in[3],
    const double T_abs_in[3][3],
    const double V[3],
    double row_rel_out[CBE_INTEGRATOR_NMOMENTS])
{
    row_rel_out[0] = m_abs;
    /* Compute NEW p_rel first; T_rel update uses NEW p_rel below. */
    double p_rel_new_a[3] = {0,0,0};
    for(int a=0;a<NUMDIMS;a++) {
        p_rel_new_a[a] = p_abs_in[a] - m_abs * V[a];
        cbe_basis_p_w(row_rel_out, a, p_rel_new_a[a]);
    }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0;a<NUMDIMS;a++) {
        const double T_rel_aa = T_abs_in[a][a]
                              - 2.0 * V[a] * p_rel_new_a[a]
                              - m_abs * V[a] * V[a];
        cbe_basis_T_w(row_rel_out, a, a, T_rel_aa);
        for(int b=a+1;b<NUMDIMS;b++) {
            const double T_rel_ab = T_abs_in[a][b]
                                  - V[a] * p_rel_new_a[b] - V[b] * p_rel_new_a[a]
                                  - m_abs * V[a] * V[b];
            cbe_basis_T_w(row_rel_out, a, b, T_rel_ab);
        }
    }
#else
    (void)T_abs_in;
#endif
}


/* --------------------------------------------------------------------------
 * SSOT conversion from stored basis-frame moments to "flux-frame" Q.
 * One helper used by BOTH the gradient pass and the per-pair reconstruction
 * so they cannot drift apart on cf_a3inv / cf_atime / NMOMENTS handling.
 *
 * Output (note: in this codebase `cf_a3inv * density_code` = physical
 * density in code units, NOT a comoving-frame density — the `a3inv` factor
 * turns the comoving-volume cell-integrated U into a physical-frame density.):
 *   Q[m][0]    = physical-frame mass density of basis m   (code units)
 *   Q[m][1..3] = physical-frame momentum density          = Q[m][0]*v_phys
 *   Q[m][4..9] = physical-frame stress density (only when NMOMENTS==10;
 *                the 3D-second-moment fence in precompiler_logic.h
 *                ensures the 3D layout [4]/[5]/[6]=diag,
 *                [7]/[8]/[9]=off-diag holds).
 *
 * "Flux-frame" = physical-frame momentum baked in (v_phys = Vel/cf_atime
 * folded into the [1..3] slots so cbe_flux_tophat_vacuum's u_out matches
 * v_phys directly). NOT a generic U/V conversion. */
KOKKOS_INLINE_FUNCTION
void cbe_build_flux_frame_Q_from_stored_moments(
    const double U[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double Vel_code[3],
    double V_i,
    double cf_a3inv,
    double cf_atime,
    double Q_out[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS])
{
    double inv_V = 1.0 / DMAX(V_i, MIN_REAL_NUMBER);
    double inv_a = 1.0 / cf_atime;
    for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
        for(int k=0; k<CBE_INTEGRATOR_NMOMENTS; k++) {
            Q_out[m][k] = U[m][k] * inv_V * cf_a3inv;
        }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        /* Relative-frame stress boost. Storage convention:
         *    p_rel = m * (v_phys - V)
         *    T_rel = m * <(v_phys - V) (v_phys - V)>
         * where V is the mesh-generating-point velocity. The flux solver
         * consumes raw absolute-frame T_abs, related by
         *    T_abs = T_rel + V (x) p_rel + p_rel (x) V + m * V (x) V
         *
         * V here is the SAME Vel_code[a]*inv_a used by the momentum boost
         * below — the FLUX-frame velocity scaling matches the existing
         * momentum-slot convention. p_rel comes from Q_out[m][1..NUMDIMS]
         * BEFORE the momentum boost runs (captured here). The stress boost
         * MUST precede the momentum boost: it reads p_rel from the still-
         * relative-frame momentum slots. */
        {
            const double rho = Q_out[m][0];
            double V_act[3]    = {0,0,0};
            double p_rel_act[3] = {0,0,0};
            for(int a=0; a<NUMDIMS; a++) {
                V_act[a]    = Vel_code[a] * inv_a;
                p_rel_act[a] = cbe_basis_p_r(Q_out[m], a);   /* pre-boost */
            }
            for(int a=0; a<NUMDIMS; a++) {
                /* Diagonal: T_abs_aa = T_rel_aa + 2*V_a*p_rel_a + rho*V_a*V_a. */
                cbe_basis_T_w(Q_out[m], a, a,
                    cbe_basis_T_r(Q_out[m], a, a)
                    + 2.0 * V_act[a] * p_rel_act[a]
                    + rho * V_act[a] * V_act[a]);
                for(int b=a+1; b<NUMDIMS; b++) {
                    /* Off-diagonal: T_abs_ab = T_rel_ab + V_a*p_rel_b
                     *                        + V_b*p_rel_a + rho*V_a*V_b. */
                    cbe_basis_T_w(Q_out[m], a, b,
                        cbe_basis_T_r(Q_out[m], a, b)
                        + V_act[a] * p_rel_act[b]
                        + V_act[b] * p_rel_act[a]
                        + rho * V_act[a] * V_act[b]);
                }
            }
        }
#endif
        /* Velocity boost on the momentum slots [1..NUMDIMS]. cbe_basis_p_a
         * is a no-op for a >= NUMDIMS so the active-only contract holds for
         * any NUMDIMS. Must run AFTER the stress boost above (which reads
         * the still-relative-frame momentum slots). */
        for(int a=0; a<NUMDIMS; a++) {
            cbe_basis_p_a(Q_out[m], a, Q_out[m][0] * Vel_code[a] * inv_a);
        }
    }
}


/* Symmetric 3x3 SSOT realizability machinery.
 *
 * The three helpers below provide the canonical PSD predicate + projection
 * for CBE basis rows, shared by the cone limiter, the face-state
 * Q-clamp, and the cell-state drift-kick repair. Storage convention: raw
 * moment slots [4..9] hold T_kl = m * R_kl = m * (S_kl + v_k v_l) (raw
 * second-moment density). Central stress S_kl = R_kl - v_k v_l is built
 * locally by the raw-row wrapper before projection.
 *
 * The raw-row wrapper does the raw <-> central conversion correctly (a
 * subtlety: passing raw R slots to a central-S API applies an effective
 * floor CBE_SPD_RELATIVE_FLOOR * trace(raw R) ~ 1e-12 * m * (trace(S) +
 * v_bulk^2), a hidden bulk-KE floor on stress). The core projection takes
 * an explicit eigenvalue_floor argument so each caller picks its own policy.
 */


/* SSOT 3x3 symmetric PSD test via all principal minors. Returns true iff
 * every principal minor of M is >= -eps_tol * scale, where scale is the
 * natural amplitude for each minor order (trace for diagonal, trace^2 for
 * 2x2, trace^3 for det). eps_tol = 0 gives strict Sylvester-style PSD;
 * eps_tol > 0 absorbs FP-scale near-zero negative roundoff. Cone limiter
 * uses 0 (strict); cell-state repair uses 1e-12. */
KOKKOS_INLINE_FUNCTION
bool cbe_sym3x3_all_principal_minors_nonneg(const double M[3][3], double eps_tol)
{
    const double trace_pos = DMAX(M[0][0] + M[1][1] + M[2][2], MIN_REAL_NUMBER);
    const double eps_diag  = eps_tol * trace_pos;
    const double eps_2x2   = eps_tol * trace_pos * trace_pos;
    const double eps_det   = eps_tol * trace_pos * trace_pos * trace_pos;
    if(M[0][0] < -eps_diag) return false;
    if(M[1][1] < -eps_diag) return false;
    if(M[2][2] < -eps_diag) return false;
    if(M[0][0]*M[1][1] - M[0][1]*M[0][1] < -eps_2x2) return false;
    if(M[0][0]*M[2][2] - M[0][2]*M[0][2] < -eps_2x2) return false;
    if(M[1][1]*M[2][2] - M[1][2]*M[1][2] < -eps_2x2) return false;
    const double det = M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[1][2])
                     - M[0][1] * (M[0][1]*M[2][2] - M[1][2]*M[0][2])
                     + M[0][2] * (M[0][1]*M[1][2] - M[1][1]*M[0][2]);
    return (det >= -eps_det);
}


/* SSOT realizability predicate for a CBE basis row in raw-slot layout.
 * A row is realizable iff m > 0 AND (NMOMENTS >= 7) M = m*T - p p^T is PSD,
 * where T_kl = U_row[4..9] are the raw second-moment density slots and
 * the PSD test uses cbe_sym3x3_all_principal_minors_nonneg with eps_tol.
 * (Note the M expression: M = m * R_kl * m - p_k p_l = T_kl * m - p_k p_l
 *  with T_kl = m * R_kl already, so M = T_kl * m - p p^T as written.)
 *
 * Consumers: cone limiter (strict, eps_tol = 0); cell-state repair
 * (FP-scale, eps_tol = CBE_REPAIR_EPS_M_REL). */
KOKKOS_INLINE_FUNCTION
bool cbe_basis_row_is_realizable(
    const double Q_row[CBE_INTEGRATOR_NMOMENTS], double eps_tol)
{
    if(!(Q_row[0] > 0) || !isfinite(Q_row[0])) return false;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    {
        /* Active-dim realizability: M = T * m - p p^T on the upper-left
         * NUMDIMS x NUMDIMS block (zero-padded outside). T loaded via the
         * stress helper so the slot layout is correct for any NUMDIMS;
         * the predicate checks principal minors only of the active block. */
        const double m = Q_row[0];
        double p[3]; cbe_basis_p_load_3(Q_row, p);
        double T[3][3]; cbe_basis_T_load_active(Q_row, T);
        double M[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for(int a=0; a<NUMDIMS; a++) {
            M[a][a] = T[a][a] * m - p[a]*p[a];
            for(int b=a+1; b<NUMDIMS; b++) {
                const double off = T[a][b] * m - p[a]*p[b];
                M[a][b] = off; M[b][a] = off;
            }
        }
        if(!cbe_symNxN_active_principal_minors_nonneg(M, eps_tol)) return false;
    }
#endif
    return true;
}


/* Realizability test on a PRIMITIVE row (ρ, v_k, S_kl). Equivalent to the
 * moment-row test cbe_basis_row_is_realizable via the identity
 *     ρ·T − p·p = ρ²·S
 * so checking ρ > floor + S PSD on primitives is exactly checking
 * ρ > 0 + (m·T − p·p) PSD on moments — no cancellation between p² and ρ·T.
 *
 * Active dim: NMOMENTS=2 (cold) tests only ρ > floor. NMOMENTS=3 (1D
 * SECONDMOMENT) tests ρ > floor AND S_xx >= 0. NMOMENTS=10 (3D) tests
 * ρ > floor AND S_3x3 PSD via the active principal-minor chain.
 *
 * Used by the primitive-gradient swap's row-scalar cone limiter. The
 * ρ-edge is strict `>` to match the existing active-density predicate
 * cbe_clamp_face_Q uses (see cbe_rho_active_floor()). Distinct from
 * cbe_basis_row_is_realizable which still operates on moment rows
 * (drift-kick cell-state repair, free-slot bookkeeping). The stress
 * block is loaded via cbe_basis_T_load_active using the same slot
 * indexing as moment rows — the primitive slot layout was designed
 * to match by construction (slot 1+NUMDIMS+... stores S_kl rather
 * than T_kl, with identical (k,l) packing). */
KOKKOS_INLINE_FUNCTION
bool cbe_primitive_row_is_realizable(
    const double prim_row[CBE_INTEGRATOR_NMOMENTS], double eps_tol)
{
    if(!(prim_row[0] > cbe_rho_active_floor()) || !isfinite(prim_row[0])) return false;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double S[3][3];
    cbe_basis_T_load_active(prim_row, S);
    if(!cbe_symNxN_active_principal_minors_nonneg(S, eps_tol)) return false;
#else
    (void)eps_tol;
#endif
    return true;
}


/* CORE: symmetric 3x3 PSD projection of a CENTRAL stress tensor. Operates
 * in-place on the 6-vector S = [Sxx, Syy, Szz, Sxy, Sxz, Syz] (NOT a raw
 * moment row — use cbe_basis_row_project_central_stress_to_PSD for that).
 * Returns true iff S was modified.
 *
 * Algorithm: 6 fixed Jacobi sweeps to diagonalize, clamp eigenvalues to
 * max(lambda_k, eigenvalue_floor), reconstruct S = V * diag(lambda) * V^T.
 *
 * eigenvalue_floor: caller-supplied. Pass 0 for true PSD projection
 * (cell-state repair); pass CBE_SPD_RELATIVE_FLOOR * trace_central_S_pos
 * for the face-clamp floor policy (computed from trace(central S), NOT
 * trace(raw R)).
 *
 * Degenerate-input fallback: fires ONLY when input is non-finite (NaN/Inf
 * in any S slot). trace_before <= 0 with finite entries is a legitimate
 * central state to project (cold flows naturally give central S near or
 * below zero) — Jacobi handles it, no isotropic-MIN_REAL injection. */
KOKKOS_INLINE_FUNCTION
bool cbe_project_central_S_to_PSD(double S[6], double eigenvalue_floor,
                                  double *dT_out)
{
    /* Layout: [0]=Sxx, [1]=Syy, [2]=Szz, [3]=Sxy, [4]=Sxz, [5]=Syz. */
    double trace_before = S[0] + S[1] + S[2];

    /* Degenerate-input fallback (non-finite only). Off-diagonal NaN check
     * is load-bearing (Jacobi propagates NaN through eigenvalues
     * silently). */
    bool all_finite = true;
    for(int q = 0; q < 6; q++) { if(!isfinite(S[q])) { all_finite = false; break; } }
    if(!all_finite || !isfinite(trace_before)) {
        S[0] = MIN_REAL_NUMBER; S[1] = MIN_REAL_NUMBER; S[2] = MIN_REAL_NUMBER;
        S[3] = 0.0; S[4] = 0.0; S[5] = 0.0;
        if(dT_out) *dT_out = 3.0 * MIN_REAL_NUMBER;
        return true;
    }

    /* Full 3x3 symmetric matrix A + identity eigenvector matrix V. */
    double A[3][3] = {{S[0], S[3], S[4]},
                      {S[3], S[1], S[5]},
                      {S[4], S[5], S[2]}};
    double V[3][3] = {{1.0, 0.0, 0.0},
                      {0.0, 1.0, 0.0},
                      {0.0, 0.0, 1.0}};

    /* 6 fixed Jacobi sweeps over the three off-diagonal positions. */
    for(int sweep = 0; sweep < 6; sweep++) {
        for(int pq = 0; pq < 3; pq++) {
            int p, q;
            if     (pq == 0) { p = 0; q = 1; }
            else if(pq == 1) { p = 0; q = 2; }
            else             { p = 1; q = 2; }
            double Apq = A[p][q];
            if(fabs(Apq) < 1.0e-300) continue;
            double theta = (A[q][q] - A[p][p]) / (2.0 * Apq);
            double t;
            if(fabs(theta) > 1.0e15) {
                t = 0.5 / theta;
            } else {
                double sgn = (theta >= 0) ? 1.0 : -1.0;
                t = sgn / (fabs(theta) + sqrt(theta*theta + 1.0));
            }
            double c   = 1.0 / sqrt(1.0 + t*t);
            double s   = t * c;
            double tau = s / (1.0 + c);
            /* Rotate A[p][p], A[q][q], A[p][q]. */
            A[p][p] -= t * Apq;
            A[q][q] += t * Apq;
            A[p][q] = 0.0; A[q][p] = 0.0;
            /* Rotate the remaining row/column r (the one not in {p,q}). */
            for(int r = 0; r < 3; r++) {
                if(r == p || r == q) continue;
                double Arp = A[r][p];
                double Arq = A[r][q];
                A[r][p] = Arp - s * (Arq + tau * Arp);
                A[r][q] = Arq + s * (Arp - tau * Arq);
                A[p][r] = A[r][p];
                A[q][r] = A[r][q];
            }
            /* Rotate eigenvector matrix V. */
            for(int r = 0; r < 3; r++) {
                double Vrp = V[r][p];
                double Vrq = V[r][q];
                V[r][p] = Vrp - s * (Vrq + tau * Vrp);
                V[r][q] = Vrq + s * (Vrp - tau * Vrq);
            }
        }
    }

    /* Eigenvalues = diagonal of converged A. Floor with caller-supplied
     * eigenvalue_floor (was: CBE_SPD_RELATIVE_FLOOR * trace_pos baked-in). */
    double lambda[3] = {A[0][0], A[1][1], A[2][2]};
    bool floored = false;
    for(int k = 0; k < 3; k++) {
        if(lambda[k] < eigenvalue_floor) { lambda[k] = eigenvalue_floor; floored = true; }
    }
    if(!floored) {
        /* Already SPD within floor — no modification, no diagnostic increment. */
        if(dT_out) *dT_out = 0.0;
        return false;
    }

    /* Reconstruct S = V * diag(lambda) * V^T (symmetric, store upper). */
    double M[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int i = 0; i < 3; i++) {
        for(int j = i; j < 3; j++) {
            double sum = 0.0;
            for(int k = 0; k < 3; k++) sum += lambda[k] * V[i][k] * V[j][k];
            M[i][j] = sum; M[j][i] = sum;
        }
    }
    S[0] = M[0][0]; S[1] = M[1][1]; S[2] = M[2][2];
    S[3] = M[0][1]; S[4] = M[0][2]; S[5] = M[1][2];

    /* dT enforcement: eigenvalue floor only adds to the trace by
     * construction (>= 0); roundoff in V*diag(lambda)*V^T may produce a
     * tiny negative residual when the floor fires on only one eigenvalue;
     * clamp to 0 to preserve the col-8 invariant. */
    double dT = (M[0][0] + M[1][1] + M[2][2]) - trace_before;
    if(!isfinite(dT) || dT < 0.0) dT = 0.0;
    if(dT_out) *dT_out = dT;
    return true;
}


/* RAW-ROW WRAPPER: takes a CBE raw moment row, converts to central S,
 * projects via cbe_project_central_S_to_PSD, rebuilds raw slots.
 *
 * U_row layout: slots [0..9] = (m, p_x, p_y, p_z, T_xx, T_yy, T_zz,
 * T_xy, T_xz, T_yz) with T_kl = m * R_kl raw second-moment density.
 * Central stress S_kl = R_kl - v_k v_l = T_kl / m - v_k v_l where v = p/m.
 * On modification, rebuilds T_kl_new = m * (S_kl_repaired + v_k v_l).
 *
 * Named explicitly around "raw moment row" to prevent passing raw slots to
 * a central-S-only API. dT_out is the raw-slot trace delta = m * (central-dT),
 * preserving the dP_sum / dT_sum accumulator semantic wired into
 * cbe_diagnostics col-8.
 *
 * SECONDMOMENT-gated; operates on the active D x D stress block via the
 * active-dimensional helpers. Inactive rows/cols stay at zero; PSD
 * projection only floors active eigenvalues. */
KOKKOS_INLINE_FUNCTION
bool cbe_basis_row_project_central_stress_to_PSD(
    double U_row[CBE_INTEGRATOR_NMOMENTS],
    double eigenvalue_floor,
    double *dT_out)
{
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    /* Mass-positivity guard: m <= 0 means v = p/m undefined; caller is
     * responsible for rho-clamping upstream (cbe_clamp_face_Q does so via
     * Q_face[m][0] <= MIN_REAL_NUMBER row-zeroing; the cell-state repair's
     * corrupt-cell guard does so at Step 0). Treat as no-op here. */
    if(!(U_row[0] > 0) || !isfinite(U_row[0])) {
        if(dT_out) *dT_out = 0;
        return false;
    }
    const double m     = U_row[0];
    const double inv_m = 1.0 / m;
    double v[3]; cbe_basis_v_load_3(U_row, v);

    /* Central stress S = T/m - v v^T on the active block. Inactive entries
     * are exactly zero on entry to the projection (cbe_basis_T_load_active_scaled
     * zeros above NUMDIMS, and the v[k]*v[l] correction is only applied
     * on the active loop). */
    double S[3][3];
    cbe_basis_T_load_active_scaled(U_row, inv_m, S);
    for(int a=0; a<NUMDIMS; a++) {
        S[a][a] -= v[a]*v[a];
        for(int b=a+1; b<NUMDIMS; b++) {
            const double off = S[a][b] - v[a]*v[b];
            S[a][b] = off; S[b][a] = off;
        }
    }

    double dT_central = 0;
    bool modified = cbe_project_central_S_active_to_PSD(S, eigenvalue_floor, &dT_central);

    if(modified) {
        /* Rebuild raw slots: T_kl = m * (S_kl + v_k v_l) on active block. */
        for(int a=0; a<NUMDIMS; a++) {
            U_row[cbe_T_idx(a, a)] = m * (S[a][a] + v[a]*v[a]);
            for(int b=a+1; b<NUMDIMS; b++)
                U_row[cbe_T_idx(a, b)] = m * (S[a][b] + v[a]*v[b]);
        }
    }
    if(dT_out) *dT_out = dT_central * m;   /* raw-slot trace delta */
    return modified;
#else
    /* SECONDMOMENT not active: nothing to project. */
    (void)U_row; (void)eigenvalue_floor;
    if(dT_out) *dT_out = 0;
    return false;
#endif
}


/* Density-only face-Q clamp + counter. For each
 * basis with Q_face[m][0] <= MIN_REAL_NUMBER, zero the ENTIRE basis row.
 * Rationale: leaving nonzero momentum/stress at zero density would crash
 * cbe_flux_tophat_vacuum's inv_rho = 1/moments[0] divide. Zeroing the row
 * marks the basis inactive at this face -- the downstream cbe_face_K_and_vn
 * helper gives it K=0, v_n=0 so it contributes nothing to the residual or
 * the flux loop.
 *
 * Also performs face-state SPD enforcement on the stress block of each
 * rho-active basis row via the raw-row wrapper
 * cbe_basis_row_project_central_stress_to_PSD, which correctly converts
 * raw <-> central. The eigenvalue floor is computed from trace(central S)
 * — NOT trace(raw R) which secretly includes the bulk-KE contribution
 * v_bulk^2 and would act as a hidden physical floor. In cold flows the
 * floor magnitude is thus much smaller (1e-12 * trace_central_S vs
 * 1e-12 * (trace_central_S + v_bulk^2) when bulk kinetic dominates).
 *
 * SECONDMOMENT-gated: the stress block runs in any dimension under
 * SECONDMOMENT; trace and projection helpers loop over active NUMDIMS so
 * 1D/2D warm physics participates correctly. */
KOKKOS_INLINE_FUNCTION
void cbe_clamp_face_Q(
    double Qface[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    long long *rho_clamp_count,
    long long *S_clamp_count)
{
#if defined(CBE_INTEGRATOR_STRICT_FACE_SPD_GUARD) && defined(CBE_INTEGRATOR_SECONDMOMENT)
    /* Strict post-clamp face-state realizability guard. Default OFF;
     * enabled by adding CBE_INTEGRATOR_STRICT_FACE_SPD_GUARD to Config.sh.
     * Only meaningful under SECONDMOMENT -- without stress storage,
     * realizability reduces to mass-positivity which the rho>MIN_REAL_NUMBER
     * branch already enforces. Surfaces post-clamp non-realizability that
     * survives the rho-zero branch + the eigenvalue-floor projection --
     * with the cone limiter upstream this should NEVER fire at
     * non-roundoff amplitude. printf-only + once-per-call latch (the print
     * is device-safe; fflush/endrun are intentionally NOT used here, they
     * are host-only functions). */
    bool printed_failure = false;
#endif
    for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
        if(Qface[m][0] <= MIN_REAL_NUMBER) {
            for(int k=0; k<CBE_INTEGRATOR_NMOMENTS; k++) Qface[m][k] = 0.0;
            if(rho_clamp_count) (*rho_clamp_count)++;
        }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        else {
            /* Eigenvalue floor from active central-S trace. Active block
             * only -- the prior 3D-hardcoded sum [4]+[5]+[6] silently set
             * a floor proportional to inactive (y/z) diagonals in low-D.
             * Mass-positivity guarded above (Qface[m][0] > MIN_REAL_NUMBER). */
            const double m_face = Qface[m][0];
            const double inv_m  = 1.0 / m_face;
            double v[3]; cbe_basis_v_load_3(Qface[m], v);
            double v_dot_v = 0.0;
            for(int a=0; a<NUMDIMS; a++) v_dot_v += v[a]*v[a];
            const double trace_R_active  = cbe_basis_T_trace_active(Qface[m], inv_m);
            const double trace_S_central = trace_R_active - v_dot_v;
            const double trace_S_pos     = DMAX(trace_S_central, MIN_REAL_NUMBER);
            const double eigenvalue_floor = CBE_SPD_RELATIVE_FLOOR * trace_S_pos;
            double dT_dump = 0.0;
            if(cbe_basis_row_project_central_stress_to_PSD(Qface[m], eigenvalue_floor,
                                                            &dT_dump)
               && S_clamp_count) {
                (*S_clamp_count)++;
            }
        }
#endif
#if defined(CBE_INTEGRATOR_STRICT_FACE_SPD_GUARD) && defined(CBE_INTEGRATOR_SECONDMOMENT)
        /* Two-tier tolerant check on rho-active rows (rho-clamped rows were
         * zeroed above and trivially pass). Strict eps=0 then tolerant eps=
         * CBE_REPAIR_EPS_M_REL. Only fail-of-both is reported (filters FP
         * eigenvalue noise). */
        if(!printed_failure && Qface[m][0] > MIN_REAL_NUMBER) {
            if(!cbe_basis_row_is_realizable(Qface[m], /*eps_tol=*/ 0.0)
               && !cbe_basis_row_is_realizable(Qface[m], CBE_REPAIR_EPS_M_REL)) {
                const double m_face_g = Qface[m][0];
                const double inv_m_g  = 1.0 / m_face_g;
                double v_g[3]; cbe_basis_v_load_3(Qface[m], v_g);
                double v_dot_v_g = 0.0;
                for(int a=0; a<NUMDIMS; a++) v_dot_v_g += v_g[a]*v_g[a];
                const double trR_g  = cbe_basis_T_trace_active(Qface[m], inv_m_g);
                const double trSc_g = trR_g - v_dot_v_g;
                /* Device-safe printf; no fflush/endrun (host-only fns) here. */
                printf("[CBE strict face guard] post-clamp non-realizable: "
                       "basis=%d m=%.3e trace_S_central=%.3e v.v=%.3e\n",
                       m, m_face_g, trSc_g, v_dot_v_g);
                printed_failure = true;
            }
        }
#endif
    }
}


/* Guarded face-Q -> (K[m], v_alpha_n[m]) construction for the root-find +
 * theta gate. On rows with Qface[m][0] > eps, this is
 * the standard v = momentum/density projection onto Ahat; on clamped-inactive
 * rows (Qface[m][0]==0 after cbe_clamp_face_Q), it returns K=0, v_n=0 so
 * the residual function and the flux loop both naturally skip the basis. */
KOKKOS_INLINE_FUNCTION
void cbe_face_K_and_vn_from_Q(
    const double Qface[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double Ahat[3],
    double K[CBE_INTEGRATOR_NBASIS],
    double v_alpha_n[CBE_INTEGRATOR_NBASIS],
    double c_x[CBE_INTEGRATOR_NBASIS])
{
    /* Per-basis face state for the vacuum mass-flux residual + flux
     * loop. K = rho (basis density on the face), v_alpha_n = v . Ahat
     * (normal velocity), c_x = sqrt(gamma_e * Ahat.S.Ahat) (normal
     * stress speed, gamma_e = 3). Inactive rows (Qface[m][0] <=
     * MIN_REAL_NUMBER post cbe_clamp_face_Q) zero all three outputs so
     * downstream consumers naturally skip them. c_x is computed via the
     * SSOT helper cbe_face_normal_stress_speed_from_Qrow so the wave-speed
     * definition matches what cbe_flux_tophat_vacuum and the residual use. */
    for(int m=0; m<CBE_INTEGRATOR_NBASIS; m++) {
        if(Qface[m][0] > MIN_REAL_NUMBER) {
            double v_basis[3]; cbe_basis_v_load_3(Qface[m], v_basis);
            v_alpha_n[m] = v_basis[0]*Ahat[0] + v_basis[1]*Ahat[1] + v_basis[2]*Ahat[2];
            K[m]         = Qface[m][0];
            c_x[m]       = cbe_face_normal_stress_speed_from_Qrow(Qface[m], Ahat);
        } else {
            v_alpha_n[m] = 0.0;
            K[m]         = 0.0;
            c_x[m]       = 0.0;
        }
    }
}


/* Branched top-hat vacuum one-sided flux. Per-basis flux on one side of a
 * face, in the SOURCE-side outward frame defined by Area_outward.
 *
 * CONTRACT (caller responsibility, NOT verified inside; device code):
 *   - moments[]      Upwind flux-frame Q (cell or face-reconstructed)
 *                    of the SOURCE basis. Layout: [m, p_x, p_y, p_z,
 *                    R_xx, R_yy, R_zz, R_xy, R_xz, R_yz] for NMOMENTS=10;
 *                    [m, p_x, p_y, p_z, R_xx, R_yy, R_zz] for NMOMENTS=7;
 *                    [m, p_x, p_y, p_z] for NMOMENTS=4. Slots [4..9] are
 *                    RAW second moments R_kl = <v_k v_l>, NOT central
 *                    stress S_kl. Central S = R - v⊗v is built locally.
 *   - vface[3]       Face velocity in code units (typically v_F_normal
 *                    times the canonical face unit normal; tangential
 *                    components are silently ignored by this function
 *                    because they cancel in the n_hat contraction).
 *   - Area_outward   SOURCE-side outward face area vector (NOT
 *                    normalized). For a pair sharing one face:
 *                       i-side call: Area_outward = +Face_Area_Vec
 *                       j-side call: Area_outward = -Face_Area_Vec
 *                    The function does NOT infer orientation from
 *                    sign(vsig) — misuse silently inverts the physics.
 *   - fluxes[]       Output, in Area_outward frame, area-integrated
 *                    (units: [Q] × velocity × area).
 *
 * RETURNS: signal speed (|u_out| + c_x) * |Area_outward|. The full
 * physical wave speed. Caller uses this for CFL and wakeup criteria.
 *
 * BRANCHING: with u_out = (<v> - vface) · n_out and
 * c_x = sqrt(gamma_e * n_out · S · n_out), gamma_e = 3. The per-basis flux
 * coefficients (F_m, pstress, gomega) come from cbe_face_flux_scalars, which
 * selects the scheme via CBE_INTEGRATOR_RP_GAUSSIAN (exact one-sided Gaussian)
 * vs the compact-support top-hat (#else). Top-hat fan:
 *     u_out >= c_x       -> F0 (cold supersonic outflow)
 *     -c_x < u_out < c_x -> warm fan (continuous into F0 at +c_x, into 0 at -c_x)
 *     u_out <= -c_x      -> 0  (no outgoing flux from this basis)
 * Cold limit S->0: c_x->0, warm fan collapses, F0 recovers exactly.
 *
 * DEFENSIVE GUARD: if rho is non-positive or non-finite or A_norm_sq <= 0
 * the function zeros the flux and returns 0. Device-safe (no endrun);
 * caller's clamp helpers should already filter these cases. */

/* Signal-speed prefactor on the normal velocity dispersion in the CBE timestep
 * estimate: vsig = |u_out| + CBE_SIGNAL_SIGMA_PREFAC * sigma_x, with
 * sigma_x = c_x/sqrt(3) = sqrt(n.S.n). Scheme-dependent: the top-hat has a
 * compact support with exact maximum signal speed |u_out| + c_x, i.e. prefactor
 * sqrt(3); the Gaussian tail is formally unbounded and uses a safety factor 3
 * (>= sqrt3). This affects ONLY the returned vsig estimate, not the flux physics. */
#if defined(CBE_INTEGRATOR_RP_GAUSSIAN)
static const double CBE_SIGNAL_SIGMA_PREFAC = 3.0;                  /* Gaussian: safety factor */
#else
static const double CBE_SIGNAL_SIGMA_PREFAC = 1.7320508075688772;  /* top-hat: exact |u|+c_x (sqrt3) */
#endif
KOKKOS_INLINE_FUNCTION
double cbe_flux_tophat_vacuum(const double moments[CBE_INTEGRATOR_NMOMENTS],
                            const double vface[3],
                            const double Area_outward[3],
                            double fluxes[CBE_INTEGRATOR_NMOMENTS])
{
    const double A_norm_sq = Area_outward[0]*Area_outward[0]
                           + Area_outward[1]*Area_outward[1]
                           + Area_outward[2]*Area_outward[2];
    if(!(moments[0] > MIN_REAL_NUMBER) || !isfinite(moments[0]) || !(A_norm_sq > 0)) {
        for(int k=0; k<CBE_INTEGRATOR_NMOMENTS; k++) fluxes[k] = 0;
        return 0;
    }
    const double A_norm  = sqrt(A_norm_sq);
    const double inv_A   = 1.0 / A_norm;
    const double n_hat[3] = { Area_outward[0]*inv_A, Area_outward[1]*inv_A, Area_outward[2]*inv_A };

    const double rho     = moments[0];
    const double inv_rho = 1.0 / rho;
    /* Always-3-vector velocity. Missing components zero-padded so the
     * v.n_hat dot product collapses cleanly to NUMDIMS terms when n_hat
     * has only active components (canonical 1D/2D face normals are
     * axis-aligned; inactive n_hat slots would contribute zero against
     * zero v). */
    double v[3]; cbe_basis_v_load_3(moments, v);
    const double v_n     = v[0]*n_hat[0]    + v[1]*n_hat[1]    + v[2]*n_hat[2];
    const double vF_n    = vface[0]*n_hat[0] + vface[1]*n_hat[1] + vface[2]*n_hat[2];
    const double u_out   = v_n - vF_n;

    /* Central stress contracted with n_hat: S_n_k = (R_kl - v_k v_l) n_hat_l.
     * Active-block only (k,l < NUMDIMS). With the dim-aware layout via
     * cbe_basis_T_load_active_scaled, R inactive entries are exactly zero
     * and v inactive components are zero, so the loop bound here is
     * cosmetic in 3D but load-bearing in 1D/2D where there is no
     * physical y/z stress to contract. c_x from the SSOT scalar helper. */
    double S_n[3] = {0};
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    {
        double R[3][3];
        cbe_basis_T_load_active_scaled(moments, inv_rho, R);
        for(int k=0; k<NUMDIMS; k++) {
            double s = 0;
            for(int l=0; l<NUMDIMS; l++) s += (R[k][l] - v[k]*v[l]) * n_hat[l];
            S_n[k] = s;
        }
    }
#endif
    const double c_x = cbe_face_normal_stress_speed_from_Qrow(moments, n_hat);

    /* Mass slot + scheme coefficients via the SSOT helper -- bit-identical to
     * what the v_F root-find residual sums over (it calls the same helper), so
     * basis-summed F_m at v_F == 0 exactly implies cell-summed mass
     * conservation, and the deposited flux can never disagree with the root
     * solve on the scheme (top-hat vs Gaussian, set by CBE_INTEGRATOR_RP_GAUSSIAN).
     * F = 0 vacuum branch returns short-circuit (no stress / momentum work). */
    const CbeFaceFluxScalars fs = cbe_face_flux_scalars(rho, u_out, c_x);
    const double F_m_per_area = fs.F_m;
    if(F_m_per_area == 0) {
        for(int k=0; k<CBE_INTEGRATOR_NMOMENTS; k++) fluxes[k] = 0;
        return (fabs(u_out) + CBE_SIGNAL_SIGMA_PREFAC * (c_x * (1.0/1.7320508075688772))) * A_norm;
    }
    fluxes[0] = F_m_per_area * A_norm;

    /* Pressure (pstress) and intrinsic-stress (gomega) flux coefficients from
     * the same helper. Top-hat: pstress = rho*zeta', gomega = rho*omega'/sigma_x.
     * Gaussian: pstress = rho*Phi, gomega = rho*phi/sigma_x (the S_n⊗S_n term). */
    const double pstress_coef = fs.pstress;
    const double gomega_coef  = fs.gomega;

    /* Momentum: F_p_k = v_k * F_m + pstress_coef * |A| * S_n_k. The flux row
     * has only NUMDIMS momentum slots in low-D builds; cbe_basis_p_w silently
     * drops writes to non-existent y/z slots. The S_n[k] math stays 3-vector
     * so the rotation-invariant cold limit (S=0 → fluxes[k+1] = v[k]*F_m)
     * collapses cleanly to the existing 1D code path with v_y=v_z=0. */
    for(int k=0; k<3; k++) {
        cbe_basis_p_w(fluxes, k, v[k] * fluxes[0] + pstress_coef * A_norm * S_n[k]);
    }

    /* Stress: F_T_ab = R_ab * F_m + pstress_coef * |A| * (v_a S_n_b + S_n_a v_b)
     * + gomega_coef * |A| * S_n_a S_n_b. The S_n⊗S_n term (gomega) is the exact
     * one-sided third-moment / heat-flux piece of the Gaussian flux (the top-hat
     * scheme uses its own compact-support gomega). Active block (a,b < NUMDIMS); the v_a S_n_b + S_n_a v_b
     * sum collapses to 2 v_a S_n_a when a==b. The R_ab * F_m piece carries RAW
     * R; the cross piece uses central S contracted with n_hat (S_n). Writes via
     * cbe_basis_T_w so inactive slots are silently untouched -- no OOB. */
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int a=0; a<NUMDIMS; a++) {
        const double R_aa = cbe_basis_T_r(moments, a, a) * inv_rho;
        cbe_basis_T_w(fluxes, a, a,
            R_aa * fluxes[0] + pstress_coef * A_norm * 2.0 * v[a] * S_n[a]
            + gomega_coef * A_norm * S_n[a] * S_n[a]);
        for(int b=a+1; b<NUMDIMS; b++) {
            const double R_ab = cbe_basis_T_r(moments, a, b) * inv_rho;
            cbe_basis_T_w(fluxes, a, b,
                R_ab * fluxes[0] + pstress_coef * A_norm * (v[a]*S_n[b] + S_n[a]*v[b])
                + gomega_coef * A_norm * S_n[a] * S_n[b]);
        }
    }
#endif

    return (fabs(u_out) + CBE_SIGNAL_SIGMA_PREFAC * (c_x * (1.0/1.7320508075688772))) * A_norm;
}


/* Conservative cell-state realizability repair.
 *
 * Single SSOT operator inside do_cbe_drift_kick_kernel. Operates on
 * RELATIVE-FRAME storage (the canonical pi.CBE_basis_moments layout after
 * the absolute round-trip).
 *
 * Conservation contract (preserved to roundoff except where explicitly
 * lossy-with-bounded-magnitude):
 *   Sum_b m_b               preserved via donation arithmetic
 *   Sum_b p_rel_b[a]        preserved (active a)
 *   Sum_b T_rel_b[a,b]      preserved (active a <= b, FULL active tensor)
 *   each row leaves realizable: m_b > 0 AND central S PSD in active dims
 *   no stochastic perturbations
 *   no global abort
 *
 * Residual sign convention (used throughout Step 5+):
 *   residual = target_cell_sum_BEFORE_repair - cell_sum_AFTER_local_repairs
 *   donor   += residual_share   =>   cell_sum_after == target_cell_sum
 *
 * Step outline:
 *   0  guard NaN/Inf, M0 <= 0 -> no-op return (upstream NaN bug class).
 *   1  snapshot target cell invariants (M0, P0, T0); V_bulk = P0/M0.
 *   2  identify empty bases (m_b < m_thresh).
 *   3  reset empties to comoving-inert in RELATIVE storage:
 *        m = m_thresh; p_rel = 0; T_rel = 0.
 *      (do NOT write m*V_bulk into stored slots; storage is
 *      relative-frame. An inert basis comoving with the cell bulk has
 *      p_rel = 0 and T_rel = 0 by definition.)
 *      Special case: if ALL bases are empty, skip reset entirely
 *      (would over-mass by NBASIS * m_thresh - M0).
 *   4  row-local PSD projection on non-empty bases. NOT cell-conservative;
 *      trace-deltas flow into Step 5 residual.
 *   5  compute residual vs Step 1 target; precompute cell scales.
 *   6  donate residual:
 *        (a) single-donor most-massive realizable; if donor stays
 *            realizable -> done.
 *        (b) multi-donor weighted by m_b at f=1; if all realizable -> done.
 *        (c) continuous bisection on f in [0,1] for the largest realizable
 *            multi-donor fraction; apply at f_lo; unabsorbed = (1-f_lo)*res.
 *   7  bounded acceptance + deterministic mass-positive proportional fallback:
 *        if |unabsorbed| <= EPS_CONS * cell_scale -> accept FP-scale loss
 *        else distribute the FULL unabsorbed residual proportionally over
 *          all bases by current mass fraction (w_b = m_b / M_cur, where
 *          M_cur = Sum m_b at Step 7 entry = M_new + f_eff*res_M). Adds
 *          exactly unab_M to the cell sum -> Sum m_b restored to M0
 *          (similarly for p_rel and T_rel). Mass positivity guaranteed
 *          since m_b_new = m_b * M0 / M_cur > 0 whenever M0 > 0 (Step 0
 *          guard) and M_cur > 0 (all m_b > 0 after Step 6).
 *   8  defensive final realizability sweep. Any non-realizable row after
 *      Step 7 gets row-local PSD projection (preserves m_b and p_b; touches
 *      stress only). Trace-delta bounded by NBASIS * CBE_SPD_RELATIVE_FLOOR *
 *      max_trace_central_S; FP-scale on smooth inputs (empirical).
 *
 * dT_out (nullable): accumulates Sum trace-delta from all PSD projections
 * (Step 4 + Step 8). Same contract as the prior legacy chain so downstream
 * cbe_diagnostics col-8 semantics carry through. */
KOKKOS_INLINE_FUNCTION
static void cbe_repair_basis_states_core(
    double (&moments)[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    double *dT_out)
{
    double dT_local = 0.0;

    /* ----- Step 0: input guard. ----- */
    {
        bool any_nonfinite = false;
        double M0_probe = 0.0;
        for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                if(!isfinite(moments[j][k])) { any_nonfinite = true; break; }
            }
            if(any_nonfinite) break;
            M0_probe += moments[j][0];
        }
        if(any_nonfinite || !(M0_probe > 0)) {
            if(dT_out) *dT_out = 0.0;
            return;
        }
    }

    /* ----- Step 1: snapshot target cell invariants. ----- */
    double M0 = 0.0;
    double P0[3] = {0.0, 0.0, 0.0};
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        M0 += moments[j][0];
        for(int a = 0; a < NUMDIMS; a++) P0[a] += cbe_basis_p_r(moments[j], a);
    }
    const double V_bulk[3] = { P0[0] / M0, P0[1] / M0, P0[2] / M0 };
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double T0[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        for(int a = 0; a < NUMDIMS; a++) {
            for(int b = a; b < NUMDIMS; b++) {
                T0[a][b] += cbe_basis_T_r(moments[j], a, b);
            }
        }
    }
#endif

    /* ----- Step 2: identify empties. ----- */
    const double m_thresh = DMAX(CBE_REPAIR_EPS_M_ABS,
                                  CBE_REPAIR_EPS_M_REL * M0 / (double)CBE_INTEGRATOR_NBASIS);
    bool is_empty[CBE_INTEGRATOR_NBASIS];
    int empty_count = 0;
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        is_empty[j] = (moments[j][0] < m_thresh);
        if(is_empty[j]) empty_count++;
    }

    /* ----- Step 3: reset empties to comoving-inert RELATIVE storage.
     * Inert basis comoving with V_bulk in absolute frame -> relative
     * storage: m = m_thresh, p_rel = 0, T_rel = 0. The cell-sum changes
     * caused by these assignments flow into Step 5 residual. Special-case:
     * if every basis is empty, skip reset (NBASIS*m_thresh > M0
     * over-masses the cell). ----- */
    if(empty_count > 0 && empty_count < CBE_INTEGRATOR_NBASIS) {
        for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
            if(!is_empty[j]) continue;
            moments[j][0] = m_thresh;
            for(int a = 0; a < NUMDIMS; a++) cbe_basis_p_w(moments[j], a, 0.0);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
            for(int a = 0; a < NUMDIMS; a++) {
                for(int b = a; b < NUMDIMS; b++) {
                    cbe_basis_T_w(moments[j], a, b, 0.0);
                }
            }
#endif
        }
    }

    /* ----- Step 4: row-local PSD projection on the (non-empty) goods.
     * Projection is NOT cell-conservative; the trace-deltas it introduces
     * flow into the Step 5 residual which the donor mechanism absorbs. */
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        if(is_empty[j]) continue;
        const double m_b = moments[j][0];
        if(!(m_b > 0)) continue;
        double v[3]; cbe_basis_v_load_3(moments[j], v);
        double v_dot_v = 0.0;
        for(int a = 0; a < NUMDIMS; a++) v_dot_v += v[a] * v[a];
        const double trace_R_active  = cbe_basis_T_trace_active(
            moments[j], 1.0 / m_b);
        const double trace_S_central = trace_R_active - v_dot_v;
        const double trace_S_pos     = DMAX(trace_S_central, MIN_REAL_NUMBER);
        const double eigenvalue_floor = CBE_SPD_RELATIVE_FLOOR * trace_S_pos;
        double dT_basis = 0.0;
        (void)cbe_basis_row_project_central_stress_to_PSD(
            moments[j], eigenvalue_floor, &dT_basis);
        dT_local += dT_basis;
    }
#endif

    /* ----- Step 5: residual = target - current_after_local_repairs.
     * Donor application later: donor += residual_share -> cell_sum
     * returns to target. ----- */
    double M_new = 0.0;
    double P_new[3] = {0.0, 0.0, 0.0};
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        M_new += moments[j][0];
        for(int a = 0; a < NUMDIMS; a++) P_new[a] += cbe_basis_p_r(moments[j], a);
    }
    double res_M    = M0 - M_new;
    double res_P[3] = { P0[0] - P_new[0], P0[1] - P_new[1], P0[2] - P_new[2] };
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double T_new[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        for(int a = 0; a < NUMDIMS; a++) {
            for(int b = a; b < NUMDIMS; b++) {
                T_new[a][b] += cbe_basis_T_r(moments[j], a, b);
            }
        }
    }
    double res_T[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int a = 0; a < NUMDIMS; a++) {
        for(int b = a; b < NUMDIMS; b++) {
            res_T[a][b] = T0[a][b] - T_new[a][b];
        }
    }
#endif

    /* Cell scales for FP-bounded acceptance test in Step 7. */
    double V_bulk_norm = 0.0;
    for(int a = 0; a < NUMDIMS; a++) V_bulk_norm += V_bulk[a] * V_bulk[a];
    V_bulk_norm = sqrt(V_bulk_norm);
    double max_abs_P0 = 0.0;
    for(int a = 0; a < NUMDIMS; a++) {
        const double ap = fabs(P0[a]); if(ap > max_abs_P0) max_abs_P0 = ap;
    }
    const double cell_scale_M = M0;
    const double cell_scale_P = DMAX(M0 * V_bulk_norm, max_abs_P0);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    double max_abs_T0 = 0.0;
    for(int a = 0; a < NUMDIMS; a++) {
        for(int b = a; b < NUMDIMS; b++) {
            const double at = fabs(T0[a][b]); if(at > max_abs_T0) max_abs_T0 = at;
        }
    }
    const double cell_scale_T = DMAX(M0 * V_bulk_norm * V_bulk_norm, max_abs_T0);
#endif

    /* ----- Step 6: donate residual.
     * (a) single-donor most-massive realizable
     * (b) multi-donor weighted by m_b at f=1
     * (c) continuous bisection on f in [0,1] ----- */
    bool   is_candidate[CBE_INTEGRATOR_NBASIS];
    double cand_m[CBE_INTEGRATOR_NBASIS];
    double cand_m_sum = 0.0;
    int    sd_idx = -1; double sd_m = -1.0;
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        is_candidate[j] = cbe_basis_row_is_realizable(
            moments[j], CBE_REPAIR_EPS_M_REL);
        cand_m[j]   = is_candidate[j] ? moments[j][0] : 0.0;
        cand_m_sum += cand_m[j];
        if(is_candidate[j] && moments[j][0] > sd_m) {
            sd_m = moments[j][0]; sd_idx = j;
        }
    }

    /* saved[]: per-row pre-donation state for revert / bisection restore. */
    double saved[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) saved[j][k] = moments[j][k];
    }

    bool single_ok = false;
    if(sd_idx >= 0 && cand_m_sum > 0) {
        moments[sd_idx][0] += res_M;
        for(int a = 0; a < NUMDIMS; a++)
            cbe_basis_p_a(moments[sd_idx], a, res_P[a]);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        for(int a = 0; a < NUMDIMS; a++) {
            for(int b = a; b < NUMDIMS; b++) {
                cbe_basis_T_w(moments[sd_idx], a, b,
                    cbe_basis_T_r(moments[sd_idx], a, b) + res_T[a][b]);
            }
        }
#endif
        single_ok = cbe_basis_row_is_realizable(
            moments[sd_idx], CBE_REPAIR_EPS_M_REL);
        if(!single_ok) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                moments[sd_idx][k] = saved[sd_idx][k];
        }
    }

    double f_eff = single_ok ? 1.0 : 0.0;

    if(!single_ok && cand_m_sum > 0) {
        /* Bisection bracket: f_lo always realizable (start at 0),
         * f_hi might not be. Initial probe at f=1. */
        double f_lo = 0.0, f_hi = 1.0;

        /* apply_at_f(f): restore all candidate rows from saved[], then add
         * (cand_m[j] / cand_m_sum) * f * residual to each candidate. */
        for(int probe = 0; probe < 2 + 40; probe++) {
            double f_try;
            if(probe == 0)      f_try = 1.0;
            else if(probe == 1) { /* if f=1 worked, stop */
                if(f_lo > 0.0) break;
                f_try = 0.5;
            } else {
                f_try = 0.5 * (f_lo + f_hi);
                if((f_hi - f_lo) < 1e-14) break;
            }
            /* Restore candidates. */
            for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
                if(!is_candidate[j]) continue;
                for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                    moments[j][k] = saved[j][k];
            }
            /* Apply weighted share at f_try. */
            for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
                if(!is_candidate[j]) continue;
                const double w = cand_m[j] / cand_m_sum;
                moments[j][0] += w * f_try * res_M;
                for(int a = 0; a < NUMDIMS; a++)
                    cbe_basis_p_a(moments[j], a, w * f_try * res_P[a]);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
                for(int a = 0; a < NUMDIMS; a++) {
                    for(int b = a; b < NUMDIMS; b++) {
                        cbe_basis_T_w(moments[j], a, b,
                            cbe_basis_T_r(moments[j], a, b) + w * f_try * res_T[a][b]);
                    }
                }
#endif
            }
            /* Check realizability of all candidates. */
            bool all_ok = true;
            for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
                if(!is_candidate[j]) continue;
                if(!cbe_basis_row_is_realizable(moments[j], CBE_REPAIR_EPS_M_REL)) {
                    all_ok = false; break;
                }
            }
            if(all_ok) {
                f_lo = f_try;
                if(f_try >= 1.0 - 1e-15) break;
            } else {
                f_hi = f_try;
            }
        }
        /* Final: ensure state reflects f_lo (the largest realizable f). */
        for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
            if(!is_candidate[j]) continue;
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                moments[j][k] = saved[j][k];
            const double w = cand_m[j] / cand_m_sum;
            moments[j][0] += w * f_lo * res_M;
            for(int a = 0; a < NUMDIMS; a++)
                cbe_basis_p_a(moments[j], a, w * f_lo * res_P[a]);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
            for(int a = 0; a < NUMDIMS; a++) {
                for(int b = a; b < NUMDIMS; b++) {
                    cbe_basis_T_w(moments[j], a, b,
                        cbe_basis_T_r(moments[j], a, b) + w * f_lo * res_T[a][b]);
                }
            }
#endif
        }
        f_eff = f_lo;
    }

    /* ----- Step 7: bounded acceptance + deterministic conservative fallback. */
    if(f_eff < 1.0 - 1e-15) {
        const double unabs_f = 1.0 - f_eff;
        const double unab_M    = unabs_f * res_M;
        const double unab_P[3] = { unabs_f*res_P[0], unabs_f*res_P[1], unabs_f*res_P[2] };
        double abs_unab_P = 0.0;
        for(int a = 0; a < NUMDIMS; a++) abs_unab_P += fabs(unab_P[a]);
        double abs_unab_T = 0.0;
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        double unab_T[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for(int a = 0; a < NUMDIMS; a++) {
            for(int b = a; b < NUMDIMS; b++) {
                unab_T[a][b] = unabs_f * res_T[a][b];
                abs_unab_T  += fabs(unab_T[a][b]);
            }
        }
#endif
        const bool M_fp = (fabs(unab_M) <= CBE_REPAIR_EPS_CONS * cell_scale_M);
        const bool P_fp = (abs_unab_P  <= CBE_REPAIR_EPS_CONS * cell_scale_P);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        const bool T_fp = (abs_unab_T  <= CBE_REPAIR_EPS_CONS * cell_scale_T);
#else
        const bool T_fp = true;
#endif
        if(!(M_fp && P_fp && T_fp)) {
            /* Deterministic mass-positive proportional fallback. Avoids
             * concentrating the FULL unabsorbed residual on a single donor,
             * which (a) could leave m_donor <= 0 for large negative unab_M,
             * and (b) would rely on a row-local PSD projection of unbounded
             * magnitude to clean up T.
             *
             * The denominator must be the CURRENT cell mass sum (post-Step-6
             * donation at f_eff), NOT
             * M_new (the pre-donation sum). At Step 7 entry the cell holds
             *   M_cur = M_new + f_eff * res_M
             * so applying (m_b / M_cur) * unab_M gives Sum w_j = 1 exactly
             * and the proportional distribution adds exactly unab_M to the
             * cell sum -> M_after = M_cur + unab_M = M_new + res_M = M0,
             * the target. Using M_new (the pre-donation sum) instead would
             * scale the donation by M_cur/M_new and break conservation when
             * f_eff > 0.
             *
             * Proportional distribution per-row:
             *   m_b   += (m_b / M_cur) * unab_M
             *   p_b   += (m_b / M_cur) * unab_P
             *   T_b   += (m_b / M_cur) * unab_T
             * preserves Sum m_b, Sum p_b, Sum T_b EXACTLY. Mass positivity:
             *   m_b_new = m_b * (M_cur + unab_M) / M_cur = m_b * M0 / M_cur
             * is positive whenever M0 > 0 (Step 0 guard) AND M_cur > 0
             * (Step 6 produces all candidate m_b > 0 by realizability, and
             * Step 3 reset gives empties m_b = m_thresh > 0). Per-row
             * realizability after proportional shift is NOT guaranteed -- the
             * defensive Step 8 sweep below row-locally projects to PSD if
             * needed. Step 8 projection touches only stress slots (preserves
             * m_b and p_b row-locally); trace-delta is bounded by
             *   NBASIS * CBE_SPD_RELATIVE_FLOOR * max_row_trace_central_S
             * and is FP-scale on smooth inputs (verified empirically by the
             * cold/warm Mac smoke; a non-FP magnitude signals upstream
             * pathological cell state surfaced via the dT_out path). */
            double M_cur = 0.0;
            for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) M_cur += moments[j][0];
            if(M_cur > 0) {
                for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
                    const double w_j = moments[j][0] / M_cur;
                    moments[j][0] += w_j * unab_M;
                    for(int a = 0; a < NUMDIMS; a++)
                        cbe_basis_p_a(moments[j], a, w_j * unab_P[a]);
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
                    for(int a = 0; a < NUMDIMS; a++) {
                        for(int b = a; b < NUMDIMS; b++) {
                            cbe_basis_T_w(moments[j], a, b,
                                cbe_basis_T_r(moments[j], a, b) + w_j * unab_T[a][b]);
                        }
                    }
#endif
                }
            }
            (void)abs_unab_P; (void)abs_unab_T;   /* used only in fp-test above */
        }
    }

    /* ----- Step 8: defensive final realizability sweep. Should not fire
     * post Step 7 (donor was already projected). Row-local PSD projection
     * preserves m_b individually; trace-delta accumulated into dT_local. */
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
    for(int j = 0; j < CBE_INTEGRATOR_NBASIS; j++) {
        if(cbe_basis_row_is_realizable(moments[j], CBE_REPAIR_EPS_M_REL))
            continue;
        const double m_b = moments[j][0];
        if(!(m_b > 0)) continue;
        double dT_last = 0.0;
        (void)cbe_basis_row_project_central_stress_to_PSD(
            moments[j], 0.0, &dT_last);
        dT_local += dT_last;
    }
#endif

    if(dT_out) *dT_out = dT_local;
}

/* Thin particle-object wrapper preserving the original conserved-path call
 * signature: forwards the particle's conserved basis moments to the
 * array-level core above. The predictor (later) calls the core directly on
 * the predicted (_pred) storage so conserved state is never touched. */
KOKKOS_INLINE_FUNCTION
static void cbe_repair_cell_conservative_basis_states(
    struct particle_data& pi,
    double *dT_out)
{
    cbe_repair_basis_states_core(pi.CBE_basis_moments, dT_out);
}


/* Donor-trial helper for the aggregate outflow-budget limiter.
 * Tests whether ALL candidate donors stay row-local realizable when fraction
 * f of the excess is redistributed from basis a. A file-scope
 * KOKKOS_INLINE_FUNCTION (rather than a captured lambda inside
 * cbe_apply_basis_outflow_budget) is device-portable under CUDA/HIP, where a
 * captured C++ lambda inside another KOKKOS_INLINE_FUNCTION is a portability
 * hazard.
 *
 * Inputs:
 *   pi          particle (read-only access to moments + moments_dt + Vel)
 *   dt          drift-kick step length
 *   V_old[3]    canonical frame for the absolute-round-trip (== pi.Vel; EXACT
 *               for row-local realizability by frame-invariance)
 *   is_cand[]   per-basis candidate flag set by the limiter loop
 *   cap_b[]     per-basis capacity weights (mass-floor-deducted)
 *   C_sum       Sum of capacities (> 0 on entry, checked by caller)
 *   excess[]    saved-outflow rates per moment slot (positive for slot 0)
 *   f           trial fraction in [0, 1]
 *
 * Returns true iff every candidate donor's projected post-update state
 * (existing rate + limiter correction * dt, frame-converted via V_old) is
 * mass-positive AND row-local realizable. */
KOKKOS_INLINE_FUNCTION
static bool cbe_outbudget_trial_donors_realizable(
    const struct particle_data& pi,
    double dt,
    const double V_old[3],
    const bool   is_cand[CBE_INTEGRATOR_NBASIS],
    const double cap_b[CBE_INTEGRATOR_NBASIS],
    double C_sum,
    const double excess[CBE_INTEGRATOR_NMOMENTS],
    double f)
{
    for(int b = 0; b < CBE_INTEGRATOR_NBASIS; b++) {
        if(!is_cand[b]) continue;
        const double w_b = cap_b[b] / C_sum;
        double m_b_abs;
        double p_abs[3];
        double T_abs[3][3];
        cbe_relative_row_to_absolute(pi.CBE_basis_moments[b], V_old,
                                      &m_b_abs, p_abs, T_abs);
        double rate_corr[CBE_INTEGRATOR_NMOMENTS];
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
            rate_corr[k] = pi.CBE_basis_moments_dt[b][k] - w_b * f * excess[k];
        }
        m_b_abs += dt * rate_corr[0];
        if(!(m_b_abs > 0)) return false;
        for(int aa = 0; aa < NUMDIMS; aa++) {
            p_abs[aa] += dt * cbe_basis_p_r(rate_corr, aa);
        }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        for(int aa = 0; aa < NUMDIMS; aa++) {
            T_abs[aa][aa] += dt * cbe_basis_T_r(rate_corr, aa, aa);
            for(int bb = aa+1; bb < NUMDIMS; bb++) {
                const double dT_ab = dt * cbe_basis_T_r(rate_corr, aa, bb);
                T_abs[aa][bb] += dT_ab; T_abs[bb][aa] += dT_ab;
            }
        }
#endif
        double row_proj[CBE_INTEGRATOR_NMOMENTS];
        cbe_absolute_to_relative_row(m_b_abs, p_abs, T_abs, V_old, row_proj);
        if(!cbe_basis_row_is_realizable(row_proj, CBE_REPAIR_EPS_M_REL)) return false;
    }
    return true;
}


/* Aggregate per-basis outflow-budget limiter. For each basis a whose dt*outflow_rate
 * exceeds f_out * m_a, scale the outflow rates by c = f_out*m_a/dt/Udot_out0
 * and redistribute the saved (1-c)*Udot_out among capacity-weighted donor
 * bases (b != a). Cell-summed dt-accumulator is preserved exactly by
 * construction (every excess[k] added to basis a is subtracted from donors
 * weighted to sum to excess[k]).
 *
 * Inputs read:
 *   pi.CBE_basis_moments[b][0]     current mass per basis
 *   pi.CBE_basis_out_rate_dt[a][k] cumulative outgoing-only flux (stored
 *                                  NEGATIVE; -ledger -> positive out rate)
 *   pi.CBE_basis_moments_dt[b][k]  full net flux rate (existing) -- used by
 *                                  the projected-donor realizability check
 *                                  (the donor check MUST include donor's
 *                                  existing rate, not just the limiter
 *                                  correction)
 *   pi.Vel                         V_old for the frame round-trip
 *
 * Inputs modified:
 *   pi.CBE_basis_moments_dt[][]    a += excess[k]; donors -= w_b * f * excess[k]
 *
 * Frame-consistency for the donor PSD check (EXACT, not
 * approximate): row-local realizability `M = m*T - p p^T PSD` is
 * frame-invariant under affine velocity transforms when the absolute<->relative
 * conversion uses the SAME V on both sides. Using V_old here is therefore
 * the right choice -- not an approximation. The drift-kick will later use
 * V_new = MMV from the post-update MMV; the V_new drift does not affect
 * row-local realizability per the same invariance.
 *
 * Projected donor state INCLUDES donor's existing rate:
 *   U_abs_proj[b] = U_abs_old[b] + dt * (CBE_basis_moments_dt[b] + apply_to_b)
 * NOT just dt * apply_to_b. This catches the failure where a donor is
 * current-state realizable but already near overdrain from its own flux
 * and the limiter's added load tips it past the realizability edge.
 *
 * Bisection on f: if full donation (f=1) breaks any donor's projected
 * realizability, bisect f down to the largest fraction that keeps all
 * donors realizable. The (1-f)*excess REMAINDER is left uncapped on basis
 * a's net rate (a loses more mass than dm_allowed would ideally permit but
 * still less than uncapped); reported via diag->n_donor_limited. NO silent
 * pretense that the cap fully protected basis a.
 *
 * Failure modes (printf/counter only -- NO endrun, NO new diagnostic
 * columns):
 *   NO_DONOR        : no eligible candidate basis (NBASIS=1 effective or
 *                     all bases at floor). Basis a stays UNCAPPED.
 *   DONOR_LIMITED   : bisection found f<1; partial cap applied + remainder
 *                     uncapped. Counted.
 *
 * Diag is nullable. Returns max-severity status code observed.
 *
 * This helper deliberately avoids any new endrun codes. */
KOKKOS_INLINE_FUNCTION
static int cbe_apply_basis_outflow_budget(
    struct particle_data& pi,
    double dt,
    struct CbeOutflowBudgetDiag *diag)
{
    int status = CBE_OUTBUDGET_OK;
    if(!(dt > 0)) return status;

    const double V_old[3] = { pi.Vel[0], pi.Vel[1], pi.Vel[2] };
    const double f_out    = CBE_OUTFLOW_BUDGET_FRAC;
    const double m_floor  = CBE_OUTFLOW_BUDGET_MIN_DONOR_MASS;

    /* Donor projected-realizability test is in the standalone
     * cbe_outbudget_trial_donors_realizable helper above (file-scope
     * KOKKOS_INLINE_FUNCTION; device-portable). */

    for(int a = 0; a < CBE_INTEGRATOR_NBASIS; a++) {
        const double m_a       = pi.CBE_basis_moments[a][0];
        const double Udot_out0 = -pi.CBE_basis_out_rate_dt[a][0]; /* positive outflow rate */
        if(!(m_a > 0) || !(Udot_out0 > 0) || !isfinite(m_a) || !isfinite(Udot_out0)) continue;

        const double dm_out     = dt * Udot_out0;
        const double dm_allowed = f_out * m_a;
        const double pre_of     = dm_out / m_a;
        if(diag && pre_of > diag->max_pre_out_frac_seen) diag->max_pre_out_frac_seen = pre_of;
        if(dm_out <= dm_allowed) continue;   /* basis a within budget */

        const double c = dm_allowed / dm_out;            /* c in (0, 1) */
        const double one_minus_c = 1.0 - c;
        double excess[CBE_INTEGRATOR_NMOMENTS];
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
            excess[k] = one_minus_c * (-pi.CBE_basis_out_rate_dt[a][k]);
        }

        /* Build candidate-donor set with capacity weights. Capacity is the
         * conservative v1 form `max(0, m_b - m_floor)` -- documented as not
         * mathematically optimal but the projected-state bisection below is
         * the safety rail. */
        bool   is_cand[CBE_INTEGRATOR_NBASIS];
        double cap_b[CBE_INTEGRATOR_NBASIS];
        double C_sum = 0.0;
        for(int b = 0; b < CBE_INTEGRATOR_NBASIS; b++) {
            is_cand[b] = false; cap_b[b] = 0.0;
            if(b == a) continue;
            const double m_b = pi.CBE_basis_moments[b][0];
            if(!(m_b > m_floor) || !isfinite(m_b)) continue;
            /* current-state quick filter; the real check is the projected one below */
            if(!cbe_basis_row_is_realizable(pi.CBE_basis_moments[b], CBE_REPAIR_EPS_M_REL)) continue;
            const double cap = m_b - m_floor;
            if(!(cap > 0)) continue;
            is_cand[b] = true; cap_b[b] = cap; C_sum += cap;
        }
        if(!(C_sum > 0)) {
            if(diag) diag->n_donor_failures++;
            status = (status > CBE_OUTBUDGET_NO_DONOR) ? status : CBE_OUTBUDGET_NO_DONOR;
            continue;   /* basis a uncapped */
        }

        /* Find largest f in [0,1] with all donors projected-realizable. The
         * trial helper is file-scope KOKKOS_INLINE_FUNCTION (above) for
         * device portability. */
        double f_lo;
        if(cbe_outbudget_trial_donors_realizable(pi, dt, V_old, is_cand, cap_b, C_sum, excess, 1.0)) {
            f_lo = 1.0;
        } else {
            f_lo = 0.0;
            double f_hi = 1.0;
            for(int iter = 0; iter < 40; iter++) {
                const double f_mid = 0.5 * (f_lo + f_hi);
                if(cbe_outbudget_trial_donors_realizable(pi, dt, V_old, is_cand, cap_b, C_sum, excess, f_mid))
                    f_lo = f_mid;
                else
                    f_hi = f_mid;
                if((f_hi - f_lo) < 1e-14) break;
            }
        }

        /* set at f=f_lo. The (1-f_lo)*excess REMAINDER stays on basis a's
         * net rate uncapped -- explicitly NOT a hidden full cap. */
        for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
            pi.CBE_basis_moments_dt[a][k] += f_lo * excess[k];
            for(int b = 0; b < CBE_INTEGRATOR_NBASIS; b++) {
                if(!is_cand[b]) continue;
                const double w_b = cap_b[b] / C_sum;
                pi.CBE_basis_moments_dt[b][k] -= w_b * f_lo * excess[k];
            }
        }

        if(f_lo >= 1.0 - 1e-15) {
            if(diag) diag->n_capped_bases++;
            status = (status >= CBE_OUTBUDGET_APPLIED) ? status : CBE_OUTBUDGET_APPLIED;
        } else {
            if(diag) { diag->n_donor_limited++; diag->n_donor_failures++; }
            status = (status > CBE_OUTBUDGET_DONOR_LIMITED) ? status : CBE_OUTBUDGET_DONOR_LIMITED;
        }
    }

    return status;
}


/* SSOT frame round-trip used by the conserved CBE drift-kick (and, later,
 * the predictor). Encapsulates kick Steps 1-5: read RELATIVE basis rows in
 * the V_old frame, convert to ABSOLUTE, apply the absolute-frame moment rate
 * scaled by `coeff` (= nfac*dt for the conserved kick), derive V_new from the
 * updated mass-weighted-mean basis velocity, and convert back to RELATIVE
 * storage in the V_new frame. `rate` is an ABSOLUTE-frame rate array (e.g.
 * CBE_basis_moments_dt). In-place use (src_rows == dst_rows) is safe: every
 * source row is read into locals before any destination row is written.
 * V_new is computed 3-wide (inactive dims default to 0, or V_old in the
 * m<=0 fallback); the caller decides which components to write. */
KOKKOS_INLINE_FUNCTION
static void cbe_apply_moment_rate_framecorrect(
    const double src_rows[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    const double V_old[3],
    const double rate[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    double coeff,
    double dst_rows[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS],
    double V_new[3])
{
    /* Step 1+2: relative storage (V_old frame) -> absolute per basis. */
    double m_b_arr[CBE_INTEGRATOR_NBASIS];
    double p_abs_arr[CBE_INTEGRATOR_NBASIS][3];
    double T_abs_arr[CBE_INTEGRATOR_NBASIS][3][3];
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
        cbe_relative_row_to_absolute(src_rows[j], V_old,
                                     &m_b_arr[j], p_abs_arr[j], T_abs_arr[j]);
    }

    /* Step 3: apply coeff * absolute-frame rate. */
    double m_new_total = 0;
    double p_abs_new_total[3] = {0,0,0};
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
        m_b_arr[j] += coeff * rate[j][0];
        m_new_total += m_b_arr[j];
        for(int a=0;a<NUMDIMS;a++) {
            p_abs_arr[j][a] += coeff * cbe_basis_p_r(rate[j], a);
            p_abs_new_total[a] += p_abs_arr[j][a];
        }
#if defined(CBE_INTEGRATOR_SECONDMOMENT)
        for(int a=0;a<NUMDIMS;a++) {
            T_abs_arr[j][a][a] += coeff * cbe_basis_T_r(rate[j], a, a);
            for(int b=a+1;b<NUMDIMS;b++) {
                const double dT_ab = coeff * cbe_basis_T_r(rate[j], a, b);
                T_abs_arr[j][a][b] += dT_ab; T_abs_arr[j][b][a] += dT_ab;
            }
        }
#endif
    }

    /* Step 4: V_new = mass-weighted-mean basis velocity (defensive copy of
     * V_old when total mass is non-positive). */
    V_new[0] = 0; V_new[1] = 0; V_new[2] = 0;
    if(m_new_total > 0) {
        for(int a=0;a<NUMDIMS;a++) V_new[a] = p_abs_new_total[a] / m_new_total;
    } else {
        for(int a=0;a<3;a++) V_new[a] = V_old[a];
    }

    /* Step 5: absolute -> relative storage in the V_new frame. */
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
        cbe_absolute_to_relative_row(m_b_arr[j], p_abs_arr[j], T_abs_arr[j],
                                     V_new, dst_rows[j]);
    }
}

/* GPU-callable per-particle drift-kick update for the CBE integrator.
 * Mirrors do_cbe_drift_kick() in cbe_integrator.cc but takes an explicit
 * particle ref so it runs in both CPU and GPU (Kokkos) contexts.
 *
 * Cell-state realizability is handled by a single call to
 * cbe_repair_cell_conservative_basis_states, following the absolute
 * round-trip (nfac + frame conversion + V_new derivation).
 *
 * *dT_out (nullable) receives the sum-over-bases trace-delta from PSD
 * projections inside the repair helper; caller passes nullptr to skip
 * diagnostic bookkeeping. The relative basis momenta sum to zero after the
 * kick; the particle's own momentum changes by Mass*(V_new-V_old), which is
 * added to pi.dp. */
KOKKOS_INLINE_FUNCTION
static void do_cbe_drift_kick_kernel(struct particle_data& pi, double dt,
                                     double *dT_out)
{
    /* Finite absolute-update round-trip.
     *
     * Per particle, per step:
     *   1. Read RELATIVE storage in current pi.Vel frame (V_old).
     *   2. Convert each basis row to ABSOLUTE (m, p_abs, T_abs) via the
     *      SSOT cbe_relative_row_to_absolute helper.
     *   3. Apply dt * dt_accumulator IN ABSOLUTE FRAME (the flux deposited
     *      the rates there; postgrav now leaves them in absolute except for
     *      a mass-closure safety net).
     *   4. Derive V_new from updated mass-weighted-mean basis velocity:
     *         V_new[k] = Σ_b p_abs_new_b[k] / Σ_b m_new_b
     *      This makes Σ p_rel = 0 exact by construction after step 5 and
     *      eliminates the need to inject a_cbe into pi.GravAccel.
     *   5. Convert ABSOLUTE → RELATIVE storage in V_new frame via the SSOT
     *      cbe_absolute_to_relative_row helper.
     *   6. Set pi.Vel = V_new (CBE-induced bulk velocity is now derived,
     *      NOT routed through GravAccel; gravity kicks are independent).
     *   7. Conservative cell-state realizability repair via SSOT helper
     *      cbe_repair_cell_conservative_basis_states.
     *
     * The nfac rate-limit (capping per-basis |dm|/m at 0.75) is preserved
     * from the legacy operator as a robustness safeguard. The legacy
     * mass + momentum closure corrections (per-basis subtraction of
     * m_b/M_total · dmom_tot[k]) are subsumed by step 4: deriving V_new
     * from the updated MMV gives Σ p_rel = 0 by construction; mass closure
     * is already enforced upstream by postgrav.
     *
     * dT_out (nullable) accumulates the sum-over-bases trace-delta from
     * PSD projections inside the Step 7 repair helper (for
     * cbe_diagnostics.txt col-8). The repair preserves Σ p_rel (projection
     * only touches stress slots), so the relative basis momenta still sum to
     * zero; the bulk momentum change is recorded in step 6. */

    /* Aggregate outflow-budget limiter.
     * Runs BEFORE nfac so subsequent steps operate on the capped rates.
     * No-op on the cold/warm 1D smoke tests we validate against (basis
     * outflows are well within f_out*m_a); fires only in the scratch-basis
     * runaway regime it was designed for. */
    {
        struct CbeOutflowBudgetDiag obdiag_local{};
        int obstatus = cbe_apply_basis_outflow_budget(pi, dt, &obdiag_local);
#ifdef CBE_INTEGRATOR_OUTBUDGET_VERBOSE
        /* Opt-in verbose smoke verification: default OFF.
         * Enable by adding CBE_INTEGRATOR_OUTBUDGET_VERBOSE to Config.sh and
         * grep stdout for "[CBE outbudget]" -- zero hits on clean smoke =
         * limiter did not fire. Device-safe printf (no fflush/endrun). */
        if(obdiag_local.n_capped_bases > 0 || obdiag_local.n_donor_limited > 0
                                            || obdiag_local.n_donor_failures > 0) {
            printf("[CBE outbudget] ID=%lld caps=%lld donor_lim=%lld donor_fail=%lld "
                   "max_pre_out_frac=%.3e status=%d\n",
                   (long long)pi.ID,
                   (long long)obdiag_local.n_capped_bases,
                   (long long)obdiag_local.n_donor_limited,
                   (long long)obdiag_local.n_donor_failures,
                   obdiag_local.max_pre_out_frac_seen, obstatus);
        }
#else
        (void)obstatus;
#endif
    }

    /* Step 0 (preserved): nfac rate-limit on per-basis mass change. */
    double dmass_tot = 0;
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) dmass_tot += dt * pi.CBE_basis_moments_dt[j][0];
    const double minv = (pi.Mass > 0) ? 1.0 / pi.Mass : 0.0;
    double biggest_dm = 1.e10;
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
        const double m_b = pi.CBE_basis_moments[j][0];
        if(!(m_b > 0)) continue;
        const double q = (dt*pi.CBE_basis_moments_dt[j][0] - m_b*minv*dmass_tot)
                       / (m_b * (1.0 + minv*dmass_tot));
        if(!isnan(q) && q < biggest_dm) biggest_dm = q;
    }
    double nfac = 1.0;
    if(biggest_dm < -0.75) nfac = -0.75 / biggest_dm;

    /* Steps 1-5 via the SSOT frame round-trip helper: read RELATIVE storage
     * in the V_old=pi.Vel frame, convert to ABSOLUTE, apply the nfac-limited
     * absolute-frame rate (coeff = nfac*dt), derive V_new from the updated
     * mass-weighted-mean basis velocity, and convert back to RELATIVE storage
     * in the V_new frame. In-place (src==dst==pi.CBE_basis_moments) is safe. */
    const double V_old[3] = { pi.Vel[0], pi.Vel[1], pi.Vel[2] };
    double V_new[3];
    cbe_apply_moment_rate_framecorrect(pi.CBE_basis_moments, V_old,
                                       pi.CBE_basis_moments_dt, nfac * dt,
                                       pi.CBE_basis_moments, V_new);

    /* Step 6: pi.Vel <- V_new (CBE-induced bulk velocity, derived not
     * GravAccel-injected). Gravity kicks against pi.GravAccel remain
     * independent and apply only EXTERNAL gravitational acceleration.
     * The momentum change goes into pi.dp, as every kick records it, so the
     * tree nodes' mean velocities follow this particle. */
    for(int a=0;a<NUMDIMS;a++) {pi.dp[a] += pi.Mass * (V_new[a] - V_old[a]); pi.Vel[a] = V_new[a];}

    /* Step 7: conservative cell-state realizability repair via the single
     * SSOT helper defined above. Operates on the now-frame-consistent
     * relative storage. */
    cbe_repair_cell_conservative_basis_states(pi, dT_out);

    /* Predictor reset: after the conserved kick, snap the predicted state to
     * the new conserved state; subsequent drifts advance pred from here and
     * the flux reads pred. Writes ONLY the derived _pred fields (conserved
     * state untouched). NOTE: pi.Vel here is the current post-kick CBE frame —
     * do_the_kick has already applied gravity before this CBE kick, and V_new
     * was derived from that gravity-inclusive frame. So CBE_VelPred starts from
     * the current gravity-inclusive frame; the predictor adds only the gravity
     * increment over later drift intervals (no double counting). */
    for(int jp=0; jp<CBE_INTEGRATOR_NBASIS; jp++)
        for(int kp=0; kp<CBE_INTEGRATOR_NMOMENTS; kp++)
            pi.CBE_basis_moments_pred[jp][kp] = pi.CBE_basis_moments[jp][kp];
    for(int ap=0; ap<3; ap++) pi.CBE_VelPred[ap] = pi.Vel[ap];
}

/* Predicted-state drift for the adaptive-timestep CBE predictor. Advances the
 * derived (CBE_basis_moments_pred, CBE_VelPred) over a drift interval so an
 * active short-dt neighbor reading these sees an estimate of where this cell
 * WILL be, instead of its stale begin-of-step conserved reservoir (the cause
 * of the adaptive over-siphon). Mirrors do_cbe_drift_kick_kernel on the _pred
 * storage: gravity enters via the FRAME velocity V_old_pred (matching
 * do_the_kick's GravAccel*dt_gravkick update of pi.Vel, which the conserved
 * CBE kick then reads as V_old), then the SSOT frame round-trip applies the
 * last-active CBE moment rate over dt_drift, then a realizability repair on the
 * predicted rows. Writes ONLY the _pred fields; conserved state untouched.
 * dt_gravkick / dt_gravkick_pm come from the same get_gravkick_factor
 * machinery as the gas VelPred prediction (predict.cc). */
KOKKOS_INLINE_FUNCTION
static void do_cbe_predict_drift_kernel(struct particle_data& pi, double dt_drift,
                                        double dt_gravkick, double dt_gravkick_pm)
{
    if(!(pi.Mass > 0)) return;

    /* Fold the gravity drift into the frame velocity (NOT the relative
     * moments): a uniform velocity shift leaves each basis's central stress
     * invariant and only advects the mean. GravAccel*dt_gravkick (no cf_a2inv,
     * matching do_the_kick); PM long-range piece by analogy to gas VelPred. */
    double V_old_pred[3];
    for(int a=0;a<3;a++) V_old_pred[a] = pi.CBE_VelPred[a] + pi.GravAccel[a]*dt_gravkick;
#ifdef PMGRID
    for(int a=0;a<3;a++) V_old_pred[a] += pi.GravPM[a]*dt_gravkick_pm;
#else
    (void)dt_gravkick_pm;
#endif

    double V_new_pred[3];
    double pred_coeff = dt_drift;
    cbe_apply_moment_rate_framecorrect(pi.CBE_basis_moments_pred, V_old_pred,
                                       pi.CBE_basis_moments_dt, pred_coeff,
                                       pi.CBE_basis_moments_pred, V_new_pred);
    for(int a=0;a<3;a++) pi.CBE_VelPred[a] = V_new_pred[a];

    /* Realizability repair on the PREDICTED rows only (array-level core, never
     * the conserved-state wrapper). dT is local/discarded — predicted-repair
     * diagnostics stay out of cbe_diagnostics until the flux consumes pred. */
    double dT_pred = 0.0;
    cbe_repair_basis_states_core(pi.CBE_basis_moments_pred, &dT_pred);
    (void)dT_pred;
}

/* GPU-callable per-particle post-gravity finalization for the CBE integrator.
 * This routine is minimal — it ONLY enforces basis-
 * mass closure (Σ_α dm_α/dt = 0 to FP) as a safety net. CBE bulk velocity
 * is NOT injected into pi.GravAccel anymore; the absolute round-trip in
 * do_cbe_drift_kick_kernel handles frame conversion and derives V_new from
 * the post-update mass-weighted-mean basis velocity.
 *
 * Moment ordering (unchanged): 0, x, y, z, xx, yy, zz, xy, xz, yz.
 *
 * Operations (pure i-side):
 *   - Sum dmom_tot[0] = Σ_α dm_α/dt across bases.
 *   - Subtract that residual from each basis's dm_α/dt weighted by m_α/M
 *     so the cell-summed mass rate is exactly zero (upstream MFM should
 *     deliver this to FP already; this is the safety net).
 *   - Leave all other dt-accumulator slots (momentum, stress) untouched —
 *     they stay in absolute frame for the drift-kick round-trip to consume.
 */
KOKKOS_INLINE_FUNCTION
static void do_cbe_postgravity_kernel(struct particle_data& pi)
{
    /* The OLD postgrav (a) injected a_cbe into
     * pi.GravAccel — routing CBE bulk velocity through the gravity kick —
     * and (b) converted dt-accumulator's per-basis dp_abs/dT_abs to relative
     * frame using a continuous-time differential formula. Both were
     * structural bugs: (a) double-counted CBE bulk velocity, and (b) lacked
     * the O(dt²) finite-step terms needed for V-changing-during-step accuracy.
     *
     * The new design moves the frame conversion into a finite absolute
     * round-trip inside do_cbe_drift_kick_kernel (see SSOT helpers
     * cbe_relative_row_to_absolute / cbe_absolute_to_relative_row above).
     * The flux already deposits dt-accumulator in absolute frame; postgrav
     * leaves those rates alone. CBE-induced bulk velocity is now derived
     * from the post-update mass-weighted-mean basis velocity inside the
     * drift-kick operator and written directly to pi.Vel. pi.GravAccel
     * carries ONLY external (real) gravitational acceleration.
     *
     * The only postgrav job left is mass-closure: enforce Σ_α dm_α/dt = 0
     * exactly as a safety net (upstream MFM should already deliver this
     * to machine eps). This is frame-invariant and is identical in absolute
     * and relative formulations. */
    double dmom_tot_mass = 0;
    for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
        dmom_tot_mass += pi.CBE_basis_moments_dt[j][0];
    }
    if(pi.Mass > 0) {
        const double m_inv = 1.0 / pi.Mass;
        for(int j=0;j<CBE_INTEGRATOR_NBASIS;j++) {
            pi.CBE_basis_moments_dt[j][0] -= pi.CBE_basis_moments[j][0] * (m_inv * dmom_tot_mass);
        }
    }
}

#endif /* CBE_INTEGRATOR */

#endif /* CBE_INTEGRATOR_FUNCTIONS_H */
