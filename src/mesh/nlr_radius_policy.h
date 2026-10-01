/* nlr_radius_policy.h — SSOT for the j-side symmetric-pair-search radius
 * under a Spec's `radius_policy`.
 *
 * Background: Mode A and Mode B both need a uniform answer to "what radius
 * sources should symmetric pair search use for particle j, given THIS Spec's
 * policy?" — and that answer depends on physics (a sink-loop Spec wants
 * P[j].KernelRadius for non-gas; an AGS-force Spec wants P[j].AGS_KernelRadius;
 * sink_feed wants max(KernelRadius, ForceSoftening); etc.).  This file is the
 * single canonical implementation of that decision.
 *
 * The corresponding NODE-level Mode B tree-prune (mode_b_node_symmetric_radius
 * in mesh/mode_b_local_walker.cc) does NOT use this helper; it reads the
 * conservative-upper-bound Extnodes[no].hmax_per_type[t] bands directly, which
 * are seeded via force_hmax_per_type_particle_radius() in gravity/forcetree.cc
 * to dominate every leaf-policy-selectable source for that type.  Node-prune
 * over-opens; this leaf-level helper is the exact policy filter.  Together
 * the two ensure tree walk reaches every leaf the leaf predicate admits.
 *
 * Callable form:
 *
 *   nlr_symmetric_radius_from_fields(type, kr, ags_kr, fs, policy)
 *     KOKKOS_INLINE_FUNCTION (device-callable via the macro fallback pattern
 *     below).  Takes the relevant per-particle fields by value, so a device
 *     kernel can call it inside a Kokkos lambda without dereferencing host
 *     structs from device.
 *     Callers MUST pass 0 for ags_kernel_radius when
 *     AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE is undefined (the field doesn't
 *     exist on P[]); pass (double)P[j].AGS_KernelRadius otherwise.
 *     Host call sites (e.g. mode_b_neighbor_symmetric_radius) read the
 *     fields with the right #ifdef gating and forward to this single form.
 *
 * KOKKOS_INLINE_FUNCTION fallback: per feedback_gpu §B.4b, this header must
 * never #include <Kokkos_Core.hpp> unconditionally.  GPU TUs that include
 * this file have already pulled Kokkos upstream and get the real
 * `inline __device__ __host__ __forceinline__` expansion; host-only TUs get
 * plain `inline` without any CUDA header coupling.
 */

#ifndef NLR_RADIUS_POLICY_H
#define NLR_RADIUS_POLICY_H

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

/* The per-particle wrappers below read KernelRadius / AGS_KernelRadius /
 * ForceSoftening / Type off `const struct particle_data&`.  particle_data.h
 * lacks an include guard, so we forward-declare here and require each TU
 * including this header to have `allvars.h` (which pulls particle_data.h
 * exactly once) visible BEFORE the include of nlr_radius_policy.h.  Every
 * existing caller already follows this order — sfc_tiles.cc / gpu_neighbor_list.cc
 * / mode_b_local_walker.cc / forcetree.cc / ghost_exchange.cc all include
 * allvars.h first.  Spec _loop.cc files inherit the order via mode_b_local_walker.h.
 *
 * The compile-flag gate on AGS_KernelRadius lives ONLY in the wrappers below;
 * callers MUST go through them and MUST NOT add their own
 * `#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE` for radius-policy use. */
struct particle_data;

/* Type-aware "symmetric radius" policy bitmask.
 *
 * For SYMMETRIC pair search, the test is r < max(h_i, h_j).  h_j depends on
 * j's type AND on which physics the Spec implements:
 *   - Gas-only pair physics (hydro density/gradient/hydro_force, mechfb, etc.):
 *     h_j = P[j].KernelRadius for Type==0, undefined for non-gas (mask excludes).
 *     Use MODE_B_RADIUS_DEFAULT = MODE_B_RADIUS_GAS_KERNEL.
 *   - AGS-style pair physics (ags_force, cbe_grad, dm_fuzzy_dmgrad-style):
 *     h_j = P[j].AGS_KernelRadius for the type-mask, possibly including gas.
 *     Use MODE_B_RADIUS_GAS_AGS | MODE_B_RADIUS_NONGAS_AGS as appropriate.
 *   - Sink/grain pair physics (sink_env1/env2/swk, grain_rt_gas):
 *     h_j = P[j].KernelRadius for non-gas (the sink's neighbor-search radius is
 *     the relevant pair-physics extent — sink reading gas density, sink-sink
 *     merge candidates, etc.).  Use MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL.
 *   - sink_feed binary-merge:
 *     h_eff_j = DMAX(P[j].KernelRadius, P[j].ForceSoftening).
 *     Use MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL | MODE_B_RADIUS_FORCE_SOFTENING.
 *
 * Bit assignments are stable; do NOT renumber.  Adding new bits goes at
 * (1u << 5) onward.
 */
typedef unsigned int mode_b_radius_policy_t;
#define MODE_B_RADIUS_GAS_KERNEL        (1u << 0)  /* Type==0:  P[j].KernelRadius */
#define MODE_B_RADIUS_GAS_AGS           (1u << 1)  /* Type==0:  P[j].AGS_KernelRadius */
#define MODE_B_RADIUS_NONGAS_KERNEL     (1u << 2)  /* Type!=0:  P[j].KernelRadius */
#define MODE_B_RADIUS_NONGAS_AGS        (1u << 3)  /* Type!=0:  P[j].AGS_KernelRadius */
#define MODE_B_RADIUS_FORCE_SOFTENING   (1u << 4)  /* any type: P[j].ForceSoftening */
#define MODE_B_RADIUS_DEFAULT           (MODE_B_RADIUS_GAS_KERNEL)

/* Legacy aggregation policy for sfc_tiles + ghost_exchange tile cache.
 * Returns P[j].KernelRadius for every type — byte-equivalent to the
 * pre-policy-threading code paths that read P[j].KernelRadius unconditionally.
 * This is the DEFAULT for build_sfc_tiles / gpu_spatial_index_build /
 * gpu_ngb_list_build so non-runner (ghost-exchange) callers see no behavior
 * change.  Runner Mode A passes Spec::radius_policy explicitly instead. */
#define MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES \
    (MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_NONGAS_KERNEL)

/* Conservative-union policy: every leaf-policy source contributes.  Used
 * exclusively by force_hmax_per_type_particle_radius (Mode B node-prune band
 * over the local force-tree, where over-opening a node is bounded perf cost
 * and is local to the host walker).  Do NOT thread this through Mode A
 * SIDX / sfc_tiles / ghost_exchange — it would inflate cached tile bands and
 * cause real BVH/import work explosions on sink/AGS configurations. */
#define MODE_B_RADIUS_ALL_SOURCES \
    (MODE_B_RADIUS_GAS_KERNEL | MODE_B_RADIUS_GAS_AGS \
     | MODE_B_RADIUS_NONGAS_KERNEL | MODE_B_RADIUS_NONGAS_AGS \
     | MODE_B_RADIUS_FORCE_SOFTENING)

/* Field-based form — device-callable.  Returns the j-side symmetric-pair-
 * search radius under `policy`, given j's type and the relevant radius fields.
 * Returns 0 when no policy bit selects a contributing source for j's type
 * (SYMMETRIC then collapses to ONEWAY for that j, which is correct). */
KOKKOS_INLINE_FUNCTION
double nlr_symmetric_radius_from_fields(int type,
                                        double kernel_radius,
                                        double ags_kernel_radius,
                                        double force_softening,
                                        mode_b_radius_policy_t policy)
{
    double h = 0.0;
    if(type == 0) {
        if((policy & MODE_B_RADIUS_GAS_KERNEL) && kernel_radius     > h) h = kernel_radius;
        if((policy & MODE_B_RADIUS_GAS_AGS)    && ags_kernel_radius > h) h = ags_kernel_radius;
    } else {
        if((policy & MODE_B_RADIUS_NONGAS_KERNEL) && kernel_radius     > h) h = kernel_radius;
        if((policy & MODE_B_RADIUS_NONGAS_AGS)    && ags_kernel_radius > h) h = ags_kernel_radius;
    }
    if((policy & MODE_B_RADIUS_FORCE_SOFTENING) && force_softening > h) h = force_softening;
    return h;
}

/* Per-particle wrapper — SSOT for AGS_KernelRadius field gating.
 *
 * Reads KernelRadius / AGS_KernelRadius (when defined) / ForceSoftening / Type
 * off the particle_data reference and forwards to nlr_symmetric_radius_from_fields.
 * All host AND device callers that resolve a per-particle pair-search reach
 * MUST go through this wrapper (or its capped sibling below) — never
 * duplicate the `#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE` guard in caller
 * code.  This is the single compile-flag-dependent decision site for radius
 * policy.  No DMIN(All.MaxKernelRadius) cap is applied; sfc_tiles / Mode A
 * compact_xyzh / Mode B leaf helper all want the uncapped reach. */
KOKKOS_INLINE_FUNCTION
double nlr_particle_symmetric_radius(const struct particle_data &p,
                                     mode_b_radius_policy_t policy)
{
    double ags_kr = 0.0;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    ags_kr = (double)p.AGS_KernelRadius;
#endif
    return nlr_symmetric_radius_from_fields((int)p.Type,
                                            (double)p.KernelRadius,
                                            ags_kr,
                                            (double)p.ForceSoftening,
                                            policy);
}

/* Per-particle wrapper, kernel-radii capped at max_kernel_radius (legacy
 * Extnodes invariant).  ForceSoftening is NOT capped — leaf-policy may admit
 * by FS so the node-prune band must dominate FS unconditionally.
 *
 * Used by force_hmax_per_type_particle_radius (with policy =
 * MODE_B_RADIUS_ALL_SOURCES, max_kernel_radius = All.MaxKernelRadius) to seed
 * the Mode B force-tree per-type bands.  Other callers should prefer the
 * uncapped form unless they explicitly need the legacy MaxKernelRadius cap. */
KOKKOS_INLINE_FUNCTION
double nlr_particle_symmetric_radius_capped(const struct particle_data &p,
                                            mode_b_radius_policy_t policy,
                                            double max_kernel_radius)
{
    double kr = (double)p.KernelRadius;
    if(kr > max_kernel_radius) kr = max_kernel_radius;
    double ags_kr = 0.0;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    ags_kr = (double)p.AGS_KernelRadius;
    if(ags_kr > max_kernel_radius) ags_kr = max_kernel_radius;
#endif
    return nlr_symmetric_radius_from_fields((int)p.Type,
                                            kr,
                                            ags_kr,
                                            (double)p.ForceSoftening,
                                            policy);
}

/* The largest symmetric radius a particle can have under this policy once drifted to the current time,
 * however far behind it is.  A drift (drift_particle_impl) rescales KernelRadius and AGS_KernelRadius by
 * at most `growth` (kernel_radius_drift_max_growth_factor) and then raises them to their floors --
 * `kernel_floor` (All.MinKernelRadius) and `ags_minsoft` (ags_return_minsoft_P); a particle the AGS
 * density does not solve for (ags_density_solves_for_P) gets its softening radius
 * (ForceSoftening_KernelRadius_P); ForceSoftening is not changed.  The AGS inputs are read only when the
 * policy has an AGS bit and AGS is compiled (pass 0 otherwise).  Every bound on a lagging particle's reach
 * goes through this one rule, on host or device. */
KOKKOS_INLINE_FUNCTION
double nlr_symmetric_radius_after_drift_from_fields(int type, double kernel_radius, double force_softening,
                                                    int ags_solves, double ags_kernel_radius,
                                                    double ags_minsoft, double softening_kernel_radius,
                                                    double growth, double kernel_floor,
                                                    mode_b_radius_policy_t policy)
{
    double kr = kernel_radius * growth;
    if(kr < kernel_floor) {kr = kernel_floor;}
    double ags_kr = 0.0;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(policy & (MODE_B_RADIUS_GAS_AGS | MODE_B_RADIUS_NONGAS_AGS)) {
        if(ags_solves) {
            ags_kr = ags_kernel_radius * growth;
            if(ags_kr < ags_minsoft) {ags_kr = ags_minsoft;}
        } else {ags_kr = softening_kernel_radius;}
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        if((1 & ADAPTIVE_GRAVSOFT_FORALL) && (type == 0) && (kr > ags_kr)) {ags_kr = kr;}   /* the gas drift sets it to KernelRadius */
#endif
    }
#else
    (void)ags_solves; (void)ags_kernel_radius; (void)ags_minsoft; (void)softening_kernel_radius;
#endif
    return nlr_symmetric_radius_from_fields(type, kr, ags_kr, force_softening, policy);
}

/* The same for particle j, its inputs read from P_arr, with the growth and floor passed in.  Reads only.
 * One body for host and device.  It is a template so that it can live here, where neither proto.h (which
 * has no include guard) nor ags_functions.h is in scope: its calls take P_arr, so the accessors are found
 * where it is used, and every caller has both. */
template <class ParticleArray>
KOKKOS_INLINE_FUNCTION
double nlr_particle_symmetric_radius_after_drift_P(int j, ParticleArray *P_arr, double growth, double kernel_floor,
                                                   mode_b_radius_policy_t radius_policy)
{
    int ags_solves = 0;
    double ags_kernel_radius = 0.0, ags_minsoft = 0.0, softening_kernel_radius = 0.0;
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(radius_policy & (MODE_B_RADIUS_GAS_AGS | MODE_B_RADIUS_NONGAS_AGS)) {
        ags_solves = ags_density_solves_for_P(j, P_arr);
        ags_kernel_radius = (double)P_arr[j].AGS_KernelRadius;
        ags_minsoft = ags_return_minsoft_P(j, P_arr);
        softening_kernel_radius = ForceSoftening_KernelRadius_P(j, P_arr);
    }
#endif
    return nlr_symmetric_radius_after_drift_from_fields((int)P_arr[j].Type, (double)P_arr[j].KernelRadius,
                                                        (double)P_arr[j].ForceSoftening, ags_solves, ags_kernel_radius,
                                                        ags_minsoft, softening_kernel_radius, growth, kernel_floor,
                                                        radius_policy);
}

/* The host form, with this step's growth and floor.  Defined in mode_b_local_walker.cc. */
double nlr_particle_symmetric_radius_after_drift(int j, struct particle_data *P_arr,
                                                 mode_b_radius_policy_t radius_policy);

#endif /* NLR_RADIUS_POLICY_H */
