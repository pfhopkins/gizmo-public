#ifndef GRAVTREE_MOMENT_SOURCES_H
#define GRAVTREE_MOMENT_SOURCES_H

/*! \file gravtree_moment_sources.h
 *  \brief Single source of truth for the per-particle gravity-tree source INPUT
 *         physics gates: RT source luminosity, sink bolometric luminosity +
 *         emission-angle vector, and cosmic-ray injection rate.
 *
 *  These per-particle quantities feed the GPU gravity walk and the node-moment
 *  build.  The gates used to be duplicated in each precompute venue (the primary
 *  GPU walk in gpu_gravtree.cc, the node-moment refresh in gpu_moment_refresh.cc,
 *  and the LET leaf synthesis in let_pack.cc); this header is the one place the
 *  gating + payload packing lives.
 *
 *  ONE canonical helper, decoration keyed on the including TU:
 *   - host / default TUs: static inline; calls the host wrappers sink_lum_bol /
 *     rt_get_source_luminosity[_chimes] / cr_get_source_injection_rate (external
 *     symbols via proto.h).  This is the eager precompute path (moment refresh,
 *     LET pack, and the dense prefill).
 *   - a TU that exposes the device-callable source cores (defines
 *     GRAVTREE_SOURCE_DEVICE_TU && GRAVTREE_SOURCE_LAZY_SUPPORTED -- only
 *     gpu_gravtree.cc): KOKKOS_INLINE_FUNCTION; calls sink_lum_bol_core and the
 *     inline device dispatchers, so the SAME helper evaluates the source payload
 *     on-device at each local particle-open with no dense O(NumPart) staging.
 *     No forked gate/packing logic.
 *
 *  A validity flag per section (rt_active / bh_active) tells the caller the array
 *  payloads were filled; the caller uses them only when the flag is set
 *  (otherwise the payload fields are indeterminate).  No physics change.
 *
 *  Include this AFTER allvars.h and proto.h; in a device TU also AFTER the source
 *  *_functions.h cores (so the inline device bodies are visible).  proto.h and
 *  sinks/sink.h have no include guards, so we do not re-include them here;
 *  sink_lum_bol() is forward-declared below for the host path.
 */

#ifdef SINK_PHOTONMOMENTUM
double sink_lum_bol(double mdot, double mass, long pindex);   /* canonical decl in sinks/sink.h (unguarded) */
#endif

/* Per-particle source inputs.  Field types match the venue storage / wire types
 * (MyFloat / Vec3<MyFloat> / double) so the venue copy-out is a plain assignment
 * with no extra conversion.  rt_active / bh_active flag whether the array
 * payloads were filled; the caller copies those out only when the flag is set
 * (otherwise the payload fields are indeterminate and must not be read). */
struct gravtree_source_inputs_t {
#ifdef RT_USE_GRAVTREE
    int     rt_active;                       /* nonzero if the particle is an active RT source */
    MyFloat src_lum[N_RT_FREQ_BINS];         /* valid only when rt_active */
#ifdef CHIMES_STELLAR_FLUXES
    double  src_lum_G0[CHIMES_LOCAL_UV_NBINS];  /* valid only when rt_active */
    double  src_lum_ion[CHIMES_LOCAL_UV_NBINS]; /* valid only when rt_active */
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    int           bh_active;                 /* nonzero if a live emitting sink (gates bh_lum/bh_angle) */
    MyFloat       bh_lum;                    /* bolometric sink luminosity; valid only when bh_active */
    Vec3<MyFloat> bh_angle;                  /* emission-angle vector;        valid only when bh_active */
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyFloat cr_inject;                       /* cosmic-ray source injection rate (0 for non-sources) */
#endif
};

/* Decoration + sink-luminosity binding keyed on whether this TU opens the device
 * source cores.  Device TU (gpu_gravtree.cc): KOKKOS_INLINE_FUNCTION + the device
 * core sink_lum_bol_core (host+device-callable).  Host/default TUs: static inline +
 * the host wrapper sink_lum_bol.  RT/CR dispatchers keep their public names in both
 * cases (device-inline where the bodies are open, host-external otherwise). */
#if defined(GRAVTREE_SOURCE_DEVICE_TU) && defined(GRAVTREE_SOURCE_LAZY_SUPPORTED)
#define GRAVTREE_FILL_DECOR                       KOKKOS_INLINE_FUNCTION
#define GRAVTREE_FILL_SINK_LUM_BOL(mdot,mass,p,P) sink_lum_bol_core((mdot),(mass),(p),(P))
#else
#define GRAVTREE_FILL_DECOR                       static inline
#define GRAVTREE_FILL_SINK_LUM_BOL(mdot,mass,p,P) sink_lum_bol((mdot),(mass),(p))
#endif

/* Fill *out using the EXACT gates the live gravity precomputes use.  Only the
 * cheap validity scalars are initialised unconditionally; array payloads are
 * written only inside an active section, and the caller mirrors that by using
 * them only when the flag is set.  No physics change relative to the legacy
 * per-venue loops.
 *
 * These payloads are computed from each source's stored state as it stands when this
 * is called: nothing here predicts a drifted field forward.  A node aggregate is filled
 * at a tree build or moment refresh and then stays fixed until the next one, as in
 * legacy (the node drift moves centres, never payloads), while a leaf opened later
 * evaluates this afresh -- so time-dependent payloads such as a star's age-dependent
 * luminosity can differ between a node and its own leaves by the interval since the
 * build.  The payloads that depend on drift-changed gas state are the gas emission
 * terms under RADTRANSFER (e.g. free-free, from the gas energy, density and mass).
 * They are read stale for an inactive gas source by contract: they feed only the
 * approximate long-range transport methods (OTVET, LEBRON), which accept somewhat
 * stale luminosities by construction.  A use that needs them current must predict
 * them itself rather than read them from here. */
GRAVTREE_FILL_DECOR void
gravtree_fill_particle_source_inputs(int p, struct particle_data *P_arr,
                                     struct gas_cell_data *CellP_arr,
                                     struct gravtree_source_inputs_t *out)
{
#ifdef RT_USE_GRAVTREE
    out->rt_active = 0;
    if(P_arr[p].Mass > 0) {
        double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
        double lum_G0[CHIMES_LOCAL_UV_NBINS], lum_ion[CHIMES_LOCAL_UV_NBINS];
        out->rt_active = rt_get_source_luminosity_chimes(p, 1, lum, lum_G0, lum_ion, P_arr, CellP_arr);
#else
        out->rt_active = rt_get_source_luminosity(p, 1, lum, P_arr, CellP_arr);
#endif
        if(out->rt_active) {
            int kf;
            for(kf = 0; kf < N_RT_FREQ_BINS; kf++) {out->src_lum[kf] = (MyFloat) lum[kf];}
#ifdef CHIMES_STELLAR_FLUXES
            for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++) {
                out->src_lum_G0[kf]  = lum_G0[kf];
                out->src_lum_ion[kf] = lum_ion[kf];
            }
#endif
        }
    }
#endif /* RT_USE_GRAVTREE */

#ifdef SINK_PHOTONMOMENTUM
    out->bh_active = 0;
    if(P_arr[p].Type == 5 && P_arr[p].Mass > 0 &&
       P_arr[p].DensityAroundParticle > 0 && P_arr[p].Sink_Mdot > 0) {
        out->bh_active = 1;     /* gated in regardless of the luminosity value (which may be 0) */
        out->bh_lum    = (MyFloat) GRAVTREE_FILL_SINK_LUM_BOL(P_arr[p].Sink_Mdot, P_arr[p].Sink_Mass, p, P_arr);
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
        out->bh_angle = P_arr[p].Sink_Specific_AngMom;
#else
        out->bh_angle = P_arr[p].GradRho;
#endif
    }
#endif /* SINK_PHOTONMOMENTUM */

#ifdef COSMIC_RAY_SUBGRID_LEBRON
    /* nonzero only for Type 4/5 sources; gas + degenerate (Mass<=0) slots return 0.
     * Mirror the legacy CPU tree, which evaluates every particle ungated. */
    out->cr_inject = 0;
    if(P_arr[p].Mass > 0) {
        out->cr_inject = (MyFloat) cr_get_source_injection_rate(p, P_arr, CellP_arr);
    }
#endif /* COSMIC_RAY_SUBGRID_LEBRON */
}

#undef GRAVTREE_FILL_DECOR
#undef GRAVTREE_FILL_SINK_LUM_BOL

#endif /* GRAVTREE_MOMENT_SOURCES_H */
