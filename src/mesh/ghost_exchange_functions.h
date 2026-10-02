/* ghost_exchange_functions.h -- shared neighbour-accept predicates for discovery.
 *
 * Two predicates, one for points and one for extended objects.  Both take the
 * RAW separation and perform the box wrap THEMSELVES, using only the canonical
 * macro family (NEAREST_XYZ / NGB_PERIODIC_BOX_LONG_{X,Y,Z}).  That family is
 * the single source of truth for box wrapping in this code and in legacy, it
 * takes all three coordinates together, and it is the only thing that knows
 * about periodic / long-box / shearing / reflecting / outflow boundaries.
 *
 * ⛔ NEVER hand-roll a wrap (`if(dx > 0.5*L) dx -= L;`) in a discovery
 * predicate, and never add a periodic_flags/box_sizes-style per-axis API.  A
 * per-axis interface CANNOT express a shearing box, where wrapping in x forces
 * a time-dependent shift in y, so it silently under-includes neighbours near
 * the radial boundary and produces wrong physics with no diagnostic.  Wrap
 * first, then test: the walkers and discovery runners then need to know
 * nothing whatsoever about box geometry.
 *
 * SCOPE: geometry ONLY.  No supply-type mask, no source/destination
 * bookkeeping, no dedup, no send-set emission -- those are caller eligibility
 * concerns and stay outside these helpers (different callers share the accept
 * but differ on supply eligibility).
 *
 * Host- and device-callable.  GX_EXCHANGE_INLINE READS (never defines/mutates)
 * KOKKOS_INLINE_FUNCTION: a GPU TU that included Kokkos upstream gets the real
 * device decoration; a plain-mpicxx host TU (ghost_exchange.cc) gets `inline`
 * and never pulls <Kokkos_Core.hpp>.
 *
 * REQUIRES the including TU to have allvars.h (hence macros.h) visible FIRST,
 * for the wrap macros and the box-geometry globals they read.  Every current
 * includer already does; this mirrors the same requirement nlr_radius_policy.h
 * documents.
 */
#ifndef GIZMO_GHOST_EXCHANGE_FUNCTIONS_H
#define GIZMO_GHOST_EXCHANGE_FUNCTIONS_H

#include "neighbor_list.h"   /* NGB_SEARCH_ONEWAY / NGB_SEARCH_SYMMETRIC */

#ifdef KOKKOS_INLINE_FUNCTION
#define GX_EXCHANGE_INLINE KOKKOS_INLINE_FUNCTION
#else
#define GX_EXCHANGE_INLINE inline
#endif

/* Point-pair accept.  Takes the RAW (unwrapped) separation and does the box
 * wrap itself through the canonical macros, so a caller cannot pass an
 * unwrapped delta into a predicate expecting a wrapped one.  Returns 1 if j is
 * a neighbour of the query, else 0.
 *   d*0             RAW separation, query minus neighbour, same order at every
 *                   call site (never mix the delta order within one test)
 *   h_q             query search radius
 *   h_j             neighbour supply radius (ignored for ONEWAY)
 *   search_mode     NGB_SEARCH_ONEWAY or NGB_SEARCH_SYMMETRIC */
GX_EXCHANGE_INLINE
int gx_pair_accept_wrap_and_test(double dx0, double dy0, double dz0,
                                 double h_q, double h_j, int search_mode)
{
    double xtmp = 0; (void)xtmp;   /* required by the NGB_PERIODIC_BOX_LONG_* macros */
    double adx = NGB_PERIODIC_BOX_LONG_X(dx0, dy0, dz0, 1);
    double ady = NGB_PERIODIC_BOX_LONG_Y(dx0, dy0, dz0, 1);
    double adz = NGB_PERIODIC_BOX_LONG_Z(dx0, dy0, dz0, 1);
    double r2 = adx*adx + ady*ady + adz*adz;
    double thresh;
    if(search_mode == NGB_SEARCH_ONEWAY) { thresh = h_q; }
    else { thresh = (h_q > h_j) ? h_q : h_j; }
    return r2 < thresh * thresh;
}

/* Conservative overlap between a search sphere of radius R and an EXTENDED
 * object (tree node, SFC tile, or the Minkowski sum of two boxes) given as a
 * centre separation plus per-axis half-widths.  Takes the RAW separation and
 * wraps it here, through the canonical macros only.
 *
 * This is a port of the legacy node check (system/ngb_codeblock_checknode.h):
 * per-axis reject at R + halfwidth, then a radial reject against the object's
 * circumscribing sphere.  For a cube of side len the circumradius is
 * 0.866*len, so the cube case is identical to legacy's
 * 0.5*len + CUBE_EDGEFACTOR_1*len; the per-axis form generalises it to the
 * non-cubic tiles and box-box sums the modern code also needs.
 *
 * SHEARING BOXES: an extended object can STRADDLE the radial wrap boundary --
 * part of its extent wraps in y, part does not -- so evaluating the y image at
 * the centre alone is wrong.  Legacy's answer, reproduced here: when the x
 * extent crosses boxHalf_X, evaluate the y image at BOTH x edges and take the
 * minimum, which brackets the wrapped and unwrapped regimes.  A point has zero
 * extent and therefore never straddles, which is why the pair accept above
 * needs no such branch.
 *
 * Requires allvars.h (and hence macros.h) to be included by the TU first, for
 * the wrap macros and boxHalf_X. */
/* THE WRAP, for extended objects.  This is the ONLY place the straddle case is
 * handled; the tests below and at the call sites consume its output and do no
 * wrapping of their own.  `hwx` is the object's x half-extent (for a box pair,
 * the SUM of the two half-extents, since that pair is a point against their
 * Minkowski sum). */
GX_EXCHANGE_INLINE
void gx_extended_wrapped_separation(double dx0, double dy0, double dz0, double hwx,
                                    double *adx_out, double *ady_out, double *adz_out)
{
    double xtmp = 0; (void)xtmp;   /* required by the NGB_PERIODIC_BOX_LONG_* macros */
    *adx_out = NGB_PERIODIC_BOX_LONG_X(dx0, dy0, dz0, -1);
    *adz_out = NGB_PERIODIC_BOX_LONG_Z(dx0, dy0, dz0, -1);
#if (BOX_SHEARING > 1)
    {
        const double dx_abs = fabs(dx0);
        if((dx_abs + hwx < boxHalf_X) || (dx_abs - hwx > boxHalf_X)) {
            /* wholly on one side of the wrap: the centre image is the object's image */
            *ady_out = NGB_PERIODIC_BOX_LONG_Y(dx0, dy0, dz0, -1);
        } else {
            /* Straddles the wrap.  The y image depends on x only through a
             * three-way selector (x below -boxHalf_X, between, or above
             * +boxHalf_X), so sampling one x in every region the object spans
             * gives its exact minimum image, and the minimum is what keeps a
             * reachable neighbour from being pruned.
             *
             * The two x edges cover the outer regions.  A THIRD sample is
             * needed because an object wide enough to contain the whole middle
             * region has both edges outside it, and the edges would then miss
             * the unshifted regime entirely: measured, that overestimates the
             * minimum by up to half a box (worst 0.4995 over 30k random
             * configurations with hwx >= boxHalf_X), i.e. exactly the silent
             * under-inclusion this predicate exists to prevent.  Clamping 0
             * into the x extent lands in the middle region whenever the object
             * reaches it, and duplicates an edge otherwise.
             *
             * Tree nodes have hwx << boxHalf_X and never reach this case, so
             * the legacy two-edge behaviour is unchanged; the third sample is
             * for the tile PAIR, whose hwx is a SUM of half-widths — a case
             * legacy never had. */
            const double x_lo = dx0 - hwx, x_hi = dx0 + hwx;
            const double x_mid = (0.0 < x_lo) ? x_lo : ((0.0 > x_hi) ? x_hi : 0.0);
            const double ady_m = NGB_PERIODIC_BOX_LONG_Y(x_lo,  dy0, dz0, -1);
            const double ady_p = NGB_PERIODIC_BOX_LONG_Y(x_hi,  dy0, dz0, -1);
            const double ady_c = NGB_PERIODIC_BOX_LONG_Y(x_mid, dy0, dz0, -1);
            double amin = (ady_m < ady_p) ? ady_m : ady_p;
            if(ady_c < amin) { amin = ady_c; }
            *ady_out = amin;
        }
    }
#else
    (void)hwx;
    *ady_out = NGB_PERIODIC_BOX_LONG_Y(dx0, dy0, dz0, -1);
#endif
}

/* Sphere-vs-extended-object overlap: the legacy node check
 * (ngb_codeblock_checknode.h) — per-axis reject at R + halfwidth, then a radial
 * reject against the object's circumscribing sphere.  For a cube of side len
 * the circumradius is 0.866*len, so the cube case is identical to legacy's
 * 0.5*len + CUBE_EDGEFACTOR_1*len; the per-axis form generalises it to the
 * non-cubic tiles the modern code also needs.
 *
 * ⛔ DO NOT "tighten" this to the exact clamped per-axis gap.  It looks
 * strictly better (exact point-to-box distance, no sqrt) and it is NOT
 * behaviour-preserving: measured, it moved evrard np=2 off bitwise
 * (Velocities max_rel 1.8e-10).  The tree/BVH search radii are evidently
 * relying on this bound's slack somewhere, so removing it changes which nodes
 * open.  Any such change is a separate, measured perf row — never a rider on a
 * correctness commit. */
GX_EXCHANGE_INLINE
int gx_extended_overlap_wrap_and_test(double dx0, double dy0, double dz0,
                                      double hwx, double hwy, double hwz,
                                      double R)
{
    double adx, ady, adz;
    gx_extended_wrapped_separation(dx0, dy0, dz0, hwx, &adx, &ady, &adz);
    if(adx > R + hwx) { return 0; }
    if(adz > R + hwz) { return 0; }
    if(ady > R + hwy) { return 0; }
    const double r_circ = sqrt(hwx*hwx + hwy*hwy + hwz*hwz);
    const double r_max = R + r_circ;
    return (adx*adx + ady*ady + adz*adz) < r_max * r_max;
}

/* Box-vs-box overlap within a search radius: a point against the Minkowski sum
 * of the two boxes.  Keeps the EXACT clamped per-axis gap this tile pass has
 * always used — it gates ghost-import volume, where over-opening costs real
 * transport, and unlike the node test above it was never on the looser bound.
 * Same wrap, same straddle handling; only the acceptance geometry differs. */
GX_EXCHANGE_INLINE
int gx_boxpair_overlap_wrap_and_test(double dx0, double dy0, double dz0,
                                     double hwx, double hwy, double hwz,
                                     double R, double R2)
{
    double a[3], hw[3];
    gx_extended_wrapped_separation(dx0, dy0, dz0, hwx, &a[0], &a[1], &a[2]);
    hw[0] = hwx; hw[1] = hwy; hw[2] = hwz;
    double dist2 = 0;
    for(int k = 0; k < 3; k++) {
        const double gap = a[k] - hw[k];
        if(gap <= 0) { continue; }        /* boxes overlap on this axis */
        if(gap > R)  { return 0; }
        dist2 += gap * gap;
    }
    return dist2 < R2;
}

#endif /* GIZMO_GHOST_EXCHANGE_FUNCTIONS_H */
