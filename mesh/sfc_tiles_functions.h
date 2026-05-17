/* sfc_tiles_functions.h — GPU-callable SFC-tile neighbor search functions.
 *
 * Contains: bbox_overlaps_sphere(), check_tile_particles_gpu(),
 *   search_neighbors_sfc_gpu().
 *
 * These are the per-particle search functions extracted from sfc_tiles.cc.
 * On CPU: sfc_tiles.cc includes this with KOKKOS_INLINE_FUNCTION redefined
 *   to empty (or just uses its own static copies directly).
 * On GPU: the Kokkos kernel includes this with KOKKOS_INLINE_FUNCTION as
 *   __device__ __host__ inline, making these callable from parallel_for.
 *
 * The functions take explicit pointers to tiles/BVH/pool arrays (no globals)
 * and explicit periodicity parameters, making them pure functions suitable
 * for GPU execution.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef SFC_TILES_FUNCTIONS_H
#define SFC_TILES_FUNCTIONS_H

#include "sfc_tiles.h"

/* Check if a search sphere overlaps an axis-aligned bounding box (with periodic wrapping).
 * Returns 1 if overlap, 0 otherwise.
 * Uses center+halfwidth approach for correct periodic distance to AABB.
 * periodic_flags[k]: 1 if axis k is periodic, 0 otherwise.
 * box_sizes[k]: box size along axis k (only used if periodic). */
KOKKOS_INLINE_FUNCTION
int bbox_overlaps_sphere_gpu(const double box_lo[3], const double box_hi[3],
                             const double pos[3], double search_r, double search_r2,
                             const int periodic_flags[3], const double box_sizes[3])
{
    double dist2 = 0;
    for(int k = 0; k < 3; k++)
    {
        double bsize = box_sizes[k];

        /* Distance from pos to center of bbox, with periodic wrapping */
        double center = 0.5 * (box_lo[k] + box_hi[k]);
        double halfwidth = 0.5 * (box_hi[k] - box_lo[k]);
        double dx = fabs(pos[k] - center);
        if(periodic_flags[k] && dx > 0.5 * bsize) dx = bsize - dx;

        /* Gap = distance from pos to nearest edge of bbox */
        double gap = dx - halfwidth;
        if(gap <= 0) continue; /* pos is inside bbox on this axis */

        if(gap > search_r) return 0; /* early rejection */
        dist2 += gap * gap;
    }
    return (dist2 < search_r2);
}


/* Process a single tile's particles against an arbitrary source position pos_i.
 * Returns neighbor count. Uses the canonical neighbor periodic macros so
 * shearing/long/reflected/outflow boundary behavior matches the CPU SFC
 * neighbor search. The source identity is purely positional — the function
 * does NOT filter j == anything; callers wanting "skip self" must do so
 * after consuming the CSR list. (Existing callers don't filter; this matches
 * the legacy ngb_treefind_* semantic that returns self when source is a
 * particle and the same particle is in the search pool.) */
/* compact_xyzh[j*4+0..3] = x,y,z,h for particle j (float, ~32MB for 2M parts).
 * Using compact float array instead of P[j].Pos/KernelRadius reduces BVH
 * traversal memory traffic by ~12× (32 vs 400 bytes/particle), allowing the
 * array to fit in GPU L2 cache and eliminate random-access cache misses. */
/* j_radius_scale: multiplier applied to the j-side kernel radius in SYMMETRIC
 * mode (1.0 = legacy behavior). Scaled-symmetric callers (TURB_DIFF_DYNAMIC
 * wide-filter loops) pass All.TurbDynamicDiffFac so the pair reach becomes
 * max(h_i, fac*h_j) — see OPEN_3d_difffilter_design.md §3. Applied at query
 * time; the cached compact_xyzh / SIDX hmax stay keyed on raw radii. */
KOKKOS_INLINE_FUNCTION
int check_tile_particles_gpu(const float *compact_xyzh, const double pos_i[3], double h_i, double h2_i,
                             double j_radius_scale,
                             sfc_tile_t *tile, int *pool, int search_mode,
                             int *store_neighbors, int count, int max_store,
                             const double box_sizes[3], const double box_halves[3],
                             const int periodic_flags[3],
                             /* Optional per-active counters; pass nullptr for fast-path. Compiler
                              * eliminates the null-checks via constant prop when null literal is passed. */
                             int *cnt_candidates_tested = nullptr,
                             int *cnt_candidates_accepted = nullptr,
                             /* Per-type supply filter (per-type hmax change). When supply_mask == 0x3f
                              * (all types) and P_gpu == nullptr, this is a no-op: every particle
                              * passes the type filter. Caller passes supply_mask narrower than ALL
                              * to skip non-supply types at leaf-accept time. P_gpu must be non-null
                              * if supply_mask is narrower than ALL. */
                             const struct particle_data *P_gpu = nullptr,
                             unsigned int supply_mask = ((1u << 6) - 1u))
{
    (void)box_sizes;
    (void)box_halves;
    (void)periodic_flags;
    MyDouble xtmp = 0; /* required by NGB_PERIODIC_BOX_LONG_* macros */
    for(int s = 0; s < tile->count; s++)
    {
        if(cnt_candidates_tested) (*cnt_candidates_tested)++;
        int j = pool[tile->first + s];
        /* Per-type supply filter at leaf. No-op when supply_mask == 0x3f (default).
         * Constant-propagated away when caller passes default args. */
        if(P_gpu) {
            int pt = (int)P_gpu[j].Type;
            if(pt < 0 || pt >= 6) continue;
            if((supply_mask & (1u << (unsigned)pt)) == 0u) continue;
        }
        double dx_raw = pos_i[0] - (double)compact_xyzh[j*4+0];
        double dy_raw = pos_i[1] - (double)compact_xyzh[j*4+1];
        double dz_raw = pos_i[2] - (double)compact_xyzh[j*4+2];
        double adx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
        double ady = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
        double adz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);

        double pair_search_r2;
        if(search_mode == NGB_SEARCH_ONEWAY) {
            pair_search_r2 = h2_i;
        } else {
            double h_j = (double)compact_xyzh[j*4+3] * j_radius_scale;
            double h_max = (h_i > h_j) ? h_i : h_j;
            pair_search_r2 = h_max * h_max;
        }

        if(adx > h_i && (search_mode == NGB_SEARCH_ONEWAY || adx * adx > pair_search_r2)) continue;
        double r2 = adx * adx + ady * ady + adz * adz;
        if(r2 < pair_search_r2) {
            /* Bounded write: count past max_store still increments (so caller can
             * detect overflow), but the write is suppressed. Used by the fused
             * single-pass build with a per-particle scratchpad of stride max_store. */
            if(store_neighbors && count < max_store) store_neighbors[count] = j;
            count++;
            if(cnt_candidates_accepted) (*cnt_candidates_accepted)++;
        }
    }
    return count;
}


/* Search for neighbors of an arbitrary source position pos_i using BVH
 * traversal over tiles. Iterative stack-based traversal (GPU-friendly,
 * no recursion). Returns count; if store_neighbors != NULL, writes indices
 * there. Decoupled from any specific P[] index so the same routine serves
 * particle-based sources (gpu_ngb_list_build's default mode) and
 * arbitrary-position sources (e.g. TURB_DRIVING_SPECTRUMGRID grid cells). */
/* j_radius_scale: see check_tile_particles_gpu — scales the j-side radius
 * (per-particle h_j at the leaf, per-node hmax at the BVH opener) in
 * SYMMETRIC mode. 1.0 = legacy. */
KOKKOS_INLINE_FUNCTION
int search_neighbors_sfc_gpu(const float *compact_xyzh, const double pos_i[3], double h_i,
                             double j_radius_scale,
                             sfc_tile_t *tiles, int ntiles,
                             int *pool, int search_mode,
                             tile_bvh_node_t *bvh, int bvh_root,
                             int *store_neighbors, int max_store,
                             const int periodic_flags[3],
                             const double box_sizes[3],
                             const double box_halves[3],
                             /* Optional per-active counters for diagnostics. Pass nullptr (default) for
                              * fast path. Compiler eliminates null-checks via constant prop. */
                             int *cnt_nodes_visited = nullptr,
                             int *cnt_tiles_visited = nullptr,
                             int *cnt_candidates_tested = nullptr,
                             int *cnt_candidates_accepted = nullptr,
                             /* Per-type hmax narrowing (per-type hmax change). Caller passes a
                              * narrow supply_mask + P_gpu to make the opener compute search_r from
                              * max-over-supply-types of node->hmax_by_type[t] (not the scalar hmax).
                              * Default (P_gpu=nullptr, supply_mask=0x3f) preserves the legacy
                              * scalar-hmax opener behavior — used by existing call sites that walk
                              * a tree whose pool is already supply-mask filtered. */
                             const struct particle_data *P_gpu = nullptr,
                             unsigned int supply_mask = ((1u << 6) - 1u))
{
    (void)ntiles;
    double h2_i = h_i * h_i;
    int count = 0;

    /* Stack-based iterative BVH traversal */
    int stack[TILE_BVH_STACK_SIZE];
    int sp = 0;
    stack[sp++] = bvh_root;

    while(sp > 0)
    {
        int node_idx = stack[--sp];
        tile_bvh_node_t *node = &bvh[node_idx];
        if(cnt_nodes_visited) (*cnt_nodes_visited)++;

        /* Compute search radius for this node. Per-type hmax filter when P_gpu given;
         * else scalar hmax (legacy behavior). For supply_mask = 0x3f the per-type
         * branch evaluates max over all 6 types == scalar hmax → identical result. */
        double node_hmax_eff;
        if(P_gpu) {
            node_hmax_eff = 0;
            for(int t = 0; t < 6; t++) {
                if((supply_mask & (1u << (unsigned)t)) == 0u) continue;
                if(node->hmax_by_type[t] > node_hmax_eff) node_hmax_eff = node->hmax_by_type[t];
            }
        } else {
            node_hmax_eff = node->hmax;
        }
        /* Scale the j-side node radius in SYMMETRIC mode (no-op when
         * j_radius_scale == 1.0; constant-propagated for default callers). */
        node_hmax_eff *= j_radius_scale;
        double search_r = (search_mode == NGB_SEARCH_ONEWAY) ? h_i : ((h_i > node_hmax_eff) ? h_i : node_hmax_eff);
        double search_r2 = search_r * search_r;

        /* Check if node's bbox (expanded by search_r) overlaps particle i */
        if(!bbox_overlaps_sphere_gpu(node->lo, node->hi, pos_i, search_r, search_r2,
                                     periodic_flags, box_sizes)) continue;

        /* Check if leaf */
        if(node->left < 0)
        {
            /* Leaf node: process the tile's particles */
            int tile_idx = -(node->left + 1);
            if(cnt_tiles_visited) (*cnt_tiles_visited)++;
            count = check_tile_particles_gpu(compact_xyzh, pos_i, h_i, h2_i, j_radius_scale, &tiles[tile_idx], pool, search_mode,
                                             store_neighbors, count, max_store,
                                             box_sizes, box_halves, periodic_flags,
                                             cnt_candidates_tested, cnt_candidates_accepted,
                                             P_gpu, supply_mask);
        }
        else
        {
            /* Internal node: push children onto stack */
            if(sp + 2 > TILE_BVH_STACK_SIZE) {
                /* Stack overflow — should not happen with STACK_SIZE=64 */
                break;
            }
            stack[sp++] = node->left;
            stack[sp++] = node->right;
        }
    }
    return count;
}

#endif /* SFC_TILES_FUNCTIONS_H */
