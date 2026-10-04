/* sfc_tiles_functions.h — GPU-callable SFC-tile neighbor search functions.
 *
 * Contains: sfc_box_may_reach(), check_tile_particles_gpu(),
 *   bvh_walk_tiles(), search_neighbors_sfc_gpu().
 *
 * These are the per-particle search functions over the tiles and BVH that
 * sfc_tiles.cc builds. They are the only traversal of that index; sfc_tiles.cc
 * itself constructs it and does not search it.
 * On GPU: the Kokkos kernel includes this with KOKKOS_INLINE_FUNCTION as
 *   __device__ __host__ inline, making these callable from parallel_for.
 *
 * The functions take explicit pointers to tiles/BVH/pool arrays.  Box wrapping
 * is NOT passed in: it belongs to the canonical macro family
 * (NEAREST_XYZ / NGB_PERIODIC_BOX_LONG_*), which these call through the shared
 * predicates in ghost_exchange_functions.h.  A per-axis periodicity argument
 * cannot express a shearing box and is why this path once silently
 * under-included neighbours there.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef SFC_TILES_FUNCTIONS_H
#define SFC_TILES_FUNCTIONS_H

#include "sfc_tiles.h"
#include "ghost_exchange_functions.h"  /* gx_extended_overlap_wrap_and_test: the canonical-wrap SSOT */

/* Can a query of reach R at pos reach anything a node or tile holds at the time of the search?
 * At the index's reference time its members lie in [lo, hi]; an undilated drift interval D later they
 * lie in [lo + u_min D - rho D, hi + u_max D + rho D] (sfc_tiles.h).  An exact frame (D = 0, nothing
 * moved) reads the box as built.  Otherwise the moved box carries a rounding allowance, so a member the
 * exact test accepts is never pruned by arithmetic.
 * Wrapping is the canonical macro family's job (ghost_exchange_functions.h); this takes no box-geometry
 * arguments.  Half-width rounded UP so the conversion from lo/hi stays conservative under FP. */
KOKKOS_INLINE_FUNCTION
int sfc_box_may_reach(const double box_lo[3], const double box_hi[3],
                      const double u_min[3], const double u_max[3], double rho,
                      const struct sfc_walk_frame &frame, const double pos[3], double R)
{
    double c[3], hw[3];
    for(int k = 0; k < 3; k++)
    {
        c[k]  = 0.5 * (box_lo[k] + box_hi[k]);
        double up = box_hi[k] - c[k], dn = c[k] - box_lo[k];
        hw[k] = (up > dn) ? up : dn;
    }
    if(!frame.exact) {
        double scale = R;
        for(int k = 0; k < 3; k++) {
            c[k]  += 0.5 * (u_min[k] + u_max[k]) * frame.D;
            hw[k] += (0.5 * (u_max[k] - u_min[k]) + rho) * frame.D;
            if(hw[k] > scale) {scale = hw[k];}
        }
        const double slack = motion_envelope_test_slack(c, pos, scale);
        for(int k = 0; k < 3; k++) {hw[k] += slack;}
    }
    return gx_extended_overlap_wrap_and_test(c[0] - pos[0], c[1] - pos[1], c[2] - pos[2],
                                             hw[0], hw[1], hw[2], R);
}

/* Process a single tile's particles against an arbitrary source position pos_i.
 * Returns neighbor count. Uses the canonical neighbor periodic macros so
 * shearing/long/reflected/outflow boundary behavior matches the CPU SFC
 * neighbor search. The source identity is purely positional — the function
 * does NOT filter j == anything; callers wanting "skip self" must do so
 * after consuming the CSR list. (Existing callers don't filter; this matches
 * the legacy ngb_treefind_* semantic that returns self when source is a
 * particle and the same particle is in the search pool.) */
/* rows[slot*SIDX_ROW_WIDTH ...] = x,y,z, reach, drifted reach for the member in pool slot `slot` (see
 * sfc_tiles.h; DOUBLE positions: float ABSOLUTE positions are invalid for GIZMO's ~1e11 dynamic range and
 * must not decide neighbour inclusion).  Under an exact frame this is the exact test and the list it
 * builds is exact.  Otherwise each member is tested as the box it can occupy now -- its row moved by its
 * tile's velocity range and residual, widened by its tile's largest half-width at the reference time --
 * against the reach it can have now; that list is a superset, which the list builder drifts and trims to
 * the exact test. */
/* j_radius_scale: multiplier applied to the j-side kernel radius in SYMMETRIC
 * mode (1.0 = legacy behavior). Scaled-symmetric callers (TURB_DIFF_DYNAMIC
 * wide-filter loops) pass All.TurbDynamicDiffFac so the pair reach becomes
 * max(h_i, fac*h_j). Applied at query time; the stored rows and bands stay
 * keyed on raw radii. */
KOKKOS_INLINE_FUNCTION
int check_tile_particles_gpu(const double *rows, const double pos_i[3], double h_i, double h2_i,
                             double j_radius_scale,
                             const sfc_tile_t *tile, const int *pool, int search_mode,
                             const struct sfc_walk_frame &frame,
                             int *store_neighbors, int count, int max_store,
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
    MyDouble xtmp = 0; /* required by NGB_PERIODIC_BOX_LONG_* macros */
    for(int s = 0; s < tile->count; s++)
    {
        if(cnt_candidates_tested) (*cnt_candidates_tested)++;
        const int slot = tile->first + s;
        int j = pool[slot];
        if(j < 0) continue;   /* a member that left the pool in place since the build (gpu_sidx_notify_member_lost) */
        /* Per-type supply filter at leaf. No-op when supply_mask == 0x3f (default).
         * Constant-propagated away when caller passes default args. */
        if(P_gpu) {
            int pt = (int)P_gpu[j].Type;
            if(pt < 0 || pt >= 6) continue;
            if((supply_mask & (1u << (unsigned)pt)) == 0u) continue;
        }
        const double *row = &rows[(size_t)slot * SIDX_ROW_WIDTH];
        int accept;
        if(frame.exact) {
            double dx_raw = pos_i[0] - row[0];
            double dy_raw = pos_i[1] - row[1];
            double dz_raw = pos_i[2] - row[2];
            double adx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
            double ady = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
            double adz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);

            double pair_search_r2;
            if(search_mode == NGB_SEARCH_ONEWAY) {
                pair_search_r2 = h2_i;
            } else {
                double h_j = row[3] * j_radius_scale;
                double h_max = (h_i > h_j) ? h_i : h_j;
                pair_search_r2 = h_max * h_max;
            }

            if(adx > h_i && (search_mode == NGB_SEARCH_ONEWAY || adx * adx > pair_search_r2)) continue;
            double r2 = adx * adx + ady * ady + adz * adz;
            accept = (r2 < pair_search_r2);
        } else {
            double reach = h_i;
            if(search_mode != NGB_SEARCH_ONEWAY) {
                const double h_j = (frame.reach_current ? row[3] : row[4]) * j_radius_scale;
                if(h_j > reach) {reach = h_j;}
            }
            double c[3], w[3], scale = reach;
            for(int k = 0; k < 3; k++) {
                c[k] = row[k] + 0.5 * (tile->u_min[k] + tile->u_max[k]) * frame.D;
                w[k] = (0.5 * (tile->u_max[k] - tile->u_min[k]) + tile->rho) * frame.D + tile->hw;
                if(w[k] > scale) {scale = w[k];}
            }
            const double slack = motion_envelope_test_slack(c, pos_i, scale);
            accept = gx_boxpair_overlap_wrap_and_test(c[0] - pos_i[0], c[1] - pos_i[1], c[2] - pos_i[2],
                                                      w[0] + slack, w[1] + slack, w[2] + slack,
                                                      reach, reach * reach);
        }
        if(accept) {
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


/* Walk the tile BVH for one query and hand every overlapping tile to a visitor.
 *
 * This is the one traversal of the tile index.  It decides which tiles a query
 * reaches -- opening nodes on the query reach (widened by the node's supply
 * reach in SYMMETRIC mode), wrapping through the canonical box macros, stepping
 * an explicit stack so it runs on the device -- and nothing else.  What happens
 * at a tile is the visitor's business: the neighbour-list build stores candidate
 * indices, a fused loop evaluates its pair kernel on the tile's particles.
 * Keeping those apart lets a new consumer reuse the opener instead of copying
 * it, and a copied opener is where the wrap convention and the opening rule
 * silently diverge.
 *
 * `visit(tile_index)` is called once per overlapping tile, in traversal order.
 * A negative root is an EMPTY index (no tiles) and walks nothing: a query into
 * an empty pool has no neighbours, and must not read a node that was never
 * built.
 *
 * Every box is read at the time of the search through `frame` (sfc_box_may_reach)
 * and every band is a drifted reach, so an index whose members have drifted
 * or been kicked since it was built still finds everything.  The walk writes
 * nothing to the index. */
template <class TileVisitor>
KOKKOS_INLINE_FUNCTION
void bvh_walk_tiles(const double pos_i[3], double h_i, double j_radius_scale,
                    int search_mode,
                    const tile_bvh_node_t *bvh, int bvh_root,
                    const struct particle_data *P_gpu, unsigned int supply_mask,
                    int *cnt_nodes_visited, int *cnt_tiles_visited,
                    TileVisitor &visit, const struct sfc_walk_frame &frame)
{
    if(bvh_root < 0) {return;}

    int stack[TILE_BVH_STACK_SIZE];
    int sp = 0;
    stack[sp++] = bvh_root;

    while(sp > 0)
    {
        int node_idx = stack[--sp];
        const tile_bvh_node_t *node = &bvh[node_idx];
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
        /* The drifted band, scaled on the j side in SYMMETRIC mode
         * (no-op when j_radius_scale == 1.0; constant-propagated for default callers). */
        node_hmax_eff *= j_radius_scale;
        double search_r = (search_mode == NGB_SEARCH_ONEWAY) ? h_i : ((h_i > node_hmax_eff) ? h_i : node_hmax_eff);

        /* Check if node's box, as it can be now, is within search_r of particle i */
        if(!sfc_box_may_reach(node->lo, node->hi, node->u_min, node->u_max, node->rho, frame, pos_i, search_r)) continue;

        if(node->left < 0)
        {
            /* Leaf node: the visitor takes the tile */
            int tile_idx = -(node->left + 1);
            if(cnt_tiles_visited) (*cnt_tiles_visited)++;
            visit(tile_idx);
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
}

/* The neighbour-list visitor: test every particle of an overlapping tile against
 * the query and store the indices that pass, bounded by max_store.  This is the
 * candidate-list form the CSR build consumes. */
struct TileCandidateStore {
    const double *rows;
    const double *pos_i;
    double h_i, h2_i, j_radius_scale;
    const sfc_tile_t *tiles;
    const int *pool;
    int search_mode;
    struct sfc_walk_frame frame;
    int *store_neighbors;
    int max_store;
    int *cnt_candidates_tested;
    int *cnt_candidates_accepted;
    const struct particle_data *P_gpu;
    unsigned int supply_mask;
    int count;

    KOKKOS_INLINE_FUNCTION
    void operator()(int tile_idx)
    {
        count = check_tile_particles_gpu(rows, pos_i, h_i, h2_i, j_radius_scale, &tiles[tile_idx], pool, search_mode,
                                         frame, store_neighbors, count, max_store,
                                         cnt_candidates_tested, cnt_candidates_accepted,
                                         P_gpu, supply_mask);
    }
};

/* Search for neighbors of an arbitrary source position pos_i using BVH
 * traversal over tiles. Returns count; if store_neighbors != NULL, writes indices
 * there. Decoupled from any specific P[] index so the same routine serves
 * particle-based sources and arbitrary-position sources (e.g.
 * TURB_DRIVING_SPECTRUMGRID grid cells). */
/* j_radius_scale: see check_tile_particles_gpu — scales the j-side radius
 * (per-particle h_j at the leaf, per-node hmax at the BVH opener) in
 * SYMMETRIC mode. 1.0 = legacy. */
KOKKOS_INLINE_FUNCTION
int search_neighbors_sfc_gpu(const double *rows, const double pos_i[3], double h_i,
                             double j_radius_scale,
                             const sfc_tile_t *tiles, int ntiles,
                             const int *pool, int search_mode,
                             const tile_bvh_node_t *bvh, int bvh_root,
                             const struct sfc_walk_frame &frame,
                             int *store_neighbors, int max_store,
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
    TileCandidateStore store{rows, pos_i, h_i, h_i * h_i, j_radius_scale,
                             tiles, pool, search_mode, frame, store_neighbors, max_store,
                             cnt_candidates_tested, cnt_candidates_accepted,
                             P_gpu, supply_mask, 0};
    bvh_walk_tiles(pos_i, h_i, j_radius_scale, search_mode, bvh, bvh_root,
                   P_gpu, supply_mask, cnt_nodes_visited, cnt_tiles_visited, store, frame);
    return store.count;
}

#endif /* SFC_TILES_FUNCTIONS_H */
