/* sfc_tiles.cc — SFC-ordered tile construction (the spatial index itself).
 *
 * Particles are already Peano-Hilbert sorted from domain decomposition.
 * We group them into tiles of ~TILE_TARGET_SIZE particles and compute
 * per-tile bounding boxes and hmax values, then a BVH over those tiles.
 * The traversal that consumes them — tile-level overlap test followed by
 * pairwise distance checks within a tile — lives in sfc_tiles_functions.h.
 *
 * Membership is type mask plus positive mass, and the per-tile fold takes its
 * reach from nlr_particle_symmetric_radius and each member's position at the
 * index's reference time from particle_motion_envelope.  Periodicity and box
 * size enter only at traversal time, which is why the traversal is not in this
 * file.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "sfc_tiles.h"


/* TILE_PERIODIC_X/Y/Z defined in sfc_tiles.h */


/* Supply-pool membership. The single definition — every path that selects
 * pool members tests through this. */
static inline int sfc_pool_member(const struct particle_data *p, int type_bitmask)
{
    if(!((1 << p->Type) & type_bitmask)) return 0;
    if(p->Mass <= 0) return 0;
    return 1;
}

/* Open a tile covering pool slots from `first`, with nothing folded in yet.
 *
 * The box starts INVERTED rather than zeroed. The BVH unions child boxes with
 * min/max, so a zeroed box would drag every ancestor out to the coordinate
 * origin and stop the opener pruning that whole subtree; inverted bounds are
 * neutral under the union and always fail the sphere-overlap gap test, which
 * is what a tile with nothing live in it needs. Folding the first member then
 * sets lo=hi=its position, so no separate "seeded" case is required.  The
 * velocity range starts inverted for the same reason. */
static inline void sfc_tile_begin(sfc_tile_t *tile, int first)
{
    tile->first = first;
    tile->count = 0;
    tile->bvh_leaf = -1;
    tile->hmax = 0;
    for(int t = 0; t < TILE_NUM_PTYPES; t++) tile->hmax_by_type[t] = 0;
    for(int k = 0; k < 3; k++) {
        tile->lo[k] = MAX_REAL_NUMBER; tile->hi[k] = -MAX_REAL_NUMBER;
        tile->u_min[k] = MAX_REAL_NUMBER; tile->u_max[k] = -MAX_REAL_NUMBER;
    }
    tile->rho = 0;
    tile->hw = 0;
}

/* Fold one live member into a tile and write its row. The tiling rule lives here
 * and nowhere else.
 *
 * The row is where the drift puts the member at the reference time ti_ref,
 * which it has not necessarily reached; the tile box holds that position widened
 * by how far from it the member can be (its half-width).  A clock the drift cannot
 * advance to ti_ref, or a motion that is not finite, is a state no search can
 * bound, and the fold refuses it (returns 1).  A member the envelope will not
 * predict because a reflecting or outflow wall may turn it round is bounded
 * instead by its speed from where it stands.
 *
 * The row holds the member's reach at its own clock (the SSOT per-particle reach under radius_policy)
 * and its drifted reach; hmax aggregates the drifted reach. */
static inline int sfc_tile_fold(sfc_tile_t *tile, struct particle_data *P,
                                const struct gas_cell_data *cells, int j,
                                mode_b_radius_policy_t radius_policy,
                                integertime ti_ref, const struct DriftKickTableView *tables,
                                double row[SIDX_ROW_WIDTH], int *all_current)
{
    const struct particle_data *p = &P[j];
    const integertime ti_j = p->Ti_current;
    if(ti_j < 0 || ti_j > ti_ref) {return 1;}
    double center[3], hw = 0.0;
    const int motion = particle_motion_envelope(j, P, cells, ti_ref, tables, center, &hw);
    if(motion == PARTICLE_MOTION_UNBOUNDED) {
        const double dl = motion_bound_widening(particle_motion_speed_bound(j, P, cells), ti_j, ti_ref, tables);
        if(!motion_bound_widening_is_valid(dl)) {return 1;}
        for(int k = 0; k < 3; k++) {
            center[k] = (double)p->Pos[k];
            if(!(center[k] - center[k] == 0.0)) {return 1;}   /* NaN or Inf, fast-math safe */
        }
        hw = 0.5 * dl + motion_envelope_rounding_floor(center, 0.5 * dl);
    }
    if(motion != PARTICLE_MOTION_CURRENT) {*all_current = 0;}
    double u_lo[3], u_hi[3], rho = 0.0;
    if(sfc_member_motion_range(j, P, cells, u_lo, u_hi, &rho)) {return 1;}
    for(int k = 0; k < 3; k++) {
        if(center[k] - hw < tile->lo[k]) tile->lo[k] = center[k] - hw;
        if(center[k] + hw > tile->hi[k]) tile->hi[k] = center[k] + hw;
        if(u_lo[k] < tile->u_min[k]) tile->u_min[k] = u_lo[k];
        if(u_hi[k] > tile->u_max[k]) tile->u_max[k] = u_hi[k];
    }
    if(rho > tile->rho) tile->rho = rho;
    if(hw > tile->hw) tile->hw = hw;
    const double hj = nlr_particle_symmetric_radius(*p, radius_policy);
    double hd = nlr_particle_symmetric_radius_after_drift(j, P, radius_policy);
    if(hd < hj) hd = hj;
    if(hd > tile->hmax) tile->hmax = hd;
    int tj = (int)p->Type;
    if(tj >= 0 && tj < TILE_NUM_PTYPES && hd > tile->hmax_by_type[tj])
        tile->hmax_by_type[tj] = hd;
    row[0] = center[0]; row[1] = center[1]; row[2] = center[2]; row[3] = hj; row[4] = hd;
    return 0;
}

int build_sfc_supply_pool(struct particle_data *P, int num_total,
                          int type_bitmask, int **pool_indices_out)
{
    int num_pool = 0;
    for(int i = 0; i < num_total; i++) {
        if(!sfc_pool_member(&P[i], type_bitmask)) continue;
        num_pool++;
    }
    if(!pool_indices_out) return num_pool;

    int *pool = (int *) mymalloc("sfc_pool", (num_pool > 0 ? num_pool : 1) * sizeof(int));
    int p = 0;
    for(int i = 0; i < num_total; i++) {
        if(!sfc_pool_member(&P[i], type_bitmask)) continue;
        pool[p++] = i;
    }
    *pool_indices_out = pool;
    return num_pool;
}

int build_sfc_tiles(struct particle_data *P, int num_total,
                    int type_bitmask, int target_tile_size,
                    sfc_tile_t **tiles_out, int **pool_indices_out,
                    int *num_pool_out, double **rows_out,
                    integertime ti_ref, const struct DriftKickTableView *tables,
                    int *all_current_out,
                    mode_b_radius_policy_t radius_policy)
{
    /* Deriving the pool and tiling it are the same walk over P[], so do them in
     * ONE pass and avoid several full streams over a large particle struct: this
     * is a bandwidth-bound loop over every particle, and each extra stream is
     * another cache line per particle. Membership and the tiling rule are the
     * shared helpers above, so this states neither of them a second time.
     *
     * The arrays are sized to their upper bounds because the member count is
     * not known until the pass ends. That costs a transient int[num_total] and
     * rows for a full pool; the persistent copies callers keep are cut to the
     * exact counts returned here. Tiles and rows are mymalloc'd above the pool,
     * so the caller frees rows, then tiles, then the pool. */
    int pool_capacity = (num_total > 0) ? num_total : 1;
    int tile_capacity = (num_total + target_tile_size - 1) / target_tile_size;
    if(tile_capacity < 1) tile_capacity = 1;

    int *pool = (int *) mymalloc("sfc_pool", pool_capacity * sizeof(int));
    sfc_tile_t *tiles = (sfc_tile_t *) mymalloc("sfc_tiles", tile_capacity * sizeof(sfc_tile_t));
    double *rows = (double *) mymalloc("sfc_rows", (size_t) pool_capacity * SIDX_ROW_WIDTH * sizeof(double));

    int num_pool = 0, ntiles = 0, all_current = 1;
    for(int i = 0; i < num_total; i++)
    {
        if(!sfc_pool_member(&P[i], type_bitmask)) continue;
        /* A member starting a fresh tile opens it; tiles cover consecutive runs
         * of target_tile_size pool slots, so `first` is the slot it opens at. */
        if((num_pool % target_tile_size) == 0) sfc_tile_begin(&tiles[ntiles++], num_pool);
        tiles[ntiles - 1].count++;
        if(sfc_tile_fold(&tiles[ntiles - 1], P, CellP, i, radius_policy, ti_ref, tables,
                         &rows[(size_t) num_pool * SIDX_ROW_WIDTH], &all_current)) {
            myfree(rows); myfree(tiles); myfree(pool);
            return -1;
        }
        pool[num_pool++] = i;
    }
    /* An empty pool still publishes one (empty, inverted-box) tile so the BVH
     * always has a root to build over. */
    if(ntiles == 0) { sfc_tile_begin(&tiles[0], 0); ntiles = 1; }

    *tiles_out = tiles;
    *pool_indices_out = pool;
    *num_pool_out = num_pool;
    *rows_out = rows;
    *all_current_out = all_current;
    return ntiles;
}


void free_sfc_tiles(sfc_tile_t *tiles, int *pool_indices)
{
    /* Free in reverse mymalloc order: tiles allocated after pool */
    if(tiles) myfree(tiles);
    if(pool_indices) myfree(pool_indices);
}


/* ================================================================
   BVH over SFC tiles — recursive midpoint subdivision.
   Since tiles are SFC-sorted, consecutive tiles are spatially local,
   so midpoint subdivision produces a spatially balanced tree.
   ================================================================ */

/* Recursive builder. Returns the index of the created node in bvh[].
   *next_node is the next free slot. */
static int build_bvh_recursive(sfc_tile_t *tiles, int tile_start, int tile_end,
                               tile_bvh_node_t *bvh, int *next_node)
{
    if(tile_end - tile_start == 1)
    {
        /* Leaf: create a node pointing to a single tile */
        int idx = (*next_node)++;
        int t = tile_start;
        for(int k = 0; k < 3; k++) {
            bvh[idx].lo[k] = tiles[t].lo[k]; bvh[idx].hi[k] = tiles[t].hi[k];
            bvh[idx].u_min[k] = tiles[t].u_min[k]; bvh[idx].u_max[k] = tiles[t].u_max[k];
        }
        bvh[idx].hmax = tiles[t].hmax;
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) bvh[idx].hmax_by_type[tt] = tiles[t].hmax_by_type[tt];
        bvh[idx].rho = tiles[t].rho;
        bvh[idx].left = -(t + 1);   /* negative encoding: leaf = -(tile_index + 1) */
        bvh[idx].right = -(t + 1);  /* same tile for both (signals leaf) */
        bvh[idx].parent = -1;
        tiles[t].bvh_leaf = idx;
        return idx;
    }

    /* Internal node: split at midpoint.  Children are built first so the
     * node's box, reach and velocity range can be their union; that also puts
     * every node's index above its children's. */
    int mid = (tile_start + tile_end) / 2;
    int left_idx = build_bvh_recursive(tiles, tile_start, mid, bvh, next_node);
    int right_idx = build_bvh_recursive(tiles, mid, tile_end, bvh, next_node);

    int idx = (*next_node)++;
    for(int k = 0; k < 3; k++) {
        bvh[idx].lo[k] = DMIN(bvh[left_idx].lo[k], bvh[right_idx].lo[k]);
        bvh[idx].hi[k] = DMAX(bvh[left_idx].hi[k], bvh[right_idx].hi[k]);
        bvh[idx].u_min[k] = DMIN(bvh[left_idx].u_min[k], bvh[right_idx].u_min[k]);
        bvh[idx].u_max[k] = DMAX(bvh[left_idx].u_max[k], bvh[right_idx].u_max[k]);
    }
    bvh[idx].hmax = DMAX(bvh[left_idx].hmax, bvh[right_idx].hmax);
    for(int tt = 0; tt < TILE_NUM_PTYPES; tt++)
        bvh[idx].hmax_by_type[tt] = DMAX(bvh[left_idx].hmax_by_type[tt], bvh[right_idx].hmax_by_type[tt]);
    bvh[idx].rho = DMAX(bvh[left_idx].rho, bvh[right_idx].rho);
    bvh[idx].left = left_idx;
    bvh[idx].right = right_idx;
    bvh[idx].parent = -1;
    bvh[left_idx].parent = idx;
    bvh[right_idx].parent = idx;
    return idx;
}

int build_tile_bvh(sfc_tile_t *tiles, int ntiles, tile_bvh_node_t **bvh_out)
{
    if(ntiles <= 0) { *bvh_out = NULL; return 0; }

    /* A binary tree with ntiles leaves has (2*ntiles - 1) total nodes */
    int max_nodes = 2 * ntiles - 1;
    tile_bvh_node_t *bvh = (tile_bvh_node_t *) mymalloc("tile_bvh", max_nodes * sizeof(tile_bvh_node_t));

    int next_node = 0;
    build_bvh_recursive(tiles, 0, ntiles, bvh, &next_node);

    *bvh_out = bvh;
    return next_node; /* root is at index (next_node - 1) */
}
