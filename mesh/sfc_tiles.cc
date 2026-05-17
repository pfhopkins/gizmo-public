/* sfc_tiles.cc — SFC-ordered tile construction and tile-based neighbor search.
 *
 * Particles are already Peano-Hilbert sorted from domain decomposition.
 * We group them into tiles of ~TILE_TARGET_SIZE particles and compute
 * per-tile bounding boxes and hmax values. Neighbor search then checks
 * tile-level overlap before doing pairwise distance checks within tiles.
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


int build_sfc_tiles(struct particle_data *P, int num_total,
                    int type_bitmask, int target_tile_size,
                    sfc_tile_t **tiles_out, int **pool_indices_out,
                    int *num_pool_out)
{
    /* Step 1: Build pool index array (particles matching type/mass filter, in P[] order).
       Since P[] is already SFC-sorted, pool_indices preserves SFC order. */
    int num_pool = 0;
    int i;
    for(i = 0; i < num_total; i++) {
        if(!((1 << P[i].Type) & type_bitmask)) continue;
        if(P[i].Mass <= 0) continue;
        num_pool++;
    }

    int *pool = (int *) mymalloc("sfc_pool", (num_pool > 0 ? num_pool : 1) * sizeof(int));
    int p = 0;
    for(i = 0; i < num_total; i++) {
        if(!((1 << P[i].Type) & type_bitmask)) continue;
        if(P[i].Mass <= 0) continue;
        pool[p++] = i;
    }

    /* Step 2: Compute number of tiles */
    int ntiles = (num_pool + target_tile_size - 1) / target_tile_size;
    if(ntiles < 1) ntiles = 1;

    sfc_tile_t *tiles = (sfc_tile_t *) mymalloc("sfc_tiles", ntiles * sizeof(sfc_tile_t));

    /* Step 3: Build tiles — single pass over pool */
    int t;
    for(t = 0; t < ntiles; t++)
    {
        int start = t * target_tile_size;
        int count = target_tile_size;
        if(start + count > num_pool) count = num_pool - start;

        tiles[t].first = start;
        tiles[t].count = count;
        tiles[t].hmax = 0;
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) tiles[t].hmax_by_type[tt] = 0;

        if(count <= 0) {
            for(int k = 0; k < 3; k++) { tiles[t].lo[k] = 0; tiles[t].hi[k] = 0; }
            continue;
        }

        /* Initialize bbox from first particle */
        int j0 = pool[start];
        tiles[t].lo[0] = tiles[t].hi[0] = P[j0].Pos[0];
        tiles[t].lo[1] = tiles[t].hi[1] = P[j0].Pos[1];
        tiles[t].lo[2] = tiles[t].hi[2] = P[j0].Pos[2];
        if(P[j0].KernelRadius > tiles[t].hmax) tiles[t].hmax = P[j0].KernelRadius;
        {
            int t0 = (int)P[j0].Type;
            if(t0 >= 0 && t0 < TILE_NUM_PTYPES &&
               P[j0].KernelRadius > tiles[t].hmax_by_type[t0])
                tiles[t].hmax_by_type[t0] = P[j0].KernelRadius;
        }

        /* Expand bbox with remaining particles */
        int s;
        for(s = 1; s < count; s++)
        {
            int j = pool[start + s];
            int k;
            for(k = 0; k < 3; k++) {
                if(P[j].Pos[k] < tiles[t].lo[k]) tiles[t].lo[k] = P[j].Pos[k];
                if(P[j].Pos[k] > tiles[t].hi[k]) tiles[t].hi[k] = P[j].Pos[k];
            }
            if(P[j].KernelRadius > tiles[t].hmax) tiles[t].hmax = P[j].KernelRadius;
            int tj = (int)P[j].Type;
            if(tj >= 0 && tj < TILE_NUM_PTYPES &&
               P[j].KernelRadius > tiles[t].hmax_by_type[tj])
                tiles[t].hmax_by_type[tj] = P[j].KernelRadius;
        }
    }

    *tiles_out = tiles;
    *pool_indices_out = pool;
    *num_pool_out = num_pool;
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
        for(int k = 0; k < 3; k++) { bvh[idx].lo[k] = tiles[t].lo[k]; bvh[idx].hi[k] = tiles[t].hi[k]; }
        bvh[idx].hmax = tiles[t].hmax;
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) bvh[idx].hmax_by_type[tt] = tiles[t].hmax_by_type[tt];
        bvh[idx].left = -(t + 1);   /* negative encoding: leaf = -(tile_index + 1) */
        bvh[idx].right = -(t + 1);  /* same tile for both (signals leaf) */
        return idx;
    }

    /* Internal node: split at midpoint */
    int mid = (tile_start + tile_end) / 2;

    /* Reserve this node's slot first (ensures parent index > children for root-last ordering) */
    /* Actually, build children first, then this node, so we can compute union bbox */
    int left_idx = build_bvh_recursive(tiles, tile_start, mid, bvh, next_node);
    int right_idx = build_bvh_recursive(tiles, mid, tile_end, bvh, next_node);

    int idx = (*next_node)++;
    for(int k = 0; k < 3; k++) {
        bvh[idx].lo[k] = DMIN(bvh[left_idx].lo[k], bvh[right_idx].lo[k]);
        bvh[idx].hi[k] = DMAX(bvh[left_idx].hi[k], bvh[right_idx].hi[k]);
    }
    bvh[idx].hmax = DMAX(bvh[left_idx].hmax, bvh[right_idx].hmax);
    for(int tt = 0; tt < TILE_NUM_PTYPES; tt++)
        bvh[idx].hmax_by_type[tt] = DMAX(bvh[left_idx].hmax_by_type[tt], bvh[right_idx].hmax_by_type[tt]);
    bvh[idx].left = left_idx;
    bvh[idx].right = right_idx;
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


/* Check if a search sphere overlaps an axis-aligned bounding box (with periodic wrapping).
 * Returns 1 if overlap, 0 otherwise.
 * Uses center+halfwidth approach for correct periodic distance to AABB. */
static int bbox_overlaps_sphere(double box_lo[3], double box_hi[3],
                                double pos[3], double search_r, double search_r2)
{
    double dist2 = 0;
    int k;
    for(k = 0; k < 3; k++)
    {
        int is_periodic = (k == 0) ? TILE_PERIODIC_X : ((k == 1) ? TILE_PERIODIC_Y : TILE_PERIODIC_Z);
        double bsize = (k == 0) ? boxSize_X : ((k == 1) ? boxSize_Y : boxSize_Z);

        /* Distance from pos to center of bbox, with periodic wrapping */
        double center = 0.5 * (box_lo[k] + box_hi[k]);
        double halfwidth = 0.5 * (box_hi[k] - box_lo[k]);
        double dx = fabs(pos[k] - center);
        if(is_periodic && dx > 0.5 * bsize) dx = bsize - dx;

        /* Gap = distance from pos to nearest edge of bbox */
        double gap = dx - halfwidth;
        if(gap <= 0) continue; /* pos is inside bbox on this axis */

        if(gap > search_r) return 0; /* early rejection */
        dist2 += gap * gap;
    }
    return (dist2 < search_r2);
}


/* Process a single tile's particles against particle i. Returns neighbor count. */
static inline int check_tile_particles(struct particle_data *P, int i, double h_i, double h2_i,
                                       sfc_tile_t *tile, int *pool, int search_mode,
                                       int *store_neighbors, int count)
{
    MyDouble xtmp = 0; /* required by NGB_PERIODIC_BOX_LONG_* macros */
    for(int s = 0; s < tile->count; s++)
    {
        int j = pool[tile->first + s];
        double dx_raw = P[i].Pos[0] - P[j].Pos[0];
        double dy_raw = P[i].Pos[1] - P[j].Pos[1];
        double dz_raw = P[i].Pos[2] - P[j].Pos[2];
        double adx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
        double ady = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
        double adz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);

        double pair_search_r2;
        if(search_mode == NGB_SEARCH_ONEWAY) {
            pair_search_r2 = h2_i;
        } else {
            double h_max = DMAX(h_i, P[j].KernelRadius);
            pair_search_r2 = h_max * h_max;
        }

        if(adx > h_i && (search_mode == NGB_SEARCH_ONEWAY || adx * adx > pair_search_r2)) continue;
        double r2 = adx * adx + ady * ady + adz * adz;
        if(r2 < pair_search_r2) {
            if(store_neighbors) store_neighbors[count] = j;
            count++;
        }
    }
    return count;
}


/* Search for neighbors of particle i using BVH traversal over tiles.
 * Iterative stack-based traversal (GPU-friendly, no recursion).
 * Returns count; if store_neighbors != NULL, writes indices there. */
static int search_neighbors_sfc(struct particle_data *P, int i, double h_i,
                                sfc_tile_t *tiles, int ntiles,
                                int *pool, int search_mode,
                                tile_bvh_node_t *bvh, int bvh_root,
                                int *store_neighbors)
{
    double h2_i = h_i * h_i;
    double pos_i[3] = {P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]};
    int count = 0;

    /* Stack-based iterative BVH traversal */
    int stack[TILE_BVH_STACK_SIZE];
    int sp = 0;
    stack[sp++] = bvh_root;

    while(sp > 0)
    {
        int node_idx = stack[--sp];
        tile_bvh_node_t *node = &bvh[node_idx];

        /* Compute search radius for this node: for SYMMETRIC, use node's subtree hmax */
        double search_r = (search_mode == NGB_SEARCH_ONEWAY) ? h_i : DMAX(h_i, node->hmax);
        double search_r2 = search_r * search_r;

        /* Check if node's bbox (expanded by search_r) overlaps particle i */
        if(!bbox_overlaps_sphere(node->lo, node->hi, pos_i, search_r, search_r2)) continue;

        /* Check if leaf */
        if(node->left < 0)
        {
            /* Leaf node: process the tile's particles */
            int tile_idx = -(node->left + 1);
            count = check_tile_particles(P, i, h_i, h2_i, &tiles[tile_idx], pool, search_mode,
                                         store_neighbors, count);
        }
        else
        {
            /* Internal node: push children onto stack */
            if(sp + 2 > TILE_BVH_STACK_SIZE) {
                /* Stack overflow — fall back to processing both children inline.
                   This shouldn't happen with STACK_SIZE=64 (supports 2^32 tiles). */
                endrun(77701);
            }
            stack[sp++] = node->left;
            stack[sp++] = node->right;
        }
    }
    return count;
}


void build_neighbor_list_sfc(struct particle_data *P, struct gas_cell_data *CellP,
                             int num_total, int *active_indices, int num_active,
                             int search_mode, int type_bitmask,
                             neighbor_list_t *out)
{
    out->num_active = num_active;
    double t_sfc_start = my_second();

    /* Build SFC tiles */
    sfc_tile_t *tiles;
    int *pool;
    int num_pool;
    int ntiles = build_sfc_tiles(P, num_total, type_bitmask, TILE_TARGET_SIZE,
                                 &tiles, &pool, &num_pool);
    double t_sfc_tiles = timediff(t_sfc_start, my_second());

    if(ThisTask == 0) {
        double max_h = 0;
        int t;
        for(t = 0; t < ntiles; t++) if(tiles[t].hmax > max_h) max_h = tiles[t].hmax;
        PRINT_STATUS("SFC-tile neighbor finder: %d tiles (%d particles/tile avg), max_h=%.4g, pool=%d particles",
                     ntiles, num_pool / (ntiles > 0 ? ntiles : 1), max_h, num_pool);
    }

    if(num_pool == 0 || ntiles == 0) {
        out->total_pairs = 0;
        myfree(tiles); myfree(pool);
        out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
        memset(out->offsets, 0, (size_t)(num_active + 1) * sizeof(int64_t));
        out->neighbors = NULL;
        return;
    }

    /* Build BVH over tiles for O(log ntiles) spatial pruning */
    tile_bvh_node_t *bvh;
    int bvh_nnodes = build_tile_bvh(tiles, ntiles, &bvh);
    int bvh_root = bvh_nnodes - 1;
    double t_sfc_bvh_end = my_second();

    /* Pass 1: count neighbors per active particle */
    int *counts = (int *) mymalloc("sfc_counts", num_active * sizeof(int));
    int a;
    for(a = 0; a < num_active; a++) {
        int i = active_indices[a];
        counts[a] = search_neighbors_sfc(P, i, P[i].KernelRadius,
                        tiles, ntiles, pool, search_mode, bvh, bvh_root, NULL);
    }

    int64_t total_pairs = 0;
    for(a = 0; a < num_active; a++) total_pairs += counts[a];
    double t_sfc_pass1_end = my_second();

    /* Free pass-1 temporaries in reverse mymalloc order before allocating output.
       Stack is: pool, tiles, bvh, counts (top). Output arrays must be at stack base. */
    myfree(counts);
    myfree(bvh);
    myfree(tiles);
    myfree(pool);

    /* Allocate output CSR (these persist after return) */
    out->total_pairs = total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(total_pairs > 0 ? total_pairs : 1) * sizeof(int));

    /* Pass 2: rebuild tiles + BVH and fill neighbor indices */
    ntiles = build_sfc_tiles(P, num_total, type_bitmask, TILE_TARGET_SIZE,
                             &tiles, &pool, &num_pool);
    bvh_nnodes = build_tile_bvh(tiles, ntiles, &bvh);
    bvh_root = bvh_nnodes - 1;

    out->offsets[0] = 0;
    for(a = 0; a < num_active; a++) {
        int i = active_indices[a];
        int n = search_neighbors_sfc(P, i, P[i].KernelRadius,
                        tiles, ntiles, pool, search_mode, bvh, bvh_root,
                        &out->neighbors[out->offsets[a]]);
        out->offsets[a + 1] = out->offsets[a] + n;
    }

    /* Cleanup pass-2 temporaries (reverse order) */
    myfree(bvh);
    myfree(tiles);
    myfree(pool);

    /* timing variables available for debugging if needed, but not printed in production */
}
