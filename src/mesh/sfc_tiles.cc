/* sfc_tiles.cc — the ghost supply pool: the index members of P[], in P[] order.
 *
 * Membership is sfc_pool_member (sfc_tiles.h), shared with the spatial index build in
 * gpu_neighbor_list.cc; the tiles, rows and BVH are built there.
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
