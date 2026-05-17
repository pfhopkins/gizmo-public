/* neighbor_list.cc — cell-list neighbor finder for GPU/CPU neighbor iteration.
 *
 * Builds a uniform grid cell-list over all pool particles, then for each
 * active particle searches neighboring cells to construct a CSR neighbor list.
 *
 * Periodic boundaries: uses the same NGB_PERIODIC_BOX_LONG_X/Y/Z macros as
 * the tree-walk neighbor finder, ensuring identical distance calculations.
 *
 * Shearing box: expands the Y cell stencil near the X boundary to account for
 * the shear-dependent position offset (Shearing_Box_Pos_Offset).
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
#include "neighbor_list.h"


/* Axis periodicity flags (matching macros.h logic) */
#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_X) && !defined(BOX_OUTFLOW_X)
#define CELLLIST_PERIODIC_X 1
#else
#define CELLLIST_PERIODIC_X 0
#endif

#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_Y) && !defined(BOX_OUTFLOW_Y)
#define CELLLIST_PERIODIC_Y 1
#else
#define CELLLIST_PERIODIC_Y 0
#endif

#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_Z) && !defined(BOX_OUTFLOW_Z)
#define CELLLIST_PERIODIC_Z 1
#else
#define CELLLIST_PERIODIC_Z 0
#endif


/* Per-active-particle neighbor search over the cell grid.
 * If store_neighbors is non-NULL, writes neighbor indices there; otherwise just counts.
 * Returns the number of neighbors found for this particle. */
static int search_neighbors_for_particle(
    struct particle_data *P, int i, double h_i,
    int *sorted, int *cell_start, int ncells[3], double cell_size[3], double box_lo[3],
    int search_mode, int shear_offset_cells,
    int *store_neighbors)
{
    double h2_i = h_i * h_i;
    MyDouble xtmp = 0; /* required by NGB_PERIODIC_BOX_LONG_* macros */

    int ci_x = (int)((P[i].Pos[0] - box_lo[0]) / cell_size[0]);
    int ci_y = (int)((P[i].Pos[1] - box_lo[1]) / cell_size[1]);
    int ci_z = (int)((P[i].Pos[2] - box_lo[2]) / cell_size[2]);
    if(ci_x < 0) ci_x = 0; if(ci_x >= ncells[0]) ci_x = ncells[0]-1;
    if(ci_y < 0) ci_y = 0; if(ci_y >= ncells[1]) ci_y = ncells[1]-1;
    if(ci_z < 0) ci_z = 0; if(ci_z >= ncells[2]) ci_z = ncells[2]-1;

    int count = 0;
    int ny = ncells[1], nz = ncells[2];

    for(int dcx = -1; dcx <= 1; dcx++)
    {
        int ncx = ci_x + dcx;
        int x_wrapped = 0;
#if CELLLIST_PERIODIC_X
        if(ncx < 0)            { ncx += ncells[0]; x_wrapped = 1; }
        else if(ncx >= ncells[0]) { ncx -= ncells[0]; x_wrapped = 1; }
#else
        if(ncx < 0 || ncx >= ncells[0]) continue;
#endif

        /* Y stencil: normally +-1, expanded near X boundary for shearing box */
        int dy_lo = -1, dy_hi = 1;
#ifdef BOX_SHEARING
        if(x_wrapped && shear_offset_cells > 0) {
            dy_lo = -1 - shear_offset_cells;
            dy_hi =  1 + shear_offset_cells;
        }
#endif

        for(int dcy = dy_lo; dcy <= dy_hi; dcy++)
        {
            int ncy = ci_y + dcy;
#if CELLLIST_PERIODIC_Y
            ncy = ((ncy % ny) + ny) % ny;
#else
            if(ncy < 0 || ncy >= ny) continue;
#endif

            for(int dcz = -1; dcz <= 1; dcz++)
            {
                int ncz = ci_z + dcz;
#if CELLLIST_PERIODIC_Z
                ncz = ((ncz % nz) + nz) % nz;
#else
                if(ncz < 0 || ncz >= nz) continue;
#endif

                int nc = ncx * ny * nz + ncy * nz + ncz;
                for(int s = cell_start[nc]; s < cell_start[nc + 1]; s++)
                {
                    int j = sorted[s];

                    /* Distance check using the same macros as the tree walk.
                       NGB_PERIODIC_BOX_LONG_* return absolute (unsigned) distance
                       to nearest periodic image. Arguments are raw signed differences.
                       The Y macro uses the X-difference for shearing box offset. */
                    double dx_raw = P[i].Pos[0] - P[j].Pos[0];
                    double dy_raw = P[i].Pos[1] - P[j].Pos[1];
                    double dz_raw = P[i].Pos[2] - P[j].Pos[2];
                    double adx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
                    double ady = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
                    double adz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);

                    double search_r2;
                    if(search_mode == NGB_SEARCH_ONEWAY) {
                        search_r2 = h2_i;
                    } else {
                        double h_max = DMAX(h_i, P[j].KernelRadius);
                        search_r2 = h_max * h_max;
                    }

                    /* Early axis-aligned rejection (same pattern as ngb_codeblock) */
                    if(adx > h_i && (search_mode == NGB_SEARCH_ONEWAY || adx * adx > search_r2)) continue;
                    double r2 = adx * adx + ady * ady + adz * adz;
                    if(r2 < search_r2) {
                        if(store_neighbors) store_neighbors[count] = j;
                        count++;
                    }
                }
            }
        }
    }
    return count;
}


void build_neighbor_list(struct particle_data *P, struct gas_cell_data *CellP,
                         int num_total, int *active_indices, int num_active,
                         int search_mode, int type_bitmask,
                         neighbor_list_t *out)
{
    out->num_active = num_active;

    /* Step 1: Find max kernel radius among pool particles */
    double max_h = 0;
    int num_pool = 0;
    for(int i = 0; i < num_total; i++) {
        if(!((1 << P[i].Type) & type_bitmask)) continue;
        if(P[i].Mass <= 0) continue;
        if(P[i].KernelRadius > max_h) max_h = P[i].KernelRadius;
        num_pool++;
    }

    if(max_h <= 0 || num_pool == 0) {
        out->total_pairs = 0;
        out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
        memset(out->offsets, 0, (size_t)(num_active + 1) * sizeof(int64_t));
        out->neighbors = NULL;
        return;
    }

    /* Step 2: Determine cell grid covering the domain */
    double box_lo[3], box_hi[3], cell_size[3];
    int ncells[3];

    /* For periodic axes, use the full box. For non-periodic, compute bounding box. */
    for(int d = 0; d < 3; d++) {
        int is_periodic = (d == 0) ? CELLLIST_PERIODIC_X : ((d == 1) ? CELLLIST_PERIODIC_Y : CELLLIST_PERIODIC_Z);
        double bsize = (d == 0) ? boxSize_X : ((d == 1) ? boxSize_Y : boxSize_Z);
        if(is_periodic) {
            box_lo[d] = 0;
            box_hi[d] = bsize;
        } else {
            box_lo[d] =  MAX_REAL_NUMBER;
            box_hi[d] = -MAX_REAL_NUMBER;
            for(int i = 0; i < num_total; i++) {
                if(!((1 << P[i].Type) & type_bitmask) || P[i].Mass <= 0) continue;
                if(P[i].Pos[d] < box_lo[d]) box_lo[d] = P[i].Pos[d];
                if(P[i].Pos[d] > box_hi[d]) box_hi[d] = P[i].Pos[d];
            }
            box_lo[d] -= max_h;
            box_hi[d] += max_h;
        }

        double extent = box_hi[d] - box_lo[d];
        ncells[d] = (int)(extent / max_h);
        if(ncells[d] < 3) ncells[d] = 3;       /* need at least 3 for +-1 stencil to make sense */
        if(ncells[d] > 512) ncells[d] = 512;    /* cap memory: 512^3 = 134M cells max */
        cell_size[d] = extent / ncells[d];
    }

    int ncells_total = ncells[0] * ncells[1] * ncells[2];
    if(ThisTask == 0) {
        PRINT_STATUS("Cell-list neighbor finder: grid %d x %d x %d = %d cells, cell_size=(%.4g, %.4g, %.4g), max_h=%.4g, pool=%d particles",
                     ncells[0], ncells[1], ncells[2], ncells_total,
                     cell_size[0], cell_size[1], cell_size[2], max_h, num_pool);
    }

    /* Step 3: Assign pool particles to cells (counting sort) */
    int *pcell = (int *) mymalloc("pcell", num_total * sizeof(int));
    int *cell_count = (int *) mymalloc("cell_count", ncells_total * sizeof(int));
    memset(cell_count, 0, ncells_total * sizeof(int));

    int ny = ncells[1], nz = ncells[2];
    for(int i = 0; i < num_total; i++) {
        if(!((1 << P[i].Type) & type_bitmask) || P[i].Mass <= 0) {
            pcell[i] = -1;
            continue;
        }
        int cx = (int)((P[i].Pos[0] - box_lo[0]) / cell_size[0]);
        int cy = (int)((P[i].Pos[1] - box_lo[1]) / cell_size[1]);
        int cz = (int)((P[i].Pos[2] - box_lo[2]) / cell_size[2]);
        if(cx < 0) cx = 0; if(cx >= ncells[0]) cx = ncells[0]-1;
        if(cy < 0) cy = 0; if(cy >= ncells[1]) cy = ncells[1]-1;
        if(cz < 0) cz = 0; if(cz >= ncells[2]) cz = ncells[2]-1;
        pcell[i] = cx * ny * nz + cy * nz + cz;
        cell_count[pcell[i]]++;
    }

    /* Step 4: Prefix sum for cell_start */
    int *cell_start = (int *) mymalloc("cell_start", (ncells_total + 1) * sizeof(int));
    cell_start[0] = 0;
    for(int c = 0; c < ncells_total; c++)
        cell_start[c + 1] = cell_start[c] + cell_count[c];

    /* Step 5: Scatter particle indices into sorted order */
    int *sorted = (int *) mymalloc("sorted", (num_pool > 0 ? num_pool : 1) * sizeof(int));
    memset(cell_count, 0, ncells_total * sizeof(int)); /* reuse as fill counter */
    for(int i = 0; i < num_total; i++) {
        if(pcell[i] < 0) continue;
        int c = pcell[i];
        sorted[cell_start[c] + cell_count[c]] = i;
        cell_count[c]++;
    }

    /* Step 6: Shearing box Y-offset in cells */
    int shear_offset_cells = 0;
#ifdef BOX_SHEARING
    shear_offset_cells = (int)(fabs(Shearing_Box_Pos_Offset) / cell_size[1]) + 2;
#endif

    /* Step 7: Pass 1 — count neighbors per active particle */
    int *counts = (int *) mymalloc("ngb_counts", num_active * sizeof(int));
    for(int a = 0; a < num_active; a++) {
        int i = active_indices[a];
        counts[a] = search_neighbors_for_particle(P, i, P[i].KernelRadius,
                        sorted, cell_start, ncells, cell_size, box_lo,
                        search_mode, shear_offset_cells, NULL);
    }

    /* Sum counts to get total pairs, then free all temporaries before allocating output */
    int64_t total_pairs = 0;
    for(int a = 0; a < num_active; a++) total_pairs += counts[a];

    /* Free all temporaries in reverse mymalloc order.
       Stack is: pcell, cell_count, cell_start, sorted, counts (top) */
    myfree(counts);
    myfree(sorted);
    myfree(cell_start);
    myfree(cell_count);
    myfree(pcell);

    /* Allocate output CSR arrays (offsets first, then neighbors — persist after return) */
    out->total_pairs = total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(total_pairs > 0 ? total_pairs : 1) * sizeof(int));

    /* Step 8: Pass 2 — rebuild cell grid and fill neighbor indices.
       We rebuild because the pass-1 temporaries were freed above to satisfy
       mymalloc stack ordering (output arrays must be at base of stack). */
    pcell = (int *) mymalloc("pcell2", num_total * sizeof(int));
    cell_count = (int *) mymalloc("cell_count2", ncells_total * sizeof(int));
    memset(cell_count, 0, ncells_total * sizeof(int));
    for(int i = 0; i < num_total; i++) {
        if(!((1 << P[i].Type) & type_bitmask) || P[i].Mass <= 0) { pcell[i] = -1; continue; }
        int cx = (int)((P[i].Pos[0] - box_lo[0]) / cell_size[0]);
        int cy = (int)((P[i].Pos[1] - box_lo[1]) / cell_size[1]);
        int cz = (int)((P[i].Pos[2] - box_lo[2]) / cell_size[2]);
        if(cx < 0) cx = 0; if(cx >= ncells[0]) cx = ncells[0]-1;
        if(cy < 0) cy = 0; if(cy >= ncells[1]) cy = ncells[1]-1;
        if(cz < 0) cz = 0; if(cz >= ncells[2]) cz = ncells[2]-1;
        pcell[i] = cx * ny * nz + cy * nz + cz;
        cell_count[pcell[i]]++;
    }
    cell_start = (int *) mymalloc("cell_start2", (ncells_total + 1) * sizeof(int));
    cell_start[0] = 0;
    for(int c = 0; c < ncells_total; c++) cell_start[c + 1] = cell_start[c] + cell_count[c];
    sorted = (int *) mymalloc("sorted2", (num_pool > 0 ? num_pool : 1) * sizeof(int));
    memset(cell_count, 0, ncells_total * sizeof(int));
    for(int i = 0; i < num_total; i++) {
        if(pcell[i] < 0) continue;
        int c = pcell[i];
        sorted[cell_start[c] + cell_count[c]] = i;
        cell_count[c]++;
    }

    /* Fill CSR: for each active particle, search and store neighbors */
    out->offsets[0] = 0;
    for(int a = 0; a < num_active; a++) {
        int i = active_indices[a];
        int n = search_neighbors_for_particle(P, i, P[i].KernelRadius,
                        sorted, cell_start, ncells, cell_size, box_lo,
                        search_mode, shear_offset_cells,
                        &out->neighbors[out->offsets[a]]);
        out->offsets[a + 1] = out->offsets[a] + n;
    }

    /* Cleanup pass-2 temporaries (reverse order) */
    myfree(sorted);
    myfree(cell_start);
    myfree(cell_count);
    myfree(pcell);
}


void free_neighbor_list(neighbor_list_t *list)
{
    /* Free in reverse mymalloc order: neighbors was allocated after offsets */
    if(list->neighbors) myfree(list->neighbors);
    if(list->offsets) myfree(list->offsets);
    list->neighbors = NULL;
    list->offsets = NULL;
    list->total_pairs = 0;
}

/* Global symmetric neighbor list: defined here, extern'd in neighbor_list.h */
neighbor_list_t gizmo_sym_neighbor_list = {NULL, NULL, 0, 0};
int *gizmo_sym_active_indices = NULL;
int gizmo_sym_num_active = 0;
int gizmo_sym_num_active_global = 0; /* Allreduce'd in gizmo_gradients_prep_symlist; reused by refresh */

void gizmo_sym_neighbor_list_free(void) {
    if(gizmo_sym_neighbor_list.offsets) { free_neighbor_list(&gizmo_sym_neighbor_list); }
    if(gizmo_sym_active_indices) { myfree(gizmo_sym_active_indices); gizmo_sym_active_indices = NULL; }
    gizmo_sym_num_active = 0;
}
