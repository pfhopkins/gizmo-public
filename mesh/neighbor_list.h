/* neighbor_list.h — cell-list neighbor finder producing CSR neighbor lists.
 *
 * Builds a compressed sparse row (CSR) neighbor list for a set of active
 * particles, using a uniform grid cell-list for spatial indexing.
 *
 * For particle a (0-indexed into the active list), its neighbors are:
 *   neighbors[offsets[a]] .. neighbors[offsets[a+1]-1]
 * Each neighbor index refers to the global P[]/CellP[] arrays.
 *
 * Handles: periodic boundaries, shearing-box Y-offsets, one-way vs symmetric
 * search, particle type filtering via bitmask.
 *
 * NOTE: This is a placeholder spatial index. For production use (especially
 * zoom-in simulations with non-uniform h), this should be upgraded to an
 * SFC-tile finder (Cornerstone/SPH-EXA pattern) which handles varying
 * resolution naturally. The CSR output interface is identical — drop-in
 * replacement.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H

#include <stdint.h>

/* CSR-format neighbor list.
 * offsets / total_pairs are 64-bit: a large symmetric search (e.g. m11i +
 * TURB_DIFF_DYNAMIC) can exceed 2^31 neighbor pairs. neighbors stays int *:
 * its *values* are particle indices into P[]/CellP[], bounded by num_total
 * (< 2^31); only the array *length* is 64-bit. */
struct neighbor_list_t {
    int64_t *offsets;       /* [num_active+1] CSR row pointers (cumulative pair counts) */
    int *neighbors;         /* [total_pairs] neighbor particle indices into P[]/CellP[] */
    int num_active;         /* number of active particles (rows in CSR) */
    int64_t total_pairs;    /* total number of neighbor pairs stored */
};

/* Search modes */
#define NGB_SEARCH_ONEWAY   0  /* r_ij < h_i only (density) */
#define NGB_SEARCH_SYMMETRIC 1 /* r_ij < max(h_i, h_j) (gradients, hydro) */

/* Build a CSR neighbor list for the given active particles.
 *
 * P, CellP:         particle arrays (including ghosts at indices >= NumPart)
 * num_total:         total number of particles including ghosts
 * active_indices:    indices of particles to find neighbors for
 * num_active:        number of active particles
 * search_mode:       NGB_SEARCH_ONEWAY or NGB_SEARCH_SYMMETRIC
 * type_bitmask:      which P[j].Type to include as neighbors (e.g., 1 for gas-only)
 * out:               output neighbor list (caller must call free_neighbor_list when done)
 *
 * Memory: uses mymalloc for out->offsets and out->neighbors. Caller must call
 * free_neighbor_list() before any myfree of allocations made before this call.
 */
void build_neighbor_list(struct particle_data *P, struct gas_cell_data *CellP,
                         int num_total, int *active_indices, int num_active,
                         int search_mode, int type_bitmask,
                         neighbor_list_t *out);

void free_neighbor_list(neighbor_list_t *list);

/* Global symmetric neighbor list cached between density and hydro_force.
 * Built after density converges, freed after hydro_force completes.
 * Used by gradients and hydro force. */
extern neighbor_list_t gizmo_sym_neighbor_list;
extern int *gizmo_sym_active_indices;
extern int gizmo_sym_num_active;
extern int gizmo_sym_num_active_global;
void gizmo_sym_neighbor_list_free(void);

#endif /* NEIGHBOR_LIST_H */
