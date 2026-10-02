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
#include "../declarations/constants.h"   /* NODELISTLENGTH */
#include "../core/timestep_functions.h" /* DriftKickTableView, held by value below */

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

/* Fixed-size export envelope carried from a requesting rank to a supply rank:
 * one query plus the start nodes the sender's walk reached on that peer.  A
 * (query,peer) needing more than NODELISTLENGTH nodes SPLITS into multiple
 * envelopes; the receiver walks each independently and dedups through the
 * matched bitmap, so a split costs an extra envelope and nothing else.
 *
 * Lives here rather than beside the host exchange code because both the host
 * receiver walk and the device receiver traversal consume it, and the device
 * one compiles in a different translation unit. */
struct gx_export_envelope_t {
    double pos[3];
    double h;
    int    n_nodes;
    int    nodes[NODELISTLENGTH];
    int    _pad;
};

/* Lives here, in a header that pulls in no Kokkos, rather than beside the device
 * traversal that consumes it.  It is plain data -- pointers and index bounds --
 * and the runner declares one as a member, so the header carrying it is included
 * by host translation units compiled without a device compiler.  The traversal
 * header cannot serve that role: it includes <Kokkos_Core.hpp> deliberately and
 * unconditionally, because the walk calls Kokkos itself, and a host unit that
 * picks that up fails with the CUDA setup header's "__CUDACC__ not defined"
 * error rather than anything that names the real cause.
 *
 */
/* The tree as the device sees it: the mirrored node arrays plus the boundaries
 * that separate the three index classes a walk can encounter.
 *
 * An index below `local_particle_slots` is a particle this rank owns; one from
 * there up to `particle_slots` is an imported ghost, which the walk reaches but
 * never reports.  One at or above `node_base` and below `pseudo_start` is a
 * node, of which those at or above `foreign_base` are imported subtrees holding
 * no local particles.  One at or above `pseudo_start` is a pseudo-particle
 * standing for another rank's subtree.  Anything in the gap between the
 * particle slots and the node base belongs to no class at all and means the
 * tree is malformed. */
struct GxDeviceTreeView {
    const Vec3<MyFloat> *node_center;
    const MyFloat       *node_len;
    const int           *node_sibling;
    const int           *node_nextnode;
    const unsigned int  *node_bitflags;
    const int           *nextnode_aux;
    int                  node_base;
    int                  particle_slots;
    int                  local_particle_slots = -1;   /* owned locals; ghosts sit above this */
    int                  node_capacity;
    int                  foreign_base;
    int                  pseudo_start;

    /* WIDEN-ON-OPEN (landing 4).  A device walk cannot take a lock, so it cannot
     * drift a node it reaches.  Instead it widens the node's own opening bound by
     * how far that node could have moved since the mirror was written:
     *     len_effective = len + TREE_NODE_WIDENING_DELTA(vmax, dt(node_ti -> now))
     * which is the SAME expression the sweep and force_drift_node apply -- the walk
     * just evaluates it lazily, for the ~4k nodes it visits, instead of eagerly for
     * ~1.4M.  Over-widening is harmless (over-inclusion, re-gated by the pair
     * kernel); UNDER-widening is silent under-inclusion, which is why `node_vmax`
     * must never be staler than `node_len` (a running max only grows) and
     * `node_ti` must never be FRESHER than `node_len`.
     *
     * Null disables widening and the walk opens on the stored length alone -- the
     * pre-landing-4 behaviour, which is correct only when something else has
     * certified the geometry current. */
    const MyGravFloat   *node_vmax = nullptr;
    const integertime   *node_ti   = nullptr;
    integertime          ti_now    = 0;
    /* BY VALUE, not by pointer. The view is captured into the device walk, and the
     * struct itself lives in host .bss -- only the table ARRAYS it names are in
     * shared space. A pointer to it would be dereferenced on device and is illegal
     * on CUDA/HIP. `drift_tables_ok` says whether it was filled. */
    struct DriftKickTableView drift_tables{};
    int                       drift_tables_ok = 0;
};


/* The distinct set of local particles a fused walk reaches, so that only those
 * have to be brought current rather than the whole rank.
 *
 * A fused walk evaluates the pair kernel at the leaf it reaches and cannot drift
 * a stale particle when it gets there -- that needs a lock.  The way out is to
 * reach the leaves twice: once to record which ones the walk touches, then a
 * drift of just those, then the evaluation.  Both passes traverse the same tree
 * with the same queries, and drifting a particle changes nothing the traversal
 * reads, so the second pass reaches exactly the set the first recorded.
 *
 * `seen` is a generation stamp rather than a set of flags, which is what keeps
 * the cost proportional to what the walk touches: it is allocated once, and the
 * generation is bumped instead of the array being cleared.  Clearing it per call
 * is the one thing that would put this back at O(local particles).
 *
 * The generation advances once per CALL, not once per pass, and the distinction
 * is load-bearing: a particle advanced for an earlier pass is still current for
 * the later ones, so carrying its claim across them is what stops the self walk
 * and every peer round re-examining the same particle.  Only the append cursor
 * is reset between passes.
 *
 * The stamp is indexed by the same local-particle slot the traversal reports, so
 * a slot can be claimed at most once per generation and the compacted list can
 * never be longer than the number of slots.  That is what removes the overflow
 * case rather than handling it.
 *
 * Lives in this header, not beside the traversal, for the reason the tree view
 * does: plain data, read by host units that no device compiler ever sees. */
/* Which phase may claim right now.  One generation deliberately spans several passes of one
 * fused call while the cursor resets at each consume, so "whose generation is this" cannot be
 * inferred from the cursor and has to be named.  NONE is the retired state: every claim in it is
 * out of phase, which is what makes a claim from a kernel that outlived its epoch visible. */
enum gx_touched_owner_t {
    GX_TOUCHED_OWNER_NONE       = 0,
    GX_TOUCHED_OWNER_FUSED_WALK = 1    /* the Mode-D fused walk's discovery passes */
};

struct GxTouchedSet {
    unsigned int *seen     = nullptr;  /* [capacity] generation stamps, never cleared */
    int          *list     = nullptr;  /* [capacity] compacted distinct local indices */
    int          *counter  = nullptr;  /* [1] append cursor for the current pass */
    int           capacity = 0;        /* owned local particle slots at the last ensure */
    unsigned int  gen      = 0;        /* this call's generation */
    int           owner    = GX_TOUCHED_OWNER_NONE;  /* constant within an epoch, like gen */
};


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
