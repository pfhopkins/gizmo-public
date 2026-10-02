/* gpu_topology_finalize.h
 *
 * Final stage of GPU tree-build: replaces the sibling / father / Father[]
 * outputs of force_update_node_recursive with three GPU passes.
 *
 *   gpu_topology_finalize_father()    Writes soa->father[k] for every
 *                                     internal node and Father[i] for every
 *                                     particle child in one kernel that
 *                                     scans soa->suns_backup.
 *
 *   gpu_topology_finalize_sibling()   Writes soa->sibling[k] for every
 *                                     internal node by walking up the father
 *                                     chain to find the next occupied
 *                                     sibling-or-aunt slot.  Pre-condition:
 *                                     finalize_father has run.
 *
 *   gpu_topology_writeback_d_to_aos() Pushes soa->sibling and soa->father
 *                                     back into Nodes_base[k].u.d.{sibling,
 *                                     father} for legacy CPU walks.  Clobbers
 *                                     u.suns via the union — safe because
 *                                     soa->suns_backup is the SoA truth.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GIZMO_GPU_TOPOLOGY_FINALIZE_H
#define GIZMO_GPU_TOPOLOGY_FINALIZE_H


#ifdef __cplusplus
extern "C" {
#endif

/* Father pass.  Writes:
 *   - soa->father[k] for k in [0, n) (n = Numnodestree).  Index k=0 (root)
 *     gets -1.
 *   - Father[i] for i in [0, TreeParticleSlots) (UVM-mirrored host array).
 * Pre-condition: soa->suns_backup populated for [0, n).  Returns 0 on
 * success. */
int gpu_topology_finalize_father(int n);

/* Sibling pass.  Writes soa->sibling[k] for k in [0, n).
 * Pre-condition: finalize_father has run (reads soa->father).  Returns 0
 * on success. */
int gpu_topology_finalize_sibling(int n);

/* AoS writeback for d.sibling/d.father.  Host loop over
 * [0, n).  Idempotent.  Returns 0 on success. */
int gpu_topology_writeback_d_to_aos(int n);

/* GPU kernel that resets per-node ephemeral fields (GravCost,
 * Ti_current, dp, Ti_lastkicked, Flag, optional payloads) after topology
 * finalize and before moment_refresh accumulates fresh moments.  Operates
 * on the absolute Nodes[] index range [TreeNodeIndexBase .. TreeNodeIndexBase+n).  On the CPU
 * path FUNR does this work inline; on the GPU path FUNR is retired
 * so this kernel takes its place.  Returns 0 on success. */
int gpu_node_reset_ephemeral(int n);

/* Father[] is UVM (SharedSpace) on GPU build so the GPU father
 * kernel can write into it directly and host readers page-fault on touch.
 * These wrappers keep the kokkos_malloc/free out of forcetree.cc (non-GPU TU
 * has no Kokkos include). */
int  *gpu_father_alloc(int tree_particle_slots);   /* returns NULL on failure */
void  gpu_father_free(int *p);

/* Extend the same UVM-allocator pattern to Nodes_base, Extnodes_base
 * and Nextnode: all four tree-storage arrays live in SharedSpace,
 * so GPU kernels can read/write them directly and the host writeback
 * loops over node-count are replaced with single-pass GPU kernels.
 *
 * NODE / extNODE struct definitions live in forcetree.h (which this TU
 * includes); the helpers take byte counts so callers stay struct-agnostic. */
void *gpu_tree_alloc_bytes(size_t bytes, const char *label);  /* generic SharedSpace alloc; label names the buffer */
void  gpu_tree_free_bytes(void *p);        /* matching free */

/*! Put a shared-space block into the state its release is cheapest from; call it
 *  immediately before releasing one.  gpu_tree_free_bytes and the gravity-tree mirror
 *  already do; gpu_father_free does NOT, because the measurement that justifies this
 *  covered the node arrays and the mirror and did not cover Father[].  Blocks whose
 *  length was never recorded are left alone. */
void  gizmo_gpu_prepare_shared_for_free(void *ptr);
void  gizmo_gpu_shared_track(void *ptr, size_t bytes);
/*! Tell the driver which side a shared block belongs to: 0 none / 1 host / 2 device, the
 *  three values the particle-storage knob already uses.  An unhinted managed block settles
 *  where it is FIRST TOUCHED, so one the host fills and the device reads stays host-resident
 *  and is fetched a page at a time for the whole of every device read of it. */
void  gizmo_gpu_place_shared(void *p, size_t nbytes, int placement, const char *what);

#ifdef __cplusplus
}
#endif


#endif /* GIZMO_GPU_TOPOLOGY_FINALIZE_H */
