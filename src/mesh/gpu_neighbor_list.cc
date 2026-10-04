/* gpu_neighbor_list.cc — GPU-accelerated neighbor list construction.
 *
 * Extracted from hydro/density_gpu.cc for reuse by any code that needs
 * a GPU-built CSR neighbor list (density, symmetric list, future loops).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <Kokkos_Core.hpp>

/* GPU All mirror: per-TU managed pointer to shared UVM allocation. */
#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../system/gpu_particles_arena.h"
#include "../core/timestep_functions.h"   /* DriftKickTableView for the walk widening */


/* The drift-factor interpolator the walk's widening uses.
 *
 * Refreshed once per call rather than per node: it is the same table the node
 * sweep mirrors (gpu_force_drift.cc), and a walk that built its own would be a
 * second interpolator answering the same question. Null on failure, which simply
 * disables widening -- the walk then opens on the stored length, and the caller's
 * certification is what makes that legal. */

static struct DriftKickTableView g_walk_drift_tables;
static double *g_walk_drift_storage = NULL;
static int g_walk_drift_tables_ok = 0;

/* Refresh once per preparation. Shares drift_kick_table_mirror_refresh with the
   particle drift and the node sweep -- same table, same units, one owner. */
static void gx_walk_drift_tables_refresh(void)
{
    g_walk_drift_tables_ok =
        (drift_kick_table_mirror_refresh(&g_walk_drift_storage, &g_walk_drift_tables) == 0);
}
#include "../core/proto.h"
#include "../declarations/gpu_error_check.h"

#include "sfc_tiles.h"
#include "sfc_tiles_functions.h"
#include "gpu_neighbor_list.h"
#include "gpu_dirty_tracker.h"
#include "neighbor_list.h"
#include "ghost_writeback.h"  /* ghost_write_detector_resnapshot_after_lazy_drift, ghost_get_num_local */
#include "ghost_exchange_functions.h"   /* the shared accept predicates (canonical box wrap) */
#include "device_tree_walk.h"          /* the one device tree traversal */
#include "../gravity/forcetree.h"       /* BITFLAG_TOPLEVEL, force_host_lazy_drift_ti */
#include "../gravity/gpu_gravity_tree.h" /* node SoA + drift certification */
#include "../gravity/gpu_morton_functions.h" /* Morton128 keys for the index build */
#include "../gravity/ags_functions.h"        /* the AGS inputs of the after-drift reach rule */
#include "../core/predict_functions.h"       /* box_wrap_position_to_primary_image */
#include <Kokkos_Sort.hpp>
#include <type_traits>
#include <cstdarg>
#include <stdexcept>

/* TILE_PERIODIC_X/Y/Z defined in sfc_tiles.h (included via gpu_neighbor_list.h) */

/* The kept indexes (gpu_neighbor_list.h): the gas index, whose owned segment is kept across sync points,
 * and the all-types index, released at every sync point.  Both are built lazily by the first list build
 * that uses them. */
static gpu_spatial_index_t g_step_sidx{};
static gpu_spatial_index_t g_step_sidx_alltypes{};


gpu_spatial_index_t *gpu_step_sidx_ptr(void) { return &g_step_sidx; }
gpu_spatial_index_t *gpu_step_sidx_alltypes_ptr(void) { return &g_step_sidx_alltypes; }


/* A radius was written (gpu_neighbor_list.h).  Each kept owned segment registers its source range
 * with the dirty tracker, so a mark reaches every segment that holds the particle; the ghost-exchange
 * supply cache holds membership only and needs none. */
void gizmo_mark_kernel_radius_dirty_indices(const int *indices, int n)
{
    if(!indices || n <= 0) return;
    gpu_dirty_tracker_mark_indices(indices, n);
}
void gizmo_mark_kernel_radius_dirty_range(int start, int end)
{
    if(end <= start) return;
    gpu_dirty_tracker_mark_range(start, end);
}

/* The owned particles' epoch: bumped whenever an index's rows over the owned
 * particles stop describing them in a way the particle count cannot show (see
 * gpu_sidx_notify_owned_changed). gpu_ngb_list_build stamps it into each index at
 * build and requires it to still match before reusing it. The imported particles
 * need no counter here: the ghost exchange already records whether a pool is live
 * and which import produced it (ghost_pool_is_live, ghost_provenance_epoch). */
static uint64_t g_sidx_owned_epoch = 0;

void gpu_sidx_notify_owned_changed(void)
{
    g_sidx_owned_epoch++;
}

/* Particles that stopped being gas in place, counted; see gpu_sidx_notify_member_lost. */
static uint64_t g_sidx_member_loss_epoch = 0;

void gpu_sidx_notify_member_lost(int particle)
{
    g_sidx_member_loss_epoch++;
    gpu_dirty_tracker_mark_indices(&particle, 1);   /* the gas index retires its slot at its next maintenance */
}

/* The order peano_hilbert_order last left the rank's particles in, recorded with the layout change that
 * followed it (gpu_sidx_notify_owned_changed): it holds while the owned epoch, the time and the particle
 * count are the ones recorded, and the next layout change of any kind advances the epoch. */
static struct {int valid; uint64_t owned_epoch; integertime ti; int num_part;} g_sidx_particles_ordered = {0, 0, 0, 0};

void gpu_sidx_notify_owned_reordered(void)
{
    gpu_sidx_notify_owned_changed();
    g_sidx_particles_ordered.valid = 1;
    g_sidx_particles_ordered.owned_epoch = g_sidx_owned_epoch;
    g_sidx_particles_ordered.ti = gizmo_host_ti_current();
    g_sidx_particles_ordered.num_part = NumPart;
}

/* Whether the particles [0, owned_end) are still in the order that decomposition left them in, at this
 * time.  That decomposition drifted every particle to this time, wrapped them into the box and measured
 * the extent around them before ordering them, so each is current and inside the key extent: a segment
 * built over them now needs no keys and no sort.  Order is a performance property of the tiles, never a
 * correctness one. */
static int sidx_particles_still_ordered(int owned_end, integertime t_now)
{
    return g_sidx_particles_ordered.valid && g_sidx_particles_ordered.owned_epoch == g_sidx_owned_epoch &&
           g_sidx_particles_ordered.ti == t_now && g_sidx_particles_ordered.num_part == owned_end;
}


/* A new sync point.  The gas index is KEPT: it describes its members as of its reference time, and
 * every walk reads it at the time of the search (sfc_tiles.h), so a drift needs nothing here.  The
 * all-types index is not kept -- only the gas index has its bounds raised when its members are kicked
 * or their motion is written -- so it is released, and the first sink call of the sync point rebuilds it. */
void gpu_step_sidx_invalidate(void)
{
    gpu_spatial_index_free(&g_step_sidx_alltypes);
}

void gpu_step_sidx_invalidate_full(void)
{
    gpu_spatial_index_free(&g_step_sidx_alltypes);
    gpu_spatial_index_free(&g_step_sidx);
}


/* How far a tile's box can grow per unit drift interval, relative to its own size or reach.  Once that
 * growth times the interval since the index was built passes SIDX_MAX_BOX_GROWTH, searching the kept
 * index costs more than rebuilding it.  Performance only: a value that lags a raise delays a rebuild and
 * never admits a wrong answer. */
static constexpr double SIDX_MAX_BOX_GROWTH = 1.0;

/* The search size, as a fraction of the segment's source particles (the rank's own particles it is built
 * over, owned_end; the 1% was measured with this denominator), above which a segment whose members are all
 * current is rebuilt rather than read in its drifted frame.  A drifted-frame walk costs several times an
 * exact one per query (a motion envelope per node and per candidate), while a build costs per member;
 * measured, the two balance near a search of one percent of the source range.  An internal performance
 * heuristic: either choice gives the same list. */
static constexpr double SIDX_FRESH_SEARCH_FRACTION = 0.01;

KOKKOS_INLINE_FUNCTION
double sidx_tile_looseness(const sfc_tile_t &tile)
{
    double w = 0;
    for(int k = 0; k < 3; k++) {
        const double growth = (tile.u_max[k] - tile.u_min[k]) + 2.0 * tile.rho;
        if(!(growth > 0)) {continue;}
        double size = tile.hi[k] - tile.lo[k];
        if(tile.hmax > size) {size = tile.hmax;}
        const double r = (size > 0) ? growth / size : MAX_REAL_NUMBER;
        if(r > w) {w = r;}
    }
    return w;
}

/* Release one segment.  Its own dirty-tracker registration goes with it, so releasing one segment never
 * touches another's. */
static void sidx_segment_free(gpu_index_segment_t *seg)
{
    /* Ordering, not cleanup. kokkos_free does not synchronize, so releasing a
     * device allocation while a kernel may still be reading it is a
     * use-after-free. The fence lives HERE, with the release, so that no caller
     * can omit it -- this function is reached from the step loop, the
     * decomposition boundary, the ghost cleanup and the cached-index staleness
     * guard, and putting the rule in any one of those leaves the next caller free
     * to reintroduce the hazard. Skipped when there is nothing device-side to release. */
    if(seg->d_compact_xyzh || seg->d_pool || seg->d_bvh || seg->d_tiles || seg->d_slot_of || seg->d_level_nodes ||
       seg->d_shear_folds || seg->looseness) {
        Kokkos::fence();
    }
    if(seg->d_compact_xyzh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_compact_xyzh);}
    if(seg->d_pool) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_pool);}
    if(seg->d_bvh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_bvh);}
    if(seg->d_tiles) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_tiles);}
    if(seg->d_slot_of) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_slot_of);}
    if(seg->d_level_nodes) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_level_nodes);}
    if(seg->d_shear_folds) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(seg->d_shear_folds);}
    if(seg->h_level_offsets) {Kokkos::kokkos_free<Kokkos::HostSpace>(seg->h_level_offsets);}
    if(seg->looseness) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(seg->looseness);}
    if(seg->dirty_handle >= 0) {gpu_dirty_tracker_unregister(seg->dirty_handle);}
    *seg = gpu_index_segment_t{};
}

void gpu_spatial_index_free(gpu_spatial_index_t *idx)
{
    sidx_segment_free(&idx->ghost);
    sidx_segment_free(&idx->owned);
    idx->cache_tbm = -1;
    idx->cache_radius_policy = MODE_B_RADIUS_DEFAULT;
}

void gpu_sidx_ghost_pool_cleanup(void)
{
    sidx_segment_free(&g_step_sidx.ghost);
    gpu_spatial_index_free(&g_step_sidx_alltypes);
}


/* Exhaustion is reported by returning NULL, so the checks below are live code and
 * the run stops cleanly instead of aborting mid-flight. A zero-byte request is
 * treated as nothing to allocate, which the callers here read as a refusal. */
static void *ngl_alloc_shared(size_t bytes, const char *label)
{
    if(bytes == 0) {return NULL;}
    return gizmo_gpu_alloc_shared(bytes, label);
}
static void *ngl_alloc_device(size_t bytes, const char *label)
{
    if(bytes == 0) {return NULL;}
    return gizmo_gpu_alloc_device(bytes, label);
}


/* ---- Keeping the kept index's bounds true (sfc_tiles.h) ----------------------------------------
 * A member's velocity range is re-read when its velocity changes, and its reach when its radius does;
 * the tile and the nodes above it are widened to cover the new value.  Bounds only ever widen.  A
 * raise costs the members it is handed and the paths from their tiles to the root -- or, when those
 * paths together cost more than one sweep of the index, one sweep of it, level by level. */

/* What one member contributes to the bounds of its tile and of the nodes above it. */
struct SidxRaise {
    double u_min[3], u_max[3];
    double rho;
    double hmax;
    double hmax_by_type[TILE_NUM_PTYPES];
};

/* Widen the bounds of a tile or node to cover `m` (a member, a tile or a node).  Called concurrently
 * -- members of one tile, walks up shared ancestors -- so each field is raised atomically; the plain
 * read first skips the atomic when the bound already covers the value.  Returns 1 if a bound moved. */
template <class Box, class Src>
KOKKOS_INLINE_FUNCTION
int sidx_widen(Box *b, const Src &m)
{
    int moved = 0;
    for(int k = 0; k < 3; k++) {
        if(m.u_min[k] < b->u_min[k] && Kokkos::atomic_fetch_min(&b->u_min[k], m.u_min[k]) > m.u_min[k]) {moved = 1;}
        if(m.u_max[k] > b->u_max[k] && Kokkos::atomic_fetch_max(&b->u_max[k], m.u_max[k]) < m.u_max[k]) {moved = 1;}
    }
    if(m.rho  > b->rho  && Kokkos::atomic_fetch_max(&b->rho,  m.rho)  < m.rho)  {moved = 1;}
    if(m.hmax > b->hmax && Kokkos::atomic_fetch_max(&b->hmax, m.hmax) < m.hmax) {moved = 1;}
    for(int t = 0; t < TILE_NUM_PTYPES; t++) {
        if(m.hmax_by_type[t] > b->hmax_by_type[t] &&
           Kokkos::atomic_fetch_max(&b->hmax_by_type[t], m.hmax_by_type[t]) < m.hmax_by_type[t]) {moved = 1;}
    }
    return moved;
}

/* The same widening by a single writer: each bound is raised in place, with no atomics. */
template <class Box, class Src>
KOKKOS_INLINE_FUNCTION
void sidx_widen_plain(Box *b, const Src &m)
{
    for(int k = 0; k < 3; k++) {
        if(m.u_min[k] < b->u_min[k]) {b->u_min[k] = m.u_min[k];}
        if(m.u_max[k] > b->u_max[k]) {b->u_max[k] = m.u_max[k];}
    }
    if(m.rho  > b->rho)  {b->rho  = m.rho;}
    if(m.hmax > b->hmax) {b->hmax = m.hmax;}
    for(int t = 0; t < TILE_NUM_PTYPES; t++) {if(m.hmax_by_type[t] > b->hmax_by_type[t]) {b->hmax_by_type[t] = m.hmax_by_type[t];}}
}

/* Re-establish, level by level from the leaves, that every node covers its children.  Each node of a
 * level has one writer, and its children were finished at the level before, so the union is plain. */
static void sidx_widen_all_levels(gpu_index_segment_t *seg)
{
    const sfc_tile_t *tiles = seg->d_tiles;
    tile_bvh_node_t *bvh = seg->d_bvh;
    const int *level_nodes = seg->d_level_nodes;
    for(int L = 0; L < seg->nlevels; L++) {
        /* each level reads the one below it, so it waits for it */
        Kokkos::parallel_for("sidx_widen_level",
                             Kokkos::RangePolicy<>(seg->h_level_offsets[L], seg->h_level_offsets[L + 1]),
                             KOKKOS_LAMBDA(int q) {
            tile_bvh_node_t *node = &bvh[level_nodes[q]];
            if(node->left < 0) {sidx_widen_plain(node, tiles[-(node->left + 1)]);}
            else {sidx_widen_plain(node, bvh[node->left]); sidx_widen_plain(node, bvh[node->right]);}
        });
        Kokkos::fence();
    }
    gizmo_gpu_check_last_error("sidx_widen_level", seg->nlevels);
}

/* ======================================================================================================
 * Building a segment of the index.
 *
 * One builder for every route, templated on where it runs: the host (its working space in the memory
 * arena, then one copy of the finished index to the device) or the device (built in place).  The steps:
 *   1. the members of the source range, compacted in particle order by a prefix sum;
 *   2. one key per member: its position at the reference time, in its primary-box image, on the Morton
 *      curve, with members outside the key extent in a class of their own after the rest;
 *   3. members and keys sorted together (only the 16-byte key and the member index move);
 *   4. one team per tile folds its members, recomputing each from the particle data, into the rows, the
 *      membership pool and the tile bounds;
 *   5. the BVH: its shape depends only on the number of tiles (midpoint split of the tile range), and its
 *      bounds are the union of its children, level by level from the leaves.
 * Two positions per member, for two jobs.  The ROW is where the drift puts it, unwrapped -- the position
 * the pair test reads, so an exact frame's list is exact.  The tile and BVH boxes hold its primary-box
 * IMAGE, so a tile stays compact however many of its members have crossed a box face since the last
 * wrapping; the walk tests boxes and rows through the canonical wrap either way.  Members outside the
 * key extent start at a fresh tile, so they never widen the tile of a member inside it.
 * ====================================================================================================== */

/* Where member j is at the reference time ti_ref, unwrapped: the drift's position (particle_motion_envelope)
 * and how far from it the member can be.  A member the envelope will not predict because a reflecting or
 * outflow wall may turn it round is bounded instead by its speed from where it stands.  Returns 1 when its
 * clock cannot be advanced to ti_ref or its position is not finite: no search can bound it. */
KOKKOS_INLINE_FUNCTION
int sidx_member_position(int j, const struct particle_data *P, const struct gas_cell_data *cells, integertime ti_ref,
                         const struct DriftKickTableView &tables, double center[3], double *hw, int *current)
{
    const integertime ti_j = P[j].Ti_current;
    if(ti_j < 0 || ti_j > ti_ref) {return 1;}
    const int motion = particle_motion_envelope(j, P, cells, ti_ref, &tables, center, hw);
    if(motion == PARTICLE_MOTION_UNBOUNDED) {
        const double dl = motion_bound_widening(particle_motion_speed_bound(j, P, cells), ti_j, ti_ref, &tables);
        if(!motion_bound_widening_is_valid(dl)) {return 1;}
        for(int k = 0; k < 3; k++) {center[k] = (double)P[j].Pos[k];}
        *hw = 0.5 * dl + motion_envelope_rounding_floor(center, 0.5 * dl);
    }
    for(int k = 0; k < 3; k++) {if(!(center[k] - center[k] == 0.0)) {return 1;}}   /* NaN or Inf, fast-math safe */
    *current = (motion == PARTICLE_MOTION_CURRENT);
    return 0;
}

/* The primary-box image of a member's unwrapped position, and the half-width its box needs there: its own,
 * plus the rounding of the fold, at the largest scale the fold worked at -- the unwrapped and the folded
 * coordinates, the box lengths and the shearing offset -- never at the folded result alone, which can be
 * tiny where the fold rounded at the box length. */
KOKKOS_INLINE_FUNCTION
void sidx_member_image(const double center[3], double hw, double image[3], double *box_hw, int *folds_up, int *folds_down)
{
    Vec3<MyDouble> img = {{center[0], center[1], center[2]}};
    const struct box_wrap_x_folds folds = box_wrap_position_to_primary_image(img);
    *folds_up = folds.up_from_below; *folds_down = folds.down_from_above;
    int moved = 0;
    for(int k = 0; k < 3; k++) {image[k] = img[k]; if(image[k] != center[k]) {moved = 1;}}
    *box_hw = hw;
    if(moved) {
        double scale = hw;
        for(int k = 0; k < 3; k++) {if(fabs(center[k]) > scale) {scale = fabs(center[k]);}}
#if defined(BOX_PERIODIC)
        if(boxSize_X > scale) {scale = boxSize_X;}
        if(boxSize_Y > scale) {scale = boxSize_Y;}
        if(boxSize_Z > scale) {scale = boxSize_Z;}
#endif
#if defined(BOX_SHEARING) && (BOX_SHEARING > 1)
        if(fabs(Shearing_Box_Pos_Offset) > scale) {scale = fabs(Shearing_Box_Pos_Offset);}
#endif
        *box_hw += motion_envelope_rounding_floor(image, scale);
    }
}

/* A member's velocity range must hold its motion in both frames the index reads it in: unwrapped, where the
 * leaf moves its row, and in its primary-box image, where its tile's box moves.  The two differ only in a
 * shearing box with a position offset (BOX_SHEARING > 1): there the image of a member folded across an x face
 * moves, along the shearing coordinate, at its velocity shifted by the shearing velocity offset once per fold,
 * up folds first -- the shift the box wrapping gives the velocity itself (do_box_wrapping).  So the range is
 * widened to cover the shifted velocity as well. */
KOKKOS_INLINE_FUNCTION
void sidx_image_velocity_union(double u_lo[3], double u_hi[3], int folds_up, int folds_down)
{
#if defined(BOX_SHEARING) && (BOX_SHEARING > 1)
    const int k = BOX_SHEARING_PHI_COORDINATE;
    double lo = u_lo[k], hi = u_hi[k];
    for(int f = 0; f < folds_up; f++) {lo -= Shearing_Box_Vel_Offset; hi -= Shearing_Box_Vel_Offset;}
    for(int f = 0; f < folds_down; f++) {lo += Shearing_Box_Vel_Offset; hi += Shearing_Box_Vel_Offset;}
    if(lo < u_lo[k]) {u_lo[k] = lo;}
    if(hi > u_hi[k]) {u_hi[k] = hi;}
#else
    (void)u_lo; (void)u_hi; (void)folds_up; (void)folds_down;
#endif
}

/* What the index holds about one member at its reference time ti_ref. */
struct SidxMember {
    double center[3];       /* the row: unwrapped, where the drift puts it at ti_ref */
    double hw;              /* how far from there it can be */
    double image[3];        /* its primary-box image, which its tile's box holds */
    double box_hw;          /* the half-width its box needs at the image */
    double u_lo[3], u_hi[3];/* the velocity range it can advance at, in both frames (sidx_image_velocity_union) */
    double rho;             /* its residual speed */
    double r, r_drifted;    /* its reach at its own clock, and the most it can be once drifted */
    int type;
    int current;            /* already at ti_ref */
    int folds_up, folds_down;
};

/* Member j's reach at its own clock, and the most it can be once drifted (never below the reach). */
KOKKOS_INLINE_FUNCTION
void sidx_member_reaches(int j, struct particle_data *P, mode_b_radius_policy_t policy, double growth, double kernel_floor,
                         double *r, double *r_drifted)
{
    *r = nlr_particle_symmetric_radius(P[j], policy);
    *r_drifted = nlr_particle_symmetric_radius_after_drift_P(j, P, growth, kernel_floor, policy);
    if(*r_drifted < *r) {*r_drifted = *r;}
}

/* What the index holds about a member that is read only at its reference time: an imported particle, which
 * is at that time when it arrives and whose segment never outlives it.  Its position (and the half-width the
 * envelope gave it, zero for a particle already at that time), its primary-box image with the rounding of
 * the fold, its reaches and its type -- and no motion: a search at the reference time moves nothing. */
KOKKOS_INLINE_FUNCTION
void sidx_describe_exact_member(const double center[3], double hw, int current, double r, double r_drifted, int type,
                                struct SidxMember &m)
{
    for(int k = 0; k < 3; k++) {m.center[k] = center[k]; m.u_lo[k] = 0.0; m.u_hi[k] = 0.0;}
    m.hw = hw;
    m.rho = 0.0;
    sidx_member_image(m.center, m.hw, m.image, &m.box_hw, &m.folds_up, &m.folds_down);
    m.r = r;
    m.r_drifted = (r_drifted < r) ? r : r_drifted;
    m.type = type;
    m.current = current;
}

/* The one rule for what the index holds about member j; exact: it is read only at the reference time
 * (sidx_describe_exact_member).  Returns 1 when no search can bound it. */
KOKKOS_INLINE_FUNCTION
int sidx_describe_member(int j, struct particle_data *P, const struct gas_cell_data *cells,
                         mode_b_radius_policy_t policy, integertime ti_ref,
                         const struct DriftKickTableView &tables, double growth, double kernel_floor, int exact,
                         struct SidxMember &m)
{
    double center[3], hw = 0.0, r, r_drifted;
    int current = 0;
    if(sidx_member_position(j, P, cells, ti_ref, tables, center, &hw, &current)) {return 1;}
    /* A member whose motion cannot be bounded is refused whether or not the index reads its motion. */
    if(sfc_member_motion_range(j, P, cells, m.u_lo, m.u_hi, &m.rho)) {return 1;}
    sidx_member_reaches(j, P, policy, growth, kernel_floor, &r, &r_drifted);
    if(exact) {sidx_describe_exact_member(center, hw, current, r, r_drifted, (int)P[j].Type, m); return 0;}
    for(int k = 0; k < 3; k++) {m.center[k] = center[k];}
    m.hw = hw; m.current = current;
    sidx_member_image(m.center, m.hw, m.image, &m.box_hw, &m.folds_up, &m.folds_down);
    sidx_image_velocity_union(m.u_lo, m.u_hi, m.folds_up, m.folds_down);
    m.r = r; m.r_drifted = r_drifted;
    m.type = (int)P[j].Type;
    return 0;
}

/* Where an index's members come from: the particles [base, base + count) of P[], read where they live.  A
 * member is named by its ordinal o in that range -- the slot map is indexed by ordinal -- and by its global
 * particle index base + o everywhere else (the pool, the neighbour lists).  Particles below owned_end are the
 * rank's own.  exact: the members are read only at the reference time (the imported particles). */
struct SidxParticleSource {
    struct particle_data *P;
    const struct gas_cell_data *cells;
    int base, count, owned_end, type_bitmask, exact;
    KOKKOS_INLINE_FUNCTION int global(int o) const {return base + o;}
    KOKKOS_INLINE_FUNCTION int is_member(int o) const {return sfc_pool_member(&P[base + o], type_bitmask);}
    KOKKOS_INLINE_FUNCTION int is_owned(int o) const {return base + o < owned_end;}
    KOKKOS_INLINE_FUNCTION int position(int o, integertime ti_ref, const struct DriftKickTableView &tables,
                                        double center[3], double *hw, int *current) const
    {return sidx_member_position(base + o, P, cells, ti_ref, tables, center, hw, current);}
    KOKKOS_INLINE_FUNCTION int describe(int o, mode_b_radius_policy_t policy, integertime ti_ref,
                                        const struct DriftKickTableView &tables, double growth, double kernel_floor,
                                        struct SidxMember &m) const
    {return sidx_describe_member(base + o, P, cells, policy, ti_ref, tables, growth, kernel_floor, exact, m);}
};

/* The most tiles num_members members can fill: the members inside the key extent fill tiles from the first,
 * and those outside start a fresh one. */
static int sidx_tiles_bound(int num_members)
{
    return num_members / TILE_TARGET_SIZE + (num_members % TILE_TARGET_SIZE != 0) + 1;
}

/* A BVH over at most INT_MAX tiles has at most this many levels (sidx_bvh_levels). */
static constexpr int SIDX_BVH_MAX_LEVELS = 32;

/* The levels of the BVH over ntiles tiles (sidx_bvh_shape: a midpoint split, so ceil(log2(ntiles)) + 1). */
static int sidx_bvh_levels(int ntiles)
{
    int levels = 1;
    while(((long long)1 << (levels - 1)) < (long long)ntiles) {levels++;}
    return levels;
}

/* The compact records of count imported particles (SidxRecordSource), at least one: rows, types, states. */
static void sidx_record_bytes(int count, size_t *rows, size_t *types, size_t *states)
{
    const size_t n = (size_t)(count > 0 ? count : 1);
    *rows = n * SIDX_ROW_WIDTH * sizeof(double); *types = n * sizeof(signed char); *states = n * sizeof(signed char);
}

/* Every allocation one index build makes, in bytes: what the segment keeps (on the device), the build's
 * working space, and the compact records of a build that reads the imported particles away from the particle
 * arrays.  The builder takes each size from here, so what is projected is what is taken.  Sized for ntiles
 * tiles: sidx_tiles_bound(num_members) before the members are sorted into tiles, the exact count after. */
struct SidxBuildMemoryPlan {
    int fits;   /* 0: the counts exceed what one segment can index, and no size below is set */
    int ntiles, nnodes, nslots, nlevels;
    /* kept by the segment: what a walk reads, then what keeps a maintained segment's bounds true */
    size_t tiles, bvh, pool, rows, level_nodes, slot_of, shear_folds, level_offsets, looseness;
    /* the build's working space; bvh_shape is on the host on either route (sidx_bvh_shape) */
    size_t scan, members, keys, leaf_of_tile, bvh_shape;
    /* one record per imported particle: its row, its type, its state */
    size_t record_rows, record_types, record_states;
};

static struct SidxBuildMemoryPlan sidx_build_memory_plan(int num_source, int num_members, int ntiles,
                                                         int maintained, int presorted)
{
    struct SidxBuildMemoryPlan m = {};
    const size_t n_source = (size_t)(num_source > 0 ? num_source : 1);
    const size_t n_mem = (size_t)(num_members > 0 ? num_members : 1);
    /* Node and slot counts are ints throughout the index; every size below is such a count, or the source
     * count, times a structure size, so none overflows once these fit. */
    if(ntiles < 1 || 2LL * ntiles - 1 > INT_MAX || (long long)ntiles * TILE_TARGET_SIZE > INT_MAX) {return m;}
    m.fits = 1;
    m.ntiles = ntiles; m.nnodes = 2 * ntiles - 1; m.nslots = ntiles * TILE_TARGET_SIZE; m.nlevels = sidx_bvh_levels(ntiles);
    m.tiles = (size_t)m.ntiles * sizeof(sfc_tile_t);
    m.bvh = (size_t)m.nnodes * sizeof(tile_bvh_node_t);
    m.pool = (size_t)m.nslots * sizeof(int);
    m.rows = (size_t)m.nslots * SIDX_ROW_WIDTH * sizeof(double);
    m.level_nodes = (size_t)m.nnodes * sizeof(int);
    if(maintained) {
        m.slot_of = n_source * sizeof(int);
#if defined(BOX_SHEARING) && (BOX_SHEARING > 1)
        m.shear_folds = (size_t)m.nslots * 2 * sizeof(int);
#endif
        m.level_offsets = (size_t)(m.nlevels + 1) * sizeof(int);
        m.looseness = sizeof(double);
    }
    m.scan = ((size_t)num_source + 1) * sizeof(int);
    m.members = n_mem * sizeof(int);
    m.keys = presorted ? 0 : n_mem * sizeof(Morton128);
    m.leaf_of_tile = (size_t)m.ntiles * sizeof(int);
    m.bvh_shape = (size_t)m.nnodes * sizeof(tile_bvh_node_t) + ((size_t)m.ntiles + 2 * (size_t)m.nnodes + (size_t)m.nlevels) * sizeof(int);
    sidx_record_bytes(num_source, &m.record_rows, &m.record_types, &m.record_states);
    return m;
}

/* The imported particles [base, base + count) as compact records, for a build that runs where the particle
 * arrays are not.  Per particle a row (position at the reference time, reach, drifted reach), its type (-1: not
 * a member under the index's mask) and its state, all packed on the host by the particle source's own rule
 * (sidx_describe_member, exact): so the build refuses a member exactly when it would refuse the particle, and
 * reports a member not at the reference time as the particle build does.  Exact members only. */
enum { SIDX_RECORD_CURRENT = 0, SIDX_RECORD_NOT_CURRENT = 1, SIDX_RECORD_REFUSED = 2 };
struct SidxRecordSource {
    const double *row;
    const signed char *type;
    const signed char *state;
    int base, count;
    KOKKOS_INLINE_FUNCTION int global(int o) const {return base + o;}
    KOKKOS_INLINE_FUNCTION int is_member(int o) const {return type[o] >= 0;}
    KOKKOS_INLINE_FUNCTION int is_owned(int) const {return 0;}
    KOKKOS_INLINE_FUNCTION int position(int o, integertime, const struct DriftKickTableView &,
                                        double center[3], double *hw, int *current) const
    {
        if(state[o] == SIDX_RECORD_REFUSED) {return 1;}
        for(int k = 0; k < 3; k++) {center[k] = row[(size_t)o * SIDX_ROW_WIDTH + k];}
        *hw = 0.0; *current = (state[o] == SIDX_RECORD_CURRENT);
        return 0;
    }
    KOKKOS_INLINE_FUNCTION int describe(int o, mode_b_radius_policy_t, integertime ti_ref,
                                        const struct DriftKickTableView &tables, double, double,
                                        struct SidxMember &m) const
    {
        double center[3], hw; int current;
        if(position(o, ti_ref, tables, center, &hw, &current)) {return 1;}
        const double *x = &row[(size_t)o * SIDX_ROW_WIDTH];
        sidx_describe_exact_member(center, hw, current, x[3], x[4], (int)type[o], m);
        return 0;
    }
};

/* The sort's own working space, per member sorted (measured with the OpenMP backend).  An allowance for
 * projecting a build's peak, and the size reported, as an estimate, when the sort is refused its memory. */
static constexpr size_t SIDX_SORT_BYTES_PER_MEMBER = 20;

/* What a build reports besides the index itself. */
struct SidxBuildReport {
    int refused;          /* a member no search can bound: the index is left unbuilt */
    int all_current;      /* every member was at ti_ref: the rows are positions, not predictions */
    int n_outside;        /* members outside the key extent */
    int n_outside_owned;  /* of those, the rank's own particles */
    size_t bytes_failed;  /* the request that could not be had, when one could not */
    int failed_memory;    /* where that request was made (SidxMemory) */
    int failed_estimated; /* bytes_failed is an estimate (the sort sizes its own request) */
    int failed_kept;      /* that request was for the arrays the segment keeps, which every route places on the device */
    char failure[160];    /* what failed, when it was not a request of known size */
};

/* How a build ended.  An allocation refused is a request of the build's own that an allocator positively
 * declined: the arena's preflight, or an allocator that returns null.  That is not proof the memory was
 * exhausted (the device allocator returns null for any error it meets), only that the request was not had.
 * Anything else that goes wrong is a failure. */
enum SidxBuildStatus {SIDX_BUILT = 0, SIDX_REFUSED_MEMBER = 1, SIDX_ALLOCATION_REFUSED = 2, SIDX_FAILED = 3};

/* Where a build's request is made. */
enum SidxMemory {SIDX_MEM_DEVICE = 0, SIDX_MEM_ARENA = 1, SIDX_MEM_HOST = 2, SIDX_MEM_SHARED = 3};
static const char *sidx_memory_name(int memory)
{
    return memory == SIDX_MEM_ARENA ? "the memory arena" : memory == SIDX_MEM_HOST ? "host memory"
         : memory == SIDX_MEM_SHARED ? "shared memory" : "device memory";
}

/* Thrown by the build when a request of its own is refused; the build then returns SIDX_ALLOCATION_REFUSED. */
struct SidxAllocationRefused {
    size_t bytes; int memory; int estimated = 0;
    /* constructed, not brace-initialised: nvcc's front end fails an internal assertion on a two-member aggregate
     * thrown here */
    SidxAllocationRefused(size_t b, int mem, int est = 0) : bytes(b), memory(mem), estimated(est) {}
};

/* The build's working space, released when the build returns, on every path.  On the host it comes from
 * the memory arena (so the memory ledger and the arena sizing see it, as they saw the host build this
 * replaces) and is released in reverse; on the device, from the device allocator. */
template <bool kStageOnHost>
struct SidxScratch {
    static constexpr int CAPACITY = 16;   /* more than the build ever holds at once */
    void *held[CAPACITY];
    int n = 0;
    void *take(size_t bytes, const char *label) {
        if(n >= CAPACITY) {throw std::runtime_error("index build working space: more blocks than it holds");}
        if(bytes == 0) {bytes = 1;}
        void *p;
        /* Every request asks the arena first, so a build that does not fit stops cleanly instead of failing
         * inside the allocator; what it already holds is given back on the way out. */
        if(kStageOnHost && !gizmo_alloc_fits_this_rank(gizmo_mymalloc_rounded_size(bytes), 1)) {throw SidxAllocationRefused(bytes, SIDX_MEM_ARENA);}
        if(kStageOnHost) {p = mymalloc(label, bytes);}
        else if(!(p = gizmo_gpu_alloc_device(bytes, label))) {throw SidxAllocationRefused(bytes, SIDX_MEM_DEVICE);}
        held[n++] = p;
        return p;
    }
    int mark() const {return n;}
    void release_to(int m) {
        if(!kStageOnHost && n > m) {Kokkos::fence();}   /* no kernel may still read what is given back (sidx_segment_free) */
        while(n > m) {
            void *p = held[--n];
            if(kStageOnHost) {myfree(p);} else {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(p);}
        }
    }
    ~SidxScratch() {release_to(0);}
};

/* The records of SidxRecordSource, on the device: pack() reads P[] on the host, staging the records in the
 * memory arena, and copies them once; the device memory is given back by release() or when the object goes.
 * pack returns 0; 1 when a request could not be had (the report says which). */
struct SidxGhostRecords {
    double *d_row = NULL;
    signed char *d_type = NULL, *d_state = NULL;
    int base = 0, count = 0;
    SidxGhostRecords() = default;
    SidxGhostRecords(const SidxGhostRecords &) = delete;
    SidxGhostRecords &operator=(const SidxGhostRecords &) = delete;
    ~SidxGhostRecords() {release();}
    int pack(const struct SidxParticleSource &src, mode_b_radius_policy_t policy, integertime ti_ref,
             const struct DriftKickTableView &tables, struct SidxBuildReport *report)
    {
        release();
        base = src.base; count = src.count;
        const size_t n = (size_t)(count > 0 ? count : 1);   /* one record per particle, at least one */
        size_t row_bytes, type_bytes, state_bytes;
        sidx_record_bytes(count, &row_bytes, &type_bytes, &state_bytes);
        try {
            SidxScratch<true> stage;   /* given back to the arena, in reverse, on every way out */
            double *h_row = (double *) stage.take(row_bytes, "ngl_sidx_ghost_stage_rows");
            signed char *h_type = (signed char *) stage.take(type_bytes, "ngl_sidx_ghost_stage_types");
            signed char *h_state = (signed char *) stage.take(state_bytes, "ngl_sidx_ghost_stage_states");
            const double growth = kernel_radius_drift_max_growth_factor(), kernel_floor = All.MinKernelRadius;
#pragma omp parallel for
            for(size_t o = 0; o < n; o++) {
                double *x = &h_row[o * SIDX_ROW_WIDTH];
                for(int k = 0; k < SIDX_ROW_WIDTH; k++) {x[k] = 0.0;}
                h_type[o] = (signed char)-1; h_state[o] = (signed char)SIDX_RECORD_REFUSED;
                const int j = src.base + (int)o;
                if((int)o >= count || !sfc_pool_member(&src.P[j], src.type_bitmask)) {continue;}
                h_type[o] = (signed char)src.P[j].Type;
                struct SidxMember m;
                if(sidx_describe_member(j, src.P, src.cells, policy, ti_ref, tables, growth, kernel_floor, 1, m)) {continue;}
                x[0] = m.center[0]; x[1] = m.center[1]; x[2] = m.center[2]; x[3] = m.r; x[4] = m.r_drifted;
                h_state[o] = (signed char)(m.current ? SIDX_RECORD_CURRENT : SIDX_RECORD_NOT_CURRENT);
            }
            using UV = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
            if(!(d_row = (double *) gizmo_gpu_alloc_device(row_bytes, "ngl_sidx_ghost_record_rows"))) {throw SidxAllocationRefused(row_bytes, SIDX_MEM_DEVICE);}
            if(!(d_type = (signed char *) gizmo_gpu_alloc_device(type_bytes, "ngl_sidx_ghost_record_types"))) {throw SidxAllocationRefused(type_bytes, SIDX_MEM_DEVICE);}
            if(!(d_state = (signed char *) gizmo_gpu_alloc_device(state_bytes, "ngl_sidx_ghost_record_states"))) {throw SidxAllocationRefused(state_bytes, SIDX_MEM_DEVICE);}
            Kokkos::deep_copy(Kokkos::View<double*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(d_row, n * SIDX_ROW_WIDTH),
                              Kokkos::View<const double*, Kokkos::HostSpace, UV>(h_row, n * SIDX_ROW_WIDTH));
            Kokkos::deep_copy(Kokkos::View<signed char*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(d_type, n),
                              Kokkos::View<const signed char*, Kokkos::HostSpace, UV>(h_type, n));
            Kokkos::deep_copy(Kokkos::View<signed char*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(d_state, n),
                              Kokkos::View<const signed char*, Kokkos::HostSpace, UV>(h_state, n));
        } catch(const SidxAllocationRefused &e) {
            report->bytes_failed = e.bytes; report->failed_memory = e.memory; report->failed_estimated = e.estimated;
            release();
            return 1;
        }
        return 0;
    }
    struct SidxRecordSource source() const {return {d_row, d_type, d_state, base, count};}
    void release()
    {
        base = 0; count = 0;   /* nothing to read, whether or not anything was allocated */
        if(!d_row && !d_type && !d_state) {return;}
        Kokkos::fence();   /* no kernel may still read them */
        if(d_row) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_row); d_row = NULL;}
        if(d_type) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_type); d_type = NULL;}
        if(d_state) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_state); d_state = NULL;}
    }
};

/* Widen a box (a tile or a node) to cover `s`, its position bounds as well as what sidx_widen_plain covers,
 * by a single writer. */
template <class Box, class Src>
KOKKOS_INLINE_FUNCTION
void sidx_union_plain(Box *b, const Src &s)
{
    for(int k = 0; k < 3; k++) {
        if(s.lo[k] < b->lo[k]) {b->lo[k] = s.lo[k];}
        if(s.hi[k] > b->hi[k]) {b->hi[k] = s.hi[k];}
    }
    sidx_widen_plain(b, s);
}

/* A tile or node that covers nothing: bounds inverted, so a union is neutral and a gap test always fails. */
template <class Box>
KOKKOS_INLINE_FUNCTION
void sidx_box_empty(Box *b)
{
    for(int k = 0; k < 3; k++) {
        b->lo[k] = MAX_REAL_NUMBER; b->hi[k] = -MAX_REAL_NUMBER;
        b->u_min[k] = MAX_REAL_NUMBER; b->u_max[k] = -MAX_REAL_NUMBER;
    }
    b->rho = 0; b->hmax = 0;
    for(int t = 0; t < TILE_NUM_PTYPES; t++) {b->hmax_by_type[t] = 0;}
}

/* A tile's bounds as its members' team builds them: each member contributes its own, the team joins them,
 * and one lane writes the result into the tile.  Also counts what the build reports about the members. */
struct SidxTileSummary {
    double lo[3], hi[3], u_min[3], u_max[3], rho, hmax, hmax_by_type[TILE_NUM_PTYPES], hw;
    int refused;       /* members no search can bound */
    int not_current;   /* members not at the reference time */
};

/* The join of tile summaries, as a Kokkos reducer: the union of the bounds, the sum of the counts. */
struct SidxTileJoin {
    using reducer = SidxTileJoin;
    using value_type = SidxTileSummary;
    using result_view_type = Kokkos::View<value_type, Kokkos::AnonymousSpace, Kokkos::MemoryUnmanaged>;
    /* The lanes exchange the summary by shuffles of int-sized pieces (sizeof / sizeof(int), truncating), so a
     * summary that is not a whole number of them would lose its tail on the device only (NlrAccumReducer). */
    static_assert(sizeof(value_type) % sizeof(int) == 0, "SidxTileSummary must be a whole number of int-sized pieces");
    value_type &value;
    KOKKOS_INLINE_FUNCTION SidxTileJoin(value_type &v) : value(v) {}
    KOKKOS_INLINE_FUNCTION void join(value_type &dst, const value_type &src) const
    {
        sidx_union_plain(&dst, src);
        if(src.hw > dst.hw) {dst.hw = src.hw;}
        dst.refused += src.refused;
        dst.not_current += src.not_current;
    }
    KOKKOS_INLINE_FUNCTION void init(value_type &v) const {sidx_box_empty(&v); v.hw = 0; v.refused = 0; v.not_current = 0;}
    KOKKOS_INLINE_FUNCTION value_type &reference() const {return value;}
    KOKKOS_INLINE_FUNCTION result_view_type view() const {return result_view_type(&value);}
    KOKKOS_INLINE_FUNCTION bool references_scalar() const {return true;}
};

/* Add one member's bounds to a tile summary: its box at its primary-box image, its velocity range and residual
 * speed, its drifted reach (overall and for its type) and its half-width. */
KOKKOS_INLINE_FUNCTION
void sidx_summary_add_member(struct SidxTileSummary &b, const struct SidxMember &m)
{
    for(int d = 0; d < 3; d++) {
        if(m.image[d] - m.box_hw < b.lo[d]) {b.lo[d] = m.image[d] - m.box_hw;}
        if(m.image[d] + m.box_hw > b.hi[d]) {b.hi[d] = m.image[d] + m.box_hw;}
        if(m.u_lo[d] < b.u_min[d]) {b.u_min[d] = m.u_lo[d];}
        if(m.u_hi[d] > b.u_max[d]) {b.u_max[d] = m.u_hi[d];}
    }
    if(m.rho > b.rho) {b.rho = m.rho;}
    if(m.r_drifted > b.hmax) {b.hmax = m.r_drifted;}
    if(m.type >= 0 && m.type < TILE_NUM_PTYPES && m.r_drifted > b.hmax_by_type[m.type]) {b.hmax_by_type[m.type] = m.r_drifted;}
    if(m.hw > b.hw) {b.hw = m.hw;}
}

/* What a build counts about its members, summed over the members of a pass. */
struct SidxBuildTally {
    int refused, not_current, outside, outside_owned;
    KOKKOS_INLINE_FUNCTION SidxBuildTally() : refused(0), not_current(0), outside(0), outside_owned(0) {}
    KOKKOS_INLINE_FUNCTION SidxBuildTally &operator+=(const SidxBuildTally &o)
    {refused += o.refused; not_current += o.not_current; outside += o.outside; outside_owned += o.outside_owned; return *this;}
};
namespace Kokkos {
template <> struct reduction_identity<SidxBuildTally> {
    KOKKOS_FORCEINLINE_FUNCTION static SidxBuildTally sum() {return SidxBuildTally();}
};
}

/* The class bit: members outside the key extent sort after every member inside it. */
static constexpr uint64_t SIDX_KEY_OUTSIDE = (uint64_t)1 << 63;

/* The extent keys are measured against: a cube from `corner` with side `len`. */
struct SidxKeyFrame {
    double corner[3];
    double len;
};

/* A member's key, from its primary-box image.  In a periodic box every image lies inside the box, so no
 * member is outside; otherwise outside means outside the extent the domain was built on, by the same
 * test the drift and the decomposition use (position_outside_domain_extent). */
KOKKOS_INLINE_FUNCTION
Morton128 sidx_member_key(const double image[3], const struct SidxKeyFrame &frame)
{
#if defined(BOX_PERIODIC)
    const int outside = 0;
#else
    const int outside = position_outside_domain_extent(image[0], image[1], image[2],
                                                       frame.corner[0], frame.corner[1], frame.corner[2], frame.len);
#endif
    uint64_t q[3];
    for(int k = 0; k < 3; k++) {
        double u = (frame.len > 0) ? (image[k] - frame.corner[k]) / frame.len + 1.0 : 1.0;
        if(!(u >= 1.0)) {u = 1.0;}
        if(!(u < 2.0))  {u = 1.9999999999999998;}
        q[k] = gpu_morton_double_to_int42(u);
    }
    Morton128 key = gpu_morton_encode128(q[0], q[1], q[2]);
    if(outside) {key.hi |= SIDX_KEY_OUTSIDE;}
    return key;
}

/* The shape of the BVH over ntiles tiles: a midpoint split of the tile range, nodes numbered children
 * first so the root is last.  Depends only on ntiles, so it is set out on the host: the nodes with their
 * links and empty bounds, which tile each leaf holds, and the nodes level by level.  Its arrays live in one
 * block (SidxBuildMemoryPlan::bvh_shape) laid out by sidx_bvh_shape; its level offsets are held here, since
 * the build reads them after the block is given back. */
struct SidxBvhShape {
    tile_bvh_node_t *nodes;
    int *leaf_of_tile, *level_nodes, *height, *cursor;
    int level_offsets[SIDX_BVH_MAX_LEVELS + 1];
    int nnodes, nlevels;
};

static int sidx_bvh_shape_recursive(SidxBvhShape &s, int t0, int t1)
{
    if(t1 - t0 == 1) {
        const int nd = s.nnodes++;
        sidx_box_empty(&s.nodes[nd]);
        s.nodes[nd].left = -(t0 + 1); s.nodes[nd].right = -(t0 + 1); s.nodes[nd].parent = -1; s.height[nd] = 0;
        s.leaf_of_tile[t0] = nd;
        return nd;
    }
    const int mid = (t0 + t1) / 2;
    const int l = sidx_bvh_shape_recursive(s, t0, mid);
    const int r = sidx_bvh_shape_recursive(s, mid, t1);
    const int nd = s.nnodes++;
    sidx_box_empty(&s.nodes[nd]);
    s.nodes[nd].left = l; s.nodes[nd].right = r; s.nodes[nd].parent = -1;
    s.nodes[l].parent = nd; s.nodes[r].parent = nd;
    s.height[nd] = 1 + std::max(s.height[l], s.height[r]);
    return nd;
}

/* Set out the shape for plan.ntiles tiles in block, of plan.bvh_shape bytes. */
static_assert(sizeof(tile_bvh_node_t) % alignof(int) == 0, "the shape's int arrays follow its nodes in one block");
static void sidx_bvh_shape(SidxBvhShape &s, const struct SidxBuildMemoryPlan &plan, void *block)
{
    const int ntiles = plan.ntiles, nnodes = plan.nnodes;
    s.nodes = (tile_bvh_node_t *)block;
    s.leaf_of_tile = (int *)(s.nodes + nnodes);
    s.level_nodes = s.leaf_of_tile + ntiles;
    s.height = s.level_nodes + nnodes;
    s.cursor = s.height + nnodes;
    s.nnodes = 0;
    sidx_bvh_shape_recursive(s, 0, ntiles);
    s.nlevels = 0;
    for(int nd = 0; nd < s.nnodes; nd++) {if(s.height[nd] + 1 > s.nlevels) {s.nlevels = s.height[nd] + 1;}}
    if(s.nnodes != nnodes || s.nlevels != plan.nlevels) {throw std::runtime_error("index build: BVH shape differs from its memory plan");}
    for(int L = 0; L <= s.nlevels; L++) {s.level_offsets[L] = 0;}
    for(int nd = 0; nd < s.nnodes; nd++) {s.level_offsets[s.height[nd] + 1]++;}
    for(int L = 0; L < s.nlevels; L++) {s.level_offsets[L + 1] += s.level_offsets[L];}
    for(int L = 0; L < s.nlevels; L++) {s.cursor[L] = s.level_offsets[L];}
    for(int nd = 0; nd < s.nnodes; nd++) {s.level_nodes[s.cursor[s.height[nd]]++] = nd;}
}

/* Device memory the segment a build keeps takes, on either route (sidx_alloc_kept allocates exactly this). */
static size_t sidx_kept_device_bytes(const struct SidxBuildMemoryPlan &m, int maintained)
{
    return m.tiles + m.bvh + m.pool + m.rows + (maintained ? m.slot_of + m.level_nodes + m.shear_folds : 0);
}

/* Allocate the arrays the segment keeps, on the device: what every segment needs to be walked, and, for a
 * maintained segment, what keeps its bounds true while it is kept.  Returns 0; 1 when one could not be had
 * (the segment is then released and the report names the request). */
static int sidx_alloc_kept(gpu_index_segment_t *seg, int maintained, const struct SidxBuildMemoryPlan &m,
                           struct SidxBuildReport *report)
{
    struct {void **dst; size_t bytes; const char *label; int maintenance;} kept[] = {
        {(void **)&seg->d_tiles, m.tiles, "ngl_sidx_dev_tiles", 0},
        {(void **)&seg->d_bvh, m.bvh, "ngl_sidx_dev_bvh", 0},
        {(void **)&seg->d_pool, m.pool, "ngl_sidx_dev_pool", 0},
        {(void **)&seg->d_compact_xyzh, m.rows, "ngl_sidx_dev_rows", 0},
        {(void **)&seg->d_slot_of, m.slot_of, "ngl_sidx_dev_slot_of", 1},
        {(void **)&seg->d_level_nodes, m.level_nodes, "ngl_sidx_dev_level_nodes", 1}};
    size_t device_bytes = 0;
    for(auto &k : kept) {
        if(k.maintenance && !maintained) {continue;}
        *k.dst = gizmo_gpu_alloc_device(k.bytes, k.label);
        if(!*k.dst) {report->bytes_failed = k.bytes; report->failed_memory = SIDX_MEM_DEVICE; report->failed_kept = 1; sidx_segment_free(seg); return 1;}
        device_bytes += k.bytes;
    }
    if(maintained && m.shear_folds > 0) {
        seg->d_shear_folds = (int *) gizmo_gpu_alloc_device(m.shear_folds, "ngl_sidx_dev_shear_folds");
        if(!seg->d_shear_folds) {report->bytes_failed = m.shear_folds; report->failed_memory = SIDX_MEM_DEVICE; report->failed_kept = 1; sidx_segment_free(seg); return 1;}
        device_bytes += m.shear_folds;
    }
    /* what the route selector projects for the kept segment is what is allocated here */
    if(device_bytes != sidx_kept_device_bytes(m, maintained)) {throw std::logic_error("index build: kept segment differs from its projection");}
    if(!maintained) {return 0;}
    seg->h_level_offsets = (int *) gizmo_gpu_alloc_host(m.level_offsets, "ngl_sidx_host_level_offsets");
    if(!seg->h_level_offsets) {report->bytes_failed = m.level_offsets; report->failed_memory = SIDX_MEM_HOST; report->failed_kept = 1; sidx_segment_free(seg); return 1;}
    seg->looseness = (double *) gizmo_gpu_alloc_shared(m.looseness, "ngl_sidx_looseness");
    if(!seg->looseness) {report->bytes_failed = m.looseness; report->failed_memory = SIDX_MEM_SHARED; report->failed_kept = 1; sidx_segment_free(seg); return 1;}
    return 0;
}

/* Build a segment over the members of src into seg.  kStageOnHost: run on the host with the working
 * space in the memory arena, then copy the segment to the device once; otherwise run on the device and
 * build it in place.  maintained: keep what raises need (the slot map, the level schedule, the shear folds,
 * the looseness); a segment that is never raised keeps only what a walk reads.  presorted: the source is
 * already in a spatially coherent order with every member current and inside the key extent
 * (sidx_particles_still_ordered), so the members are tiled in source order with no keys and no sort; the
 * fold still refuses a member no search can bound.  Returns a SidxBuildStatus (report->failure says how a
 * build failed); on anything but SIDX_BUILT the segment is left invalid and nothing the build took is held. */
template <class Exec, bool kStageOnHost, class Source>
static int sidx_build_segment(const Source &src, mode_b_radius_policy_t policy, integertime ti_ref,
                              const struct DriftKickTableView &tables, const struct SidxKeyFrame &frame,
                              int maintained, int presorted, gpu_index_segment_t *seg, struct SidxBuildReport *report)
{
    using Mem = typename std::conditional<kStageOnHost, Kokkos::HostSpace, GIZMO_KOKKOS_DEVICE_SPACE>::type;
    using UV = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
    const double growth = kernel_radius_drift_max_growth_factor(), kernel_floor = All.MinKernelRadius;
    const int num_source = src.count;
    report->refused = 0; report->all_current = 1; report->n_outside = 0; report->n_outside_owned = 0;
    SidxScratch<kStageOnHost> scratch;
    try {
        /* The member count first, so that on the host the index can be staged at its largest possible size
         * BEFORE the build's transients: the arena is a stack, and the transients (above the stage) are then
         * released before the index is allocated on the device. */
        int members_counted = 0;
        Kokkos::parallel_reduce("sidx_build_count", Kokkos::RangePolicy<Exec>(0, num_source),
                                KOKKOS_LAMBDA(int o, int &c) {if(src.is_member(o)) {c++;}}, members_counted);
        /* Sized for the most tiles the members can fill until they are sorted into tiles, exactly after. */
        const struct SidxBuildMemoryPlan bound = sidx_build_memory_plan(num_source, members_counted,
                                                                        sidx_tiles_bound(members_counted), maintained, presorted);
        if(!bound.fits) {throw std::runtime_error("index build: more members than one segment can index");}
        const size_t n_slot_of = bound.slot_of / sizeof(int);
        void *p_tiles = NULL, *p_bvh = NULL, *p_pool = NULL, *p_rows = NULL, *p_slot = NULL, *p_level = NULL, *p_folds = NULL;
        if(kStageOnHost) {
            p_tiles = scratch.take(bound.tiles, "ngl_sidx_stage_tiles");
            p_bvh = scratch.take(bound.bvh, "ngl_sidx_stage_bvh");
            p_pool = scratch.take(bound.pool, "ngl_sidx_stage_pool");
            p_rows = scratch.take(bound.rows, "ngl_sidx_stage_rows");
            if(bound.slot_of > 0) {p_slot = scratch.take(bound.slot_of, "ngl_sidx_stage_slot_of");}
            p_level = scratch.take(bound.level_nodes, "ngl_sidx_stage_level_nodes");
            if(bound.shear_folds > 0) {p_folds = scratch.take(bound.shear_folds, "ngl_sidx_stage_shear_folds");}
        }
        const int transients = scratch.mark();
        /* 1. members, compacted in source order, as ordinals */
        const size_t n_scan = bound.scan / sizeof(int);
        Kokkos::View<int*, Mem, UV> first_member((int *)scratch.take(bound.scan, "ngl_sidx_build_scan"), n_scan);
        Kokkos::parallel_scan("sidx_build_members", Kokkos::RangePolicy<Exec>(0, num_source),
                              KOKKOS_LAMBDA(int o, int &acc, const bool final) {
            if(final) {first_member(o) = acc;}
            if(src.is_member(o)) {acc++;}
            if(final && o == num_source - 1) {first_member(num_source) = acc;}
        });
        int num_members = 0;
        if(num_source > 0) {Kokkos::deep_copy(num_members, Kokkos::subview(first_member, (size_t)num_source));}
        if(num_members != members_counted) {throw std::runtime_error("index build: membership changed between its passes");}
        const size_t n_mem = bound.members / sizeof(int);
        Kokkos::View<int*, Mem, UV> member((int *)scratch.take(bound.members, "ngl_sidx_build_members"), n_mem);
        Kokkos::parallel_for("sidx_build_member_list", Kokkos::RangePolicy<Exec>(0, num_source), KOKKOS_LAMBDA(int o) {
            if(src.is_member(o)) {member(first_member(o)) = o;}
        });
        Exec().fence();
        gizmo_gpu_check_last_error("sidx_build_member_list", num_source);
        if(!presorted) {
            Kokkos::View<Morton128*, Mem, UV> key((Morton128 *)scratch.take(bound.keys, "ngl_sidx_build_keys"), n_mem);
            /* 2. keys: only the position is needed here */
            SidxBuildTally keys_tally;
            Kokkos::parallel_reduce("sidx_build_keys", Kokkos::RangePolicy<Exec>(0, num_members),
                                    KOKKOS_LAMBDA(int s, SidxBuildTally &c) {
                const int o = member(s);
                double center[3], hw = 0.0, image[3], box_hw;
                int current = 0, folds_up = 0, folds_down = 0;
                if(src.position(o, ti_ref, tables, center, &hw, &current)) {
                    c.refused++; key(s).hi = 0; key(s).lo = 0; return;
                }
                sidx_member_image(center, hw, image, &box_hw, &folds_up, &folds_down);
                key(s) = sidx_member_key(image, frame);
                if(key(s).hi & SIDX_KEY_OUTSIDE) {c.outside++; if(src.is_owned(o)) {c.outside_owned++;}}
            }, Kokkos::Sum<SidxBuildTally>(keys_tally));
            Exec().fence();
            gizmo_gpu_check_last_error("sidx_build_keys", num_members);
            report->n_outside = keys_tally.outside; report->n_outside_owned = keys_tally.outside_owned;
            if(keys_tally.refused > 0) {report->refused = 1; sidx_segment_free(seg); return SIDX_REFUSED_MEMBER;}
            /* 3. sort: only the key and the member ordinal move */
            if(num_members > 1) {
                /* The sort allocates its own working space; a refusal of it is a refused allocation like any other
                 * (a typed one, std::bad_alloc), its size known only by the allowance. */
                try {
                    Kokkos::Experimental::sort_by_key(Exec(), Kokkos::subview(key, std::make_pair(0, num_members)),
                                                      Kokkos::subview(member, std::make_pair(0, num_members)), Morton128Less{});
                    Exec().fence();
                } catch(const std::bad_alloc &) {
                    throw SidxAllocationRefused((size_t)num_members * SIDX_SORT_BYTES_PER_MEMBER,
                                                kStageOnHost ? SIDX_MEM_HOST : SIDX_MEM_DEVICE, 1);
                }
            }
        }
        /* 4. tiles: members inside the extent fill tiles from slot 0; the rest start at a fresh tile */
        const int n_inside = num_members - report->n_outside;
        const int tiles_inside = (n_inside + TILE_TARGET_SIZE - 1) / TILE_TARGET_SIZE;
        const int tiles_outside = (report->n_outside + TILE_TARGET_SIZE - 1) / TILE_TARGET_SIZE;
        const int ntiles = (tiles_inside + tiles_outside > 0) ? tiles_inside + tiles_outside : 1;
        const struct SidxBuildMemoryPlan plan = sidx_build_memory_plan(num_source, num_members, ntiles, maintained, presorted);
        if(!plan.fits) {throw std::runtime_error("index build: more members than one segment can index");}
        const int nslots = plan.nslots;
        const size_t n_folds = plan.shear_folds / sizeof(int);
        /* What the index keeps: on the device route allocated now and built in place; on the host route
         * already staged above, and copied once below. */
        if(!kStageOnHost) {
            if(sidx_alloc_kept(seg, maintained, plan, report)) {return SIDX_ALLOCATION_REFUSED;}
            p_tiles = seg->d_tiles; p_bvh = seg->d_bvh; p_pool = seg->d_pool; p_rows = seg->d_compact_xyzh;
            p_slot = seg->d_slot_of; p_folds = seg->d_shear_folds;
            /* the level schedule is kept by a maintained segment, a working array of the build otherwise */
            p_level = maintained ? (void *)seg->d_level_nodes : scratch.take(plan.level_nodes, "ngl_sidx_build_level_nodes");
        }
        Kokkos::View<sfc_tile_t*, Mem, UV> tiles((sfc_tile_t *)p_tiles, (size_t)ntiles);
        Kokkos::View<tile_bvh_node_t*, Mem, UV> bvh((tile_bvh_node_t *)p_bvh, (size_t)plan.nnodes);
        Kokkos::View<int*, Mem, UV> pool((int *)p_pool, (size_t)nslots);
        Kokkos::View<double*, Mem, UV> rows((double *)p_rows, (size_t)nslots * SIDX_ROW_WIDTH);
        Kokkos::View<int*, Mem, UV> slot_of((int *)p_slot, n_slot_of);
        Kokkos::View<int*, Mem, UV> level_nodes((int *)p_level, (size_t)plan.nnodes);
        Kokkos::View<int*, Mem, UV> shear_folds((int *)p_folds, p_folds ? n_folds : 0);
        Kokkos::deep_copy(pool, -1);
        Kokkos::deep_copy(slot_of, -1);
        /* One team per tile.  Each lane writes its member's row and pool slot and contributes its bounds; the
         * team joins them (SidxTileJoin) and one lane writes the tile, so the tile has a single writer.  The
         * counts the build reports are summed over the tiles the same way.  On the device at most one lane per
         * member (a team of TILE_TARGET_SIZE lanes); on the host one lane per team, since the tiles alone keep every thread busy there and a
         * team of host threads would only add a rendezvous per tile. */
        SidxBuildTally fold_tally;
        const auto fold = KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<Exec>::member_type &team, SidxBuildTally &c) {
            const int t = team.league_rank();
            const int k0 = (t < tiles_inside) ? t * TILE_TARGET_SIZE : n_inside + (t - tiles_inside) * TILE_TARGET_SIZE;
            const int k_end = (t < tiles_inside) ? n_inside : num_members;
            const int k1 = (k0 + TILE_TARGET_SIZE < k_end) ? k0 + TILE_TARGET_SIZE : k_end;
            struct SidxTileSummary sum;
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, k0, (k1 > k0) ? k1 : k0), [&](int k, struct SidxTileSummary &b) {
                const int o = member(k), slot = t * TILE_TARGET_SIZE + (k - k0);
                struct SidxMember m;
                if(src.describe(o, policy, ti_ref, tables, growth, kernel_floor, m)) {
                    b.refused++;   /* its position or its motion cannot be bounded */
                    return;
                }
                if(!m.current) {b.not_current++;}
                double *row = &rows((size_t)slot * SIDX_ROW_WIDTH);
                row[0] = m.center[0]; row[1] = m.center[1]; row[2] = m.center[2]; row[3] = m.r; row[4] = m.r_drifted;
                pool(slot) = src.global(o);
                if(slot_of.extent(0) > 0) {slot_of(o) = slot;}
                if(shear_folds.extent(0) > 0) {
                    shear_folds(2 * (size_t)slot) = m.folds_up; shear_folds(2 * (size_t)slot + 1) = m.folds_down;
                }
                sidx_summary_add_member(b, m);
            }, SidxTileJoin(sum));
            Kokkos::single(Kokkos::PerTeam(team), [&]() {
                sfc_tile_t *tile = &tiles(t);
                tile->first = t * TILE_TARGET_SIZE; tile->count = (k1 > k0) ? k1 - k0 : 0; tile->bvh_leaf = -1;
                sidx_box_empty(tile);
                sidx_union_plain(tile, sum);
                tile->hw = sum.hw;
                c.refused += sum.refused; c.not_current += sum.not_current;
            });
        };
        const int fold_team = Kokkos::SpaceAccessibility<Exec, Kokkos::HostSpace>::accessible ? 1 : TILE_TARGET_SIZE;
        Kokkos::parallel_reduce("sidx_build_fold", Kokkos::TeamPolicy<Exec>(ntiles, fold_team),
                                fold, Kokkos::Sum<SidxBuildTally>(fold_tally));
        Exec().fence();
        gizmo_gpu_check_last_error("sidx_build_fold", ntiles);
        /* whether every member was already at ti_ref, and any member the fold could not bound */
        report->all_current = (fold_tally.not_current == 0);
        if(fold_tally.refused > 0) {report->refused = 1; sidx_segment_free(seg); return SIDX_REFUSED_MEMBER;}
        /* 5. the BVH: its shape set out on the host, its bounds level by level from the leaves.  The shape's
         * block is taken last and given back at the end of this block, so the arena stays a stack. */
        SidxBvhShape shape;
        {
            Kokkos::View<int*, Mem, UV> leaf_of_tile((int *)scratch.take(plan.leaf_of_tile, "ngl_sidx_build_leaf_of_tile"), (size_t)ntiles);
            SidxScratch<true> host_side;
            sidx_bvh_shape(shape, plan, host_side.take(plan.bvh_shape, "ngl_sidx_build_bvh_shape"));
            Kokkos::deep_copy(bvh, Kokkos::View<const tile_bvh_node_t*, Kokkos::HostSpace, UV>(shape.nodes, (size_t)shape.nnodes));
            Kokkos::deep_copy(level_nodes, Kokkos::View<const int*, Kokkos::HostSpace, UV>(shape.level_nodes, (size_t)shape.nnodes));
            Kokkos::deep_copy(leaf_of_tile, Kokkos::View<const int*, Kokkos::HostSpace, UV>(shape.leaf_of_tile, (size_t)ntiles));
            Kokkos::parallel_for("sidx_build_tile_leaf", Kokkos::RangePolicy<Exec>(0, ntiles), KOKKOS_LAMBDA(int t) {
                tiles(t).bvh_leaf = leaf_of_tile(t);
            });
            Exec().fence();
        }
        /* One writer per node, and its children were finished at the level before: plain unions. */
        for(int L = 0; L < shape.nlevels; L++) {
            Kokkos::parallel_for("sidx_build_bvh_level",
                                 Kokkos::RangePolicy<Exec>(shape.level_offsets[L], shape.level_offsets[L + 1]),
                                 KOKKOS_LAMBDA(int q) {
                tile_bvh_node_t *node = &bvh(level_nodes(q));
                if(node->left < 0) {sidx_union_plain(node, tiles(-(node->left + 1)));}
                else {sidx_union_plain(node, bvh(node->left)); sidx_union_plain(node, bvh(node->right));}
            });
            Exec().fence();
        }
        gizmo_gpu_check_last_error("sidx_build_bvh_level", shape.nlevels);
        double looseness = 0;
        if(maintained) {
            Kokkos::parallel_reduce("sidx_build_looseness", Kokkos::RangePolicy<Exec>(0, ntiles),
                                    KOKKOS_LAMBDA(int t, double &w) {const double x = sidx_tile_looseness(tiles(t)); if(x > w) {w = x;}},
                                    Kokkos::Max<double>(looseness));
            Exec().fence();
        }
        /* The build's transients go now; on the host route, before the index is allocated on the device and
         * the stage copied to it once (the stage goes when this returns). */
        scratch.release_to(transients);
        if(kStageOnHost) {
            if(sidx_alloc_kept(seg, maintained, plan, report)) {return SIDX_ALLOCATION_REFUSED;}
            Kokkos::deep_copy(Kokkos::View<sfc_tile_t*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_tiles, (size_t)ntiles), tiles);
            Kokkos::deep_copy(Kokkos::View<tile_bvh_node_t*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_bvh, (size_t)shape.nnodes), bvh);
            Kokkos::deep_copy(Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_pool, (size_t)nslots), pool);
            Kokkos::deep_copy(Kokkos::View<double*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_compact_xyzh, (size_t)nslots * SIDX_ROW_WIDTH), rows);
            if(maintained) {
                Kokkos::deep_copy(Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_slot_of, n_slot_of), slot_of);
                Kokkos::deep_copy(Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_level_nodes, (size_t)shape.nnodes), level_nodes);
            }
            if(seg->d_shear_folds) {
                Kokkos::deep_copy(Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, UV>(seg->d_shear_folds, n_folds), shear_folds);
            }
        }
        if(maintained) {
            memcpy(seg->h_level_offsets, shape.level_offsets, (size_t)(shape.nlevels + 1) * sizeof(int));
            *seg->looseness = looseness;
            seg->nlevels = shape.nlevels;
        }
        Kokkos::fence();
        seg->ntiles = ntiles;
        seg->bvh_nnodes = shape.nnodes;
        seg->bvh_root = shape.nnodes - 1;
        seg->num_pool = nslots;
        seg->source_base = src.base;
        seg->source_count = num_source;
    } catch(const SidxAllocationRefused &e) {
        report->bytes_failed = e.bytes; report->failed_memory = e.memory; report->failed_estimated = e.estimated;
        sidx_segment_free(seg);
        return SIDX_ALLOCATION_REFUSED;
    } catch(const std::exception &e) {
        /* A failure inside the sort (whose own working space Kokkos allocates) or a copy: not known to be
         * memory, so not retried elsewhere.  Its own text says what. */
        snprintf(report->failure, sizeof(report->failure), "%s", e.what());
        sidx_segment_free(seg);
        return SIDX_FAILED;
    }
    return SIDX_BUILT;
}

/* Where a build runs: on the host, its working space in the memory arena, the segment then copied to the
 * device once; or on the device, built in place -- the rank's own particles read where they live (the
 * particle storage is shared with the device), the imported ones from compact records. */
enum SidxBuildRoute {SIDX_ROUTE_HOST = 0, SIDX_ROUTE_DEVICE = 1};

/* The shapes of build whose device cost differs: the rank's own members already in order (no keys, no sort),
 * the rank's own members sorted, and the imported members, sorted, from their compact records. */
enum SidxBuildShape {SIDX_SHAPE_OWNED_PRESORTED = 0, SIDX_SHAPE_OWNED_SORTED = 1, SIDX_SHAPE_IMPORTED = 2, SIDX_NUM_SHAPES = 3};

/* The fewest members for which a build of each shape takes the device route.  No build takes it until its
 * builder has been checked on a device; these are then set from the measured crossover of the two routes. */
static constexpr int SIDX_DEVICE_ROUTE_MIN_MEMBERS[SIDX_NUM_SHAPES] = {INT_MAX, INT_MAX, INT_MAX};

/* The share of the device's free memory a build may plan on: ranks sharing a device allocate between the query
 * and the build.  The projection does not count the particle storage, which the device route reads where it lives
 * and whose placement is decided where it is allocated (gpu_particles_arena). */
static constexpr double SIDX_DEVICE_FREE_FRACTION = 0.8;

/* At most how many members a build of src has, without reading a particle: for the gas index over all of the
 * rank's own particles, its gas cells -- the gas block, and any grains promoted to gas since it was last
 * rearranged (an import never changes N_gas) -- else the particles in the source. */
static int sidx_member_bound(const struct SidxParticleSource &src)
{
    if(src.P == P && src.type_bitmask == 1 && src.base == 0 && src.count == src.owned_end) {
        long long gas = N_gas;
#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)
        gas += Grains_promoted;
#endif
        return (gas < src.count) ? (int)gas : src.count;
    }
    return src.count;
}

/* The device route's peak in device memory: while sorting, the working space, the sort's own and the records;
 * after it, the working space, the records and the kept segment. */
static size_t sidx_device_route_peak(const struct SidxBuildMemoryPlan &m, int num_members, int maintained, int presorted,
                                     int imported)
{
    const size_t records = imported ? m.record_rows + m.record_types + m.record_states : 0;
    const size_t working = m.scan + m.members + m.keys + records;
    const size_t sorting = working + (presorted ? 0 : (size_t)num_members * SIDX_SORT_BYTES_PER_MEMBER);
    const size_t building = working + m.leaf_of_tile + (maintained ? 0 : m.level_nodes) + sidx_kept_device_bytes(m, maintained);
    return (sorting > building) ? sorting : building;
}

/* Whether bytes fit within fraction of the device's free memory: 1 when they do, or when that memory is not
 * known (a refused allocation is then the backstop). */
static int sidx_device_has_room(size_t bytes, double fraction)
{
    size_t free_bytes = 0, total_bytes = 0;
    if(!gizmo_gpu_device_memory(&free_bytes, &total_bytes)) {return 1;}
    return (double)bytes <= fraction * (double)free_bytes;
}

/* The route a build takes: the device route when there is a device apart from the host, the source is the
 * particle storage the device shares (the device route reads members where they live, so never an array of a
 * caller's own), the build is of a size at which that route is the faster, its projected peak fits the device's
 * free memory, and the drift tables are on the device; otherwise the host route.  O(1); below the size, no query
 * of the device at all.  An unknown free memory does not refuse the route: a refused allocation then sends the
 * build to the host. */
static int sidx_build_route(const struct SidxParticleSource &src, int maintained, int presorted)
{
    if(std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>::value) {return SIDX_ROUTE_HOST;}
    if(src.P != P) {return SIDX_ROUTE_HOST;}
    const int shape = src.exact ? SIDX_SHAPE_IMPORTED : (presorted ? SIDX_SHAPE_OWNED_PRESORTED : SIDX_SHAPE_OWNED_SORTED);
    const int members = sidx_member_bound(src);
    if(members < SIDX_DEVICE_ROUTE_MIN_MEMBERS[shape]) {return SIDX_ROUTE_HOST;}
    const struct SidxBuildMemoryPlan m = sidx_build_memory_plan(src.count, members, sidx_tiles_bound(members), maintained, presorted);
    if(!m.fits) {return SIDX_ROUTE_HOST;}
    if(!sidx_device_has_room(sidx_device_route_peak(m, members, maintained, presorted, src.exact), SIDX_DEVICE_FREE_FRACTION)) {
        return SIDX_ROUTE_HOST;
    }
    gx_walk_drift_tables_refresh();
    return g_walk_drift_tables_ok ? SIDX_ROUTE_DEVICE : SIDX_ROUTE_HOST;
}

/* Whether the host route has room for a build of src: its peak in the memory arena -- the staged segment, the
 * working space and the BVH shape, held together -- and, where the device's free memory is known, all of that
 * memory still covers the segment it keeps (which the host route places on the device too; no margin here,
 * since a no is final).  A projection only; every request the build makes is checked again when it is made. */
static int sidx_host_route_fits(const struct SidxParticleSource &src, int maintained, int presorted)
{
    const int members = sidx_member_bound(src);
    const struct SidxBuildMemoryPlan m = sidx_build_memory_plan(src.count, members, sidx_tiles_bound(members), maintained, presorted);
    if(!m.fits) {return 0;}
    if(!std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>::value &&
       !sidx_device_has_room(sidx_kept_device_bytes(m, maintained), 1.0)) {return 0;}
    const size_t always[] = {m.tiles, m.bvh, m.pool, m.rows, m.level_nodes, m.scan, m.members, m.leaf_of_tile, m.bvh_shape};
    const size_t when_sized[] = {m.slot_of, m.shear_folds, m.keys};
    size_t bytes = 0; int blocks = 0;
    for(size_t b : always) {bytes += gizmo_mymalloc_rounded_size(b); blocks++;}
    for(size_t b : when_sized) {if(b > 0) {bytes += gizmo_mymalloc_rounded_size(b); blocks++;}}
    return gizmo_alloc_fits_this_rank(bytes, blocks);
}

/* One attempt to build seg over src on one route.  On return nothing it took is held but, on SIDX_BUILT,
 * the segment. */
static int sidx_build_attempt(int route, const struct SidxParticleSource &src, mode_b_radius_policy_t policy,
                              integertime ti_ref, const struct DriftKickTableView &host_tables,
                              const struct SidxKeyFrame &frame, int maintained, int presorted,
                              gpu_index_segment_t *seg, struct SidxBuildReport *report)
{
    *report = SidxBuildReport{};
    if(route == SIDX_ROUTE_HOST) {
        return sidx_build_segment<Kokkos::DefaultHostExecutionSpace, true>(src, policy, ti_ref, host_tables, frame,
                                                                           maintained, presorted, seg, report);
    }
    const struct DriftKickTableView device_tables = g_walk_drift_tables;   /* refreshed when the route was chosen */
    if(!src.exact) {
        return sidx_build_segment<Kokkos::DefaultExecutionSpace, false>(src, policy, ti_ref, device_tables, frame,
                                                                        maintained, presorted, seg, report);
    }
    /* the imported particles (the exact source) */
    try {
        struct SidxGhostRecords records;
        if(records.pack(src, policy, ti_ref, host_tables, report)) {return SIDX_ALLOCATION_REFUSED;}
        return sidx_build_segment<Kokkos::DefaultExecutionSpace, false>(records.source(), policy, ti_ref, device_tables, frame,
                                                                        maintained, presorted, seg, report);
    } catch(const std::exception &e) {
        snprintf(report->failure, sizeof(report->failure), "%s", e.what());
        return SIDX_FAILED;
    }
}

/* Build seg over the members of src, described at the current time, on the route sidx_build_route chooses.
 * A device build refused for memory, and only that, is tried once on the host, when the host route has room
 * and this rank has not already asked to stop.  maintained: the segment will be kept and raised, so it keeps
 * what raises need and registers its source range with the dirty tracker.  Returns SIDX_BUILT; otherwise the
 * segment is left invalid and the controlled stop has been requested here, naming the cause. */
static int sidx_build_segment_now(const struct SidxParticleSource &src, int maintained, int presorted,
                                  mode_b_radius_policy_t radius_policy, const char *caller_label,
                                  gpu_index_segment_t *seg, struct SidxBuildReport *report)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

#if defined(BOX_PERIODIC)
    /* Hard guard on the box lengths the wrap depends on.  The device predicates
     * image through the canonical macros, which read the boxSize_ and boxHalf_
     * lengths from
     * this TU's AllDeviceMirror; if that mirror is unsynced they read zero,
     * every wrapped separation collapses, the per-axis prunes go dead and the
     * BVH silently degenerates into an exhaustive scan.  Fail loud at the call
     * site so a regression in the mirror sync can never quietly cost the index.
     * (Previously this checked the per-axis flags/box_sizes copies; those are
     * gone with that API, but the failure mode belongs to the globals and is
     * unchanged, so the guard now reads them directly.) */
    {
        const double box_len[3] = { boxSize_X, boxSize_Y, boxSize_Z };
        const int    wraps[3]   = { TILE_PERIODIC_X, TILE_PERIODIC_Y, TILE_PERIODIC_Z };
        for(int k = 0; k < 3; k++) {
            if(wraps[k] && !(box_len[k] > 0.0)) {
                printf("spatial index build: axis %d is periodic but its box length is %g "
                       "(caller='%s'). Likely cause: this TU's AllDeviceMirror not synced from "
                       "host All. Confirm gizmo_gpu_sync_all() has run for this timestep.\n",
                       k, box_len[k], caller_label ? caller_label : "?");
                fflush(stdout);
                endrun(913004);
            }
        }
    }
#endif
    const integertime ti_ref = gizmo_host_ti_current();
    const struct DriftKickTableView host_tables = drift_kick_table_view_host();
    /* The extent keys are measured against: the box when periodic (every image lies inside it), else the
     * padded extent of the last domain decomposition. */
    struct SidxKeyFrame frame;
#if defined(BOX_PERIODIC)
    frame.corner[0] = frame.corner[1] = frame.corner[2] = 0.0;
    frame.len = DMAX(boxSize_X, DMAX(boxSize_Y, boxSize_Z));
#else
    for(int k = 0; k < 3; k++) {frame.corner[k] = DomainCorner[k];}
    frame.len = DomainLen;
#endif
    const int route = sidx_build_route(src, maintained, presorted);
    int rc = sidx_build_attempt(route, src, radius_policy, ti_ref, host_tables, frame, maintained, presorted, seg, report);
    /* A device build refused an allocation is tried on the host -- unless what was refused is the segment
     * itself, which every route needs, or the host route has no room, or a stop is already asked. */
    const struct SidxBuildReport first = *report;
    const char *no_retry = NULL;
    int retried = 0;
    if(rc == SIDX_ALLOCATION_REFUSED && route == SIDX_ROUTE_DEVICE) {
        if(first.failed_kept) {no_retry = "what was refused is the index itself, which the host route needs too";}
        else if(!sidx_host_route_fits(src, maintained, presorted)) {no_retry = "the host route has no room (the memory arena, or the device for the index it keeps)";}
        else if(gizmo_controlled_stop_local_reason()) {no_retry = "a stop is already requested";}
        else {
            retried = 1;
            rc = sidx_build_attempt(SIDX_ROUTE_HOST, src, radius_policy, ti_ref, host_tables, frame, maintained, presorted, seg, report);
        }
    }
    if(rc != SIDX_BUILT) {
        char msg[520], earlier[200] = "";
        const char *who = caller_label ? caller_label : "?";
        const char *where = retried ? "on the host route (after the device route)" : (route == SIDX_ROUTE_DEVICE ? "on the device route" : "on the host route");
        if(retried) {snprintf(earlier, sizeof(earlier), "; the device route had been refused %s%.1f MB of %s", first.failed_estimated ? "about " : "",
                              (double)first.bytes_failed / (1024.0 * 1024.0), sidx_memory_name(first.failed_memory));}
        else if(no_retry) {snprintf(earlier, sizeof(earlier), "; not tried on the host route: %s", no_retry);}
        if(rc == SIDX_REFUSED_MEMBER) {
            snprintf(msg, sizeof(msg), "spatial index build (caller '%s') %s: a member's clock is outside [0, now] or its "
                     "position or velocity is not finite, so no search can bound where it is%s; spatial index left unbuilt", who, where, earlier);
            gizmo_request_controlled_stop(7740, msg, __FILE__, __LINE__, __FUNCTION__);
        } else if(rc == SIDX_ALLOCATION_REFUSED) {
            snprintf(msg, sizeof(msg), "spatial index build (caller '%s'): %s was refused %s%.1f MB of %s for %d particles%s; "
                     "spatial index left unbuilt", who, where, report->failed_estimated ? "about " : "", (double)report->bytes_failed / (1024.0 * 1024.0),
                     sidx_memory_name(report->failed_memory), src.count, earlier);
            gizmo_request_controlled_stop(7712, msg, __FILE__, __LINE__, __FUNCTION__);
        } else {
            snprintf(msg, sizeof(msg), "spatial index build (caller '%s'): the build %s failed (%s) for %d particles%s; "
                     "spatial index left unbuilt", who, where, report->failure, src.count, earlier);
            gizmo_request_controlled_stop(7712, msg, __FILE__, __LINE__, __FUNCTION__);
        }
        return rc;
    }
    /* A member of the rank's own outside the extent the domain was built on: the next step decomposes
     * (the drift latch), and meanwhile its own tile bounds it exactly. */
    if(report->n_outside_owned > 0) {DomainExtentOutgrownLocal = 1;}

    seg->ti_ref = ti_ref;
    seg->rows_are_positions = report->all_current;
    seg->rebuild_needed = 0;
    seg->valid = 1;
    /* Register the source range with the dirty tracker. The rows were written from
     * the live P[] under this index's radius policy, and nothing between there and
     * here can mutate it, so the range starts clean: the first refresh would
     * recompute values it already holds. */
    if(maintained) {
        seg->dirty_handle = gpu_dirty_tracker_register(src.base, src.count, 1);
        if(seg->dirty_handle < 0) {seg->rebuild_needed = 1;}   /* untracked, it cannot be kept: rebuilt on its next use */
    }
    return SIDX_BUILT;
}


enum { SIDX_RAISE_MOTION = 1, SIDX_RAISE_REACH = 2, SIDX_REFIT_REACH = 4 };
/* A raise record's type for a particle that left the pool in place: its slot is retired. */
static constexpr int SIDX_RAISE_LEFT_POOL = -2;
/* Raise records staged per round: a bounded buffer beside the segment, whatever the member count. */
static constexpr int SIDX_RAISE_CHUNK = 1 << 16;

/* The ratio of the two ways to widen the nodes above n touched members: walking each member's path, or
 * one sweep of the index.  The walk is taken while it costs at most this times the sweep. */
static constexpr double SIDX_PATH_WALK_FACTOR = 1.0;

/* What one particle contributes to the index.  Computed on the host, which reads P[] and CellP[] where
 * their pages live -- a device kernel reading them faults the pages across -- and applied on the device,
 * which owns the bounds.  j < 0 marks a particle that contributes nothing; type < 0, no reach. */
struct SidxRaiseRecord {
    int j, type;
    double r, r_drifted;
    double u_min[3], u_max[3], rho;
};

/* Raise the kept owned segment of idx over particles list[0..n) (host indices; the first n particles of its
 * source range when list is null).  Only an owned segment is raised: an imported particle's motion and
 * radius do not change while it is imported.
 * MOTION re-reads each member's velocity range, after its velocity changed; REACH rewrites its row's reaches
 * from its current radius and raises the reach bands, after its radius changed.  REFIT (with REACH, over the
 * whole source range: list null) recomputes the reach bands of every tile and node from its members instead
 * of raising them, after every radius changed; the boxes, motion ranges and membership stay as they are.  A
 * particle that left the pool in place (another type now; gpu_sidx_notify_member_lost) has its slot retired,
 * so no walk returns it; one with no mass (no pair kernel takes it) is skipped.
 * Returns 0; 1 when a member's motion is not finite, so nothing can bound it (the run is stopped); 2 when
 * the raise could not be staged.  Either way the segment no longer bounds its members and is marked to be
 * rebuilt by the next list build, which releases it. */
static int sidx_raise_members(gpu_spatial_index_t *idx, const int *list, int n, int what)
{
    gpu_index_segment_t *seg = &idx->owned;
    if(n <= 0 || !seg->valid) {return 0;}
    if((what & SIDX_REFIT_REACH) && (list || n != seg->source_count)) {seg->rebuild_needed = 1; return 2;}   /* a refit takes every member */
    const int source_base = seg->source_base, source_end = seg->source_base + seg->source_count;
    const int type_bitmask = idx->cache_tbm;
    const mode_b_radius_policy_t policy = idx->cache_radius_policy;
    const int refit = (what & SIDX_REFIT_REACH) != 0;
    const double touched_tiles = (n < seg->ntiles) ? (double) n : (double) seg->ntiles;
    const int walk_paths = !refit && (touched_tiles * seg->nlevels <= SIDX_PATH_WALK_FACTOR * (double)(seg->ntiles + seg->bvh_nnodes));
    /* The records are staged in chunks through one host and one device buffer, so a raise or refit over
     * every member costs a bounded amount of memory beside the segment, whatever the member count. */
    const int chunk = (n < SIDX_RAISE_CHUNK) ? n : SIDX_RAISE_CHUNK;
    std::vector<SidxRaiseRecord> rec;
    try {rec.resize((size_t)chunk);} catch(const std::bad_alloc &) {seg->rebuild_needed = 1; return 2;}
    SidxRaiseRecord *d_rec = (SidxRaiseRecord *) ngl_alloc_device((size_t) chunk * sizeof(SidxRaiseRecord), "sidx_raise_records");
    if(!d_rec) {seg->rebuild_needed = 1; return 2;}
    sfc_tile_t *tiles = seg->d_tiles;
    tile_bvh_node_t *bvh = seg->d_bvh;
    const int *slot_of = seg->d_slot_of;
    double *rows = seg->d_compact_xyzh;
    double *looseness = seg->looseness;
    const int *shear_folds = seg->d_shear_folds;
    int *pool = seg->d_pool;
    if(refit) {
        /* every member is about to contribute its reach again, so the bands start from nothing */
        Kokkos::parallel_for("sidx_refit_tiles", seg->ntiles, KOKKOS_LAMBDA(int t) {
            tiles[t].hmax = 0;
            for(int y = 0; y < TILE_NUM_PTYPES; y++) {tiles[t].hmax_by_type[y] = 0;}
        });
        Kokkos::parallel_for("sidx_refit_nodes", seg->bvh_nnodes, KOKKOS_LAMBDA(int q) {
            bvh[q].hmax = 0;
            for(int y = 0; y < TILE_NUM_PTYPES; y++) {bvh[q].hmax_by_type[y] = 0;}
        });
        Kokkos::fence();
    }
    for(int k0 = 0; k0 < n; k0 += chunk) {
        const int nk = (n - k0 < chunk) ? n - k0 : chunk;
        int fault = 0;
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(|:fault)
#endif
        for(int k = 0; k < nk; k++) {
            SidxRaiseRecord &r = rec[(size_t)k];
            const int j = list ? list[k0 + k] : source_base + k0 + k;
            r.j = -1; r.type = -1; r.r = 0; r.r_drifted = 0; r.rho = 0;
            for(int d = 0; d < 3; d++) {r.u_min[d] = MAX_REAL_NUMBER; r.u_max[d] = -MAX_REAL_NUMBER;}
            if(j < source_base || j >= source_end) {continue;}
            const int type = (int)P[j].Type;
            if(type < 0 || type >= TILE_NUM_PTYPES || !((1 << type) & type_bitmask)) {r.j = j; r.type = SIDX_RAISE_LEFT_POOL; continue;}
            if(!sfc_pool_member(&P[j], type_bitmask)) {continue;}
            if((what & SIDX_RAISE_MOTION) && sfc_member_motion_range(j, P, CellP, r.u_min, r.u_max, &r.rho)) {fault |= 1; continue;}
            if(what & SIDX_RAISE_REACH) {
                r.r = nlr_particle_symmetric_radius(P[j], policy);
                r.r_drifted = nlr_particle_symmetric_radius_after_drift(j, P, policy);
                if(r.r_drifted < r.r) {r.r_drifted = r.r;}
                r.type = type;
            }
            r.j = j;
        }
        if(fault) {
            Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_rec);
            seg->rebuild_needed = 1;
            gizmo_request_controlled_stop(7740, "a particle kept in the neighbour index has a velocity that is not finite, "
                                          "so no search can bound where it is", __FILE__, __LINE__, __FUNCTION__);
            return 1;
        }
        {
            Kokkos::View<const SidxRaiseRecord*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> hv(rec.data(), (size_t) nk);
            Kokkos::View<SidxRaiseRecord*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>> dv(d_rec, (size_t) nk);
            Kokkos::deep_copy(dv, hv);
        }
        Kokkos::parallel_for("sidx_raise_members", nk, KOKKOS_LAMBDA(int k) {
            const SidxRaiseRecord &r = d_rec[k];
            if(r.j < 0) {return;}
            const int slot = slot_of[r.j - source_base];
            if(slot < 0) {return;}
            if(r.type == SIDX_RAISE_LEFT_POOL) {pool[slot] = -1; return;}
            struct SidxRaise m;
            for(int d = 0; d < 3; d++) {m.u_min[d] = r.u_min[d]; m.u_max[d] = r.u_max[d];}
            /* the range must hold the member's motion in its tile's image frame as well (sidx_image_velocity_union) */
            if(shear_folds && m.u_min[0] <= m.u_max[0]) {sidx_image_velocity_union(m.u_min, m.u_max, shear_folds[2 * (size_t)slot], shear_folds[2 * (size_t)slot + 1]);}
            m.rho = r.rho; m.hmax = 0;
            for(int t = 0; t < TILE_NUM_PTYPES; t++) {m.hmax_by_type[t] = 0;}
            if(r.type >= 0) {
                /* the row holds the reaches themselves, which may fall; the bands only rise */
                rows[(size_t)slot * SIDX_ROW_WIDTH + 3] = r.r;
                rows[(size_t)slot * SIDX_ROW_WIDTH + 4] = r.r_drifted;
                m.hmax = r.r_drifted; m.hmax_by_type[r.type] = r.r_drifted;
            }
            sfc_tile_t *tile = &tiles[slot / TILE_TARGET_SIZE];
            /* A tile that already covers the member needs nothing above it either: whoever raised it is
             * raising its ancestors too, or they were covering it already. */
            if(!sidx_widen(tile, m)) {return;}
            if(refit) {return;}
            Kokkos::atomic_max(looseness, sidx_tile_looseness(*tile));
            if(!walk_paths) {return;}
            for(int node = tile->bvh_leaf; node >= 0 && sidx_widen(&bvh[node], m); node = bvh[node].parent) {}
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("sidx_raise_members", nk);
    }
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_rec);
    if(!walk_paths) {sidx_widen_all_levels(seg);}
    if(refit) {
        /* a band that fell makes a tile's size floor smaller, so its looseness is taken again */
        double w = 0;
        Kokkos::parallel_reduce("sidx_refit_looseness", seg->ntiles,
                                KOKKOS_LAMBDA(int t, double &x) {const double y = sidx_tile_looseness(tiles[t]); if(y > x) {x = y;}},
                                Kokkos::Max<double>(w));
        *looseness = w;
    }
    return 0;
}

/* The kept gas index follows the velocities of its members: particles idx_host[0..n) (host indices)
 * just had their velocity changed.  Called after the kick of the active set, and by every other writer
 * of a particle's velocity through gizmo_motion_bound_raise.  Costs the list and the paths above it. */
void gpu_step_sidx_raise_motion(const int *idx_host, int n)
{
    gpu_spatial_index_t *idx = &g_step_sidx;
    if(!idx->owned.valid || n <= 0 || !idx_host) {return;}
    GIZMO_GPU_ENSURE_ALL_FRESH();
    sidx_raise_members(idx, idx_host, n, SIDX_RAISE_MOTION);
}

/* Hand back a list with no pairs, after a stop has already been requested.  It
 * requests nothing itself, so the first error stays the one reported.  d_active
 * may be left null: the runner reads that as "a stop is pending, return". */
static void ngl_leave_csr_empty(gpu_neighbor_list_t *gnl, int num_active)
{
    /* The row offsets are what makes the empty list readable: every consumer walks
     * offsets[aa]..offsets[aa+1] unconditionally, and only reaches the neighbour
     * array when total_pairs > 0. So when the build fails before the offsets exist,
     * this allocates them here rather than handing back a null row index. The
     * request is (num_active+1) 8-byte slots -- orders of magnitude below the walk
     * scratchpad and the pair list whose failure brings us here -- so it is served
     * from what those releases just returned. Should even that fail, the consumer
     * reads a null row index and dies where it would have died anyway; there is no
     * smaller allocation left to fall back to. */
    if(!gnl->offsets) {
        gnl->offsets = (int64_t *) ngl_alloc_shared((size_t)(num_active + 1) * sizeof(int64_t),
                                                    "ngl_pairs_offsets");
    }
    if(gnl->offsets) {for(int aa = 0; aa <= num_active; aa++) {gnl->offsets[aa] = 0;}}
    gnl->total_pairs = 0;
}

/* Hand back the empty list under a controlled stop this call requests itself, its message formatted
 * from fmt.  The first stop requested stays the one reported. */
static void ngl_stop_empty(gpu_neighbor_list_t *gnl, int num_active, int code, int line, const char *fmt, ...)
{
    char msg[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg, sizeof(msg), fmt, args);
    va_end(args);
    gizmo_request_controlled_stop(code, msg, __FILE__, line, "gpu_ngb_list_build");
    ngl_leave_csr_empty(gnl, num_active);
}

/* The per-active and per-pair arrays are the largest transients this loop asks for --
 * the walk scratchpad alone is 512 int slots per active particle. When one cannot be
 * had, the node is out of memory at the size this loop needs. Release the transients,
 * leave a VALID EMPTY list (every active with zero neighbours) so that whatever runs
 * between here and the stop draining walks nothing rather than reading a half-built
 * list, and name the buffer that could not be had. */
static void ngl_build_leave_empty(gpu_neighbor_list_t *gnl, int num_active,
                                  double *d_radii, double *d_source_pos,
                                  int *d_scratch, int *d_counts,
                                  const char *what, size_t bytes)
{
    if(d_counts)     {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_counts);}
    if(d_scratch)    {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_scratch);}
    if(d_source_pos) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);}
    if(d_radii)      {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);}
    ngl_leave_csr_empty(gnl, num_active);
    char msg[256];
    snprintf(msg, sizeof(msg),
             "gpu_ngb_list_build: could not allocate %s (%.1f MB) for %d active particles; "
             "neighbour list left empty",
             what, (double) bytes / (1024.0 * 1024.0), num_active);
    gizmo_request_controlled_stop(7711, msg, __FILE__, __LINE__, __FUNCTION__);
}

/* Keep, in each row, exactly the candidates the loop's pair test accepts: a member of the list's types
 * within R_i of query aa for ONEWAY, within max(R_i, j_radius_scale * r_j) for
 * SYMMETRIC, r_j being the member's reach under the loop's radius policy.  Every candidate must already
 * be current.  Candidates from first_exact on (the imported particles) were found by an exact search, so
 * they are kept as they are.  Rows keep their order and their members' order.  The compaction runs in
 * place in the host copy `ngb`: each row is trimmed on its own, then the kept prefixes move down in row
 * order, which never overwrites a row that has not moved yet; the result replaces the front of the device
 * list. */
static void ngl_trim_rows_to_exact(gpu_neighbor_list_t *gnl, std::vector<int> &ngb, int num_active,
                                   const double *q_pos, const double *q_radius, double radius_factor,
                                   int search_mode, int type_bitmask, mode_b_radius_policy_t radius_policy,
                                   double j_radius_scale, const struct particle_data *Pp, int first_exact)
{
    int64_t *off = gnl->offsets;
    std::vector<int64_t> kept((size_t)num_active);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
#endif
    for(int aa = 0; aa < num_active; aa++) {
        const double R = q_radius[aa] * radius_factor;
        int64_t w = off[aa];
        for(int64_t n = off[aa]; n < off[aa + 1]; n++) {
            const int j = ngb[(size_t)n];
            if(j >= first_exact) {ngb[(size_t)w++] = j; continue;}
            const struct particle_data &pj = Pp[j];
            if(!((1 << pj.Type) & type_bitmask)) {continue;}
            const double h_j = (search_mode == NGB_SEARCH_ONEWAY) ? 0.0
                             : nlr_particle_symmetric_radius(pj, radius_policy) * j_radius_scale;
            if(!gx_pair_accept_wrap_and_test(q_pos[aa*3+0] - (double)pj.Pos[0], q_pos[aa*3+1] - (double)pj.Pos[1],
                                             q_pos[aa*3+2] - (double)pj.Pos[2], R, h_j, search_mode)) {continue;}
            ngb[(size_t)w++] = j;
        }
        kept[aa] = w - off[aa];
    }
    int64_t total = 0;
    for(int aa = 0; aa < num_active; aa++) {
        const int64_t start = off[aa];
        off[aa] = total;
        if(total != start && kept[aa] > 0) {memmove(&ngb[(size_t)total], &ngb[(size_t)start], (size_t)kept[aa] * sizeof(int));}
        total += kept[aa];
    }
    off[num_active] = total;
    gnl->total_pairs = total;
    if(total > 0) {
        Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(ngb.data(), (size_t)total);
        Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>> d(gnl->neighbors, (size_t)total);
        Kokkos::deep_copy(d, h);
    }
}

/* What a walk reads of one segment, and the frame it reads it in.  root < 0: the walk skips it. */
struct SidxSegmentView {
    const double *rows;
    const sfc_tile_t *tiles;
    int ntiles;
    const int *pool;
    const tile_bvh_node_t *bvh;
    int root;
    struct sfc_walk_frame frame;
};

static struct SidxSegmentView sidx_segment_view(const gpu_index_segment_t &seg, const struct sfc_walk_frame &frame, int walk)
{
    struct SidxSegmentView v;
    v.rows = seg.d_compact_xyzh; v.tiles = seg.d_tiles; v.ntiles = seg.ntiles; v.pool = seg.d_pool; v.bvh = seg.d_bvh;
    v.root = (walk && seg.valid) ? seg.bvh_root : -1;
    v.frame = frame;
    return v;
}

/* The one search of an index: the owned segment, then the ghost segment, each in its own frame, into one
 * output.  Returns the TRUE number of neighbours; stores at most `capacity` of them, the owned ones first
 * (the ghost search stores after them, into what capacity they left), so a list whose true count exceeds
 * the capacity is re-searched in the same order into a larger output. */
KOKKOS_INLINE_FUNCTION
int sidx_search_segments(const struct SidxSegmentView &own, const struct SidxSegmentView &ghost,
                         const double pos[3], double h, double j_radius_scale, int search_mode,
                         int *store, int capacity)
{
    const int n_own = (own.root < 0) ? 0
                    : search_neighbors_sfc_gpu(own.rows, pos, h, j_radius_scale, own.tiles, own.ntiles, own.pool,
                                               search_mode, own.bvh, own.root, own.frame, store, capacity);
    const int stored = (n_own < capacity) ? n_own : capacity;
    const int n_ghost = (ghost.root < 0) ? 0
                      : search_neighbors_sfc_gpu(ghost.rows, pos, h, j_radius_scale, ghost.tiles, ghost.ntiles, ghost.pool,
                                                 search_mode, ghost.bvh, ghost.root, ghost.frame,
                                                 store + stored, capacity - stored);
    return n_own + n_ghost;
}

void gpu_ngb_list_build(struct particle_data *P_shared, int num_total,
                        int *active_indices_host, int num_active,
                        int search_mode, int type_bitmask,
                        gpu_neighbor_list_t *gnl,
                        gpu_spatial_index_t *cached_idx,
                        double search_radius_factor,
                        const double *search_radii_host,
                        const double *source_positions_host,
                        const char *caller_label,
                        double j_kernel_radius_scale,
                        mode_b_radius_policy_t radius_policy)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    gnl->num_active = num_active;
    /* Everything this function hands back starts null, so a failure part-way through
     * leaves a list whose unbuilt parts are recognisable rather than whatever the
     * caller's struct happened to contain. Several callers declare it uninitialised. */
    gnl->d_active = NULL; gnl->offsets = NULL; gnl->neighbors = NULL; gnl->total_pairs = 0;
    double t_entry = my_second(); /* DIAG: entry */
    const double cpu_rows_child0 = CPU_ChildCharged;
    /* A cached index holds only the types it was built for. Walked under a
     * different type_bitmask it would return neighbours outside the caller's
     * mask (e.g. DM in a gas-only density walk, whose lazy drift then stops on
     * Type=1). Refuse it here with a clear message, so a Spec routed to the
     * wrong cache fails at the right layer: request the stop, hand back the
     * empty list and return without walking. The companion warning in
     * gpu_neighbor_list.h on gpu_step_sidx_alltypes_ptr documents this
     * invariant. */
    if (cached_idx && (cached_idx->owned.valid || cached_idx->ghost.valid) &&
        cached_idx->cache_tbm != type_bitmask) {
        fprintf(stderr,
            "gpu_ngb_list_build FATAL: caller='%s' type_bitmask=0x%x but cached "
            "SIDX was built with tbm=0x%x. The cached compact_xyzh / pool only "
            "contains the originally-built types — walking it under a different "
            "mask returns wrong-type neighbors and triggers downstream aborts "
            "(e.g. drift_particle 'no prediction into past allowed' on Type=1 "
            "during a gas-only density walk). Spec author: declare "
            "sidx_cache_kind matching your neighbor_type_mask, or rebuild "
            "the cache. See gpu_neighbor_list.h docstring on "
            "gpu_step_sidx_alltypes_ptr.\n",
            caller_label ? caller_label : "?", type_bitmask, cached_idx->cache_tbm);
        fflush(stderr);
        endrun(913005);
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }
    /* The same refusal for a radius_policy mismatch: a cached row's reach is the
     * build-time policy's, so walking it under another policy would use the
     * wrong leaf-side reach and miss valid pairs. */
    if (cached_idx && (cached_idx->owned.valid || cached_idx->ghost.valid) &&
        cached_idx->cache_radius_policy != radius_policy) {
        fprintf(stderr,
            "gpu_ngb_list_build FATAL: caller='%s' radius_policy=0x%x but cached "
            "SIDX was built with radius_policy=0x%x. The cached compact_xyzh[j*4+3] "
            "encodes the build-time policy's per-particle pair-search reach; "
            "walking it under a different policy returns wrong leaf-h_j values "
            "and silently misses valid pairs.  Spec author: ensure Spec::radius_policy "
            "matches the cache's owner (or declare SidxCacheKind::None to rebuild "
            "per call).\n",
            caller_label ? caller_label : "?",
            (unsigned)radius_policy, (unsigned)cached_idx->cache_radius_policy);
        fflush(stderr);
        endrun(913006);
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }
    /* Early-out: with no active particles there is nothing to search.
     * Skip the SIDX build/refresh AND all kernel launches.  Allocate 1-element
     * stubs so the caller's gpu_ngb_list_free path is well-defined (it always
     * frees neighbors/offsets/d_active). */
    if(num_active == 0) {
        gnl->d_active  = (int *) ngl_alloc_shared(sizeof(int), "ngl_pairs_active_stub");
        gnl->offsets   = (int64_t *) ngl_alloc_shared(sizeof(int64_t), "ngl_pairs_offsets_stub");
        gnl->neighbors = (int *) ngl_alloc_device(sizeof(int), "ngl_pairs_neighbors_stub");
        /* Single elements: failing these means the node has no memory left at all.
         * The empty list is the same one every other exhausted allocation here
         * hands back, so the callers need no separate case. */
        if(!gnl->d_active || !gnl->offsets || !gnl->neighbors) {
            ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL,
                                  "the empty-list placeholders", sizeof(int64_t));
            cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
            return;
        }
        gnl->offsets[0] = 0;
        gnl->total_pairs = 0;
        cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
        return;
    }

    /* Use cached spatial index if available, otherwise build fresh.
     * If caller provided a cached_idx but it's not yet built, populate it (this
     * enables persistent caching across calls — caller controls invalidation). */
    /* Value-initialise: the positional form listed fewer entries than the struct has members and
     * put them in the wrong slots, so it overrode the members that carry a non-zero default --
     * dirty_handle (-1 = not registered), cache_tbm (-1 = no build yet) and cache_radius_policy.
     * A zero dirty_handle passes the `>= 0` test that guards gpu_dirty_tracker_unregister, so an
     * index that had never registered anything would unregister handle 0, which belongs to
     * whoever did register first. */
    gpu_spatial_index_t local_idx = {};
    /* An index built for this call alone is released on every way out of it; the
     * list it produced holds no part of it. Freeing an index that was never built
     * (the cached route) does nothing. */
    struct ngl_call_index_release {
        gpu_spatial_index_t *owned;
        ~ngl_call_index_release() {gpu_spatial_index_free(owned);}
    } local_idx_release = {&local_idx};
    gpu_spatial_index_t *idx = cached_idx ? cached_idx : &local_idx;
    const int maintained = (cached_idx != NULL);   /* a kept index is raised; one built for this call is not */
    const integertime t_now = gizmo_host_ti_current();
    const struct DriftKickTableView host_tables = drift_kick_table_view_host();
    /* The particles this call searches: the rank's own, [0, owned_end), and, while an import is live, the
     * imported ones, [ghost_base, ghost_end).  A caller takes the imported particles by passing the whole
     * count, or leaves them out by passing the count of its own (hii_fb searches only the rank's own gas,
     * and the imported segment it leaves out stays as it is for the next caller).  ghost_base cannot move
     * during this call: every layout change refuses a live import.  Only the particle storage holds
     * imported particles: a caller searching an array of its own (group finding) has none. */
    const int pool_live = ghost_pool_is_live() && (P_shared == P);
    const int ghost_base = pool_live ? ghost_get_num_local() : num_total;
    const int ghost_end = pool_live ? ghost_base + ghost_get_num_ghosts() : num_total;
    const int owned_end = (num_total < ghost_base) ? num_total : ghost_base;
    const int walk_ghosts = pool_live && (num_total > ghost_base);
    if(walk_ghosts && num_total != ghost_end) {
        ngl_stop_empty(gnl, num_active, 7741, __LINE__, "gpu_ngb_list_build (caller '%s'): asked for %d of the %d imported "
                       "particles; a neighbour list covers all of them or none", caller_label ? caller_label : "?",
                       num_total - ghost_base, ghost_end - ghost_base);
        return;
    }
    /* Whether the members are already at the time of this search.  The rank's own
     * particles by the full-drift certificate (move_particles drifts only the active
     * set, so it does not advance it, which is what makes it a proof rather than a
     * convention); the imported ones by their owners having advanced them before
     * packing.  A new timestep advances All.Ti_Current, so a certificate from an
     * earlier time simply stops matching. */
    const int owned_current = (gizmo_full_drift_ti() == t_now);
    const int ghosts_current = (ghost_pool_current_ti() == t_now);
    const int pool_current = owned_current && (!walk_ghosts || ghosts_current);
    struct SidxBuildReport report;

    /* The owned segment.  Released unless it still describes the rank's own particles: the
     * same range, and the owned epoch (a change of membership or a position written outside
     * a drift, which no count can show).  A drift changes neither: the segment is read at the
     * time of the search (sfc_tiles.h), so reuse across a drift is the common path, and an
     * import or its cleanup does not touch it. */
    gpu_index_segment_t *own = &idx->owned;
    /* The all-types index holds sinks, stars and dark matter, whose positions and membership can be
     * written in place between two of its searches (sink repositioning, accretion, mergers) without
     * advancing the owned epoch; its owned segment is therefore rebuilt at every import, as the whole
     * index was before it was split. */
    const int owned_follows_import = (idx == &g_step_sidx_alltypes);
    if(own->valid && (own->source_count != owned_end || own->owned_epoch_when_built != g_sidx_owned_epoch ||
                      own->ti_ref > t_now || own->rebuild_needed ||
                      (owned_follows_import && (own->ghost_provenance_when_built != ghost_provenance_epoch() ||
                                                own->member_loss_epoch_when_built != g_sidx_member_loss_epoch)))) {
        sidx_segment_free(own);
    }
    /* A kept segment is rebuilt instead when its boxes can have grown by more than their
     * own size, or when every member is current and the search is large enough to repay a
     * fresh segment, which it reads exactly.  A small search reads the kept one in its
     * drifted frame and the trim makes its list exact, for far less than a build. */
    if(own->valid && own->ti_ref < t_now) {
        const double D_kept = get_drift_factor_impl(own->ti_ref, t_now, 1.0, &host_tables);
        const int fresh_repays = owned_current && (double)num_active >= SIDX_FRESH_SEARCH_FRACTION * (double)owned_end;
        if(fresh_repays || !(*own->looseness * D_kept <= SIDX_MAX_BOX_GROWTH)) {
            sidx_segment_free(own);
        }
    }
    /* Bring a kept segment's reaches up to date.  Every particle whose radius changed
     * since was marked dirty (gizmo_mark_kernel_radius_dirty_*); its row takes its
     * current reach and the bands above it are raised to cover it.  All-dirty refits
     * every reach and band in place: the boxes, motion ranges and membership do not
     * depend on the radii.  A refresh that cannot be staged rebuilds the segment
     * instead: a fresh segment starts current.  Per-segment state means consuming-and-
     * clearing this segment's bits leaves the other registered segments' bitsets untouched. */
    if(own->valid && own->dirty_handle >= 0) {
        const int handle = own->dirty_handle;
        if(2 * (long)gpu_dirty_tracker_popcount(handle) >= (long)own->source_count && own->source_count > 0) {
            /* most members changed: one refit, whose bands may also fall, rather than a raise per member */
            if(sidx_raise_members(idx, NULL, own->source_count, SIDX_RAISE_REACH | SIDX_REFIT_REACH) == 0) {
                gpu_dirty_tracker_clear(handle);
                /* a fallen band shrinks a tile's size floor, so the box growth is judged again */
                if(own->ti_ref < t_now && !(*own->looseness * get_drift_factor_impl(own->ti_ref, t_now, 1.0, &host_tables) <= SIDX_MAX_BOX_GROWTH)) {
                    sidx_segment_free(own);
                }
            } else {
                sidx_segment_free(own);
            }
        } else if(gpu_dirty_tracker_popcount(handle) > 0) {
            /* Drain bitset -> host list -> raise. */
            std::vector<int> dirty_host;
            dirty_host.reserve(gpu_dirty_tracker_popcount(handle));
            gpu_dirty_tracker_consume(handle,
                [](int j, void *ud){ ((std::vector<int> *)ud)->push_back(j); },
                &dirty_host);
            if(sidx_raise_members(idx, dirty_host.data(), (int)dirty_host.size(), SIDX_RAISE_REACH)) {
                sidx_segment_free(own);
            }
        }
    }
    if(!own->valid) {
        /* A build that is refused an allocation, or meets a member no search can bound, leaves
         * the segment invalid and has already asked for the stop, naming the cause.
         * There is nothing to walk, so hand back the same empty list any other
         * exhausted allocation here produces. */
        const struct SidxParticleSource src = {P_shared, CellP, 0, owned_end, owned_end, type_bitmask, 0};
        /* Only the gas index takes the order as it stands: every in-place write of a gas particle's
         * position, and every particle that becomes gas, advances the owned epoch, which is not so for the
         * other types.  (A particle that stops being gas leaves its slot to be retired, not the order.) */
        const int presorted = (P_shared == P) && (type_bitmask == 1) && sidx_particles_still_ordered(owned_end, t_now);
        if(sidx_build_segment_now(src, maintained, presorted, radius_policy, caller_label, own, &report)) {
            ngl_leave_csr_empty(gnl, num_active);
            return;
        }
        own->owned_epoch_when_built = g_sidx_owned_epoch;
        own->member_loss_epoch_when_built = g_sidx_member_loss_epoch;
        own->ghost_provenance_when_built = ghost_provenance_epoch();
        idx->cache_tbm = type_bitmask;
        idx->cache_radius_policy = radius_policy;
    }

    /* The ghost segment: the imported particles, exact.  An imported particle is current
     * when it arrives and stays so for as long as it is imported (Ti_Current advances only
     * between sync points, and its owner writes its values back, not its position, radius
     * or type), so the segment needs no raises and is released with the import
     * (gpu_sidx_ghost_pool_cleanup).  Rebuilt when the import it was built over is not
     * the live one, or the time is not the one it was built at. */
    gpu_index_segment_t *ghost = &idx->ghost;
    if(ghost->valid && (!pool_live || ghost->ghost_provenance_when_built != ghost_provenance_epoch() ||
                        ghost->source_base != ghost_base || ghost->source_count != ghost_end - ghost_base ||
                        ghost->ti_ref != t_now)) {
        sidx_segment_free(ghost);
    }
    /* Searched exactly only while its owners' certificate holds; otherwise there is no exact
     * search of the imported particles, kept or fresh. */
    if(walk_ghosts && !ghosts_current) {
        ngl_stop_empty(gnl, num_active, 7742, __LINE__, "gpu_ngb_list_build (caller '%s'): the imported particles are "
                       "not certified to be at the time of this search, so no exact search of them exists; neighbour "
                       "list left empty", caller_label ? caller_label : "?");
        return;
    }
    if(walk_ghosts && !ghost->valid) {
        /* Built aside and published only whole: a build that fails leaves no part of a
         * segment visible, and the owned segment untouched. */
        gpu_index_segment_t fresh;
        const struct SidxParticleSource src = {P_shared, CellP, ghost_base, ghost_end - ghost_base, ghost_base, type_bitmask, 1};
        if(sidx_build_segment_now(src, 0, 0, radius_policy, caller_label, &fresh, &report)) {
            ngl_leave_csr_empty(gnl, num_active);
            return;
        }
        if(!fresh.rows_are_positions) {
            sidx_segment_free(&fresh);
            ngl_stop_empty(gnl, num_active, 7742, __LINE__, "gpu_ngb_list_build (caller '%s'): an imported particle is "
                           "not at the time of this search although its import was certified current; neighbour list "
                           "left empty", caller_label ? caller_label : "?");
            return;
        }
        fresh.ghost_provenance_when_built = ghost_provenance_epoch();
        *ghost = fresh;
        idx->cache_tbm = type_bitmask;
        idx->cache_radius_policy = radius_policy;
    }

    /* How this search reads each segment (sfc_walk_frame).  The owned segment: exact when
     * nothing can have moved since its rows were written -- built at this time from members
     * that were all current then, and all current now; each member's reach is its stored
     * one when they are current (the rows were brought up to date above), otherwise its
     * drifted reach.  The ghost segment is always read exactly. */
    struct sfc_walk_frame own_frame;
    own_frame.D = (own->ti_ref < t_now) ? get_drift_factor_impl(own->ti_ref, t_now, 1.0, &host_tables) : 0.0;
    own_frame.exact = (own->ti_ref == t_now) && owned_current && own->rows_are_positions;
    own_frame.reach_current = owned_current;
    struct sfc_walk_frame ghost_frame;
    ghost_frame.D = 0.0;
    ghost_frame.exact = 1;
    ghost_frame.reach_current = 1;
    if(!(own_frame.D >= 0.0 && own_frame.D < 1.0e30)) {
        ngl_stop_empty(gnl, num_active, 7740, __LINE__, "gpu_ngb_list_build (caller '%s'): the drift interval since the "
                       "neighbour index was built is %g, which bounds nothing; neighbour list left empty",
                       caller_label ? caller_label : "?", own_frame.D);
        return;
    }
    const struct SidxSegmentView own_view = sidx_segment_view(*own, own_frame, 1);
    const struct SidxSegmentView ghost_view = sidx_segment_view(*ghost, ghost_frame, walk_ghosts);

    /* Active indices: always re-uploaded (changes per call) */
    size_t active_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(int);
    gnl->d_active = (int *) ngl_alloc_shared(active_bytes, "ngl_pairs_active");
    if(!gnl->d_active) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the active-index list", active_bytes); return;}
    memcpy(gnl->d_active, active_indices_host, num_active * sizeof(int));

    /* Each query's radius and position.  A caller may supply either (a loop with a
     * different kernel than the particle's own, e.g. KernelRadiusDM or AGS_Hsml; a
     * source not backed by a P[] entry, e.g. a grid cell).  Otherwise they are the
     * query particle's own -- its reach under the loop's policy, its stored position,
     * which must be at the time of the search (checked below) -- and never the index's
     * rows, which describe members at the index's reference time.  Positions:
     * [aa*3 + k] for axis k. */
    size_t radii_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(double);
    double *d_radii = (double *) ngl_alloc_shared(radii_bytes, "ngl_pairs_radii");
    if(!d_radii) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the per-active search radii", radii_bytes); return;}
    size_t srcpos_bytes = (size_t)((num_active > 0) ? num_active : 1) * 3 * sizeof(double);
    double *d_source_pos = (double *) ngl_alloc_shared(srcpos_bytes, "ngl_pairs_source_pos");
    if(!d_source_pos) {ngl_build_leave_empty(gnl, num_active, d_radii, NULL, NULL, NULL, "the per-active source positions", srcpos_bytes); return;}
    /* A query named by its particle index is searched from where that particle is, so it must be at the time of
     * the search: the import that brought its neighbours, the walk on the other route and the pair kernels all
     * read its stored position, and a prediction here would put the list in a frame of its own.  A source at
     * any other position passes it (source_positions_host). */
    long queries_not_current = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:queries_not_current)
#endif
    for(int aa = 0; aa < num_active; aa++) {
        const int i = active_indices_host[aa];
        d_radii[aa] = search_radii_host ? search_radii_host[aa] : nlr_particle_symmetric_radius(P_shared[i], radius_policy);
        if(source_positions_host) {
            for(int k = 0; k < 3; k++) {d_source_pos[aa*3 + k] = source_positions_host[aa*3 + k];}
        } else {
            if(P_shared[i].Ti_current != t_now) {queries_not_current++;}
            for(int k = 0; k < 3; k++) {d_source_pos[aa*3 + k] = (double)P_shared[i].Pos[k];}
        }
    }
    if(queries_not_current > 0) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
        ngl_stop_empty(gnl, num_active, 7743, __LINE__, "gpu_ngb_list_build (caller '%s'): %ld of %d queries named by "
                       "particle index are not at the time of the search, so they have no position the rest of the step "
                       "agrees on; neighbour list left empty", caller_label ? caller_label : "?", queries_not_current, num_active);
        return;
    }

    /* Allocate CSR offsets (64-bit row pointers) */
    size_t offsets_bytes = (size_t)(num_active + 1) * sizeof(int64_t);
    gnl->offsets = (int64_t *) ngl_alloc_shared(offsets_bytes, "ngl_pairs_offsets");
    if(!gnl->offsets) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, NULL, NULL, "the CSR row offsets", offsets_bytes); return;}

    /* Per-particle scratchpad for fused single-pass build. Each active particle
     * gets a fixed stride (NGL_SCRATCH_STRIDE) of int slots in d_scratch; the BVH
     * walk emits j-indices directly there, with the count tracked in d_counts.
     * After scan + compact we transcribe into the dense CSR neighbors[] array.
     * Memory: stride * num_active * 4 bytes (e.g. 256 * 2M * 4 = 2GB). */
    constexpr int NGL_SCRATCH_STRIDE = 512;
    /* size_t cast required: int * int overflows for num_active > ~4.19M (e.g. fire_m11i
     * gas-per-rank), wrapping to negative int → ~UINT64_MAX after promotion to size_t. */
    size_t na_safe = (size_t)((num_active > 0) ? num_active : 1);
    size_t scratch_bytes = na_safe * (size_t)NGL_SCRATCH_STRIDE * sizeof(int);
    int *d_scratch = (int *) ngl_alloc_device(scratch_bytes, "ngl_pairs_scratch");
    if(!d_scratch) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, NULL, NULL, "the neighbour-walk scratchpad", scratch_bytes); return;}
    size_t counts_bytes = na_safe * sizeof(int);
    int *d_counts  = (int *) ngl_alloc_device(counts_bytes, "ngl_pairs_counts");
    if(!d_counts) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, d_scratch, NULL, "the per-active neighbour counts", counts_bytes); return;}

    /* Drain any prior async GPU work before the passes below. */
    Kokkos::fence();
    /* Fused single pass: BVH walk + write neighbors into per-particle scratchpad */
    {
        int *scratch = d_scratch;
        int *counts = d_counts;
        int smode = search_mode;

        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const struct SidxSegmentView own_walk = own_view, ghost_walk = ghost_view;
        Kokkos::parallel_for("ngb_fused", num_active, KOKKOS_LAMBDA(int aa) {
            double h_i = radii[aa] * sr_fac;
            double pos_i[3] = {src_pos[aa*3+0], src_pos[aa*3+1], src_pos[aa*3+2]};
            int cnt = sidx_search_segments(own_walk, ghost_walk, pos_i, h_i, j_rad_scale, smode,
                                           &scratch[(size_t)aa * NGL_SCRATCH_STRIDE], NGL_SCRATCH_STRIDE);
            counts[aa] = cnt;
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("ngb_fused", num_active);

    }

    /* Count overflow particles (count > stride: they need a re-walk in compact phase) */
    int overflow_count = 0;
    {
        int *counts = d_counts;
        Kokkos::parallel_reduce("ngb_overflow_check", num_active,
            KOKKOS_LAMBDA(int aa, int &local) {
                if(counts[aa] > NGL_SCRATCH_STRIDE) local++;
            }, overflow_count);
        Kokkos::fence();
    }

    /* GPU exclusive prefix scan: counts → offsets, returning total.
       Counts are correct even for overflow particles (search_neighbors_sfc_gpu
       returns the true count regardless of bounded write). */
    int64_t total_ll = 0;
    {
        int *counts = d_counts;
        int64_t *offsets = gnl->offsets;
        Kokkos::parallel_scan("ngb_offsets_scan", num_active,
            KOKKOS_LAMBDA(int aa, int64_t &update, const bool final) {
                int64_t v = (int64_t)counts[aa];
                if(final) offsets[aa] = update;
                update += v;
            }, total_ll);
        Kokkos::fence();
    }
    /* Sanity guard: the CSR index (gnl->total_pairs / gnl->offsets) is 64-bit,
     * so INT_MAX is no longer a ceiling. Still abort cleanly on a nonsensical
     * total_pairs (negative => prefix-scan corruption; or a count so large the
     * neighbors allocation byte size would overflow size_t) rather than letting
     * a bad value reach kokkos_malloc / the device compact pass. */
    if(total_ll < 0 || (uint64_t)total_ll > (uint64_t)(SIZE_MAX / sizeof(int))) {
        fprintf(stderr,
            "[NGL FATAL rank=%d] CSR total_pairs implausible: caller=%s num_active=%d "
            "total_pairs=%lld  search_radius_factor=%g j_kernel_radius_scale=%g\n",
            ThisTask, caller_label ? caller_label : "?", num_active, (long long)total_ll,
            search_radius_factor, j_kernel_radius_scale);
        fflush(stderr);
        endrun(915100);
    }
    int64_t total = total_ll;
    gnl->offsets[num_active] = total;
    gnl->total_pairs = total;

    /* Allocate CSR neighbors array (length is 64-bit; element type stays int) */
    size_t neighbors_bytes = (size_t)((total > 0) ? total : 1) * sizeof(int);
    gnl->neighbors = (int *) ngl_alloc_device(neighbors_bytes, "ngl_pairs_neighbors");
    if(!gnl->neighbors) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, d_scratch, d_counts, "the CSR neighbour list", neighbors_bytes); return;}

    /* Compact: copy from per-particle scratchpad into dense CSR neighbors[]. */
    {
        int *scratch = d_scratch;
        int *counts = d_counts;
        int64_t *offsets = gnl->offsets;
        int *neighbors = gnl->neighbors;
        int smode = search_mode;
        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const struct SidxSegmentView own_walk = own_view, ghost_walk = ghost_view;
        Kokkos::parallel_for("ngb_compact", num_active, KOKKOS_LAMBDA(int aa) {
            int n = counts[aa];
            int64_t dst = offsets[aa];
            if(n <= NGL_SCRATCH_STRIDE) {
                size_t src = (size_t)aa * NGL_SCRATCH_STRIDE;
                for(int k = 0; k < n; k++) neighbors[dst + k] = scratch[src + k];
            } else {
                /* Overflow path: re-walk BVH writing directly into neighbors[] */
                double h_i = radii[aa] * sr_fac;
                double pos_i[3] = {src_pos[aa*3+0], src_pos[aa*3+1], src_pos[aa*3+2]};
                sidx_search_segments(own_walk, ghost_walk, pos_i, h_i, j_rad_scale, smode, &neighbors[dst], n);
            }
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("ngb_compact", num_active);
    }

    /* Free temporaries.  The staged query radii and positions are kept for the trim below. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_scratch);
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_counts);


    /* Two things may remain to be done to the walk's list, both on one host copy of it.
     *
     * The members it found that are behind the time of this search are drifted to it
     * (Attack C: lazy drift).  A kernel reads each neighbour's predicted state --
     * CellP[j].VelPred / Density / InternalEnergyPred / KernelRadius -- at this time,
     * and drift_particle(j, time1) provides all of it.  Not needed when the pool is
     * already current (pool_current above): every drift would return at once.
     * MEASURED on a production run: on fulldrift steps this sweep was 92.0 billion
     * pool visits with zero members behind, at ~97 s per rank.
     *
     * Then, unless the walk read the index exactly, its list is a superset -- every
     * member the pair test can accept once drifted -- and each row is trimmed to the
     * members that test accepts at their current positions and reaches. */
    const integertime time1 = t_now;
    const int need_drift = !pool_current;
    const int need_trim = !own_frame.exact;
    if(gnl->total_pairs > 0 && gnl->neighbors && (need_drift || need_trim)) {
        std::vector<int> ngb_host((size_t)gnl->total_pairs);
        gpu_ngb_copy_neighbors_to_host(gnl, ngb_host.data());
        if(need_drift) {
            /* Ghosts imported for this step were advanced to the current time by
             * their owners before being packed, so the whole imported segment is
             * already current and there is nothing to confirm per ghost. The pool's
             * stamp is compared against the time THIS call needs rather than trusted
             * on its own, so a pool carried over from an earlier time still gets
             * checked particle by particle.  The imported particles start at
             * ghost_base whatever count the caller passed. */
            /* Collect the distinct members that are behind, then advance them in one
             * threaded pass, rather than calling drift_particle once per visit.
             *
             * The drift is real per-particle work -- it runs the implicit
             * thermochemistry solve through set_eos_pressure -- so it must not be
             * strictly serial. A member appears once per PAIR, and two threads
             * testing one particle's Ti_current before either writes would advance
             * it twice, so the distinct set is established first. The stamp is
             * generation-counted and never needs clearing between calls.
             *
             * MEASURED on a production run: this leaves the loop's own rank skew at
             * a tenth of what the per-visit form generated, and that skew was being
             * absorbed by the convergence barrier downstream. */
            static std::vector<unsigned int> pool_seen;
            static unsigned int pool_seen_gen = 0;
            static std::vector<int> pool_behind;
            if((int)pool_seen.size() < num_total) {pool_seen.assign((size_t)num_total, 0u);}
            if(++pool_seen_gen == 0u) {std::fill(pool_seen.begin(), pool_seen.end(), 0u); pool_seen_gen = 1u;}
            pool_behind.clear();
            for(int64_t idx_n = 0; idx_n < gnl->total_pairs; idx_n++) {
                int j = ngb_host[idx_n];
                if(j < 0 || j >= num_total) continue;
                if(ghosts_current && j >= ghost_base) continue;
                if(pool_seen[(size_t)j] == pool_seen_gen) continue;
                pool_seen[(size_t)j] = pool_seen_gen;
                if(P[j].Ti_current != time1) {pool_behind.push_back(j);}
            }
            {
                const int n_behind = (int)pool_behind.size();
                const int *behind_idx = pool_behind.data();
                drift_particles_batch(behind_idx, n_behind, time1);
            }
            /* The drift moved the KernelRadius of exactly the members it advanced, so
             * those are the radii to mark dirty; a member it left alone kept its own. */
            gizmo_mark_kernel_radius_dirty_indices(pool_behind.data(), (int)pool_behind.size());
            /* Move detector baseline past the lazy drift's Ti_current/Pos updates
             * — those are predicted-state setup, not kernel writes that need
             * writeback. Subsequent kernel-side writes to ghost particles will
             * still be flagged by ghost_write_detector_end(). No-op when
             * GIZMO_GPU_ARENA_DEBUG is undefined or detector is inactive. */
            ghost_write_detector_resnapshot_after_lazy_drift();
        }
        if(need_trim) {
            ngl_trim_rows_to_exact(gnl, ngb_host, num_active, d_source_pos, d_radii, search_radius_factor,
                                   search_mode, type_bitmask, radius_policy, j_kernel_radius_scale, P_shared,
                                   walk_ghosts ? ghost_base : num_total);
        }
    }
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);

    /* Release this call's own index now, so its fence and frees are charged to this
     * build (the release at scope exit then finds nothing to do). */
    gpu_spatial_index_free(&local_idx);

    /* Charge list-build wall, less any kernel time already charged inside it,
     * so the two rows never overlap. */
    cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
}


void gpu_ngb_list_free(gpu_neighbor_list_t *gnl)
{
    if(gnl->neighbors) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->neighbors);
    if(gnl->offsets)   Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->offsets);
    if(gnl->d_active)  Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->d_active);
    gnl->neighbors = NULL;
    gnl->offsets = NULL;
    gnl->d_active = NULL;
    gnl->num_active = 0;
    gnl->total_pairs = 0;
}

/* Copy gnl->neighbors (DEVICE_SPACE) into a host buffer. Caller owns host_dest.
   For host-side per-source loops (radfb_local, merge_split, density.cc:868
   HYDRO_VOLUME_CORRECTIONS path, turb_powerspectra, twopoint) that index
   neighbors[] from CPU code. */
void gpu_ngb_copy_neighbors_to_host(const gpu_neighbor_list_t *gnl, int *host_dest)
{
    if(!host_dest || !gnl || gnl->total_pairs <= 0 || !gnl->neighbors) {return;}
    Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        h(host_dest, (size_t)gnl->total_pairs);
    Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        d(gnl->neighbors, (size_t)gnl->total_pairs);
    Kokkos::deep_copy(h, d);
}


/* High-level wrapper: build a symmetric neighbor list on GPU and return it
   in the mymalloc-based neighbor_list_t format expected by gradient/hydro.
   Called from accel.cc (which is NOT compiled by nvcc).
   search_radius_factor: multiplier on KernelRadius (default 1.0; >1 for TURB_DIFF_DYNAMIC). */
void gpu_build_symmetric_neighbor_list(struct particle_data *P_host, int num_total,
                                       int *active_indices, int num_active,
                                       neighbor_list_t *out,
                                       double search_radius_factor)
{
    /* Use the per-step particle arena instead of a dedicated full-NumPart memcpy.
     * The arena's fast path is a no-op when valid (e.g. when gradient/hydro
     * already populated it earlier in the step), avoiding ~2.3s of redundant
     * P-copy on small-N symlist invocations.  Pass the global CellP so the
     * arena's "valid" state remains consistent across mixed P/CellP consumers. */
    gpu_particles_arena_set_site("gpu_build_symmetric_neighbor_list");
    gpu_particles_arena_acquire(num_total, P_host, CellP);
    struct particle_data *P_shared = gpu_particles_arena_P();

    /* Build GPU CSR — share gas-only SIDX with density via the step-persistent cache.
     *
     * search_radius_factor is applied to BOTH the i-side query radius and the
     * j-side kernel radius — a genuinely symmetric scaled search. This repairs
     * the j-side under-search in the shared symlist (hydro-gradient Velocity_hat
     * wide filter under TURB_DIFF_DYNAMIC).
     *
     * RADIUS SEMANTICS: pass EXPLICIT raw per-active radii (P[i].KernelRadius)
     * so search_radius_factor multiplies the RAW kernel radius, as the runner
     * Spec path does with its explicit fac*raw radii. */
    std::vector<double> symlist_raw_radii((num_active > 0) ? (size_t)num_active : 1);
    for(int aa = 0; aa < num_active; aa++) {
        symlist_raw_radii[aa] = (double) P_shared[active_indices[aa]].KernelRadius;
    }
    gpu_neighbor_list_t gpu_nl;
    gpu_ngb_list_build(P_shared, num_total, active_indices, num_active,
                       NGB_SEARCH_SYMMETRIC, 1 /* gas only */, &gpu_nl, gpu_step_sidx_ptr(),
                       search_radius_factor, symlist_raw_radii.data(), NULL, "symlist",
                       search_radius_factor /* j_kernel_radius_scale */);


    /* Copy CSR into mymalloc neighbor_list_t */
    out->num_active = num_active;
    out->total_pairs = gpu_nl.total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(gpu_nl.total_pairs > 0 ? gpu_nl.total_pairs : 1) * sizeof(int));
    /* gpu_nl.offsets is SharedSpace (UVM) → host memcpy is fine.
     * gpu_nl.neighbors is DEVICE_SPACE (CudaSpace) → must use deep_copy, not host memcpy. */
    memcpy(out->offsets, gpu_nl.offsets, (size_t)(num_active + 1) * sizeof(int64_t));
    if(gpu_nl.total_pairs > 0) {
        Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            h_neighbors(out->neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            d_neighbors(gpu_nl.neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::deep_copy(h_neighbors, d_neighbors);
    }


    /* Free the GPU list (the kept gas index is not the list's).
     * Arena is intentionally not released — subsequent gradient/hydro callers
     * benefit from the fast-path skip. */
    gpu_ngb_list_free(&gpu_nl);

}


/* Cross-type variant: i-side is the caller's active_indices (any type), j-side
   is filtered by j_type_bitmask. Caller supplies per-active search radii (so
   e.g. KernelRadiusDM or AGS_Hsml can be used instead of P[i].KernelRadius).
   Returns a neighbor_list_t in the mymalloc format, same as the symmetric
   variant — consumers look identical. */
void gpu_build_cross_type_neighbor_list(struct particle_data *P_host, int num_total,
                                        int *i_active_indices, int num_active,
                                        const double *i_search_radii_host,
                                        int j_type_bitmask, int search_mode,
                                        neighbor_list_t *out)
{
    /* Use the per-step particle arena to avoid a redundant full-NumPart memcpy
     * (see gpu_build_symmetric_neighbor_list for rationale). */
    gpu_particles_arena_set_site("gpu_build_cross_type_neighbor_list");
    gpu_particles_arena_acquire(num_total, P_host, CellP);
    struct particle_data *P_shared = gpu_particles_arena_P();

    /* Build GPU CSR with explicit per-i radii and j-side type filter */
    gpu_neighbor_list_t gpu_nl;
    gpu_ngb_list_build(P_shared, num_total, i_active_indices, num_active,
                       search_mode, j_type_bitmask, &gpu_nl, NULL,
                       1.0 /* search_radius_factor */, i_search_radii_host, NULL, "xtype");

    /* Copy CSR into mymalloc neighbor_list_t */
    out->num_active = num_active;
    out->total_pairs = gpu_nl.total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(gpu_nl.total_pairs > 0 ? gpu_nl.total_pairs : 1) * sizeof(int));
    /* See gpu_build_symmetric_neighbor_list for why neighbors needs deep_copy. */
    memcpy(out->offsets, gpu_nl.offsets, (size_t)(num_active + 1) * sizeof(int64_t));
    if(gpu_nl.total_pairs > 0) {
        Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            h_neighbors(out->neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            d_neighbors(gpu_nl.neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::deep_copy(h_neighbors, d_neighbors);
    }

    /* Free GPU temporaries.  Arena is intentionally retained for subsequent callers. */
    gpu_ngb_list_free(&gpu_nl);
}


/* ===================================================================== */
/* Device receiver traversal for request-driven ghost discovery.          */
/* ===================================================================== */
/* The supply-rank half: each received envelope carries a peer's query plus the
 * start nodes that peer's walk reached in THIS rank's tree, and the answer is
 * the set of local particles the query admits.  This is the device form of
 * mode_b_walk_from_start_nodes: it finds the particles that may be neighbours once
 * drifted, and both backends hand those to the same accept
 * (gx_send_set_accept_rows), which drifts them and keeps the exact set.
 *
 * Node geometry comes from the SoA mirror, never the managed Nodes[]/Extnodes[]
 * arrays: streaming those from a kernel is memory-bound to the point of erasing
 * the win.  Leaf fields are staged into a compact array first, for the same
 * reason and because the walk-export path does not build the tile geometry that
 * other device searches read.
 *
 * The device walk cannot drift a stale node the way the host walk does (that
 * needs a lock).  It does not have to: it only runs when the node sweep has
 * certified the whole tree current at this time AND no host lazy drift has
 * happened since, which is the same pair of conditions the device gravity walk
 * uses.  The certification is checked immediately before launch, because the
 * discovery walks themselves arm the lazy-drift latch — a sender walk earlier in
 * the same exchange can withdraw device legality.  Ordering violations therefore
 * make the test fail and route to the host, never corrupt a result. */

/* Per-leaf staged record.  `type` is negative for a particle the receiver can never
 * send -- no mass, or not in the supply pool -- so the kernel tests one integer where
 * the host tests three conditions across two arrays.  `pos` is where the particle is
 * once drifted to the current time, and `half_width` how far from there it can be
 * along any axis (particle_motion_envelope), negative when its motion cannot be
 * bounded.  The leaf test adds its own rounding allowance (motion_envelope_test_slack),
 * so the device's candidates are a superset of what the host's exact test accepts. */
struct gx_recv_leaf_t {
    double pos[3];
    double half_width;
    int    type;
};

/* The leaf half of the receiver walk: record a locally-owned particle that may be a
 * neighbour of the query once drifted to the current time: kept if the box round the
 * position the drift will give it reaches the query sphere; one whose motion cannot be bounded is kept.  This is
 * discovery only -- the candidates go to gx_send_set_accept_rows, which drifts them
 * and makes the exact decision -- so the walk may over-include but never decides.
 *
 * Candidates past `cap` are counted but not written, so the caller can re-run the
 * row against a buffer sized to the true count. */
struct GxRecvCandidates {
    const struct gx_recv_leaf_t *leaves;
    unsigned int supply_mask;
    int         *out;
    int          cap;
    int          n_found;

    KOKKOS_INLINE_FUNCTION
    void visit(int j, double qx, double qy, double qz, double reach)
    {
        const struct gx_recv_leaf_t &lf = leaves[j];
        if(lf.type < 0) {return;}
        if(!(supply_mask & (1u << lf.type))) {return;}
        int keep = 1;
        if(lf.half_width >= 0.0) {
            const double q[3] = {qx, qy, qz};
            const double hw = lf.half_width + motion_envelope_test_slack(lf.pos, q, reach);
            keep = gx_boxpair_overlap_wrap_and_test(lf.pos[0] - qx, lf.pos[1] - qy, lf.pos[2] - qz,
                                                    hw, hw, hw, reach, reach * reach);
        }
        if(keep) {
            if(n_found < cap) {out[n_found] = j;}
            n_found++;
        }
    }
};

/* One envelope's traversal.  Returns the number of candidates, writing the first
 * `cap` of them to `out` (the caller re-runs with the true count when a row
 * overflows its scratch slot).  Mirrors mode_b_walk_impl's three index classes
 * exactly; see mesh/mode_b_local_walker.cc for the host original.
 *
 * `anomaly` reports the one state the host treats as fatal: an index in the gap
 * between the particle slots and the node base, which belongs to neither and
 * means the tree is malformed.  The host stops the run there, so the device
 * cannot simply stop walking -- that would silently truncate an envelope.  It
 * records the state and the caller reproduces the host's stop. */
KOKKOS_INLINE_FUNCTION
static int gx_recv_walk_one(const struct gx_export_envelope_t &env,
                            unsigned int supply_mask,
                            int type_mask_trusted,
                            const Vec3<MyFloat> *node_center,
                            const MyFloat *node_len,
                            const int *node_sibling,
                            const int *node_nextnode,
                            const unsigned int *node_bitflags,
                            const int *nextnode_aux,
                            const struct gx_recv_leaf_t *leaves,
                            int tree_base, int tree_slots, int node_capacity,
                            int foreign_base, int pseudo_start, int num_local,
                            int *anomaly, int *out, int cap)
{
    GxDeviceTreeView tree;
    tree.node_center    = node_center;
    tree.node_len       = node_len;
    tree.node_sibling   = node_sibling;
    tree.node_nextnode  = node_nextnode;
    tree.node_bitflags  = node_bitflags;
    tree.type_mask_trusted = type_mask_trusted;   /* read on the host; a device read of the
                                                   * host global would be silently wrong */
    tree.nextnode_aux   = nextnode_aux;
    tree.node_base      = tree_base;
    tree.particle_slots = tree_slots;
    tree.local_particle_slots = num_local;
    tree.node_capacity  = node_capacity;
    tree.foreign_base   = foreign_base;
    tree.pseudo_start   = pseudo_start;

    GxRecvCandidates emit;
    emit.leaves      = leaves;
    emit.supply_mask = supply_mask;
    emit.out         = out;
    emit.cap         = cap;
    emit.n_found     = 0;

    gx_device_tree_walk(env, tree, emit, anomaly, supply_mask);
    return emit.n_found;
}

/* Envelopes are processed in fixed-size batches, and each batch's accepted pairs
 * are emitted through one buffer of a fixed size, so neither the scratch nor the
 * output footprint is set by however many envelopes arrived or how many pairs
 * they admit.  Rows are independent and the answer is a set, so splitting
 * changes nothing about the result. */
static const long GX_RECV_BATCH  = 4096;   /* envelopes per walk launch */
static const int  GX_RECV_STRIDE = 512;    /* per-row scratch slots before a re-walk */
/* One budget serves the scratch and the emission buffer; a batch whose pairs do
 * not fit is emitted in several passes rather than abandoned, because declining
 * on a dense batch would send exactly the large-N case this exists to serve back
 * to the host. */
static const long GX_RECV_PAIR_FLOOR = GX_RECV_BATCH * (long)GX_RECV_STRIDE;

/* Exhaustion is reported by returning NULL, so the caller's NULL check decides
   what to do. */
template <class T>
static T *gx_recv_alloc(const char *label, size_t count)
{
    return (T *) gizmo_gpu_alloc_device(count * sizeof(T), label);
}


/* Describe this rank's tree to the device walk, or say why it cannot be
 * described.  Returns 0 with `out` filled, or 1 with `out` untouched.
 *
 * Every caller of the device traversal needs the same answer to the same two
 * questions -- is there a mirror, and does it cover everything a walk can reach
 * -- and the index arithmetic that separates the three index classes is the same
 * arithmetic in every case.  Written once, because a second copy would be a
 * second place for the class boundaries to be got wrong, and a walk that reads
 * one boundary wrong does not fail, it answers short.
 *
 * `local_particle_slots` is the CALLER's, and it is the one thing here that is
 * genuinely per-caller: a walk answering for this rank alone passes the owned
 * count, while one meant to see imported ghosts too passes the full slot count.
 * See mesh/device_tree_walk.h.
 *
 * What this does NOT decide is whether the geometry is CURRENT.  That is a
 * policy question -- who is expected to have drifted the nodes, and what to do
 * when nobody has -- and the answer differs by caller, so each one settles it at
 * its own site rather than inheriting a rule written for another. */
int gx_device_tree_view_build(struct GxDeviceTreeView *out, int local_particle_slots,
                              const char *caller)
{
    /* A tree mirror is required, and it exists whenever a tree does: the build
     * pipeline that fills it runs inside force_treebuild, which every
     * configuration performs because neighbour search needs the tree. */
    const struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa || !soa->center || !soa->len || !soa->sibling || !soa->nextnode ||
       !soa->bitflags || !soa->nextnode_aux ||
       All.TreeNodeIndexBase <= 0 || Numnodestree <= 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d has no usable tree mirror; answering on the host\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }

    const int tree_base     = All.TreeNodeIndexBase;
    const int tree_slots    = All.TreeParticleSlots;
    const int node_capacity = gpu_gravity_tree_capacity();
    const int foreign_base  = tree_base + MaxNodes;
    const int pseudo_start  = tree_base + MaxNodes + MaxForeignNodes;
    /* The walk may reach any foreign node that was installed, so the mirror has to
     * cover the foreign slots that have storage behind them -- AllocatedForeignNodes,
     * this rank's actual import, which is what the mirror is sized to.  NOT
     * MaxForeignNodes: that is the shared INDEX ceiling used above to place the
     * pseudo-particle region, it is the worst rank's import rather than this one's,
     * and no node is ever installed in the gap between the two, so nothing points
     * there.  Testing the ceiling would decline on every rank whose import is
     * smaller than the largest, which is nearly all of them.  A short mirror is a
     * precondition failure, not a malformed tree: decline, and the host answers --
     * but say so, because a run that quietly answered everything on the host would
     * otherwise look exactly like a run where the device did the work.  This says
     * nothing about whether the geometry is current; that is the caller's. */
    if(node_capacity < MaxNodes + AllocatedForeignNodes ||
       soa->nextnode_aux_size < tree_slots + NTopleaves) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d tree mirror covers %d nodes and %d particle links, short of the %d nodes and %d links the walk can reach; answering on the host\n",
                   caller, ThisTask, node_capacity, soa->nextnode_aux_size,
                   MaxNodes + AllocatedForeignNodes, tree_slots + NTopleaves);
            fflush(stdout);
        }
        return 1;
    }

    /* WIDEN-ON-OPEN inputs (landing 4). The walk re-bounds each node it opens from
       the pair (len, node_ti) plus vmax, instead of requiring a sweep to have
       advanced every node first. Left null if the mirror does not carry them, in
       which case the walk opens on the stored length alone. */
    out->node_vmax            = soa->vmax;
    out->node_ti              = soa->node_ti;
    out->ti_now               = All.Ti_Current;
    gx_walk_drift_tables_refresh();
    out->drift_tables_ok      = g_walk_drift_tables_ok;
    if(g_walk_drift_tables_ok) {out->drift_tables = g_walk_drift_tables;}
    out->node_center          = soa->center;
    out->node_len             = soa->len;
    out->node_sibling         = soa->sibling;
    out->node_nextnode        = soa->nextnode;
    out->node_bitflags        = soa->bitflags;
    out->nextnode_aux         = soa->nextnode_aux;
    out->node_base            = tree_base;
    out->particle_slots       = tree_slots;
    out->local_particle_slots = local_particle_slots;
    out->node_capacity        = node_capacity;
    out->foreign_base         = foreign_base;
    out->pseudo_start         = pseudo_start;
    out->type_mask_trusted    = TypePresenceMaskTrusted;
    return 0;
}

/* ---- The touched-set workspace ------------------------------------------
 *
 * Rank-local and persistent, so a call that reaches a thousand leaves pays for a
 * thousand rather than for the whole rank.  The contract is in the header; what
 * follows is the storage and the three operations that use it.
 *
 * The generation stamp is the same one the Mode A neighbour-list hook uses a few
 * hundred lines above (`pool_seen` / `pool_seen_gen`) -- allocated once, never
 * cleared, a wrap re-zeroing it -- moved onto the device because the recorder is
 * a kernel.  Nothing scans it: entries are touched only for leaves a walk
 * actually reaches, which is what keeps this inside the rule that a step with a
 * handful of active particles does no work proportional to the rank. */
static struct GxTouchedSet g_touched_set;

int gx_touched_set_ensure(int local_particle_slots)
{
    if(local_particle_slots <= 0) {return 1;}
    if(g_touched_set.capacity >= local_particle_slots && g_touched_set.seen) {return 0;}

    /* Grown, not resized in place: the stamps say which slots a PREVIOUS
     * generation claimed, and the slots have been re-indexed underneath them by
     * whatever grew the rank.  Re-zeroing and restarting the generation is the
     * honest response; carrying stamps across would let a stale one suppress a
     * particle that genuinely needs drifting. */
    /* Ordering, not cleanup -- the same rule gpu_spatial_index_free states above:
       kokkos_free does not synchronize, so releasing storage a kernel may still
       be reading is a use-after-free. The recording kernels are fenced by
       drift_and_mark before this is ever reached, but the fence belongs WITH the
       release so no future caller has to know that. Skipped when there is
       nothing to release. */
    if(g_touched_set.seen || g_touched_set.list || g_touched_set.counter) {Kokkos::fence();}
    if(g_touched_set.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.seen);}
    if(g_touched_set.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.list);}
    if(g_touched_set.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.counter);}
    g_touched_set = GxTouchedSet{};

    unsigned int *seen    = (unsigned int *) ngl_alloc_shared((size_t)local_particle_slots * sizeof(unsigned int), "touched_set_seen");
    int          *list    = (int *)          ngl_alloc_shared((size_t)local_particle_slots * sizeof(int),          "touched_set_list");
    int          *counter = (int *)          ngl_alloc_shared(sizeof(int),                                         "touched_set_counter");
    if(!seen || !list || !counter) {
        if(seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(seen);}
        if(list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(list);}
        if(counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(counter);}
        return 1;
    }
    for(int k = 0; k < local_particle_slots; k++) {seen[k] = 0u;}
    *counter = 0;

    g_touched_set.seen     = seen;
    g_touched_set.list     = list;
    g_touched_set.counter  = counter;
    g_touched_set.capacity = local_particle_slots;
    g_touched_set.gen      = 0u;
    return 0;
}

/* Lifecycle faults, counted per rank.  A caller that is permanently declining to its own safe
 * route, or a retire that found work nobody drained, is otherwise indistinguishable from the
 * mechanism simply working. */
static long long g_touched_refused_epochs  = 0;
static long long g_touched_retire_faults   = 0;
/* Whether this epoch's claims have already been waited for.  A consume fences; retire needs to
 * know so it does not fence a second time on every successful call. */
static int       g_touched_consumed_this_epoch = 0;

/* Each fault owns its latch.  One shared latch let the first report -- typically the benign
 * "another phase still holds it" refusal -- silence the serious one for the rest of the run,
 * so a recorder retired with claims nobody consumed could pass unmentioned. */
static void touched_report_once_(int *latch, const char *what)
{
    if(!latch || *latch) {return;}
    *latch = 1;
    printf("touched set: task %d: %s\n", ThisTask, what);
    fflush(stdout);
}

/* Open an epoch for `owner`.  Returns 0 when the recorder is now this owner's, nonzero when it
 * could not be handed over -- and a refusal leaves the generation, the cursor and the owner
 * exactly as they were, because the whole point is not to disturb whoever still holds them.
 *
 * THE ADMISSION CONDITION IS THAT NOBODY HOLDS THE RECORDER, and nothing else.  The cursor
 * cannot serve: a consume resets it while the generation deliberately lives on across the
 * remaining passes of the same call, so a zero cursor mid-call means "this pass is drained",
 * never "the recorder is free".  Admitting on a zero cursor would let a second owner take the
 * generation between two passes of a live call, and the first owner's next claim would then
 * arrive out of phase -- which is why the holder, not the cursor, is the question.
 *
 * ⚠ A refusal is NOT self-healing.  The generation is deliberately not advanced, so the
 * previous epoch's stamps still stand and the particles carrying them can no longer re-enter
 * the list.  Every caller therefore has to ACT on a refusal -- the fused walk declines the
 * whole call collectively -- and none of them may
 * proceed as if the epoch had opened. */
int gx_touched_set_begin_epoch_owned(int owner)
{
    if(!g_touched_set.seen || !g_touched_set.list || !g_touched_set.counter) {
        g_touched_refused_epochs++;
        return 1;
    }
    if(g_touched_set.owner != GX_TOUCHED_OWNER_NONE) {
        g_touched_refused_epochs++;
        static int latch_busy = 0;
        touched_report_once_(&latch_busy,
                             "an epoch was requested while another phase still held the recorder; "
                             "the requester takes its own fallback");
        return 1;
    }
    /* Zero is the never-claimed value, so a wrap has to skip it AND clear the
     * stamps -- otherwise a slot still carrying the old maximum would read as
     * claimed by the new generation and its particle would silently go
     * undrifted. */
    if(++g_touched_set.gen == 0u) {
        for(int k = 0; k < g_touched_set.capacity; k++) {g_touched_set.seen[k] = 0u;}
        g_touched_set.gen = 1u;
    }
    *g_touched_set.counter = 0;
    g_touched_set.owner    = owner;
    g_touched_consumed_this_epoch = 0;
    Kokkos::memory_fence();   /* the epoch is open before any claim in it is visible */
    return 0;
}

/* The fused walk's entry point.  It holds the recorder for its whole call, across every
 * discovery pass, and retires at the end of that call.  It reports rather than swallows,
 * because its caller (the walk preparation) already has the collective decline that a
 * refusal needs. */
int gx_touched_set_begin_call(void)
{
    return gx_touched_set_begin_epoch_owned(GX_TOUCHED_OWNER_FUSED_WALK);
}

/* End an epoch.  MUST be called after the owner's last consume, which is also what fences the
 * recording kernel -- retiring with a kernel still in flight would leave claims arriving into
 * an epoch that no longer exists.
 *
 * Both faults are loud.  Retiring someone else's epoch means two phases disagree about who is
 * running, and returning quietly would leave the real holder's epoch closed underneath it.  A
 * nonzero cursor means claims were recorded and never drained, i.e. particles the walk reached
 * were never drifted -- the silent staleness this recorder exists to prevent, so it is
 * reported rather than cleared quietly. */
void gx_touched_set_retire(int owner)
{
    if(g_touched_set.owner != owner) {
        g_touched_retire_faults++;
        static int latch_wrong_owner = 0;
        touched_report_once_(&latch_wrong_owner, "a phase tried to retire a recorder epoch it does not hold");
        return;
    }
    /* Wait for the recording kernel before reading the cursor -- but only when nothing else
     * already has.  A consume fences, so after one the cursor is readable and a second wait
     * would be a host-side device synchronize on the critical path of every fused call.  It is
     * the path with NO consume that the check below exists for, and there nothing has waited,
     * so reading the cursor would report the epoch clean whenever the kernel's writes were
     * merely not visible yet -- blind in exactly the case it is for.  A full fence, not a
     * memory fence: the kernel may still be running, so ordering the accesses is not enough. */
    if(!g_touched_consumed_this_epoch) {Kokkos::fence();}
    if(g_touched_set.counter && *g_touched_set.counter != 0) {
        g_touched_retire_faults++;
        static int latch_undrained = 0;
        touched_report_once_(&latch_undrained,
                             "a recorder epoch retired with claims nobody consumed; "
                             "the particles they name were never drifted");
        *g_touched_set.counter = 0;
    }
    g_touched_set.owner = GX_TOUCHED_OWNER_NONE;
    Kokkos::memory_fence();   /* the epoch is closed before the next one can open */
}

long long gx_touched_set_refused_epochs(void) {return g_touched_refused_epochs;}
long long gx_touched_set_retire_faults(void)  {return g_touched_retire_faults;}

struct GxTouchedSet gx_touched_set_view(void) {return g_touched_set;}

/* Give the workspace back.  Called once, at shutdown, before Kokkos is
 * finalized -- Kokkos must not be torn down while an allocation it is tracking
 * is still owned.  Deliberately NOT hung off the tree epoch or the domain
 * decomposition: this storage is persistent on purpose, and freeing it there
 * would turn a once-per-run allocation into per-rebuild churn, which is the
 * cost the generation stamp exists to avoid. */
void gx_touched_set_release(void)
{
    /* Ordering, not cleanup (gpu_spatial_index_free states the rule): kokkos_free
       does not synchronize. This runs on the controlled-stop path as well as the
       normal one, where a kernel may well still be in flight, so incidental
       completion is not an ownership contract. Fence once if anything is held. */
    if(g_touched_set.seen || g_touched_set.list || g_touched_set.counter) {Kokkos::fence();}
    if(g_touched_set.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.seen);}
    if(g_touched_set.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.list);}
    if(g_touched_set.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.counter);}
    g_touched_set = GxTouchedSet{};
}

/* ============================================================================
 * The motion-target set and the raise it feeds (declared in gpu_neighbor_list.h).
 * ========================================================================== */
static struct GxMotionTargetSet g_motion_targets;
static int g_motion_targets_armed = 0;

int gx_motion_target_ensure(int local_particle_slots)
{
    if(local_particle_slots <= 0) {return 1;}
    if(g_motion_targets.capacity >= local_particle_slots && g_motion_targets.seen) {return 0;}
    gx_motion_target_release();
    unsigned int *seen    = (unsigned int *) ngl_alloc_shared((size_t)local_particle_slots * sizeof(unsigned int), "motion_target_seen");
    int          *list    = (int *)          ngl_alloc_shared((size_t)local_particle_slots * sizeof(int),          "motion_target_list");
    int          *counter = (int *)          ngl_alloc_shared(sizeof(int),                                         "motion_target_counter");
    if(!seen || !list || !counter) {
        if(seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(seen);}
        if(list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(list);}
        if(counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(counter);}
        return 1;
    }
    for(int k = 0; k < local_particle_slots; k++) {seen[k] = 0u;}
    *counter = 0;
    g_motion_targets.seen = seen; g_motion_targets.list = list; g_motion_targets.counter = counter;
    g_motion_targets.capacity = local_particle_slots; g_motion_targets.gen = 0u;
    return 0;
}

void gx_motion_target_begin_call(void)
{
    if(++g_motion_targets.gen == 0u) {
        for(int k = 0; k < g_motion_targets.capacity; k++) {g_motion_targets.seen[k] = 0u;}
        g_motion_targets.gen = 1u;
    }
    if(g_motion_targets.counter) {*g_motion_targets.counter = 0;}
}

struct GxMotionTargetSet gx_motion_target_view(void) {return g_motion_targets;}

void gx_motion_target_mark_host(int j)
{
    /* Serial host callers only (the writeback apply loop, a host module's own
     * loop); the device mark is the atomic form. */
    struct GxMotionTargetSet &ts = g_motion_targets;
    if(!ts.seen || j < 0 || j >= ts.capacity) {return;}
    if(ts.seen[j] == ts.gen) {return;}
    ts.seen[j] = ts.gen;
    const int slot = (*ts.counter)++;
    if(slot < ts.capacity) {ts.list[slot] = j;}
}

void gx_motion_target_set_armed(int armed) {g_motion_targets_armed = armed;}
int  gx_motion_target_armed(void) {return g_motion_targets_armed;}

void gx_motion_target_consume(void)
{
    if(!g_motion_targets.counter || !g_motion_targets.list) {return;}
    Kokkos::fence();   /* the markers may be device kernels */
    const int claimed = *g_motion_targets.counter;
    const int n = (claimed < g_motion_targets.capacity) ? claimed : g_motion_targets.capacity;
    if(n > 0) {gizmo_motion_bound_raise(g_motion_targets.list, n);}
    *g_motion_targets.counter = 0;
}

void gx_motion_target_release(void)
{
    if(g_motion_targets.seen || g_motion_targets.list || g_motion_targets.counter) {Kokkos::fence();}
    if(g_motion_targets.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.seen);}
    if(g_motion_targets.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.list);}
    if(g_motion_targets.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.counter);}
    g_motion_targets = GxMotionTargetSet{};
}

void gizmo_motion_bound_raise(const int *idx, int n)
{
    if(n <= 0 || !idx) {return;}
    /* The tree records which top-level nodes changed, for the exchange at the
     * next tree-update phase. */
    gravity_note_motion_bound(idx, n);
    /* The kept gas neighbour index follows the same velocities. */
    gpu_step_sidx_raise_motion(idx, n);
}

void gx_touched_set_drift_and_mark(integertime time1)
{
    if(!g_touched_set.counter || !g_touched_set.list) {return;}

    /* The recorder is a device kernel writing into shared space, so its writes
     * are not visible to the loop below until it has finished.  Without this the
     * count reads as whatever was there when the launch returned. */
    Kokkos::fence();

    /* The cursor counts every claim, including any the list had no room for, so
     * it is what the append attempted rather than what the list holds.  The two
     * agree unless the generation stamp and the list disagree about how many
     * owned slots there are, which the ensure makes impossible -- but the
     * consequence of trusting the cursor if they ever did is a read past the
     * allocation, here, several kernel launches before the anomaly that records
     * the overflow is ever looked at.  So the allocation bounds the read, and the
     * anomaly stays the thing that reports it. */
    g_touched_consumed_this_epoch = 1;   /* the fence above has waited for this epoch's claims */
    const int claimed = *g_touched_set.counter;
    const int n = (claimed < g_touched_set.capacity) ? claimed : g_touched_set.capacity;
    /* First claims for THIS pass. Reported separately from the pass and call
       counts so a reduction factor is read rather than inferred from a quotient
       whose denominator has to be guessed. */
    if(n <= 0) {return;}

    /* Which of these are actually behind is drift_particles_batch's question and
     * it already answers it -- deciding it here as well would be a second place
     * owning the same test.  So it is asked to hand its compaction back, in
     * place: the recorded list becomes the advanced list, and `n_drifted` is how
     * much of it is live.  The status it returns says only that a controlled stop
     * is pending somewhere, so nothing branches on it. */
    int n_drifted = 0;
    (void) drift_particles_batch(g_touched_set.list, n, time1,
                                 g_touched_set.list, &n_drifted);

    /* drift_particle rescales KernelRadius, so a particle it ADVANCED owes a
     * dirty mark -- and only those.  Marking everything the walk recorded would
     * mark the already-current majority too: measured on a mixed-timebin vehicle,
     * only 19.3% of recorded particles were behind, so that is several times the
     * cache invalidation the work actually justifies, and it can push the dirty
     * tracker over its promote-to-all threshold for nothing. */
    if(n_drifted > 0) {gizmo_mark_kernel_radius_dirty_indices(g_touched_set.list, n_drifted);}

    *g_touched_set.counter = 0;
}


/* Put this rank into the state a fused device walk needs, and describe its tree.
 * Returns 0 with `out` filled, or 1 with the walk declined and the host to answer.
 *
 * A fused walk evaluates the pair kernel at the leaf it just reached, so unlike
 * a discovery walk it reads particle fields, and unlike the host walker it
 * cannot drift a stale one when it gets there -- that needs a lock.  So both the
 * particles and the node geometry have to be current BEFORE the launch, and this
 * is where that is arranged, once, ahead of any discovery round.
 *
 * The particle half is NOT arranged here.  It used to be a full-rank drift, on
 * the argument that discovering the reached set first would mean doing the walk
 * twice and buy nothing.  Both halves of that were wrong at small active counts:
 * a call below ten thousand actives reaches on the order of a thousand leaves
 * while the drift advanced essentially the entire local pool, and the second
 * traversal costs a small fraction of the drift it removes.  Discovery is now
 * what decides which particles are brought current, per pass, at the three
 * evaluation sites in mesh/neighbor_loop_runner.cc -- the two shapes this rank's
 * own self walk takes, and the queries its peers sent.
 *
 * What this preparation still owes the walk is everything that cannot be
 * discovered: the node geometry, which decides where the walk goes and so cannot
 * be repaired from what it reached.  It sweeps the whole tree when nothing else
 * has.  The receiver walk declines instead of sweeping when gravity is compiled
 * in, on the ground that it would drift the entire tree on behalf of a walk that
 * touches part of it; that reasoning applies here too, and the only thing
 * standing against it is that a fused walk cannot take the lock a host walk uses
 * to drift a node when it arrives.  Closing that gap is a separate piece of
 * work.  The one state a sweep cannot repair is a host lazy drift that has
 * already advanced nodes at this time: a sweep skips nodes that are current, so
 * their mirrors stay behind, and the host has to answer. */
int gx_device_fused_walk_prepare(struct GxDeviceTreeView *out, const char *caller)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    const int num_local = ghost_get_num_local();
    if(num_local <= 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d owns no local particles, so it cannot answer on the device; this call falls back for every rank\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }

    if(gx_device_tree_view_build(out, num_local, caller) != 0) {return 1;}

    /* A walk from the root dereferences the root node before it tests anything,
     * so an empty or half-built tree is not a walk that returns nothing -- it is
     * an unbounded read.  The node capacity cannot catch it, because that is the
     * mirror's allocation rather than the number of nodes actually built, so a
     * built count of zero passes every bound test in the traversal.  Test the
     * built count itself, here, where declining is still free. */
    if(Numnodestree <= 0 || out->node_base >= out->pseudo_start) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d has %d built tree nodes; a walk from the root has nothing to enter, answering on the host\n",
                   caller, ThisTask, Numnodestree);
            fflush(stdout);
        }
        return 1;
    }

    /* ---- Everything above is a CAPABILITY question, answered from state that
     * already exists: does this rank own particles, is there a mirror, does it
     * cover what a walk can reach, is there a root to enter.  Everything below
     * REPAIRS state and is the expensive half.
     *
     * The order matters and is not incidental.  A rank that cannot describe its
     * tree declines, and under the collective vote that decision pulls the WHOLE
     * call back to Mode A -- so any repair performed before the checks is work
     * paid for an answer that is then thrown away, on every rank.  Ask first,
     * repair second.
     *
     * Nothing below invalidates the view built above: the drift neither creates
     * nor destroys particles, and the node sweep writes through the mirror
     * arrays the view already points at rather than reallocating them. */

    /* The walk reads particle fields and cannot drift a stale one when it gets
     * there.  It does NOT follow that the whole rank has to be current: what the
     * walk reads is the leaves it reaches, and a call below ten thousand actives
     * reaches on the order of a thousand of them out of half a million.  So the
     * particles are brought current per discovery pass, against the set that pass
     * actually recorded, at mesh/neighbor_loop_runner.cc's three evaluation sites --
     * which is the same three-stage shape the host backend at those sites has
     * always had, and the same one move_particles, the neighbour-list hook, the
     * Mode B walker and the ghost send-list certify already use.
     *
     * Nothing is voted on here for the particles, and that is deliberate.  Every
     * rank carries the same obligation and discharges it the same way; the drift
     * enters no collective and its only failure report is that a controlled stop
     * is already pending, which is a property of the run rather than of this
     * rank.  A vote would create the divergence it was meant to prevent.
     *
     * The workspace the recording writes into IS a capability, and is arranged
     * here for exactly that reason: it can fail rank-locally, so it belongs where
     * declining is still collective. */
    if(gx_touched_set_ensure(out->local_particle_slots) != 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d could not hold the touched-set workspace for %d local particles; this call falls back for every rank\n",
                   caller, ThisTask, out->local_particle_slots);
            fflush(stdout);
        }
        return 1;
    }
    /* The epoch is a capability too, and it is refused only when another phase still holds the
     * recorder.  Declining here routes it through the same collective vote as a workspace that
     * could not be allocated, so every rank falls back together -- and a refusal must never be
     * walked past, because it leaves the previous epoch's stamps standing. */
    if(gx_touched_set_begin_call() != 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d could not open a touched-set epoch; this call falls back for every rank\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }
    /* From here the epoch is OPEN, and this function owns it until it reports success.  The
     * caller arms its own guard only once we return 0, so anything that leaves by another door
     * has to close the epoch itself -- otherwise the recorder stays held, every later call is
     * refused, and the run is pinned to its fallback while looking exactly like the mechanism
     * working.  That is the failure the caller's guard was added to prevent, and it reappears
     * here because opening and closing sat with different owners.  Holding it in a scope guard
     * puts both with whoever opened it, so a return added below cannot reintroduce it. */
    struct FusedPrepareEpoch {
        int held = 1;
        ~FusedPrepareEpoch() {if(held) {gx_touched_set_retire(GX_TOUCHED_OWNER_FUSED_WALK);}}
    } prepare_epoch;

    const int nodes_already_current = gpu_gravity_tree_nodes_current_at(All.Ti_Current) ? 1 : 0;
    /* Nothing dirtied them, so the span the census measures starts here. */
    if(!nodes_already_current) {
        /* A host lazy drift already advanced nodes at this time.  The sweep skips
         * nodes that are current, so their mirrors would stay behind it, and that
         * state cannot be repaired here.  Expected on some calls rather than
         * exceptional -- but reported, because the fallback is collective and an
         * unexplained absence of the device path is indistinguishable from a
         * device path that ran. */
        /* A host lazy drift may have advanced nodes at this time without writing
         * their mirrors, and the ordinary sweep skips exactly those.  Mode-D used
         * to answer that by forcing the variant that rewrites EVERY mirror --
         * ~1.99M of them, to serve a walk that opens a few thousand nodes.
         *
         * It no longer needs to, because the walk now carries its own answer to
         * both classes of staleness:
         *   (a) advanced-but-unmirrored -- the node dirty set recorded exactly
         *       which nodes those are, and repairing them is O(Ndirty), measured
         *       at ~9,200 per rank per span against ~1.99M mirrors;
         *   (b) never-advanced -- widen-on-open re-bounds the node from
         *       (len, node_ti, vmax), so it needs no advancement at all.
         * The two are one design: (b) is what makes (a) sufficient.
         *
         * ⛔ This is NOT "delete the sweep".  Mode-D stops FORCING one; the sweep
         * survives for its other three callers, and the sibling device receiver
         * walk already declines for this exact reason and says so in code
         * (gpu_neighbor_list.cc, "Gravity owns the sweep").  Mode-D is adopting
         * the policy its sibling already has.
         *
         * Any doubt falls back to the old behaviour: an unsafe or overflowed
         * dirty epoch, or a mirror that cannot carry the widening inputs, makes
         * gpu_node_dirty_bring_gravity_current() return nonzero and the full sweep runs.
         * Its claims are nodes the host drifted, so bringing them current is publishing them:
         * the same routine gravity uses, one consumer of the claim list. */
        int sweep_rc = 0;
        int _sweep_needed  = 1;
        /* ⛔ Widening needs its inputs. Without them the walk would open on the
         * stored length alone, which is only legal when something else certified
         * the geometry -- so a missing table or mirror DECLINES to the old sweep
         * rather than quietly narrowing the bound. */
        const int widen_armed = (out->drift_tables_ok && out->node_ti && out->node_vmax);
        if(gpu_gravity_tree_nodes_current_at(All.Ti_Current)) {
            _sweep_needed = 0;                       /* fully certified: widening not needed */
        } else if(!widen_armed) {
            _sweep_needed = 1;                       /* cannot widen -> must sweep */
        } else if(gpu_gravity_tree_oneway_safe_at(All.Ti_Current)) {
            _sweep_needed = 0;                       /* already safe to walk */
        } else if(gpu_node_dirty_bring_gravity_current(All.Ti_Current) == 0) {
            _sweep_needed = 0;                       /* O(Ndirty) repair sufficed */
        }
        if(_sweep_needed) {
            { sweep_rc = gpu_force_drift_nodes_ex(All.Ti_Current, /*refresh_mirrors_already_current=*/1); }

            if(sweep_rc != 0) {
                static int reported = 0;
                if(!reported) {
                    reported = 1;
                    printf("%s: task %d could not sweep the node geometry current; this call falls back for every rank\n",
                           caller, ThisTask);
                    fflush(stdout);
                }
                return 1;
            }
            /* The sweep has just rewritten every mirror, so whatever the lazy drift
             * had accumulated is answered and the census span restarts. ⛔ Only on
             * the path that actually swept: the repair path leaves the span alone,
             * because it answered the dirty set rather than the whole tree. */
        }
    }

    prepare_epoch.held = 0;   /* success: the epoch passes to the caller's guard, which retires it */
    return 0;
}

int gx_device_receiver_walk(const struct gx_export_envelope_t *envelopes, long n_env,
                            const int *envelope_peer,
                            unsigned int supply_mask, int search_mode,
                            mode_b_radius_policy_t radius_policy,
                            double j_radius_scale, double safety_factor,
                            const int *j_to_pool, int npart_bound,
                            int num_pool, struct ghost_send_set *send_set)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* ONEWAY only.  A symmetric search accepts a pair on the NEIGHBOUR's radius as
     * well as the query's, so it has to bound, per node, how far the particles of
     * each type below it reach -- Extnodes[].hmax_per_type[].  Those bands are host
     * AoS only.  The per-node type-presence bits the device walk does carry answer
     * whether a type is there, which is enough to skip a node but not to decide
     * whether one of its particles reaches back, so they do not close this gap.
     * Symmetric callers are answered by the host walk. */
    if(search_mode != NGB_SEARCH_ONEWAY) {return GX_RECEIVER_DECLINED;}
    if(n_env <= 0 || num_pool <= 0 || !send_set) {return GX_RECEIVER_DECLINED;}

    const int num_local = ghost_get_num_local();
    if(num_local <= 0) {return GX_RECEIVER_DECLINED;}

    /* DISPATCH FLOOR.  Staging the local leaves costs O(num_local) whatever the
     * envelopes ask for, while the traversal it enables scales with the envelope
     * count -- so on a step with little to answer the staging is pure overhead,
     * and the corridor already loses the small and middle bins to fixed per-call
     * costs of exactly this kind.  Requiring each envelope to cover a fixed
     * number of staged particles keeps the device path out of that regime.
     * Structural and rank-local: no caller identity, no tuning knob.
     * The ratio is provisional and is what the pricing arm sets. */
    const long GX_RECV_LEAVES_PER_ENVELOPE = 32;
    if(n_env * GX_RECV_LEAVES_PER_ENVELOPE < (long)num_local) {return GX_RECEIVER_DECLINED;}

    struct GxDeviceTreeView tree_view;
    if(gx_device_tree_view_build(&tree_view, num_local, "gx_device_receiver_walk") != 0) {return GX_RECEIVER_DECLINED;}
    const int tree_base     = tree_view.node_base;
    const int tree_slots    = tree_view.particle_slots;
    const int node_capacity = tree_view.node_capacity;
    const int foreign_base  = tree_view.foreign_base;
    const int pseudo_start  = tree_view.pseudo_start;

    /* The traversal reads node geometry, which drifts, so the nodes have to be
     * current before it runs -- the host walk achieves that by drifting each
     * stale node as it reaches it, under a lock, which a kernel cannot do.
     *
     * Requiring some earlier caller to have swept them is not enough: the sweep
     * lives on the gravity path, so with self-gravity disabled nothing would
     * ever perform it and this traversal could never run at all, on problems
     * whose tree and mirror are perfectly valid.  Perform the sweep here when
     * nothing else has.  It is the same work the host walk would do node by
     * node, done once for the whole tree, and it is rank-local.
     *
     * The one case it cannot repair is a host lazy drift that has already
     * advanced nodes to this time: the sweep skips nodes that are already
     * current and so would leave their mirrors behind. Then, and only then, the
     * host answers. */
    if(!gpu_gravity_tree_nodes_current_at(All.Ti_Current)) {
        /* The geometry is current by neither route -- no sweep certified it and
         * this tree was not built at this time.  A host lazy drift at this time
         * advanced nodes without their mirrors, and a sweep skips already-current
         * nodes, so that state is unrepairable here and the host answers.
         * A tree BUILT after such a drift never reaches this branch: the build
         * rewrites every node and every mirror, and no later host walk can re-arm
         * the latch at this time because force_drift_node returns early on a node
         * that is already current. */
        if(force_host_lazy_drift_ti() == All.Ti_Current) {return GX_RECEIVER_DECLINED;}
#ifdef SELFGRAVITY_OFF
        /* No gravity walk exists to sweep the nodes in this build, so without
         * this the traversal could never run at all, however valid the tree and
         * its mirror are. */
        if(gpu_force_drift_nodes(All.Ti_Current) != 0) {return GX_RECEIVER_DECLINED;}
#else
        /* Gravity owns the sweep. Sweeping here instead would drift the whole
         * tree eagerly where the walks drift only what they touch, so leave the
         * geometry alone and let the host answer. */
        return GX_RECEIVER_DECLINED;
#endif
    }

    using DevSp = GIZMO_KOKKOS_DEVICE_SPACE;

    /* One envelope's accepted pairs are distinct local particles, so a row can
     * never exceed num_local.  Sizing the emission buffer to at least that lets
     * every row fit whole, so no batch is ever abandoned for being dense and the
     * decline paths all sit before anything is written. */
    const long pair_cap = (GX_RECV_PAIR_FLOOR > (long)num_local) ? GX_RECV_PAIR_FLOOR : (long)num_local;

    /* Device buffers.  Allocated once for the whole call and reused across every
     * batch.  A failure here is reported as a decline, not an abort: the caller
     * runs the host walk, and this window holds no collectives so no other rank
     * needs to agree. */
    /* Host buffers first: std::vector throws rather than returning null, and a
     * throw escaping this function would take the rank down between the envelope
     * exchange and the caller's reduction.  Reserving them here turns that into
     * the same decline every other resource failure produces. */
    std::vector<struct gx_recv_leaf_t>     leaf_h;
    std::vector<struct gx_export_envelope_t> env_h;
    std::vector<int>     counts_h;
    std::vector<int64_t> offsets_h;
    std::vector<int>     pairs_h;
    std::vector<struct gx_candidate_row> rows_h;
    try {
        leaf_h.resize((size_t)num_local);
        env_h.resize((size_t)GX_RECV_BATCH);
        counts_h.resize((size_t)GX_RECV_BATCH);
        offsets_h.resize((size_t)GX_RECV_BATCH);
        pairs_h.resize((size_t)pair_cap);
        rows_h.resize((size_t)GX_RECV_BATCH);
    } catch(const std::bad_alloc &) {
        printf("gx_device_receiver_walk: task %d could not reserve host staging for %d local leaves; answering on the host\n",
               ThisTask, num_local);
        fflush(stdout);
        return GX_RECEIVER_DECLINED;
    }

    struct gx_recv_leaf_t *leaf_d = gx_recv_alloc<struct gx_recv_leaf_t>("gx_recv_leaf", (size_t)num_local);
    struct gx_export_envelope_t *env_d = gx_recv_alloc<struct gx_export_envelope_t>("gx_recv_env", (size_t)GX_RECV_BATCH);
    int     *scratch_d = gx_recv_alloc<int>("gx_recv_scratch", (size_t)GX_RECV_BATCH * GX_RECV_STRIDE);
    int     *counts_d  = gx_recv_alloc<int>("gx_recv_counts",  (size_t)GX_RECV_BATCH);
    int64_t *offsets_d = gx_recv_alloc<int64_t>("gx_recv_offsets", (size_t)GX_RECV_BATCH);
    int     *pairs_d   = gx_recv_alloc<int>("gx_recv_pairs",   (size_t)pair_cap);
    int     *anomaly_d = gx_recv_alloc<int>("gx_recv_anomaly", 1);

    if(!leaf_d || !env_d || !scratch_d || !counts_d || !offsets_d || !pairs_d || !anomaly_d) {
        printf("gx_device_receiver_walk: task %d could not reserve device buffers for %d local leaves; answering on the host\n",
               ThisTask, num_local);
        fflush(stdout);
        if(leaf_d)    {Kokkos::kokkos_free<DevSp>(leaf_d);}
        if(env_d)     {Kokkos::kokkos_free<DevSp>(env_d);}
        if(scratch_d) {Kokkos::kokkos_free<DevSp>(scratch_d);}
        if(counts_d)  {Kokkos::kokkos_free<DevSp>(counts_d);}
        if(offsets_d) {Kokkos::kokkos_free<DevSp>(offsets_d);}
        if(pairs_d)   {Kokkos::kokkos_free<DevSp>(pairs_d);}
        if(anomaly_d) {Kokkos::kokkos_free<DevSp>(anomaly_d);}
        return GX_RECEIVER_DECLINED;
    }

    using UmHostI   = Kokkos::View<int*,     Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmHostI64 = Kokkos::View<int64_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI    = Kokkos::View<int*,     DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI64  = Kokkos::View<int64_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    int status_staging_failed = 0;
    /* Stage the leaf fields the walk reads.  One pass over the local particles,
     * host-side, into a compact record; the AoS is never touched from device.  The
     * position a particle that is behind will be drifted to is taken here, from the host
     * tables, so the kernel reads the record instead of the particle state and the tables. */
    {
        const integertime ti_now = All.Ti_Current;
        const struct DriftKickTableView drift_tables = drift_kick_table_view_host();
        long n_unmapped = 0;
#pragma omp parallel for schedule(static) reduction(+:n_unmapped)
        for(int j = 0; j < num_local; j++) {
            struct gx_recv_leaf_t rec;
            int in_pool = 0;
            if(j_to_pool && j < npart_bound) {
                const int pp = j_to_pool[j];
                in_pool = (pp >= 0 && pp < num_pool);
            }
            /* The pool holds every particle of positive mass (gx_send_set_accept_rows), so one
               missing from it is a stale or corrupt map, stopped below rather than hidden. */
            if(P[j].Mass > 0 && !in_pool) {n_unmapped++;}
            rec.type = (P[j].Mass > 0 && in_pool) ? (int)P[j].Type : -1;
            rec.pos[0] = (double)P[j].Pos[0]; rec.pos[1] = (double)P[j].Pos[1]; rec.pos[2] = (double)P[j].Pos[2];
            rec.half_width = 0.0;
            if(rec.type >= 0) {
                double hw = 0.0;
                const int motion = particle_motion_envelope(j, P, CellP, ti_now, &drift_tables, rec.pos, &hw);
                rec.half_width = (motion == PARTICLE_MOTION_UNBOUNDED) ? -1.0 : hw;
            }
            leaf_h[(size_t)j] = rec;
        }
        if(n_unmapped > 0) {
            printf("gx_device_receiver_walk: task %d has %ld particles of positive mass with no supply-pool slot\n",
                   ThisTask, n_unmapped);
            fflush(stdout);
            gizmo_request_controlled_stop(7738, "gx_device_receiver_walk: supply map missing particles of positive mass",
                                          __FILE__, __LINE__, __FUNCTION__);
            status_staging_failed = 1;
        }
        Kokkos::View<struct gx_recv_leaf_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            lh(leaf_h.data(), (size_t)num_local);
        Kokkos::View<struct gx_recv_leaf_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            ld(leaf_d, (size_t)num_local);
        Kokkos::deep_copy(ld, lh);
    }

    {   /* clear the malformed-tree report */
        int zero = 0;
        Kokkos::View<const int, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> zh(&zero);
        Kokkos::View<int, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>> zd(anomaly_d);
        Kokkos::deep_copy(zd, zh);
    }

    const Vec3<MyFloat>  *node_center   = tree_view.node_center;
    const MyFloat        *node_len      = tree_view.node_len;
    const int            *node_sibling  = tree_view.node_sibling;
    const int            *node_nextnode = tree_view.node_nextnode;
    const unsigned int   *node_bitflags = tree_view.node_bitflags;
    const int            *nextnode_aux  = tree_view.nextnode_aux;

    int status = status_staging_failed ? GX_RECEIVER_FAILED : GX_RECEIVER_COMPLETED;

    for(long base = 0; base < n_env && status == GX_RECEIVER_COMPLETED; base += GX_RECV_BATCH) {
        const int nb = (int)((n_env - base < GX_RECV_BATCH) ? (n_env - base) : GX_RECV_BATCH);
        for(int b = 0; b < nb; b++) {env_h[(size_t)b] = envelopes[base + b];}
        {
            Kokkos::View<struct gx_export_envelope_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                eh(env_h.data(), (size_t)nb);
            Kokkos::View<struct gx_export_envelope_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                ed(env_d, (size_t)nb);
            Kokkos::deep_copy(ed, eh);
        }

        /* Walk: bounded write into the row's scratch slot, true count returned. */
        {
            const struct gx_export_envelope_t *env_v = env_d;
            const struct gx_recv_leaf_t *leaf_v = leaf_d;
            int *scratch_v = scratch_d; int *counts_v = counts_d; int *anom_v = anomaly_d;
            const int stride = GX_RECV_STRIDE;
            const int type_mask_trusted_v = TypePresenceMaskTrusted;   /* captured by value */
            Kokkos::parallel_for("gx_recv_walk", nb, KOKKOS_LAMBDA(int b) {
                counts_v[b] = gx_recv_walk_one(env_v[b], supply_mask, type_mask_trusted_v,
                                               node_center, node_len, node_sibling,
                                               node_nextnode, node_bitflags, nextnode_aux,
                                               leaf_v,
                                               tree_base, tree_slots, node_capacity,
                                               foreign_base, pseudo_start, num_local,
                                               anom_v,
                                               &scratch_v[(size_t)b * stride], stride);
            });
            Kokkos::fence();
            gizmo_gpu_check_last_error("gx_recv_walk", nb);
        }

        int64_t batch_total = 0;
        {
            int *counts_v = counts_d; int64_t *offsets_v = offsets_d;
            Kokkos::parallel_scan("gx_recv_scan", nb,
                KOKKOS_LAMBDA(int b, int64_t &update, const bool final) {
                    const int64_t v = (int64_t)counts_v[b];
                    if(final) {offsets_v[b] = update;}
                    update += v;
                }, batch_total);
            Kokkos::fence();
        }
        if(batch_total < 0) {
            printf("gx_device_receiver_walk: task %d prefix sum over %d envelopes produced a negative pair count (%lld)\n",
                   ThisTask, nb, (long long)batch_total);
            fflush(stdout);
            endrun(90001025);
            status = GX_RECEIVER_FAILED; break;
        }
        if(batch_total == 0) {continue;}

        Kokkos::deep_copy(UmHostI(counts_h.data(), (size_t)nb),  UmDevI(counts_d, (size_t)nb));
        Kokkos::deep_copy(UmHostI64(offsets_h.data(), (size_t)nb), UmDevI64(offsets_d, (size_t)nb));

        /* Emit in as many passes as the fixed output buffer needs.  Split points
         * are row boundaries chosen from the counts already in hand, so every
         * pass fits by construction and no batch is ever abandoned for being
         * dense. */
        int r0 = 0;
        while(r0 < nb && status == GX_RECEIVER_COMPLETED) {
            int r1 = r0; int64_t sub_total = 0;
            while(r1 < nb && sub_total + (int64_t)counts_h[(size_t)r1] <= pair_cap) {
                sub_total += (int64_t)counts_h[(size_t)r1];
                r1++;
            }
            if(r1 == r0) {
                /* The buffer holds num_local, and one envelope cannot admit a
                 * local particle twice -- the exported subtrees are disjoint --
                 * so a row this large means the traversal emitted a duplicate. */
                printf("gx_device_receiver_walk: task %d envelope %ld admits %d pairs against %d local particles; the traversal has emitted a duplicate\n",
                       ThisTask, base + r0, counts_h[(size_t)r0], num_local);
                fflush(stdout);
                endrun(90001026);
                status = GX_RECEIVER_FAILED; break;
            }
            const int64_t sub_base = offsets_h[(size_t)r0];
            {
                const struct gx_export_envelope_t *env_v = env_d;
                const struct gx_recv_leaf_t *leaf_v = leaf_d;
                int *scratch_v = scratch_d; int *counts_v = counts_d;
                int64_t *offsets_v = offsets_d; int *pairs_v = pairs_d; int *anom_v = anomaly_d;
                const int stride = GX_RECV_STRIDE;
                const int rr0 = r0;
                const int type_mask_trusted_v = TypePresenceMaskTrusted;   /* captured by value */
                Kokkos::parallel_for("gx_recv_compact", r1 - r0, KOKKOS_LAMBDA(int t) {
                    const int b = rr0 + t;
                    const int n = counts_v[b];
                    if(n <= 0) {return;}
                    const int64_t dst = offsets_v[b] - sub_base;
                    if(n <= stride) {
                        const size_t src = (size_t)b * stride;
                        for(int c = 0; c < n; c++) {pairs_v[dst + c] = scratch_v[src + c];}
                    } else {
                        /* Overflowed its scratch slot: re-walk straight into the
                         * final position. */
                        gx_recv_walk_one(env_v[b], supply_mask, type_mask_trusted_v,
                                         node_center, node_len, node_sibling,
                                         node_nextnode, node_bitflags, nextnode_aux,
                                         leaf_v,
                                         tree_base, tree_slots, node_capacity,
                                         foreign_base, pseudo_start, num_local,
                                         anom_v, &pairs_v[dst], n);
                    }
                });
                Kokkos::fence();
                gizmo_gpu_check_last_error("gx_recv_compact", r1 - r0);
            }

            Kokkos::deep_copy(UmHostI(pairs_h.data(), (size_t)sub_total),
                              UmDevI(pairs_d, (size_t)sub_total));
            /* Hand each row's candidates to the shared accept, which drifts them and
             * keeps the exact set.  Envelopes arrive grouped by peer in ascending
             * order, so a peer that spans a pass or a batch boundary simply carries
             * on in the next call. */
            long n_rows = 0;
            for(int b = r0; b < r1; b++) {
                const int t = envelope_peer[base + b];
                if(t < 0 || t >= NTask) {
                    /* Every envelope's sender is known; one that is not would have
                     * its candidates silently dropped, so stop rather than under-include. */
                    printf("gx_device_receiver_walk: task %d envelope %ld names sender %d, outside 0..%d\n",
                           ThisTask, base + b, t, NTask - 1);
                    fflush(stdout);
                    gizmo_request_controlled_stop(7737, "gx_device_receiver_walk: envelope with no valid sender",
                                                  __FILE__, __LINE__, __FUNCTION__);
                    status = GX_RECEIVER_FAILED;
                    break;
                }
                if(t == ThisTask) {continue;}
                const int64_t off = offsets_h[(size_t)b] - sub_base;
                struct gx_candidate_row row = {&env_h[(size_t)b], t, &pairs_h[(size_t)off], counts_h[(size_t)b]};
                rows_h[(size_t)n_rows++] = row;
            }
            if(status == GX_RECEIVER_COMPLETED &&
               gx_send_set_accept_rows(send_set, rows_h.data(), n_rows, search_mode,
                                       radius_policy, j_radius_scale, safety_factor) != 0) {
                status = GX_RECEIVER_FAILED;
            }
            r0 = r1;
        }
    }

    int anomaly = 0;
    {
        Kokkos::View<int, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ah(&anomaly);
        Kokkos::View<const int, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ad(anomaly_d);
        Kokkos::deep_copy(ah, ad);
    }

    Kokkos::kokkos_free<DevSp>(leaf_d);
    Kokkos::kokkos_free<DevSp>(env_d);
    Kokkos::kokkos_free<DevSp>(scratch_d);
    Kokkos::kokkos_free<DevSp>(counts_d);
    Kokkos::kokkos_free<DevSp>(offsets_d);
    Kokkos::kokkos_free<DevSp>(pairs_d);
    Kokkos::kokkos_free<DevSp>(anomaly_d);

    /* The host walk stops the run on this state, so reaching it here means the
     * same thing.  Falling back would only hide a malformed tree: the host walk
     * would meet it too. */
    if(anomaly) {
        printf("gx_device_receiver_walk: task %d walked into the index gap between the particle slots and the node base, or was handed an incompletely filled tree view; either way the walk cannot answer\n",
               ThisTask);
        fflush(stdout);
        endrun(90001024);
        return GX_RECEIVER_FAILED;
    }
    /* A pass that gave up reports FAILED rather than declining: what it already
     * handed to the send set cannot be taken back by a host rerun. */
    return status;
}


/* Per-TU init function: sets this TU's All_ptr to the shared UVM allocation */
