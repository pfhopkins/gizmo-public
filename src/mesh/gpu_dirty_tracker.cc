/* gpu_dirty_tracker.cc — see gpu_dirty_tracker.h for design. */

#include "gpu_dirty_tracker.h"

#include <cstdlib>
#include <cstring>
#include <cstdint>

/* Promote-to-all threshold per cache (matches the old G_DIRTY_PROMOTE_THRESHOLD
 * in gpu_neighbor_list.cc). Above this many unique dirty indices, narrow
 * refresh costs ~the same as full, plus we'd be paying the staging cost for
 * the index list. Per-cache, not global. */
static const int G_DIRTY_PROMOTE_THRESHOLD = 1 << 20; /* 1M */

/* Bitset word type. uint64_t lets us iterate set bits via __builtin_ctzll. */
typedef uint64_t bword_t;
static const int BWORD_BITS = 64;

struct dirty_cache_t {
    int valid;          /* 0 = slot free */
    int base;           /* particle index base */
    int count;          /* particle index count (bitset size = count bits) */
    bword_t *bits;      /* (count + 63) / 64 words */
    int all_dirty;      /* if set, refresh all [base, base+count) */
    int popcount_hint;  /* approximate popcount; precise on demand */
};

/* Small fixed registry. Today: gas + alltypes => 2 caches; future home/ghost
 * split => up to 4. 8 is comfortable headroom. */
#define GPU_DIRTY_MAX_CACHES 8
static struct dirty_cache_t g_caches[GPU_DIRTY_MAX_CACHES];

static inline int n_words_(int count)
{
    return (count + BWORD_BITS - 1) / BWORD_BITS;
}

static inline void clear_cache_state_(struct dirty_cache_t *c)
{
    if(c->bits && c->count > 0) {
        memset(c->bits, 0, (size_t)n_words_(c->count) * sizeof(bword_t));
    }
    c->all_dirty = 0;
    c->popcount_hint = 0;
}

gpu_dirty_handle_t gpu_dirty_tracker_register(int base, int count, int start_clean)
{
    if(count <= 0) return -1;
    for(int h = 0; h < GPU_DIRTY_MAX_CACHES; h++) {
        if(!g_caches[h].valid) {
            g_caches[h].base = base;
            g_caches[h].count = count;
            g_caches[h].bits = (bword_t *) malloc((size_t)n_words_(count) * sizeof(bword_t));
            memset(g_caches[h].bits, 0, (size_t)n_words_(count) * sizeof(bword_t));
            /* A caller that has just written every row from the live P[] passes
             * start_clean=1, so the first consume refreshes nothing. Anyone else
             * starts all-dirty, which refreshes the whole range once. */
            g_caches[h].all_dirty = start_clean ? 0 : 1;
            g_caches[h].popcount_hint = 0;
            g_caches[h].valid = 1;
            return h;
        }
    }
    return -1; /* registry full */
}

void gpu_dirty_tracker_unregister(gpu_dirty_handle_t handle)
{
    if(handle < 0 || handle >= GPU_DIRTY_MAX_CACHES) return;
    if(!g_caches[handle].valid) return;
    if(g_caches[handle].bits) { free(g_caches[handle].bits); g_caches[handle].bits = NULL; }
    g_caches[handle].valid = 0;
    g_caches[handle].count = 0;
    g_caches[handle].base = 0;
    g_caches[handle].all_dirty = 0;
    g_caches[handle].popcount_hint = 0;
}

static inline void cache_set_bit_(struct dirty_cache_t *c, int local_idx)
{
    bword_t *w = &c->bits[local_idx / BWORD_BITS];
    bword_t mask = ((bword_t)1) << (local_idx % BWORD_BITS);
    if((*w & mask) == 0) {
        *w |= mask;
        c->popcount_hint++;
        if(c->popcount_hint > G_DIRTY_PROMOTE_THRESHOLD) {
            c->all_dirty = 1;
            /* leave bits set; consume() ignores them when all_dirty */
        }
    }
}

void gpu_dirty_tracker_mark_indices(const int *indices, int n)
{
    if(!indices || n <= 0) return;
    for(int h = 0; h < GPU_DIRTY_MAX_CACHES; h++) {
        struct dirty_cache_t *c = &g_caches[h];
        if(!c->valid) continue;
        if(c->all_dirty) continue; /* already covered */
        int base = c->base, count = c->count;
        for(int k = 0; k < n; k++) {
            int j = indices[k];
            int local = j - base;
            if(local < 0 || local >= count) continue;
            cache_set_bit_(c, local);
            if(c->all_dirty) break; /* promoted mid-loop */
        }
    }
}

void gpu_dirty_tracker_mark_range(int start, int end)
{
    if(end <= start) return;
    for(int h = 0; h < GPU_DIRTY_MAX_CACHES; h++) {
        struct dirty_cache_t *c = &g_caches[h];
        if(!c->valid) continue;
        if(c->all_dirty) continue;
        int base = c->base, count = c->count;
        int lo = start - base; if(lo < 0) lo = 0;
        int hi = end   - base; if(hi > count) hi = count;
        if(hi <= lo) continue;
        int n = hi - lo;
        if(c->popcount_hint + n > G_DIRTY_PROMOTE_THRESHOLD) {
            c->all_dirty = 1;
            continue;
        }
        /* Whole-word write where possible, edge bits handled per-bit. */
        int word_lo = lo / BWORD_BITS;
        int word_hi = (hi - 1) / BWORD_BITS;
        if(word_lo == word_hi) {
            for(int b = lo; b < hi; b++) cache_set_bit_(c, b);
        } else {
            for(int b = lo; b < (word_lo + 1) * BWORD_BITS && b < hi; b++) cache_set_bit_(c, b);
            for(int w = word_lo + 1; w < word_hi; w++) {
                bword_t prev = c->bits[w];
                c->bits[w] = (bword_t)~(bword_t)0;
                /* popcount_hint update: only the bits that flipped 0->1 */
                int delta = BWORD_BITS - __builtin_popcountll(prev);
                c->popcount_hint += delta;
                if(c->popcount_hint > G_DIRTY_PROMOTE_THRESHOLD) {
                    c->all_dirty = 1;
                    break;
                }
            }
            if(!c->all_dirty) {
                for(int b = word_hi * BWORD_BITS; b < hi; b++) cache_set_bit_(c, b);
            }
        }
    }
}

void gpu_dirty_tracker_consume(gpu_dirty_handle_t handle,
                               void (*callback)(int j, void *userdata),
                               void *userdata)
{
    if(handle < 0 || handle >= GPU_DIRTY_MAX_CACHES) return;
    struct dirty_cache_t *c = &g_caches[handle];
    if(!c->valid || !callback) return;

    if(c->all_dirty) {
        int end = c->base + c->count;
        for(int j = c->base; j < end; j++) callback(j, userdata);
    } else {
        int nw = n_words_(c->count);
        for(int w = 0; w < nw; w++) {
            bword_t bits = c->bits[w];
            while(bits) {
                int bit = __builtin_ctzll(bits);
                bits &= bits - 1;
                int local = w * BWORD_BITS + bit;
                if(local < c->count) callback(c->base + local, userdata);
            }
        }
    }
    clear_cache_state_(c);
}

int gpu_dirty_tracker_popcount(gpu_dirty_handle_t handle)
{
    if(handle < 0 || handle >= GPU_DIRTY_MAX_CACHES) return 0;
    struct dirty_cache_t *c = &g_caches[handle];
    if(!c->valid) return 0;
    if(c->all_dirty) return c->count;
    return c->popcount_hint;
}

int gpu_dirty_tracker_is_all_dirty(gpu_dirty_handle_t handle)
{
    if(handle < 0 || handle >= GPU_DIRTY_MAX_CACHES) return 0;
    if(!g_caches[handle].valid) return 0;
    return g_caches[handle].all_dirty;
}
