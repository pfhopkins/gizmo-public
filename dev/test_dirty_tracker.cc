/* test_dirty_tracker.cc — semantic unit tests for src/mesh/gpu_dirty_tracker.{h,cc}.
 *
 * Build standalone (no Kokkos / no MPI / no allvars dependency), from the repo root:
 *   c++ -std=c++17 -O0 -g dev/test_dirty_tracker.cc \
 *       src/mesh/gpu_dirty_tracker.cc -o /tmp/test_dirty_tracker
 *
 * Tests cover the following cases:
 *   1. overlapping cache ranges both receive marks
 *   2. consuming one cache does not clear the other
 *   3. duplicate indices dedup
 *   4. out-of-range marks ignored
 *   5. range boundary start/end correct
 *   6. unregister/re-register does not leave stale handles
 *   7. threshold promotion behavior is per-cache, not global
 *
 * Each test prints "PASS" or "FAIL: ..."; exit code is non-zero if any failed. */

#include "../src/mesh/gpu_dirty_tracker.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static int g_fail = 0;

#define ASSERT_EQ(label, expected, actual) do { \
    long long _e = (long long)(expected), _a = (long long)(actual); \
    if(_e != _a) { fprintf(stderr, "FAIL %s: expected %lld, got %lld\n", label, _e, _a); g_fail++; } \
} while(0)

#define ASSERT_TRUE(label, cond) do { \
    if(!(cond)) { fprintf(stderr, "FAIL %s: condition false\n", label); g_fail++; } \
} while(0)

/* Helper: collect all consumed indices into a vector via consume callback. */
static void collect_cb(int j, void *ud) {
    ((std::vector<int> *)ud)->push_back(j);
}

static void test_overlap_two_caches() {
    fprintf(stderr, "test_overlap_two_caches: ");
    gpu_dirty_handle_t h1 = gpu_dirty_tracker_register(0, 100, 0);
    gpu_dirty_handle_t h2 = gpu_dirty_tracker_register(50, 100, 0); /* covers [50, 150) */
    ASSERT_TRUE("h1 valid", h1 >= 0);
    ASSERT_TRUE("h2 valid", h2 >= 0);

    /* Fresh registrations start all_dirty=1 — consume to clear. */
    std::vector<int> v;
    gpu_dirty_tracker_consume(h1, collect_cb, &v); v.clear();
    gpu_dirty_tracker_consume(h2, collect_cb, &v); v.clear();

    int idx[1] = { 75 };
    gpu_dirty_tracker_mark_indices(idx, 1);
    ASSERT_EQ("h1 popcount after mark j=75", 1, gpu_dirty_tracker_popcount(h1));
    ASSERT_EQ("h2 popcount after mark j=75", 1, gpu_dirty_tracker_popcount(h2));

    /* Consume h1; h2 must still have the bit. */
    v.clear(); gpu_dirty_tracker_consume(h1, collect_cb, &v);
    ASSERT_EQ("h1 consumed count", 1, (int)v.size());
    if(v.size() == 1) ASSERT_EQ("h1 consumed value", 75, v[0]);
    ASSERT_EQ("h1 popcount after consume", 0, gpu_dirty_tracker_popcount(h1));
    ASSERT_EQ("h2 popcount after h1 consume (should be untouched)", 1, gpu_dirty_tracker_popcount(h2));

    v.clear(); gpu_dirty_tracker_consume(h2, collect_cb, &v);
    ASSERT_EQ("h2 consumed count", 1, (int)v.size());
    if(v.size() == 1) ASSERT_EQ("h2 consumed value", 75, v[0]);

    gpu_dirty_tracker_unregister(h1);
    gpu_dirty_tracker_unregister(h2);
    fprintf(stderr, "done\n");
}

static void test_dedup() {
    fprintf(stderr, "test_dedup: ");
    gpu_dirty_handle_t h = gpu_dirty_tracker_register(0, 100, 0);
    std::vector<int> v;
    gpu_dirty_tracker_consume(h, collect_cb, &v); v.clear();

    int idx[5] = { 42, 42, 42, 42, 42 };
    gpu_dirty_tracker_mark_indices(idx, 5);
    ASSERT_EQ("popcount with 5 dup marks of j=42", 1, gpu_dirty_tracker_popcount(h));
    gpu_dirty_tracker_unregister(h);
    fprintf(stderr, "done\n");
}

static void test_out_of_range_ignored() {
    fprintf(stderr, "test_out_of_range_ignored: ");
    gpu_dirty_handle_t h = gpu_dirty_tracker_register(0, 100, 0); /* covers [0, 100) */
    std::vector<int> v;
    gpu_dirty_tracker_consume(h, collect_cb, &v); v.clear();

    int idx[3] = { 200, -1, 99 };
    gpu_dirty_tracker_mark_indices(idx, 3);
    ASSERT_EQ("popcount with 200, -1, 99: only 99 covered", 1, gpu_dirty_tracker_popcount(h));
    gpu_dirty_tracker_unregister(h);
    fprintf(stderr, "done\n");
}

static void test_range_boundaries() {
    fprintf(stderr, "test_range_boundaries: ");
    gpu_dirty_handle_t h = gpu_dirty_tracker_register(50, 50, 0); /* covers [50, 100) */
    std::vector<int> v;
    gpu_dirty_tracker_consume(h, collect_cb, &v); v.clear();

    /* Inclusive 50 to exclusive 100. */
    gpu_dirty_tracker_mark_range(0, 100); /* should mark only 50..100 in this cache */
    ASSERT_EQ("popcount marking [0,100) into cache [50,100)", 50, gpu_dirty_tracker_popcount(h));

    v.clear(); gpu_dirty_tracker_consume(h, collect_cb, &v);
    ASSERT_EQ("consumed count = 50", 50, (int)v.size());
    if(!v.empty()) {
        ASSERT_EQ("first consumed value", 50, v[0]);
        ASSERT_EQ("last consumed value", 99, v[(int)v.size() - 1]);
    }

    gpu_dirty_tracker_unregister(h);
    fprintf(stderr, "done\n");
}


static void test_unregister_re_register() {
    fprintf(stderr, "test_unregister_re_register: ");
    gpu_dirty_handle_t h1 = gpu_dirty_tracker_register(0, 100, 0);
    int idx[1] = { 50 };
    std::vector<int> v;
    gpu_dirty_tracker_consume(h1, collect_cb, &v); v.clear();

    gpu_dirty_tracker_mark_indices(idx, 1);
    ASSERT_EQ("h1 marked", 1, gpu_dirty_tracker_popcount(h1));

    gpu_dirty_tracker_unregister(h1);

    /* Mark again — no live caches. Should be a clean no-op. */
    gpu_dirty_tracker_mark_indices(idx, 1);

    /* Re-register on same range. */
    gpu_dirty_handle_t h2 = gpu_dirty_tracker_register(0, 100, 0);
    ASSERT_TRUE("re-register valid", h2 >= 0);
    /* Fresh registration: starts all_dirty (newborn cache); not stale from h1. */
    v.clear(); gpu_dirty_tracker_consume(h2, collect_cb, &v);
    ASSERT_EQ("re-registered consume size = count (newborn all_dirty)", 100, (int)v.size());
    ASSERT_TRUE("post-consume not all_dirty", !gpu_dirty_tracker_is_all_dirty(h2));
    ASSERT_EQ("post-consume popcount", 0, gpu_dirty_tracker_popcount(h2));

    gpu_dirty_tracker_unregister(h2);
    fprintf(stderr, "done\n");
}

static void test_promote_per_cache() {
    fprintf(stderr, "test_promote_per_cache: ");
    /* Need range > 1M (G_DIRTY_PROMOTE_THRESHOLD) to trigger promote. */
    const int N = 2 * (1 << 20); /* 2M */
    gpu_dirty_handle_t hbig   = gpu_dirty_tracker_register(0, N, 0);
    gpu_dirty_handle_t hsmall = gpu_dirty_tracker_register(0, 100, 0);
    std::vector<int> v;
    gpu_dirty_tracker_consume(hbig, collect_cb, &v); v.clear();
    gpu_dirty_tracker_consume(hsmall, collect_cb, &v); v.clear();

    /* mark_range(0, N) covers all of hbig's [0,N) and only [0,100) of hsmall. */
    gpu_dirty_tracker_mark_range(0, N);
    ASSERT_TRUE("hbig promoted to all_dirty", gpu_dirty_tracker_is_all_dirty(hbig));
    ASSERT_TRUE("hsmall NOT promoted (only 100 marks)", !gpu_dirty_tracker_is_all_dirty(hsmall));
    ASSERT_EQ("hsmall popcount = 100", 100, gpu_dirty_tracker_popcount(hsmall));

    gpu_dirty_tracker_unregister(hbig);
    gpu_dirty_tracker_unregister(hsmall);
    fprintf(stderr, "done\n");
}

int main() {
    test_overlap_two_caches();
    test_dedup();
    test_out_of_range_ignored();
    test_range_boundaries();
    test_unregister_re_register();
    test_promote_per_cache();
    if(g_fail) {
        fprintf(stderr, "\n=== FAIL: %d assertion(s) ===\n", g_fail);
        return 1;
    }
    fprintf(stderr, "\n=== ALL PASS ===\n");
    return 0;
}
