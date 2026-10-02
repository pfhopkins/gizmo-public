/* gpu_dispatch_templates.h — inline helpers for GPU kernel dispatch.
 *
 * gizmo_gpu_kernel_launch wraps the `Kokkos::parallel_for + Kokkos::fence +
 * gizmo_gpu_check_last_error` trio that appears after every Kokkos::parallel_for
 * in the GPU TUs.  Condenses 3 lines per call site to 1 and provides a single
 * place to add profiling / logging hooks later.
 *
 * Usage:
 *     gizmo_gpu_kernel_launch("density_kernel", N, KOKKOS_LAMBDA(int aa) {
 *         // per-particle work
 *     });
 *
 * For batched dispatches that want to log the batch range on error:
 *     gizmo_gpu_kernel_launch("cooling_loop", batch_n, KOKKOS_LAMBDA(int j) {
 *         do_the_cooling_for_particle(j, kp, kc);
 *     }, batch_start);
 *
 * Include AFTER <Kokkos_Core.hpp>.
 */
#ifndef GIZMO_GPU_DISPATCH_TEMPLATES_H
#define GIZMO_GPU_DISPATCH_TEMPLATES_H


#include <utility>
#include <type_traits>
#include "gpu_error_check.h"

template<typename F>
inline void gizmo_gpu_kernel_launch(const char *tag, int N, F&& f, int batch_start = -1)
{
    Kokkos::parallel_for(tag, N, std::forward<F>(f));
    Kokkos::fence();
    gizmo_gpu_check_last_error(tag, N, batch_start);
}


/* Team-policy sibling of gizmo_gpu_kernel_launch, for kernels whose per-item
 * work is itself parallel (a neighbor-list row, a node's children, a frequency
 * bin sweep). `league` teams of `team_size` lanes each; the functor receives the
 * team member, takes its item index from member.league_rank(), and divides the
 * item's work with TeamThreadRange.
 *
 * vector_length is fixed at 1, which is a portability requirement rather than a
 * preference: Kokkos caps vector_length at the backend warp size -- 32 on CUDA,
 * 64 on HIP -- and SILENTLY clamps an over-request, so a width carried on the
 * vector dimension cannot express a 64-lane division on NVIDIA. Carrying it on
 * the team dimension expresses both widths on both backends, and makes the
 * nested reduction single-level. It also keeps the item's memory stream
 * coalesced: with vector_length 1 the team rank is the fastest-varying thread
 * index, so lane t of a team reads element t of the item's contiguous run.
 *
 * gizmo_gpu_team_size_max reports the largest LAUNCHABLE team size for a
 * functor. It is an occupancy/register-derived legality bound and nothing more
 * -- in particular it does NOT shrink when the functor's reduction value grows,
 * so it must not be used as a memory or spill budget.
 */
template<typename F>
inline void gizmo_gpu_team_kernel_launch(const char *tag, int league, int team_size, F&& f)
{
    Kokkos::TeamPolicy<> policy(league, team_size, 1);
    Kokkos::parallel_for(tag, policy, std::forward<F>(f));
    Kokkos::fence();
    gizmo_gpu_check_last_error(tag, league, -1);
}

template<typename F>
inline int gizmo_gpu_team_size_max(const F& f)
{
    Kokkos::TeamPolicy<> probe(1, 1, 1);
    return probe.team_size_max(f, Kokkos::ParallelForTag());
}

/* How many work items the default execution space can have in flight at once:
 * the hardware's lane count, as Kokkos reports it for the backend in use.
 *
 * It answers one question -- are there more lanes than there are items to give
 * them -- which is what decides whether a loop should divide ONE item's work
 * among several lanes instead of giving each lane its own. Keyed on the machine
 * rather than on a tuned item count, so the same code makes the right choice on
 * a laptop's thread pool and on a GPU three orders of magnitude wider. */
inline int gizmo_gpu_lane_count()
{
    const int n = (int) Kokkos::DefaultExecutionSpace().concurrency();
    return (n > 0) ? n : 1;
}

/* True when the default execution space runs on the host, where a kernel's
 * items are already spread across the available threads and there are no SIMT
 * lanes left to divide an item's work among. Team-per-item policies are a
 * device idiom; on a host backend the flat launch above is the right shape. */
inline constexpr bool gizmo_gpu_default_space_is_host()
{
    return std::is_same<Kokkos::DefaultExecutionSpace::memory_space,
                        Kokkos::HostSpace>::value;
}


#endif /* GIZMO_GPU_DISPATCH_TEMPLATES_H */
