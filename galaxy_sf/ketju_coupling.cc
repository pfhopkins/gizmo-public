/* KETJU regularized integrator coupling for GIZMO
 *
 * Implements the coupling between GIZMO's gravity tree and the MSTAR
 * regularized N-body integrator. Based on the GADGET4-KETJU architecture
 * (Mannerkoski+ 2022) but adapted for GIZMO's data structures and extended
 * to keep particles individually visible in the tree for radiation physics.
 *
 * Architecture:
 *   1. limit_timesteps(): force chain particles to shared timebin
 *   2. find_regions(): detect chain regions around massive particles
 *   3. run_integration(): collect particles, subtract tree force, run MSTAR
 *   4. set_final_velocities(): swap in true velocities after drift
 *   5. finish_step(): clear region data at end of step
 */

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <vector>
#include <set>
#include <unordered_map>
#include <queue>
#include <utility>
#include <cfloat>
#include <queue>
#include <functional>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "ketju_coupling.h"

#ifdef KETJU_REGULARIZATION

/* PN classification for KETJU integrator (0 = PN/BH, 1 = non-PN/star, -1 = not Type 5).
 * Derived from ProtoStellarStage: stage 7 (relic/BH) → 0, stage <7 (star) → 1. */
#define KETJU_PN_CLASS(i) ((P[(i)].ProtoStellarStage == 7) ? 0 : 1)

/* Cost tracking for load balancing (keyed by first center ID per region) */
static std::unordered_map<MyIDType, double> region_previous_cost;

/* ============================================================
 *  Lightweight particle data for MPI communication
 * ============================================================ */
struct ketju_mpi_particle {
    MyIDType ID;
    int Type;
    int Task;
    int Index;
    double Mass;
    double Pos[3];
    double Vel[3];
#ifdef SINGLE_STAR_SINK_DYNAMICS
    short int KetjuPNClass;
    double Spin[3];
#endif
};

/* ============================================================
 *  Per-particle extra data tracked alongside integrator arrays
 * ============================================================ */
struct ketju_extra_data {
    MyIDType ID;
    int Task;
    int Index;
    int Type;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    short int KetjuPNClass;
#endif
};

/* ============================================================
 *  MPI task group for parallel KETJU computation
 *  Based on public GADGET4-KETJU (Mannerkoski+ 2023)
 * ============================================================ */
struct KetjuTaskGroup {
    MPI_Group group;
    MPI_Comm comm;
    int rank;      /* rank within this group, MPI_UNDEFINED if not a member */
    int size;
    int root;      /* rank of root within the group */
    int root_sim;  /* rank of root in MPI_COMM_WORLD */

    KetjuTaskGroup() : group(MPI_GROUP_NULL), comm(MPI_COMM_NULL), rank(MPI_UNDEFINED), size(0), root(MPI_UNDEFINED), root_sim(MPI_UNDEFINED) {}

    void init(const std::vector<int> &task_sim_indices, int tag = 0) {
        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group_incl(world_group, task_sim_indices.size(), task_sim_indices.data(), &group);
        MPI_Group_size(group, &size);
        MPI_Group_rank(group, &rank);
        MPI_Comm_create_group(MPI_COMM_WORLD, group, tag, &comm);
        MPI_Group_free(&world_group);
        root = 0;
        root_sim = task_sim_indices[0];
    }

    void free_comms() {
        if(comm != MPI_COMM_NULL) { MPI_Comm_free(&comm); comm = MPI_COMM_NULL; }
        if(group != MPI_GROUP_NULL) { MPI_Group_free(&group); group = MPI_GROUP_NULL; }
        rank = MPI_UNDEFINED; size = 0; root = MPI_UNDEFINED; root_sim = MPI_UNDEFINED;
    }

    /* Set both groups to have the same root task (from intersection) */
    void set_common_root(KetjuTaskGroup &other) {
        if(group == MPI_GROUP_NULL || other.group == MPI_GROUP_NULL) return;
        MPI_Group world_group, intersection;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group_intersection(group, other.group, &intersection);
        if(intersection != MPI_GROUP_EMPTY) {
            int root_inter = 0;
            MPI_Group_translate_ranks(intersection, 1, &root_inter, group, &root);
            MPI_Group_translate_ranks(intersection, 1, &root_inter, other.group, &other.root);
            MPI_Group_translate_ranks(group, 1, &root, world_group, &root_sim);
            MPI_Group_translate_ranks(other.group, 1, &other.root, world_group, &other.root_sim);
        }
        MPI_Group_free(&intersection);
        MPI_Group_free(&world_group);
    }

    int is_member() const { return rank != MPI_UNDEFINED; }
    int is_root() const { return rank == root; }
};

/* Compute info for scheduling regions across tasks */
struct region_compute_info {
    int compute_sequence_position; /* which sequential run this region belongs to */
    int first_task_index;          /* first MPI task in compute group */
    int final_task_index;          /* last MPI task in compute group */
};

/* ============================================================
 *  A single KETJU region
 * ============================================================ */
struct KetjuRegion {
    std::vector<ketju_mpi_particle> centers;      /* BH/massive-star region centers */
    std::vector<ketju_mpi_particle> all_particles; /* all chain members (gathered on root) */
    std::set<int> local_member_indices;            /* local P[] indices in this region */
    int total_particle_count;
    int num_pn_particles;                          /* count of PN-enabled (SMBH) particles */

    /* MPI task groups for parallel computation */
    KetjuTaskGroup affected_tasks;  /* tasks that hold particles in this region */
    KetjuTaskGroup compute_tasks;   /* tasks that run the integrator */
    region_compute_info compute_info;
    std::vector<int> affected_sim_indices;       /* world ranks of affected tasks */
    std::vector<int> particle_counts_on_affected; /* particle count per affected task */

    /* integrator state (only allocated on compute root) */
    struct ketju_system *integrator;
    std::vector<ketju_extra_data> extra_data;

    /* CoM frame data (stored at setup, used in scatter) */
    double com_pos[3];
    double com_vel[3];
    integertime ti_step;     /* integer timestep for this region */

    KetjuRegion() : total_particle_count(0), num_pn_particles(0), integrator(NULL), ti_step(0) {
        for(int j = 0; j < 3; j++) { com_pos[j] = 0; com_vel[j] = 0; }
        compute_info.compute_sequence_position = 0;
        compute_info.first_task_index = 0;
        compute_info.final_task_index = 0;
    }
    ~KetjuRegion() {
        if(integrator) { ketju_free_system(integrator); free(integrator); }
        affected_tasks.free_comms();
        compute_tasks.free_comms();
    }
};

/* ============================================================
 *  Loop scheduling for parallel N² force computation
 *  Distributes the upper triangle of the N×N pair loop across
 *  num_proc tasks. Returns the start index for edge_index-th block.
 *  Based on public GADGET4-KETJU (Mannerkoski+ 2023).
 * ============================================================ */
static int loop_scheduling_block_edge(int Nloop, int num_proc, int edge_index)
{
    if(Nloop <= 0 || num_proc <= 0) return 0;
    if(num_proc >= Nloop) num_proc = Nloop;
    if(edge_index <= 0) return 0;
    if(edge_index >= num_proc) return Nloop - 1;
    double P = (double)num_proc, N = (double)Nloop;
    return (int)ceil(N - 0.5 - sqrt(N * (N - 1) * (P - edge_index) / P + 0.25));
}

/* ============================================================
 *  Module-level state (persists within a timestep)
 * ============================================================ */
static std::vector<KetjuRegion> ActiveRegions;
static std::set<int> AllKetjuParticleIndices; /* local indices of all particles in any region */

/* ============================================================
 *  Region persistence: cached regions across steps
 *  Communicators are preserved when membership is unchanged.
 * ============================================================ */
struct CachedRegionKey {
    MyIDType ID;
    int Task;
    int Index; /* local P[] index — stable between domain decompositions */
    bool operator==(const CachedRegionKey &o) const { return ID == o.ID && Task == o.Task && Index == o.Index; }
    bool operator<(const CachedRegionKey &o) const {
        if(ID != o.ID) return ID < o.ID;
        if(Task != o.Task) return Task < o.Task;
        return Index < o.Index;
    }
};

struct CachedRegionInfo {
    std::vector<CachedRegionKey> sorted_particle_keys; /* sorted (ID,Task) for staleness check */
    KetjuTaskGroup affected_tasks;
    KetjuTaskGroup compute_tasks;
    std::vector<int> affected_sim_indices;
    std::vector<int> particle_counts_on_affected;
    region_compute_info compute_info;
    int total_particle_count;

    CachedRegionInfo() : total_particle_count(0) {}
    ~CachedRegionInfo() {
        affected_tasks.free_comms();
        compute_tasks.free_comms();
    }
    /* move-only to protect MPI handles */
    CachedRegionInfo(CachedRegionInfo &&o) noexcept
        : sorted_particle_keys(std::move(o.sorted_particle_keys)),
          affected_sim_indices(std::move(o.affected_sim_indices)),
          particle_counts_on_affected(std::move(o.particle_counts_on_affected)),
          compute_info(o.compute_info),
          total_particle_count(o.total_particle_count)
    {
        affected_tasks = o.affected_tasks; o.affected_tasks.comm = MPI_COMM_NULL; o.affected_tasks.group = MPI_GROUP_NULL;
        compute_tasks = o.compute_tasks; o.compute_tasks.comm = MPI_COMM_NULL; o.compute_tasks.group = MPI_GROUP_NULL;
    }
    CachedRegionInfo &operator=(CachedRegionInfo &&o) noexcept {
        if(this != &o) {
            affected_tasks.free_comms(); compute_tasks.free_comms();
            sorted_particle_keys = std::move(o.sorted_particle_keys);
            affected_sim_indices = std::move(o.affected_sim_indices);
            particle_counts_on_affected = std::move(o.particle_counts_on_affected);
            compute_info = o.compute_info; total_particle_count = o.total_particle_count;
            affected_tasks = o.affected_tasks; o.affected_tasks.comm = MPI_COMM_NULL; o.affected_tasks.group = MPI_GROUP_NULL;
            compute_tasks = o.compute_tasks; o.compute_tasks.comm = MPI_COMM_NULL; o.compute_tasks.group = MPI_GROUP_NULL;
        }
        return *this;
    }
    CachedRegionInfo(const CachedRegionInfo &) = delete;
    CachedRegionInfo &operator=(const CachedRegionInfo &) = delete;
};

static std::vector<CachedRegionInfo> CachedRegions;
static int KetjuRegionsStale = 1; /* 1 = force rebuild (after domain decomp or first call) */

/* Cache for gather_chain_centers result — shared between ketju_limit_timesteps
 * and ketju_find_regions to avoid a redundant MPI_Allgatherv on COMM_WORLD */
static std::vector<ketju_mpi_particle> CachedChainCenters;
static int CachedChainCentersValid = 0;

/* ============================================================
 *  Cost-based scheduling (Phase B)
 *  Based on public GADGET4-KETJU (Mannerkoski+ 2023).
 *  Estimates per-region cost and distributes compute tasks.
 * ============================================================ */
static double region_cost_estimate(int region_index)
{
    if(ActiveRegions[region_index].centers.empty()) return 1.0;
    MyIDType key = ActiveRegions[region_index].centers[0].ID;
    auto it = region_previous_cost.find(key);
    if(it != region_previous_cost.end()) return it->second;
    /* no previous cost — estimate from N² */
    double N = ActiveRegions[region_index].total_particle_count;
    return 20.0 * N * N;
}

/* Number of sequential rounds chosen by allocate_compute_tasks_for_regions.
 * Each region is assigned a round index (compute_sequence_position) and
 * ketju_run_integration iterates over rounds in outer loop, so regions in
 * the same round run in parallel while tasks move to the next round only
 * after all regions in the current round finish. Matches the multi-round
 * scheduler used by GADGET4-KETJU / gadget_tnt_new for better load balance. */
static int g_ketju_num_sequential_runs = 0;

/* Bucket regions into `num_runs` rounds using a longest-processing-time (LPT)
 * heuristic: repeatedly place the most-expensive unassigned region into the
 * least-loaded bucket. Each bucket is capped at NTask regions so every
 * region gets at least one task. */
static std::vector<std::vector<int>> lpt_bucket_regions(int num_runs)
{
    int n_regions = (int)ActiveRegions.size();
    typedef std::pair<double, int> cost_pair;   /* (cost, index) */

    /* max-heap over regions by cost */
    std::priority_queue<cost_pair> region_pq;
    for(int r = 0; r < n_regions; r++) region_pq.push({region_cost_estimate(r), r});

    /* min-heap over buckets by accumulated cost */
    std::priority_queue<cost_pair, std::vector<cost_pair>, std::greater<cost_pair>> bucket_pq;
    for(int b = 0; b < num_runs; b++) bucket_pq.push({0.0, b});

    std::vector<std::vector<int>> buckets(num_runs);
    while(!region_pq.empty() && !bucket_pq.empty()) {
        double rc = region_pq.top().first;
        int    r  = region_pq.top().second;
        region_pq.pop();
        double bc = bucket_pq.top().first;
        int    b  = bucket_pq.top().second;
        bucket_pq.pop();
        buckets[b].push_back(r);
        /* only push this bucket back if it can still fit more regions */
        if((int)buckets[b].size() < NTask)
            bucket_pq.push({bc + rc, b});
    }
    return buckets;
}

/* Estimate wall-time for a given bucketing: sum over buckets of the
 * max (region_cost / tasks_allocated_to_region). Scaling is cut off at
 * N_particles / particles_per_task_scaling_limit because the integrator
 * stops getting faster past that point. */
static double estimate_total_time_for(const std::vector<std::vector<int>>& buckets)
{
    const int ppts_limit = 7; /* see Mannerkoski+2022 scaling data */
    double total = 0;
    for(const auto& bucket : buckets) {
        if(bucket.empty()) continue;
        double bucket_cost = 0;
        for(int r : bucket) bucket_cost += region_cost_estimate(r);
        if(bucket_cost <= 0) bucket_cost = 1.0;

        double bucket_max = 0;
        for(int r : bucket) {
            double frac = region_cost_estimate(r) / bucket_cost;
            int n_tasks = DMAX(1, (int)(frac * NTask + 0.5));
            int npart = ActiveRegions[r].total_particle_count;
            int max_useful = DMAX(1, npart / ppts_limit);
            if(n_tasks > max_useful) n_tasks = max_useful;
            double t_reg = region_cost_estimate(r) / (double)n_tasks;
            if(t_reg > bucket_max) bucket_max = t_reg;
        }
        total += bucket_max;
    }
    return total;
}

/* Assign MPI task ranges to the regions in one bucket, proportional to cost. */
static void assign_tasks_for_bucket(const std::vector<int>& bucket_regions, int seq_pos)
{
    if(bucket_regions.empty()) return;
    double bucket_cost = 0;
    for(int r : bucket_regions) bucket_cost += region_cost_estimate(r);
    if(bucket_cost <= 0) bucket_cost = 1.0;

    int min_per_task = (int)DMAX(All.KetjuMinParticlesPerTask, 2);
    int tasks_used = 0;
    for(size_t i = 0; i < bucket_regions.size(); i++) {
        int r = bucket_regions[i];
        double frac = region_cost_estimate(r) / bucket_cost;
        int n_tasks = DMAX(1, (int)(frac * NTask + 0.5));

        /* respect minimum particles per task */
        if(ActiveRegions[r].total_particle_count > 0) {
            int max_tasks = ActiveRegions[r].total_particle_count / min_per_task;
            if(max_tasks < 1) max_tasks = 1;
            if(n_tasks > max_tasks) n_tasks = max_tasks;
        }

        /* don't exceed available tasks in this bucket. If the bucket is
         * already full (tasks_used >= NTask) we still need to assign this
         * region somewhere valid — co-locate it on the last task. Without
         * this clamp the indices go out of bounds (e.g. range=[4,4] on
         * NTask=4) and MPI_Group_incl/MPI_Comm_create_group return a
         * corrupt group → segfault in MPI_Group_intersection later. */
        if(tasks_used >= NTask) {
            ActiveRegions[r].compute_info.first_task_index          = NTask - 1;
            ActiveRegions[r].compute_info.final_task_index          = NTask - 1;
            ActiveRegions[r].compute_info.compute_sequence_position = seq_pos;
            continue;
        }
        if(tasks_used + n_tasks > NTask) n_tasks = NTask - tasks_used;
        if(n_tasks < 1) n_tasks = 1;

        ActiveRegions[r].compute_info.first_task_index          = tasks_used;
        ActiveRegions[r].compute_info.final_task_index          = tasks_used + n_tasks - 1;
        ActiveRegions[r].compute_info.compute_sequence_position = seq_pos;
        tasks_used += n_tasks;
    }
}

static void allocate_compute_tasks_for_regions(void)
{
    int n_regions = (int)ActiveRegions.size();
    if(n_regions == 0) { g_ketju_num_sequential_runs = 0; return; }

    /* Sweep num_runs from the minimum (limited by NTask per bucket) upward
     * and pick the one with the smallest estimated total wall-time. We
     * exit early when adding another round stops improving the estimate
     * (the cost surface is empirically monotone past the optimum). */
    int min_runs = (n_regions + NTask - 1) / NTask;
    if(min_runs < 1) min_runs = 1;

    double best_time = DBL_MAX;
    std::vector<std::vector<int>> best_buckets;
    for(int num_runs = min_runs; num_runs <= n_regions; num_runs++) {
        auto buckets = lpt_bucket_regions(num_runs);
        double t = estimate_total_time_for(buckets);
        if(t < best_time) {
            best_time   = t;
            best_buckets = std::move(buckets);
        } else {
            break;
        }
    }

    for(size_t s = 0; s < best_buckets.size(); s++) {
        assign_tasks_for_bucket(best_buckets[s], (int)s);
    }
    g_ketju_num_sequential_runs = (int)best_buckets.size();
}

static void update_region_costs(void)
{
    /* After integration, store actual cost for next step. The cost is only
     * known on each region's compute_root (it has the integrator), but the
     * map MUST be globally consistent — otherwise next step's
     * region_cost_estimate / allocate_compute_tasks_for_regions returns
     * different compute_info.first/final_task_index on different tasks,
     * which then produces different compute groups and deadlocks
     * MPI_Comm_create_group inside compute_tasks.init.
     *
     * Pack one cost per region into a dense array indexed by region, then
     * MPI_Allreduce(MAX). Non-roots contribute 0 → the compute_root's value
     * wins. */
    int n_regions = (int)ActiveRegions.size();
    if(n_regions == 0) return;

    std::vector<double> local_costs(n_regions, 0.0);
    for(int r = 0; r < n_regions; r++) {
        if(!ActiveRegions[r].compute_tasks.is_root()) continue;
        if(!ActiveRegions[r].integrator) continue;
        if(ActiveRegions[r].centers.empty()) continue;

        local_costs[r] = (double)(ActiveRegions[r].integrator->perf->successful_steps +
                                  ActiveRegions[r].integrator->perf->failed_steps) *
                         ActiveRegions[r].total_particle_count * ActiveRegions[r].total_particle_count;
    }

    std::vector<double> global_costs(n_regions, 0.0);
    MPI_Allreduce(local_costs.data(), global_costs.data(), n_regions, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    for(int r = 0; r < n_regions; r++) {
        if(ActiveRegions[r].centers.empty()) continue;
        if(global_costs[r] <= 0) continue;
        MyIDType key = ActiveRegions[r].centers[0].ID;
        region_previous_cost[key] = global_costs[r];
    }
}

/* ============================================================
 *  Helper: is this particle a valid chain center?
 * ============================================================ */
static int is_chain_center(int i)
{
    /* STARFORGE: all sinks are Type 5; ProtoStellarStage selects role.
     * Stages 0-6 (proto / MS / pre-collapse) use KetjuMinStarMass.
     * Stage 7 (relic/BH) uses KetjuMinBHMass. */
    if(P[i].Type == 5 && P[i].ProtoStellarStage <  7 && All.KetjuMinStarMass > 0) {
        if(P[i].Mass >= All.KetjuMinStarMass) return 1;
    }
    if(P[i].Type == 5 && P[i].ProtoStellarStage == 7 && All.KetjuMinBHMass > 0) {
        if(P[i].Mass >= All.KetjuMinBHMass) return 1;
    }
    return 0;
}

/* ============================================================
 *  Helper: is this particle eligible for chain membership?
 * ============================================================ */
static int is_chain_eligible(int i)
{
    /* STARFORGE: every Type 5 sink is a chain-member candidate (regardless of stage) */
    if(P[i].Type == 5) return 1;
    return 0;
}

/* ============================================================
 *  Helper: is this particle PN-enabled? (BH only)
 * ============================================================ */
static int is_pn_particle(ketju_mpi_particle *p)
{
#ifdef SINGLE_STAR_SINK_DYNAMICS
    if(p->Type == 5 && p->KetjuPNClass == 0) return 1; /* relic BH: PN */
#ifdef KETJU_PN_COMPACT_OBJECTS
    if(p->Type == 5 && p->KetjuPNClass == 1) return 1; /* opt-in: stellar-mass non-relic sinks treated as PN */
#endif
#endif
    return 0;
}

/* ============================================================
 *  Helper: distance between two positions (with box wrapping)
 * ============================================================ */
static double particle_distance(double *pos1, double *pos2)
{
    double dx, dy, dz;
    dx = pos1[0] - pos2[0];
    dy = pos1[1] - pos2[1];
    dz = pos1[2] - pos2[2];
#ifdef BOX_PERIODIC
    NEAREST_XYZ(dx, dy, dz, -1);
#endif
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* ============================================================
 *  Helper: GIZMO softening kernel (same as gravtree.cc)
 * ============================================================ */
static double softened_force_factor(double r, double h)
{
    /* Returns G*m/r^2 factor with softening: equivalent to 1/r^2 for r>h */
    if(h <= 0 || r >= h) return 1.0 / (r * r * r);
    double h_inv = 1.0 / h;
    double u = r * h_inv;
    if(u < 0.5) {
        return h_inv * h_inv * h_inv * (32.0/3.0 + u*u*(32.0*u - 38.4));
    } else {
        return h_inv * h_inv * h_inv * (-1.0/(30.0*u*u*u) + 64.0/3.0 + u*u*(-48.0 + u*(38.4 - 32.0/3.0*u)));
    }
}

/* ============================================================
 *  Parse PN terms string (same logic as GADGET4-KETJU)
 * ============================================================ */
static int parse_pn_terms(void)
{
    char terms[20];
    strncpy(terms, All.KetjuPNTerms, sizeof(terms));
    terms[sizeof(terms)-1] = '\0';
    for(int i = 0; terms[i]; i++) terms[i] = tolower(terms[i]);

    if(strcmp(terms, "all") == 0) return KETJU_PN_ALL;
    if(strcmp(terms, "no_spin") == 0) return KETJU_PN_DYNAMIC_ALL;
    if(strcmp(terms, "none") == 0) return KETJU_PN_NONE;

    int flags = KETJU_PN_NONE;
    for(int i = 0; terms[i]; i++) {
        switch(terms[i]) {
            case '2': flags |= KETJU_PN_1_0_ACC; break;
            case '4': flags |= KETJU_PN_2_0_ACC; break;
            case '5': flags |= KETJU_PN_2_5_ACC; break;
            case '6': flags |= KETJU_PN_3_0_ACC; break;
            case '7': flags |= KETJU_PN_3_5_ACC; break;
            case 's': flags |= KETJU_PN_SPIN_ALL; break;
            case 'c': flags |= KETJU_PN_THREEBODY; break;
            default:
                if(ThisTask == 0) printf("KETJU: Warning: unknown PN term character '%c'\n", terms[i]);
                break;
        }
    }
    return flags;
}

/* ============================================================
 *  PHASE 2: Find chain centers and build regions
 * ============================================================ */

/* Gather all chain centers across MPI tasks */
static std::vector<ketju_mpi_particle> gather_chain_centers(void)
{
    /* count local centers */
    int n_local = 0;
    for(int i = 0; i < NumPart; i++) {
        if(is_chain_center(i)) n_local++;
    }

    /* fill local array */
    std::vector<ketju_mpi_particle> local_centers(n_local);
    int k = 0;
    for(int i = 0; i < NumPart; i++) {
        if(is_chain_center(i)) {
            local_centers[k].ID = P[i].ID;
            local_centers[k].Type = P[i].Type;
            local_centers[k].Task = ThisTask;
            local_centers[k].Index = i;
            local_centers[k].Mass = P[i].Mass;
            for(int j = 0; j < 3; j++) {
                local_centers[k].Pos[j] = P[i].Pos[j];
                local_centers[k].Vel[j] = P[i].Vel[j];
            }
#ifdef SINGLE_STAR_SINK_DYNAMICS
            local_centers[k].KetjuPNClass = (P[i].Type == 5) ? KETJU_PN_CLASS(i) : -1;
            for(int j = 0; j < 3; j++) local_centers[k].Spin[j] = (P[i].Type == 5) ? P[i].KetjuSpin[j] : 0;
#endif
            k++;
        }
    }

    /* gather counts */
    std::vector<int> counts(NTask), displs(NTask);
    MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int n_total = 0;
    for(int t = 0; t < NTask; t++) {
        displs[t] = n_total;
        n_total += counts[t];
    }
    if(n_total == 0) return {};

    /* gather data — use raw bytes since we don't have a registered MPI type */
    int local_bytes = n_local * sizeof(ketju_mpi_particle);
    std::vector<int> byte_counts(NTask), byte_displs(NTask);
    for(int t = 0; t < NTask; t++) {
        byte_counts[t] = counts[t] * sizeof(ketju_mpi_particle);
        byte_displs[t] = displs[t] * sizeof(ketju_mpi_particle);
    }

    std::vector<ketju_mpi_particle> all_centers(n_total);
    MPI_Allgatherv(local_centers.data(), local_bytes, MPI_BYTE,
                   all_centers.data(), byte_counts.data(), byte_displs.data(),
                   MPI_BYTE, MPI_COMM_WORLD);

    return all_centers;
}

/* Merge overlapping regions: two centers within 2*r_region get merged */
static std::vector<std::vector<ketju_mpi_particle>> merge_overlapping_regions(
    const std::vector<ketju_mpi_particle> &centers, double merge_radius)
{
    int n = centers.size();
    std::vector<int> group(n);
    for(int i = 0; i < n; i++) group[i] = i;

    /* union-find */
    auto find = [&](int x) {
        while(group[x] != x) { group[x] = group[group[x]]; x = group[x]; }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if(a != b) group[a] = b;
    };

    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            double dist = particle_distance(
                const_cast<double*>(centers[i].Pos),
                const_cast<double*>(centers[j].Pos));
            if(dist < merge_radius) unite(i, j);
        }
    }

    /* group centers by their root */
    std::vector<std::vector<ketju_mpi_particle>> regions;
    std::vector<int> root_to_region(n, -1);
    for(int i = 0; i < n; i++) {
        int r = find(i);
        if(root_to_region[r] < 0) {
            root_to_region[r] = regions.size();
            regions.push_back({});
        }
        regions[root_to_region[r]].push_back(centers[i]);
    }

    return regions;
}

/* Find local particles within r_region of any center in a region */
static std::set<int> find_local_members(const std::vector<ketju_mpi_particle> &centers, double radius)
{
    std::set<int> members;
    for(int i = 0; i < NumPart; i++) {
        if(!is_chain_eligible(i)) continue;
        for(size_t c = 0; c < centers.size(); c++) {
            double dist = particle_distance(P[i].Pos, const_cast<double*>(centers[c].Pos));
            if(dist < radius) { members.insert(i); break; }
        }
    }
    return members;
}

/* ============================================================
 *  PHASE 3: MPI gather and integrator setup
 * ============================================================ */

static void setup_integrator(KetjuRegion &reg)
{
    /* gather particle counts per task */
    int n_local = reg.local_member_indices.size();
    std::vector<int> counts(NTask);
    MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    reg.total_particle_count = 0;
    for(int t = 0; t < NTask; t++) reg.total_particle_count += counts[t];
    if(reg.total_particle_count < 2) return; /* need at least 2 particles */

    /* create affected_tasks communicator: only tasks with particles */
    {
        std::vector<int> affected_indices;
        for(int t = 0; t < NTask; t++) { if(counts[t] > 0) affected_indices.push_back(t); }
        if(affected_indices.empty()) return;
        reg.affected_tasks.init(affected_indices, 0);
        reg.affected_sim_indices = affected_indices;
        reg.particle_counts_on_affected.resize(affected_indices.size());
        for(size_t i = 0; i < affected_indices.size(); i++)
            reg.particle_counts_on_affected[i] = counts[affected_indices[i]];
    }

    /* create compute_tasks communicator: contiguous range from cost-based allocator.
     * Scatter uses MPI_Scatterv on affected_tasks, so compute/affected groups can
     * safely differ — data is forwarded via Send/Recv between roots if needed. */
    {
        std::vector<int> compute_indices;
        for(int t = reg.compute_info.first_task_index; t <= reg.compute_info.final_task_index; t++)
            compute_indices.push_back(t);
        reg.compute_tasks.init(compute_indices, 1);
        reg.compute_tasks.set_common_root(reg.affected_tasks);
    }

    /* fill local particle data */
    std::vector<ketju_mpi_particle> local_parts(n_local);
    int k = 0;
    for(int idx : reg.local_member_indices) {
        ketju_mpi_particle &mp = local_parts[k];
        mp.ID = P[idx].ID;
        mp.Type = P[idx].Type;
        mp.Task = ThisTask;
        mp.Index = idx;
        mp.Mass = P[idx].Mass;
        for(int j = 0; j < 3; j++) {
            mp.Pos[j] = P[idx].Pos[j];
            mp.Vel[j] = P[idx].Vel[j];
        }
#ifdef SINGLE_STAR_SINK_DYNAMICS
        mp.KetjuPNClass = (P[idx].Type == 5) ? KETJU_PN_CLASS(idx) : -1;
        for(int j = 0; j < 3; j++) mp.Spin[j] = (P[idx].Type == 5) ? P[idx].KetjuSpin[j] : 0;
#endif
        k++;
    }

    /* gather all particles on affected root */
    if(reg.affected_tasks.is_member()) {
        int n_aff = reg.affected_tasks.size;
        std::vector<int> aff_counts(n_aff), byte_counts(n_aff), byte_displs(n_aff);
        MPI_Allgather(&n_local, 1, MPI_INT, aff_counts.data(), 1, MPI_INT, reg.affected_tasks.comm);
        int total_bytes = 0;
        for(int t = 0; t < n_aff; t++) {
            byte_counts[t] = aff_counts[t] * sizeof(ketju_mpi_particle);
            byte_displs[t] = total_bytes;
            total_bytes += byte_counts[t];
        }
        if(reg.affected_tasks.is_root()) {
            reg.all_particles.resize(reg.total_particle_count);
        }
        MPI_Gatherv(local_parts.data(), n_local * sizeof(ketju_mpi_particle), MPI_BYTE,
                    reg.all_particles.data(), byte_counts.data(), byte_displs.data(),
                    MPI_BYTE, reg.affected_tasks.root, reg.affected_tasks.comm);
    }

    /* if compute root differs from affected root, forward data */
    if(reg.affected_tasks.root_sim != reg.compute_tasks.root_sim) {
        if(ThisTask == reg.affected_tasks.root_sim) {
            MPI_Send(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                     MPI_BYTE, reg.compute_tasks.root_sim, 999, MPI_COMM_WORLD);
        }
        if(ThisTask == reg.compute_tasks.root_sim) {
            reg.all_particles.resize(reg.total_particle_count);
            MPI_Recv(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                     MPI_BYTE, reg.affected_tasks.root_sim, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    /* broadcast particle data to all compute tasks */
    if(reg.compute_tasks.is_member()) {
        if(!reg.compute_tasks.is_root()) reg.all_particles.resize(reg.total_particle_count);
        MPI_Bcast(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                  MPI_BYTE, reg.compute_tasks.root, reg.compute_tasks.comm);
    }

    /* set up integrator on compute tasks */
    if(reg.compute_tasks.is_member()) {
        /* count PN particles */
        reg.num_pn_particles = 0;
        for(int i = 0; i < reg.total_particle_count; i++) {
            if(is_pn_particle(&reg.all_particles[i])) reg.num_pn_particles++;
        }

        /* sort: PN particles first */
        std::stable_sort(reg.all_particles.begin(), reg.all_particles.end(),
            [](const ketju_mpi_particle &a, const ketju_mpi_particle &b) {
                int a_pn = 0, b_pn = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
                if(a.Type == 5 && a.KetjuPNClass == 0) a_pn = 1;
                if(b.Type == 5 && b.KetjuPNClass == 0) b_pn = 1;
#endif
#ifdef KETJU_PN_COMPACT_OBJECTS
                if(a.Type == 5 && a.KetjuPNClass == 1) a_pn = 1;
                if(b.Type == 5 && b.KetjuPNClass == 1) b_pn = 1;
#endif
                return a_pn > b_pn;
            });

        /* compute CoM position and velocity */
        double total_mass = 0;
        for(int j = 0; j < 3; j++) { reg.com_pos[j] = 0; reg.com_vel[j] = 0; }
        for(int i = 0; i < reg.total_particle_count; i++) {
            for(int j = 0; j < 3; j++) {
                reg.com_pos[j] += reg.all_particles[i].Mass * reg.all_particles[i].Pos[j];
                reg.com_vel[j] += reg.all_particles[i].Mass * reg.all_particles[i].Vel[j];
            }
            total_mass += reg.all_particles[i].Mass;
        }
        for(int j = 0; j < 3; j++) { reg.com_pos[j] /= total_mass; reg.com_vel[j] /= total_mass; }

        /* allocate integrator — use compute_tasks.comm for parallel integration */
        reg.integrator = (struct ketju_system *)calloc(1, sizeof(struct ketju_system));
        int n_other = reg.total_particle_count - reg.num_pn_particles;
        ketju_create_system(reg.integrator, reg.num_pn_particles, n_other, reg.compute_tasks.comm);

        /* set units */
        reg.integrator->constants->G = All.G;
#if defined(C_LIGHT_CODE)
        reg.integrator->constants->c = C_LIGHT_CODE;
#else
        reg.integrator->constants->c = 2.9979e10 / All.UnitVelocity_in_cm_per_s;
#endif

        /* set options */
        reg.integrator->options->PN_flags = parse_pn_terms();
        reg.integrator->options->gbs_relative_tolerance = All.KetjuIntegrationTolerance;
        reg.integrator->options->enable_bh_mergers = (All.KetjuEnableBHMergerKicks >= 0);
        reg.integrator->options->enable_bh_merger_kicks = All.KetjuEnableBHMergerKicks;
        if(All.KetjuMaxStepCount > 0) reg.integrator->options->max_step_count = All.KetjuMaxStepCount;
        if(All.KetjuUseStarStarSoftening) {
            double h = All.ForceSoftening[4] * All.cf_atime;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            if(All.ForceSoftening[5] > 0) h = DMIN(h, All.ForceSoftening[5] * All.cf_atime);
#endif
            reg.integrator->options->star_star_softening = h;
        }

        /* fill particle data (relative to CoM, converted to physical coordinates) */
        reg.extra_data.resize(reg.total_particle_count);
        struct ketju_system_physical_state *ps = reg.integrator->physical_state;
        ps->time = 0;
        double a = All.cf_atime; /* scale factor: 1 for non-cosmo */
        for(int i = 0; i < reg.total_particle_count; i++) {
            ps->mass[i] = reg.all_particles[i].Mass;
            for(int j = 0; j < 3; j++) {
                /* comoving -> physical: pos_phys = (pos_comov - com_comov) * a */
                ps->pos[i][j] = (reg.all_particles[i].Pos[j] - reg.com_pos[j]) * a;
                /* GIZMO Vel stores canonical momentum p = a^2 * dx_comov/dt_phys.
                 * Peculiar velocity = Vel / a. Physical velocity = v_pec + H*r_phys.
                 * MSTAR needs physical velocity (includes Hubble flow).
                 * In non-cosmo: a=1, Vel is already physical, H=0. */
                double v_pec_rel = (reg.all_particles[i].Vel[j] - reg.com_vel[j]) / a;
                ps->vel[i][j] = v_pec_rel;
                if(All.ComovingIntegrationOn)
                    ps->vel[i][j] += hubble_function(a) * ps->pos[i][j];
            }
            reg.extra_data[i].ID = reg.all_particles[i].ID;
            reg.extra_data[i].Task = reg.all_particles[i].Task;
            reg.extra_data[i].Index = reg.all_particles[i].Index;
            reg.extra_data[i].Type = reg.all_particles[i].Type;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            reg.extra_data[i].KetjuPNClass = reg.all_particles[i].KetjuPNClass;
            if(i < reg.num_pn_particles) {
                for(int j = 0; j < 3; j++) ps->spin[i][j] = reg.all_particles[i].Spin[j];
            }
#endif
        }
        reg.integrator->particle_extra_data = reg.extra_data.data();
        reg.integrator->particle_extra_data_elem_size = sizeof(ketju_extra_data);

        if(reg.compute_tasks.is_root() && ThisTask == 0) {
            printf("KETJU: Region with %d particles (%d PN), compute tasks %d-%d\n",
                   reg.total_particle_count, reg.num_pn_particles,
                   reg.compute_info.first_task_index, reg.compute_info.final_task_index);
        }
    }

    /* broadcast CoM data to all tasks that need it (affected group for scatter) */
    if(reg.affected_tasks.is_member()) {
        MPI_Bcast(reg.com_pos, 3, MPI_DOUBLE, reg.affected_tasks.root, reg.affected_tasks.comm);
        MPI_Bcast(reg.com_vel, 3, MPI_DOUBLE, reg.affected_tasks.root, reg.affected_tasks.comm);
    }
}

/* Setup integrator reusing existing communicators (region persistence).
 * Skips MPI_Comm_create_group — the expensive part — but still gathers
 * particle data and builds the integrator state for this step. */
static void setup_integrator_reuse(KetjuRegion &reg)
{
    if(reg.total_particle_count < 2) return;

    int n_local = reg.local_member_indices.size();

    /* fill local particle data */
    std::vector<ketju_mpi_particle> local_parts(n_local);
    int k = 0;
    for(int idx : reg.local_member_indices) {
        ketju_mpi_particle &mp = local_parts[k];
        mp.ID = P[idx].ID;
        mp.Type = P[idx].Type;
        mp.Task = ThisTask;
        mp.Index = idx;
        mp.Mass = P[idx].Mass;
        for(int j = 0; j < 3; j++) {
            mp.Pos[j] = P[idx].Pos[j];
            mp.Vel[j] = P[idx].Vel[j];
        }
#ifdef SINGLE_STAR_SINK_DYNAMICS
        mp.KetjuPNClass = (P[idx].Type == 5) ? KETJU_PN_CLASS(idx) : -1;
        for(int j = 0; j < 3; j++) mp.Spin[j] = (P[idx].Type == 5) ? P[idx].KetjuSpin[j] : 0;
#endif
        k++;
    }

    /* gather all particles on affected root (reusing existing communicator) */
    if(reg.affected_tasks.is_member()) {
        int n_aff = reg.affected_tasks.size;
        std::vector<int> aff_counts(n_aff), byte_counts(n_aff), byte_displs(n_aff);
        MPI_Allgather(&n_local, 1, MPI_INT, aff_counts.data(), 1, MPI_INT, reg.affected_tasks.comm);
        int total_bytes = 0;
        for(int t = 0; t < n_aff; t++) {
            byte_counts[t] = aff_counts[t] * sizeof(ketju_mpi_particle);
            byte_displs[t] = total_bytes;
            total_bytes += byte_counts[t];
        }
        if(reg.affected_tasks.is_root()) {
            reg.all_particles.resize(reg.total_particle_count);
        }
        MPI_Gatherv(local_parts.data(), n_local * sizeof(ketju_mpi_particle), MPI_BYTE,
                    reg.all_particles.data(), byte_counts.data(), byte_displs.data(),
                    MPI_BYTE, reg.affected_tasks.root, reg.affected_tasks.comm);
    }

    /* forward from affected root to compute root if they differ */
    if(reg.affected_tasks.root_sim != reg.compute_tasks.root_sim) {
        if(ThisTask == reg.affected_tasks.root_sim) {
            MPI_Send(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                     MPI_BYTE, reg.compute_tasks.root_sim, 999, MPI_COMM_WORLD);
        }
        if(ThisTask == reg.compute_tasks.root_sim) {
            reg.all_particles.resize(reg.total_particle_count);
            MPI_Recv(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                     MPI_BYTE, reg.affected_tasks.root_sim, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    /* broadcast to all compute tasks */
    if(reg.compute_tasks.is_member()) {
        if(!reg.compute_tasks.is_root()) reg.all_particles.resize(reg.total_particle_count);
        MPI_Bcast(reg.all_particles.data(), reg.total_particle_count * sizeof(ketju_mpi_particle),
                  MPI_BYTE, reg.compute_tasks.root, reg.compute_tasks.comm);
    }

    /* set up integrator on compute tasks (same as setup_integrator from here) */
    if(reg.compute_tasks.is_member()) {
        reg.num_pn_particles = 0;
        for(int i = 0; i < reg.total_particle_count; i++) {
            if(is_pn_particle(&reg.all_particles[i])) reg.num_pn_particles++;
        }

        std::stable_sort(reg.all_particles.begin(), reg.all_particles.end(),
            [](const ketju_mpi_particle &a, const ketju_mpi_particle &b) {
                int a_pn = 0, b_pn = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
                if(a.Type == 5 && a.KetjuPNClass == 0) a_pn = 1;
                if(b.Type == 5 && b.KetjuPNClass == 0) b_pn = 1;
#endif
#ifdef KETJU_PN_COMPACT_OBJECTS
                if(a.Type == 5 && a.KetjuPNClass == 1) a_pn = 1;
                if(b.Type == 5 && b.KetjuPNClass == 1) b_pn = 1;
#endif
                return a_pn > b_pn;
            });

        double total_mass = 0;
        for(int j = 0; j < 3; j++) { reg.com_pos[j] = 0; reg.com_vel[j] = 0; }
        for(int i = 0; i < reg.total_particle_count; i++) {
            for(int j = 0; j < 3; j++) {
                reg.com_pos[j] += reg.all_particles[i].Mass * reg.all_particles[i].Pos[j];
                reg.com_vel[j] += reg.all_particles[i].Mass * reg.all_particles[i].Vel[j];
            }
            total_mass += reg.all_particles[i].Mass;
        }
        for(int j = 0; j < 3; j++) { reg.com_pos[j] /= total_mass; reg.com_vel[j] /= total_mass; }

        reg.integrator = (struct ketju_system *)calloc(1, sizeof(struct ketju_system));
        int n_other = reg.total_particle_count - reg.num_pn_particles;
        ketju_create_system(reg.integrator, reg.num_pn_particles, n_other, reg.compute_tasks.comm);

        reg.integrator->constants->G = All.G;
#if defined(C_LIGHT_CODE)
        reg.integrator->constants->c = C_LIGHT_CODE;
#else
        reg.integrator->constants->c = 2.9979e10 / All.UnitVelocity_in_cm_per_s;
#endif
        reg.integrator->options->PN_flags = parse_pn_terms();
        reg.integrator->options->gbs_relative_tolerance = All.KetjuIntegrationTolerance;
        reg.integrator->options->enable_bh_mergers = (All.KetjuEnableBHMergerKicks >= 0);
        reg.integrator->options->enable_bh_merger_kicks = All.KetjuEnableBHMergerKicks;
        if(All.KetjuMaxStepCount > 0) reg.integrator->options->max_step_count = All.KetjuMaxStepCount;
        if(All.KetjuUseStarStarSoftening) {
            double h = All.ForceSoftening[4] * All.cf_atime;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            if(All.ForceSoftening[5] > 0) h = DMIN(h, All.ForceSoftening[5] * All.cf_atime);
#endif
            reg.integrator->options->star_star_softening = h;
        }

        reg.extra_data.resize(reg.total_particle_count);
        struct ketju_system_physical_state *ps = reg.integrator->physical_state;
        ps->time = 0;
        double a = All.cf_atime;
        for(int i = 0; i < reg.total_particle_count; i++) {
            ps->mass[i] = reg.all_particles[i].Mass;
            for(int j = 0; j < 3; j++) {
                ps->pos[i][j] = (reg.all_particles[i].Pos[j] - reg.com_pos[j]) * a;
                double v_pec_rel = (reg.all_particles[i].Vel[j] - reg.com_vel[j]) / a;
                ps->vel[i][j] = v_pec_rel;
                if(All.ComovingIntegrationOn)
                    ps->vel[i][j] += hubble_function(a) * ps->pos[i][j];
            }
            reg.extra_data[i].ID = reg.all_particles[i].ID;
            reg.extra_data[i].Task = reg.all_particles[i].Task;
            reg.extra_data[i].Index = reg.all_particles[i].Index;
            reg.extra_data[i].Type = reg.all_particles[i].Type;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            reg.extra_data[i].KetjuPNClass = reg.all_particles[i].KetjuPNClass;
            if(i < reg.num_pn_particles) {
                for(int j = 0; j < 3; j++) ps->spin[i][j] = reg.all_particles[i].Spin[j];
            }
#endif
        }
        reg.integrator->particle_extra_data = reg.extra_data.data();
        reg.integrator->particle_extra_data_elem_size = sizeof(ketju_extra_data);

        if(reg.compute_tasks.is_root() && ThisTask == 0) {
            printf("KETJU: Region with %d particles (%d PN), compute tasks %d-%d [cached comms]\n",
                   reg.total_particle_count, reg.num_pn_particles,
                   reg.compute_info.first_task_index, reg.compute_info.final_task_index);
        }
    }

    if(reg.affected_tasks.is_member()) {
        MPI_Bcast(reg.com_pos, 3, MPI_DOUBLE, reg.affected_tasks.root, reg.affected_tasks.comm);
        MPI_Bcast(reg.com_vel, 3, MPI_DOUBLE, reg.affected_tasks.root, reg.affected_tasks.comm);
    }
}

/* ============================================================
 *  PHASE 4: Negative half-step kick and integration
 * ============================================================ */

/* Subtract softened pairwise gravity between chain members from their velocities.
 * This undoes the tree force contribution so MSTAR can replace it with exact forces. */
static void do_negative_halfstep_kick(KetjuRegion &reg, double kick_factor)
{
    if(!reg.compute_tasks.is_member() || !reg.integrator) return;

    int n = reg.integrator->num_particles;
    struct ketju_system_physical_state *ps = reg.integrator->physical_state;

    /* softening: use the star softening (comoving -> physical) */
    double h = All.ForceSoftening[4] * All.cf_atime;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    if(All.ForceSoftening[5] > 0) h = DMIN(h, All.ForceSoftening[5] * All.cf_atime);
#endif

    /* distribute N² loop across compute tasks */
    int loop_start = loop_scheduling_block_edge(n, reg.compute_tasks.size, reg.compute_tasks.rank);
    int loop_end   = loop_scheduling_block_edge(n, reg.compute_tasks.size, reg.compute_tasks.rank + 1);

    std::vector<double> dv(3 * n, 0.0);

    for(int i = loop_start; i < loop_end; i++) {
        for(int j = i + 1; j < n; j++) {
            double dr[3], r2 = 0;
            for(int k = 0; k < 3; k++) {
                dr[k] = ps->pos[i][k] - ps->pos[j][k];
                r2 += dr[k] * dr[k];
            }
            double r = sqrt(r2);
            if(r == 0) continue;

            double fac = kick_factor * All.G * softened_force_factor(r, h);
            for(int k = 0; k < 3; k++) {
                dv[3*i + k] += ps->mass[j] * fac * dr[k];
                dv[3*j + k] -= ps->mass[i] * fac * dr[k];
            }
        }
    }

    /* collect results across compute tasks and apply */
    MPI_Allreduce(MPI_IN_PLACE, dv.data(), 3 * n, MPI_DOUBLE, MPI_SUM, reg.compute_tasks.comm);

    for(int i = 0; i < n; i++) {
        for(int k = 0; k < 3; k++) {
            ps->vel[i][k] += dv[3*i + k];
        }
    }
}

/* ============================================================
 *  Expand tight binaries: if the tightest binary has a period
 *  much shorter than the timestep and is an outlier, double its
 *  semi-major axis to prevent integrator stalling.
 *  Based on the public GADGET4-KETJU algorithm (Mannerkoski+ 2023).
 * ============================================================ */
static void expand_tight_binaries(KetjuRegion &reg, double dt_physical)
{
    if(!reg.compute_tasks.is_root() || !reg.integrator) return;
    if(All.KetjuExpandBinariesFactor <= 0) return; /* disabled */

    int n = reg.integrator->num_particles;
    struct ketju_system_physical_state *ps = reg.integrator->physical_state;

    /* find the tightest bound binary */
    int best_i = -1, best_j = -1;
    double min_period = 1e30;
    int n_tight = 0; /* count of binaries with period < 2*min_period */

    for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
            double dr[3], dv[3], r2 = 0, v2 = 0;
            for(int k = 0; k < 3; k++) {
                dr[k] = ps->pos[i][k] - ps->pos[j][k];
                dv[k] = ps->vel[i][k] - ps->vel[j][k];
                r2 += dr[k] * dr[k]; v2 += dv[k] * dv[k];
            }
            double GM = All.G * (ps->mass[i] + ps->mass[j]);
            double E = 0.5 * v2 - GM / sqrt(r2);
            if(E >= 0) continue; /* unbound */
            double a = -GM / (2.0 * E);
            double period = 2.0 * M_PI * a * sqrt(a / GM);
            if(period < min_period) {
                min_period = period; best_i = i; best_j = j;
            }
        }
    }

    if(best_i < 0) return; /* no bound binaries */
    if(min_period > All.KetjuExpandBinariesFactor * dt_physical) return; /* not tight enough */

    /* count binaries with period < 2*min_period */
    for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
            double dr[3], dv[3], r2 = 0, v2 = 0;
            for(int k = 0; k < 3; k++) {
                dr[k] = ps->pos[i][k] - ps->pos[j][k];
                dv[k] = ps->vel[i][k] - ps->vel[j][k];
                r2 += dr[k] * dr[k]; v2 += dv[k] * dv[k];
            }
            double GM = All.G * (ps->mass[i] + ps->mass[j]);
            double E = 0.5 * v2 - GM / sqrt(r2);
            if(E >= 0) continue;
            double a = -GM / (2.0 * E);
            double period = 2.0 * M_PI * a * sqrt(a / GM);
            if(period < 2.0 * min_period) n_tight++;
        }
    }

    /* only expand if it's an outlier (< max(N/100, 5) tight binaries) */
    if(n_tight > DMAX(n / 100, 5)) return;

    /* double the semi-major axis: scale relative velocity by 1/sqrt(2)
     * so a_new = 2*a_old (via virial theorem: E_new = E_old/2) */
    double scale = 1.0 / sqrt(2.0);
    double com_vel[3];
    double mtot = ps->mass[best_i] + ps->mass[best_j];
    for(int k = 0; k < 3; k++) {
        com_vel[k] = (ps->mass[best_i] * ps->vel[best_i][k] + ps->mass[best_j] * ps->vel[best_j][k]) / mtot;
    }
    for(int k = 0; k < 3; k++) {
        ps->vel[best_i][k] = com_vel[k] + scale * (ps->vel[best_i][k] - com_vel[k]);
        ps->vel[best_j][k] = com_vel[k] + scale * (ps->vel[best_j][k] - com_vel[k]);
    }

    printf("KETJU: expanded tight binary (period=%g -> %g, dt=%g) in region with %d particles\n",
           min_period, min_period * 2.0 * sqrt(2.0), dt_physical, n);
}

/* Run the full integration step for one region.
 * Uses proper cosmological kick factors for comoving integration
 * (based on public GADGET4-KETJU, Mannerkoski+ 2023). */
static void integrate_region(KetjuRegion &reg, double dt_physical)
{
    if(dt_physical <= 0) return;

    /* Compute kick factors: in comoving integration, the half-step kick
     * factor is the gravkick integral (dt/a²), not just 0.5*dt_physical.
     * The scale factor correction accounts for the tree force being computed
     * in comoving coordinates while KETJU works in physical coordinates. */
    double first_halfkick, second_halfkick;
    if(All.ComovingIntegrationOn) {
        integertime t0 = All.Ti_Current;
        integertime t1 = t0 + reg.ti_step;
        integertime tmid = t0 + reg.ti_step / 2;
        double a0 = All.TimeBegin * exp(t0 * All.Timebase_interval);
        double a1 = All.TimeBegin * exp(t1 * All.Timebase_interval);
        first_halfkick  = get_gravkick_factor(t0, tmid, -1, 0) * a0;
        second_halfkick = get_gravkick_factor(tmid, t1, -1, 0) * a1;
    } else {
        first_halfkick  = 0.5 * dt_physical;
        second_halfkick = 0.5 * dt_physical;
    }

    /* subtract first half of tree force */
    do_negative_halfstep_kick(reg, first_halfkick);

    /* expand tight binaries before integration to prevent stalling.
     * Runs on compute root only — broadcast updated velocities to all compute tasks. */
    expand_tight_binaries(reg, dt_physical);
    if(reg.compute_tasks.is_member() && reg.integrator && reg.compute_tasks.size > 1) {
        int n = reg.integrator->num_particles;
        MPI_Bcast(reg.integrator->physical_state->vel, 3 * n, MPI_DOUBLE,
                  reg.compute_tasks.root, reg.compute_tasks.comm);
    }

    /* run MSTAR integrator (all compute tasks participate via compute_tasks.comm) */
    if(reg.compute_tasks.is_member() && reg.integrator) {
        ketju_run_integrator(reg.integrator, dt_physical);

        int n_steps = reg.integrator->perf->successful_steps + reg.integrator->perf->failed_steps;
        double dE = reg.integrator->perf->relative_energy_error;

        if(reg.compute_tasks.is_root())
        printf("KETJU [task %d]: Integrated %d particles for dt=%g, %d steps (%d failed), dE/E=%g\n",
               ThisTask, reg.total_particle_count, dt_physical,
               reg.integrator->perf->successful_steps,
               reg.integrator->perf->failed_steps, dE);

        if(All.KetjuMaxStepCount > 0 && n_steps >= All.KetjuMaxStepCount) {
            printf("KETJU WARNING: integration hit max step count (%d) — results may be inaccurate!\n",
                   All.KetjuMaxStepCount);
        }
        if(fabs(dE) > 1e-4) {
            printf("KETJU WARNING: large energy error dE/E=%g in region with %d particles\n",
                   dE, reg.total_particle_count);
        }
    }

    /* subtract second half of tree force */
    do_negative_halfstep_kick(reg, second_halfkick);
}

/* ============================================================
 *  PHASE 5a: Merger handling
 *
 *  After MSTAR integration, merged particles have mass=0.
 *  Detect these, transfer metals to the survivor (nearest
 *  surviving particle), and mark merged particles for removal.
 *  Controlled by KETJU_MERGE_STARS (Type 4+4) and
 *  KETJU_MERGE_BH (Type 5+5, Type 5+4).
 * ============================================================ */

/* Merger action descriptor — determined on affected root, broadcast to all affected tasks */
struct ketju_merger_action {
    MyIDType id_merged, id_survivor;
    int type_merged, type_survivor;
    int task_merged, task_survivor;
    int idx_merged, idx_survivor;
    double M_merged, M_survivor_old, M_survivor_new;
};

#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
static void handle_mergers(KetjuRegion &reg,
                           const std::vector<double> &final_abs_pos,
                           const std::vector<double> &final_mass,
                           const std::vector<double> &original_mass)
{
    if(!reg.affected_tasks.is_member()) return;

    /* ---- Phase 1: affected root determines merger pairs ---- */
    int n_mergers = 0;
    std::vector<ketju_merger_action> actions;

    if(reg.affected_tasks.is_root()) {
        int n = reg.total_particle_count;
        for(int i = 0; i < n; i++) {
            if(final_mass[i] > 0) continue;

            int type_merged = reg.extra_data[i].Type;

            /* find nearest surviving particle */
            int survivor = -1;
            double min_dist2 = 1e60;
            for(int j = 0; j < n; j++) {
                if(j == i || final_mass[j] <= 0) continue;
                double d2 = 0;
                for(int k = 0; k < 3; k++) {
                    double dx = final_abs_pos[3*i + k] - final_abs_pos[3*j + k];
                    d2 += dx * dx;
                }
                if(d2 < min_dist2) { min_dist2 = d2; survivor = j; }
            }
            if(survivor < 0) continue;

            int type_survivor = reg.extra_data[survivor].Type;

            int merge_allowed = 0;
#ifdef KETJU_MERGE_BH
            if(type_merged == 5 || type_survivor == 5) merge_allowed = 1;
#endif
#ifdef KETJU_MERGE_STARS
            if(type_merged == 4 && type_survivor == 4) merge_allowed = 1;
#endif
            if(!merge_allowed) continue;

            double M_merged = original_mass[i];
            double M_survivor_old = original_mass[survivor];
            double M_survivor_new = final_mass[survivor];
            if(M_survivor_new <= 0 || M_merged <= 0) continue;

            ketju_merger_action act;
            act.id_merged      = reg.extra_data[i].ID;
            act.id_survivor    = reg.extra_data[survivor].ID;
            act.type_merged    = type_merged;
            act.type_survivor  = type_survivor;
            act.task_merged    = reg.extra_data[i].Task;
            act.task_survivor  = reg.extra_data[survivor].Task;
            act.idx_merged     = reg.extra_data[i].Index;
            act.idx_survivor   = reg.extra_data[survivor].Index;
            act.M_merged       = M_merged;
            act.M_survivor_old = M_survivor_old;
            act.M_survivor_new = M_survivor_new;
            actions.push_back(act);

            printf("KETJU MERGER: ID %llu (Type %d, M=%g) merged into ID %llu (Type %d, M=%g->%g)\n",
                   (unsigned long long)act.id_merged, type_merged, M_merged,
                   (unsigned long long)act.id_survivor, type_survivor,
                   M_survivor_old, M_survivor_new);
        }
        n_mergers = actions.size();
    }

    /* ---- Phase 2: broadcast merger list to all affected tasks ---- */
    MPI_Bcast(&n_mergers, 1, MPI_INT, reg.affected_tasks.root, reg.affected_tasks.comm);
    if(n_mergers == 0) return;

    if(!reg.affected_tasks.is_root()) actions.resize(n_mergers);
    MPI_Bcast(actions.data(), n_mergers * sizeof(ketju_merger_action), MPI_BYTE,
              reg.affected_tasks.root, reg.affected_tasks.comm);

    /* ---- Phase 3: each task executes its local part ---- */
#ifdef METALS
    const int buf_size = NUM_METAL_SPECIES;
#else
    const int buf_size = 0;
#endif

    for(int m = 0; m < n_mergers; m++) {
        ketju_merger_action &act = actions[m];

#ifdef METALS
        if(buf_size > 0) {
            std::vector<double> merged_metals(buf_size, 0.0);

            if(act.task_merged == ThisTask) {
                for(int k = 0; k < NUM_METAL_SPECIES; k++)
                    merged_metals[k] = P[act.idx_merged].Metallicity[k];
            }

            /* point-to-point transfer of metals */
            if(act.task_merged != act.task_survivor) {
                if(act.task_merged == ThisTask)
                    MPI_Send(merged_metals.data(), buf_size, MPI_DOUBLE,
                             act.task_survivor, 997, MPI_COMM_WORLD);
                if(act.task_survivor == ThisTask)
                    MPI_Recv(merged_metals.data(), buf_size, MPI_DOUBLE,
                             act.task_merged, 997, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if(act.task_survivor == ThisTask) {
                for(int k = 0; k < NUM_METAL_SPECIES; k++)
                    P[act.idx_survivor].Metallicity[k] =
                        (act.M_survivor_old * P[act.idx_survivor].Metallicity[k] +
                         act.M_merged * merged_metals[k]) / act.M_survivor_new;
            }
        }
#endif

        /* mark merged particle for removal */
        if(act.task_merged == ThisTask) {
            P[act.idx_merged].Mass = 0;
        }
    }
}
#endif /* KETJU_MERGE_STARS || KETJU_MERGE_BH */

/* ============================================================
 *  PHASE 5b: Scatter results and apply velocity trick
 *
 *  After MSTAR integration, positions are in CoM-relative frame.
 *  Convert back to GIZMO frame and set up velocity trick so that
 *  GIZMO's standard drift lands particles at the correct positions.
 *
 *  Velocity trick:
 *    Vel = com_vel + (r_new - r_old) / dt_drift
 *  so that after drift: Pos += Vel * dt_drift
 *    = com_vel*dt_drift + r_new - r_old
 *    = (com_pos + com_vel*dt_drift) + r_new - (com_pos + r_old)
 *    = com_pos_drifted + r_new - Pos_current
 *  => new Pos = com_pos_drifted + r_new  (correct!)
 *
 *  True velocity (KetjuFinalVel) = com_vel + v_integrator
 *  is swapped in after the drift by ketju_set_final_velocities().
 * ============================================================ */

/* Per-particle data packed for MPI_Scatterv (one struct per particle) */
struct ketju_scatter_particle {
    MyIDType ID;
    int Task;
    int Index;
    int Type;
    double Mass;
    double Pos[3];       /* absolute comoving position */
    double Vel[3];       /* CoM-relative peculiar velocity */
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
    double OriginalMass;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    short int KetjuPNClass;
    double Spin[3];      /* only meaningful for PN particles */
    int is_pn;           /* 1 if this is a PN particle (index < n_pn) */
#endif
};

/* comparison for sorting scatter particles by Task (for MPI_Scatterv) */
static int scatter_particle_task_cmp(const void *a, const void *b)
{
    const ketju_scatter_particle *pa = (const ketju_scatter_particle *)a;
    const ketju_scatter_particle *pb = (const ketju_scatter_particle *)b;
    return (pa->Task > pb->Task) - (pa->Task < pb->Task);
}

static void scatter_results(KetjuRegion &reg)
{
    int n = reg.total_particle_count;
    bool is_affected     = reg.affected_tasks.is_member();
    bool is_compute_root = (ThisTask == reg.compute_tasks.root_sim);
    bool needs_forward   = (reg.compute_tasks.root_sim != reg.affected_tasks.root_sim);

    /* compute_root must participate to pack and Send the result to affected_root,
     * even if it is not itself a member of the affected_tasks group. */
    if(!is_affected && !(is_compute_root && needs_forward)) return;

    if(is_affected) {
        MPI_Bcast(&n, 1, MPI_INT, reg.affected_tasks.root, reg.affected_tasks.comm);
    }
    if(n < 2) return;

    int n_aff = reg.affected_tasks.size;
    std::vector<ketju_scatter_particle> scatter_buf;

    /* ---- Phase 1: compute root packs results into scatter_buf ---- */
    /* compute end-of-step scale factor for position/velocity back-conversion */
    double a_end = All.cf_atime; /* 1 for non-cosmo */
    double hubble_end = 0; /* H(a_end); 0 for non-cosmo */
    if(All.ComovingIntegrationOn) {
        a_end = All.cf_atime * exp(reg.ti_step * All.Timebase_interval);
        hubble_end = hubble_function(a_end);
    }

    if(reg.compute_tasks.is_root() && reg.integrator) {
        struct ketju_system_physical_state *ps = reg.integrator->physical_state;
        scatter_buf.resize(n);

        for(int i = 0; i < n; i++) {
            ketju_scatter_particle &sp = scatter_buf[i];
            sp.ID    = reg.extra_data[i].ID;
            sp.Task  = reg.extra_data[i].Task;
            sp.Index = reg.extra_data[i].Index;
            sp.Type  = reg.extra_data[i].Type;
            sp.Mass  = ps->mass[i];
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
            sp.OriginalMass = reg.all_particles[i].Mass;
#endif
            for(int j = 0; j < 3; j++) {
                /* physical -> comoving: use a_end since positions are at end of step */
                sp.Pos[j] = reg.com_pos[j] + ps->pos[i][j] / a_end;
                /* physical velocity -> canonical momentum: remove Hubble flow, multiply by a_end.
                 * sp.Vel stores CoM-relative canonical momentum (com_vel added later in velocity trick).
                 * In non-cosmo: a_end=1, hubble_end=0, so sp.Vel = ps->vel (unchanged). */
                sp.Vel[j] = (ps->vel[i][j] - hubble_end * ps->pos[i][j]) * a_end;
            }
#ifdef SINGLE_STAR_SINK_DYNAMICS
            sp.KetjuPNClass = reg.extra_data[i].KetjuPNClass;
            sp.is_pn = (i < reg.num_pn_particles) ? 1 : 0;
            for(int j = 0; j < 3; j++)
                sp.Spin[j] = (i < reg.num_pn_particles) ? ps->spin[i][j] : 0;
#endif
        }

        /* log merger info from integrator */
        if(reg.integrator->num_mergers > 0) {
            for(int m = 0; m < reg.integrator->num_mergers; m++) {
                struct ketju_bh_merger_info *info = &reg.integrator->merger_infos[m];
                printf("KETJU BH MERGER [region]: m1=%g m2=%g -> m_rem=%g, chi1=%g chi2=%g -> chi_rem=%g, v_kick=%g\n",
                       info->m1, info->m2, info->m_remnant,
                       info->chi1, info->chi2, info->chi_remnant, info->v_kick);
            }
        }
    }

    /* ---- Phase 2: forward from compute root to affected root if they differ ---- */
    if(reg.compute_tasks.root_sim != reg.affected_tasks.root_sim) {
        if(ThisTask == reg.compute_tasks.root_sim) {
            MPI_Send(scatter_buf.data(), n * sizeof(ketju_scatter_particle),
                     MPI_BYTE, reg.affected_tasks.root_sim, 998, MPI_COMM_WORLD);
        }
        if(ThisTask == reg.affected_tasks.root_sim) {
            scatter_buf.resize(n);
            MPI_Recv(scatter_buf.data(), n * sizeof(ketju_scatter_particle),
                     MPI_BYTE, reg.compute_tasks.root_sim, 998, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    /* compute_root that is not in affected_tasks has finished its job (Send done) */
    if(!is_affected) return;

    /* ---- Phase 3: affected root sorts by Task, Scattervs to affected tasks ---- */

    /* build map: world task -> affected rank */
    std::unordered_map<int, int> task_to_aff_rank;
    for(size_t i = 0; i < reg.affected_sim_indices.size(); i++)
        task_to_aff_rank[reg.affected_sim_indices[i]] = (int)i;

    std::vector<int> sendcounts(n_aff, 0);
    std::vector<int> byte_sendcounts(n_aff, 0), byte_displs(n_aff, 0);

    if(reg.affected_tasks.is_root()) {
        /* sort particles by Task */
        qsort(scatter_buf.data(), n, sizeof(ketju_scatter_particle), scatter_particle_task_cmp);

        /* compute sendcounts per affected rank */
        for(int i = 0; i < n; i++) {
            auto it = task_to_aff_rank.find(scatter_buf[i].Task);
            if(it != task_to_aff_rank.end())
                sendcounts[it->second]++;
        }
        for(int t = 0; t < n_aff; t++) byte_sendcounts[t] = sendcounts[t] * sizeof(ketju_scatter_particle);
        for(int t = 1; t < n_aff; t++) byte_displs[t] = byte_displs[t-1] + byte_sendcounts[t-1];
    }

    /* scatter sendcounts so each task knows how many particles it receives */
    int my_count = 0;
    MPI_Scatter(sendcounts.data(), 1, MPI_INT, &my_count, 1, MPI_INT,
                reg.affected_tasks.root, reg.affected_tasks.comm);

    /* scatter particle data */
    std::vector<ketju_scatter_particle> my_particles(my_count);
    MPI_Scatterv(scatter_buf.data(), byte_sendcounts.data(), byte_displs.data(), MPI_BYTE,
                 my_particles.data(), my_count * sizeof(ketju_scatter_particle), MPI_BYTE,
                 reg.affected_tasks.root, reg.affected_tasks.comm);

    /* ---- Phase 4: each task applies velocity trick to its particles ---- */
    /* also collect full particle data on affected root for stellar collision check + mergers */
    std::vector<double> final_abs_pos, final_rel_vel, final_mass;
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
    std::vector<double> original_mass;
#endif

    /* affected root keeps the sorted scatter_buf for collision/merger checks */
    if(reg.affected_tasks.is_root()) {
        final_abs_pos.resize(3 * n);
        final_rel_vel.resize(3 * n);
        final_mass.resize(n);
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
        original_mass.resize(n);
#endif
        /* also rebuild extra_data and all_particles on root from scatter_buf for collision/merger code */
        reg.extra_data.resize(n);
        reg.all_particles.resize(n);
        for(int i = 0; i < n; i++) {
            final_abs_pos[3*i]   = scatter_buf[i].Pos[0];
            final_abs_pos[3*i+1] = scatter_buf[i].Pos[1];
            final_abs_pos[3*i+2] = scatter_buf[i].Pos[2];
            final_rel_vel[3*i]   = scatter_buf[i].Vel[0];
            final_rel_vel[3*i+1] = scatter_buf[i].Vel[1];
            final_rel_vel[3*i+2] = scatter_buf[i].Vel[2];
            final_mass[i] = scatter_buf[i].Mass;
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
            original_mass[i] = scatter_buf[i].OriginalMass;
#endif
            reg.extra_data[i].ID    = scatter_buf[i].ID;
            reg.extra_data[i].Task  = scatter_buf[i].Task;
            reg.extra_data[i].Index = scatter_buf[i].Index;
            reg.extra_data[i].Type  = scatter_buf[i].Type;
            reg.all_particles[i].Mass = scatter_buf[i].Mass;
        }
    }

    for(int i = 0; i < my_count; i++) {
        ketju_scatter_particle &sp = my_particles[i];
        int idx = sp.Index;

        /* sanity check */
        if(P[idx].ID != sp.ID) {
            printf("KETJU ERROR: particle index mismatch! ID %llu vs %llu on task %d\n",
                   (unsigned long long)P[idx].ID, (unsigned long long)sp.ID, ThisTask);
            endrun(667);
        }

        /* update mass */
        P[idx].Mass = sp.Mass;

        /* skip velocity trick for merged particles (mass=0, will be removed) */
        if(sp.Mass <= 0) {
            P[idx].KetjuIntegrated = 1;
            continue;
        }

        /* compute per-particle drift factor: what drift_particle() will use */
        double dt_drift = get_drift_factor(P[idx].Ti_current, P[idx].Ti_current + reg.ti_step, idx, 0);

        if(dt_drift > 0) {
            for(int j = 0; j < 3; j++) {
                /* desired final position = final_abs_pos + CoM displacement during drift */
                double desired_pos = sp.Pos[j] + reg.com_vel[j] * dt_drift;
                double delta_pos = desired_pos - P[idx].Pos[j];
#ifdef BOX_PERIODIC
                double box = (j == 0) ? boxSize_X : (j == 1) ? boxSize_Y : boxSize_Z;
                if(box > 0) { while(delta_pos > 0.5 * box) delta_pos -= box; while(delta_pos < -0.5 * box) delta_pos += box; }
#endif
                P[idx].Vel[j] = delta_pos / dt_drift;
                P[idx].KetjuFinalVel[j] = sp.Vel[j] + reg.com_vel[j];
            }
        } else {
            /* no drift — set position and velocity directly */
            for(int j = 0; j < 3; j++) {
                P[idx].Pos[j] = sp.Pos[j];
                P[idx].Vel[j] = sp.Vel[j] + reg.com_vel[j];
                P[idx].KetjuFinalVel[j] = P[idx].Vel[j];
            }
        }

#ifdef SINGLE_STAR_SINK_DYNAMICS
        /* update spin for PN particles */
        if(sp.is_pn && P[idx].Type == 5) {
            for(int j = 0; j < 3; j++) P[idx].KetjuSpin[j] = sp.Spin[j];
        }
#endif
        P[idx].KetjuIntegrated = 1;
    }

    /* ---- Phase 5: geometric coalescence check for Type 5 sinks ----
     *  All sinks have ProtoStellarRadius_inSolar populated by the protostellar
     *  evolution: protostars (stage 1-5) get their physical radius, MS stars
     *  (stage 6) get the MS radius, relics (stage 7) get 6 * 2GM/c^2 (an ISCO-
     *  like criterion). Coalescence is therefore detected when (R_i + R_j) >
     *  separation, regardless of whether the pair is star-star, star-BH, or
     *  BH-BH. handle_mergers() (gated on KETJU_MERGE_BH) accepts any Type 5
     *  pair and transfers metals to the survivor. */
#if defined(KETJU_MERGE_BH) && defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
    {
        const double R_SUN_IN_CM = 6.957e10;
        std::vector<double> my_radii(my_count, 0.0);
        for(int i = 0; i < my_count; i++) {
            if(my_particles[i].Type != 5) continue;
            int idx = my_particles[i].Index;
            if(P[idx].ProtoStellarRadius_inSolar > 0)
                my_radii[i] = P[idx].ProtoStellarRadius_inSolar * R_SUN_IN_CM;
        }

        std::vector<double> sink_radius_cm;
        std::vector<int> recv_counts(n_aff, 0), recv_displs(n_aff, 0);
        MPI_Gather(&my_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                   reg.affected_tasks.root, reg.affected_tasks.comm);
        if(reg.affected_tasks.is_root()) {
            for(int t = 1; t < n_aff; t++) recv_displs[t] = recv_displs[t-1] + recv_counts[t-1];
            sink_radius_cm.resize(n, 0.0);
        }
        MPI_Gatherv(my_radii.data(), my_count, MPI_DOUBLE,
                     sink_radius_cm.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                     reg.affected_tasks.root, reg.affected_tasks.comm);

        if(reg.affected_tasks.is_root()) {
            for(int i = 0; i < n; i++) {
                if(final_mass[i] <= 0 || reg.extra_data[i].Type != 5) continue;
                if(sink_radius_cm[i] <= 0) continue;
                for(int j = i + 1; j < n; j++) {
                    if(final_mass[j] <= 0 || reg.extra_data[j].Type != 5) continue;
                    if(sink_radius_cm[j] <= 0) continue;
                    /* If both are PN-enabled BHs AND KETJU is integrating PN terms,
                     * MSTAR's BH-merger handler (ketju_calculate_bh_merger_remnant_properties)
                     * will handle the coalescence with full GR remnant properties
                     * (mass + spin + recoil kick) at its 12·2GM/c² threshold.
                     * Skip here to avoid double-merging the same BH-BH pair.
                     * When PN is off (KetjuPNTerms="none"), MSTAR is dormant and
                     * we keep the geometric check active for BH-BH too. */
                    if(reg.extra_data[i].KetjuPNClass == 0 && reg.extra_data[j].KetjuPNClass == 0
                       && strcmp(All.KetjuPNTerms, "none") != 0) continue;
                    double d2 = 0;
                    for(int k = 0; k < 3; k++) {
                        double dx = final_abs_pos[3*i + k] - final_abs_pos[3*j + k];
                        d2 += dx * dx;
                    }
                    double sep_cm = sqrt(d2) * UNIT_LENGTH_IN_CGS * All.cf_atime;
                    if(sep_cm < sink_radius_cm[i] + sink_radius_cm[j]) {
                        int victim = (final_mass[i] < final_mass[j]) ? i : j;
                        int surv = (victim == i) ? j : i;
                        final_mass[surv] += final_mass[victim];
                        final_mass[victim] = 0;
                        double Mi = original_mass[i] * UNIT_MASS_IN_SOLAR;
                        double Mj = original_mass[j] * UNIT_MASS_IN_SOLAR;
                        printf("KETJU SINK COALESCENCE: ID %llu (%.2f Msun) + ID %llu (%.2f Msun) at sep=%.2e cm (R1+R2=%.2e cm)\n",
                            (unsigned long long)reg.extra_data[i].ID, Mi,
                            (unsigned long long)reg.extra_data[j].ID, Mj,
                            sep_cm, sink_radius_cm[i] + sink_radius_cm[j]);
                    }
                }
            }
        }
    }
#endif

    /* handle mergers (integrator mergers + stellar collisions) — merger decisions
     * are made on affected root then broadcast to all affected tasks for execution */
#if defined(KETJU_MERGE_STARS) || defined(KETJU_MERGE_BH)
    handle_mergers(reg, final_abs_pos, final_mass, original_mass);
#endif
}

/* ============================================================
 *  PUBLIC INTERFACE FUNCTIONS
 * ============================================================ */

/* ============================================================
 *  Helper: move a particle between timebin linked lists
 * ============================================================ */
static void move_particle_timebin(int i, int old_bin, int new_bin)
{
    if(old_bin == new_bin) return;

    TimeBinCount[old_bin]--;
    if(PrevInTimeBin[i] >= 0)
        NextInTimeBin[PrevInTimeBin[i]] = NextInTimeBin[i];
    else
        FirstInTimeBin[old_bin] = NextInTimeBin[i];
    if(NextInTimeBin[i] >= 0)
        PrevInTimeBin[NextInTimeBin[i]] = PrevInTimeBin[i];
    else
        LastInTimeBin[old_bin] = PrevInTimeBin[i];

    if(FirstInTimeBin[new_bin] >= 0)
        PrevInTimeBin[FirstInTimeBin[new_bin]] = i;
    NextInTimeBin[i] = FirstInTimeBin[new_bin];
    PrevInTimeBin[i] = -1;
    FirstInTimeBin[new_bin] = i;
    if(LastInTimeBin[new_bin] < 0)
        LastInTimeBin[new_bin] = i;
    TimeBinCount[new_bin]++;

    P[i].TimeBin = new_bin;
}

/* PHASE 6: Synchronize chain particles within each region.
 *
 * SUBCYCLING: instead of forcing all chain particles globally to the minimum
 * timebin (which drags the entire simulation down), we synchronize particles
 * within each region to the region's MAXIMUM active timebin. This is the
 * normal hydro/gravity step for the slowest-evolving particle in the region.
 * MSTAR handles internal substeps via its adaptive GBS integrator (η=10⁻⁷),
 * so tight binaries are resolved internally without shrinking the global step.
 *
 * All particles in a region must be on the same timebin for the velocity trick
 * to work (they all drift by the same dt after KETJU sets their velocities). */
void ketju_limit_timesteps(void)
{
    CachedChainCentersValid = 0;

    if(All.KetjuRegionRadius <= 0) return;
    if(All.KetjuMinBHMass <= 0 && All.KetjuMinStarMass <= 0) return;

    /* gather chain center positions (cache for ketju_find_regions) */
    std::vector<ketju_mpi_particle> centers = gather_chain_centers();
    CachedChainCenters = centers;
    CachedChainCentersValid = 1;
    if(centers.empty()) return;

    /* merge overlapping centers into regions */
    double merge_radius = 2.0 * All.KetjuRegionRadius;
    std::vector<std::vector<ketju_mpi_particle>> region_centers = merge_overlapping_regions(centers, merge_radius);

    int n_moved_local = 0;

    for(size_t r = 0; r < region_centers.size(); r++) {
        /* find local members of this region */
        std::set<int> local_members = find_local_members(region_centers[r], All.KetjuRegionRadius);

        /* find local MAXIMUM active timebin (largest = longest step) */
        int local_max_bin = 0;
        for(int idx : local_members) {
            if(P[idx].TimeBin > 0 && TimeBinActive[P[idx].TimeBin] && P[idx].TimeBin > local_max_bin)
                local_max_bin = P[idx].TimeBin;
        }

        /* global maximum across all tasks for this region */
        int global_max_bin;
        MPI_Allreduce(&local_max_bin, &global_max_bin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(global_max_bin <= 0) continue;

        /* move all local chain particles in this region to the max timebin
         * (particles on shorter timebins get promoted to the longer step) */
        for(int idx : local_members) {
            if(P[idx].TimeBin < global_max_bin) {
                move_particle_timebin(idx, P[idx].TimeBin, global_max_bin);
                n_moved_local++;
            }
        }
    }

    int n_moved_global;
    MPI_Reduce(&n_moved_local, &n_moved_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0 && n_moved_global > 0) {
        printf("KETJU SUBCYCLE: Promoted %d particles to longer timebins (MSTAR handles internal accuracy)\n",
               n_moved_global);
    }
}

/* Collect sorted (ID,Task) keys for a region across all tasks (for staleness check).
 * Uses (ID,Task) pairs instead of bare IDs to handle duplicate star IDs correctly. */
static std::vector<CachedRegionKey> gather_sorted_region_keys(const std::set<int> &local_members)
{
    int n_local = local_members.size();
    std::vector<CachedRegionKey> local_keys(n_local);
    int k = 0;
    for(int idx : local_members) { local_keys[k].ID = P[idx].ID; local_keys[k].Task = ThisTask; local_keys[k].Index = idx; k++; }

    std::vector<int> counts(NTask), displs(NTask);
    MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int n_total = 0;
    for(int t = 0; t < NTask; t++) { displs[t] = n_total; n_total += counts[t]; }

    std::vector<CachedRegionKey> all_keys(n_total);
    std::vector<int> byte_counts(NTask), byte_displs(NTask);
    for(int t = 0; t < NTask; t++) {
        byte_counts[t] = counts[t] * sizeof(CachedRegionKey);
        byte_displs[t] = displs[t] * sizeof(CachedRegionKey);
    }
    MPI_Allgatherv(local_keys.data(), n_local * sizeof(CachedRegionKey), MPI_BYTE,
                   all_keys.data(), byte_counts.data(), byte_displs.data(),
                   MPI_BYTE, MPI_COMM_WORLD);
    std::sort(all_keys.begin(), all_keys.end());
    return all_keys;
}

void ketju_find_regions(void)
{
    /* clear previous state — KetjuIntegrated is cleared here (not in finish_step)
     * so the flag survives across the step boundary for guard checks in kicks/predict/gravtree */
#ifdef HERMITE_INTEGRATION
    /* Mark particles that were KETJU-integrated last step as "Hermite history stale".
     * If they remain in a chain this step they never reach the stale check (in-region
     * guard fires first); if they exit, the next eligible_for_hermite() call falls
     * back to KDK once before re-engaging Hermite with fresh OldPos/Acc/Jerk. */
    for(int i = 0; i < NumPart; i++) {if(P[i].KetjuIntegrated) {P[i].HermiteHistoryStale = 1;}}
#endif
    for(int i = 0; i < NumPart; i++) {P[i].KetjuIntegrated = 0;}
    ActiveRegions.clear();
    AllKetjuParticleIndices.clear();

    /* validate parameters */
    if(All.KetjuRegionRadius <= 0) return;
    if(All.KetjuMinBHMass <= 0 && All.KetjuMinStarMass <= 0) return;

    /* reuse chain centers from ketju_limit_timesteps if available, else gather fresh */
    std::vector<ketju_mpi_particle> centers;
    if(CachedChainCentersValid) {
        centers = std::move(CachedChainCenters);
        CachedChainCentersValid = 0;
    } else {
        centers = gather_chain_centers();
    }
    if(centers.empty()) { CachedRegions.clear(); return; }

    /* merge overlapping regions */
    double merge_radius = 2.0 * All.KetjuRegionRadius;
    std::vector<std::vector<ketju_mpi_particle>> region_centers = merge_overlapping_regions(centers, merge_radius);

    /* build regions (first pass: find members, collect IDs, check staleness) */
    int n_regions = region_centers.size();
    int n_cached = CachedRegions.size();
    int any_changed = KetjuRegionsStale || (n_regions != n_cached);

    /* per-region: collect sorted (ID,Task) keys and compare to cache */
    std::vector<std::vector<CachedRegionKey>> new_keys(n_regions);
    std::vector<int> region_changed(n_regions, 1); /* 1 = needs rebuild */

    for(int r = 0; r < n_regions; r++) {
        KetjuRegion reg;
        reg.centers = region_centers[r];
        reg.local_member_indices = find_local_members(reg.centers, All.KetjuRegionRadius);

        int n_local = reg.local_member_indices.size();
        MPI_Allreduce(&n_local, &reg.total_particle_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        /* skip regions with fewer than 2 particles — nothing to integrate */
        if(reg.total_particle_count < 2) {
            ActiveRegions.push_back(std::move(reg)); /* placeholder to keep indexing aligned with cache */
            new_keys[r] = {};
            continue;
        }

        for(int idx : reg.local_member_indices)
            AllKetjuParticleIndices.insert(idx);

        ActiveRegions.push_back(std::move(reg));

        /* collect sorted (ID,Task) keys for staleness comparison */
        new_keys[r] = gather_sorted_region_keys(ActiveRegions[r].local_member_indices);

        /* check if this region matches a cached one */
        if(!any_changed && r < n_cached && new_keys[r] == CachedRegions[r].sorted_particle_keys) {
            region_changed[r] = 0; /* unchanged — can reuse communicators */
        }
    }

    /* allocate compute tasks across regions based on cost estimates */
    allocate_compute_tasks_for_regions();

    /* second pass: set up integrators — reuse cached communicators where possible */
    int n_reused = 0;
    for(int r = 0; r < n_regions; r++) {
        if(!region_changed[r]) {
            /* reuse cached communicators — transfer ownership to ActiveRegion */
            ActiveRegions[r].affected_tasks = CachedRegions[r].affected_tasks;
            CachedRegions[r].affected_tasks.comm = MPI_COMM_NULL;
            CachedRegions[r].affected_tasks.group = MPI_GROUP_NULL;
            ActiveRegions[r].compute_tasks = CachedRegions[r].compute_tasks;
            CachedRegions[r].compute_tasks.comm = MPI_COMM_NULL;
            CachedRegions[r].compute_tasks.group = MPI_GROUP_NULL;
            ActiveRegions[r].affected_sim_indices = std::move(CachedRegions[r].affected_sim_indices);
            ActiveRegions[r].particle_counts_on_affected = std::move(CachedRegions[r].particle_counts_on_affected);
            ActiveRegions[r].compute_info = CachedRegions[r].compute_info;
            /* still need to gather particles and set up integrator state, but skip MPI_Comm_create */
            setup_integrator_reuse(ActiveRegions[r]);
            n_reused++;
        } else {
            /* full rebuild */
            setup_integrator(ActiveRegions[r]);
        }
    }

    /* clear old cache (frees any remaining communicators not transferred) */
    CachedRegions.clear();
    KetjuRegionsStale = 0;

    if(ThisTask == 0 && !ActiveRegions.empty()) {
        printf("KETJU: Found %d region(s) from %d center(s) at t=%g (%d reused from cache)\n",
               n_regions, (int)centers.size(), All.Time, n_reused);
        for(int r = 0; r < n_regions; r++) {
            printf("KETJU:   region %d: %d particles, %d centers%s\n",
                   r, ActiveRegions[r].total_particle_count, (int)ActiveRegions[r].centers.size(),
                   region_changed[r] ? "" : " [cached]");
        }
    }
}

void ketju_run_integration(void)
{
    if(ActiveRegions.empty()) return;

    /* Outer loop over sequential rounds. Within a round, regions run in
     * parallel on disjoint task ranges; rounds run back-to-back so that
     * expensive regions can use all tasks first before smaller regions. */
    int n_rounds = g_ketju_num_sequential_runs > 0 ? g_ketju_num_sequential_runs : 1;
    for(int seq = 0; seq < n_rounds; seq++) {
    for(size_t r = 0; r < ActiveRegions.size(); r++) {
        KetjuRegion &reg = ActiveRegions[r];
        if(reg.compute_info.compute_sequence_position != seq) continue;
        if(reg.total_particle_count < 2) continue;

        /* Subcycling: use the MAXIMUM active timebin among chain members.
         * This is the normal hydro/gravity step for these particles.
         * MSTAR internally substeps to whatever accuracy it needs (GBS adaptive).
         * Previously we used the MINIMUM, which dragged the whole simulation down. */
        integertime local_max_ti = 0;
        for(int idx : reg.local_member_indices) {
            if(!TimeBinActive[P[idx].TimeBin]) continue;
            integertime ti_i = GET_PARTICLE_INTEGERTIME(idx);
            if(ti_i > local_max_ti) local_max_ti = ti_i;
        }

        /* global maximum across all tasks (all tasks participate since region spans tasks) */
        integertime global_max_ti;
        MPI_Allreduce(&local_max_ti, &global_max_ti, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        if(global_max_ti <= 0 || global_max_ti > TIMEBASE) continue;

        reg.ti_step = global_max_ti;

        /* convert to physical time for MSTAR integrator */
        double dt_physical;
        if(All.ComovingIntegrationOn) {
            integertime t0 = All.Ti_Current;
            integertime t1 = t0 + reg.ti_step;
            double a_mid = All.cf_atime * exp(0.5 * reg.ti_step * All.Timebase_interval);
            dt_physical = get_gravkick_factor(t0, t1, -1, 0) * a_mid;
        } else {
            dt_physical = reg.ti_step * All.Timebase_interval;
        }

        if(dt_physical <= 0 || !isfinite(dt_physical)) continue;

        /* integrate and scatter results */
        integrate_region(reg, dt_physical);
        scatter_results(reg);
    }
    } /* end for(seq) */

    /* update cost estimates for next step's load balancing */
    update_region_costs();

    if(ThisTask == 0 && !ActiveRegions.empty()) {
        int total_mergers = 0;
        for(size_t r = 0; r < ActiveRegions.size(); r++) {
            if(ActiveRegions[r].integrator && ActiveRegions[r].integrator->num_mergers > 0)
                total_mergers += ActiveRegions[r].integrator->num_mergers;
        }
        if(total_mergers > 0)
            printf("KETJU: %d merger(s) this step\n", total_mergers);
    }
}

void ketju_set_final_velocities(void)
{
    /* called after drift: swap in true physical velocities and correct dp[] momentum.
     * The velocity trick set Vel = delta_pos/dt_drift for the drift.
     * Now we replace it with the true MSTAR velocity and correct dp[]. */
    for(int i = 0; i < NumPart; i++) {
        if(P[i].KetjuIntegrated) {
            for(int j = 0; j < 3; j++) {
                /* dp correction: the difference between true velocity and trick velocity,
                 * converted to momentum. This ensures the kick integrator is consistent. */
                double dv = P[i].KetjuFinalVel[j] - P[i].Vel[j];
                P[i].dp[j] += P[i].Mass * dv;
                P[i].Vel[j] = P[i].KetjuFinalVel[j];
            }
        }
    }
}

void ketju_finish_step(void)
{
    /* NOTE: KetjuIntegrated flags are NOT cleared here — they persist until
     * ketju_find_regions() at the start of the next step, so that guard checks
     * in kicks.cc/predict.cc/gravtree.cc can see which particles were KETJU-integrated */

    /* cache region communicators for reuse next step */
    CachedRegions.clear();
    for(size_t r = 0; r < ActiveRegions.size(); r++) {
        CachedRegionInfo ci;
        /* Collect sorted (ID,Task,Index) keys for staleness check next step.
         * MUST be globally consistent on every task — using ActiveRegions[r].all_particles
         * is wrong because it is only populated on the affected_root after Gatherv,
         * leaving the cached keys empty on every other task and causing split decisions
         * (some tasks setup_integrator_reuse, others setup_integrator) → MPI deadlock. */
        ci.sorted_particle_keys = gather_sorted_region_keys(ActiveRegions[r].local_member_indices);

        /* transfer communicator ownership — avoid free in KetjuRegion destructor */
        ci.affected_tasks = ActiveRegions[r].affected_tasks;
        ActiveRegions[r].affected_tasks.comm = MPI_COMM_NULL;
        ActiveRegions[r].affected_tasks.group = MPI_GROUP_NULL;
        ci.compute_tasks = ActiveRegions[r].compute_tasks;
        ActiveRegions[r].compute_tasks.comm = MPI_COMM_NULL;
        ActiveRegions[r].compute_tasks.group = MPI_GROUP_NULL;

        ci.affected_sim_indices = std::move(ActiveRegions[r].affected_sim_indices);
        ci.particle_counts_on_affected = std::move(ActiveRegions[r].particle_counts_on_affected);
        ci.compute_info = ActiveRegions[r].compute_info;
        ci.total_particle_count = ActiveRegions[r].total_particle_count;

        CachedRegions.push_back(std::move(ci));
    }

    /* free region data (integrators, particle arrays — but comms are now in cache) */
    ActiveRegions.clear();
    AllKetjuParticleIndices.clear();
}

int ketju_is_particle_in_region(int i)
{
    return (AllKetjuParticleIndices.find(i) != AllKetjuParticleIndices.end()) ? 1 : 0;
}

void ketju_mark_regions_stale(void)
{
    KetjuRegionsStale = 1;
}

/* ============================================================
 *  Phase C: HDF5 output for KETJU diagnostics
 *  Based on public GADGET4-KETJU (Mannerkoski+ 2023).
 *  Writes per-BH trajectories, region data, and merger events.
 *  Only Task 0 writes; data collected from compute roots.
 * ============================================================ */

struct ketju_bh_output {
    double pos[3], vel[3], spin[3];
    double mass;
    int timestep_index, n_particles, n_bh;
};

struct ketju_merger_output {
    MyIDType id1, id2;
    double m1, m2, m_remnant;
    double time, redshift;
};

static hid_t ketju_output_file = -1;
static int ketju_output_tstep_index = 0;

void ketju_open_output_file(void)
{
    if(ThisTask != 0) return;
    char buf[512];
    snprintf(buf, sizeof(buf), "%sketju_output.hdf5", All.OutputDir);
    if(RestartFlag == 0) {
        ketju_output_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        /* create groups */
        hid_t grp = H5Gcreate2(ketju_output_file, "BHs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(grp);
        /* create extendable timestep dataset */
        hsize_t dims[1] = {0}, maxdims[1] = {H5S_UNLIMITED}, chunk[1] = {64};
        hid_t space = H5Screate_simple(1, dims, maxdims);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl, 1, chunk);
        H5Dcreate2(ketju_output_file, "timesteps", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Pclose(dcpl); H5Sclose(space);
        /* create extendable merger dataset */
        hid_t merger_type = H5Tcreate(H5T_COMPOUND, sizeof(ketju_merger_output));
        H5Tinsert(merger_type, "ID1", HOFFSET(ketju_merger_output, id1), H5T_NATIVE_LLONG);
        H5Tinsert(merger_type, "ID2", HOFFSET(ketju_merger_output, id2), H5T_NATIVE_LLONG);
        H5Tinsert(merger_type, "M1", HOFFSET(ketju_merger_output, m1), H5T_NATIVE_DOUBLE);
        H5Tinsert(merger_type, "M2", HOFFSET(ketju_merger_output, m2), H5T_NATIVE_DOUBLE);
        H5Tinsert(merger_type, "M_remnant", HOFFSET(ketju_merger_output, m_remnant), H5T_NATIVE_DOUBLE);
        H5Tinsert(merger_type, "time", HOFFSET(ketju_merger_output, time), H5T_NATIVE_DOUBLE);
        H5Tinsert(merger_type, "redshift", HOFFSET(ketju_merger_output, redshift), H5T_NATIVE_DOUBLE);
        dims[0] = 0;
        space = H5Screate_simple(1, dims, maxdims);
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl, 1, chunk);
        H5Dcreate2(ketju_output_file, "mergers", merger_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Pclose(dcpl); H5Sclose(space); H5Tclose(merger_type);
    } else {
        ketju_output_file = H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);
        /* read current timestep index */
        hid_t ds = H5Dopen2(ketju_output_file, "timesteps", H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t dims[1]; H5Sget_simple_extent_dims(sp, dims, NULL);
        ketju_output_tstep_index = (int)dims[0];
        H5Sclose(sp); H5Dclose(ds);
    }
    if(ThisTask == 0 && ketju_output_file >= 0)
        printf("KETJU: Opened HDF5 output file (tstep_index=%d)\n", ketju_output_tstep_index);
}

void ketju_close_output_file(void)
{
    if(ThisTask == 0 && ketju_output_file >= 0) {
        H5Fclose(ketju_output_file);
        ketju_output_file = -1;
    }
}

static void ketju_hdf5_append_double(const char *dset_name, double value)
{
    if(ThisTask != 0 || ketju_output_file < 0) return;
    hid_t ds = H5Dopen2(ketju_output_file, dset_name, H5P_DEFAULT);
    hid_t sp = H5Dget_space(ds);
    hsize_t dims[1]; H5Sget_simple_extent_dims(sp, dims, NULL);
    hsize_t newdims[1] = {dims[0] + 1};
    H5Dset_extent(ds, newdims);
    H5Sclose(sp);
    sp = H5Dget_space(ds);
    hsize_t offset[1] = {dims[0]}, count[1] = {1};
    H5Sselect_hyperslab(sp, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t memsp = H5Screate_simple(1, count, NULL);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, memsp, sp, H5P_DEFAULT, &value);
    H5Sclose(memsp); H5Sclose(sp); H5Dclose(ds);
}

static void ketju_hdf5_write_merger(ketju_merger_output *merger)
{
    if(ThisTask != 0 || ketju_output_file < 0) return;
    hid_t ds = H5Dopen2(ketju_output_file, "mergers", H5P_DEFAULT);
    hid_t sp = H5Dget_space(ds);
    hsize_t dims[1]; H5Sget_simple_extent_dims(sp, dims, NULL);
    hsize_t newdims[1] = {dims[0] + 1};
    H5Dset_extent(ds, newdims);
    H5Sclose(sp);
    sp = H5Dget_space(ds);
    hsize_t offset[1] = {dims[0]}, count[1] = {1};
    H5Sselect_hyperslab(sp, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t memsp = H5Screate_simple(1, count, NULL);
    hid_t mtype = H5Dget_type(ds);
    H5Dwrite(ds, mtype, memsp, sp, H5P_DEFAULT, merger);
    H5Tclose(mtype); H5Sclose(memsp); H5Sclose(sp); H5Dclose(ds);
}

void ketju_write_output(void)
{
    if(ActiveRegions.empty()) return;
    if(ThisTask != 0 && ketju_output_file < 0) return; /* only task 0 has the file */

    /* write timestep */
    double phys_time = All.Time; /* scale factor in cosmo, time otherwise */
    ketju_hdf5_append_double("timesteps", phys_time);
    ketju_output_tstep_index++;
}

#endif /* KETJU_REGULARIZATION */
