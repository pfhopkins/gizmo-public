// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

// Functions for constructing and working with the MST coordinate system

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

//////////////////////////////////////////////////////////
// Path finding between particles in a constructed tree //
//////////////////////////////////////////////////////////
// These functions are intended to be called through the dispatcher function
// defined in internal.h.
// See there for explanation of arguments.

// Special cases for small Nd = max_tree_distance
bool ketju_check_relative_proximity_ND_1(int v1, int v2,
                                         const struct ketju_system *R, int *d,
                                         int path[1], int sign[1]) {

    const struct graph_vertex *vertex = R->internal_state->mst->vertex;

    int diff = abs(vertex[v1].level - vertex[v2].level);

    if (diff == 1) {
        if (vertex[v2].parent == v1) {
            *d = 1;
            sign[0] = -1;
            path[0] = vertex[v2].edgetoparent;
            return true;
        } else if (vertex[v1].parent == v2) {
            *d = 1;
            sign[0] = +1;
            path[0] = vertex[v1].edgetoparent;
            return true;
        }
    }
    return false;
}

bool ketju_check_relative_proximity_ND_2(int v1, int v2,
                                         const struct ketju_system *R, int *d,
                                         int path[2], int sign[2]) {

    const struct graph_vertex *vertex = R->internal_state->mst->vertex;

    int diff = abs(vertex[v1].level - vertex[v2].level);

    if (diff == 1) {
        if (vertex[v2].parent == v1) {
            *d = 1;
            sign[0] = -1;
            path[0] = vertex[v2].edgetoparent;
            return true;
        } else if (vertex[v1].parent == v2) {
            *d = 1;
            sign[0] = +1;
            path[0] = vertex[v1].edgetoparent;
            return true;
        }
    } else if (diff == 2) {
        // Need to check that there actually is a valid parent before checking
        // the parents parent.
        if (vertex[v2].parent >= 0 && vertex[vertex[v2].parent].parent == v1) {
            *d = 2;
            sign[0] = -1;
            sign[1] = -1;
            path[0] = vertex[v2].edgetoparent;
            path[1] = vertex[vertex[v2].parent].edgetoparent;
            return true;
        } else if (vertex[v1].parent >= 0
                   && vertex[vertex[v1].parent].parent == v2) {
            *d = 2;
            sign[0] = +1;
            sign[1] = +1;
            path[0] = vertex[v1].edgetoparent;
            path[1] = vertex[vertex[v1].parent].edgetoparent;
            return 1;
        }
    } else {
        if (vertex[v1].parent == vertex[v2].parent) {
            // only possible in diff == 0 case
            *d = 2;
            sign[0] = -1;
            sign[1] = +1;
            path[0] = vertex[v2].edgetoparent;
            path[1] = vertex[v1].edgetoparent;
            return true;
        }
    }
    return false;
}

// general case
bool ketju_check_relative_proximity_general(int v1, int v2, const int Nd,
                                            const struct ketju_system *R,
                                            int *d, int path[/*Nd*/],
                                            int sign[/*Nd*/]) {

    const struct graph_vertex *vertex = R->internal_state->mst->vertex;

    // Some placeholder invalid values to make it easier to spot if the results
    // are used wrong.
    for (int i = 0; i < Nd; i++) {
        path[i] = -1000;
        sign[i] = 0;
    }

    if (vertex[v1].level == vertex[v2].level) {

        int halfNd = (int)floor(0.5 * Nd);
        int pathend[Nd];
        memset(pathend, 0, sizeof pathend);

        for (int i = 0; i < halfNd; i++) {
            sign[i] = +1;
            path[i] = vertex[v1].edgetoparent;
            pathend[i] = vertex[v2].edgetoparent;
            v1 = vertex[v1].parent;
            v2 = vertex[v2].parent;
            if (v1 == v2) {
                *d = 2 * (i + 1);
                int c = 0;
                for (int m = i + 1; m < (*d); m++) {
                    sign[m] = -1;
                    path[m] = pathend[i - c];
                    c++;
                }
                return true;
            }
        }

    } else {
        *d = abs(vertex[v1].level - vertex[v2].level);
        if (*d <= Nd) {
            int lo, hi, thissign;
            if (vertex[v1].level > vertex[v2].level) {
                lo = v2;
                hi = v1;
                thissign = +1;
            } else {
                lo = v1;
                hi = v2;
                thissign = -1;
            }
            for (int i = 0; i < *d; i++) {
                path[i] = vertex[hi].edgetoparent;
                sign[i] = thissign;
                hi = vertex[hi].parent;
            }
            if (lo == hi) { return true; }
        }
    }

    return false;
}

/////////////////////////////////////////////////////////////////////
// Minimum spanning tree (MST) construction using Prim's algorithm //
/////////////////////////////////////////////////////////////////////

// get the vertices belonging to an edge in Prim's algorithm
static void get_v1_v2(int N, int edge, int *v1, int *v2) {

    int final_element_on_row = 0;
    int elements_on_this_row;

    for (int row = 0; row < N - 1; row++) {
        if (row == 0) {
            elements_on_this_row = N - 2 - row;
        } else {
            elements_on_this_row = N - 1 - row;
        }
        final_element_on_row += elements_on_this_row;
        if (edge <= final_element_on_row) {
            *v1 = row;
            int diff = final_element_on_row - edge;
            *v2 = N - 1 - diff;
            break;
        }
    }
}

static int cmp_weight_index2(const void *a, const void *b) {
    struct weight_index_2 *a1 = (struct weight_index_2 *)a;
    struct weight_index_2 *a2 = (struct weight_index_2 *)b;
    if ((*a1).weight > (*a2).weight)
        return 1;
    else if ((*a1).weight < (*a2).weight)
        return -1;
    else
        return 0;
}

static int cmp(const void *a, const void *b) {
    struct graph_edge *a1 = (struct graph_edge *)a;
    struct graph_edge *a2 = (struct graph_edge *)b;
    if ((*a1).weight > (*a2).weight)
        return 1;
    else if ((*a1).weight < (*a2).weight)
        return -1;
    else
        return 0;
}

static int edge_search_for_prim(struct ketju_system *R, int num_local_edges) {

    timer_start(&R->perf->tprim3);
    const struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    struct weight_index wi;

    wi.weight = DBL_MAX;
    wi.index = -1;
    for (int i = 0; i < num_local_edges; i++) {
        if (mst->vertex[mst->local_edge[i].vertex1].inMST
            != mst->vertex[mst->local_edge[i].vertex2].inMST) {
            wi.weight = mst->local_edge[i].weight;
            wi.index = mst->local_edge[i].id;
            break;
        }
    }

    timer_stop(&R->perf->tprim3);
    timer_start(&R->perf->tprim4);

    MPI_Allreduce(MPI_IN_PLACE, &wi, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  R->internal_state->task_info->main_comm);

    timer_stop(&R->perf->tprim4);

    return wi.index;
}

static void construct_MST_Prim(struct ketju_system *R) {

    const struct ketju_system_physical_state *ps = R->physical_state;
    const struct task_info *ti = R->internal_state->task_info;
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;

    timer_start(&R->perf->time_mst);

    // Construct the edges handled by this task
    int block_size, cum_num_edge, istart, jstart;
    loop_scheduling_N2(R->num_particles, ti->main_group_size,
                       ti->main_group_rank, &istart, &jstart, &block_size,
                       &cum_num_edge);

    int num_local_edges = 0;
    for (int i = istart; i < R->num_particles; i++) {

        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < R->num_particles; j++) {
            if (num_local_edges >= block_size) {
                // We've processed all edges for this task.
                // If there are more tasks than edges, block_size may also be 0.
                goto after_edge_construction; // exit both loops
            }
            double w = 0;
            for (int k = 0; k < 3; k++) {
                double ds = ps->pos[i][k] - ps->pos[j][k];
                w += ds * ds;
            }
            w = sqrt(w);

            mst->local_edge[num_local_edges].weight = w;
            mst->local_edge[num_local_edges].vertex1 = i;
            mst->local_edge[num_local_edges].vertex2 = j;
            mst->local_edge[num_local_edges].id =
                (int)(cum_num_edge + num_local_edges);

            num_local_edges++;
        }
    }
after_edge_construction:

    qsort(mst->local_edge, num_local_edges, sizeof(struct graph_edge), cmp);

    // Setup vertices
    for (int i = 0; i < R->num_particles; i++) {
        mst->vertex[i].id = i;
        mst->vertex[i].degree = 0;
        mst->vertex[i].edgetoparent = -1;
        mst->vertex[i].parent = -1;
        mst->vertex[i].level = 0;
        mst->vertex[i].inMST = false;
    }

    // Find the root of the tree = the particle closest to the CoM.
    double rmin = DBL_MAX;
    int root = -1;
    double CoM_pos[3], CoM_vel[3];
    ketju_compute_CoM(R, CoM_pos, CoM_vel);

    for (int i = 0; i < R->num_particles; i++) {
        double r2 = 0;
        for (int k = 0; k < 3; k++) {
            double ds = ps->pos[i][k] - CoM_pos[k];
            r2 += ds * ds;
        }
        if (rmin > sqrt(r2)) {
            rmin = sqrt(r2);
            root = i;
        }
    }
    mst->vertex[root].inMST = true;
    for (int k = 0; k < 3; k++) {
        mst->root_pos[k] = ps->pos[root][k];
        mst->root_vel[k] = ps->vel[root][k];
    }

    // Construct the MST
    for (int c = 0; c < R->num_particles - 1; ++c) {
        int next_edge = edge_search_for_prim(R, num_local_edges);

        int v1 = -1, v2 = -1;
        get_v1_v2(R->num_particles, next_edge, &v1, &v2);

        mst->edge_in_mst[c].vertex1 = v1;
        mst->edge_in_mst[c].vertex2 = v2;

        mst->vertex[v1].degree++;
        mst->vertex[v2].degree++;

        if (!mst->vertex[v1].inMST) {

            mst->vertex[v1].parent = v2;
            mst->vertex[v1].level = mst->vertex[v2].level + 1;
            mst->vertex[v1].inMST = true;
            mst->vertex[v1].edgetoparent = c;
            for (int k = 0; k < 3; k++) {
                mst->edge_rel_pos[c][k] = ps->pos[v1][k] - ps->pos[v2][k];
                mst->edge_rel_vel[c][k] = ps->vel[v1][k] - ps->vel[v2][k];
            }
        } else {
            mst->vertex[v2].parent = v1;
            mst->vertex[v2].level = mst->vertex[v1].level + 1;
            mst->vertex[v2].inMST = true;
            mst->vertex[v2].edgetoparent = c;
            for (int k = 0; k < 3; k++) {
                mst->edge_rel_pos[c][k] = ps->pos[v2][k] - ps->pos[v1][k];
                mst->edge_rel_vel[c][k] = ps->vel[v2][k] - ps->vel[v1][k];
            }
        }
    }
    timer_stop(&R->perf->time_mst);
}

/////////////////////////////////////////////////////////////////////////
// Approximate MST consruction using oct-tree based spatial partioning //
/////////////////////////////////////////////////////////////////////////

static double get_rmax(struct ketju_system *R) {
    double rmax2 = 0;
    double CoM_pos[3], CoM_vel[3];
    ketju_compute_CoM(R, CoM_pos, CoM_vel);
    for (int i = 0; i < R->num_particles; i++) {
        double ds;
        double r2 = 0;
        for (int k = 0; k < 3; k++) {
            ds = R->physical_state->pos[i][k] - CoM_pos[k];
            r2 += ds * ds;
        }
        if (r2 > rmax2) { rmax2 = r2; }
    }
    return sqrt(rmax2);
}

static int get_octant(int triplet[3]) {
    if (triplet[0] == 1) {
        if (triplet[1] == 1) {
            if (triplet[2] == -1) {
                return 6;
            } else {
                return 7;
            }
        } else {
            if (triplet[2] == -1) {
                return 4;
            } else {
                return 5;
            }
        }
    } else {
        if (triplet[1] == 1) {
            if (triplet[2] == -1) {
                return 2;
            } else {
                return 3;
            }
        } else {
            if (triplet[2] == -1) {
                return 0;
            } else {
                return 1;
            }
        }
    }
    return -1;
}

static void get_triplet_from_octant(int octant, int *triplet) {

    if (octant == 0) {
        triplet[0] = -1;
        triplet[1] = -1;
        triplet[2] = -1;
        return;
    }
    if (octant == 1) {
        triplet[0] = -1;
        triplet[1] = -1;
        triplet[2] = +1;
        return;
    }
    if (octant == 2) {
        triplet[0] = -1;
        triplet[1] = +1;
        triplet[2] = -1;
        return;
    }
    if (octant == 3) {
        triplet[0] = -1;
        triplet[1] = +1;
        triplet[2] = +1;
        return;
    }
    if (octant == 4) {
        triplet[0] = +1;
        triplet[1] = -1;
        triplet[2] = -1;
        return;
    }
    if (octant == 5) {
        triplet[0] = +1;
        triplet[1] = -1;
        triplet[2] = +1;
        return;
    }
    if (octant == 6) {
        triplet[0] = +1;
        triplet[1] = +1;
        triplet[2] = -1;
        return;
    }
    if (octant == 7) {
        triplet[0] = +1;
        triplet[1] = +1;
        triplet[2] = +1;
        return;
    }
    return;
}

static void list_octree_level(struct octree_node *Node, int *leaf_sizes,
                              int **leaf_particles, int *Nleafs) {

    if (Node->is_leaf) {
        leaf_sizes[*Nleafs] = Node->Npart;
        leaf_particles[*Nleafs] = calloc(Node->Npart, sizeof(int));
        for (int i = 0; i < Node->Npart; i++) {
            leaf_particles[*Nleafs][i] = Node->object_list[i];
        }
        (*Nleafs)++;
    } else {
        for (int oct = 0; oct < 8; oct++) {
            if (Node->children[oct]) {
                list_octree_level(Node->children[oct], leaf_sizes,
                                  leaf_particles, Nleafs);
            }
        }
    }
}

static void get_number_of_leafs(struct octree_node *Node, int *Nleaf) {
    if (Node->is_leaf) {
        (*Nleaf)++;
    } else {
        for (int oct = 0; oct < 8; oct++) {
            if (Node->children[oct]) {
                get_number_of_leafs(Node->children[oct], Nleaf);
            }
        }
    }
}

static void free_octree(struct octree_node *Node) {
    for (int oct = 0; oct < 8; oct++) {
        if (Node->children[oct]) { free_octree(Node->children[oct]); }
    }
    free(Node->object_list);
    free(Node);
}

static void allocate_new_octree_node(struct octree_node *Node, int oct,
                                     int *label) {

    Node->children[oct] = calloc(1, sizeof(struct octree_node));
    Node->children[oct]->parent = Node;
    Node->children[oct]->half_size = Node->half_size / 2;
    Node->children[oct]->Npart = 0;
    Node->children[oct]->level = Node->level + 1;
    Node->children[oct]->is_leaf = 0;
    Node->children[oct]->ID = *label;

    (*label)++;
}

static void octree(struct ketju_system *R, struct octree_node *Node,
                   int max_leaf_size, int *label) {

    struct octree_data *octd = R->internal_state->octree_data;

    int NpartOctant[8] = {0};
    int triplet[3];
    for (int i = 0; i < Node->Npart; i++) {
        int j = Node->object_list[i];
        for (int k = 0; k < 3; k++) {
            triplet[k] =
                R->physical_state->pos[j][k] - Node->center[k] >= 0 ? 1 : -1;
        }
        int oct = get_octant(triplet);

        if (NpartOctant[oct] == 0) {
            allocate_new_octree_node(Node, oct, label);
        }

        octd->indexlist[oct][NpartOctant[oct]] = j;

        NpartOctant[oct]++;
        Node->children[oct]->Npart++;
    }

    for (int oct = 0; oct < 8; oct++) {

        if (NpartOctant[oct] == 0) { continue; }

        Node->children[oct]->object_list =
            calloc(Node->children[oct]->Npart, sizeof(int));

        for (int i = 0; i < Node->children[oct]->Npart; i++) {
            Node->children[oct]->object_list[i] = octd->indexlist[oct][i];
        }

        get_triplet_from_octant(oct, triplet);

        for (int k = 0; k < 3; k++) {
            Node->children[oct]->center[k] =
                Node->center[k] + triplet[k] * Node->children[oct]->half_size;
        }

        for (int i = 0; i < Node->children[oct]->Npart; i++) {
            octd->particle_in_node[Node->children[oct]->object_list[i]] =
                Node->children[oct];
        }

        if (Node->children[oct]->Npart == 1) {
            Node->children[oct]->is_leaf = 1;
        }
    }

    for (int oct = 0; oct < 8; oct++) {
        if (Node->children[oct]) {
            if (Node->children[oct]->Npart > max_leaf_size) {
                octree(R, Node->children[oct], max_leaf_size, label);
            } else {
                Node->children[oct]->is_leaf = 1;
            }
        }
    }
}

static void spatialpartition(struct ketju_system *R, int max_leaf_size) {

    struct octree_data *octd = R->internal_state->octree_data;

    double rmax = get_rmax(R);

    octd->root_node = calloc(1, sizeof(struct octree_node));

    octd->root_node->ID = 0;
    octd->root_node->parent = NULL;
    octd->root_node->half_size = 1.001 * rmax;
    octd->root_node->Npart = R->num_particles;
    octd->root_node->object_list = calloc(R->num_particles, sizeof(int));
    octd->root_node->level = 0;
    for (int i = 0; i < 3; i++) {
        octd->root_node->center[i] = 0.0;
    }

    for (int i = 0; i < R->num_particles; i++) {
        octd->root_node->object_list[i] = i;
    }

    int label = 1;
    octree(R, octd->root_node, max_leaf_size, &label);
}

// Get the edge index from id with the largest length (weight) given in w.
static int get_longest_edge(double w[3], int id[3]) {

    double L = -1;
    int ind = -1;
    for (int m = 0; m < 3; m++) {
        if (w[m] > L) {
            L = w[m];
            ind = m;
        }
    }
    return id[ind];
}

static int get_edge(int Nleafs, int va, int vb) {

    int N = Nleafs - 2;
    int lo, hi;
    lo = va, hi = vb;
    if (va > vb) {
        lo = vb;
        hi = va;
    }
    int index = hi - 1;
    for (int p = 0; p < lo; p++) {
        index += N;
        N--;
    }
    return index;
}

static void solve_meta_MST(struct ketju_system *R, int Nleafs, int *Npart,
                           int **particle_indices, int *num_edge) {

    const struct task_info *ti = R->internal_state->task_info;
    const struct ketju_system_physical_state *ps = R->physical_state;
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    struct octree_data *octd = R->internal_state->octree_data;

    int num_all_edges = Nleafs * (Nleafs - 1) / 2;
    struct weight_index_2 *W =
        calloc(num_all_edges, sizeof(struct weight_index_2));
    int *in_MST = calloc(Nleafs, sizeof(int));

    int cum_num_edge = 0, istart = 0;
    int ThisBlock, jstart;
    loop_scheduling_N2(Nleafs, ti->main_group_size, ti->main_group_rank,
                       &istart, &jstart, &ThisBlock, &cum_num_edge);

    struct octree_node *Node_i, *Node_j;
    int counter = 0, ind1 = -1, ind2 = -1;

    int *ind1s = calloc(num_all_edges, sizeof(int));
    int *ind2s = calloc(num_all_edges, sizeof(int));
    double *weights = calloc(num_all_edges, sizeof(double));

    for (int i = istart; i < Nleafs; i++) {

        Node_i = octd->particle_in_node[particle_indices[i][0]];

        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < Nleafs; j++) {

            if (counter >= ThisBlock) { goto after_loops; }

            Node_j = octd->particle_in_node[particle_indices[j][0]];

            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = Node_i->center[k] - Node_j->center[k];
                r2 = ds * ds;
            }
            const double node_separation = sqrt(r2);
            const double node_size =
                1.01 * sqrt(2)
                * (Node_i->half_size + Node_j->half_size); // with safety buffer

            if (node_separation <= node_size) {

                double min_weight = DBL_MAX;
                for (int q = 0; q < Npart[i]; q++) {
                    for (int p = 0; p < Npart[j]; p++) {
                        double r2 = 0;
                        int i1 = particle_indices[i][q];
                        int i2 = particle_indices[j][p];
                        for (int k = 0; k < 3; k++) {
                            double ds = ps->pos[i1][k] - ps->pos[i2][k];
                            r2 += ds * ds;
                        }
                        double r = sqrt(r2);
                        if (r < min_weight) {
                            min_weight = r;
                            ind1 = i1;
                            ind2 = i2;
                        }
                    }
                }

            } else {

                double min_weight = DBL_MAX;
                for (int q = 0; q < Npart[i]; q++) {
                    int i1 = particle_indices[i][q];
                    double r2 = 0;
                    for (int k = 0; k < 3; k++) {
                        double ds = ps->pos[i1][k] - Node_j->center[k];
                        r2 += ds * ds;
                    }
                    double r = sqrt(r2);
                    if (r < min_weight) {
                        min_weight = r;
                        ind1 = i1;
                    }
                }
                min_weight = DBL_MAX;
                for (int p = 0; p < Npart[j]; p++) {
                    int i2 = particle_indices[j][p];
                    double r2 = 0;
                    for (int k = 0; k < 3; k++) {
                        double ds = ps->pos[i2][k] - Node_i->center[k];
                        r2 += ds * ds;
                    }
                    double r = sqrt(r2);
                    if (r < min_weight) {
                        min_weight = r;
                        ind2 = i2;
                    }
                }
            }

            r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = ps->pos[ind1][k] - ps->pos[ind2][k];
                r2 += ds * ds;
            }
            double r = sqrt(r2);

            weights[cum_num_edge + counter] = r;
            ind1s[cum_num_edge + counter] = ind1;
            ind2s[cum_num_edge + counter] = ind2;
            counter++;
        }
    }
after_loops:;
    int *recvcounts = calloc(ti->main_group_size, sizeof(int));
    int *displs = calloc(ti->main_group_size, sizeof(int));
    for (int i = 0; i < ti->main_group_size; ++i) {
        int dummy, dummy2;
        loop_scheduling_N2(Nleafs, ti->main_group_size, i, &dummy, &dummy2,
                           &recvcounts[i], &displs[i]);
    }

    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ind1s, recvcounts,
                   displs, MPI_INT, ti->main_comm);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ind2s, recvcounts,
                   displs, MPI_INT, ti->main_comm);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, weights, recvcounts,
                   displs, MPI_DOUBLE, ti->main_comm);

    free(recvcounts);
    free(displs);

    counter = 0;
    for (int i = 0; i < Nleafs; i++) {
        for (int j = i + 1; j < Nleafs; j++) {
            W[counter].weight = weights[counter];
            W[counter].id1 = ind1s[counter];
            W[counter].id2 = ind2s[counter];
            W[counter].index1 = i;
            W[counter].index2 = j;
            counter++;
        }
    }

    free(weights);
    free(ind1s);
    free(ind2s);

    // Filter out edges that clearly aren't part of the final tree, these get
    // the corresponding entry in the flag table set
    istart = 0;
    int istop = num_all_edges;

    int block = (int)floor(num_all_edges / ti->main_group_size);
    istart = ti->main_group_rank * block;
    istop = (ti->main_group_rank + 1) * block;
    if (ti->main_group_rank == ti->main_group_size - 1) istop = num_all_edges;

    bool *flag = calloc(num_all_edges, sizeof(bool));

    int Noff = 0;
    double w[3];
    int ind[3];
    for (int e1 = istart; e1 < istop; e1++) {

        int v1 = W[e1].index1;
        int v2 = W[e1].index2;
        ind[0] = e1;
        w[0] = W[e1].weight;

        if (w[0] < 0) continue;

        for (int v3 = 0; v3 < Nleafs; v3++) {

            if (w[0] < 0) continue;

            if (v3 == v1 || v3 == v2) continue;

            ind[1] = get_edge(Nleafs, v1, v3);
            ind[2] = get_edge(Nleafs, v2, v3);

            w[1] = W[ind[1]].weight;
            w[2] = W[ind[2]].weight;
            if (w[1] < 0 || w[2] < 0) continue;
            int elong = get_longest_edge(w, ind);

            flag[elong] = true;
            W[elong].weight *= -1;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, flag, num_all_edges, MPI_C_BOOL, MPI_LOR,
                  ti->main_comm);

    // How many are filtered
    for (int e1 = 0; e1 < num_all_edges; e1++) {
        if (flag[e1]) Noff++;
    }

    // Collect the remaining edges
    struct weight_index_2 *Wfilter =
        calloc(num_all_edges - Noff, sizeof(struct weight_index_2));

    counter = 0;
    for (int e1 = 0; e1 < num_all_edges; e1++) {
        if (flag[e1] == 0) {
            Wfilter[counter].weight = W[e1].weight;
            Wfilter[counter].id1 = W[e1].id1;
            Wfilter[counter].id2 = W[e1].id2;
            Wfilter[counter].index1 = W[e1].index1;
            Wfilter[counter].index2 = W[e1].index2;
            counter++;
        }
    }

    num_all_edges -= Noff;
    free(W);
    free(flag);

    qsort(Wfilter, num_all_edges, sizeof(struct weight_index_2),
          cmp_weight_index2);

    in_MST[0] = 1;

    int remaining = Nleafs - 1;
    istart = 0;
    istop = num_all_edges;

    block = (int)floor(num_all_edges / ti->main_group_size);
    istart = ti->main_group_rank * block;
    istop = (ti->main_group_rank + 1) * block;
    if (ti->main_group_rank == ti->main_group_size - 1) istop = num_all_edges;

    do {
        struct weight_index wi;
        wi.weight = DBL_MAX;
        wi.index = -1;

        for (int i = istart; i < istop; i++) {
            int i1 = Wfilter[i].index1;
            int i2 = Wfilter[i].index2;
            if (in_MST[i1] != in_MST[i2]) {
                wi.index = i;
                wi.weight = Wfilter[i].weight;
                break;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &wi, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                      ti->main_comm);

        in_MST[Wfilter[wi.index].index1] = 1;
        in_MST[Wfilter[wi.index].index2] = 1;

        mst->edge_in_mst[*num_edge].vertex1 = Wfilter[wi.index].id1;
        mst->edge_in_mst[*num_edge].vertex2 = Wfilter[wi.index].id2;
        mst->edge_in_mst[*num_edge].weight = Wfilter[wi.index].weight;
        (*num_edge)++;

        remaining--;
    } while (remaining > 0);

    free(Wfilter);
    free(in_MST);
}

static void solve_sub_MST(struct ketju_system *R, int Npart, int *Particles,
                          int *num_edge) {
    if (Npart == 1) return;

    const struct ketju_system_physical_state *ps = R->physical_state;
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;

    int num_all_edges = Npart * (Npart - 1) / 2;
    struct weight_index_2 *W =
        calloc(num_all_edges, sizeof(struct weight_index_2));
    int *in_MST = calloc(Npart, sizeof(int));

    int counter = 0;
    for (int i = 0; i < Npart; i++) {
        int ind1 = Particles[i];
        for (int j = i + 1; j < Npart; j++) {
            int ind2 = Particles[j];
            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = ps->pos[ind1][k] - ps->pos[ind2][k];
                r2 += ds * ds;
            }
            double r = sqrt(r2);
            W[counter].weight = r;
            W[counter].index1 = i;
            W[counter].index2 = j;
            W[counter].id1 = ind1;
            W[counter].id2 = ind2;
            counter++;
        }
    }

    qsort(W, num_all_edges, sizeof(struct weight_index_2), cmp_weight_index2);

    in_MST[0] = 1;

    int remaining = Npart - 1;

    do {

        for (int i = 0; i < num_all_edges; i++) {
            int i1 = W[i].index1;
            int i2 = W[i].index2;
            if (in_MST[i1] != in_MST[i2]) {
                if (in_MST[i1] == 1) {
                    in_MST[i2] = 1;
                } else {
                    in_MST[i1] = 1;
                }

                mst->edge_in_mst[*num_edge].vertex1 = W[i].id1;
                mst->edge_in_mst[*num_edge].vertex2 = W[i].id2;
                mst->edge_in_mst[*num_edge].weight = W[i].weight;
                (*num_edge)++;
                break;
            }
        }

        remaining--;
    } while (remaining > 0);

    free(W);
    free(in_MST);
}

static void swap_edge(struct ketju_system *R, int i, int j) {
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;

    int v1 = mst->edge_in_mst[i].vertex1;
    int v2 = mst->edge_in_mst[i].vertex2;
    double w = mst->edge_in_mst[i].weight;

    mst->edge_in_mst[i].vertex1 = mst->edge_in_mst[j].vertex1;
    mst->edge_in_mst[i].vertex2 = mst->edge_in_mst[j].vertex2;
    mst->edge_in_mst[i].weight = mst->edge_in_mst[j].weight;

    mst->edge_in_mst[j].vertex1 = v1;
    mst->edge_in_mst[j].vertex2 = v2;
    mst->edge_in_mst[j].weight = w;
}

static void combine_all_MST(struct ketju_system *R, int num_edge) {

    const struct ketju_system_physical_state *ps = R->physical_state;
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;

    int root = -1;
    double rmin = DBL_MAX;
    double CoM_pos[3], CoM_vel[3];
    ketju_compute_CoM(R, CoM_pos, CoM_vel);

    for (int i = 0; i < R->num_particles; i++) {
        double r2 = 0;
        for (int k = 0; k < 3; k++) {
            double ds = ps->pos[i][k] - CoM_pos[k];
            r2 += ds * ds;
        }
        double r = sqrt(r2);
        if (rmin > r) {
            rmin = r;
            root = i;
        }
    }

    int counter = 0;
    mst->vertex[root].inMST = true;
    for (int k = 0; k < 3; k++) {
        mst->root_pos[k] = ps->pos[root][k];
        mst->root_vel[k] = ps->vel[root][k];
    }

    int istop = num_edge;

    do {

        struct weight_index wi;
        wi.weight = DBL_MAX;
        wi.index = -1;

        for (int i = counter; i < istop; i++) {
            int v1 = mst->edge_in_mst[i].vertex1;
            int v2 = mst->edge_in_mst[i].vertex2;
            if (mst->vertex[v1].inMST != mst->vertex[v2].inMST) {
                if (mst->edge_in_mst[i].weight < wi.weight) {
                    wi.weight = mst->edge_in_mst[i].weight;
                    wi.index = i;
                }
            }
        }

        int new = mst->edge_in_mst[wi.index].vertex1;
        int old = mst->edge_in_mst[wi.index].vertex2;
        if (mst->vertex[new].inMST == 1) {
            old = mst->edge_in_mst[wi.index].vertex1;
            new = mst->edge_in_mst[wi.index].vertex2;
        }

        mst->vertex[new].inMST = 1;
        mst->vertex[new].parent = old;
        mst->vertex[new].level = mst->vertex[old].level + 1;
        mst->vertex[new].edgetoparent = counter;

        for (int k = 0; k < 3; k++) {
            mst->edge_rel_pos[counter][k] = ps->pos[new][k] - ps->pos[old][k];
            mst->edge_rel_vel[counter][k] = ps->vel[new][k] - ps->vel[old][k];
        }

        swap_edge(R, counter, wi.index);
        counter++;

    } while (counter < num_edge);
}

static void construct_MST_divide_and_conquer_Prim(struct ketju_system *R) {

    timer_start(&R->perf->time_mst);
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    struct octree_data *octd = R->internal_state->octree_data;

    for (int i = 0; i < R->num_particles; i++) {
        mst->vertex[i].id = i;
        mst->vertex[i].degree = 0;
        mst->vertex[i].edgetoparent = -1;
        mst->vertex[i].parent = -1;
        mst->vertex[i].level = 0;
        mst->vertex[i].inMST = 0;
    }

    int Nleafs = 0;
    int max_leaf_size = (int)floor(3 * sqrt(R->num_particles));

    timer_start(&R->perf->tprim1);
    spatialpartition(&R[0], max_leaf_size);
    get_number_of_leafs(octd->root_node, &Nleafs);
    int *leaf_sizes = calloc(Nleafs, sizeof(int));
    int **leaf_particles = calloc(Nleafs, sizeof(int *));
    timer_stop(&R->perf->tprim1);

    timer_start(&R->perf->tprim2);
    Nleafs = 0;
    list_octree_level(octd->root_node, leaf_sizes, leaf_particles, &Nleafs);
    timer_stop(&R->perf->tprim2);

    timer_start(&R->perf->tprim3);
    int num_edge = 0;
    for (int i = 0; i < Nleafs; i++) {
        solve_sub_MST(R, leaf_sizes[i], leaf_particles[i], &num_edge);
    }
    timer_stop(&R->perf->tprim3);

    timer_start(&R->perf->tprim4);
    solve_meta_MST(R, Nleafs, leaf_sizes, leaf_particles, &num_edge);
    timer_stop(&R->perf->tprim4);

    timer_start(&R->perf->tprim5);
    combine_all_MST(R, num_edge);
    timer_stop(&R->perf->tprim5);

    free_octree(octd->root_node);
    for (int i = 0; i < Nleafs; i++) {
        free(leaf_particles[i]);
    }
    free(leaf_particles);
    free(leaf_sizes);

    timer_stop(&R->perf->time_mst);
}

// main construction function

void ketju_construct_MST(struct ketju_system *R) {
    // Use the basic prim method for small systems always, since divide and
    // conquer is not useful for small systems.
    if (R->options->mst_algorithm == KETJU_MST_PRIM || R->num_particles < 10) {
        construct_MST_Prim(R);
    } else if (R->options->mst_algorithm == KETJU_MST_DIVIDE_AND_CONQUER_PRIM) {
        construct_MST_divide_and_conquer_Prim(R);
    } else {
        char msg[100];
        sprintf(msg, "Invalid mst_algorithm value %d!",
                R->options->mst_algorithm);
        terminate(msg);
    }
}
