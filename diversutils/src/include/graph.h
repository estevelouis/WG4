/*
 *      DiversUtils - Functions to measure diversity
 *
 * Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <pthread.h>
#include <stdint.h>

#ifndef ENABLE_AVX256
#define ENABLE_AVX256 0
#endif

#ifndef MST_IMPLEMENTATION_VERSION
#define MST_IMPLEMENTATION_VERSION 3
#endif

// #define GRAPH_CAPACITY_STEP 64
#define GRAPH_CAPACITY_STEP 4096 / sizeof(struct graph_node)

#ifndef FP_MODES
#define FP_MODES
enum {
	FP32,
	FP64
};
#endif

struct graph_neighbour {
	// uint64_t index;
	uint32_t index;
	float distance;
};

struct graph_node {
	struct word2vec_entry* word2vec_entry_pointer;
	struct graph_neighbour* neighbours;
	union {
		float* fp32;
		double* fp64;
	} vector;
	pthread_mutex_t mutex_local_node;
	double relative_proportion;
	uint32_t absolute_proportion;
	uint32_t capacity_neighbours;
	uint32_t num_neighbours;
	uint16_t num_dimensions;
	uint8_t already_considered;
};

struct matrix {
	union {
		float* fp32;
		double* fp64;
	} bfr;
	uint8_t* active;
	uint8_t* active_final;
	uint32_t a;
	uint32_t b;
	uint8_t fp_mode;
};

struct graph {
	struct graph_node* nodes;
	uint64_t num_nodes;
	uint64_t capacity;
	pthread_mutex_t mutex_nodes;
    pthread_mutex_t mutex_matrix;
	struct matrix dist_mat;
	int16_t num_dimensions;
	uint8_t dist_mat_must_be_freed;
};

// ---- <legacy> ----

#define MAX_ABUNDANCE +1000

#define WORD2VEC_KEY_BUFFER_SIZE 64

#define RANDOM_MODULO 65536
#define RANDOM_MIN_VALUE -128.0
#define RANDOM_MAX_VALUE +128.0

// ---- </legacy> ----

// ---- <threading> ----

struct row_thread_arg {
	float* vector;
	const struct graph* g;
	uint64_t i;
	uint64_t start_j;
	uint64_t end_j;
};

void* row_thread(void*);

struct batch_row_thread_arg {
	float* vector;
	const struct graph* g;
	uint64_t i;
	uint64_t start_j;
	uint64_t end_j;
};

struct matrix_thread_arg {
	struct matrix* m;
	struct graph* g;
	uint8_t thread_rank;
	uint8_t thread_total_count;
};

void* matrix_thread(void*);

// ---- </threading> ----


// ---- <matrix> ----

int32_t create_matrix(struct matrix* const, const uint32_t, const uint32_t, const int8_t);
int32_t distance_matrix_from_graph(const struct graph* const restrict, struct matrix* const restrict);
void distance_row_from_graph(const struct graph* const restrict, const int32_t, float* const restrict);
int32_t distance_row_from_graph_multithread(const struct graph* const, const uint64_t, float* const, const int16_t);
int32_t distance_row_batch_from_graph_multithread(const struct graph* const, const uint64_t, float* const, const int16_t, int16_t);
int32_t distance_matrix_from_graph_multithread(struct graph* const, struct matrix* const, const int16_t);
void free_matrix(struct matrix*);
int32_t stats_matrix(struct matrix*, double*, double*, double*, double*);

// ---- </matrix> ----

// ---- <other> ----

int32_t void_strcmp(const void*, const void*);

struct distance_two_nodes {
	// struct graph_node* a;
	// struct graph_node* b;
	uint32_t a;
	uint32_t b;
	// double distance;
	// int32_t usable;
	// int64_t usable;
	float distance;
	// int32_t usable;
};

// void create_distance_two_nodes(struct distance_two_nodes* restrict, struct graph_node* restrict, struct graph_node* restrict, const int8_t, const float* const);
void create_distance_two_nodes(struct distance_two_nodes* restrict distance, uint32_t a, uint32_t b, const float distance_value);
int32_t double_cmp(const void*, const void*);
// ---- </other> ----

// ---- <graph> ----

// #define GRAPH_NODE_NEIGHBOUR_STEP 64;
#define GRAPH_NODE_NEIGHBOUR_STEP 4096 / sizeof(struct graph_neighbour)

#ifndef INITIALIZE_GRAPH_NODE_WITH_RANDOM_VECTOR
#define INITIALIZE_GRAPH_NODE_WITH_RANDOM_VECTOR 0
#endif

enum {
	GRAPH_NODE_FP32,
	GRAPH_NODE_FP64
};

int32_t create_graph_node(struct graph_node* restrict, const uint16_t, const uint8_t);
int32_t request_more_neighbour_capacity_graph_node(struct graph_node*);
void free_graph_node(struct graph_node* restrict, const int8_t);
int32_t create_graph(struct graph* restrict const, const int32_t, const int16_t, const int8_t);
int32_t request_more_capacity_graph(struct graph* restrict const);
int32_t create_graph_empty(struct graph* restrict const);
void compute_graph_relative_proportions(struct graph* const);
int32_t compute_graph_dist_mat(struct graph* const, const int16_t);
// ---- </graph> ----

// ---- <word2vec> ----

struct word2vec_entry {
	float* vector;
	struct graph_node* graph_node_pointer;
	uint64_t num_occurrences;
	uint64_t graph_node_index;
    pthread_mutex_t mutex;
	char key[WORD2VEC_KEY_BUFFER_SIZE];
	uint8_t active_in_current_graph;
};

struct word2vec {
	float* vectors;
	struct word2vec_entry* keys;
	uint64_t num_vectors;
	uint16_t num_dimensions;
};

int32_t word2vec_entry_cmp(const void* restrict, const void* restrict);
int32_t load_word2vec_binary(struct word2vec* restrict, const char* restrict);
void reset_word2vec_active_in_current_graph(struct word2vec* restrict const w2v);
void free_word2vec(struct word2vec* restrict);
int32_t word2vec_key_to_index(const struct word2vec* restrict, const char* restrict);
struct word2vec_entry* word2vec_find_closest(const struct word2vec* restrict, const char* restrict);
void free_graph(struct graph* restrict);

int32_t word2vec_to_graph_fp32(struct graph*, struct word2vec*, char**, char**, int32_t, int32_t, const char * const);
// ---- </word2vec> ----

// ---- <heap> ----

struct graph_distance_heap {
	struct graph* g;
	struct distance_two_nodes* distances;
	uint64_t num_distances;
	uint64_t num_distances_memory;
};

enum {
	MIN_HEAP,
	MAX_HEAP
};

void siftdown_min_heap(struct graph_distance_heap* const restrict, const uint64_t);
int32_t siftdown(struct graph_distance_heap* const restrict, const int32_t, const int64_t);
// void heapify_min_heap(struct graph_distance_heap* restrict, const uint64_t);
void heapify_min_heap(struct graph_distance_heap* const restrict);
int32_t heapify(struct graph_distance_heap* restrict, const int32_t);
int32_t create_graph_distance_heap(struct graph_distance_heap* restrict, struct graph* restrict, const struct matrix* const);
void pop_graph_distance_min_heap(struct graph_distance_heap* restrict, const uint64_t);
void pop_graph_distance_min_heap_v2(struct graph_distance_heap* restrict, const uint64_t);
void pop_graph_distance_min_heap_v3(struct graph_distance_heap* restrict, const uint64_t);
int32_t pop_graph_distance_heap(struct graph_distance_heap*, const int32_t, const int32_t);
void free_graph_distance_heap(struct graph_distance_heap*);
// ---- </heap> ----

// ---- <minimum_spanning_tree> ----
struct minimum_spanning_tree {
	struct graph_distance_heap* heap;
	struct distance_two_nodes* distances;
	struct graph_node** nodes;
	uint64_t num_distances;
	uint64_t num_nodes;
	uint64_t num_active_distances;
	uint64_t num_active_nodes;
};

enum {
	MST_PRIMS_ALGORITHM
};

int32_t create_minimum_spanning_tree(struct minimum_spanning_tree*, struct graph_distance_heap*);
void free_minimum_spanning_tree(struct minimum_spanning_tree*);
int32_t find_minimum_acceptable_arc(struct minimum_spanning_tree*, uint64_t, double, int32_t);
int32_t calculate_minimum_spanning_tree(struct minimum_spanning_tree*, struct matrix*, int32_t);
// ---- </minimum_spanning_tree> ----

// ---- <disparities> ----
int32_t agg_mst_from_minimum_spanning_tree(struct minimum_spanning_tree*, double*);
int32_t functional_evenness_from_minimum_spanning_tree(struct minimum_spanning_tree*, double*);
int32_t functional_dispersion_from_graph(struct graph* const, double* const, const int8_t);
int32_t functional_divergence_modified_from_graph(struct graph* const, double* const, const int8_t);
int32_t pairwise_from_graph(struct graph* const, double* const, const int8_t, const struct matrix* const);
int32_t _weitzman(struct matrix*, double*);
int32_t weitzman_from_graph(struct graph* const, double* const, const int8_t);
int32_t _lexicographic(const struct matrix* const, double* const, long double* const);
int32_t lexicographic_from_graph(struct graph* const, double* const, long double* const, const int8_t, const struct matrix* const);
int32_t stirling_from_graph(struct graph*, double* restrict const, const double, const double, const int8_t, const struct matrix* restrict const m_);
int32_t ricotta_szeidl_from_graph(struct graph* const, double* const, const double, const int8_t, const struct matrix* const);
int32_t chao_et_al_functional_diversity_from_graph(struct graph* const, double* const, double* const, const double, const int8_t, const struct matrix* const);
int32_t leinster_cobbold_diversity_from_graph(struct graph* const, double* const, double* const, const double, const int8_t, const struct matrix* const);
int32_t scheiner_species_phylogenetic_functional_diversity_from_graph(struct graph* const, double* const, double* const, const double, const int8_t, const struct matrix* const);
int32_t nhc_e_q_from_graph(struct graph* const, double* const res_nhc, double* const res_e_q);
// ---- </disparities> ----

// ---- <iterative_disparities> ----
struct iterative_state_pairwise_from_graph {
	int64_t n;
	double result;
	struct graph* g;
	pthread_mutex_t mutex;
	int32_t i;
};

struct iterative_state_stirling_from_graph {
	int64_t n;
	double alpha;
	double beta;
	double result;
	struct graph* g;
	pthread_mutex_t mutex;
	int32_t i;
};

struct iterative_state_leinster_cobbold_from_graph {
	int64_t n;
	double alpha;
	double hill_number;
	double entropy;
	struct graph* g;
	pthread_mutex_t mutex;
	int32_t i;
};

struct thread_args_aggregator {
	union {
		struct iterative_state_pairwise_from_graph* const pairwise;
		struct iterative_state_stirling_from_graph* const stirling;
		struct iterative_state_leinster_cobbold_from_graph* const leinster_cobbold;
	} iter_state;
	int32_t i;
	const float* const vector;
};

int32_t create_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph* const restrict, struct graph* const g);
void iterate_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph* const restrict, const float* const);
void* iterate_iterative_state_pairwise_from_graph_thread(void*);

#if ENABLE_AVX256 == 1
void* iterate_iterative_state_pairwise_from_graph_avx_thread(void*);
#endif
void finalise_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph* const restrict);

int32_t create_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph* const restrict, struct graph* const, double, double);
void iterate_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph* const restrict, const float* const);
void* iterate_iterative_state_stirling_from_graph_thread(void*);
void finalise_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph* const restrict);

int32_t create_iterative_state_leinster_cobbold_from_graph(struct iterative_state_leinster_cobbold_from_graph* const restrict, struct graph* const, double);
void iterate_iterative_state_leinster_cobbold_from_graph(struct iterative_state_leinster_cobbold_from_graph* const restrict, const float* const);
void* iterate_iterative_state_leinster_cobbold_from_graph_thread(void*);
void finalise_iterative_state_leinster_cobbold_from_graph(struct iterative_state_leinster_cobbold_from_graph* const restrict);

// ---- </iterative_disparities> ----


#endif
