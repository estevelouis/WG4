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

#include <errno.h>
#include <immintrin.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "graph.h"

#include "cupt/parser.h"
#include "distances.h"
#include "general_constants.h"
#include "logging.h"
#include "stats.h"

const int32_t CONSTANT_RELATIVE_PROPORTION = 1;

const int32_t ENABLE_EXTREMES = 1;
const int32_t EXTREME_STEP = 50;
const double EXTREME_RATIO = 4.0;

int32_t create_matrix(struct matrix *const m, const uint32_t a, const uint32_t b, const int8_t fp_mode) {
  if (!(fp_mode == FP32 || fp_mode == FP64)) {
    perror("Unknown fp_mode in create_matrix\n");
    return 1;
  }

  m->a = a;
  m->b = b;
  size_t malloc_size = ((size_t)a) * ((size_t)b);
  switch (fp_mode) {
  case FP32:
    malloc_size *= sizeof(float);
    break;
  case FP64:
    malloc_size *= sizeof(double);
    break;
  }

  void *malloc_pointer;

  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_failure;
  }
  memset(malloc_pointer, '\0', malloc_size);
  switch (fp_mode) {
  case FP32:
    m->bfr.fp32 = (float *)malloc_pointer;
    break;
  case FP64:
    m->bfr.fp64 = (double *)malloc_pointer;
    break;
  }

  malloc_size = a * b * sizeof(int8_t);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    switch (fp_mode) {
    case FP32:
      free(m->bfr.fp32);
      break;
    case FP64:
      free(m->bfr.fp64);
      break;
    }
    goto malloc_failure;
  }
  memset(malloc_pointer, '\0', malloc_size);
  m->active = (uint8_t *)malloc_pointer;

  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    switch (fp_mode) {
    case FP32:
      free(m->bfr.fp32);
      break;
    case FP64:
      free(m->bfr.fp64);
      break;
    }
    goto malloc_failure;
  }
  memset(malloc_pointer, '\0', malloc_size);
  m->active_final = (uint8_t *)malloc_pointer;

  m->fp_mode = fp_mode;

  // m->to_free = 1;

  return 0;

malloc_failure:
  fprintf(stderr, "Failed to malloc %li bytes\n", malloc_size);
  return 1;
}

int32_t distance_matrix_from_graph(const struct graph *const restrict g, struct matrix *const restrict m) {
  if (m->a != m->b) {
    perror("m->a != m->b\n");
    return 1;
  }
  if (((uint32_t)g->num_nodes) != m->a) {
    perror("g->num_nodes != m->a\n");
    return 1;
  }

  for (uint64_t i = 0; i < m->a; i++) {
    switch (m->fp_mode) {
    case FP32:
      m->bfr.fp32[i * m->b + i] = 0.0f;
      break;
    case FP64:
      m->bfr.fp64[i * m->b + i] = 0.0;
      break;
    }
    for (uint64_t j = i + 1; j < m->b; j++) {
      switch (m->fp_mode) {
      case FP32:
#if ENABLE_AVX256 == 1
        m->bfr.fp32[i * m->b + j] =
            (float)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
        m->bfr.fp32[i * m->b + j] =
            (float)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
        m->bfr.fp32[j * m->b + i] = m->bfr.fp32[i * m->b + j];
        break;
      case FP64:
        m->bfr.fp64[i * m->b + j] =
            (double)cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
        m->bfr.fp64[j * m->b + i] = m->bfr.fp64[i * m->b + j];
        break;
      }
    }
  }

  return 0;
}

void *row_thread(void *args) {
  float *vector = (((struct row_thread_arg *)args)->vector);
  const struct graph *g = (((struct row_thread_arg *)args)->g);
  uint64_t i = ((struct row_thread_arg *)args)->i;
  uint64_t start_j = ((struct row_thread_arg *)args)->start_j;
  uint64_t end_j = ((struct row_thread_arg *)args)->end_j;

  for (uint64_t j = start_j; j < end_j; j++) {
#if ENABLE_AVX256 == 1
    vector[j] = cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
    vector[j] = cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
  }
  return NULL;
}

void distance_row_from_graph(const struct graph *const restrict g, const int32_t i, float *const restrict vector) {
  for (uint64_t j = 0; j < g->num_nodes; j++) {
#if ENABLE_AVX256 == 1
    vector[j] = cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
    vector[j] = cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
  }
}

int32_t distance_row_from_graph_multithread(const struct graph *const g, const uint64_t i, float *const vector,
                                            const int16_t num_row_threads) {
  pthread_t threads[num_row_threads];
  struct row_thread_arg args[num_row_threads];
  uint64_t start_j = 0;
  uint64_t end_j;
  for (int16_t k = 0; k < num_row_threads; k++) {
    end_j = start_j + floor(g->num_nodes / num_row_threads);
    if (((uint64_t)k) < g->num_nodes % ((uint64_t)num_row_threads)) {
      end_j++;
    }
    args[k].i = i;
    args[k].start_j = start_j;
    args[k].end_j = end_j;
    args[k].vector = vector;
    args[k].g = g;
    if (pthread_create(&(threads[k]), NULL, row_thread, &(args[k])) != 0) {
      perror("failed to create row thread\n");
      return 1;
    }
    start_j = end_j;
  }

  for (int16_t a = 0; a < num_row_threads; a++) {
    if (pthread_join(threads[a], NULL) != 0) {
      perror("failed to join matrix thread\n");
      return 1;
    }
  }

  return 0;
}

void *batch_row_thread(void *args) {
  float *vector = (((struct row_thread_arg *)args)->vector);
  const struct graph *g = (((struct row_thread_arg *)args)->g);
  uint64_t i = ((struct row_thread_arg *)args)->i;
  uint64_t start_j = ((struct row_thread_arg *)args)->start_j;
  uint64_t end_j = ((struct row_thread_arg *)args)->end_j;

  for (uint64_t j = start_j; j < end_j; j++) {
#if ENABLE_AVX256 == 1
    vector[j] = cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
    vector[j] = cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
  }
  return NULL;
}

int32_t distance_row_batch_from_graph_multithread(const struct graph *const g, const uint64_t i, float *const vector,
                                                  const int16_t num_threads, int16_t batch_size) {
  pthread_t threads[num_threads];
  struct batch_row_thread_arg args[num_threads];

  if (batch_size > num_threads) {
    perror("In distance_row_batch_from_graph_multithread, batch_size must be lower or equal to the number of threads\n");
    return 1;
  }

  if (i + batch_size >= g->num_nodes) {
    batch_size = g->num_nodes - i;
  }

  int16_t thread_index = 0;
  const int16_t step_threads_to_launch = floor(num_threads / batch_size);
  for (int16_t a = 0; a < batch_size; a++) {
    int16_t num_threads_to_launch = step_threads_to_launch;
    if (a < num_threads % batch_size) {
      num_threads_to_launch++;
    }
    uint64_t start_j = 0;
    uint64_t end_j;
    const uint64_t step_j = floor(g->num_nodes / num_threads_to_launch);
    for (int16_t k = 0; k < num_threads_to_launch; k++) {
      end_j = start_j + step_j;
      if (((uint64_t)k) < g->num_nodes % ((uint64_t)num_threads_to_launch)) {
        end_j++;
      }
      args[thread_index] = (struct batch_row_thread_arg){
          .i = i + a,
          .start_j = start_j,
          .end_j = end_j,
          .vector = &(vector[a * g->num_nodes]),
          .g = g,
      };
      if (pthread_create(&(threads[thread_index]), NULL, batch_row_thread, &(args[thread_index])) != 0) {
        perror("Failed to create batch row thread\n");
        return 1;
      }
      start_j = end_j;
      thread_index++;
    }
  }

  for (int16_t a = 0; a < num_threads; a++) {
    if (pthread_join(threads[a], NULL) != 0) {
      perror("failed to join matrix thread\n");
      return 1;
    }
  }

  return 0;
}

void *matrix_thread(void *args) {
  int32_t thread_rank = (int32_t)(((struct matrix_thread_arg *)args)->thread_rank);
  int32_t thread_total_count = (int32_t)(((struct matrix_thread_arg *)args)->thread_total_count);
  struct matrix *m = ((struct matrix_thread_arg *)args)->m;
  struct graph *g = ((struct matrix_thread_arg *)args)->g;
  for (uint64_t i = (uint64_t)thread_rank; i < m->a; i += thread_total_count) {
    switch (m->fp_mode) {
    case FP32:
      m->bfr.fp32[i * m->b + i] = 0.0f;
      break;
    case FP64:
      m->bfr.fp64[i * m->b + i] = 0.0;
      break;
    }
    for (uint64_t j = i + 1; j < m->b; j++) {
      switch (m->fp_mode) {
      case FP32:
#if ENABLE_AVX256 == 1
        m->bfr.fp32[i * m->b + j] =
            (float)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
        m->bfr.fp32[i * m->b + j] =
            (float)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
        m->bfr.fp32[j * m->b + i] = m->bfr.fp32[i * m->b + j];
        break;
      case FP64:
        m->bfr.fp64[i * m->b + j] =
            (double)cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
        m->bfr.fp64[j * m->b + i] = m->bfr.fp64[i * m->b + j];
        break;
      }
    }
  }
  return NULL;
}

int32_t distance_matrix_from_graph_multithread(struct graph *const g, struct matrix *const m,
                                               const int16_t num_matrix_threads) {
  if (m->a != m->b) {
    perror("m->a != m->b\n");
    return 1;
  }
  if (((uint32_t)g->num_nodes) != m->a) {
    perror("g->num_nodes != m->a\n");
    printf("g->num_nodes: %lu\n", g->num_nodes);
    printf("m->a: %u\n", m->a);
    return 1;
  }

  pthread_t threads[num_matrix_threads];
  struct matrix_thread_arg args[num_matrix_threads];
  for (int32_t i = 0; i < num_matrix_threads; i++) {
    args[i].thread_rank = i;
    args[i].thread_total_count = num_matrix_threads;
    args[i].m = m;
    args[i].g = g;
    if (pthread_create(&(threads[i]), NULL, matrix_thread, &(args[i])) != 0) {
      perror("failed to create matrix thread\n");
      return 1;
    }
  }

  for (int32_t i = 0; i < num_matrix_threads; i++) {
    if (pthread_join(threads[i], NULL) != 0) {
      perror("failed to join matrix thread\n");
      return 1;
    }
  }

  return 0;
}

void free_matrix(struct matrix *m) {
  // if(m->to_free == 0){return;}
  /*
  switch(m->fp_mode){
          case FP32:
                  free(m->bfr.fp32);
                  break;
          case FP64:
                  free(m->bfr.fp64);
                  break;
          default:
                  perror("Unknown m->fp_mode in free_matrix\n");
                  exit(1);
  }
  */
  free((void *)m->bfr.fp32);
  free(m->active);
  free(m->active_final);
  // m->to_free = 0;
}

int32_t stats_matrix(struct matrix *m, double *avg_p, double *std_p, double *min_p, double *max_p) {
  double sum = 0.0;
  double min;
  double max;
  if (m->a <= 0 || m->b <= 0) {
    perror("incorrect matrix dimensions\n");
    return 1;
  }
  switch (m->fp_mode) {
  case FP32:
    min = (double)m->bfr.fp32[0];
    max = (double)m->bfr.fp32[0];
    break;
  case FP64:
    min = m->bfr.fp64[0];
    max = m->bfr.fp64[0];
    break;
  default:
    perror("unknown floating point mode\n");
    return 1;
  }
  for (uint64_t i = 0; i < m->a; i++) {
    for (uint64_t j = 0; j < m->b; j++) {
      double value;
      switch (m->fp_mode) {
      case FP32:
        value = (double)m->bfr.fp32[i * m->b + j];
        break;
      case FP64:
        value = m->bfr.fp64[i * m->b + j];
        break;
      default:
        perror("unknown floating point mode\n");
        return 1;
      }
      if (value < min) {
        min = value;
      }
      if (value > max) {
        max = value;
      }
      sum += value;
    }
  }
  double avg = sum / ((double)(m->a * m->b));

  sum = 0.0;
  for (uint64_t i = 0; i < m->a; i++) {
    for (uint64_t j = 0; j < m->b; j++) {
      double value;
      switch (m->fp_mode) {
      case FP32:
        value = (double)m->bfr.fp32[i * m->b + j];
        break;
      case FP64:
        value = m->bfr.fp64[i * m->b + j];
        break;
      default:
        perror("unknown floating point mode\n");
        return 1;
      }
      sum += pow(value - avg, 2.0);
    }
  }

  double std = pow(sum / ((double)(m->a * m->b)), 0.5);

  (*avg_p) = avg;
  (*std_p) = std;
  (*min_p) = min;
  (*max_p) = max;

  return 0;
}

int32_t void_strcmp(const void *a, const void *b) {
  return strcmp((char *)a, (char *)b);
}

// int32_t create_graph_node(struct graph_node* restrict node, const int16_t num_dimensions, const int8_t fp_mode){
int32_t create_graph_node(struct graph_node *restrict node, const uint16_t num_dimensions, const uint8_t fp_mode) {
  static int32_t num_graph_node_created = 0;

  size_t malloc_size;

  // node->mutex_local_node = (PTHREAD_MUTEX_INITIALIZER);
  pthread_mutex_init(&node->mutex_local_node, NULL);

  node->already_considered = 0;

  node->num_dimensions = num_dimensions;

  node->num_neighbours = 0;
  node->capacity_neighbours = GRAPH_NODE_NEIGHBOUR_STEP;
  malloc_size = node->capacity_neighbours * sizeof(struct graph_neighbour);
  node->neighbours = (struct graph_neighbour *)malloc(malloc_size);
  if (node->neighbours == NULL) {
    fprintf(stderr, "failed to initialise node->neighbours with %li bytes\n", malloc_size);
    return 1;
  }
  memset(node->neighbours, '\0', malloc_size);

  if (INITIALIZE_GRAPH_NODE_WITH_RANDOM_VECTOR) {
    malloc_size = 0;
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
      malloc_size = num_dimensions * sizeof(float);
      break;
    case GRAPH_NODE_FP64:
      malloc_size = num_dimensions * sizeof(double);
      break;
    }
    void *malloc_pointer = malloc(malloc_size);
    if (malloc_pointer == NULL) {
      goto malloc_fail;
    }
    memset(malloc_pointer, 0, malloc_size);
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
      (*node).vector.fp32 = (float *)malloc_pointer;
      break;
    case GRAPH_NODE_FP64:
      (*node).vector.fp64 = (double *)malloc_pointer;
      break;
    }

    for (int16_t i = 0; i < num_dimensions; i++) {
      double current_rand = (double)(rand() % RANDOM_MODULO);
      current_rand /= RANDOM_MODULO;
      current_rand *= (RANDOM_MAX_VALUE - RANDOM_MIN_VALUE);
      current_rand += RANDOM_MIN_VALUE;
      if (ENABLE_EXTREMES && num_graph_node_created % EXTREME_STEP == 0) {
        current_rand *= EXTREME_RATIO;
      }
      switch (fp_mode) {
      case GRAPH_NODE_FP32:
        (*node).vector.fp32[i] = (float)current_rand;
        break;
      case GRAPH_NODE_FP64:
        (*node).vector.fp64[i] = current_rand;
        break;
      }
    }
  } else {
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
      node->vector.fp32 = NULL;
      break;
    case GRAPH_NODE_FP64:
      node->vector.fp64 = NULL;
      break;
    }
  }

  node->relative_proportion = 0.0;
  node->absolute_proportion = 0;

  num_graph_node_created++;

  return 0;

malloc_fail:
  perror("malloc failed\n");
  return 1;
}

int32_t request_more_neighbour_capacity_graph_node(struct graph_node *node) {
  int32_t new_capacity = node->capacity_neighbours + GRAPH_NODE_NEIGHBOUR_STEP;
  size_t malloc_size = new_capacity * sizeof(struct graph_neighbour);
  node->neighbours = (struct graph_neighbour *)realloc(node->neighbours, malloc_size);
  if (node->neighbours == NULL) {
    fprintf(stderr, "realloc (graph node neighbours) failed for %li bytes\n", malloc_size);
    return 1;
  }
  node->capacity_neighbours = new_capacity;

  return 0;
}

void free_graph_node(struct graph_node *restrict node, const int8_t fp_mode) {
  switch (fp_mode) {
  case GRAPH_NODE_FP32:
    free((*node).vector.fp32);
    break;
  case GRAPH_NODE_FP64:
    free((*node).vector.fp64);
    break;
  }
  free(node->neighbours);
  pthread_mutex_destroy(&(node->mutex_local_node));
}

int32_t create_graph(struct graph *restrict const g, const int32_t num_nodes, const int16_t num_dimensions,
                     const int8_t fp_mode) {
  size_t malloc_size = num_nodes * sizeof(struct graph_node);
  void *malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, 0, malloc_size);
  g->nodes = (struct graph_node *)malloc_pointer;
  g->num_nodes = num_nodes;
  g->capacity = num_nodes;
  g->num_dimensions = num_dimensions;
  // g->dist_mat.to_free = 0;
  g->dist_mat = (struct matrix){
      .fp_mode = fp_mode,
  };
  g->dist_mat_must_be_freed = 0;
  pthread_mutex_init(&(g->mutex_nodes), NULL);
  pthread_mutex_init(&(g->mutex_matrix), NULL);

  double relative_proportion_sum = 0.0;

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    int32_t res;
    double proportion;
    res = create_graph_node(&(g->nodes[i]), num_dimensions, fp_mode);
    if (res != 0) {
      goto create_graph_node_failed;
    }
    if (CONSTANT_RELATIVE_PROPORTION) {
      proportion = 1.0;
    } else {
      proportion = (double)(rand() % MAX_ABUNDANCE);
    }
    g->nodes[i].relative_proportion = proportion;
    g->nodes[i].absolute_proportion = (int32_t)proportion;

    relative_proportion_sum += proportion;
  }
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    g->nodes[i].relative_proportion /= relative_proportion_sum;
  }

  return 0;

malloc_fail:
  perror("malloc failed\n");
  printf("errno: %i\n", errno);
  goto return_failure;

create_graph_node_failed:
  perror("create_graph_node failed\n");
  free(g->nodes);
  goto return_failure;

return_failure:
  return 1;
}

int32_t request_more_capacity_graph(struct graph *restrict const g) {
  void *malloc_pointer = realloc(g->nodes, (g->capacity + GRAPH_CAPACITY_STEP) * sizeof(struct graph_node));
  if (malloc_pointer == NULL) {
    perror("failed to realloc\n");
    return 1;
  }
  g->nodes = (struct graph_node *)malloc_pointer;
  memset(&(g->nodes[g->capacity]), '\0', GRAPH_CAPACITY_STEP * sizeof(struct graph_node));
  g->capacity += GRAPH_CAPACITY_STEP;
  return 0;
}

int32_t create_graph_empty(struct graph *restrict const g) {
  g->nodes = NULL;
  g->num_nodes = 0;
  g->capacity = 0;
  g->dist_mat = (struct matrix){
      .fp_mode = FP32,
  };
  g->dist_mat_must_be_freed = 0;

  if (request_more_capacity_graph(g) != 0) {
    return 1;
  }

  return 0;
}

void compute_graph_relative_proportions(struct graph *const g) {
  uint64_t sum = 0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum += g->nodes[i].absolute_proportion;
  }
  /*
  if(sum == 0){
          perror("sum == 0; division by zero\n");
  }
  */
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    g->nodes[i].relative_proportion = ((double)g->nodes[i].absolute_proportion) / ((double)sum);
  }
}

int32_t compute_graph_dist_mat(struct graph *const g, const int16_t num_matrix_threads) {
  /*
  size_t alloc_size;

  alloc_size = (size_t) (g->num_nodes * g->num_nodes * sizeof(float));
  g->dist_mat = (float*) realloc(g->dist_mat, alloc_size);
  if(g->dist_mat == NULL){
          perror("failed to realloc\n");
          return 1;
  }
  memset(g->dist_mat, '\0', alloc_size);

  for(i = 0 ; i < g->num_nodes ; i++){
          for(j = i + 1 ; j
  }
  */
  if (distance_matrix_from_graph_multithread(g, &(g->dist_mat), num_matrix_threads) != 0) {
    perror("failed to call distance_matrix_from_graph\n");
    return 1;
  }

  g->dist_mat_must_be_freed = 1;

  return 0;
}

int32_t word2vec_entry_cmp(const void *restrict a, const void *restrict b) {
  return strcmp(((struct word2vec_entry *)a)->key, ((struct word2vec_entry *)b)->key);
}

int32_t load_word2vec_binary(struct word2vec *restrict w2v, const char *restrict path) {
  const int32_t log_bfr_size = 256;
  char log_bfr[log_bfr_size];
  // printf("Creating graph from word2vec binary %s\n", path);
  memset(log_bfr, '\0', log_bfr_size * sizeof(char));
  snprintf(log_bfr, log_bfr_size, "Creating graph from word2vec binary: %s", path);
  info_format(__FILE__, __func__, __LINE__, log_bfr);

  FILE *file_p;
  file_p = fopen(path, "r");
  if (file_p == NULL) {
    perror("failed to open word2vec file\n");
    return 1;
  }
  const int32_t first_line_buffer_size = 64;
  char bfr_read[first_line_buffer_size];
  memset(bfr_read, '\0', first_line_buffer_size);
  if (fgets(bfr_read, first_line_buffer_size, file_p) == 0) { // compliance
    perror("initial fgets didn't read anything\n");
    return 1;
  }
  char *start_num_dimensions = strchr(bfr_read, ' ');
  if (start_num_dimensions == NULL) {
    perror("cannot find number of dimensions in word2vec binary\n");
    fclose(file_p);
    return 1;
  }
  int64_t num_vectors = strtol(bfr_read, &start_num_dimensions, 10);
  start_num_dimensions++;
  char *end_num_dimensions = strchr(start_num_dimensions, '\n');
  if (end_num_dimensions == NULL) {
    perror("cannot find end of number of dimensions in word2vec binary\n");
    fclose(file_p);
    return 1;
  }
  int16_t num_dimensions = (int16_t)strtol(start_num_dimensions, &end_num_dimensions, 10);

  // printf("number of nodes: %li\nnumber of dimensions: %i\n", num_vectors, num_dimensions);
  memset(log_bfr, '\0', log_bfr_size * sizeof(char));
  snprintf(log_bfr, log_bfr_size, "Number of nodes: %li; number of dimensions: %i", num_vectors, num_dimensions);
  info_format(__FILE__, __func__, __LINE__, log_bfr);

  w2v->num_dimensions = num_dimensions;
  w2v->num_vectors = num_vectors;

  const size_t word2vec_vector_buffer_read_size = 1024;
  char vector_buffer_read[word2vec_vector_buffer_read_size];
  memset(vector_buffer_read, '\0', word2vec_vector_buffer_read_size);

  size_t malloc_size = num_vectors * num_dimensions * sizeof(float);
  void *malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, '\0', malloc_size);
  w2v->vectors = (float *)malloc_pointer;

  malloc_size = w2v->num_vectors * sizeof(struct word2vec_entry);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, '\0', malloc_size);
  w2v->keys = (struct word2vec_entry *)malloc_pointer;

  for (uint64_t i = 0; i < w2v->num_vectors; i++) {
    pthread_mutex_init(&(w2v->keys[i].mutex), NULL);
  }

  uint64_t parsing_key = 1;
  uint64_t h = 0;
  uint64_t i = 0;
  uint64_t j = 0;
  uint64_t k = 0;
  uint64_t actual_vector_buffer_read_size;
  actual_vector_buffer_read_size = (uint64_t)fread(vector_buffer_read, 1, word2vec_vector_buffer_read_size, file_p);
  j = 0;

  while (i < w2v->num_vectors) {
    if (parsing_key) {
      while (j < actual_vector_buffer_read_size) { // fread
        if (vector_buffer_read[j] == ' ') {
          parsing_key = 0;
          j++;
          break;
        }
        if (h < WORD2VEC_KEY_BUFFER_SIZE - 1) {
          w2v->keys[i].key[h] = vector_buffer_read[j];
          h++;
        }
        j++;
      }
    }

    if (j == actual_vector_buffer_read_size) { // fread
      actual_vector_buffer_read_size = fread(vector_buffer_read, 1, word2vec_vector_buffer_read_size, file_p);
      j = 0;
    }

    if (!parsing_key) {
      while (k < w2v->num_dimensions && j < actual_vector_buffer_read_size) { // fread
        if (j + sizeof(float) >= (uint64_t)actual_vector_buffer_read_size) {
          fseek(file_p, j - actual_vector_buffer_read_size, SEEK_CUR);
          actual_vector_buffer_read_size = fread(vector_buffer_read, 1, word2vec_vector_buffer_read_size, file_p);
          j = 0;
        }
        memcpy((void *)&(w2v->vectors[i * w2v->num_dimensions + k]), (void *)&(vector_buffer_read[j]), sizeof(float));
        k++;
        j += sizeof(float);
      }
      if (k == w2v->num_dimensions) {
        k = 0;
        h = 0;
        parsing_key = 1;
        w2v->keys[i].vector = &(w2v->vectors[i * w2v->num_dimensions]);
        i++;
        while (j < actual_vector_buffer_read_size && vector_buffer_read[j] == '\n') { // fread
          j++;
        }
      }
    }

    if (j >= actual_vector_buffer_read_size) {
      actual_vector_buffer_read_size = fread(vector_buffer_read, 1, word2vec_vector_buffer_read_size, file_p);
      j = 0;
    }
  }

  fclose(file_p);

  for (uint64_t i = 0; i < w2v->num_vectors; i++) {
    w2v->keys[i].active_in_current_graph = 0;
  }

  qsort((void *)w2v->keys, w2v->num_vectors, sizeof(struct word2vec_entry), word2vec_entry_cmp);

  return 0;

malloc_fail:
  perror("malloc failed\n");
  printf("errno: %i\n", errno);
  if (w2v->vectors != NULL) {
    free(w2v->vectors);
  }
  if (w2v->keys != NULL) {
    free(w2v->keys);
  }
  fclose(file_p);
  goto return_failure;

return_failure:
  return 1;
}

void free_word2vec(struct word2vec *restrict w2v) {
  for (uint64_t i = 0; i < w2v->num_vectors; i++) {
    pthread_mutex_destroy(&(w2v->keys[i].mutex));
  }
  free(w2v->vectors);
  free(w2v->keys);
}

int32_t word2vec_key_to_index(const struct word2vec *restrict w2v, const char *restrict key) {
  int32_t lower_bound = 0;
  int32_t higher_bound = w2v->num_vectors;
  while (lower_bound <= higher_bound) {
    int32_t middle_index = lower_bound + floor((higher_bound - lower_bound) / 2.0);
    int32_t cmp_res = strcmp(w2v->keys[middle_index].key, key);
    if (cmp_res == 0) {
      return middle_index;
    } else if (cmp_res < 0) {
      lower_bound = middle_index + 1;
    } else {
      higher_bound = middle_index - 1;
    }
  }
  return -1;
}

struct word2vec_entry *word2vec_find_closest(const struct word2vec *restrict w2v, const char *restrict target) {
  const int64_t index_target = word2vec_key_to_index(w2v, target);
  if (index_target == -1) {
    printf("cannot find target word\n");
    return NULL;
  }
  uint64_t index_best = 0;
  int32_t found_starting_point = 0;
  while (index_best < w2v->num_vectors) {
    if (strcmp(w2v->keys[index_best].key, target) != 0) {
      found_starting_point = 1;
      break;
    }
    if (index_best == UINT64_MAX) {
      perror("UINT64 overflow\n");
      exit(1);
    }
    index_best++;
  }
  if (!found_starting_point) {
    printf("cannot find starting point\n");
    return NULL;
  }
  float distance_best;

#if ENABLE_AVX256 == 1
  distance_best = cosine_distance_fp32_avx(w2v->keys[index_best].vector, w2v->keys[index_target].vector, w2v->num_dimensions);
#else
  distance_best = cosine_distance_fp32(w2v->keys[index_best].vector, w2v->keys[index_target].vector, w2v->num_dimensions);
#endif
  for (uint64_t i = 1; i < w2v->num_vectors; i++) {
    if (i == (uint64_t)index_target) {
      continue;
    }
    float local_distance;
#if ENABLE_AVX256 == 1
    local_distance = cosine_distance_fp32_avx(w2v->keys[i].vector, w2v->keys[index_target].vector, w2v->num_dimensions);
#else
    local_distance = cosine_distance_fp32(w2v->keys[i].vector, w2v->keys[index_target].vector, w2v->num_dimensions);
#endif
    if (local_distance < distance_best) {
      distance_best = local_distance;
      index_best = i;
    }
  }
  return &(w2v->keys[index_best]);
}

void free_graph(struct graph *restrict g) {
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    // free_graph_node(&(g->nodes[i]), FP32); // not needed as pointers to w2v entries
    free(g->nodes[i].neighbours);
  }
  free(g->nodes);
  pthread_mutex_destroy(&(g->mutex_nodes));
  // if(g->dist_mat != NULL){free(g->dist_mat);}
  if (g->dist_mat_must_be_freed) {
    free_matrix(&(g->dist_mat));
    g->dist_mat_must_be_freed = 0;
  }
  pthread_mutex_destroy(&(g->mutex_matrix));
}

void create_distance_two_nodes(struct distance_two_nodes *restrict distance, struct graph_node *restrict a,
                               struct graph_node *restrict b, const int8_t fp_mode, const float *const distance_value) {
  if (a < b) {
    (*distance).a = a;
    (*distance).b = b;
  } else {
    (*distance).a = b;
    (*distance).b = a;
  }

  if (distance_value == NULL) {
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
#if ENABLE_AVX256 == 1
      distance->distance =
          cosine_distance_fp32_avx(distance->a->vector.fp32, distance->b->vector.fp32, (int32_t)a->num_dimensions);
#else
      distance->distance = cosine_distance_fp32(distance->a->vector.fp32, distance->b->vector.fp32, (int32_t)a->num_dimensions);
#endif
      break;
    case GRAPH_NODE_FP64:
      distance->distance = cosine_distance(distance->a->vector.fp64, distance->b->vector.fp64, (int32_t)a->num_dimensions);
      break;
    }
  } else {
    distance->distance = (*distance_value);
  }
  (*distance).usable = 1;
}

void siftdown_min_heap(struct graph_distance_heap *restrict heap, const uint64_t current_index) {
  uint64_t child_a_index;
  uint64_t child_b_index;
  if ((UINT64_MAX - 2) / 2 >= current_index) {
    child_a_index = (current_index * 2 + 1);
    child_b_index = child_a_index + 1;
  } else {
    perror("UINT64 overflow\n");
    exit(1);
  }
  const int32_t reachable_a = child_a_index < heap->num_distances;
  const int32_t reachable_b = child_b_index < heap->num_distances;

  if (reachable_a && heap->distances[child_a_index].distance < heap->distances[current_index].distance) { // better_a
    if (reachable_b && heap->distances[child_b_index].distance < heap->distances[current_index].distance &&
        !(heap->distances[child_a_index].distance < heap->distances[child_b_index].distance)) { // better_b && !a_better_than_b
      struct distance_two_nodes swap = heap->distances[current_index];
      heap->distances[current_index] = heap->distances[child_b_index];
      heap->distances[child_b_index] = swap;
      siftdown_min_heap(heap, child_b_index);
    } else {
      struct distance_two_nodes swap = heap->distances[current_index];
      heap->distances[current_index] = heap->distances[child_a_index];
      heap->distances[child_a_index] = swap;
      siftdown_min_heap(heap, child_a_index);
    }
  } else if (reachable_b && heap->distances[child_b_index].distance < heap->distances[current_index].distance) { // better_b
    struct distance_two_nodes swap = heap->distances[current_index];
    heap->distances[current_index] = heap->distances[child_b_index];
    heap->distances[child_b_index] = swap;
    siftdown_min_heap(heap, child_b_index);
  } else if (reachable_a && reachable_b &&
             heap->distances[child_b_index].distance <
                 heap->distances[child_a_index].distance) { // modification to potentially improve the speed of pop_graph by
                                                            // always having the lowest of the two childrens on the left
    struct distance_two_nodes swap = heap->distances[child_a_index];
    heap->distances[child_a_index] = heap->distances[child_b_index];
    heap->distances[child_b_index] = swap;
    siftdown_min_heap(heap, child_a_index);
    siftdown_min_heap(heap, child_b_index);
  }
}

int32_t siftdown(struct graph_distance_heap *restrict heap, const int32_t direction, const int64_t current_index) {
  if (direction == MIN_HEAP) {
    siftdown_min_heap(heap, current_index);
  } else if (direction == MAX_HEAP) {
    /* ... */
  } else {
    return 1;
  }

  return 0;
}

void heapify_min_heap(struct graph_distance_heap *restrict heap, const uint64_t current_index) {
  uint64_t child_a_index;
  uint64_t child_b_index;
  if ((UINT64_MAX - 2) / 2 >= current_index) {
    child_a_index = (current_index * 2 + 1);
    child_b_index = child_a_index + 1;
  } else {
    perror("UINT64 overflow\n");
    exit(1);
  }

  if (child_a_index < heap->num_distances) {
    heapify_min_heap(heap, child_a_index);
  }
  if (child_b_index < heap->num_distances) {
    heapify_min_heap(heap, child_b_index);
  }

  siftdown_min_heap(heap, current_index);
}

int32_t heapify(struct graph_distance_heap *restrict heap, const int32_t direction, const uint64_t current_index) {
  uint64_t child_a_index;
  uint64_t child_b_index;
  if ((UINT64_MAX - 2) / 2 >= current_index) {
    child_a_index = (current_index * 2 + 1);
    child_b_index = child_a_index + 1;
  } else {
    perror("UINT64 overflow\n");
    return 1;
  }

  if (direction == MIN_HEAP) {
    if (child_a_index < heap->num_distances) {
      heapify_min_heap(heap, child_a_index);
    }
    if (child_b_index < heap->num_distances) {
      heapify_min_heap(heap, child_b_index);
    }
    siftdown_min_heap(heap, current_index);
  } else if (direction == MAX_HEAP) {
    /* ... */
  } else {
    return 1;
  }
  return 0;
}

int32_t create_graph_distance_heap(struct graph_distance_heap *restrict heap, struct graph *restrict g, const int8_t fp_mode,
                                   const struct matrix *const m_) {
  if (m_ == NULL) {
    perror("for efficiency reasons, pass an actual matrix to create_graph_distance_heap instead of NULL\n");
  }

  (*heap).g = g;

  int64_t num_distances = 0;
  for (uint64_t i = (*g).num_nodes - 1; i >= 1; i--) {
    num_distances += (int64_t)i;
  }

  size_t malloc_size;
  void *malloc_pointer;
  malloc_size = num_distances * sizeof(struct distance_two_nodes);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, 0, malloc_size);
  (*heap).distances = (struct distance_two_nodes *)malloc_pointer;
  (*heap).num_distances = num_distances;

  uint64_t distance_index = 0;
  int64_t m = 0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = i + 1; j < g->num_nodes; j++) {
      if (m_ == NULL) {
        create_distance_two_nodes(&(heap->distances[distance_index]), &(g->nodes[i]), &(g->nodes[j]), fp_mode, NULL);
      } else {
        create_distance_two_nodes(&(heap->distances[distance_index]), &(g->nodes[i]), &(g->nodes[j]), fp_mode,
                                  &(m_->bfr.fp32[i * m_->b + j]));
      }
      if (distance_index == UINT64_MAX) {
        perror("UINT64 overflow\n");
        goto return_failure;
      }
      distance_index++;
      m++;
    }
  }
  if (distance_index != (*heap).num_distances) {
    perror("distance_index != (*heap).num_distances\n");
    return 1;
  }

  heapify(heap, MIN_HEAP, 0);

  return 0;

malloc_fail:
  perror("malloc failed\n");
  goto return_failure;

return_failure:
  return 1;
}

void pop_graph_distance_min_heap(struct graph_distance_heap *restrict heap, const uint64_t current_index) {
  uint64_t child_a_index;
  uint64_t child_b_index;
  if ((UINT64_MAX - 2) / 2 >= current_index) {
    child_a_index = (current_index * 2 + 1);
    child_b_index = child_a_index + 1;
  } else {
    perror("UINT64 overflow\n");
    exit(1);
  }
  int32_t reachable_a = child_a_index < (*heap).num_distances;
  int32_t reachable_b = child_b_index < (*heap).num_distances;
  struct distance_two_nodes swap;

  heap->distances[current_index].usable = 0;

  if (reachable_a) {
    reachable_a &= (*heap).distances[child_a_index].usable;
  }
  if (reachable_b) {
    reachable_b &= (*heap).distances[child_b_index].usable;
  }

  if (reachable_a) {
    if (reachable_b) {
      if ((*heap).distances[child_a_index].distance <= (*heap).distances[child_b_index].distance) { // a_better_than_b
        swap = (*heap).distances[current_index];
        (*heap).distances[current_index] = (*heap).distances[child_a_index];
        (*heap).distances[child_a_index] = swap;

        pop_graph_distance_min_heap(heap, child_a_index);

      } else {
        swap = (*heap).distances[current_index];
        (*heap).distances[current_index] = (*heap).distances[child_b_index];
        (*heap).distances[child_b_index] = swap;

        pop_graph_distance_min_heap(heap, child_b_index);
      }
    } else {
      swap = (*heap).distances[current_index];
      (*heap).distances[current_index] = (*heap).distances[child_a_index];
      (*heap).distances[child_a_index] = swap;

      pop_graph_distance_min_heap(heap, child_a_index);
    }
  }
}

int32_t pop_graph_distance_heap(struct graph_distance_heap *heap, const int32_t direction, const int32_t current_index) {
  if (direction == MIN_HEAP) {
    pop_graph_distance_min_heap(heap, current_index);
  } else if (direction == MAX_HEAP) {
    /* ... */
  } else {
    return 1;
  }

  return 0;
}

void free_graph_distance_heap(struct graph_distance_heap *heap) {
  free((*heap).distances);
}

int32_t create_minimum_spanning_tree(struct minimum_spanning_tree *mst, struct graph_distance_heap *heap) {
  size_t malloc_size;
  void *malloc_pointer;

  int32_t num_nodes;
  int32_t num_distances;

  (*mst).heap = heap;

  num_nodes = (*((*heap).g)).num_nodes;
  num_distances = num_nodes - 1;

  malloc_size = num_distances * sizeof(struct distance_two_nodes);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, 0, malloc_size);
  (*mst).distances = (struct distance_two_nodes *)malloc_pointer;
  (*mst).num_distances = num_distances;
  (*mst).num_active_distances = 0;

  malloc_size = num_nodes * sizeof(struct graph_node *);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, 0, malloc_size);
  (*mst).nodes = (struct graph_node **)malloc_pointer;
  (*mst).num_nodes = num_nodes;
  (*mst).num_active_nodes = 0;

  return 0;

malloc_fail:
  perror("malloc failed\n");
  return 1;
}

void free_minimum_spanning_tree(struct minimum_spanning_tree *mst) {
  free((*mst).distances);
  free((*mst).nodes);
}

int32_t find_minimum_acceptable_arc(struct minimum_spanning_tree *mst, uint64_t current_index, double distance_to_beat,
                                    int32_t consider_distance_to_beat) {
  if (consider_distance_to_beat && (*((*mst).heap)).distances[current_index].distance > distance_to_beat) {
    return -1;
  }

  int32_t current_node_usable = 1;

  int32_t a_already_present = mst->heap->distances[current_index].a->already_considered;
  int32_t b_already_present = mst->heap->distances[current_index].b->already_considered;

  current_node_usable = a_already_present ^ b_already_present;
  if (current_node_usable) {
    return current_index;
  }

  uint64_t child_a_index;
  uint64_t child_b_index;
  if ((UINT64_MAX - 2) / 2 >= current_index) {
    child_a_index = (current_index * 2 + 1);
    child_b_index = child_a_index + 1;
  } else {
    perror("UINT64 overflow\n");
    exit(1);
  }
  int32_t reachable_a = child_a_index < (*((*mst).heap)).num_distances;
  int32_t reachable_b = child_b_index < (*((*mst).heap)).num_distances;

  if (reachable_a) {
    reachable_a &= (*((*mst).heap)).distances[child_a_index].usable;
  }
  if (reachable_b) {
    reachable_b &= (*((*mst).heap)).distances[child_b_index].usable;
  }

  int32_t result_a = -1;
  int32_t result_b = -1;
  if (reachable_a) {
    result_a = find_minimum_acceptable_arc(mst, child_a_index, distance_to_beat, consider_distance_to_beat);
    if (result_a != -1) {
      double new_distance = (*((*mst).heap)).distances[result_a].distance;
      if (!consider_distance_to_beat || new_distance < distance_to_beat) {
        distance_to_beat = new_distance;
        consider_distance_to_beat = 1;
      }
    }
  }
  if (reachable_b) {
    result_b = find_minimum_acceptable_arc(mst, child_b_index, distance_to_beat, consider_distance_to_beat);
  }

  if (result_a == -1) {
    return result_b;
  } else {
    if (result_b == -1) {
      return result_a;
    } else {
      if ((*((*mst).heap)).distances[result_a].distance < (*((*mst).heap)).distances[result_b].distance) {
        return result_a;
      } else {
        return result_b;
      }
    }
  }
}

int32_t calculate_minimum_spanning_tree(struct minimum_spanning_tree *mst, struct matrix *m_in, int32_t method) {
  int64_t num_added_arc = 0;
  for (uint64_t i = 0; i < mst->heap->g->num_nodes; i++) {
    mst->heap->g->nodes[i].already_considered = 0;
  }

  if (method == MST_PRIMS_ALGORITHM) {
    if ((*mst).num_active_nodes == 0) {
      (*mst).distances[(*mst).num_active_distances] = (*((*mst).heap)).distances[0];
      (*mst).num_active_distances++;
      (*mst).nodes[(*mst).num_active_nodes] = (*((*mst).heap)).distances[0].a;
      (*mst).num_active_nodes++;
      (*mst).nodes[(*mst).num_active_nodes] = (*((*mst).heap)).distances[0].b;
      (*mst).num_active_nodes++;

      mst->heap->distances[0].a->already_considered = 1;
      mst->heap->distances[0].b->already_considered = 1;

      pop_graph_distance_heap((*mst).heap, MIN_HEAP, 0);
    }

    while ((*mst).num_active_nodes < (*mst).num_nodes) {
      int32_t a_already_present = 0;
      int32_t b_already_present = 0;

      int64_t index_of_arc_to_add = find_minimum_acceptable_arc(mst, 0, 0.0, 0);
      num_added_arc++;

      if (index_of_arc_to_add < 0 || ((uint64_t)index_of_arc_to_add) >= (*((*mst).heap)).num_distances) {
        break;
      }

      a_already_present = mst->heap->distances[index_of_arc_to_add].a->already_considered;
      b_already_present = mst->heap->distances[index_of_arc_to_add].b->already_considered;

      if (!(a_already_present ^ b_already_present)) {
        printf("a_already_present: %i, b_already_present: %i\n", a_already_present, b_already_present);
      }
      (*mst).distances[(*mst).num_active_distances] = (*((*mst).heap)).distances[index_of_arc_to_add];
      (*mst).num_active_distances++;
      if (!a_already_present) {
        (*mst).nodes[(*mst).num_active_nodes] = (*((*mst).heap)).distances[index_of_arc_to_add].a;
        (*mst).num_active_nodes++;

        mst->heap->distances[index_of_arc_to_add].a->already_considered = 1;
      }
      if (!b_already_present) {
        (*mst).nodes[(*mst).num_active_nodes] = (*((*mst).heap)).distances[index_of_arc_to_add].b;
        (*mst).num_active_nodes++;

        mst->heap->distances[index_of_arc_to_add].b->already_considered = 1;
      }
      pop_graph_distance_heap((*mst).heap, MIN_HEAP, index_of_arc_to_add);
    }
  }

  for (uint64_t i = 0; i < mst->num_distances; i++) {
    int64_t index_a = (int64_t)(mst->distances[i].a - mst->heap->g->nodes);
    int64_t index_b = (int64_t)(mst->distances[i].b - mst->heap->g->nodes);

    m_in->active[index_a * m_in->b + index_b] = 1;
    m_in->active[index_b * m_in->b + index_a] = 1;
    // assumes symmetry
    switch (m_in->fp_mode) {
    case FP32:
      m_in->bfr.fp32[index_a * m_in->b + index_b] = (float)mst->heap->distances[i].distance;
      m_in->bfr.fp32[index_b * m_in->b + index_a] = (float)mst->heap->distances[i].distance;
      break;
    case FP64:
      m_in->bfr.fp64[index_a * m_in->b + index_b] = mst->heap->distances[i].distance;
      m_in->bfr.fp64[index_b * m_in->b + index_a] = mst->heap->distances[i].distance;
      break;
    }
  }

  if ((*mst).num_active_nodes != (*mst).num_nodes || (*mst).num_active_distances != (*mst).num_distances) {
    printf("num_active_nodes: %lu, num_nodes: %lu, num_active_distances: %lu, num_distances: %lu\n", (*mst).num_active_nodes,
           (*mst).num_nodes, (*mst).num_active_distances, (*mst).num_distances);

    perror("improper number of nodes or distances when calculating minimum spanning tree\n");
  }
  return 0;
}

int32_t functional_evenness_from_minimum_spanning_tree(struct minimum_spanning_tree *mst, double *result_buffer) {
  // see Villéger et al. (2008)

  void *malloc_pointer;
  size_t malloc_size;
  double ew_sum;
  double *all_ew;
  double upper_sum;
  double one_over_n_minus_one;
  double result;

  malloc_size = (*mst).num_active_distances * sizeof(double);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, 0, malloc_size);
  all_ew = (double *)malloc_pointer;

  ew_sum = 0.0;
  one_over_n_minus_one = (1.0 / ((*mst).num_active_nodes - 1.0));

  for (uint64_t i = 0; i < (*mst).num_active_distances; i++) {
    double weight_a = (*(mst->distances[i].a)).relative_proportion;
    double weight_b = (*(mst->distances[i].b)).relative_proportion;
    all_ew[i] = mst->distances[i].distance / (weight_a + weight_b);
    ew_sum += all_ew[i];
  }

  upper_sum = 0.0;

  for (uint64_t i = 0; i < mst->num_active_distances; i++) {
    double pew = all_ew[i] / ew_sum;
    if (pew < one_over_n_minus_one) {
      upper_sum += pew;
    } else {
      upper_sum += one_over_n_minus_one;
    }
  }

  result = (upper_sum - one_over_n_minus_one) / (1.0 - one_over_n_minus_one);
  (*result_buffer) = result;

  return 0;

malloc_fail:
  perror("malloc failed\n");
  return 1;
}

int32_t functional_dispersion_from_graph(struct graph *const g, double *const result_buffer, const int8_t fp_mode) {
  // see Laliberté & Legendre (2010)

  struct graph_node centroid;
  int32_t res;
  double result;
  double sum_relative_proportion;

  res = create_graph_node(&centroid, (*g).nodes[0].num_dimensions, fp_mode);
  if (res != 0) {
    goto create_graph_node_failed;
  }

  centroid.num_dimensions = g->nodes[0].num_dimensions;
  size_t malloc_size;
  switch (fp_mode) {
  case FP32:
    if (centroid.vector.fp32 == NULL) {
      malloc_size = g->nodes[0].num_dimensions * sizeof(float);
      centroid.vector.fp32 = (float *)malloc(malloc_size);
      if (centroid.vector.fp32 == NULL) {
        perror("malloc failed\n");
        return 1;
      }
      memset(centroid.vector.fp32, '\0', malloc_size);
    }
    break;
  case FP64:
    if (centroid.vector.fp64 == NULL) {
      malloc_size = g->nodes[0].num_dimensions * sizeof(double);
      centroid.vector.fp64 = (double *)malloc(malloc_size);
      if (centroid.vector.fp64 == NULL) {
        perror("malloc failed\n");
        return 1;
      }
      memset(centroid.vector.fp64, '\0', malloc_size);
    }
    break;
  }

  for (uint64_t i = 0; i < (*g).num_nodes; i++) {
    for (uint64_t j = 0; j < (*g).nodes[i].num_dimensions; j++) {
      switch (fp_mode) {
      case GRAPH_NODE_FP32:
        centroid.vector.fp32[j] += (*g).nodes[i].vector.fp32[j] * (*g).nodes[i].relative_proportion;
        break;
      case GRAPH_NODE_FP64:
        centroid.vector.fp64[j] += (*g).nodes[i].vector.fp64[j] * (*g).nodes[i].relative_proportion;
        break;
      }
    }
  }

  result = 0.0;
  sum_relative_proportion = 0.0;
  for (uint64_t i = 0; i < (*g).num_nodes; i++) {
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
#if ENABLE_AVX256 == 1
      result += cosine_distance_fp32_avx(g->nodes[i].vector.fp32, centroid.vector.fp32, g->nodes[i].num_dimensions) *
                g->nodes[i].relative_proportion;
#else
      result += cosine_distance_fp32(g->nodes[i].vector.fp32, centroid.vector.fp32, g->nodes[i].num_dimensions) *
                g->nodes[i].relative_proportion;
#endif
      break;
    case GRAPH_NODE_FP64:
      result += cosine_distance(g->nodes[i].vector.fp64, centroid.vector.fp64, g->nodes[i].num_dimensions) *
                g->nodes[i].relative_proportion;
      break;
    }
    sum_relative_proportion += (*g).nodes[i].relative_proportion;
  }
  result /= sum_relative_proportion;
  (*result_buffer) = result;

  free(centroid.vector.fp32);

  return 0;

create_graph_node_failed:
  perror("created_graph_node_failed\n");
  return 1;
}

int32_t functional_divergence_modified_from_graph(struct graph *const g, double *const result_buffer, const int8_t fp_mode) {
  // see Villéger et al. (2008)
  // modifications: compute centroids based on all points instead of convex hull; cosine distance instead of euclidean

  struct graph_node centroid;
  int32_t res;

  res = create_graph_node(&centroid, g->nodes[0].num_dimensions, fp_mode);
  if (res != 0) {
    perror("created_graph_node_failed\n");
    return 1;
  }

  centroid.num_dimensions = g->nodes[0].num_dimensions;
  size_t malloc_size;
  switch (fp_mode) {
  case FP32:
    if (centroid.vector.fp32 == NULL) {
      malloc_size = g->nodes[0].num_dimensions * sizeof(float);
      centroid.vector.fp32 = (float *)malloc(malloc_size);
      if (centroid.vector.fp32 == NULL) {
        perror("malloc failed\n");
        return 1;
      }
      memset(centroid.vector.fp32, '\0', malloc_size);
    }
    break;
  case FP64:
    if (centroid.vector.fp64 == NULL) {
      malloc_size = g->nodes[0].num_dimensions * sizeof(double);
      centroid.vector.fp64 = (double *)malloc(malloc_size);
      if (centroid.vector.fp64 == NULL) {
        perror("malloc failed\n");
        return 1;
      }
      memset(centroid.vector.fp64, '\0', malloc_size);
    }
    break;
  }

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = 0; j < g->nodes[i].num_dimensions; j++) {
      switch (fp_mode) {
      case GRAPH_NODE_FP32:
        centroid.vector.fp32[j] += g->nodes[i].vector.fp32[j] * g->nodes[i].relative_proportion;
        break;
      case GRAPH_NODE_FP64:
        centroid.vector.fp64[j] += g->nodes[i].vector.fp64[j] * g->nodes[i].relative_proportion;
        break;
      }
    }
  }
  double distances_to_centroid[g->num_nodes];
  memset(distances_to_centroid, '\0', g->num_nodes * sizeof(double));
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    switch (fp_mode) {
    case GRAPH_NODE_FP32:
#if ENABLE_AVX256 == 1
      distances_to_centroid[i] =
          (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, centroid.vector.fp32, g->nodes[i].num_dimensions) *
          g->nodes[i].relative_proportion;
#else
      distances_to_centroid[i] =
          (double)cosine_distance_fp32(g->nodes[i].vector.fp32, centroid.vector.fp32, g->nodes[i].num_dimensions) *
          g->nodes[i].relative_proportion;
#endif
      break;
    case GRAPH_NODE_FP64:
      distances_to_centroid[i] = cosine_distance(g->nodes[i].vector.fp64, centroid.vector.fp64, g->nodes[i].num_dimensions) *
                                 g->nodes[i].relative_proportion;
      break;
    }
  }
  double sum_distances_to_centroid = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum_distances_to_centroid += distances_to_centroid[i];
  }
  double avg_distances_to_centroid = sum_distances_to_centroid / g->num_nodes;
  double weighted_deviance = 0.0;
  double weighted_deviance_abs = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    double local_deviance = distances_to_centroid[i] - avg_distances_to_centroid;
    double local_deviance_abs = local_deviance;
    if (local_deviance_abs < 0.0) {
      local_deviance_abs *= -1.0;
    }
    weighted_deviance += g->nodes[i].relative_proportion * local_deviance;
    weighted_deviance_abs += g->nodes[i].relative_proportion * local_deviance_abs;
  }

  double functional_divergence =
      (weighted_deviance + avg_distances_to_centroid) / (weighted_deviance_abs + avg_distances_to_centroid);
  (*result_buffer) = functional_divergence;

  free(centroid.vector.fp32);

  return 0;
}

/* ==== ITERATIVE FUNCTIONS ==== */

/* ---- PAIRWISE ---- */

int32_t create_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph *const restrict iter_state,
                                                   struct graph *const g) {
  iter_state->i = 0; // single-threaded version uses this one, but multi-threaded version uses other i
  iter_state->n = (g->num_nodes * (g->num_nodes - 1)) / 2;
  iter_state->result = 0.0;
  iter_state->g = g;
  pthread_mutex_init(&(iter_state->mutex), NULL);
  return 0;
}

void iterate_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph *const restrict iter_state,
                                                 const float *const vector) {
  for (uint64_t j = (uint64_t)(iter_state->i + 1); j < (uint64_t)iter_state->g->num_nodes; j++) {
    iter_state->result += (double)vector[j];
  }
  iter_state->i++;
}

void *iterate_iterative_state_pairwise_from_graph_thread(void *args) {
  struct iterative_state_pairwise_from_graph *iter_state = ((struct thread_args_aggregator *)args)->iter_state.pairwise;
  const float *const vector = ((struct thread_args_aggregator *)args)->vector;
  uint64_t j = ((struct thread_args_aggregator *)args)->i + 1;

  float sum = 0.0f;

  const uint64_t n = (uint64_t)iter_state->g->num_nodes;
  while (j < n) {
    sum += vector[j];
    j++;
  }

  pthread_mutex_lock(&(iter_state->mutex));
  iter_state->result += sum;
  pthread_mutex_unlock(&(iter_state->mutex));

  return NULL;
}

#if ENABLE_AVX256 == 1
void *iterate_iterative_state_pairwise_from_graph_avx_thread(void *args) {
  struct iterative_state_pairwise_from_graph *iter_state = ((struct thread_args_aggregator *)args)->iter_state.pairwise;
  const float *const vector = ((struct thread_args_aggregator *)args)->vector;
  uint64_t j = ((struct thread_args_aggregator *)args)->i + 1;

  __m256 avx_sum = _mm256_setzero_ps();

  const uint64_t n = (uint64_t)iter_state->g->num_nodes;
  while (j < n) {
    avx_sum = _mm256_add_ps(avx_sum, _mm256_loadu_ps(&(vector[j])));
    j += 8;
  }

  float sum = 0.0;
  float local_vec[8];
  _mm256_storeu_ps(local_vec, avx_sum);
  for (int32_t k = 0; k < 8; k++) {
    sum += local_vec[k];
  }

  while (j < n) {
    sum += vector[j];
    j++;
  }

  pthread_mutex_lock(&(iter_state->mutex));
  iter_state->result += sum;
  pthread_mutex_unlock(&(iter_state->mutex));

  return NULL;
}
#endif

void finalise_iterative_state_pairwise_from_graph(struct iterative_state_pairwise_from_graph *const restrict iter_state) {
  iter_state->result /= (double)iter_state->n;
  pthread_mutex_destroy(&(iter_state->mutex));
}

/* ---- STIRLING ---- */

int32_t create_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph *const restrict iter_state,
                                                   struct graph *const g, double alpha, double beta) {
  iter_state->i = 0;
  iter_state->n = (g->num_nodes * (g->num_nodes - 1)) / 2;
  iter_state->alpha = alpha;
  iter_state->beta = beta;
  iter_state->result = 0.0;
  iter_state->g = g;
  pthread_mutex_init(&(iter_state->mutex), NULL);
  return 0;
}

void iterate_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph *const restrict iter_state,
                                                 const float *const vector) {
  for (uint64_t j = 0; j < iter_state->g->num_nodes; j++) {
    iter_state->result +=
        pow((double)vector[j], iter_state->alpha) *
        pow(iter_state->g->nodes[iter_state->i].relative_proportion * iter_state->g->nodes[j].relative_proportion,
            iter_state->beta);
  }
  iter_state->i++;
}

void *iterate_iterative_state_stirling_from_graph_thread(void *args) {
  struct iterative_state_stirling_from_graph *iter_state = ((struct thread_args_aggregator *)args)->iter_state.stirling;
  const float *const vector = ((struct thread_args_aggregator *)args)->vector;
  uint64_t j = ((struct thread_args_aggregator *)args)->i + 1;

  double sum = 0.0;

  const uint64_t n = (uint64_t)iter_state->g->num_nodes;
  while (j < n) {
    sum += pow((double)vector[j], iter_state->alpha) *
           pow(iter_state->g->nodes[iter_state->i].relative_proportion * iter_state->g->nodes[j].relative_proportion,
               iter_state->beta);
    j++;
  }

  pthread_mutex_lock(&(iter_state->mutex));
  iter_state->result += sum;
  pthread_mutex_unlock(&(iter_state->mutex));

  return NULL;
}

void finalise_iterative_state_stirling_from_graph(struct iterative_state_stirling_from_graph *const restrict iter_state) {
  pthread_mutex_destroy(&(iter_state->mutex));
}

/* ---- LEINSTER_COBBOLD ---- */

int32_t create_iterative_state_leinster_cobbold_from_graph(
    struct iterative_state_leinster_cobbold_from_graph *const restrict iter_state, struct graph *const g, double alpha) {
  iter_state->i = 0;
  iter_state->n = (g->num_nodes * (g->num_nodes - 1)) / 2;
  iter_state->alpha = alpha;
  if (alpha != 1.0) {
    iter_state->hill_number = 0.0;
  } else {
    iter_state->hill_number = 1.0;
  }
  iter_state->entropy = 0.0;
  iter_state->g = g;
  pthread_mutex_init(&(iter_state->mutex), NULL);
  return 0;
}

void iterate_iterative_state_leinster_cobbold_from_graph(
    struct iterative_state_leinster_cobbold_from_graph *const restrict iter_state, const float *const vector) {
  const double u = 1.0;

  double local_agg = 0.0;
  for (uint64_t j = 0; j < iter_state->g->num_nodes; j++) {
    double distance = (double)vector[j];
    double similarity = 1.0 - distance;
    local_agg += iter_state->g->nodes[j].relative_proportion * pow(E, -u * similarity);
  }
  if (iter_state->alpha != 1.0) {
    iter_state->hill_number += pow(local_agg, iter_state->alpha - 1.0);
  } else {
    iter_state->hill_number *= pow(local_agg, iter_state->g->nodes[iter_state->i].relative_proportion);
  }
  iter_state->i++;
}

void *iterate_iterative_state_leinster_cobbold_from_graph_thread(void *args) {
  struct iterative_state_leinster_cobbold_from_graph *iter_state =
      ((struct thread_args_aggregator *)args)->iter_state.leinster_cobbold;
  const float *const vector = ((struct thread_args_aggregator *)args)->vector;
  // uint64_t j = ((struct thread_args_aggregator*) args)->i + 1;
  uint64_t j = 0;

  const double u = 1.0;

  double local_agg = 0.0;

  const uint64_t n = (uint64_t)iter_state->g->num_nodes;
  while (j < n) {
    double distance = (double)vector[j];
    double similarity = 1.0 - distance;
    local_agg += iter_state->g->nodes[j].relative_proportion * pow(E, -u * similarity);
    j++;
  }

  pthread_mutex_lock(&(iter_state->mutex));
  if (iter_state->alpha != 1.0) {
    iter_state->hill_number += pow(local_agg, iter_state->alpha - 1.0);
  } else {
    iter_state->hill_number *= pow(local_agg, iter_state->g->nodes[iter_state->i].relative_proportion);
  }
  pthread_mutex_unlock(&(iter_state->mutex));

  return NULL;
}

void finalise_iterative_state_leinster_cobbold_from_graph(
    struct iterative_state_leinster_cobbold_from_graph *const restrict iter_state) {
  const double LOGARITHMIC_BASE = E;

  if (iter_state->alpha != 1.0) {
    iter_state->hill_number = pow(iter_state->hill_number, 1.0 / (1.0 - iter_state->alpha));
  } else {
    iter_state->hill_number = pow(iter_state->hill_number, -1.0);
  }

  iter_state->entropy = log(iter_state->hill_number) / log(LOGARITHMIC_BASE);

  pthread_mutex_destroy(&(iter_state->mutex));
}

// ================

int32_t pairwise_from_graph(struct graph *const g, double *const result_buffer, const int8_t fp_mode,
                            const struct matrix *const m_) {
  // see Mouchet et al. (2010) referencing Walker, Kinzig & Langridge (1999)

  double result;

  result = 0.0;
  int64_t m = 0;
  int64_t n = (g->num_nodes * (g->num_nodes - 1)) / 2;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = i + 1; j < g->num_nodes; j++) {
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          result += m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          result += m_->bfr.fp64[i * m_->b + j];
          break;
        }
        m++;
        continue;
      }
      switch (fp_mode) {
      case GRAPH_NODE_FP32:
#if ENABLE_AVX256 == 1
        result += cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
        result += cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
        break;
      case GRAPH_NODE_FP64:
        result += cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
        break;
      }
      m++;
    }
  }

  result /= n;

  (*result_buffer) = result;

  return 0;
}

int32_t word2vec_to_graph_fp32(struct graph *g, struct word2vec *w2v, char **cupt_paths, char **cupt_paths_true_positives,
                               int32_t num_cupt_paths, int32_t ud_column) {
  for (uint64_t i = 0; i < w2v->num_vectors; i++) {
    w2v->keys[i].active_in_current_graph = 0;
    w2v->keys[i].num_occurrences = 0;
  }
  int32_t num_nodes = 0;
  for (int32_t i = 0; i < num_cupt_paths; i++) {
    if (cupt_paths_true_positives == NULL) {
      printf("processing %s\n", cupt_paths[i]);
    } else {
      printf("processing %s (true positives: %s)\n", cupt_paths[i], cupt_paths_true_positives[i]);
    }
    struct cupt_sentence_iterator csi;
    struct cupt_sentence_iterator csi_tp;
    int32_t err = create_cupt_sentence_iterator(&csi, cupt_paths[i]);
    if (err != 0) {
      perror("failed to call create_cupt_sentence_iterator\n");
      return 1;
    }
    if (cupt_paths_true_positives != NULL) {
      err = create_cupt_sentence_iterator(&csi_tp, cupt_paths_true_positives[i]);
      if (err != 0) {
        perror("failed to call create_cupt_sentence_iterator\n");
        free_cupt_sentence_iterator(&csi);
        return 1;
      }
    }
    err = iterate_cupt_sentence_iterator(&csi);
    if (err != 0) {
      perror("failed to call iterate_cupt_sentence_iterator\n");
      free_cupt_sentence_iterator(&csi);
      if (cupt_paths_true_positives != NULL) {
        free_cupt_sentence_iterator(&csi_tp);
      }
      return 1;
    }
    if (cupt_paths_true_positives != NULL) {
      err = iterate_cupt_sentence_iterator(&csi_tp);
      if (err != 0) {
        perror("failed to call iterate_cupt_sentence_iterator\n");
        free_cupt_sentence_iterator(&csi);
        free_cupt_sentence_iterator(&csi_tp);
        return 1;
      }
    }
    while (!csi.file_is_done) {
      const int32_t max_mwe = 32;
      const int32_t max_tokens_per_mwe = 32;
      const int32_t size_token_mwe = 32;
      size_t mwe_bfr_size = max_mwe * max_tokens_per_mwe * size_token_mwe;
      char mwe[mwe_bfr_size];
      memset(mwe, '\0', mwe_bfr_size);
      char mwe_tp[mwe_bfr_size];
      memset(mwe_tp, '\0', mwe_bfr_size);
      int32_t mwe_lengths[max_mwe];
      memset(mwe_lengths, '\0', max_mwe * sizeof(int32_t));
      int32_t mwe_tp_lengths[max_mwe];
      memset(mwe_tp_lengths, '\0', max_mwe * sizeof(int32_t));
      int8_t mwe_correct_span[max_mwe];
      memset(mwe_correct_span, 1, max_mwe * sizeof(int8_t));
      for (int32_t j = 0; j < csi.current_sentence.num_tokens; j++) {
        if (ud_column == UD_MWE) {
          if (strcmp(csi.current_sentence.tokens[j].mwe, "") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "_") == 0 ||
              strcmp(csi.current_sentence.tokens[j].mwe, "-") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "*") == 0) {
            continue;
          }

          char *strtok_placeholder = NULL;
          strtok_placeholder = strtok(csi.current_sentence.tokens[j].mwe, ";");
          while (strtok_placeholder != NULL) {
            int64_t mwe_num = strtol(strtok_placeholder, NULL, 10);
            if (mwe_num < max_mwe) {
              size_t bytes_to_cpy = strlen(csi.current_sentence.tokens[j].lemma);
              if (bytes_to_cpy > size_token_mwe - 1) {
                bytes_to_cpy = size_token_mwe - 1;
              }
              memcpy(&(mwe[mwe_num * max_tokens_per_mwe * size_token_mwe + mwe_lengths[mwe_num] * size_token_mwe]),
                     csi.current_sentence.tokens[j].lemma, bytes_to_cpy);
              mwe_lengths[mwe_num]++;
            }

            if (cupt_paths_true_positives != NULL && (strcmp(csi_tp.current_sentence.tokens[j].mwe, "") == 0 ||
                                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "_") == 0 ||
                                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "-") == 0 ||
                                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "*") == 0)) {
              mwe_correct_span[mwe_num] = 0;
            }

            strtok_placeholder = strtok(NULL, ";");
          }

          continue;
        }
        int32_t index;
        switch (ud_column) {
        case UD_FORM:
          index = word2vec_key_to_index(w2v, csi.current_sentence.tokens[j].form);
          break;
        case UD_LEMMA:
          index = word2vec_key_to_index(w2v, csi.current_sentence.tokens[j].lemma);
          break;
        default:
          perror("unknown UD column\n");
          return 1;
        }
        if (index != -1) {
          if (w2v->keys[index].active_in_current_graph == 0) {
            w2v->keys[index].active_in_current_graph = 1;
            num_nodes++;
          }
          w2v->keys[index].num_occurrences++;
        }
      }

      for (int32_t k = 0; k < max_mwe; k++) {
        if (mwe_lengths[k] == 0) {
          continue;
        }

        if (cupt_paths_true_positives != NULL && !(mwe_correct_span[k])) {
          continue;
        }

        qsort(&(mwe[k * max_tokens_per_mwe * size_token_mwe]), mwe_lengths[k], size_token_mwe, void_strcmp);

        const int32_t bfr_size = 512;
        int32_t index_bfr = 0;
        char bfr[bfr_size];
        memset(bfr, '\0', bfr_size);
        memcpy(bfr + index_bfr, "_MWE_", 5);
        index_bfr += 5;
        size_t bytes_to_cpy = 0;
        for (int32_t m = 0; m < mwe_lengths[k]; m++) {
          bytes_to_cpy = strlen(&(mwe[k * max_tokens_per_mwe * size_token_mwe + m * size_token_mwe]));
          if (m < mwe_lengths[k] - 1) {
            bytes_to_cpy++;
          }
          if (bytes_to_cpy > ((size_t)bfr_size - 1 - index_bfr)) {
            bytes_to_cpy = bfr_size - 1 - index_bfr;
          }
          if (m < mwe_lengths[k] - 1) {
            memcpy(bfr + index_bfr, &(mwe[k * max_tokens_per_mwe * size_token_mwe + m * size_token_mwe]), bytes_to_cpy - 1);
            index_bfr += bytes_to_cpy - 1;
            bfr[index_bfr] = '_';
            index_bfr++;
          } else {
            memcpy(bfr + index_bfr, &(mwe[k * max_tokens_per_mwe * size_token_mwe + m * size_token_mwe]), bytes_to_cpy);
            index_bfr += bytes_to_cpy;
          }
        }

        int32_t index = word2vec_key_to_index(w2v, bfr);
        if (index != -1) {
          if (w2v->keys[index].active_in_current_graph == 0) {
            w2v->keys[index].active_in_current_graph = 1;
            num_nodes++;
          }
          w2v->keys[index].num_occurrences++;
        }
      }

      err = iterate_cupt_sentence_iterator(&csi);
      if (err != 0) {
        perror("failed to call iterate_cupt_sentence_iterator\n");
        free_cupt_sentence_iterator(&csi);
        if (cupt_paths_true_positives != NULL) {
          free_cupt_sentence_iterator(&csi_tp);
        }
        return 1;
      }
      if (cupt_paths_true_positives != NULL) {
        err = iterate_cupt_sentence_iterator(&csi_tp);
        if (err != 0) {
          perror("failed to call iterate_cupt_sentence_iterator\n");
          free_cupt_sentence_iterator(&csi);
          free_cupt_sentence_iterator(&csi_tp);
          return 1;
        }
      }
    }
    free_cupt_sentence_iterator(&csi);
    if (cupt_paths_true_positives != NULL) {
      free_cupt_sentence_iterator(&csi_tp);
    }
  }

  size_t malloc_size = num_nodes * sizeof(struct graph_node);
  void *malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    goto malloc_fail;
  }
  memset(malloc_pointer, '\0', malloc_size);
  g->nodes = (struct graph_node *)malloc_pointer;
  g->num_nodes = num_nodes;

  double relative_proportion_sum = 0.0;

  uint64_t i = 0; // changes for iterative update -> disabled
  uint64_t j = 0;
  while (i < g->num_nodes && j < w2v->num_vectors) {
    if (w2v->keys[j].active_in_current_graph == 1) {
      g->nodes[i].num_dimensions = w2v->num_dimensions;
      g->nodes[i].vector.fp32 = w2v->keys[j].vector;
      g->nodes[i].relative_proportion = (double)w2v->keys[j].num_occurrences;
      g->nodes[i].absolute_proportion = (int32_t)w2v->keys[j].num_occurrences;
      relative_proportion_sum += (double)w2v->keys[j].num_occurrences;
      i++;
    }
    j++;
  }

  for (uint64_t k = 0; k < g->num_nodes; k++) {
    g->nodes[k].relative_proportion /= relative_proportion_sum;
  }

  return 0;

malloc_fail:
  perror("malloc failed\n");
  printf("errno: %i\n", errno);
  goto return_failure;

return_failure:
  return 1;
}

int32_t _weitzman(struct matrix *m, double *res) {
  int32_t argmin_dim1_index = 0;
  int32_t argmin_dim2_index = 1;
  int32_t found_a_value = 0;

  for (uint64_t i = 0; i < m->a; i++) {
    for (uint64_t j = i + 1; j < m->b; j++) {
      int32_t index_in_buffer = (i * m->b) + j;
      if (m->active[index_in_buffer] == 0) {
        continue;
      }
      switch (m->fp_mode) {
      case FP32:
        if (found_a_value == 0 || m->bfr.fp32[index_in_buffer] < m->bfr.fp32[(argmin_dim1_index * m->b) + argmin_dim2_index]) {
          argmin_dim1_index = i;
          argmin_dim2_index = j;
        }
        break;
      case FP64:
        if (found_a_value == 0 || m->bfr.fp64[index_in_buffer] < m->bfr.fp64[(argmin_dim1_index * m->b) + argmin_dim2_index]) {
          argmin_dim1_index = i;
          argmin_dim2_index = j;
        }
        break;
      }
      found_a_value = 1;
    }
  }
  if (!found_a_value) {
    return 0;
  }

  double local_res_a = -1.0;
  double local_res_b = -1.0;

  size_t malloc_size_a = m->a;
  size_t malloc_size_b = m->b;
  void *malloc_pointer;
  malloc_pointer = malloc(malloc_size_a);
  if (malloc_pointer == NULL) {
    goto malloc_failure;
  }
  int8_t *a_mem = (int8_t *)malloc_pointer;

  for (int32_t i = 0; i < (int32_t)m->a; i++) {
    int32_t index_in_buffer = (i * m->b) + argmin_dim2_index;
    a_mem[i] = m->active[index_in_buffer];
    m->active[index_in_buffer] = 0;
  }

  int32_t err = _weitzman(m, &local_res_a);
  if (err != 0) {
    goto _weitzman_failure;
  }
  for (int32_t i = 0; i < (int32_t)m->a; i++) {
    int32_t index_in_buffer = (i * m->b) + argmin_dim2_index;
    m->active[index_in_buffer] = a_mem[i];
  }
  free(a_mem);

  malloc_pointer = malloc(malloc_size_b);
  if (malloc_pointer == NULL) {
    goto malloc_failure;
  }
  int8_t *b_mem = (int8_t *)malloc_pointer;

  for (uint64_t j = 0; j < m->b; j++) {
    int32_t index_in_buffer = (argmin_dim1_index * m->b) + j;
    b_mem[j] = m->active[index_in_buffer];
    m->active[index_in_buffer] = 0;
  }
  err = _weitzman(m, &local_res_b);
  if (err != 0) {
    goto _weitzman_failure;
  }
  for (uint64_t j = 0; j < m->b; j++) {
    int32_t index_in_buffer = (argmin_dim1_index * m->b) + j;
    m->active[index_in_buffer] = b_mem[j];
  }
  free(b_mem);

  double result = 0.0;
  switch (m->fp_mode) {
  case FP32:
    result = (double)m->bfr.fp32[(argmin_dim1_index * m->b) + argmin_dim2_index];
    break;
  case FP64:
    result = (double)m->bfr.fp64[(argmin_dim1_index * m->b) + argmin_dim2_index];
    break;
  }
  if (local_res_a > local_res_b) {
    result += local_res_a;
  } else {
    result += local_res_b;
  }
  (*res) = result;

  return 0;

malloc_failure:
  perror("failed to malloc\n");
  return 1;

_weitzman_failure:
  perror("failed to call _weitzman\n");
  return 1;
}

int32_t weitzman_from_graph(struct graph *const g, double *const res, const int8_t fp_mode) {
  struct matrix m;
  int32_t err = create_matrix(&m, g->num_nodes, g->num_nodes, fp_mode);
  if (err != 0) {
    perror("failed to call create_matrix\n");
    return 1;
  }
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = i + 1; j < g->num_nodes; j++) {
      int32_t index_in_buffer = (i * m.b) + j;
      double distance = 0.0;
      switch (fp_mode) {
      case FP32:
#if ENABLE_AVX256 == 1
        distance =
            (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
        distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
        m.bfr.fp32[index_in_buffer] = (float)distance;
        break;
      case FP64:
        distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
        m.bfr.fp64[index_in_buffer] = distance;
        break;
      }
      m.active[index_in_buffer] = 1;
    }
  }

  err = _weitzman(&m, res);
  if (err != 0) {
    free_matrix(&m);
    perror("failed to call _weitzman\n");
    return 1;
  }

  free_matrix(&m);
  return 0;
}

int32_t double_cmp(const void *a, const void *b) {
  double a_ = *((double *)a);
  double b_ = *((double *)b);
  if (a_ < b_) {
    return -1;
  } else if (a_ > b_) {
    return 1;
  } else {
    return 0;
  }
}

int32_t _lexicographic(const struct matrix *const m, double *const res, long double *const res_hybrid) {
  long double _m = (long double)m->a;
  // long double c_m = (long double) (pow(M_PI, (_m / 2.0)) / lgamma((_m / 2.0) + 1.0));
  long double c_m = (long double)(pow(PI, (_m / 2.0)) / lgamma((_m / 2.0) + 1.0));

  int32_t argmin_dim1_index = 0;
  int32_t argmin_dim2_index = 1;
  int32_t found_a_value = 0;

  void *malloc_pointer;

  size_t malloc_size_champion;
  malloc_size_champion = m->b * sizeof(double);
  malloc_pointer = malloc(malloc_size_champion);
  if (malloc_pointer == NULL) {
    goto malloc_failure;
  }
  memset(malloc_pointer, '\0', malloc_size_champion);
  double *is_min_champion = (double *)malloc_pointer;
  int32_t length_champion = 0;
  int32_t i_champion = -1;
  int32_t j_champion = -1;

  size_t malloc_size_challenger;
  malloc_size_challenger = m->b * sizeof(double);
  malloc_pointer = malloc(malloc_size_challenger);
  if (malloc_pointer == NULL) {
    goto malloc_failure;
  }
  memset(malloc_pointer, '\0', malloc_size_challenger);
  double *is_min_challenger = (double *)malloc_pointer;
  int32_t length_challenger = 0;

  for (uint64_t i = 0; i < m->a; i++) {
    if (m->active[i * m->b + i] == 0) {
      continue;
    }
    memset(is_min_challenger, '\0', m->b * sizeof(double));
    length_challenger = 0;

    int32_t local_i_challenger = -1;
    int32_t local_j_challenger = -1;
    double local_d_challenger = -1.0;

    for (uint64_t j = 0; j < m->b; j++) { // necessary if we want proper challengers?
      if (i == j) {
        continue;
      }
      if (m->active[j * m->b + j] == 0) {
        continue;
      }
      int32_t index_in_buffer = (i * m->b) + j;
      if (m->active[index_in_buffer] == 0) {
        continue;
      }
      double distance = -1.0;
      switch (m->fp_mode) {
      case FP32:
        distance = (double)m->bfr.fp32[index_in_buffer];
        break;
      case FP64:
        distance = m->bfr.fp64[index_in_buffer];
        break;
      }

      if (length_challenger == 0 || distance < local_d_challenger) {
        local_i_challenger = i;
        local_j_challenger = j;
        local_d_challenger = distance;
      }

      is_min_challenger[length_challenger] = distance;
      length_challenger++;
    }

    if (length_challenger == 0) {
      continue;
    }
    qsort(is_min_challenger, length_challenger, sizeof(double), double_cmp);
    if (found_a_value == 0) {
      memcpy(is_min_champion, is_min_challenger, m->b * sizeof(double));
      length_champion = length_challenger;
      i_champion = local_i_challenger;
      j_champion = local_j_challenger;
      found_a_value = 1;
      continue;
    }

    for (int32_t k = 0; k < length_challenger && k < length_champion; k++) {
      if (is_min_challenger[k] == is_min_champion[k]) {
        continue;
      } else if (is_min_challenger[k] < is_min_champion[k]) {
        memcpy(is_min_champion, is_min_challenger, m->b * sizeof(double));
        length_champion = length_challenger;
        i_champion = local_i_challenger;
        j_champion = local_j_challenger;
        break;
      } else {
        break;
      }
    }
  }
  if (!found_a_value) {
    return 0;
  }

  argmin_dim1_index = i_champion;
  argmin_dim2_index = j_champion;

  m->active_final[(argmin_dim1_index * m->b) + argmin_dim2_index] = 1;

  double local_res_a = 0.0;
  long double local_res_hybrid_a = 0.0;

  // changing to dim1 as it is the point itself, not the closest point in the rest of graph
  for (uint64_t j = 0; j < m->a; j++) {
    m->active[argmin_dim1_index * m->b + j] = 0;
    m->active[j * m->b + argmin_dim1_index] = 0; // assuming square matrix
  }

  int32_t err = _lexicographic(m, &local_res_a, &local_res_hybrid_a);
  if (err != 0) {
    goto _lexicographic_failure;
  }

  double result = 0.0;
  switch (m->fp_mode) {
  case FP32:
    result = (double)m->bfr.fp32[(argmin_dim1_index * m->b) + argmin_dim2_index];
    break;
  case FP64:
    result = (double)m->bfr.fp64[(argmin_dim1_index * m->b) + argmin_dim2_index];
    break;
  }
  long double result_hybrid = c_m * ((long double)pow(result, (double)_m));
  result += local_res_a;
  result_hybrid += local_res_hybrid_a;

  (*res) = result;
  (*res_hybrid) = result_hybrid;

  free(is_min_challenger);
  free(is_min_champion);

  return 0;

malloc_failure:
  perror("failed to malloc\n");
  return 1;

_lexicographic_failure:
  perror("failed to call _lexicographic\n");
  return 1;
}

int32_t lexicographic_from_graph(struct graph *const g, double *const res, long double *const res_hybrid, const int8_t fp_mode,
                                 const struct matrix *const m_) {
  if (m_ != NULL) {
    memset(m_->active, 1, m_->a * m_->b * sizeof(int8_t));
    memset(m_->active_final, '\0', m_->a * m_->b * sizeof(int8_t));
    if (_lexicographic(m_, res, res_hybrid) != 0) {
      perror("failed to call _lexicographic with matrix given as argument\n");
      return 1;
    }
    return 0;
  } else {
    struct matrix m;
    int32_t err = create_matrix(&m, g->num_nodes, g->num_nodes, fp_mode);
    if (err != 0) {
      perror("failed to call create_matrix\n");
      return 1;
    }

    for (uint64_t i = 0; i < g->num_nodes; i++) {
      for (uint64_t j = i + 1; j < g->num_nodes; j++) {
        int32_t index_in_buffer = (i * m.b) + j;
        double distance = 0.0;
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance =
              (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
          distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
          m.bfr.fp32[index_in_buffer] = (float)distance;
          break;
        case FP64:
          distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          m.bfr.fp64[index_in_buffer] = distance;
          break;
        }
        m.active[index_in_buffer] = 1;
        m.active_final[index_in_buffer] = 0;
      }
    }

    err = _lexicographic(&m, res, res_hybrid);
    if (err != 0) {
      free_matrix(&m);
      perror("failed to call _lexicographic\n");
      return 1;
    }

    free_matrix(&m);
  }
  return 0;
}

int32_t stirling_from_graph(struct graph *g, double *restrict const result, const double alpha_arg, const double beta_arg,
                            const int8_t fp_mode, const struct matrix *restrict const m_) {
  double local_result = 0.0;

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      if (i == j) {
        continue;
      }
      double distance = 0.0;
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          distance = (double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = m_->bfr.fp64[i * m_->b + j];
          break;
        }
      } else {
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance =
              (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
          distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
          break;
        case FP64:
          distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          break;
        }
      }
      double proportion_product = g->nodes[i].relative_proportion * g->nodes[j].relative_proportion;
      local_result += pow(distance, alpha_arg) * pow(proportion_product, beta_arg);
    }
  }

  (*result) = local_result;

  return 0;
}

int32_t ricotta_szeidl_from_graph(struct graph *const g, double *const result, const double alpha_arg, const int8_t fp_mode,
                                  const struct matrix *const m_) {
  double local_result = 0.0;

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    double local_neg_sum = 1.0;
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      if (i == j) {
        continue;
      }
      double distance = 0.0;
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          distance = (double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = m_->bfr.fp64[i * m_->b + j];
          break;
        }
      } else {
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance =
              (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
          distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
          break;
        case FP64:
          distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          break;
        }
      }
      if (distance < 0.0) {
        printf("distance is negative: %f\n", distance);
      }
      local_neg_sum -= distance * g->nodes[j].relative_proportion;
    }
    double product = g->nodes[i].relative_proportion * pow(local_neg_sum, alpha_arg - 1.0);
    if (isnan(product)) {
      printf("product is nan; relative_proportion: %f, pow: %f, local_neg_sum: %f, alpha[k] - 1.0: %f\r",
             g->nodes[i].relative_proportion, pow(local_neg_sum, alpha_arg - 1.0), local_neg_sum, alpha_arg - 1.0);
      continue;
    }
    local_result += product;
  }

  local_result = (1.0 - local_result) / (alpha_arg - 1.0);
  (*result) = local_result;
  return 0;
}

int32_t chao_et_al_functional_diversity_from_graph(struct graph *const g, double *const div_result, double *const hill_result,
                                                   const double alpha, const int8_t fp_mode, const struct matrix *const m_) {
  // see Chao et al. (2014)

  const double LOGARITHMIC_BASE = E;

  double rao_q = 0.0;
  double *matrix = NULL;
  if (m_ == NULL) {
    void *malloc_pointer = NULL;
    size_t malloc_size = g->num_nodes * g->num_nodes * sizeof(double);
    malloc_pointer = malloc(malloc_size);
    if (malloc_pointer == NULL) {
      perror("failed to malloc\n");
      return 1;
    }
    memset(malloc_pointer, '\0', malloc_size);
    matrix = (double *)malloc_pointer;
  }

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      double distance = 0.0;
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          distance = (double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = m_->bfr.fp64[i * m_->b + j];
          break;
        default:
          perror("unknown FP mode\n");
          return 1;
        }
      } else {
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance =
              (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
          distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
          break;
        case FP64:
          distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          break;
        default:
          perror("unknown FP mode\n");
          return 1;
        }
      }
      if (m_ == NULL) {
        matrix[(i * g->num_nodes) + j] = distance;
      }
      rao_q += distance * g->nodes[i].relative_proportion * g->nodes[j].relative_proportion;
    }
  }

  double diversity = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      double distance;
      if (m_ == NULL) {
        distance = matrix[(i * g->num_nodes) + j];
      } else {
        switch (m_->fp_mode) {
        case FP32:
          distance = (double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = m_->bfr.fp64[i * m_->b + j];
          break;
        default:
          perror("unknown FP mode\n");
          return 1;
        }
      }
      if (alpha != 1.0) {
        diversity += distance * pow((g->nodes[i].relative_proportion * g->nodes[j].relative_proportion) / rao_q, alpha);
      } else {
        double product_ratio = (g->nodes[i].relative_proportion * g->nodes[j].relative_proportion) / rao_q;
        diversity += distance * product_ratio * (log(product_ratio) / log(LOGARITHMIC_BASE));
      }
    }
  }
  if (alpha != 1.0) {
    diversity = pow(diversity, 1.0 / (1.0 - alpha));
  } else {
    diversity *= -1.0;
    diversity = pow(LOGARITHMIC_BASE, diversity);
  }

  double hill_number = pow(diversity / rao_q, 0.5);

  (*div_result) = diversity;
  (*hill_result) = hill_number;

  free(matrix);

  return 0;
}

int32_t leinster_cobbold_diversity_from_graph(struct graph *const g, double *const div_result, double *const hill_result,
                                              const double alpha, const int8_t fp_mode, const struct matrix *const m_) {
  // see Leinster & Cobbold (2012)

  const double LOGARITHMIC_BASE = E;

  const double u = 1.0;

  double hill_number;
  if (alpha != 1.0) {
    hill_number = 0.0;
  } else {
    hill_number = 1.0;
  }

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    double local_agg = 0.0;
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      double distance = 0.0;
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          distance = (double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = m_->bfr.fp64[i * m_->b + j];
          break;
        }
      } else {
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance =
              (double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#else
          distance = (double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif
          break;
        case FP64:
          distance = cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          break;
        }
      }

      double similarity = 1.0 - distance;

      local_agg += g->nodes[j].relative_proportion * pow(E, -u * similarity);
    }
    if (alpha != 1.0) {
      hill_number += pow(local_agg, alpha - 1.0);
    } else {
      hill_number *= pow(local_agg, g->nodes[i].relative_proportion);
    }
  }
  if (alpha != 1.0) {
    hill_number = pow(hill_number, 1.0 / (1.0 - alpha));
  } else {
    hill_number = pow(hill_number, -1.0);
  }

  double entropy = log(hill_number) / log(LOGARITHMIC_BASE);

  (*div_result) = entropy;
  (*hill_result) = hill_number;

  return 0;
}

int32_t scheiner_species_phylogenetic_functional_diversity_from_graph(struct graph *const g, double *const div_result,
                                                                      double *const hill_result, const double alpha,
                                                                      const int8_t fp_mode, const struct matrix *const m_) {
  // see Scheiner (2012)

  const long double LOGARITHMIC_BASE = E;

  const int8_t USE_LOGARITHM = 0;

  void *malloc_pointer = NULL;
  size_t malloc_size = g->num_nodes * sizeof(long double);
  malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    perror("failed to malloc\n");
    return 1;
  }
  memset(malloc_pointer, '\0', malloc_size);
  long double *vector = (long double *)malloc_pointer;
  long double norm_sum = 0.0;

  long double m = (long double)g->nodes[0].num_dimensions;
  // long double c_m = (long double) (pow(M_PI, (m / 2.0)) / lgamma((m / 2.0) + 1.0));
  long double c_m = (long double)(pow(PI, (m / 2.0)) / lgamma((m / 2.0) + 1.0));
  // long double c_m_ln = log(c_m);

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    int32_t index_min_distance = -1;
    long double min_distance = -1.0;
    for (uint64_t j = 0; j < g->num_nodes; j++) {
      if (i == j) {
        continue;
      }
      long double distance = 0.0;
      if (m_ != NULL) {
        switch (m_->fp_mode) {
        case FP32:
          distance = (long double)m_->bfr.fp32[i * m_->b + j];
          break;
        case FP64:
          distance = (long double)m_->bfr.fp64[i * m_->b + j];
          break;
        }
      } else {
        // the paper makes use of euclidean distance
        switch (fp_mode) {
        case FP32:
#if ENABLE_AVX256 == 1
          distance = (long double)cosine_distance_fp32_avx(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32,
                                                           g->nodes[i].num_dimensions);
#else
          distance =
              (long double)cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->nodes[i].num_dimensions);
#endif

          // distance = (long double) minkowski_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, 2.0,
          // g->nodes[i].num_dimensions);
          break;
        case FP64:
          distance = (long double)cosine_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, g->nodes[i].num_dimensions);
          // distance = minkowski_distance(g->nodes[i].vector.fp64, g->nodes[j].vector.fp64, 2.0, g->nodes[i].num_dimensions);
          break;
        }
      }

      if (index_min_distance == -1 || distance < min_distance) {
        index_min_distance = j;
        min_distance = distance;
      }
    }

    // long double v_i = c_m * pow(min_distance, m);
    long double d_pow = (long double)pow(min_distance, m);
    // long double d_pow_ln = log(d_pow);
    // long double v_i_ln = c_m_ln + d_pow_ln;
    // long double v_i = pow(M_E, v_i_ln);
    long double v_i = (long double)(c_m * d_pow);

    long double abundance = (long double)g->nodes[i].absolute_proportion;
    // long double abundance_ln = log(abundance);

    long double phylogenetic_distance = 1.0; // !
    // long double phylogenetic_distance_ln = log(phylogenetic_distance);

    // matrix[(i * g->num_nodes) + j] = abundance * phylogenetic_distance * v_i;
    // vector[i] = abundance * phylogenetic_distance * v_i;
    // vector[i] = abundance_ln + phylogenetic_distance_ln + v_i_ln;
    long double local_agg;
    /*
    if(USE_LOGARITHM){
            local_agg = abundance_ln + phylogenetic_distance_ln + v_i_ln;
    } else {
            local_agg = abundance * phylogenetic_distance * v_i;
    }
    */
    local_agg = abundance * phylogenetic_distance * v_i;
    vector[i] = local_agg;
    norm_sum += local_agg;
  }

  long double diversity = 0.0;
  // modification, not like Chao et al. (2014) (?)
  long double hill_number = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    long double distance = vector[i];
    long double product_ratio = distance / norm_sum;
    if (USE_LOGARITHM) {
      product_ratio = pow(LOGARITHMIC_BASE, product_ratio); // change
    }
    if (alpha != 1.0) {
      hill_number += pow(product_ratio, alpha);
    } else {
      diversity += product_ratio * (log(product_ratio) / log(LOGARITHMIC_BASE));
    }
  }
  if (alpha != 1.0) {
    hill_number = pow(hill_number, 1.0 / (1.0 - alpha));
    diversity = pow(hill_number, 1.0 / alpha);
  } else {
    diversity *= -1.0;
    hill_number = pow(LOGARITHMIC_BASE, diversity);
  }

  (*div_result) = (double)diversity;
  (*hill_result) = (double)hill_number;

  free(vector);

  return 0;
}

int32_t nhc_e_q_from_graph(struct graph *const g, double *const res_nhc, double *const res_e_q) {
  size_t alloc_size = g->num_nodes * sizeof(double);
  double *proportions = malloc(alloc_size);
  if (proportions == NULL) {
    perror("malloc failed\n");
    return 1;
  }
  for (uint32_t i = 0; i < g->num_nodes; i++) {
    proportions[i] = (double)g->nodes[i].absolute_proportion; // abundances, not relative
  }

  qsort(proportions, g->num_nodes, sizeof(double), double_cmp);

  // reverse order
  for (uint32_t i = 0; i < floor(((double)g->num_nodes) / 2.0); i++) {
    double placeholder = proportions[i];
    proportions[i] = proportions[g->num_nodes - 1 - i];
    proportions[g->num_nodes - 1 - i] = placeholder;
  }

  for (uint32_t i = 0; i < g->num_nodes; i++) {
    proportions[i] = log(proportions[i]);
  }

  double diff_min_max = 100.0;
  const double update_divider = 10.0;
  double min_b = -diff_min_max;
  double max_b = 0.0;
  double min_b_prime = -diff_min_max;
  double max_b_prime = diff_min_max;
  const uint32_t total_subiter = 32;
  const uint32_t total_iter = 8;
  double least_mse = 0.0;
  double best_b = 0.0;
  double least_mse_prime = 0.0;
  double best_b_prime = 0.0;
  for (uint32_t current_iter = 0; current_iter < total_iter; current_iter++) {
    least_mse = 0.0;
    best_b = 0.0;
    least_mse_prime = 0.0;
    best_b_prime = 0.0;
    for (uint32_t current_subiter = 0; current_subiter < total_subiter; current_subiter++) {
      double b = min_b + (max_b - min_b) * (((double)current_subiter) / ((double)(total_subiter - 1)));
      double mse = 0.0;
      double b_prime = min_b_prime + (max_b_prime - min_b_prime) * (((double)current_subiter) / ((double)(total_subiter - 1)));
      double mse_prime = 0.0;
      for (uint32_t i = 0; i < g->num_nodes; i++) {
        double rank = ((double)i) + 1.0;
        double real_value = proportions[i] / rank; // log of abundance on rank of abundance
        double prediction = b * rank;
        mse += pow(real_value - prediction, 2.0);

        double rank_scaled = rank / ((double)g->num_nodes);
        double real_value_prime = rank_scaled / proportions[i]; // scaled rank over log of abundance
        double prediction_prime = b_prime * rank;
        mse_prime += pow(real_value_prime - prediction_prime, 2.0);
      }
      if (current_subiter == 0 || mse < least_mse) {
        least_mse = mse;
        best_b = b;
      }
      if (current_subiter == 0 || mse_prime < least_mse_prime) {
        least_mse_prime = mse_prime;
        best_b_prime = b_prime;
      }
    }
    diff_min_max /= update_divider;
    min_b = best_b - diff_min_max / 2.0;
    max_b = best_b + diff_min_max / 2.0;
    min_b_prime = best_b_prime - diff_min_max / 2.0;
    max_b_prime = best_b_prime + diff_min_max / 2.0;
  }

  (*res_nhc) = best_b;
  (*res_e_q) = -(2.0 / PI) * atan(best_b_prime);

  free(proportions);

  return 0;
}
