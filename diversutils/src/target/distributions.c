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

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "distributions.h"

#ifndef FP_MODES
#define FP_MODES
enum { FP32, FP64 };
#endif

int32_t create_zipfian_distribution(struct zipfian_distribution *z, double s, uint32_t n, const int8_t fp_mode) {
  size_t malloc_size = n;
  switch (fp_mode) {
  case FP32:
    malloc_size *= 4;
    break;
  case FP64:
    malloc_size *= 8;
    break;
  }
  void *malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    perror("failed to malloc\n");
    return 1;
  }
  memset(malloc_pointer, '\0', malloc_size);
  switch (fp_mode) {
  case FP32:
    z->vector.fp32 = (float *)malloc_pointer;
    break;
  case FP64:
    z->vector.fp64 = (double *)malloc_pointer;
    break;
  }

  double p_sum = 0.0;
  for (uint32_t i = 1; i <= n; i++) {
    double p = pow(((double)i), -s);
    switch (fp_mode) {
    case FP32:
      z->vector.fp32[i - 1] = p;
      break;
    case FP64:
      z->vector.fp64[i - 1] = p;
      break;
    }
    p_sum += p;
  }
  for (uint32_t i = 0; i < n; i++) {
    switch (fp_mode) {
    case FP32:
      z->vector.fp32[i] /= p_sum;
      break;
    case FP64:
      z->vector.fp64[i] /= p_sum;
      break;
    }
  }

  z->s = s;
  z->n = n;

  return 0;
}

void free_zipfian_distribution(struct zipfian_distribution *z, const int8_t fp_mode) {
  switch (fp_mode) {
  case FP32:
    free(z->vector.fp32);
    break;
  case FP64:
    free(z->vector.fp64);
    break;
  }
}

double mean_squared_error(double *v, double *w, uint32_t n) {
  double sum = 0.0;
  for (uint32_t i = 0; i < n; i++) {
    sum += pow(w[i] - v[i], 2.0);
  }
  return (sum / ((double)n));
}

int32_t double_cmp_reverse(const void *a, const void *b) {
  double a_ = (*((double *)a));
  double b_ = (*((double *)b));
  if (a_ > b_) {
    return -1;
  } else if (a_ < b_) {
    return 1;
  }
  return 0;
}

int32_t zipfian_fit(double *v, uint32_t n, double *result) {
  const int32_t num_precision_levels = 8;
  const int32_t num_iter = 32;
  double lower_bound = 0.0;
  double upper_bound = 10.0;
  double best_s = -1.0;
  double best_mse = -1.0;
  double *v_copy;
  size_t alloc_size;
  /*
  time_t t_start, t_end;

  t_start = time(NULL);
  */

  alloc_size = n * sizeof(double);
  v_copy = malloc(alloc_size);
  if (v_copy == NULL) {
    perror("failed to malloc\n");
    return 1;
  }
  memcpy(v_copy, v, alloc_size);

  qsort((void *)v_copy, n, sizeof(double), double_cmp_reverse);

  for (int32_t i = 0; i < num_precision_levels; i++) {
    double window = (upper_bound - lower_bound);
    double step = window / (((double)num_iter) - 1.0);
    best_s = -1.0;
    best_mse = -1.0;
    for (int32_t j = 0; j < num_iter; j++) {
      double s = lower_bound + step * ((double)j);
      struct zipfian_distribution z;
      int32_t err = create_zipfian_distribution(&z, s, n, FP64);
      if (err != 0) {
        perror("failed to call create_zipfian_distribution\n");
        free(v_copy);
        return 1;
      }

      double mse = mean_squared_error(v_copy, z.vector.fp64, n);
      if (j == 0 || mse < best_mse) {
        best_s = s;
        best_mse = mse;
      }

      free_zipfian_distribution(&z, FP64);
    }

    lower_bound = best_s - ((window / 10.0) / 2.0);
    upper_bound = best_s + ((window / 10.0) / 2.0);
    if (lower_bound < 0.0) {
      lower_bound = 0.0;
    }
  }

  (*result) = best_s;

  free(v_copy);

  /*
  t_end = time(NULL);
  printf("Zipfian curvature computed in %lus\n", t_end - t_start);
  */

  return 0;
}

int32_t zipfian_fit_from_graph(struct graph *g, double *result) {
  size_t malloc_size = g->num_nodes * sizeof(double);
  void *malloc_pointer = malloc(malloc_size);
  if (malloc_pointer == NULL) {
    perror("failed to malloc\n");
    return 1;
  }
  memset(malloc_pointer, '\0', malloc_size);
  double *series = (double *)malloc_pointer;

  for (uint64_t i = 0; i < g->num_nodes; i++) {
    series[i] = g->nodes[i].relative_proportion;
  }

  int32_t err = zipfian_fit(series, g->num_nodes, result);
  if (err != 0) {
    perror("failed to call zipfian_fit\n");
    free(series);
    return 1;
  }

  free(series);

  return 0;
}
