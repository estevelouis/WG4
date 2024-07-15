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

#ifndef _POSIX_C_SOURCE
// for clock_gettime in time.h: https://stackoverflow.com/questions/69145941/why-does-clock-gettime-not-compile-when-using-c99
// #define _POSIX_C_SOURCE >= 199309L
#define _POSIX_C_SOURCE 199309L
#endif

#include <time.h>

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cpu.h"
#include "dfunctions.h"
#include "distributions.h"
#include "graph.h"
#include "logging.h"
#include "measurement.h"
#include "sorted_array/array.h"
#include "stats.h"

int32_t time_ns_delta(int64_t *const delta) {
  static struct timespec static_ts = {0};

  struct timespec new_ts;

  if (clock_gettime(CLOCK_MONOTONIC, &new_ts) != 0) {
    perror("Failed to call clock_gettime\n");
    return 1;
  }

  if (delta != NULL) {
    (*delta) = (new_ts.tv_sec - static_ts.tv_sec) * 1000000000 + (new_ts.tv_nsec - static_ts.tv_nsec);
  }

  static_ts = new_ts;

  return 0;
}

int32_t virtual_memory_consumption(int64_t *const res) {
  FILE *f;
  const int32_t bfr_size = 2048;
  char bfr[bfr_size];
  char *strtok_pointer;

  if (res == NULL) {
    perror("Cannot pass NULL pointer to virtual_memory_consumption\n");
    return 1;
  }

  memset(bfr, '\0', bfr_size * sizeof(char));

  f = fopen("/proc/self/stat", "r");

  if (f == NULL) {
    perror("Failed to open /proc/self/stat\n");
    return 1;
  }

  if (fgets(bfr, bfr_size, f) == 0) {
    perror("No bytes could be read from /proc/self/stat\n");
    fclose(f);
    return 1;
  }

  strtok_pointer = strtok(bfr, " ");
  for (int32_t i = 0; i < 23; i++) {
    strtok_pointer = strtok(NULL, " ");
  }
  (*res) = strtol(strtok_pointer, NULL, 10);

  fclose(f);

  return 0;
}

void timing_and_memory(FILE *f_timing_ptr, FILE *f_memory_ptr, const uint8_t enable_output_timing,
                       const uint8_t enable_output_memory) {
  int64_t ns_delta, virtual_mem;
  if (enable_output_timing) {
    if (time_ns_delta(&ns_delta) != 0) {
      perror("Failed to call time_ns_delta\n");
      exit(1);
    } else {
      fprintf(f_timing_ptr, "\t%li", ns_delta);
    }
  }
  if (enable_output_memory) {
    if (virtual_memory_consumption(&virtual_mem) != 0) {
      perror("Failed to call virtual_memory_consumption\n");
      exit(1);
    } else {
      fprintf(f_memory_ptr, "\t%li", virtual_mem);
    }
  }
}

int32_t wrap_diversity_1r_0a(struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr,
                             FILE *f_timing_ptr, FILE *f_memory_ptr,
                             int32_t (*df)(struct graph *const, double *const, const int8_t, const struct matrix *const),
                             const uint8_t enable_timings, const uint8_t enable_output_timing,
                             const uint8_t enable_output_memory) {
  double res1;
  time_t t, delta_t;
  const int32_t bfr_size = 64;
  char bfr[bfr_size];
  t = time(NULL);
  if (df(g, &res1, fp_mode, m) != 0) {
    perror("Failed to call diversity function in wrap_diversity_1r_0a\n");
    return EXIT_FAILURE;
  }
  delta_t = time(NULL) - t;
  // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
  memset(bfr, '\0', bfr_size * sizeof(char));
  snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
  if (enable_timings) {
    info_format(__FILE__, __func__, __LINE__, bfr);
  }
  fprintf(f_ptr, "\t%.10e", res1);
  timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
  return 0;
}

int32_t wrap_diversity_1r_0a_no_matrix(struct graph *const g, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr,
                                       FILE *f_memory_ptr, int32_t (*df)(struct graph *const, double *const, const int8_t),
                                       const uint8_t enable_timings, const uint8_t enable_output_timing,
                                       const uint8_t enable_output_memory) {
  double res1;
  time_t t, delta_t;
  const int32_t bfr_size = 64;
  char bfr[bfr_size];
  t = time(NULL);
  if (df(g, &res1, fp_mode) != 0) {
    perror("Failed to call diversity function in wrap_diversity_1r_0a\n");
    return EXIT_FAILURE;
  }
  delta_t = time(NULL) - t;
  // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
  memset(bfr, '\0', bfr_size * sizeof(char));
  snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
  if (enable_timings) {
    info_format(__FILE__, __func__, __LINE__, bfr);
  }
  fprintf(f_ptr, "\t%.10e", res1);
  timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
  return 0;
}

int32_t wrap_diversity_2r_1a(struct graph *const g, struct matrix *const m, const double alpha, const int8_t fp_mode,
                             FILE *f_ptr, FILE *f_timing_ptr, FILE *f_memory_ptr,
                             int32_t (*df)(struct graph *const, double *const, double *const, double, const int8_t,
                                           const struct matrix *const),
                             const uint8_t enable_timings, const uint8_t enable_output_timing,
                             const uint8_t enable_output_memory) {
  double res1, res2;
  time_t t, delta_t;
  const int32_t bfr_size = 64;
  char bfr[bfr_size];
  t = time(NULL);
  if (df(g, &res1, &res2, alpha, fp_mode, m) != 0) {
    perror("Failed to call diversity function in wrap_diversity_2r_1a\n");
    return EXIT_FAILURE;
  }
  delta_t = time(NULL) - t;
  // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
  memset(bfr, '\0', bfr_size * sizeof(char));
  snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
  if (enable_timings) {
    info_format(__FILE__, __func__, __LINE__, bfr);
  }
  fprintf(f_ptr, "\t%.10e\t%.10e", res1, res2);
  timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
  return 0;
}

int32_t
wrap_diversity_2r_0a(struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr,
                     FILE *f_memory_ptr,
                     int32_t (*df)(struct graph *const, double *const, double *const, const int8_t, const struct matrix *const),
                     const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory) {
  double res1, res2;
  time_t t, delta_t;
  const int32_t bfr_size = 64;
  char bfr[bfr_size];
  t = time(NULL);
  if (df(g, &res1, &res2, fp_mode, m) != 0) {
    perror("Failed to call diversity function in wrap_diversity_2r_0a\n");
    return EXIT_FAILURE;
  }
  delta_t = time(NULL) - t;
  // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
  memset(bfr, '\0', bfr_size * sizeof(char));
  snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
  if (enable_timings) {
    info_format(__FILE__, __func__, __LINE__, bfr);
  }
  fprintf(f_ptr, "\t%.10e\t%.10e", res1, res2);
  timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
  return 0;
}

int32_t wrap_diversity_2r_0a_long_double_alt(
    struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr, FILE *f_memory_ptr,
    int32_t (*df)(struct graph *const, double *const, long double *const, const int8_t, const struct matrix *const),
    const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory) {
  double res1;
  long double res2;
  time_t t, delta_t;
  const int32_t bfr_size = 64;
  char bfr[bfr_size];
  t = time(NULL);
  if (df(g, &res1, &res2, fp_mode, m) != 0) {
    perror("Failed to call diversity function in wrap_diversity_2r_0a\n");
    return EXIT_FAILURE;
  }
  delta_t = time(NULL) - t;
  // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
  memset(bfr, '\0', bfr_size * sizeof(char));
  snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
  if (enable_timings) {
    info_format(__FILE__, __func__, __LINE__, bfr);
  }
  fprintf(f_ptr, "\t%.10e\t%.10Le", res1, res2);
  timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
  return 0;
}

/*
int32_t wrap_diversity_1r_0a(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE*
f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double*, const int8_t, struct matrix* const), const uint8_t
enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){ double res1; time_t t, delta_t; const
int32_t bfr_size = 64; char bfr[bfr_size]; t = time(NULL); if(df(g, &res1, fp_mode, m) != 0){perror("Failed to call diversity
function in wrap_diversity_1r_0a\n"); return EXIT_FAILURE;} delta_t = time(NULL) - t;
        // if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
        memset(bfr, '\0', bfr_size * sizeof(char));
        snprintf(bfr, bfr_size, "Computed df in %lis", delta_t);
        if(enable_timings){info_format(__FILE__, __func__, __LINE__, bfr);}
        fprintf(f_ptr, "\t%.10e", res1);
        timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
        return 0;
}
*/

int32_t apply_diversity_functions_to_graph(const uint64_t i, struct measurement_configuration *const mcfg,
                                           struct measurement_structure_references *const sref,
                                           struct measurement_mutables *const mmut) {
  const uint8_t enable_distance_computation =
      mcfg->enable.disparity_functions &&
      (mcfg->enable.stirling || mcfg->enable.ricotta_szeidl || mcfg->enable.pairwise ||
       mcfg->enable.chao_et_al_functional_diversity || mcfg->enable.scheiner_species_phylogenetic_functional_diversity ||
       mcfg->enable.leinster_cobbold_diversity || mcfg->enable.lexicographic || mcfg->enable.functional_evenness ||
       mcfg->enable.functional_dispersion || mcfg->enable.functional_divergence_modified);

  int32_t err;
  // if(enable_iterative_distance_computation){
  if (mcfg->threading.enable_iterative_distance_computation) {
    size_t local_malloc_size = sref->g->num_nodes * sizeof(float) * mcfg->threading.row_generation_batch_size;
    float *vector_batch = (float *)malloc(local_malloc_size);
    if (vector_batch == NULL) {
      goto malloc_fail;
    }
    memset(vector_batch, '\0', local_malloc_size);

    local_malloc_size = sref->g->num_nodes * sizeof(uint8_t) * mcfg->threading.row_generation_batch_size;
    uint8_t *used = (uint8_t *)malloc(local_malloc_size);
    if (used == NULL) {
      goto malloc_fail;
    }
    memset(used, '\0', local_malloc_size);

    struct iterative_state_pairwise_from_graph iter_state_pairwise;
    struct iterative_state_stirling_from_graph iter_state_stirling;
    struct iterative_state_leinster_cobbold_from_graph iter_state_leinster_cobbold;
    if (create_iterative_state_pairwise_from_graph(&iter_state_pairwise, sref->g) != 0) {
      perror("failled to call create_iterative_state_pairwise_from_graph\n");
      return 1;
    }
    if (create_iterative_state_stirling_from_graph(&iter_state_stirling, sref->g, mcfg->div_param.stirling_alpha,
                                                   mcfg->div_param.stirling_beta) != 0) {
      perror("failled to call create_iterative_state_stirling_from_graph\n");
      return 1;
    }
    if (create_iterative_state_leinster_cobbold_from_graph(&iter_state_leinster_cobbold, sref->g,
                                                           mcfg->div_param.leinster_cobbold_diversity_alpha) != 0) {
      perror("failled to call create_iterative_state_leinster_cobbold_from_graph\n");
      return 1;
    }

    int32_t i_index = 0;
    double sum = 0.0;
    for (uint64_t h = 0; h < sref->g->num_nodes; h += mcfg->threading.row_generation_batch_size) {
      if (mcfg->threading.enable_multithreaded_row_generation) {
        distance_row_batch_from_graph_multithread(sref->g, i_index, vector_batch, mcfg->threading.num_row_threads,
                                                  mcfg->threading.row_generation_batch_size);
      } else {
        perror("single-threaded version of row computation has been disabled\n");
        return 1;
      }

      /* // DO NOT REMOVE
      // for(uint64_t m = h ; m < sref->g->num_nodes ; m++){
      for(uint64_t m = 0 ; m < row_generation_batch_size && h + m < sref->g->num_nodes ; m++){
              for(uint64_t p = i_index + 1 ; p < sref->g->num_nodes ; p++){
                      sum += (double) vector_batch[m * sref->g->num_nodes + p];
              }

              if(mcfg->enable.stirling){iterate_iterative_state_stirling_from_graph(&iter_state_stirling, &(vector_batch[m *
      sref->g->num_nodes]));} if(mcfg->enable.pairwise){iterate_iterative_state_pairwise_from_graph(&iter_state_pairwise,
      &(vector_batch[m * sref->g->num_nodes]));}

              i_index = h + m + 1;
              iter_state_stirling.i = i_index;
              iter_state_pairwise.i = i_index;
      }
      */

      uint64_t actual_row_generation_batch_size = mcfg->threading.row_generation_batch_size;
      if (h + mcfg->threading.row_generation_batch_size >= sref->g->num_nodes) {
        // if(sref->g->num_nodes - h < row_generation_batch_size){
        actual_row_generation_batch_size = sref->g->num_nodes - h;
      }

      if (mcfg->enable.pairwise) {
        pthread_t agg_threads[mcfg->threading.row_generation_batch_size]; // still have full buffer to make it writeable?
        struct thread_args_aggregator agg_thread_args[mcfg->threading.row_generation_batch_size];
        i_index = 0; // ?
        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          struct thread_args_aggregator local_args = (struct thread_args_aggregator){
              .iter_state.pairwise = &iter_state_pairwise,
              .i = i_index,
              .vector = &(vector_batch[m * sref->g->num_nodes]),
          };
          memcpy(&(agg_thread_args[m]), &local_args, sizeof(struct thread_args_aggregator));

#if ENABLE_AVX256 == 1
          if (pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_pairwise_from_graph_avx_thread,
                             &(agg_thread_args[m])) != 0) {
#else
          if (pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_pairwise_from_graph_thread,
                             &(agg_thread_args[m])) != 0) {
#endif
            perror("Failed to call pthread_create\n");
            return 1;
          }

          i_index = h + m + 1;             // !
          iter_state_pairwise.i = i_index; // !
        }

        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          pthread_join(agg_threads[m], NULL);
        }
      }
      if (mcfg->enable.stirling) {
        pthread_t agg_threads[mcfg->threading.row_generation_batch_size];
        struct thread_args_aggregator agg_thread_args[mcfg->threading.row_generation_batch_size];
        i_index = 0; // ?
        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          struct thread_args_aggregator local_args = (struct thread_args_aggregator){
              .iter_state.stirling = &iter_state_stirling,
              .i = i_index,
              .vector = &(vector_batch[m * sref->g->num_nodes]),
          };
          memcpy(&(agg_thread_args[m]), &local_args, sizeof(struct thread_args_aggregator));

#if ENABLE_AVX256 == 1
          // if(pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_stirling_from_graph_avx_thread,
          // &(agg_thread_args[m])) != 0){
          perror("iterate_iterative_state_stirling_from_graph_avx_thread not implemented\n");
          return 1;
#else
          if (pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_stirling_from_graph_thread,
                             &(agg_thread_args[m])) != 0) {
            perror("Failed to call pthread_create\n");
            return 1;
          }
#endif

          i_index = h + m + 1;             // !
          iter_state_stirling.i = i_index; // !
        }

        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          pthread_join(agg_threads[m], NULL);
        }
      }
      if (mcfg->enable.leinster_cobbold_diversity) {
        pthread_t agg_threads[mcfg->threading.row_generation_batch_size];
        struct thread_args_aggregator agg_thread_args[mcfg->threading.row_generation_batch_size];
        i_index = 0; // ?
        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          struct thread_args_aggregator local_args = (struct thread_args_aggregator){
              .iter_state.leinster_cobbold = &iter_state_leinster_cobbold,
              .i = i_index,
              .vector = &(vector_batch[m * sref->g->num_nodes]),
          };
          memcpy(&(agg_thread_args[m]), &local_args, sizeof(struct thread_args_aggregator));

#if ENABLE_AVX256 == 1
          // if(pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_leinster_cobbold_from_graph_avx_thread,
          // &(agg_thread_args[m])) != 0){
          perror("iterate_iterative_state_leinster_cobbold_from_graph_avx_thread not implemented\n");
          return 1;
#else
          if (pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_leinster_cobbold_from_graph_thread,
                             &(agg_thread_args[m])) != 0) {
            perror("Failed to call pthread_create\n");
            return 1;
          }
#endif

          i_index = h + m + 1;                     // !
          iter_state_leinster_cobbold.i = i_index; // !
        }

        for (uint64_t m = 0; m < actual_row_generation_batch_size; m++) {
          pthread_join(agg_threads[m], NULL);
        }
      }
    }
    free(vector_batch);
    free(used);

    finalise_iterative_state_pairwise_from_graph(&iter_state_pairwise);
    finalise_iterative_state_stirling_from_graph(&iter_state_stirling);
    finalise_iterative_state_leinster_cobbold_from_graph(&iter_state_leinster_cobbold);

    double mu_dist = sum / ((double)(sref->g->num_nodes * (sref->g->num_nodes - 1) / 2));

    // fprintf(mcfg->io.f_ptr, "%lu\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%c", i+1, mmut->sentence.num,
    // mmut->sentence.num_all, mmut->document.num, mcfg->io.w2v_path,
    // sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes, mu_dist,
    // '?'); // recomputing sigma dist would be expensive // DO NOT REMOVE
    fprintf(mcfg->io.f_ptr, "%lu\t%li\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%c", i + 1, mmut->sentence.num_containing_mwe,
            mmut->sentence.num_containing_mwe_tp_only, mmut->sentence.num_all, mmut->document.num_all, mcfg->io.w2v_path,
            sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes,
            mu_dist, '?'); // recomputing sigma dist would be expensive

    /*
    if(mcfg->enable.pairwise){printf("[log] [end iter] pairwise: %f\n", iter_state_pairwise.result); fprintf(mcfg->io.f_ptr,
    "\t%.10e", iter_state_pairwise.result);} if(mcfg->enable.stirling){printf("[log] [end iter] stirling: %f\n",
    iter_state_stirling.result); fprintf(mcfg->io.f_ptr, "\t%.10e", iter_state_stirling.result);}
    if(mcfg->enable.leinster_cobbold_diversity){printf("[log] [end iter] leinster_cobbold: %f, %f\n",
    iter_state_leinster_cobbold.entropy, iter_state_leinster_cobbold.hill_number); fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e",
    iter_state_leinster_cobbold.entropy, iter_state_leinster_cobbold.hill_number);}
    */

    const int32_t log_bfr_size = 128;
    char log_bfr[log_bfr_size];

    if (mcfg->enable.pairwise) {
      memset(log_bfr, '\0', log_bfr_size);
      snprintf(log_bfr, log_bfr_size, "[end iter] pairwise: %f", iter_state_pairwise.result);
      info_format(__FILE__, __func__, __LINE__, log_bfr);
      fprintf(mcfg->io.f_ptr, "\t%.10e", iter_state_pairwise.result);
    }
    if (mcfg->enable.stirling) {
      memset(log_bfr, '\0', log_bfr_size);
      snprintf(log_bfr, log_bfr_size, "[end iter] stirling: %f", iter_state_stirling.result);
      info_format(__FILE__, __func__, __LINE__, log_bfr);
      fprintf(mcfg->io.f_ptr, "\t%.10e", iter_state_stirling.result);
    }
    if (mcfg->enable.leinster_cobbold_diversity) {
      memset(log_bfr, '\0', log_bfr_size);
      snprintf(log_bfr, log_bfr_size, "[end iter] leinster_cobbold: %f, %f", iter_state_leinster_cobbold.entropy,
               iter_state_leinster_cobbold.hill_number);
      info_format(__FILE__, __func__, __LINE__, log_bfr);
      fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e", iter_state_leinster_cobbold.entropy, iter_state_leinster_cobbold.hill_number);
    }

    fprintf(mcfg->io.f_ptr, "\n");
  } else {
    int64_t ns_delta, virtual_mem;

    if (mcfg->io.enable_output_timing) {
      // fprintf(mcfg->io.f_timing_ptr, "%lu\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu", i+1, mmut->sentence.num,
      // mmut->sentence.num_all, mmut->document.num, mcfg->io.w2v_path,
      // sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes); // DO NOT
      // REMOVE
      fprintf(mcfg->io.f_timing_ptr, "%lu\t%li\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu", i + 1, mmut->sentence.num_containing_mwe,
              mmut->sentence.num_containing_mwe_tp_only, mmut->sentence.num_all, mmut->document.num_all, mcfg->io.w2v_path,
              sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes);
      if (time_ns_delta(NULL) != 0) {
        goto time_ns_delta_failure;
      }
    }
    if (mcfg->io.enable_output_memory) {
      // fprintf(mcfg->io.f_memory_ptr, "%lu\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu", i+1, mmut->sentence.num,
      // mmut->sentence.num_all, mmut->document.num, mcfg->io.w2v_path,
      // sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes); // DO NOT
      // REMOVE
      fprintf(mcfg->io.f_memory_ptr, "%lu\t%li\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu", i + 1, mmut->sentence.num_containing_mwe,
              mmut->sentence.num_containing_mwe_tp_only, mmut->sentence.num_all, mmut->document.num_all, mcfg->io.w2v_path,
              sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes);
    }

    struct matrix m_mst = {
        .fp_mode = FP64,
    };
    if (mcfg->enable.functional_evenness) {
      if (create_matrix(&m_mst, sref->g->num_nodes, sref->g->num_nodes, FP64) != 0) {
        perror("failed to call create_matrix for MST\n");
        return 1;
      }
      memset(m_mst.active, '\0', sref->g->num_nodes * sref->g->num_nodes * sizeof(uint8_t));
      memset(m_mst.active_final, '\0', sref->g->num_nodes * sref->g->num_nodes * sizeof(uint8_t));
      if (mcfg->io.enable_output_timing) {
        if (time_ns_delta(&ns_delta) != 0) {
          goto time_ns_delta_failure;
        } else {
          fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
        }
      }
      if (mcfg->io.enable_output_memory) {
        if (virtual_memory_consumption(&virtual_mem) != 0) {
          goto virtual_memory_consumption_failure;
        } else {
          fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
        }
      }
    }

    time_t t, delta_t;

    struct matrix m = {
        .fp_mode = FP32,
    };
    if (enable_distance_computation) {
      if (create_matrix(&m, (uint32_t)sref->g->num_nodes, (uint32_t)sref->g->num_nodes, FP32) != 0) {
        perror("failed to call create_matrix\n");
        return 1;
      }
      t = time(NULL);
      if (mcfg->threading.enable_multithreaded_matrix_generation) {
        if (distance_matrix_from_graph_multithread(sref->g, &m, mcfg->threading.num_matrix_threads) != 0) {
          perror("failed to call distance_matrix_from_graph_multithread\n");
          return 1;
        }
      } else {
        // if(distance_matrix_from_graph(sref->g, ANGULAR_MINKOWSKI_DISTANCE_ORDER, &m) != 0){
        if (distance_matrix_from_graph(sref->g, &m) != 0) {
          perror("failed to call distance_matrix_from_graph\n");
          return 1;
        }
      }
      if (mcfg->io.enable_output_timing) {
        if (time_ns_delta(&ns_delta) != 0) {
          goto time_ns_delta_failure;
        } else {
          fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
        }
      }
      if (mcfg->io.enable_output_memory) {
        if (virtual_memory_consumption(&virtual_mem) != 0) {
          goto virtual_memory_consumption_failure;
        } else {
          fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
        }
      }
      delta_t = time(NULL) - t;
      if (mcfg->io.enable_timings) {
        printf("[log] [time] Computed matrix in %lis\n", delta_t);
      }
    }

    if (mmut->mst_initialised) {
      free_graph_distance_heap(sref->heap);
      if (mcfg->enable.functional_evenness) {
        free_minimum_spanning_tree(sref->mst);
      }
    }

    struct graph_distance_heap local_heap;

    double mu_dist = NAN;
    double sigma_dist = NAN;
    float mu_dist_fp32 = NAN;
    float sigma_dist_fp32 = NAN;
    if (enable_distance_computation && mcfg->enable.functional_evenness) {
      err = create_graph_distance_heap(&local_heap, sref->g, GRAPH_NODE_FP32, &m);
      if (err != 0) {
        perror("failed to call create_graph_distance_sref->heap\n");
        return EXIT_FAILURE;
      }

      (*sref->heap) = local_heap;

      if (mcfg->io.enable_output_timing) {
        if (time_ns_delta(&ns_delta) != 0) {
          goto time_ns_delta_failure;
        } else {
          fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
        }
      }
      if (mcfg->io.enable_output_memory) {
        if (virtual_memory_consumption(&virtual_mem) != 0) {
          goto virtual_memory_consumption_failure;
        } else {
          fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
        }
      }
    }

    if (enable_distance_computation) {
      switch (m.fp_mode) {
      case FP32:
        avg_and_std_fp32(m.bfr.fp32, m.a * m.b, &mu_dist_fp32, &sigma_dist_fp32);
        mu_dist = (double)mu_dist_fp32;
        sigma_dist = (double)sigma_dist_fp32;
        break;
      case FP64:
        avg_and_std_fp64(m.bfr.fp64, m.a * m.b, &mu_dist, &sigma_dist);
        break;
      default:
        perror("unknown FP mode\n");
        return 1;
      }
      if (mcfg->io.enable_output_timing) {
        if (time_ns_delta(&ns_delta) != 0) {
          goto time_ns_delta_failure;
        } else {
          fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
        }
      }
      if (mcfg->io.enable_output_memory) {
        if (virtual_memory_consumption(&virtual_mem) != 0) {
          goto virtual_memory_consumption_failure;
        } else {
          fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
        }
      }
    }

    // fprintf(mcfg->io.f_ptr, "%lu\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%.10e", i+1, mmut->sentence.num,
    // mmut->sentence.num_all, mmut->document.num, mcfg->io.w2v_path,
    // sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes, mu_dist,
    // sigma_dist); // DO NOT REMOVE
    fprintf(mcfg->io.f_ptr, "%lu\t%li\t%li\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%.10e", i + 1,
            mmut->sentence.num_containing_mwe, mmut->sentence.num_containing_mwe_tp_only, mmut->sentence.num_all,
            mmut->document.num_all, mcfg->io.w2v_path,
            sref->sorted_array_discarded_because_not_in_vector_database->num_elements, mmut->best_s, sref->g->num_nodes,
            mu_dist, sigma_dist);

    if (mcfg->enable.disparity_functions) {
      if (mcfg->io.enable_output_timing) {
        if (time_ns_delta(NULL) != 0) {
          goto time_ns_delta_failure;
        }
      }

      if (mcfg->enable.functional_evenness) {
        t = time(NULL);
        struct minimum_spanning_tree local_mst;

        err = create_minimum_spanning_tree(&local_mst, &local_heap);
        if (err != 0) {
          perror("failed to call create_minimum_spanning_tree\n");
          return EXIT_FAILURE;
        }

        err = calculate_minimum_spanning_tree(&local_mst, &m_mst, MST_PRIMS_ALGORITHM);
        if (err != 0) {
          perror("failed to call calculate_minimum_spanning_tree\n");
          return EXIT_FAILURE;
        }

        (*sref->mst) = local_mst;
        mmut->mst_initialised = 1;

        delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed MST in %lis\n", delta_t);
        }
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      mmut->prev_num_nodes = sref->g->num_nodes;

      if (mcfg->enable.stirling) {
        double stirling;
        t = time(NULL);
        err = stirling_from_graph(sref->g, &stirling, mcfg->div_param.stirling_alpha, mcfg->div_param.stirling_beta,
                                  GRAPH_NODE_FP32, &m);
        if (err != 0) {
          perror("failed to call stirling_from_graph\n");
          return EXIT_FAILURE;
        }
        delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Stirling in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", stirling);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }

      if (mcfg->enable.ricotta_szeidl) {
        double ricotta_szeidl;
        t = time(NULL);
        err = ricotta_szeidl_from_graph(sref->g, &ricotta_szeidl, mcfg->div_param.ricotta_szeidl_alpha, GRAPH_NODE_FP32, &m);
        if (err != 0) {
          perror("failed to call ricotta_szeidl_from_graph\n");
          return EXIT_FAILURE;
        }
        delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Ricotta-Szeidl in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", ricotta_szeidl);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }

      if (mcfg->enable.pairwise) {
        if (wrap_diversity_1r_0a(sref->g, &m, GRAPH_NODE_FP32, mcfg->io.f_ptr, mcfg->io.f_timing_ptr, mcfg->io.f_memory_ptr,
                                 pairwise_from_graph, mcfg->io.enable_timings, mcfg->io.enable_output_timing,
                                 mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.chao_et_al_functional_diversity) {
        if (wrap_diversity_2r_1a(sref->g, &m, mcfg->div_param.chao_et_al_functional_diversity_alpha, GRAPH_NODE_FP32,
                                 mcfg->io.f_ptr, mcfg->io.f_timing_ptr, mcfg->io.f_memory_ptr,
                                 chao_et_al_functional_diversity_from_graph, mcfg->io.enable_timings,
                                 mcfg->io.enable_output_timing, mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.scheiner_species_phylogenetic_functional_diversity) {
        if (wrap_diversity_2r_1a(sref->g, &m, mcfg->div_param.scheiner_species_phylogenetic_functional_diversity_alpha,
                                 GRAPH_NODE_FP32, mcfg->io.f_ptr, mcfg->io.f_timing_ptr, mcfg->io.f_memory_ptr,
                                 scheiner_species_phylogenetic_functional_diversity_from_graph, mcfg->io.enable_timings,
                                 mcfg->io.enable_output_timing, mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.leinster_cobbold_diversity) {
        if (wrap_diversity_2r_1a(sref->g, &m, mcfg->div_param.leinster_cobbold_diversity_alpha, GRAPH_NODE_FP32, mcfg->io.f_ptr,
                                 mcfg->io.f_timing_ptr, mcfg->io.f_memory_ptr, leinster_cobbold_diversity_from_graph,
                                 mcfg->io.enable_timings, mcfg->io.enable_output_timing, mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.lexicographic) {
        if (wrap_diversity_2r_0a_long_double_alt(sref->g, &m, GRAPH_NODE_FP32, mcfg->io.f_ptr, mcfg->io.f_timing_ptr,
                                                 mcfg->io.f_memory_ptr, lexicographic_from_graph, mcfg->io.enable_timings,
                                                 mcfg->io.enable_output_timing, mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.functional_evenness) {
        if (!mmut->mst_initialised) {
          perror("sref->mst not initialised\n");
          return 1;
        }
        double functional_evenness;
        t = time(NULL);
        err = functional_evenness_from_minimum_spanning_tree(sref->mst, &functional_evenness);
        if (err != 0) {
          perror("failed to call functional_evenness_from_minimum_spanning_tree\n");
          return EXIT_FAILURE;
        }
        delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed FEve in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", functional_evenness);
        timing_and_memory(mcfg->io.f_timing_ptr, mcfg->io.f_memory_ptr, mcfg->io.enable_output_timing,
                          mcfg->io.enable_output_memory);
      }

      if (mcfg->enable.functional_dispersion) {
        if (wrap_diversity_1r_0a_no_matrix(sref->g, GRAPH_NODE_FP32, mcfg->io.f_ptr, mcfg->io.f_timing_ptr,
                                           mcfg->io.f_memory_ptr, functional_dispersion_from_graph, mcfg->io.enable_timings,
                                           mcfg->io.enable_output_timing, mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }

      if (mcfg->enable.functional_divergence_modified) {
        if (wrap_diversity_1r_0a_no_matrix(sref->g, GRAPH_NODE_FP32, mcfg->io.f_ptr, mcfg->io.f_timing_ptr,
                                           mcfg->io.f_memory_ptr, functional_divergence_modified_from_graph,
                                           mcfg->io.enable_timings, mcfg->io.enable_output_timing,
                                           mcfg->io.enable_output_memory) != 0) {
          return 1;
        }
      }
    }

    if (mcfg->enable.non_disparity_functions) {
      struct cpu_info local_cpu_info = {0};
      if (get_cpu_info(&local_cpu_info) != 0) {
        perror("Failed to call get_cpu_info\n");
        return 1;
      }
      if (mcfg->enable.shannon_weaver_entropy) {
        double res_entropy;
        double res_hill_number;
        time_t t = time(NULL);

#if ENABLE_NON_DISPARITY_MULTITHREADING == 1
        if (non_disparity_multithread(entropy_shannon_weaver_transform_proportion, agregate_local_add, NULL, 0.0, sref->g, 0.0,
                                      0.0, &res_entropy, &res_hill_number,
                                      (int32_t)local_cpu_info.cardinality_virtual_cores) != 0) {
          perror("Failed to call non_disparity_multithread\n");
          return 1;
        }
        res_hill_number = entropy_shannon_weaver_to_hill_number(res_entropy);
#else
        shannon_weaver_entropy_from_graph(sref->g, &res_entropy, &res_hill_number);
#endif

        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed SW entropy in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.good_entropy) {
        double res;
        time_t t = time(NULL);

#if ENABLE_NON_DISPARITY_MULTITHREADING == 1
        if (non_disparity_multithread(entropy_good_transform_proportion, agregate_local_add, NULL, 0.0, sref->g,
                                      mcfg->div_param.good_alpha, mcfg->div_param.good_beta, &res, NULL,
                                      (int32_t)local_cpu_info.cardinality_virtual_cores) != 0) {
          perror("Failed to call non_disparity_multithread\n");
          return 1;
        }
#else
        good_entropy_from_graph(sref->g, &res, mcfg->div_param.good_alpha, mcfg->div_param.good_beta);
#endif

        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Good entropy in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.renyi_entropy) {
        double res_entropy;
        double res_hill_number;
        time_t t = time(NULL);

#if ENABLE_NON_DISPARITY_MULTITHREADING == 1
        if (non_disparity_multithread(entropy_renyi_transform_proportion, agregate_local_add, entropy_renyi_finalise, 0.0,
                                      sref->g, mcfg->div_param.renyi_alpha, 0.0, &res_entropy, &res_hill_number,
                                      (int32_t)local_cpu_info.cardinality_virtual_cores) != 0) {
          perror("Failed to call non_disparity_multithread\n");
          return 1;
        }
        res_hill_number = entropy_renyi_to_hill_number(res_entropy);
#else
        renyi_entropy_from_graph(sref->g, &res_entropy, &res_hill_number, mcfg->div_param.renyi_alpha);
#endif

        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Renyi entropy in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.patil_taillie_entropy) {
        double res_entropy;
        double res_hill_number;
        time_t t = time(NULL);
        patil_taillie_entropy_from_graph(sref->g, &res_entropy, &res_hill_number, mcfg->div_param.patil_taillie_alpha);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Patil-Taillie entropy in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.q_logarithmic_entropy) {
        double res_entropy;
        double res_hill_number;
        time_t t = time(NULL);
        q_logarithmic_entropy_from_graph(sref->g, &res_entropy, &res_hill_number, mcfg->div_param.q_logarithmic_q);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed q-logarithmic entropy in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.simpson_index) {
        double res;
        time_t t = time(NULL);
        simpson_index_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Simpson index in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.simpson_dominance_index) {
        double res;
        time_t t = time(NULL);
        simpson_dominance_index_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Simpson dominance index in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.hill_number_standard) {
        double res;
        time_t t = time(NULL);
        hill_number_standard_from_graph(sref->g, &res, mcfg->div_param.hill_number_standard_alpha);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Hill number (standard) in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.hill_evenness) {
        double res;
        time_t t = time(NULL);
        hill_evenness_from_graph(sref->g, &res, mcfg->div_param.hill_evenness_alpha, mcfg->div_param.hill_evenness_beta);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Hill evenness in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.berger_parker_index) {
        double res;
        time_t t = time(NULL);
        berger_parker_index_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Berger Parker index in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.junge1994_page22) {
        double res;
        time_t t = time(NULL);
        junge1994_page22_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Junge 1994 p22 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.brillouin_diversity) {
        double res;
        time_t t = time(NULL);
        brillouin_diversity_from_graph(sref->g, &res);
        // printf("res: %f\n", res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed Brillouin diversity in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.mcintosh_index) {
        double res;
        time_t t = time(NULL);
        mcintosh_index_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed McIntosh index in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_entropy_over_log_n_species_pielou1975) {
        double res;
        time_t t = time(NULL);
        sw_entropy_over_log_n_species_pielou1975_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) entropy over log n species Pielou 1975 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_heip) {
        double res;
        time_t t = time(NULL);
        sw_e_heip_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E Heip in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_one_minus_d) {
        double res;
        time_t t = time(NULL);
        sw_e_one_minus_D_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E one minus D in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_one_over_ln_d_williams1964) {
        double res;
        time_t t = time(NULL);
        sw_e_one_over_D_williams1964_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E one over ln D Williams 1964 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_minus_ln_d_pielou1977) {
        double res;
        time_t t = time(NULL);
        sw_e_minus_ln_D_pielou1977_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E minus ln D Pielou 1977 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_f_2_1_alatalo1981) {
        double res;
        time_t t = time(NULL);
        sw_f_2_1_alatalo1981_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) F_2_1 Alatalo 1981 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_g_2_1_molinari1989) {
        double res;
        time_t t = time(NULL);
        sw_g_2_1_molinari1989_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) G_2_1 Molinari 1989 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_bulla1994) {
        double res;
        time_t t = time(NULL);
        sw_e_bulla1994_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E Bulla 1994 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_o_bulla1994) {
        double res;
        time_t t = time(NULL);
        sw_o_bulla1994_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) O bulla 1994 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_mci_pielou1969) {
        double res;
        time_t t = time(NULL);
        sw_e_mci_pielou1969_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E MCI Pielou 1969 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_prime_camargo1993) {
        double res;
        time_t t = time(NULL);
        if (mcfg->threading.enable_sw_e_prime_camargo1993_multithreading) {
          sw_e_prime_camargo1993_from_graph_multithread(sref->g, &res, mcfg->threading.num_matrix_threads);
        } else {
          sw_e_prime_camargo1993_from_graph(sref->g, &res);
        }
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E prime Camargo 1993 in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
      if (mcfg->enable.sw_e_var_smith_and_wilson1996_original) {
        double res;
        time_t t = time(NULL);
        sw_e_var_smith_and_wilson1996_original_from_graph(sref->g, &res);
        time_t delta_t = time(NULL) - t;
        if (mcfg->io.enable_timings) {
          printf("[log] [time] Computed (SW) E var Smith and Wilson 1996 original in %lis\n", delta_t);
        }
        fprintf(mcfg->io.f_ptr, "\t%.10e", res);
        if (mcfg->io.enable_output_timing) {
          if (time_ns_delta(&ns_delta) != 0) {
            goto time_ns_delta_failure;
          } else {
            fprintf(mcfg->io.f_timing_ptr, "\t%li", ns_delta);
          }
        }
        if (mcfg->io.enable_output_memory) {
          if (virtual_memory_consumption(&virtual_mem) != 0) {
            goto virtual_memory_consumption_failure;
          } else {
            fprintf(mcfg->io.f_memory_ptr, "\t%li", virtual_mem);
          }
        }
      }
    }

    fprintf(mcfg->io.f_ptr, "\n");
    if (mcfg->io.enable_output_timing) {
      fprintf(mcfg->io.f_timing_ptr, "\n");
    }
    if (mcfg->io.enable_output_memory) {
      fprintf(mcfg->io.f_memory_ptr, "\n");
    }

    if (enable_distance_computation) {
      free_matrix(&m);
      if (mcfg->enable.functional_evenness) {
        free_matrix(&m_mst);
      }
    }
  }

  return 0;

time_ns_delta_failure:
  perror("Failed to call time_ns_delta\n");
  return 1;

virtual_memory_consumption_failure:
  perror("Failed to call virtual_memory_consumption\n");
  return 1;

malloc_fail:
  perror("malloc failed\n");
  return 1;
}
