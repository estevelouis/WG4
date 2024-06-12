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

#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <stdint.h>
#include <stdio.h>

#include "graph.h"
#include "sorted_array/array.h"

struct measurement_diversity_parameters {
  const double stirling_alpha;
  const double stirling_beta;
  const double ricotta_szeidl_alpha;
  const double chao_et_al_functional_diversity_alpha;
  const double scheiner_species_phylogenetic_functional_diversity_alpha;
  const double leinster_cobbold_diversity_alpha;
  const double good_alpha;
  const double good_beta;
  const double renyi_alpha;
  const double patil_taillie_alpha;
  const double q_logarithmic_q;
  const double hill_number_standard_alpha;
  const double hill_evenness_alpha;
  const double hill_evenness_beta;
};

struct measurement_diversity_enabler {
  const uint8_t stirling;
  const uint8_t ricotta_szeidl;
  const uint8_t pairwise;
  const uint8_t lexicographic;
  const uint8_t chao_et_al_functional_diversity;
  const uint8_t scheiner_species_phylogenetic_functional_diversity;
  const uint8_t leinster_cobbold_diversity;
  const uint8_t functional_evenness;
  const uint8_t functional_dispersion;
  const uint8_t functional_divergence_modified;
  const uint8_t non_disparity_functions;
  const uint8_t disparity_functions;
  const uint8_t shannon_weaver_entropy;
  const uint8_t good_entropy;
  const uint8_t renyi_entropy;
  const uint8_t patil_taillie_entropy;
  const uint8_t q_logarithmic_entropy;
  const uint8_t simpson_index;
  const uint8_t simpson_dominance_index;
  const uint8_t hill_number_standard;
  const uint8_t hill_evenness;
  const uint8_t berger_parker_index;
  const uint8_t junge1994_page22;
  const uint8_t brillouin_diversity;
  const uint8_t mcintosh_index;
  const uint8_t sw_entropy_over_log_n_species_pielou1975;
  const uint8_t sw_e_heip;
  const uint8_t sw_e_one_minus_d;
  const uint8_t sw_e_one_over_ln_d_williams1964;
  const uint8_t sw_e_minus_ln_d_pielou1977;
  const uint8_t sw_f_2_1_alatalo1981;
  const uint8_t sw_g_2_1_molinari1989;
  const uint8_t sw_e_bulla1994;
  const uint8_t sw_o_bulla1994;
  const uint8_t sw_e_mci_pielou1969;
  const uint8_t sw_e_prime_camargo1993;
  const uint8_t sw_e_var_smith_and_wilson1996_original;
};

struct measurement_io {
  const char *const w2v_path;
  const char *const jsonl_content_key;
  const char *const input_path;
  const char *const output_path;
  const char *const output_path_timing;
  const char *const output_path_memory;
  const char *const udpipe_model_path;
  FILE *f_ptr;
  FILE *f_timing_ptr;
  FILE *f_memory_ptr;
  const uint8_t enable_timings;
  const uint8_t enable_output_timing;
  const uint8_t enable_output_memory;
  const uint8_t force_timing_and_memory_to_output_path;
};

struct measurement_threading {
  const int32_t num_row_threads;
  const int32_t num_matrix_threads;
  const int32_t num_file_reading_threads;
  const uint8_t enable_multithreaded_matrix_generation;
  const uint8_t enable_iterative_distance_computation;
  const uint8_t enable_multithreaded_row_generation;
  const int8_t row_generation_batch_size;
};

struct measurement_step {
  const uint8_t enable_count_recompute_step;
  const uint8_t use_log10;
  const uint64_t recompute_step;     // const?
  const double recompute_step_log10; // const?
};

struct measurement_step_parameters {
  struct measurement_step sentence;
  struct measurement_step document;
};

// !
struct measurement_configuration {
  const uint32_t target_column;
  const char *const jsonl_content_key;
  const struct measurement_diversity_parameters div_param;
  const struct measurement_diversity_enabler enable;
  struct measurement_io io;
  const struct measurement_threading threading;
  const struct measurement_step_parameters steps; // const?
};

// !
struct measurement_structure_references {
  struct graph *const g;
  struct minimum_spanning_tree *const mst;
  struct graph_distance_heap *const heap;
  struct word2vec *const w2v;
  struct sorted_array *const sorted_array_discarded_because_not_in_vector_database;
};

struct measurement_mutable_counters {
  uint64_t num;
  uint64_t num_all;
  uint64_t count_target;
  double stacked_log;
};

// !
struct measurement_mutables {
  double best_s;
  double prev_best_s;
  int64_t prev_num_nodes;
  struct measurement_mutable_counters sentence;
  struct measurement_mutable_counters document;
  uint8_t mst_initialised;
  pthread_mutex_t mutex;
};

// !
struct measurement_file_thread {
  const int32_t i;
  const char *const filename;
  struct measurement_configuration *const mcfg;
  struct measurement_structure_references *const sref;
  struct measurement_mutables *const mmut;
};

int32_t time_ns_delta(int64_t *const delta);

int32_t virtual_memory_consumption(int64_t *const res);

void timing_and_memory(FILE *f_timing_ptr, FILE *f_memory_ptr, const uint8_t enable_output_timing,
                       const uint8_t enable_output_memory);

int32_t wrap_diversity_1r_0a(struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr,
                             FILE *f_timing_ptr, FILE *f_memory_ptr,
                             int32_t (*df)(struct graph *const, double *const, const int8_t, const struct matrix *const),
                             const uint8_t enable_timings, const uint8_t enable_output_timing,
                             const uint8_t enable_output_memory);
int32_t wrap_diversity_1r_0a_no_matrix(struct graph *const g, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr,
                                       FILE *f_memory_ptr, int32_t (*df)(struct graph *const, double *const, const int8_t),
                                       const uint8_t enable_timings, const uint8_t enable_output_timing,
                                       const uint8_t enable_output_memory);
int32_t wrap_diversity_2r_1a(struct graph *const g, struct matrix *const m, const double alpha, const int8_t fp_mode,
                             FILE *f_ptr, FILE *f_timing_ptr, FILE *f_memory_ptr,
                             int32_t (*df)(struct graph *const, double *const, double *const, double, const int8_t,
                                           const struct matrix *const),
                             const uint8_t enable_timings, const uint8_t enable_output_timing,
                             const uint8_t enable_output_memory);
int32_t
wrap_diversity_2r_0a(struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr,
                     FILE *f_memory_ptr,
                     int32_t (*df)(struct graph *const, double *const, double *const, const int8_t, const struct matrix *const),
                     const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory);
int32_t wrap_diversity_2r_0a_long_double_alt(
    struct graph *const g, struct matrix *const m, const int8_t fp_mode, FILE *f_ptr, FILE *f_timing_ptr, FILE *f_memory_ptr,
    int32_t (*df)(struct graph *const, double *const, long double *const, const int8_t, const struct matrix *const),
    const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory);
/* int32_t wrap_diversity_1r_0a(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE*
 * f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double*, const int8_t, struct matrix* const), const
 * uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory); */

/*
int32_t apply_diversity_functions_to_graph(struct graph* g, struct minimum_spanning_tree* mst, struct graph_distance_heap* heap,
FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int64_t* previous_g_num_nodes_p, int64_t* num_sentences_p, int64_t*
num_all_sentences_p, int64_t* num_documents_p, double* best_s, int8_t* mst_initialised, uint64_t i, struct sorted_array*
sorted_array_discarded_because_not_in_vector_database, char* w2v_path, uint32_t num_row_threads, uint32_t num_matrix_threads,
        uint8_t enable_stirling,
        uint8_t enable_ricotta_szeidl,
        uint8_t enable_pairwise,
        uint8_t enable_lexicographic,
        uint8_t enable_chao_et_al_functional_diversity,
        uint8_t enable_scheiner_species_phylogenetic_functional_diversity,
        uint8_t enable_leinster_cobbold_diversity,
        uint8_t enable_multithreaded_matrix_generation,
        uint8_t enable_timings,
        uint8_t enable_iterative_distance_computation,
        uint8_t enable_multithreaded_row_generation,
        uint8_t row_generation_batch_size,
        uint8_t enable_output_timing,
        uint8_t enable_output_memory,
        uint8_t enable_functional_evenness,
        uint8_t enable_functional_dispersion,
        uint8_t enable_functional_divergence_modified,
        uint8_t enable_non_disparity_functions,
        uint8_t enable_disparity_functions,
        uint8_t enable_shannon_weaver_entropy,
        uint8_t enable_good_entropy,
        uint8_t enable_renyi_entropy,
        uint8_t enable_patil_taillie_entropy,
        uint8_t enable_q_logarithmic_entropy,
        uint8_t enable_simpson_index,
        uint8_t enable_simpson_dominance_index,
        uint8_t enable_hill_number_standard,
        uint8_t enable_hill_evenness,
        uint8_t enable_berger_parker_index,
        uint8_t enable_junge1994_page22,
        uint8_t enable_brillouin_diversity,
        uint8_t enable_mcintosh_index,
        uint8_t enable_sw_entropy_over_log_n_species_pielou1975,
        uint8_t enable_sw_e_heip,
        uint8_t enable_sw_e_one_minus_d,
        uint8_t enable_sw_e_one_over_ln_d_williams1964,
        uint8_t enable_sw_e_minus_ln_d_pielou1977,
        uint8_t enable_sw_f_2_1_alatalo1981,
        uint8_t enable_sw_g_2_1_molinari1989,
        uint8_t enable_sw_e_bulla1994,
        uint8_t enable_sw_o_bulla1994,
        uint8_t enable_sw_e_mci_pielou1969,
        uint8_t enable_sw_e_prime_camargo1993,
        uint8_t enable_sw_e_var_smith_and_wilson1996_original,
        double stirling_alpha,
        double stirling_beta,
        double ricotta_szeidl_alpha,
        double chao_et_al_functional_diversity_alpha,
        double scheiner_species_phylogenetic_functional_diversity_alpha,
        double leinster_cobbold_diversity_alpha,
        double good_alpha,
        double good_beta,
        double renyi_alpha,
        double patil_taillie_alpha,
        double q_logarithmic_q,
        double hill_number_standard_alpha,
        double hill_evenness_alpha,
        double hill_evenness_beta
);
*/
// int32_t apply_diversity_functions_to_graph(struct measurement_configuration * const mcfg, struct
// measurement_structure_references * const sref, struct measurement_mutables * const mmut);
int32_t apply_diversity_functions_to_graph(const uint64_t i, struct measurement_configuration *const mcfg,
                                           struct measurement_structure_references *const sref,
                                           struct measurement_mutables *const mmut);

// int32_t measurement(struct measurement_configuration * const mcfg);

#endif
