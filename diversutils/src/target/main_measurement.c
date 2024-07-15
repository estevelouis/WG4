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

#include <sys/stat.h>
#include <sys/types.h> // ?
#include <unistd.h>    // ?

#include "cpu.h"

#include "dfunctions.h"
#include "distributions.h"
#include "graph.h"
#include "logging.h"
#include "measurement.h"
#include "sorted_array/array.h"
#include "stats.h"

#include "cupt/load.h"
#include "cupt/parser.h"
#include "jsonl/load.h"
#include "jsonl/parser.h"

#include "udpipe/interface/cinterface.h"

#include "macroconfig.h"

#include "filter.h"

const uint8_t ENABLE_DISTANCE_COMPUTATION =
    ENABLE_DISPARITY_FUNCTIONS &&
    (ENABLE_STIRLING || ENABLE_RICOTTA_SZEIDL || ENABLE_PAIRWISE || ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY ||
     ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY || ENABLE_LEINSTER_COBBOLD_DIVERSITY || ENABLE_LEXICOGRAPHIC ||
     ENABLE_FUNCTIONAL_EVENNESS || ENABLE_FUNCTIONAL_DISPERSION || ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED);

enum { CONLLU_COLUMN_UPOS, CONLLU_COLUMN_DEPREL };
const int32_t CONLLU_COLUMNS_TO_ADD[2] = {CONLLU_COLUMN_UPOS, CONLLU_COLUMN_DEPREL};
// const int32_t NUM_CONLLU_COLUMNS_TO_ADD = (int32_t) (sizeof(CONLLU_COLUMNS_TO_ADD) / sizeof(int32_t));
const int32_t NUM_CONLLU_COLUMNS_TO_ADD = 0;

const int32_t CONLLU_ADD_FORM = 0;

double stacked_sentence_count_log10 = SENTENCE_COUNT_RECOMPUTE_STEP_LOG10;
int64_t stacked_sentence_count_target;
double stacked_document_count_log10 = DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10;
int64_t stacked_document_count_target;

int32_t read_list_simulation_files(const char *const path_list_simulation_files, char ***const resulting_paths,
                                   int32_t *num_input_paths) {
  FILE *f_paths_ptr = fopen(path_list_simulation_files, "r");
  if (f_paths_ptr == NULL) {
    fprintf(stderr, "cannot open %s\n", path_list_simulation_files);
    return 1;
  }
  // char* bfr[MAX_FILES];
  char **bfr;
  const size_t bfr_step = 32;
  size_t bfr_capacity = bfr_step;
  size_t bfr_alloc_size = bfr_capacity * sizeof(char *);
  bfr = malloc(bfr_alloc_size);
  if (bfr == NULL) {
    goto failure_alloc;
  }
  memset(bfr, '\0', bfr_alloc_size);

  char bfr_read[BFR_SIZE];
  memset(bfr_read, '\0', BFR_SIZE);

  int32_t num_files = 0;
  while (fgets(bfr_read, BFR_SIZE - 1, f_paths_ptr)) {
    if (bfr_read[0] == '#') {
      continue;
    }

    size_t bytes_to_cpy = strlen(bfr_read);
    if (bfr_read[bytes_to_cpy - 1] == '\n') {
      bytes_to_cpy--;
    }
    void *local_malloc_p = malloc(bytes_to_cpy + 1);
    if (local_malloc_p == NULL) {
      perror("malloc failed\n");
      return 1;
    }
    memset(local_malloc_p, '\0', bytes_to_cpy + 1);
    memcpy(local_malloc_p, bfr_read, bytes_to_cpy);
    bfr[num_files] = (char *)local_malloc_p;

    num_files++;

    if (((size_t)num_files) >= bfr_capacity) {
      size_t new_capacity = bfr_capacity + bfr_step;
      bfr_alloc_size = new_capacity * sizeof(char *);
      bfr = realloc(bfr, bfr_alloc_size);
      if (bfr == NULL) {
        goto failure_alloc;
      }
      memset(&bfr[bfr_capacity], '\0', bfr_step * sizeof(char *));
      bfr_capacity = new_capacity;
    }
  }

  fclose(f_paths_ptr);

  *resulting_paths = bfr;
  *num_input_paths = num_files;

  return 0;

failure_alloc:
  perror("alloc failed\n");
  return 1;
}

int32_t measurement(struct measurement_configuration *const mcfg) {
  const int32_t log_bfr_size = 512;
  char log_bfr[512];

  const uint8_t enable_distance_computation =
      mcfg->enable.disparity_functions &&
      (mcfg->enable.stirling || mcfg->enable.ricotta_szeidl || mcfg->enable.pairwise ||
       mcfg->enable.chao_et_al_functional_diversity || mcfg->enable.scheiner_species_phylogenetic_functional_diversity ||
       mcfg->enable.leinster_cobbold_diversity || mcfg->enable.lexicographic || mcfg->enable.functional_evenness ||
       mcfg->enable.functional_dispersion || mcfg->enable.functional_divergence_modified);

  stacked_sentence_count_target = log(stacked_sentence_count_log10) / log(10.0);
  stacked_document_count_target = log(stacked_document_count_log10) / log(10.0);

  struct sorted_array sorted_array_discarded_because_not_in_vector_database;
  if (create_sorted_array(&sorted_array_discarded_because_not_in_vector_database, 0,
                          sizeof(struct sorted_array_str_int_element), sorted_array_str_int_cmp) != 0) {
    perror("failed to call create_sorted_array\n");
    return 1;
  }

  struct graph g = {0};
  if (create_graph_empty(&g) != 0) {
    perror("failed to call create_graph_empty\n");
    return 1;
  }
  struct graph_distance_heap heap;
  struct minimum_spanning_tree mst;

  memset(&heap, '\0', sizeof(struct graph_distance_heap));
  memset(&mst, '\0', sizeof(struct minimum_spanning_tree));

  int32_t err;

  struct word2vec w2v;
  err = load_word2vec_binary(&w2v, mcfg->io.w2v_path);
  if (err != 0) {
    fprintf(stderr, "failed to call load_word2vec binary: %s\n", mcfg->io.w2v_path);
    return 1;
  }

  for (uint64_t i = 0; i < w2v.num_vectors; i++) {
    w2v.keys[i].active_in_current_graph = 0;
    w2v.keys[i].num_occurrences = 0;
    w2v.keys[i].graph_node_pointer = NULL; // !
  }

  int32_t num_input_paths = 0;
  int32_t num_input_paths_true_positive = 0;
  char **input_paths = NULL;
  char **input_paths_tp = NULL; // !

  // <--
  if (read_list_simulation_files(mcfg->io.input_path, &input_paths, &num_input_paths) != 0) {
    goto failure_read_list_simulation_files;
  }
  if (mcfg->io.input_path_tp != NULL) {
    if (read_list_simulation_files(mcfg->io.input_path_tp, &input_paths_tp, &num_input_paths_true_positive) != 0) {
      goto failure_read_list_simulation_files;
    }
  }

  mcfg->io.f_ptr = fopen(mcfg->io.output_path, "w");
  if (mcfg->io.f_ptr == NULL) {
    fprintf(stderr, "Failed to open file: %s\n", mcfg->io.output_path);
    return EXIT_FAILURE;
  }
  // fprintf(mcfg->io.f_ptr,
  // "num_active_files\tnum_active_sentences\tnum_all_sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
  // // DO NOT REMOVE
  fprintf(mcfg->io.f_ptr, "num_active_files\tnum_sentences_containing_mwe\tnum_sentences_containing_mwe_tp_only\tnum_all_"
                          "sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");

  if (mcfg->enable.disparity_functions) {
    if (mcfg->enable.stirling) {
      fprintf(mcfg->io.f_ptr, "\tstirling_alpha%.10e_beta%.10e", mcfg->div_param.stirling_alpha, mcfg->div_param.stirling_beta);
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.ricotta_szeidl) {
      fprintf(mcfg->io.f_ptr, "\tricotta_szeidl_alpha%.10e", mcfg->div_param.ricotta_szeidl_alpha);
    }
    if (mcfg->enable.pairwise) {
      fprintf(mcfg->io.f_ptr, "\tpairwise");
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.chao_et_al_functional_diversity) {
      fprintf(mcfg->io.f_ptr, "\tchao_et_al_functional_diversity_alpha%.10e\tchao_et_al_functional_hill_number_alpha%.10e",
              mcfg->div_param.chao_et_al_functional_diversity_alpha, mcfg->div_param.chao_et_al_functional_diversity_alpha);
    }
    if (!mcfg->threading.enable_iterative_distance_computation &&
        mcfg->enable.scheiner_species_phylogenetic_functional_diversity) {
      fprintf(mcfg->io.f_ptr,
              "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e\tscheiner_species_phylogenetic_functional_hill_"
              "number_alpha%.10e",
              mcfg->div_param.scheiner_species_phylogenetic_functional_diversity_alpha,
              mcfg->div_param.scheiner_species_phylogenetic_functional_diversity_alpha);
    }
    if (mcfg->enable.leinster_cobbold_diversity) {
      fprintf(mcfg->io.f_ptr, "\tleinster_cobbold_diversity_alpha%.10e\tleinster_cobbold_hill_number_alpha%.10e",
              mcfg->div_param.leinster_cobbold_diversity_alpha, mcfg->div_param.leinster_cobbold_diversity_alpha);
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.lexicographic) {
      fprintf(mcfg->io.f_ptr, "\tlexicographic\tlexicographic_hybrid_scheiner");
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_evenness) {
      fprintf(mcfg->io.f_ptr, "\tfunctional_evenness");
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_dispersion) {
      fprintf(mcfg->io.f_ptr, "\tfunctional_dispersion");
    }
    if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_divergence_modified) {
      fprintf(mcfg->io.f_ptr, "\tfunctional_divergence_modified");
    }
  }
  if (mcfg->enable.non_disparity_functions) {
    if (mcfg->enable.shannon_weaver_entropy) {
      fprintf(mcfg->io.f_ptr, "\tshannon_weaver_entropy\tshannon_weaver_hill_number");
    }
    if (mcfg->enable.good_entropy) {
      fprintf(mcfg->io.f_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", mcfg->div_param.good_alpha, mcfg->div_param.good_beta);
    }
    if (mcfg->enable.renyi_entropy) {
      fprintf(mcfg->io.f_ptr, "\trenyi_entropy_alpha%.4e\trenyi_hill_number_alpha%.4e", mcfg->div_param.renyi_alpha,
              mcfg->div_param.renyi_alpha);
    }
    if (mcfg->enable.patil_taillie_entropy) {
      fprintf(mcfg->io.f_ptr, "\tpatil_taillie_entropy_alpha%.4e\tpatil_taillie_hill_number_alpha%.4e",
              mcfg->div_param.patil_taillie_alpha, mcfg->div_param.patil_taillie_alpha);
    }
    if (mcfg->enable.q_logarithmic_entropy) {
      fprintf(mcfg->io.f_ptr, "\tq_logarithmic_entropy_alpha%.4e\tq_logarithmic_hill_number_alpha%.4e",
              mcfg->div_param.q_logarithmic_q, mcfg->div_param.q_logarithmic_q);
    }
    if (mcfg->enable.simpson_index) {
      fprintf(mcfg->io.f_ptr, "\tsimpson_index");
    }
    if (mcfg->enable.simpson_dominance_index) {
      fprintf(mcfg->io.f_ptr, "\tsimpson_dominance_index");
    }
    if (mcfg->enable.hill_number_standard) {
      fprintf(mcfg->io.f_ptr, "\thill_number_standard_alpha%.4e", mcfg->div_param.hill_number_standard_alpha);
    }
    if (mcfg->enable.hill_evenness) {
      fprintf(mcfg->io.f_ptr, "\thill_evenness_alpha%.4e_beta%.4e", mcfg->div_param.hill_evenness_alpha,
              mcfg->div_param.hill_evenness_beta);
    }
    if (mcfg->enable.berger_parker_index) {
      fprintf(mcfg->io.f_ptr, "\tberger_parker_index");
    }
    if (mcfg->enable.junge1994_page22) {
      fprintf(mcfg->io.f_ptr, "\tjunge1994_page22");
    }
    if (mcfg->enable.brillouin_diversity) {
      fprintf(mcfg->io.f_ptr, "\tbrillouin_diversity");
    }
    if (mcfg->enable.mcintosh_index) {
      fprintf(mcfg->io.f_ptr, "\tmcintosh_index");
    }
    if (mcfg->enable.sw_entropy_over_log_n_species_pielou1975) {
      fprintf(mcfg->io.f_ptr, "\tsw_entropy_over_log_n_species_pielou1975");
    }
    if (mcfg->enable.sw_e_heip) {
      fprintf(mcfg->io.f_ptr, "\tsw_e_heip");
    }
    if (mcfg->enable.sw_e_one_minus_d) {
      fprintf(mcfg->io.f_ptr, "\te_one_minus_D");
    }
    if (mcfg->enable.sw_e_one_over_ln_d_williams1964) {
      fprintf(mcfg->io.f_ptr, "\te_one_over_ln_D_williams1964");
    }
    if (mcfg->enable.sw_e_minus_ln_d_pielou1977) {
      fprintf(mcfg->io.f_ptr, "\te_minus_ln_D_pielou1977");
    }
    if (mcfg->enable.sw_f_2_1_alatalo1981) {
      fprintf(mcfg->io.f_ptr, "\tsw_f_2_1_alatalo1981");
    }
    if (mcfg->enable.sw_g_2_1_molinari1989) {
      fprintf(mcfg->io.f_ptr, "\tsw_g_2_1_molinari1989");
    }
    if (mcfg->enable.sw_e_bulla1994) {
      fprintf(mcfg->io.f_ptr, "\tsw_e_bulla1994");
    }
    if (mcfg->enable.sw_o_bulla1994) {
      fprintf(mcfg->io.f_ptr, "\tsw_o_bulla1994");
    }
    if (mcfg->enable.sw_e_mci_pielou1969) {
      fprintf(mcfg->io.f_ptr, "\tsw_e_mci_pielou1969");
    }
    if (mcfg->enable.sw_e_prime_camargo1993) {
      fprintf(mcfg->io.f_ptr, "\tsw_e_prime_camargo1993");
    }
    if (mcfg->enable.sw_e_var_smith_and_wilson1996_original) {
      fprintf(mcfg->io.f_ptr, "\tsw_e_var_smith_and_wilson1996_original");
    }
  }
  fprintf(mcfg->io.f_ptr, "\n");

  mcfg->io.f_timing_ptr = NULL;
  if (mcfg->io.enable_output_timing) {
    mcfg->io.f_timing_ptr = fopen(mcfg->io.output_path_timing, "w");
    if (mcfg->io.f_timing_ptr == NULL) {
      fprintf(stderr, "Failed to open file: %s\n", mcfg->io.output_path_timing);
      return EXIT_FAILURE;
    }
    // fprintf(mcfg->io.f_timing_ptr,
    // "num_active_files\tnum_active_sentences\tnum_all_sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
    fprintf(mcfg->io.f_timing_ptr,
            "num_active_files\tnum_active_sentences\tnum_all_sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn");

    if (mcfg->enable.disparity_functions) {
      // ----
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_timing_ptr, "\tm_mst_creation");
      }
      if (enable_distance_computation) {
        fprintf(mcfg->io.f_timing_ptr, "\tdist_matrix_computation");
      }
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_timing_ptr, "\tdist_heap_computation");
      }
      if (enable_distance_computation) {
        fprintf(mcfg->io.f_timing_ptr, "\tdist_stat_computation");
      }
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_timing_ptr, "\tmst_computation");
      }
      // ----

      if (mcfg->enable.stirling) {
        fprintf(mcfg->io.f_timing_ptr, "\tstirling_alpha%.10e_beta%.10e", mcfg->div_param.stirling_alpha,
                mcfg->div_param.stirling_beta);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.ricotta_szeidl) {
        fprintf(mcfg->io.f_timing_ptr, "\tricotta_szeidl_alpha%.10e", mcfg->div_param.ricotta_szeidl_alpha);
      }
      if (mcfg->enable.pairwise) {
        fprintf(mcfg->io.f_timing_ptr, "\tpairwise");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.chao_et_al_functional_diversity) {
        fprintf(mcfg->io.f_timing_ptr, "\tchao_et_al_functional_diversity_alpha%.10e",
                mcfg->div_param.chao_et_al_functional_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation &&
          mcfg->enable.scheiner_species_phylogenetic_functional_diversity) {
        fprintf(mcfg->io.f_timing_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e",
                mcfg->div_param.scheiner_species_phylogenetic_functional_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.leinster_cobbold_diversity) {
        fprintf(mcfg->io.f_timing_ptr, "\tleinster_cobbold_diversity_alpha%.10e",
                mcfg->div_param.leinster_cobbold_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.lexicographic) {
        fprintf(mcfg->io.f_timing_ptr, "\tlexicographic");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_timing_ptr, "\tfunctional_evenness");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_dispersion) {
        fprintf(mcfg->io.f_timing_ptr, "\tfunctional_dispersion");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_divergence_modified) {
        fprintf(mcfg->io.f_timing_ptr, "\tfunctional_divergence_modified");
      }
    }
    if (mcfg->enable.non_disparity_functions) {
      if (mcfg->enable.shannon_weaver_entropy) {
        fprintf(mcfg->io.f_timing_ptr, "\tshannon_weaver_entropy");
      }
      if (mcfg->enable.good_entropy) {
        fprintf(mcfg->io.f_timing_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", mcfg->div_param.good_alpha,
                mcfg->div_param.good_beta);
      }
      if (mcfg->enable.renyi_entropy) {
        fprintf(mcfg->io.f_timing_ptr, "\trenyi_entropy_alpha%.4e", mcfg->div_param.renyi_alpha);
      }
      if (mcfg->enable.patil_taillie_entropy) {
        fprintf(mcfg->io.f_timing_ptr, "\tpatil_taillie_entropy_alpha%.4e", mcfg->div_param.patil_taillie_alpha);
      }
      if (mcfg->enable.q_logarithmic_entropy) {
        fprintf(mcfg->io.f_timing_ptr, "\tq_logarithmic_entropy_alpha%.4e", mcfg->div_param.q_logarithmic_q);
      }
      if (mcfg->enable.simpson_index) {
        fprintf(mcfg->io.f_timing_ptr, "\tsimpson_index");
      }
      if (mcfg->enable.simpson_dominance_index) {
        fprintf(mcfg->io.f_timing_ptr, "\tsimpson_dominance_index");
      }
      if (mcfg->enable.hill_number_standard) {
        fprintf(mcfg->io.f_timing_ptr, "\thill_number_standard_alpha%.4e", mcfg->div_param.hill_number_standard_alpha);
      }
      if (mcfg->enable.hill_evenness) {
        fprintf(mcfg->io.f_timing_ptr, "\thill_evenness_alpha%.4e_beta%.4e", mcfg->div_param.hill_evenness_alpha,
                mcfg->div_param.hill_evenness_beta);
      }
      if (mcfg->enable.berger_parker_index) {
        fprintf(mcfg->io.f_timing_ptr, "\tberger_parker_index");
      }
      if (mcfg->enable.junge1994_page22) {
        fprintf(mcfg->io.f_timing_ptr, "\tjunge1994_page22");
      }
      if (mcfg->enable.brillouin_diversity) {
        fprintf(mcfg->io.f_timing_ptr, "\tbrillouin_diversity");
      }
      if (mcfg->enable.mcintosh_index) {
        fprintf(mcfg->io.f_timing_ptr, "\tmcintosh_index");
      }
      if (mcfg->enable.sw_entropy_over_log_n_species_pielou1975) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_entropy_over_log_n_species_pielou1975");
      }
      if (mcfg->enable.sw_e_heip) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_e_heip");
      }
      if (mcfg->enable.sw_e_one_minus_d) {
        fprintf(mcfg->io.f_timing_ptr, "\te_one_minus_D");
      }
      if (mcfg->enable.sw_e_one_over_ln_d_williams1964) {
        fprintf(mcfg->io.f_timing_ptr, "\te_one_over_ln_D_williams1964");
      }
      if (mcfg->enable.sw_e_minus_ln_d_pielou1977) {
        fprintf(mcfg->io.f_timing_ptr, "\te_minus_ln_D_pielou1977");
      }
      if (mcfg->enable.sw_f_2_1_alatalo1981) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_f_2_1_alatalo1981");
      }
      if (mcfg->enable.sw_g_2_1_molinari1989) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_g_2_1_molinari1989");
      }
      if (mcfg->enable.sw_e_bulla1994) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_e_bulla1994");
      }
      if (mcfg->enable.sw_o_bulla1994) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_o_bulla1994");
      }
      if (mcfg->enable.sw_e_mci_pielou1969) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_e_mci_pielou1969");
      }
      if (mcfg->enable.sw_e_prime_camargo1993) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_e_prime_camargo1993");
      }
      if (mcfg->enable.sw_e_var_smith_and_wilson1996_original) {
        fprintf(mcfg->io.f_timing_ptr, "\tsw_e_var_smith_and_wilson1996_original");
      }
    }
    fprintf(mcfg->io.f_timing_ptr, "\n");
  }

  mcfg->io.f_memory_ptr = NULL;
  if (mcfg->io.enable_output_memory) {
    mcfg->io.f_memory_ptr = fopen(mcfg->io.output_path_memory, "w");
    if (mcfg->io.f_memory_ptr == NULL) {
      fprintf(stderr, "Failed to open file: %s\n", mcfg->io.output_path_memory);
      return EXIT_FAILURE;
    }
    // fprintf(mcfg->io.f_memory_ptr,
    // "num_active_files\tnum_active_sentences\tnum_all_sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
    fprintf(mcfg->io.f_memory_ptr,
            "num_active_files\tnum_active_sentences\tnum_all_sentences\tnum_documents\tw2v\tnum_discarded_types\ts\tn");

    if (mcfg->enable.disparity_functions) {
      // ----
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_memory_ptr, "\tm_mst_creation");
      }
      if (enable_distance_computation) {
        fprintf(mcfg->io.f_memory_ptr, "\tdist_matrix_computation");
      }
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_memory_ptr, "\tdist_heap_computation");
      }
      if (enable_distance_computation) {
        fprintf(mcfg->io.f_memory_ptr, "\tdist_stat_computation");
      }
      if (mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_memory_ptr, "\tmst_computation");
      }
      // ----

      if (mcfg->enable.stirling) {
        fprintf(mcfg->io.f_memory_ptr, "\tstirling_alpha%.10e_beta%.10e", mcfg->div_param.stirling_alpha,
                mcfg->div_param.stirling_beta);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.ricotta_szeidl) {
        fprintf(mcfg->io.f_memory_ptr, "\tricotta_szeidl_alpha%.10e", mcfg->div_param.ricotta_szeidl_alpha);
      }
      if (mcfg->enable.pairwise) {
        fprintf(mcfg->io.f_memory_ptr, "\tpairwise");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.chao_et_al_functional_diversity) {
        fprintf(mcfg->io.f_memory_ptr, "\tchao_et_al_functional_diversity_alpha%.10e",
                mcfg->div_param.chao_et_al_functional_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation &&
          mcfg->enable.scheiner_species_phylogenetic_functional_diversity) {
        fprintf(mcfg->io.f_memory_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e",
                mcfg->div_param.scheiner_species_phylogenetic_functional_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.leinster_cobbold_diversity) {
        fprintf(mcfg->io.f_memory_ptr, "\tleinster_cobbold_diversity_alpha%.10e",
                mcfg->div_param.leinster_cobbold_diversity_alpha);
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.lexicographic) {
        fprintf(mcfg->io.f_memory_ptr, "\tlexicographic");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_evenness) {
        fprintf(mcfg->io.f_memory_ptr, "\tfunctional_evenness");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_dispersion) {
        fprintf(mcfg->io.f_memory_ptr, "\tfunctional_dispersion");
      }
      if (!mcfg->threading.enable_iterative_distance_computation && mcfg->enable.functional_divergence_modified) {
        fprintf(mcfg->io.f_memory_ptr, "\tfunctional_divergence_modified");
      }
    }
    if (mcfg->enable.non_disparity_functions) {
      if (mcfg->enable.shannon_weaver_entropy) {
        fprintf(mcfg->io.f_memory_ptr, "\tshannon_weaver_entropy");
      }
      if (mcfg->enable.good_entropy) {
        fprintf(mcfg->io.f_memory_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", mcfg->div_param.good_alpha,
                mcfg->div_param.good_beta);
      }
      if (mcfg->enable.renyi_entropy) {
        fprintf(mcfg->io.f_memory_ptr, "\trenyi_entropy_alpha%.4e", mcfg->div_param.renyi_alpha);
      }
      if (mcfg->enable.patil_taillie_entropy) {
        fprintf(mcfg->io.f_memory_ptr, "\tpatil_taillie_entropy_alpha%.4e", mcfg->div_param.patil_taillie_alpha);
      }
      if (mcfg->enable.q_logarithmic_entropy) {
        fprintf(mcfg->io.f_memory_ptr, "\tq_logarithmic_entropy_alpha%.4e", mcfg->div_param.q_logarithmic_q);
      }
      if (mcfg->enable.simpson_index) {
        fprintf(mcfg->io.f_memory_ptr, "\tsimpson_index");
      }
      if (mcfg->enable.simpson_dominance_index) {
        fprintf(mcfg->io.f_memory_ptr, "\tsimpson_dominance_index");
      }
      if (mcfg->enable.hill_number_standard) {
        fprintf(mcfg->io.f_memory_ptr, "\thill_number_standard_alpha%.4e", mcfg->div_param.hill_number_standard_alpha);
      }
      if (mcfg->enable.hill_evenness) {
        fprintf(mcfg->io.f_memory_ptr, "\thill_evenness_alpha%.4e_beta%.4e", mcfg->div_param.hill_evenness_alpha,
                mcfg->div_param.hill_evenness_beta);
      }
      if (mcfg->enable.berger_parker_index) {
        fprintf(mcfg->io.f_memory_ptr, "\tberger_parker_index");
      }
      if (mcfg->enable.junge1994_page22) {
        fprintf(mcfg->io.f_memory_ptr, "\tjunge1994_page22");
      }
      if (mcfg->enable.brillouin_diversity) {
        fprintf(mcfg->io.f_memory_ptr, "\tbrillouin_diversity");
      }
      if (mcfg->enable.mcintosh_index) {
        fprintf(mcfg->io.f_memory_ptr, "\tmcintosh_index");
      }
      if (mcfg->enable.sw_entropy_over_log_n_species_pielou1975) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_entropy_over_log_n_species_pielou1975");
      }
      if (mcfg->enable.sw_e_heip) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_e_heip");
      }
      if (mcfg->enable.sw_e_one_minus_d) {
        fprintf(mcfg->io.f_memory_ptr, "\te_one_minus_D");
      }
      if (mcfg->enable.sw_e_one_over_ln_d_williams1964) {
        fprintf(mcfg->io.f_memory_ptr, "\te_one_over_ln_D_williams1964");
      }
      if (mcfg->enable.sw_e_minus_ln_d_pielou1977) {
        fprintf(mcfg->io.f_memory_ptr, "\te_minus_ln_D_pielou1977");
      }
      if (mcfg->enable.sw_f_2_1_alatalo1981) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_f_2_1_alatalo1981");
      }
      if (mcfg->enable.sw_g_2_1_molinari1989) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_g_2_1_molinari1989");
      }
      if (mcfg->enable.sw_e_bulla1994) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_e_bulla1994");
      }
      if (mcfg->enable.sw_o_bulla1994) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_o_bulla1994");
      }
      if (mcfg->enable.sw_e_mci_pielou1969) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_e_mci_pielou1969");
      }
      if (mcfg->enable.sw_e_prime_camargo1993) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_e_prime_camargo1993");
      }
      if (mcfg->enable.sw_e_var_smith_and_wilson1996_original) {
        fprintf(mcfg->io.f_memory_ptr, "\tsw_e_var_smith_and_wilson1996_original");
      }
    }
    fprintf(mcfg->io.f_memory_ptr, "\n");
  }

  struct measurement_structure_references sref = {
      .g = &g,
      .mst = &mst,
      .heap = &heap,
      .w2v = &w2v,
      .sorted_array_discarded_because_not_in_vector_database = &sorted_array_discarded_because_not_in_vector_database,
  };

  struct measurement_mutables mmut = {
      .best_s = -1.0,
      .prev_best_s = -1.0,
      .prev_num_nodes = 0,
      .mst_initialised = 0,
      .sentence =
          (struct measurement_mutable_counters){
              .num_containing_mwe = 0,
              .num_containing_mwe_tp_only = 0,
              .num_all = 0,
              .count_target = 1,
              .stacked_log = 0.0,
          },
      .document =
          (struct measurement_mutable_counters){
              .num_all = 0,
              .count_target = 1,
              .stacked_log = 0.0,
          },
  };
  if (pthread_mutex_init(&(mmut.mutex), NULL) != 0) {
    perror("Failed to call pthread_mutex_init for mmut\n");
    return 1;
  }

  int64_t num_sentences = 0;
  int64_t max_num_sentences = 10000;
  const int8_t use_max_num_sentences = 0;

  int64_t num_documents = 0;
  int64_t max_num_documents = 10000;
  const uint8_t use_max_num_documents = 0;

  int64_t num_all_sentences = 0; // even ignored sentences
  int64_t max_num_all_sentences = 50000;
  const int8_t use_max_num_all_sentences = 0;

  enum { CUPT, JSONL };

  // pthread_mutex_t misc_var_mutex;
  if (pthread_mutex_init(&mmut.mutex, NULL) != 0) {
    perror("Failed to call pthread_mutex_init for misc_var_mutex\n");
    return 1;
  }

  size_t alloc_size_file_thread = sizeof(pthread_t) * mcfg->threading.num_file_reading_threads;
  pthread_t *const threads_file_reading = malloc(alloc_size_file_thread);
  if (threads_file_reading == NULL) {
    perror("malloc failed\n");
    return 1;
  }
  memset(threads_file_reading, '\0', alloc_size_file_thread);

  size_t alloc_size_file_thread_args = sizeof(struct measurement_file_thread) * mcfg->threading.num_file_reading_threads;
  struct measurement_file_thread *const thread_args_file_reading = malloc(alloc_size_file_thread_args);
  if (thread_args_file_reading == NULL) {
    perror("malloc failed\n");
    return 1;
  }
  memset(thread_args_file_reading, '\0', alloc_size_file_thread_args);

  int32_t i;

  for (i = 0; i < num_input_paths && ((!use_max_num_all_sentences) || num_all_sentences < max_num_all_sentences) &&
              ((!use_max_num_sentences) || num_sentences < max_num_sentences) &&
              ((!use_max_num_documents) || num_documents < max_num_documents);
       i++) {
    int32_t thread_index = i % mcfg->threading.num_file_reading_threads;
    if (i >= mcfg->threading.num_file_reading_threads) {
      if (pthread_join(threads_file_reading[thread_index], NULL) != 0) {
        fprintf(stderr, "Failed to pthread_join for file of index %i; thread_index: %i\n", i, thread_index);
        return 1;
      }
      memset(&(threads_file_reading[thread_index]), '\0', sizeof(pthread_t));
    }

    num_documents++;
    /* // DO NOT REMOVE
    if(input_paths_tp == NULL){
            printf("processing %s\n", input_paths[i]);
    } else {
            printf("processing %s (true positives: %s)\n", input_paths[i], input_paths_tp[i]);
    }
    */
    size_t input_path_len = strlen(input_paths[i]);
    int32_t current_file_format = 0;
    if (strcmp(&(input_paths[i][input_path_len - 5]), ".cupt") == 0) {
      current_file_format = CUPT;
      memset(log_bfr, '\0', log_bfr_size * sizeof(char));
      snprintf(log_bfr, log_bfr_size, "CUPT: %s", input_paths[i]);
      info_format(__FILE__, __func__, __LINE__, log_bfr);
    } else if (strcmp(&(input_paths[i][input_path_len - 6]), ".jsonl") == 0) {
      current_file_format = JSONL;
      memset(log_bfr, '\0', log_bfr_size * sizeof(char));
      snprintf(log_bfr, log_bfr_size, "JSONL: %s", input_paths[i]);
      info_format(__FILE__, __func__, __LINE__, log_bfr);
    } else {
      printf("unknown file type: %s\n", input_paths[i]);
      return 1;
    }

    struct measurement_file_thread mft = {
        .i = i,
        .filename = input_paths[i],
        .mcfg = mcfg,
        .sref = &sref,
        .mmut = &mmut,
    };
    if (mcfg->io.input_path_tp != NULL) {
      mft.filename_tp = input_paths_tp[i];
    } else {
      mft.filename_tp = NULL;
    }
    memcpy(&(thread_args_file_reading[thread_index]), &mft, sizeof(struct measurement_file_thread));
    if (current_file_format == CUPT) {
      if (pthread_create(&(threads_file_reading[thread_index]), NULL, cupt_to_graph_thread,
                         &(thread_args_file_reading[thread_index])) != 0) {
        perror("failed to call pthread_create\n");
        return 1;
      }
    } else if (current_file_format == JSONL) {
      if (pthread_create(&(threads_file_reading[thread_index]), NULL, jsonl_to_graph_thread,
                         &(thread_args_file_reading[thread_index])) != 0) {
        perror("failed to call pthread_create\n");
        return 1;
      }
    }
  }

  for (int32_t j = 0; j < i && j < mcfg->threading.num_file_reading_threads; j++) {
    if (pthread_join(threads_file_reading[j], NULL) != 0) {
      fprintf(stderr, "Failed to pthread_join for file of index %i; thread_index: %i\n", i, j);
      return 1;
    }
  }

  memset(log_bfr, '\0', log_bfr_size);
  snprintf(log_bfr, log_bfr_size, "Running final snapshot with %li nodes", g.num_nodes);
  info_format(__FILE__, __func__, __LINE__, log_bfr);

  int32_t return_status = 0;

  pthread_mutex_lock(&g.mutex_nodes);
  compute_graph_relative_proportions(&g);
  if (zipfian_fit_from_graph(&g, &mmut.best_s) != 0) {
    perror("Failed to call zipfian_fit_from_graph\n");
    return_status = 1;
  } else {
    if (apply_diversity_functions_to_graph(i, mcfg, &sref, &mmut) != 0) {
      perror("Failed to call apply_diversity_functions_to_graph\n");
      return_status = 1;
    }
  }

  pthread_mutex_unlock(&g.mutex_nodes);

  pthread_mutex_destroy(&mmut.mutex);

  free_graph(&g);
  free_graph_distance_heap(&heap);
  free_minimum_spanning_tree(&mst);

  free_sorted_array(&sorted_array_discarded_because_not_in_vector_database);

  free_word2vec(&w2v);

  free(threads_file_reading);
  free(thread_args_file_reading);

  // for(int32_t i = 0 ; i < num_files ; i++){
  for (int32_t i = 0; i < num_input_paths; i++) {
    // free(bfr[i]);
    free(input_paths[i]);
    if (mcfg->io.input_path_tp != NULL) {
      free(input_paths_tp[i]);
    }
  }
  free(input_paths);
  if (mcfg->io.input_path_tp != NULL) {
    free(input_paths_tp);
  }

  fclose(mcfg->io.f_ptr);
  if (mcfg->io.f_timing_ptr != NULL) {
    fclose(mcfg->io.f_timing_ptr);
  }
  if (mcfg->io.f_memory_ptr != NULL) {
    fclose(mcfg->io.f_memory_ptr);
  }

  if (pthread_mutex_destroy(&(mmut.mutex)) != 0) {
    perror("Failed to call pthread_mutex_destroy for mmut\n");
    return_status = 1;
  }

  return return_status;

failure_read_list_simulation_files:
  perror("failed to call read_list_simulation_files\n");
  return 1;
}

int32_t main(int32_t argc, char **argv) {
  char *argv_w2v_path = NULL;
  char *argv_jsonl_content_key = NULL;
  char *argv_input_path = NULL;
  char *argv_input_path_tp = NULL;
  char *argv_output_path = NULL;
  char *argv_output_path_timing = NULL;
  char *argv_output_path_memory = NULL;
  char *argv_udpipe_model_path = NULL;
  // uint32_t argv_tokenization_method = TOKENIZATION_METHOD;
  uint32_t argv_target_column = TARGET_COLUMN;
  int32_t argv_num_row_threads = NUM_ROW_THREADS;
  int32_t argv_num_matrix_threads = NUM_MATRIX_THREADS;
  int32_t argv_num_file_reading_threads = NUM_FILE_READING_THREADS;
  uint8_t argv_enable_token_utf8_normalisation = ENABLE_TOKEN_UTF8_NORMALISATION;
  uint8_t argv_enable_stirling = ENABLE_STIRLING;
  uint8_t argv_enable_ricotta_szeidl = ENABLE_RICOTTA_SZEIDL;
  uint8_t argv_enable_pairwise = ENABLE_PAIRWISE;
  uint8_t argv_enable_lexicographic = ENABLE_LEXICOGRAPHIC;
  uint8_t argv_enable_chao_et_al_functional_diversity = ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY;
  uint8_t argv_enable_scheiner_species_phylogenetic_functional_diversity =
      ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY;
  uint8_t argv_enable_leinster_cobbold_diversity = ENABLE_LEINSTER_COBBOLD_DIVERSITY;
  uint8_t argv_enable_multithreaded_matrix_generation = ENABLE_MULTITHREADED_MATRIX_GENERATION;
  uint8_t argv_enable_timings = ENABLE_TIMINGS;
  uint8_t argv_enable_iterative_distance_computation = ENABLE_ITERATIVE_DISTANCE_COMPUTATION;
  uint8_t argv_enable_multithreaded_row_generation = ENABLE_MULTITHREADED_ROW_GENERATION;
  int8_t argv_row_generation_batch_size = ROW_GENERATION_BATCH_SIZE;
  uint8_t argv_enable_sentence_count_recompute_step = ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP;
  uint8_t argv_sentence_recompute_step_use_log10 = SENTENCE_RECOMPUTE_STEP_USE_LOG10;
  uint8_t argv_enable_document_count_recompute_step = ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP;
  uint8_t argv_document_recompute_step_use_log10 = DOCUMENT_RECOMPUTE_STEP_USE_LOG10;
  uint8_t argv_enable_output_timing = ENABLE_OUTPUT_TIMING;
  uint8_t argv_enable_output_memory = ENABLE_OUTPUT_MEMORY;
  uint64_t argv_sentence_count_recompute_step = SENTENCE_COUNT_RECOMPUTE_STEP;
  double argv_sentence_count_recompute_step_log10 = SENTENCE_COUNT_RECOMPUTE_STEP_LOG10;
  uint64_t argv_document_count_recompute_step = DOCUMENT_COUNT_RECOMPUTE_STEP;
  double argv_document_count_recompute_step_log10 = DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10;
  uint8_t argv_enable_functional_evenness = ENABLE_FUNCTIONAL_EVENNESS;
  uint8_t argv_enable_functional_dispersion = ENABLE_FUNCTIONAL_DISPERSION;
  uint8_t argv_enable_functional_divergence_modified = ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED;
  uint8_t argv_enable_non_disparity_functions = ENABLE_NON_DISPARITY_FUNCTIONS;
  uint8_t argv_enable_disparity_functions = ENABLE_DISPARITY_FUNCTIONS;
  uint8_t argv_enable_shannon_weaver_entropy = ENABLE_SHANNON_WEAVER_ENTROPY;
  uint8_t argv_enable_good_entropy = ENABLE_GOOD_ENTROPY;
  uint8_t argv_enable_renyi_entropy = ENABLE_RENYI_ENTROPY;
  uint8_t argv_enable_patil_taillie_entropy = ENABLE_PATIL_TAILLIE_ENTROPY;
  uint8_t argv_enable_q_logarithmic_entropy = ENABLE_Q_LOGARITHMIC_ENTROPY;
  uint8_t argv_enable_simpson_index = ENABLE_SIMPSON_INDEX;
  uint8_t argv_enable_simpson_dominance_index = ENABLE_SIMPSON_DOMINANCE_INDEX;
  uint8_t argv_enable_hill_number_standard = ENABLE_HILL_NUMBER_STANDARD;
  uint8_t argv_enable_hill_evenness = ENABLE_HILL_EVENNESS;
  uint8_t argv_enable_berger_parker_index = ENABLE_BERGER_PARKER_INDEX;
  uint8_t argv_enable_junge1994_page22 = ENABLE_JUNGE1994_PAGE22;
  uint8_t argv_enable_brillouin_diversity = ENABLE_BRILLOUIN_DIVERSITY;
  uint8_t argv_enable_mcintosh_index = ENABLE_MCINTOSH_INDEX;
  uint8_t argv_enable_sw_entropy_over_log_n_species_pielou1975 = ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975;
  uint8_t argv_enable_sw_e_heip = ENABLE_SW_E_HEIP;
  uint8_t argv_enable_sw_e_one_minus_d = ENABLE_SW_E_ONE_MINUS_D;
  uint8_t argv_enable_sw_e_one_over_ln_d_williams1964 = ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964;
  uint8_t argv_enable_sw_e_minus_ln_d_pielou1977 = ENABLE_SW_E_MINUS_LN_D_PIELOU1977;
  uint8_t argv_enable_sw_f_2_1_alatalo1981 = ENABLE_SW_F_2_1_ALATALO1981;
  uint8_t argv_enable_sw_g_2_1_molinari1989 = ENABLE_SW_G_2_1_MOLINARI1989;
  uint8_t argv_enable_sw_e_bulla1994 = ENABLE_SW_E_BULLA1994;
  uint8_t argv_enable_sw_o_bulla1994 = ENABLE_SW_O_BULLA1994;
  uint8_t argv_enable_sw_e_mci_pielou1969 = ENABLE_SW_E_MCI_PIELOU1969;
  uint8_t argv_enable_sw_e_prime_camargo1993 = ENABLE_SW_E_PRIME_CAMARGO1993;
  uint8_t argv_enable_sw_e_prime_camargo1993_multithreading = ENABLE_SW_E_PRIME_CAMARGO1993_MULTITHREADING;
  uint8_t argv_enable_sw_e_var_smith_and_wilson1996_original = ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL;
  double argv_stirling_alpha = STIRLING_ALPHA;
  double argv_stirling_beta = STIRLING_BETA;
  double argv_ricotta_szeidl_alpha = RICOTTA_SZEIDL_ALPHA;
  double argv_chao_et_al_functional_diversity_alpha = CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA;
  double argv_scheiner_species_phylogenetic_functional_diversity_alpha =
      SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA;
  double argv_leinster_cobbold_diversity_alpha = LEINSTER_COBBOLD_DIVERSITY_ALPHA;
  double argv_good_alpha = GOOD_ALPHA;
  double argv_good_beta = GOOD_BETA;
  double argv_renyi_alpha = RENYI_ALPHA;
  double argv_patil_taillie_alpha = PATIL_TAILLIE_ALPHA;
  double argv_q_logarithmic_q = Q_LOGARITHMIC_Q;
  double argv_hill_number_standard_alpha = HILL_NUMBER_STANDARD_ALPHA;
  double argv_hill_evenness_alpha = HILL_EVENNESS_ALPHA;
  double argv_hill_evenness_beta = HILL_EVENNESS_BETA;

  uint8_t argv_force_timing_and_memory_to_output_path = 0;

  for (int32_t i = 1; i < argc; i++) {
    if (strncmp(argv[i], "--w2v_path=", 11) == 0) {
      argv_w2v_path = argv[i] + 11;
    } else if (strncmp(argv[i], "--target_column=", 16) == 0) {
      // argv_target_column = argv[i] + 16;
      if (strcmp(argv[i] + 16, "UD_MWE") == 0) {
        argv_target_column = UD_MWE;
      } else if (strcmp(argv[i] + 16, "UD_FORM") == 0) {
        argv_target_column = UD_FORM;
      } else if (strcmp(argv[i] + 16, "UD_LEMMA") == 0) {
        argv_target_column = UD_LEMMA;
      } else {
        fprintf(stderr, "Unknown target_column: %s\n", argv[i] + 16);
        return 1;
      }
    } else if (strncmp(argv[i], "--num_row_threads=", 18) == 0) {
      argv_num_row_threads = (int32_t)strtol(argv[i] + 18, NULL, 10);
    } else if (strncmp(argv[i], "--num_matrix_threads=", 21) == 0) {
      argv_num_matrix_threads = (int32_t)strtol(argv[i] + 21, NULL, 10);
    } else if (strncmp(argv[i], "--num_file_reading_threads=", 27) == 0) {
      argv_num_file_reading_threads = (int32_t)strtol(argv[i] + 27, NULL, 10);
    } else if (strncmp(argv[i], "--jsonl_content_key=", 20) == 0) {
      argv_jsonl_content_key = argv[i] + 20;
    } else if (strncmp(argv[i], "--input_path=", 13) == 0) {
      argv_input_path = argv[i] + 13;
    } else if (strncmp(argv[i], "--input_path_tp=", 16) == 0) {
      argv_input_path_tp = argv[i] + 16;
    } else if (strncmp(argv[i], "--output_path=", 14) == 0) {
      argv_output_path = argv[i] + 14;
    } else if (strncmp(argv[i], "--output_path_timing=", 21) == 0) {
      argv_output_path_timing = argv[i] + 21;
    } else if (strncmp(argv[i], "--output_path_memory=", 21) == 0) {
      argv_output_path_memory = argv[i] + 21;
    } else if (strncmp(argv[i], "--udpipe_model_path=", 20) == 0) {
      argv_udpipe_model_path = argv[i] + 20;
    } else if (strncmp(argv[i], "--enable_multithreaded_matrix_generation=", 41) == 0) {
      argv_enable_multithreaded_matrix_generation = (argv[i][41] == '1');
    } else if (strncmp(argv[i], "--enable_timings=", 17) == 0) {
      argv_enable_timings = (argv[i][17] == '1');
    } else if (strncmp(argv[i], "--enable_iterative_distance_computation=", 40) == 0) {
      argv_enable_iterative_distance_computation = (argv[i][40] == '1');
    } else if (strncmp(argv[i], "--enable_sentence_count_recompute_step=", 39) == 0) {
      argv_enable_sentence_count_recompute_step = (argv[i][39] == '1');
    } else if (strncmp(argv[i], "--enable_document_count_recompute_step=", 39) == 0) {
      argv_enable_document_count_recompute_step = (argv[i][39] == '1');
    } else if (strncmp(argv[i], "--sentence_recompute_step_use_log10=", 36) == 0) {
      argv_sentence_recompute_step_use_log10 = (argv[i][36] == '1');
    } else if (strncmp(argv[i], "--document_recompute_step_use_log10=", 36) == 0) {
      argv_document_recompute_step_use_log10 = (argv[i][36] == '1');
    } else if (strncmp(argv[i], "--enable_output_timing=", 23) == 0) {
      argv_enable_output_timing = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_output_memory=", 23) == 0) {
      argv_enable_output_memory = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_token_utf8_normalisation=", 34) == 0) {
      argv_enable_token_utf8_normalisation = (argv[i][34] == '1');
    } else if (strncmp(argv[i], "--enable_stirling=", 18) == 0) {
      argv_enable_stirling = (argv[i][18] == '1');
    } else if (strncmp(argv[i], "--enable_ricotta_szeidl=", 24) == 0) {
      argv_enable_ricotta_szeidl = (argv[i][24] == '1');
    } else if (strncmp(argv[i], "--enable_pairwise=", 18) == 0) {
      argv_enable_pairwise = (argv[i][18] == '1');
    } else if (strncmp(argv[i], "--enable_lexicographic=", 23) == 0) {
      argv_enable_lexicographic = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_chao_et_al_functional_diversity=", 41) == 0) {
      argv_enable_chao_et_al_functional_diversity = (argv[i][41] == '1');
    } else if (strncmp(argv[i], "--enable_scheiner_species_phylogenetic_functional_diversity=", 60) == 0) {
      argv_enable_scheiner_species_phylogenetic_functional_diversity = (argv[i][60] == '1');
    } else if (strncmp(argv[i], "--enable_leinster_cobbold_diversity=", 36) == 0) {
      argv_enable_leinster_cobbold_diversity = (argv[i][36] == '1');
    } else if (strncmp(argv[i], "--enable_functional_evenness=", 29) == 0) {
      argv_enable_functional_evenness = (argv[i][29] == '1');
    } else if (strncmp(argv[i], "--enable_functional_dispersion=", 31) == 0) {
      argv_enable_functional_dispersion = (argv[i][31] == '1');
    } else if (strncmp(argv[i], "--enable_functional_divergence_modified=", 40) == 0) {
      argv_enable_functional_divergence_modified = (argv[i][40] == '1');
    } else if (strncmp(argv[i], "--enable_non_disparity_functions=", 33) == 0) {
      argv_enable_non_disparity_functions = (argv[i][33] == '1');
    } else if (strncmp(argv[i], "--enable_disparity_functions=", 29) == 0) {
      argv_enable_disparity_functions = (argv[i][29] == '1');
    } else if (strncmp(argv[i], "--enable_shannon_weaver_entropy=", 32) == 0) {
      argv_enable_shannon_weaver_entropy = (argv[i][32] == '1');
    } else if (strncmp(argv[i], "--enable_good_entropy=", 22) == 0) {
      argv_enable_good_entropy = (argv[i][22] == '1');
    } else if (strncmp(argv[i], "--enable_renyi_entropy=", 23) == 0) {
      argv_enable_renyi_entropy = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_patil_taillie_entropy=", 31) == 0) {
      argv_enable_patil_taillie_entropy = (argv[i][31] == '1');
    } else if (strncmp(argv[i], "--enable_q_logarithmic_entropy=", 31) == 0) {
      argv_enable_q_logarithmic_entropy = (argv[i][31] == '1');
    } else if (strncmp(argv[i], "--enable_simpson_index=", 23) == 0) {
      argv_enable_simpson_index = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_simpson_dominance_index=", 33) == 0) {
      argv_enable_simpson_dominance_index = (argv[i][33] == '1');
    } else if (strncmp(argv[i], "--enable_hill_number_standard=", 30) == 0) {
      argv_enable_hill_number_standard = (argv[i][30] == '1');
    } else if (strncmp(argv[i], "--enable_hill_evenness=", 23) == 0) {
      argv_enable_hill_evenness = (argv[i][23] == '1');
    } else if (strncmp(argv[i], "--enable_berger_parker_index=", 29) == 0) {
      argv_enable_berger_parker_index = (argv[i][29] == '1');
    } else if (strncmp(argv[i], "--enable_junge1994_page22=", 26) == 0) {
      argv_enable_junge1994_page22 = (argv[i][26] == '1');
    } else if (strncmp(argv[i], "--enable_brillouin_diversity=", 29) == 0) {
      argv_enable_brillouin_diversity = (argv[i][29] == '1');
    } else if (strncmp(argv[i], "--enable_mcintosh_index=", 24) == 0) {
      argv_enable_mcintosh_index = (argv[i][24] == '1');
    } else if (strncmp(argv[i], "--enable_sw_entropy_over_log_n_species_pielou1975=", 50) == 0) {
      argv_enable_sw_entropy_over_log_n_species_pielou1975 = (argv[i][50] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_heip=", 19) == 0) {
      argv_enable_sw_e_heip = (argv[i][19] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_one_over_ln_d_williams1964=", 41) == 0) {
      argv_enable_sw_e_one_over_ln_d_williams1964 = (argv[i][41] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_minus_ln_d_pielou1977=", 36) == 0) {
      argv_enable_sw_e_minus_ln_d_pielou1977 = (argv[i][36] == '1');
    } else if (strncmp(argv[i], "--enable_sw_f_2_1_alatalo1981=", 30) == 0) {
      argv_enable_sw_f_2_1_alatalo1981 = (argv[i][30] == '1');
    } else if (strncmp(argv[i], "--enable_sw_g_2_1_molinari1989=", 31) == 0) {
      argv_enable_sw_g_2_1_molinari1989 = (argv[i][31] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_bulla1994=", 24) == 0) {
      argv_enable_sw_e_bulla1994 = (argv[i][24] == '1');
    } else if (strncmp(argv[i], "--enable_sw_o_bulla1994=", 24) == 0) {
      argv_enable_sw_o_bulla1994 = (argv[i][24] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_mci_pielou1969=", 29) == 0) {
      argv_enable_sw_e_mci_pielou1969 = (argv[i][29] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_prime_camargo1993=", 32) == 0) {
      argv_enable_sw_e_prime_camargo1993 = (argv[i][32] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_prime_camargo1993_multithreading=", 47) == 0) {
      argv_enable_sw_e_prime_camargo1993_multithreading = (argv[i][47] == '1');
    } else if (strncmp(argv[i], "--enable_sw_e_var_smith_and_wilson1996_original=", 48) == 0) {
      argv_enable_sw_e_var_smith_and_wilson1996_original = (argv[i][48] == '1');
    } else if (strncmp(argv[i], "--stirling_alpha=", 17) == 0) {
      argv_stirling_alpha = strtod(argv[i] + 17, NULL);
    } else if (strncmp(argv[i], "--stirling_beta=", 16) == 0) {
      argv_stirling_beta = strtod(argv[i] + 16, NULL);
    } else if (strncmp(argv[i], "--ricotta_szeidl_alpha=", 23) == 0) {
      argv_ricotta_szeidl_alpha = strtod(argv[i] + 23, NULL);
    } else if (strncmp(argv[i], "--chao_et_al_functional_diversity_alpha=", 40) == 0) {
      argv_chao_et_al_functional_diversity_alpha = strtod(argv[i] + 40, NULL);
    } else if (strncmp(argv[i], "--scheiner_species_phylogenetic_functional_diversity_alpha=", 59) == 0) {
      argv_scheiner_species_phylogenetic_functional_diversity_alpha = strtod(argv[i] + 59, NULL);
    } else if (strncmp(argv[i], "--leinster_cobbold_diversity_alpha=", 35) == 0) {
      argv_leinster_cobbold_diversity_alpha = strtod(argv[i] + 35, NULL);
    } else if (strncmp(argv[i], "--renyi_alpha=", 14) == 0) {
      argv_renyi_alpha = strtod(argv[i] + 14, NULL);
    } else if (strncmp(argv[i], "--patil_taillie_alpha=", 22) == 0) {
      argv_patil_taillie_alpha = strtod(argv[i] + 22, NULL);
    } else if (strncmp(argv[i], "--hill_number_standard_alpha=", 29) == 0) {
      argv_hill_number_standard_alpha = strtod(argv[i] + 29, NULL);
    } else if (strncmp(argv[i], "--hill_evenness_alpha=", 22) == 0) {
      argv_hill_evenness_alpha = strtod(argv[i] + 22, NULL);
    } else if (strncmp(argv[i], "--hill_evenness_beta=", 21) == 0) {
      argv_hill_evenness_beta = strtod(argv[i] + 21, NULL);
    } else if (strncmp(argv[i], "--good_alpha=", 13) == 0) {
      argv_good_alpha = strtod(argv[i] + 13, NULL);
    } else if (strncmp(argv[i], "--good_beta=", 12) == 0) {
      argv_good_beta = strtod(argv[i] + 12, NULL);
    } else if (strncmp(argv[i], "--row_generation_batch_size=", 28) == 0) {
      argv_row_generation_batch_size = (int8_t)strtol(argv[i] + 28, NULL, 10);
    } else if (strncmp(argv[i], "--sentence_count_recompute_step=", 32) == 0) {
      argv_sentence_count_recompute_step = strtol(argv[i] + 32, NULL, 10);
    } else if (strncmp(argv[i], "--sentence_count_recompute_step_log10=", 38) == 0) {
      argv_sentence_count_recompute_step_log10 = strtod(argv[i] + 38, NULL);
    } else if (strncmp(argv[i], "--document_count_recompute_step=", 32) == 0) {
      argv_document_count_recompute_step = strtol(argv[i] + 32, NULL, 10);
    } else if (strncmp(argv[i], "--document_count_recompute_step_log10=", 38) == 0) {
      argv_document_count_recompute_step_log10 = strtod(argv[i] + 38, NULL);
    } else if (strncmp(argv[i], "--force_timing_and_memory_to_output_path=", 41) == 0) {
      argv_force_timing_and_memory_to_output_path = (argv[i][41] == '1');
    }
    // else if(strncmp(argv[i], "--tokenization_method=", 22) == 0){argv_tokenization_method = strtol(argv[i] + 22, NULL, 10);}
    else {
      fprintf(stderr, "Unknown argument: %s\n", argv[i]);
      return 1;
    }
  }

  if (argv_w2v_path == NULL) {
    argv_w2v_path = W2V_PATH;
  }
  if (argv_jsonl_content_key == NULL) {
    argv_jsonl_content_key = JSONL_CONTENT_KEY;
  }
  if (argv_input_path == NULL) {
    argv_input_path = INPUT_PATH;
  }
  if (argv_output_path == NULL) {
    argv_output_path = OUTPUT_PATH;
  }
  if (argv_output_path_timing == NULL) {
    argv_output_path_timing = OUTPUT_PATH_TIMING;
  }
  if (argv_output_path_memory == NULL) {
    argv_output_path_memory = OUTPUT_PATH_MEMORY;
  }
  if (argv_udpipe_model_path == NULL && TOKENIZATION_METHOD == 2) {
    fprintf(stderr, "Setting TOKENIZATION_METHOD to 2 (embedded UDPipe) requires defining --udpipe_model_path=...\n");
    return 1;
  }

  if (argv_force_timing_and_memory_to_output_path) {
    size_t path_len, delta, alloc_size, len_suffix;
    char *ptr_last_slash;

    ptr_last_slash = strrchr(argv_output_path, '/');
    if (ptr_last_slash == NULL) {
      ptr_last_slash = argv_output_path - 1;
    }
    // ptr_last_slash++;

    delta = ((size_t)(ptr_last_slash - argv_output_path)) + 1;
    len_suffix = strlen("measurement_output_timing.tsv");
    path_len = delta + len_suffix;
    alloc_size = path_len + 1;
    argv_output_path_timing = malloc(alloc_size);
    if (argv_output_path_timing == NULL) {
      goto malloc_fail;
    }
    memset(argv_output_path_timing, '\0', alloc_size);
    memcpy(argv_output_path_timing, argv_output_path, delta);
    memcpy(argv_output_path_timing + delta, "measurement_output_timing.tsv", len_suffix);

    len_suffix = strlen("measurement_output_memory.tsv");
    path_len = delta + len_suffix;
    alloc_size = path_len + 1;
    argv_output_path_memory = malloc(alloc_size);
    if (argv_output_path_memory == NULL) {
      free(argv_output_path_timing);
      goto malloc_fail;
    }
    memset(argv_output_path_memory, '\0', alloc_size);
    memcpy(argv_output_path_memory, argv_output_path, delta);
    memcpy(argv_output_path_memory + delta, "measurement_output_memory.tsv", len_suffix);
  }

  if (argv_num_row_threads == -1 || argv_num_matrix_threads == -1 || argv_row_generation_batch_size == -1 ||
      argv_num_file_reading_threads == -1) {
    struct cpu_info local_cpu_info = {0};
    if (get_cpu_info(&local_cpu_info) != 0) {
      perror("Failed to call get_cpu_info\n");
      return 1;
    }
    const int32_t log_bfr_size = 256;
    char log_bfr[256];
    memset(log_bfr, '\0', log_bfr_size);
    snprintf(log_bfr, log_bfr_size, "Number of virtual cores: %u", local_cpu_info.cardinality_virtual_cores);
    info_format(__FILE__, __func__, __LINE__, log_bfr);

    memset(log_bfr, '\0', log_bfr_size);
    if (local_cpu_info.avx256_capable) {
      snprintf(log_bfr, log_bfr_size, "CPU can do AVX256 instructions");
    } else {
      snprintf(log_bfr, log_bfr_size, "CPU cannot do AVX256 instructions");
    }
    info_format(__FILE__, __func__, __LINE__, log_bfr);

    memset(log_bfr, '\0', log_bfr_size);
    if (local_cpu_info.avx512_capable) {
      snprintf(log_bfr, log_bfr_size, "CPU can do AVX512 instructions");
    } else {
      snprintf(log_bfr, log_bfr_size, "CPU cannot do AVX512 instructions");
    }
    info_format(__FILE__, __func__, __LINE__, log_bfr);

    if (argv_num_row_threads == -1) {
      argv_num_row_threads = (int32_t)local_cpu_info.cardinality_virtual_cores;
    }
    if (argv_num_matrix_threads == -1) {
      argv_num_matrix_threads = (int32_t)local_cpu_info.cardinality_virtual_cores;
    }
    if (argv_row_generation_batch_size == -1) {
      argv_row_generation_batch_size = (int32_t)local_cpu_info.cardinality_virtual_cores;
    }
    if (argv_num_file_reading_threads == -1) {
      argv_num_file_reading_threads = (int32_t)local_cpu_info.cardinality_virtual_cores;
    }
  }

  printf("w2v_path: %s\n", argv_w2v_path);
  printf("jsonl_content_key: %s\n", argv_jsonl_content_key);
  printf("input_path: %s\n", argv_input_path);
  if (argv_input_path_tp != NULL) {
    printf("input_path_tp: %s\n", argv_input_path_tp);
  }
  printf("output_path: %s\n", argv_output_path);
  printf("force_timing_and_memory_to_output_path: %u\n", argv_force_timing_and_memory_to_output_path);
  printf("output_path_timing: %s\n", argv_output_path_timing);
  printf("output_path_memory: %s\n", argv_output_path_memory);
#if TOKENIZATION_METHOD == 2
  printf("udpipe_model_path: %s\n", argv_udpipe_model_path);
#endif

  printf("target_column: %u\n", argv_target_column);
  printf("enable_token_utf8_normalisation: %u\n", argv_enable_token_utf8_normalisation);
  printf("num_row_threads: %i\n", argv_num_row_threads);
  printf("num_matrix_threads: %i\n", argv_num_matrix_threads);
  printf("num_file_reading_threads: %i\n", argv_num_file_reading_threads);

  printf("enable_multithreaded_matrix_generation: %u\n", argv_enable_multithreaded_matrix_generation);
  printf("enable_timings: %u\n", argv_enable_timings);
  printf("enable_iterative_distance_computation: %u\n", argv_enable_iterative_distance_computation);
  printf("enable_multithreaded_row_generation: %u\n", argv_enable_multithreaded_row_generation);
  printf("row_generation_batch_size: %i\n", argv_row_generation_batch_size);
  printf("enable_sw_e_prime_camargo1993_multithreading: %u\n", argv_enable_sw_e_prime_camargo1993_multithreading);
  printf("sentence_count_recompute_step: %lu\n", argv_sentence_count_recompute_step);
  printf("enable_sentence_count_recompute_step: %u\n", argv_enable_sentence_count_recompute_step);
  printf("sentence_recompute_step_use_log10: %u\n", argv_sentence_recompute_step_use_log10);
  printf("sentence_count_recompute_step_log10: %f\n", argv_sentence_count_recompute_step_log10);
  printf("document_count_recompute_step: %lu\n", argv_document_count_recompute_step);
  printf("enable_document_count_recompute_step: %u\n", argv_enable_document_count_recompute_step);
  printf("document_recompute_step_use_log10: %u\n", argv_document_recompute_step_use_log10);
  printf("document_count_recompute_step_log10: %f\n", argv_document_count_recompute_step_log10);
  printf("enable_output_timing: %u\n", argv_enable_output_timing);
  printf("enable_output_memory: %u\n", argv_enable_output_memory);

  printf("enable_stirling: %u\n", argv_enable_stirling);
  printf("enable_ricotta_szeidl: %u\n", argv_enable_ricotta_szeidl);
  printf("enable_pairwise: %u\n", argv_enable_pairwise);
  printf("enable_lexicographic: %u\n", argv_enable_lexicographic);
  printf("enable_chao_et_al_functional_diversity: %u\n", argv_enable_chao_et_al_functional_diversity);
  printf("enable_scheiner_species_phylogenetic_functional_diversity: %u\n",
         argv_enable_scheiner_species_phylogenetic_functional_diversity);
  printf("enable_leinster_cobbold_diversity: %u\n", argv_enable_leinster_cobbold_diversity);
  printf("enable_functional_evenness: %u\n", argv_enable_functional_evenness);
  printf("enable_functional_dispersion: %u\n", argv_enable_functional_dispersion);
  printf("enable_functional_divergence_modified: %u\n", argv_enable_functional_divergence_modified);
  printf("enable_non_disparity_functions: %u\n", argv_enable_non_disparity_functions);
  printf("enable_disparity_functions: %u\n", argv_enable_disparity_functions);
  printf("enable_shannon_weaver_entropy: %u\n", argv_enable_shannon_weaver_entropy);
  printf("enable_good_entropy: %u\n", argv_enable_good_entropy);
  printf("enable_renyi_entropy: %u\n", argv_enable_renyi_entropy);
  printf("enable_patil_taillie_entropy: %u\n", argv_enable_patil_taillie_entropy);
  printf("enable_q_logarithmic_entropy: %u\n", argv_enable_q_logarithmic_entropy);
  printf("enable_simpson_index: %u\n", argv_enable_simpson_index);
  printf("enable_simpson_dominance_index: %u\n", argv_enable_simpson_dominance_index);
  printf("enable_hill_number_standard: %u\n", argv_enable_hill_number_standard);
  printf("enable_hill_evenness: %u\n", argv_enable_hill_evenness);
  printf("enable_berger_parker_index: %u\n", argv_enable_berger_parker_index);
  printf("enable_junge1994_page22: %u\n", argv_enable_junge1994_page22);
  printf("enable_brillouin_diversity: %u\n", argv_enable_brillouin_diversity);
  printf("enable_mcintosh_index: %u\n", argv_enable_mcintosh_index);
  printf("enable_sw_entropy_over_log_n_species_pielou1975: %u\n", argv_enable_sw_entropy_over_log_n_species_pielou1975);
  printf("enable_sw_e_heip: %u\n", argv_enable_sw_e_heip);
  printf("enable_sw_e_one_minus_d: %u\n", argv_enable_sw_e_one_minus_d);
  printf("enable_sw_e_one_over_ln_d_williams1964: %u\n", argv_enable_sw_e_one_over_ln_d_williams1964);
  printf("enable_sw_e_minus_ln_d_pielou1977: %u\n", argv_enable_sw_e_minus_ln_d_pielou1977);
  printf("enable_sw_f_2_1_alatalo1981: %u\n", argv_enable_sw_f_2_1_alatalo1981);
  printf("enable_sw_g_2_1_molinari1989: %u\n", argv_enable_sw_g_2_1_molinari1989);
  printf("enable_sw_e_bulla1994: %u\n", argv_enable_sw_e_bulla1994);
  printf("enable_sw_o_bulla1994: %u\n", argv_enable_sw_o_bulla1994);
  printf("enable_sw_e_mci_pielou1969: %u\n", argv_enable_sw_e_mci_pielou1969);
  printf("enable_sw_e_prime_camargo1993: %u\n", argv_enable_sw_e_prime_camargo1993);
  printf("enable_sw_e_var_smith_and_wilson1996_original: %u\n", argv_enable_sw_e_var_smith_and_wilson1996_original);

  printf("stirling_alpha: %f\n", argv_stirling_alpha);
  printf("stirling_beta: %f\n", argv_stirling_beta);
  printf("ricotta_szeidl_alpha: %f\n", argv_ricotta_szeidl_alpha);
  printf("chao_et_al_functional_diversity_alpha: %f\n", argv_chao_et_al_functional_diversity_alpha);
  printf("scheiner_species_phylogenetic_functional_diversity_alpha: %f\n",
         argv_scheiner_species_phylogenetic_functional_diversity_alpha);
  printf("leinster_cobbold_diversity_alpha: %f\n", argv_leinster_cobbold_diversity_alpha);
  printf("good_alpha: %f\n", argv_good_alpha);
  printf("good_beta: %f\n", argv_good_beta);
  printf("renyi_alpha: %f\n", argv_renyi_alpha);
  printf("patil_taillie_alpha: %f\n", argv_patil_taillie_alpha);
  printf("q_logarithmic_q: %f\n", argv_q_logarithmic_q);
  printf("hill_number_standard_alpha: %f\n", argv_hill_number_standard_alpha);
  printf("hill_evenness_alpha: %f\n", argv_hill_evenness_alpha);
  printf("hill_evenness_beta: %f\n", argv_hill_evenness_beta);

#if TOKENIZATION_METHOD == 2
  ensure_proper_udpipe_pipeline_size();
#elif TOKENIZATION_METHOD == 1
  jsonl_init_tokenization();
#endif

#if ENABLE_FILTER == 1
  if (filter_ready() != 0) {
    perror("Failed to call filter_ready\n");
    return 1;
  }
#endif

  int32_t err;

#if TOKENIZATION_METHOD == 2
  // udpipe_pipeline_create_global("/home/esteve/Documents/thesis/other_repos/udpipe/sandbox_models/english-ewt-ud-2.5-191206.udpipe",
  // "tokenizer", "none", "none", "vertical"); // unsure about "none" for parser
  udpipe_pipeline_create_global(argv_udpipe_model_path, "tokenizer", "none", "none",
                                "vertical"); // unsure about "none" for parser
// udpipe_pipeline_print_global_info();
#endif

  struct measurement_configuration mcfg = {
      .target_column = argv_target_column,
      .enable_token_utf8_normalisation = argv_enable_token_utf8_normalisation,
      .jsonl_content_key = argv_jsonl_content_key,
      .div_param =
          (struct measurement_diversity_parameters){
              .stirling_alpha = argv_stirling_alpha,
              .stirling_beta = argv_stirling_beta,
              .ricotta_szeidl_alpha = argv_ricotta_szeidl_alpha,
              .chao_et_al_functional_diversity_alpha = argv_chao_et_al_functional_diversity_alpha,
              .scheiner_species_phylogenetic_functional_diversity_alpha =
                  argv_scheiner_species_phylogenetic_functional_diversity_alpha,
              .leinster_cobbold_diversity_alpha = argv_leinster_cobbold_diversity_alpha,
              .good_alpha = argv_good_alpha,
              .good_beta = argv_good_beta,
              .renyi_alpha = argv_renyi_alpha,
              .patil_taillie_alpha = argv_patil_taillie_alpha,
              .q_logarithmic_q = argv_q_logarithmic_q,
              .hill_number_standard_alpha = argv_hill_number_standard_alpha,
              .hill_evenness_alpha = argv_hill_evenness_alpha,
              .hill_evenness_beta = argv_hill_evenness_beta,
          },
      .enable =
          (struct measurement_diversity_enabler){
              .stirling = argv_enable_stirling,
              .ricotta_szeidl = argv_enable_ricotta_szeidl,
              .pairwise = argv_enable_pairwise,
              .lexicographic = argv_enable_lexicographic,
              .chao_et_al_functional_diversity = argv_enable_chao_et_al_functional_diversity,
              .scheiner_species_phylogenetic_functional_diversity =
                  argv_enable_scheiner_species_phylogenetic_functional_diversity,
              .leinster_cobbold_diversity = argv_enable_leinster_cobbold_diversity,
              .functional_evenness = argv_enable_functional_evenness,
              .functional_dispersion = argv_enable_functional_dispersion,
              .functional_divergence_modified = argv_enable_functional_divergence_modified,
              .non_disparity_functions = argv_enable_non_disparity_functions,
              .disparity_functions = argv_enable_disparity_functions,
              .shannon_weaver_entropy = argv_enable_shannon_weaver_entropy,
              .good_entropy = argv_enable_good_entropy,
              .renyi_entropy = argv_enable_renyi_entropy,
              .patil_taillie_entropy = argv_enable_patil_taillie_entropy,
              .q_logarithmic_entropy = argv_enable_q_logarithmic_entropy,
              .simpson_index = argv_enable_simpson_index,
              .simpson_dominance_index = argv_enable_simpson_dominance_index,
              .hill_number_standard = argv_enable_hill_number_standard,
              .hill_evenness = argv_enable_hill_evenness,
              .berger_parker_index = argv_enable_berger_parker_index,
              .junge1994_page22 = argv_enable_junge1994_page22,
              .brillouin_diversity = argv_enable_brillouin_diversity,
              .mcintosh_index = argv_enable_mcintosh_index,
              .sw_entropy_over_log_n_species_pielou1975 = argv_enable_sw_entropy_over_log_n_species_pielou1975,
              .sw_e_heip = argv_enable_sw_e_heip,
              .sw_e_one_minus_d = argv_enable_sw_e_one_minus_d,
              .sw_e_one_over_ln_d_williams1964 = argv_enable_sw_e_one_over_ln_d_williams1964,
              .sw_e_minus_ln_d_pielou1977 = argv_enable_sw_e_minus_ln_d_pielou1977,
              .sw_f_2_1_alatalo1981 = argv_enable_sw_f_2_1_alatalo1981,
              .sw_g_2_1_molinari1989 = argv_enable_sw_g_2_1_molinari1989,
              .sw_e_bulla1994 = argv_enable_sw_e_bulla1994,
              .sw_o_bulla1994 = argv_enable_sw_o_bulla1994,
              .sw_e_mci_pielou1969 = argv_enable_sw_e_mci_pielou1969,
              .sw_e_prime_camargo1993 = argv_enable_sw_e_prime_camargo1993,
              .sw_e_var_smith_and_wilson1996_original = argv_enable_sw_e_var_smith_and_wilson1996_original,
          },
      .io =
          (struct measurement_io){
              .w2v_path = argv_w2v_path,
              .jsonl_content_key = argv_jsonl_content_key,
              .input_path = argv_input_path,
              .input_path_tp = argv_input_path_tp,
              .output_path = argv_output_path,
              .output_path_timing = argv_output_path_timing,
              .output_path_memory = argv_output_path_memory,
              .udpipe_model_path = argv_udpipe_model_path,
              .f_ptr = NULL,
              .f_timing_ptr = NULL,
              .f_memory_ptr = NULL,
              .enable_timings = argv_enable_timings,
              .enable_output_timing = argv_enable_output_timing,
              .enable_output_memory = argv_enable_output_memory,
              .force_timing_and_memory_to_output_path = argv_force_timing_and_memory_to_output_path,
          },
      .threading =
          (struct measurement_threading){
              .num_row_threads = argv_num_row_threads,
              .num_matrix_threads = argv_num_matrix_threads,
              .num_file_reading_threads = argv_num_file_reading_threads,
              .enable_multithreaded_matrix_generation = argv_enable_multithreaded_matrix_generation,
              .enable_iterative_distance_computation = argv_enable_iterative_distance_computation,
              .enable_multithreaded_row_generation = argv_enable_multithreaded_row_generation,
              .row_generation_batch_size = argv_row_generation_batch_size,
              .enable_sw_e_prime_camargo1993_multithreading = argv_enable_sw_e_prime_camargo1993_multithreading,
          },
      .steps =
          (struct measurement_step_parameters){
              .sentence =
                  (struct measurement_step){
                      .enable_count_recompute_step = argv_enable_sentence_count_recompute_step,
                      .use_log10 = argv_sentence_recompute_step_use_log10,
                      .recompute_step = argv_sentence_count_recompute_step,
                      .recompute_step_log10 = argv_sentence_count_recompute_step_log10,
                  },
              .document =
                  (struct measurement_step){
                      .enable_count_recompute_step = argv_enable_document_count_recompute_step,
                      .use_log10 = argv_document_recompute_step_use_log10,
                      .recompute_step = argv_document_count_recompute_step,
                      .recompute_step_log10 = argv_document_count_recompute_step_log10,
                  },
          },
  };

  err = measurement(&mcfg);
  if (err != 0) {
    perror("failed to call measurement\n");
    if (argv_force_timing_and_memory_to_output_path) {
      free(argv_output_path_timing);
      free(argv_output_path_memory);
    }
    return 1;
  }

  if (argv_force_timing_and_memory_to_output_path) {
    free(argv_output_path_timing);
    free(argv_output_path_memory);
  }

#if ENABLE_FILTER == 1
  filter_release();
#endif

  return 0;

malloc_fail:
  perror("malloc failed\n");
  return 1;
}
