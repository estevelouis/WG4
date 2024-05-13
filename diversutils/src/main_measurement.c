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

#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "graph.h"
#include "distributions.h"
#include "stats.h"
#include "dfunctions.h"
#include "sorted_array/array.h"

#include "jsonl/parser.h"

#define BFR_SIZE 256
#define MAX_FILES 1024

#ifndef ENABLE_STIRLING
#define ENABLE_STIRLING 1
#endif
#ifndef ENABLE_RICOTTA_SZEIDL
#define ENABLE_RICOTTA_SZEIDL 1
#endif
#ifndef ENABLE_PAIRWISE
#define ENABLE_PAIRWISE 1
#endif
#ifndef ENABLE_LEXICOGRAPHIC
#define ENABLE_LEXICOGRAPHIC 0
#endif
#ifndef ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY
#define ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY 1
#endif
#ifndef ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY
#define ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY 1
#endif
#ifndef ENABLE_LEINSTER_COBBOLD_DIVERSITY
#define ENABLE_LEINSTER_COBBOLD_DIVERSITY 0
#endif
#ifndef ENABLE_FUNCTIONAL_EVENNESS
#define ENABLE_FUNCTIONAL_EVENNESS 0
#endif
#ifndef ENABLE_FUNCTIONAL_DISPERSION
#define ENABLE_FUNCTIONAL_DISPERSION 1
#endif
#ifndef ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED
#define ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED 1
#endif

#ifndef ENABLE_NON_DISPARITY_FUNCTIONS
#define ENABLE_NON_DISPARITY_FUNCTIONS 1
#endif

#ifndef ENABLE_SHANNON_WEAVER_ENTROPY
#define ENABLE_SHANNON_WEAVER_ENTROPY 1
#endif
#ifndef ENABLE_GOOD_ENTROPY
#define ENABLE_GOOD_ENTROPY 1
#endif
#ifndef ENABLE_RENYI_ENTROPY
#define ENABLE_RENYI_ENTROPY 1
#endif
#ifndef ENABLE_PATIL_TAILLIE_ENTROPY
#define ENABLE_PATIL_TAILLIE_ENTROPY 1
#endif
#ifndef ENABLE_Q_LOGARITHMIC_ENTROPY
#define ENABLE_Q_LOGARITHMIC_ENTROPY 1
#endif
#ifndef ENABLE_SIMPSON_INDEX
#define ENABLE_SIMPSON_INDEX 1
#endif
#ifndef ENABLE_SIMPSON_DOMINANCE_INDEX
#define ENABLE_SIMPSON_DOMINANCE_INDEX 1
#endif
#ifndef ENABLE_HILL_NUMBER_STANDARD
#define ENABLE_HILL_NUMBER_STANDARD 1
#endif
#ifndef ENABLE_HILL_EVENNESS
#define ENABLE_HILL_EVENNESS 1
#endif
#ifndef ENABLE_BERGER_PARKER_INDEX
#define ENABLE_BERGER_PARKER_INDEX 1
#endif
#ifndef ENABLE_JUNGE1994_PAGE22
#define ENABLE_JUNGE1994_PAGE22 1
#endif
#ifndef ENABLE_BRILLOUIN_DIVERSITY
#define ENABLE_BRILLOUIN_DIVERSITY 1
#endif
#ifndef ENABLE_MCINTOSH_INDEX
#define ENABLE_MCINTOSH_INDEX 1
#endif
#ifndef ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975
#define ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975 1
#endif
#ifndef ENABLE_SW_E_HEIP
#define ENABLE_SW_E_HEIP 1
#endif
#ifndef ENABLE_SW_E_ONE_MINUS_D
#define ENABLE_SW_E_ONE_MINUS_D 1
#endif
#ifndef ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964
#define ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964 1
#endif
#ifndef ENABLE_SW_E_MINUS_LN_D_PIELOU1977
#define ENABLE_SW_E_MINUS_LN_D_PIELOU1977 1
#endif
#ifndef ENABLE_SW_F_2_1_ALATALO1981
#define ENABLE_SW_F_2_1_ALATALO1981 1
#endif
#ifndef ENABLE_SW_G_2_1_MOLINARI1989
#define ENABLE_SW_G_2_1_MOLINARI1989 1
#endif
#ifndef ENABLE_SW_E_BULLA1994
#define ENABLE_SW_E_BULLA1994 1
#endif
#ifndef ENABLE_SW_O_BULLA1994
#define ENABLE_SW_O_BULLA1994 1
#endif
#ifndef ENABLE_SW_E_MCI_PIELOU1969
#define ENABLE_SW_E_MCI_PIELOU1969 1
#endif
#ifndef ENABLE_SW_E_PRIME_CAMARGO1993
#define ENABLE_SW_E_PRIME_CAMARGO1993 1
#endif
#ifndef ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL
#define ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL 1
#endif


#define GOOD_ALPHA 2.0
#define GOOD_BETA 2.0
#define RENYI_ALPHA 2.0
#define PATIL_TAILLIE_ALPHA 1.0
#define Q_LOGARITHMIC_Q 2.0
#define HILL_NUMBER_STANDARD_ALPHA 2.0
#define HILL_EVENNESS_ALPHA 2.0
#define HILL_EVENNESS_BETA 1.0

#define STIRLING_ALPHA 1.0
#define STIRLING_BETA 2.0
#define RICOTTA_SZEIDL_ALPHA 2.0
#define CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA 2.0
#define SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA 2.0
#define LEINSTER_COBBOLD_DIVERSITY_ALPHA 2.0

#ifndef ENABLE_MULTITHREADED_MATRIX_GENERATION
#define ENABLE_MULTITHREADED_MATRIX_GENERATION 1
#endif

#ifndef ENABLE_TIMINGS
#define ENABLE_TIMINGS 1
#endif

#ifndef ENABLE_ITERATIVE_DISTANCE_COMPUTATION
#define ENABLE_ITERATIVE_DISTANCE_COMPUTATION 0
#endif
#ifndef ENABLE_MULTITHREADED_ROW_GENERATION
#define ENABLE_MULTITHREADED_ROW_GENERATION 1
#endif
#ifndef ROW_GENERATION_BATCH_SIZE
#define ROW_GENERATION_BATCH_SIZE 16
#endif

#ifndef SENTENCE_COUNT_RECOMPUTE_STEP
#define SENTENCE_COUNT_RECOMPUTE_STEP 1
#endif
#ifndef ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP
#define ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP 0
#endif
#ifndef SENTENCE_RECOMPUTE_STEP_USE_LOG10
#define SENTENCE_RECOMPUTE_STEP_USE_LOG10 0
#endif
#ifndef SENTENCE_COUNT_RECOMPUTE_STEP_LOG10
#define SENTENCE_COUNT_RECOMPUTE_STEP_LOG10 0.1
#endif

#ifndef DOCUMENT_COUNT_RECOMPUTE_STEP
#define DOCUMENT_COUNT_RECOMPUTE_STEP 1
#endif
#ifndef ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP
#define ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP 0
#endif
#ifndef DOCUMENT_RECOMPUTE_STEP_USE_LOG10
#define DOCUMENT_RECOMPUTE_STEP_USE_LOG10 0
#endif
#ifndef DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10
#define DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10 0.1
#endif

#ifndef INPUT_PATH
#define INPUT_PATH "measurement_files.txt"
#endif
#ifndef OUTPUT_PATH
#define OUTPUT_PATH "measurement_output.tsv"
#endif
#ifndef OUTPUT_PATH_TIMING
#define OUTPUT_PATH_TIMING "measurement_output_timing.tsv"
#endif
#ifndef OUTPUT_PATH_MEMORY
#define OUTPUT_PATH_MEMORY "measurement_output_memory.tsv"
#endif

#ifndef ENABLE_OUTPUT_TIMING
#define ENABLE_OUTPUT_TIMING 1
#endif
#ifndef ENABLE_OUTPUT_MEMORY
#define ENABLE_OUTPUT_MEMORY 1
#endif


#ifndef TARGET_COLUMN
#define TARGET_COLUMN UD_MWE
#endif

#ifndef W2V_PATH
#define W2V_PATH "/home/esteve/Documents/thesis/other_repos/word2vec/bin/MWE_S2S_IT_11GB_100d_skip-gram.bin"
#endif

#ifndef JSONL_CONTENT_KEY
#define JSONL_CONTENT_KEY "text"
#endif

const uint8_t ENABLE_DISTANCE_COMPUTATION = ENABLE_DISPARITY_FUNCTIONS && (ENABLE_STIRLING || ENABLE_RICOTTA_SZEIDL || ENABLE_PAIRWISE || ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY || ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY || ENABLE_LEINSTER_COBBOLD_DIVERSITY || ENABLE_LEXICOGRAPHIC || ENABLE_FUNCTIONAL_EVENNESS || ENABLE_FUNCTIONAL_DISPERSION || ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED);

enum {
	CONLLU_COLUMN_UPOS,
	CONLLU_COLUMN_DEPREL
};
const int32_t CONLLU_COLUMNS_TO_ADD[2] = {CONLLU_COLUMN_UPOS, CONLLU_COLUMN_DEPREL};
// const int32_t NUM_CONLLU_COLUMNS_TO_ADD = (int32_t) (sizeof(CONLLU_COLUMNS_TO_ADD) / sizeof(int32_t));
const int32_t NUM_CONLLU_COLUMNS_TO_ADD = 0;

const int32_t CONLLU_ADD_FORM = 0;

double stacked_sentence_count_log10 = SENTENCE_COUNT_RECOMPUTE_STEP_LOG10;
uint64_t stacked_sentence_count_target;
double stacked_document_count_log10 = DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10;
uint64_t stacked_document_count_target;

static int32_t time_ns_delta(int64_t* const delta){
	static struct timespec static_ts = (struct timespec) {};

	struct timespec new_ts;

	if(clock_gettime(CLOCK_MONOTONIC, &new_ts) != 0){
		perror("Failed to call clock_gettime\n");
		return 1;
	}

	if(delta != NULL){
		(*delta) = (new_ts.tv_sec - static_ts.tv_sec) * 1000000000 + (new_ts.tv_nsec - static_ts.tv_nsec);
	}

	static_ts = new_ts;

	return 0;
}

int32_t virtual_memory_consumption(int64_t* const res){
	FILE* f;
	const int32_t bfr_size = 2048;
	char bfr[bfr_size];
	char* strtok_pointer;

	if(res == NULL){
		perror("Cannot pass NULL pointer to virtual_memory_consumption\n");
		return 1;
	}

	memset(bfr, '\0', bfr_size * sizeof(char));

	f = fopen("/proc/self/stat", "r");

	if(f == NULL){
		perror("Failed to open /proc/self/stat\n");
		return 1;
	}

	if(fgets(bfr, bfr_size, f) == 0){
		perror("No bytes could be read from /proc/self/stat\n");
		fclose(f);
		return 1;
	}

	strtok_pointer = strtok(bfr, " ");
	for(int32_t i = 0 ; i < 23 ; i++){
		strtok_pointer = strtok(NULL, " ");
	}
	(*res) = strtol(strtok_pointer, NULL, 10);

	fclose(f);

	return 0;
}

void timing_and_memory(FILE* f_timing_ptr, FILE* f_memory_ptr, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	int64_t ns_delta, virtual_mem;
	if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){perror("Failed to call time_ns_delta\n"); exit(1);} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
	if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){perror("Failed to call virtual_memory_consumption\n"); exit(1);} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
}

int32_t wrap_diversity_1r_0a(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double* const, const int8_t, const struct matrix* const), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, fp_mode, m) != 0){perror("Failed to call diversity function in wrap_diversity_1r_0a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e", res1);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}

int32_t wrap_diversity_1r_0a_no_matrix(struct graph* const g, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double* const, const int8_t), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, fp_mode) != 0){perror("Failed to call diversity function in wrap_diversity_1r_0a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e", res1);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}

int32_t wrap_diversity_2r_1a(struct graph* const g, struct matrix* const m, const double alpha, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double* const, double* const, double, const int8_t, const struct matrix* const), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1, res2;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, &res2, alpha, fp_mode, m) != 0){perror("Failed to call diversity function in wrap_diversity_2r_1a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e\t%.10e", res1, res2);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}

int32_t wrap_diversity_2r_0a(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double* const, double* const, const int8_t, const struct matrix* const), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1, res2;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, &res2, fp_mode, m) != 0){perror("Failed to call diversity function in wrap_diversity_2r_0a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e\t%.10e", res1, res2);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}

int32_t wrap_diversity_2r_0a_long_double_alt(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double* const, long double* const, const int8_t, const struct matrix* const), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1;
	long double res2;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, &res2, fp_mode, m) != 0){perror("Failed to call diversity function in wrap_diversity_2r_0a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(enable_timings){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e\t%.10Le", res1, res2);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}

/*
int32_t wrap_diversity_1r_0a(struct graph* const g, struct matrix* const m, const int8_t fp_mode, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int32_t (*df)(struct graph* const, double*, const int8_t, struct matrix* const), const uint8_t enable_timings, const uint8_t enable_output_timing, const uint8_t enable_output_memory){
	double res1;
	time_t t, delta_t;
	t = time(NULL);
	if(df(g, &res1, fp_mode, m) != 0){perror("Failed to call diversity function in wrap_diversity_1r_0a\n"); return EXIT_FAILURE;}
	delta_t = time(NULL) - t;
	if(ENABLE_TIMINGS){printf("[log] [time] Computed df in %lis\n", delta_t);}
	fprintf(f_ptr, "\t%.10e", res1);
	timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
	return 0;
}
*/

int32_t apply_diversity_functions_to_graph(struct graph* g, struct minimum_spanning_tree* mst, struct graph_distance_heap* heap, FILE* f_ptr, FILE* f_timing_ptr, FILE* f_memory_ptr, int64_t* previous_g_num_nodes_p, int64_t* num_sentences_p, int64_t* num_all_sentences_p, double* best_s, int8_t* mst_initialised, uint64_t i, struct sorted_array* sorted_array_discarded_because_not_in_vector_database,
// 	uint32_t target_column,
	uint32_t num_row_threads,
	uint32_t num_matrix_threads,
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
// 	uint8_t enable_sentence_count_recompute_step,
// 	uint8_t sentence_recompute_step_use_log10,
// 	uint8_t enable_document_count_recompute_step,
// 	uint8_t document_recompute_step_use_log10,
 	uint8_t enable_output_timing,
 	uint8_t enable_output_memory,
// 	uint64_t sentence_count_recompute_step,
// 	double sentence_count_recompute_step_log10,
// 	uint64_t document_count_recompute_step,
// 	double document_count_recompute_step_log10,
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
){
	const uint8_t enable_distance_computation = enable_disparity_functions && (enable_stirling || enable_ricotta_szeidl || enable_pairwise || enable_chao_et_al_functional_diversity || enable_scheiner_species_phylogenetic_functional_diversity || enable_leinster_cobbold_diversity || enable_lexicographic || enable_functional_evenness || enable_functional_dispersion || enable_functional_divergence_modified);

	int32_t err;
	if(enable_iterative_distance_computation){
		size_t local_malloc_size = g->num_nodes * sizeof(float) * row_generation_batch_size;
		float* vector_batch = (float*) malloc(local_malloc_size);
		if(vector_batch == NULL){goto malloc_fail;}
		memset(vector_batch, '\0', local_malloc_size);

		local_malloc_size = g->num_nodes * sizeof(uint8_t) * row_generation_batch_size;
		uint8_t* used = (uint8_t*) malloc(local_malloc_size);
		if(used == NULL){goto malloc_fail;}
		memset(used, '\0', local_malloc_size);

		struct iterative_state_stirling_from_graph iter_state_stirling;
		struct iterative_state_pairwise_from_graph iter_state_pairwise;
		if(create_iterative_state_stirling_from_graph(&iter_state_stirling, g, stirling_alpha, stirling_beta) != 0){
			perror("failled to call create_iterateive_state_stirling_from_graph\n");
			return 1;
		}
		if(create_iterative_state_pairwise_from_graph(&iter_state_pairwise, g) != 0){
			perror("failled to call create_iterateive_state_pairwise_from_graph\n");
			return 1;
		}

		int32_t i_index = 0;
		double sum = 0.0;
		for(uint64_t h = 0 ; h < g->num_nodes ; h += row_generation_batch_size){
			if(enable_multithreaded_row_generation){
				distance_row_batch_from_graph_multithread(g, i_index, vector_batch, num_row_threads, row_generation_batch_size);
			} else {
				perror("single-threaded version of row computation has been disabled\n");
				return 1;
			}

			/* // DO NOT REMOVE
			// for(uint64_t m = h ; m < g->num_nodes ; m++){
			for(uint64_t m = 0 ; m < row_generation_batch_size && h + m < g->num_nodes ; m++){
				for(uint64_t p = i_index + 1 ; p < g->num_nodes ; p++){
					sum += (double) vector_batch[m * g->num_nodes + p];
				}
	
				if(enable_stirling){iterate_iterative_state_stirling_from_graph(&iter_state_stirling, &(vector_batch[m * g->num_nodes]));}
				if(enable_pairwise){iterate_iterative_state_pairwise_from_graph(&iter_state_pairwise, &(vector_batch[m * g->num_nodes]));}
	
				i_index = h + m + 1;
				iter_state_stirling.i = i_index;
				iter_state_pairwise.i = i_index;
			}
			*/

			uint64_t actual_row_generation_batch_size = row_generation_batch_size;
			if(h + row_generation_batch_size >= g->num_nodes){
				actual_row_generation_batch_size = g->num_nodes - h;
			}

			pthread_t agg_threads[row_generation_batch_size]; // still have full buffer to make it writeable?
			struct thread_args_aggregator agg_thread_args[row_generation_batch_size];
			if(enable_pairwise){
				for(uint64_t m = 0 ; m < actual_row_generation_batch_size ; m++){
					struct thread_args_aggregator local_args = (struct thread_args_aggregator) {
						.iter_state.pairwise = &iter_state_pairwise,
						.i = i_index,
						.vector = &(vector_batch[m * g->num_nodes]),
					};
					memcpy(&(agg_thread_args[m]), &local_args, sizeof(struct thread_args_aggregator));

					if(pthread_create(&(agg_threads[m]), NULL, iterate_iterative_state_pairwise_from_graph_avx_thread, &(agg_thread_args[m])) != 0){
						perror("Failed to call pthread_create\n");
						return 1;
					}

					i_index = h + m + 1; // !
					iter_state_stirling.i = i_index; // !
					iter_state_pairwise.i = i_index; // !
				}

				for(uint64_t m = 0 ; m < actual_row_generation_batch_size ; m++){
					pthread_join(agg_threads[m], NULL);
				}
			}
		}
		free(vector_batch);
		free(used);


		// finalise_iterative_state_stirling_from_graph(&iter_state_stirling);
		finalise_iterative_state_pairwise_from_graph(&iter_state_pairwise);

		double mu_dist = sum / ((double) (g->num_nodes * (g->num_nodes - 1) / 2));

		fprintf(f_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%c", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes, mu_dist, '?'); // recomputing sigma dist would be expensive

		if(enable_stirling){printf("[log] [end iter] stirling: %f\n", iter_state_stirling.result);}
		fprintf(f_ptr, "\t%.10e", iter_state_stirling.result);
		if(enable_pairwise){printf("[log] [end iter] pairwise: %f\n", iter_state_pairwise.result);}
		fprintf(f_ptr, "\t%.10e", iter_state_pairwise.result);

		fprintf(f_ptr, "\n");
	} else {
		int64_t ns_delta, virtual_mem;

		if(enable_output_timing){
			fprintf(f_timing_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes);
			if(time_ns_delta(NULL) != 0){goto time_ns_delta_failure;}
		}
		if(enable_output_memory){
			fprintf(f_memory_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes);
		}

		struct matrix m_mst;
		if(enable_functional_evenness){
			if(create_matrix(&m_mst, g->num_nodes, g->num_nodes, FP64) != 0){
				perror("failed to call create_matrix for MST\n");
				return 1;
			}
			memset(m_mst.active, '\0', g->num_nodes * g->num_nodes * sizeof(uint8_t));
			memset(m_mst.active_final, '\0', g->num_nodes * g->num_nodes * sizeof(uint8_t));
			if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
			if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
		}

		time_t t, delta_t;

		struct matrix m;
		if(enable_distance_computation){
			if(create_matrix(&m, (uint32_t) g->num_nodes, (uint32_t) g->num_nodes, FP32) != 0){
				perror("failed to call create_matrix\n");
				return 1;
			}
			t = time(NULL);
			if(enable_multithreaded_matrix_generation){
				if(distance_matrix_from_graph_multithread(g, &m, num_matrix_threads) != 0){
					perror("failed to call distance_matrix_from_graph_multithread\n");
					return 1;
				}
			} else {
				// if(distance_matrix_from_graph(g, ANGULAR_MINKOWSKI_DISTANCE_ORDER, &m) != 0){
				if(distance_matrix_from_graph(g, &m) != 0){
					perror("failed to call distance_matrix_from_graph\n");
					return 1;
				}
			}
			if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
			if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			delta_t = time(NULL) - t;
			if(enable_timings){
				printf("[log] [time] Computed matrix in %lis\n", delta_t);
			}
		}

		if((*mst_initialised)){
			free_graph_distance_heap(heap);
			if(enable_functional_evenness){
				free_minimum_spanning_tree(mst);
			}
		}

		struct graph_distance_heap local_heap;

		double mu_dist = NAN;
		double sigma_dist = NAN;
		float mu_dist_fp32 = NAN;
		float sigma_dist_fp32 = NAN;
		if(enable_functional_evenness){
			err = create_graph_distance_heap(&local_heap, g, GRAPH_NODE_FP32, &m);
			if(err != 0){
				perror("failed to call create_graph_distance_heap\n");
				return EXIT_FAILURE;
			}
	
			(*heap) = local_heap;

			if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
			if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
		}

		if(enable_distance_computation){
			switch(m.fp_mode){
				case FP32:
					avg_and_std_fp32(m.bfr.fp32, m.a * m.b, &mu_dist_fp32, &sigma_dist_fp32);
					mu_dist = (double) mu_dist_fp32;
					sigma_dist = (double) sigma_dist_fp32;
					break;
				case FP64:
					avg_and_std_fp64(m.bfr.fp64, m.a * m.b, &mu_dist, &sigma_dist);
					break;
				default:
					perror("unknown FP mode\n");
					return 1;
			}
			if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
			if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
		}

		fprintf(f_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%.10e", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes, mu_dist, sigma_dist);

		if(enable_disparity_functions){
			if(enable_output_timing){if(time_ns_delta(NULL) != 0){goto time_ns_delta_failure;}}
	
			if(enable_functional_evenness){
				t = time(NULL);
				struct minimum_spanning_tree local_mst;
	
				err = create_minimum_spanning_tree(&local_mst, &local_heap);
				if(err != 0){
					perror("failed to call create_minimum_spanning_tree\n");
					return EXIT_FAILURE;
				}
	
				err = calculate_minimum_spanning_tree(&local_mst, &m_mst, MST_PRIMS_ALGORITHM);
				if(err != 0){
					perror("failed to call calculate_minimum_spanning_tree\n");
					return EXIT_FAILURE;	
				}
				
				(*mst) = local_mst;
				(*mst_initialised) = 1;
	
				delta_t = time(NULL) - t;
				if(enable_timings){
					printf("[log] [time] Computed MST in %lis\n", delta_t);
				}
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			(*previous_g_num_nodes_p) = g->num_nodes;
	
			if(enable_stirling){
				double stirling;
				t = time(NULL);
				err = stirling_from_graph(g, &stirling, stirling_alpha, stirling_beta, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call stirling_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(enable_timings){
					printf("[log] [time] Computed Stirling in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", stirling);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
	
			if(enable_ricotta_szeidl){
				double ricotta_szeidl;
				t = time(NULL);
				err = ricotta_szeidl_from_graph(g, &ricotta_szeidl, ricotta_szeidl_alpha, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call ricotta_szeidl_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(enable_timings){
					printf("[log] [time] Computed Ricotta-Szeidl in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", ricotta_szeidl);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
	
			if(enable_pairwise){
				if(wrap_diversity_1r_0a(g, &m, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, pairwise_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_chao_et_al_functional_diversity){
				if(wrap_diversity_2r_1a(g, &m, chao_et_al_functional_diversity_alpha, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, chao_et_al_functional_diversity_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_scheiner_species_phylogenetic_functional_diversity){
				if(wrap_diversity_2r_1a(g, &m, scheiner_species_phylogenetic_functional_diversity_alpha, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, scheiner_species_phylogenetic_functional_diversity_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_leinster_cobbold_diversity){
				if(wrap_diversity_2r_1a(g, &m, leinster_cobbold_diversity_alpha, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, leinster_cobbold_diversity_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_lexicographic){
				if(wrap_diversity_2r_0a_long_double_alt(g, &m, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, lexicographic_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_functional_evenness){
				if(!(*mst_initialised)){
					perror("mst not initialised\n");
					return 1;
				}
				double functional_evenness;
				t = time(NULL);
				err = functional_evenness_from_minimum_spanning_tree(mst, &functional_evenness);
				if(err != 0){
					perror("failed to call functional_evenness_from_minimum_spanning_tree\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(enable_timings){
					printf("[log] [time] Computed FEve in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", functional_evenness);
				timing_and_memory(f_timing_ptr, f_memory_ptr, enable_output_timing, enable_output_memory);
			}
	
			if(enable_functional_dispersion){
				if(wrap_diversity_1r_0a_no_matrix(g, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, functional_dispersion_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
	
			if(enable_functional_divergence_modified){
				if(wrap_diversity_1r_0a_no_matrix(g, GRAPH_NODE_FP32, f_ptr, f_timing_ptr, f_memory_ptr, functional_divergence_modified_from_graph, enable_timings, enable_output_timing, enable_output_memory) != 0){return 1;}
			}
		}

		if(enable_non_disparity_functions){
			if(enable_shannon_weaver_entropy){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				shannon_weaver_entropy_from_graph(g, &res_entropy, &res_hill_number);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed SW entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_good_entropy){
				double res;
				time_t t = time(NULL);
				good_entropy_from_graph(g, &res, good_alpha, good_beta);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Good entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_renyi_entropy){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				renyi_entropy_from_graph(g, &res_entropy, &res_hill_number, renyi_alpha);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Renyi entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_patil_taillie_entropy){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				patil_taillie_entropy_from_graph(g, &res_entropy, &res_hill_number, patil_taillie_alpha);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Patil-Taillie entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_q_logarithmic_entropy){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				q_logarithmic_entropy_from_graph(g, &res_entropy, &res_hill_number, q_logarithmic_q);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed q-logarithmic entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_simpson_index){
				double res;
				time_t t = time(NULL);
				simpson_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Simpson index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_simpson_dominance_index){
				double res;
				time_t t = time(NULL);
				simpson_dominance_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Simpson dominance index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_hill_number_standard){
				double res;
				time_t t = time(NULL);
				hill_number_standard_from_graph(g, &res, hill_number_standard_alpha);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Hill number (standard) in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_hill_evenness){
				double res;
				time_t t = time(NULL);
				hill_evenness_from_graph(g, &res, hill_evenness_alpha, hill_evenness_beta);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Hill evenness in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_berger_parker_index){
				double res;
				time_t t = time(NULL);
				berger_parker_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Berger Parker index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_junge1994_page22){
				double res;
				time_t t = time(NULL);
				junge1994_page22_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Junge 1994 p22 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_brillouin_diversity){
				double res;
				time_t t = time(NULL);
				brillouin_diversity_from_graph(g, &res);
				// printf("res: %f\n", res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed Brillouin diversity in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_mcintosh_index){
				double res;
				time_t t = time(NULL);
				mcintosh_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed McIntosh index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_entropy_over_log_n_species_pielou1975){
				double res;
				time_t t = time(NULL);
				sw_entropy_over_log_n_species_pielou1975_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) entropy over log n species Pielou 1975 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_heip){
				double res;
				time_t t = time(NULL);
				sw_e_heip_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E Heip in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_one_minus_d){
				double res;
				time_t t = time(NULL);
				sw_e_one_minus_D_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E one minus D in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_one_over_ln_d_williams1964){
				double res;
				time_t t = time(NULL);
				sw_e_one_over_D_williams1964_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E one over ln D Williams 1964 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_minus_ln_d_pielou1977){
				double res;
				time_t t = time(NULL);
				sw_e_minus_ln_D_pielou1977_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E minus ln D Pielou 1977 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_f_2_1_alatalo1981){
				double res;
				time_t t = time(NULL);
				sw_f_2_1_alatalo1981_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) F_2_1 Alatalo 1981 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_g_2_1_molinari1989){
				double res;
				time_t t = time(NULL);
				sw_g_2_1_molinari1989_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) G_2_1 Molinari 1989 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_bulla1994){
				double res;
				time_t t = time(NULL);
				sw_e_bulla1994_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E Bulla 1994 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_o_bulla1994){
				double res;
				time_t t = time(NULL);
				sw_o_bulla1994_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) O bulla 1994 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_mci_pielou1969){
				double res;
				time_t t = time(NULL);
				sw_e_mci_pielou1969_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E MCI Pielou 1969 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_prime_camargo1993){
				double res;
				time_t t = time(NULL);
				sw_e_prime_camargo1993_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E prime Camargo 1993 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
			if(enable_sw_e_var_smith_and_wilson1996_original){
				double res;
				time_t t = time(NULL);
				sw_e_var_smith_and_wilson1996_original_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(enable_timings){printf("[log] [time] Computed (SW) E var Smith and Wilson 1996 original in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
				if(enable_output_timing){if(time_ns_delta(&ns_delta) != 0){goto time_ns_delta_failure;} else {fprintf(f_timing_ptr, "\t%li", ns_delta);}}
				if(enable_output_memory){if(virtual_memory_consumption(&virtual_mem) != 0){goto virtual_memory_consumption_failure;} else {fprintf(f_memory_ptr, "\t%li", virtual_mem);}}
			}
		}

		fprintf(f_ptr, "\n");
		if(enable_output_timing){fprintf(f_timing_ptr, "\n");}
		if(enable_output_memory){fprintf(f_memory_ptr, "\n");}

		if(enable_distance_computation){
			free_matrix(&m);
			if(enable_functional_evenness){
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

int32_t measurement(
	struct word2vec* w2v,
	char* jsonl_content_key,
	char* input_path,
	char* output_path,
	char* output_path_timing,
	char* output_path_memory,
	uint32_t target_column,
	uint32_t num_row_threads,
	uint32_t num_matrix_threads,
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
	uint8_t enable_sentence_count_recompute_step,
	uint8_t sentence_recompute_step_use_log10,
	uint8_t enable_document_count_recompute_step,
	uint8_t document_recompute_step_use_log10,
	uint8_t enable_output_timing,
	uint8_t enable_output_memory,
	uint64_t sentence_count_recompute_step,
	double sentence_count_recompute_step_log10,
	uint64_t document_count_recompute_step,
	double document_count_recompute_step_log10,
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
){
	const uint8_t enable_distance_computation = enable_disparity_functions && (enable_stirling || enable_ricotta_szeidl || enable_pairwise || enable_chao_et_al_functional_diversity || enable_scheiner_species_phylogenetic_functional_diversity || enable_leinster_cobbold_diversity || enable_lexicographic || enable_functional_evenness || enable_functional_dispersion || enable_functional_divergence_modified);

	stacked_sentence_count_target = log(stacked_sentence_count_log10) / log(10.0);
	stacked_document_count_target = log(stacked_document_count_log10) / log(10.0);

	struct sorted_array sorted_array_discarded_because_not_in_vector_database;
	if(create_sorted_array(&sorted_array_discarded_because_not_in_vector_database, 0, sizeof(struct sorted_array_str_int_element), sorted_array_str_int_cmp) != 0){
		perror("failed to call create_sorted_array\n");
		return 1;
	}

	struct graph g;
	if(create_graph_empty(&g) != 0){
		perror("failed to call create_graph_empty\n");
		return 1;
	}
	struct graph_distance_heap heap;
	struct minimum_spanning_tree mst;

	memset(&heap, '\0', sizeof(struct graph_distance_heap));
	memset(&mst, '\0', sizeof(struct minimum_spanning_tree));

	int32_t zipfian_n = 0;

	double previous_best_s = -1.0;

	for(uint64_t i = 0 ; i < w2v->num_vectors ; i++){
		w2v->keys[i].active_in_current_graph = 0;
		w2v->keys[i].num_occurrences = 0;
		w2v->keys[i].graph_node_pointer = NULL; // !
	}

	FILE* f_paths_ptr = fopen(input_path, "r");
	if(f_paths_ptr == NULL){
		fprintf(stderr, "cannot open %s\n", input_path);
		return 1;
	}
	char* bfr[MAX_FILES];
	char bfr_read[BFR_SIZE];
	memset(bfr_read, '\0', BFR_SIZE);

	int32_t num_files = 0;
	while(fgets(bfr_read, BFR_SIZE - 1, f_paths_ptr)){
		if(bfr_read[0] == '#'){
			continue;
		}

		size_t bytes_to_cpy = strlen(bfr_read);
		if(bfr_read[bytes_to_cpy - 1] == '\n'){
			bytes_to_cpy--;
		}
		void* local_malloc_p = malloc(bytes_to_cpy + 1);
		if(local_malloc_p == NULL){
			perror("malloc failed\n");
			return 1;
		}
		memset(local_malloc_p, '\0', bytes_to_cpy + 1);
		memcpy(local_malloc_p, bfr_read, bytes_to_cpy);
		bfr[num_files] = (char*) local_malloc_p;

		num_files++;
	}
	int32_t num_input_paths = num_files;
	char** input_paths = bfr;
	char** input_paths_true_positives = NULL;

	FILE* f_ptr = fopen(output_path, "w");
	if(f_ptr == NULL){fprintf(stderr, "Failed to open file: %s\n", output_path); return EXIT_FAILURE;}
	fprintf(f_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");

	if(enable_disparity_functions){
		if(enable_stirling){fprintf(f_ptr, "\tstirling_alpha%.10e_beta%.10e", stirling_alpha, stirling_beta);}
		if(!enable_iterative_distance_computation && enable_ricotta_szeidl){fprintf(f_ptr, "\tricotta_szeidl_alpha%.10e", ricotta_szeidl_alpha);}
		if(enable_pairwise){fprintf(f_ptr, "\tpairwise");}
		if(!enable_iterative_distance_computation && enable_chao_et_al_functional_diversity){fprintf(f_ptr, "\tchao_et_al_functional_diversity_alpha%.10e\tchao_et_al_functional_hill_number_alpha%.10e", chao_et_al_functional_diversity_alpha, chao_et_al_functional_diversity_alpha);}
		if(!enable_iterative_distance_computation && enable_scheiner_species_phylogenetic_functional_diversity){fprintf(f_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e\tscheiner_species_phylogenetic_functional_hill_number_alpha%.10e", scheiner_species_phylogenetic_functional_diversity_alpha, scheiner_species_phylogenetic_functional_diversity_alpha);}
		if(!enable_iterative_distance_computation && enable_leinster_cobbold_diversity){fprintf(f_ptr, "\tleinster_cobbold_diversity_alpha%.10e\tleinster_cobbold_hill_number_alpha%.10e", leinster_cobbold_diversity_alpha, leinster_cobbold_diversity_alpha);}
		if(!enable_iterative_distance_computation && enable_lexicographic){fprintf(f_ptr, "\tlexicographic\tlexicographic_hybrid_scheiner");}
		if(!enable_iterative_distance_computation && enable_functional_evenness){fprintf(f_ptr, "\tfunctional_evenness");}
		if(!enable_iterative_distance_computation && enable_functional_dispersion){fprintf(f_ptr, "\tfunctional_dispersion");}
		if(!enable_iterative_distance_computation && enable_functional_divergence_modified){fprintf(f_ptr, "\tfunctional_divergence_modified");}
	}
	if(enable_non_disparity_functions){
		if(enable_shannon_weaver_entropy){fprintf(f_ptr, "\tshannon_weaver_entropy\tshannon_weaver_hill_number");}
		if(enable_good_entropy){fprintf(f_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", good_alpha, good_beta);}
		if(enable_renyi_entropy){fprintf(f_ptr, "\trenyi_entropy_alpha%.4e\trenyi_hill_number_alpha%.4e", renyi_alpha, renyi_alpha);}
		if(enable_patil_taillie_entropy){fprintf(f_ptr, "\tpatil_taillie_entropy_alpha%.4e\tpatil_taillie_hill_number_alpha%.4e", patil_taillie_alpha, patil_taillie_alpha);}
		if(enable_q_logarithmic_entropy){fprintf(f_ptr, "\tq_logarithmic_entropy_alpha%.4e\tq_logarithmic_hill_number_alpha%.4e", q_logarithmic_q, q_logarithmic_q);}
		if(enable_simpson_index){fprintf(f_ptr, "\tsimpson_index");}
		if(enable_simpson_dominance_index){fprintf(f_ptr, "\tsimpson_dominance_index");}
		if(enable_hill_number_standard){fprintf(f_ptr, "\thill_number_standard_alpha%.4e", hill_number_standard_alpha);}
		if(enable_hill_evenness){fprintf(f_ptr, "\thill_evenness_alpha%.4e_beta%.4e", hill_evenness_alpha, hill_evenness_beta);}
		if(enable_berger_parker_index){fprintf(f_ptr, "\tberger_parker_index");}
		if(enable_junge1994_page22){fprintf(f_ptr, "\tjunge1994_page22");}
		if(enable_brillouin_diversity){fprintf(f_ptr, "\tbrillouin_diversity");}
		if(enable_mcintosh_index){fprintf(f_ptr, "\tmcintosh_index");}
		if(enable_sw_entropy_over_log_n_species_pielou1975){fprintf(f_ptr, "\tsw_entropy_over_log_n_species_pielou1975");}
		if(enable_sw_e_heip){fprintf(f_ptr, "\tsw_e_heip");}
		if(enable_sw_e_one_minus_d){fprintf(f_ptr, "\te_one_minus_D");}
		if(enable_sw_e_one_over_ln_d_williams1964){fprintf(f_ptr, "\te_one_over_ln_D_williams1964");}
		if(enable_sw_e_minus_ln_d_pielou1977){fprintf(f_ptr, "\te_minus_ln_D_pielou1977");}
		if(enable_sw_f_2_1_alatalo1981){fprintf(f_ptr, "\tsw_f_2_1_alatalo1981");}
		if(enable_sw_g_2_1_molinari1989){fprintf(f_ptr, "\tsw_g_2_1_molinari1989");}
		if(enable_sw_e_bulla1994){fprintf(f_ptr, "\tsw_e_bulla1994");}
		if(enable_sw_o_bulla1994){fprintf(f_ptr, "\tsw_o_bulla1994");}
		if(enable_sw_e_mci_pielou1969){fprintf(f_ptr, "\tsw_e_mci_pielou1969");}
		if(enable_sw_e_prime_camargo1993){fprintf(f_ptr, "\tsw_e_prime_camargo1993");}
		if(enable_sw_e_var_smith_and_wilson1996_original){fprintf(f_ptr, "\tsw_e_var_smith_and_wilson1996_original");}
	}
	fprintf(f_ptr, "\n");

	FILE* f_timing_ptr = NULL;
	if(enable_output_timing){
		f_timing_ptr = fopen(output_path_timing, "w");
		if(f_timing_ptr == NULL){fprintf(stderr, "Failed to open file: %s\n", output_path_timing); return EXIT_FAILURE;}
		// fprintf(f_timing_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
		fprintf(f_timing_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn");

		if(enable_disparity_functions){
			// ----
			if(enable_functional_evenness){fprintf(f_timing_ptr, "\tm_mst_creation");}
			if(enable_distance_computation){fprintf(f_timing_ptr, "\tdist_matrix_computation");}
			if(enable_functional_evenness){fprintf(f_timing_ptr, "\tdist_heap_computation");}
			if(enable_distance_computation){fprintf(f_timing_ptr, "\tdist_stat_computation");}
			if(enable_functional_evenness){fprintf(f_timing_ptr, "\tmst_computation");}
			// ----

			if(enable_stirling){fprintf(f_timing_ptr, "\tstirling_alpha%.10e_beta%.10e", stirling_alpha, stirling_beta);}
			if(!enable_iterative_distance_computation && enable_ricotta_szeidl){fprintf(f_timing_ptr, "\tricotta_szeidl_alpha%.10e", ricotta_szeidl_alpha);}
			if(enable_pairwise){fprintf(f_timing_ptr, "\tpairwise");}
			if(!enable_iterative_distance_computation && enable_chao_et_al_functional_diversity){fprintf(f_timing_ptr, "\tchao_et_al_functional_diversity_alpha%.10e", chao_et_al_functional_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_scheiner_species_phylogenetic_functional_diversity){fprintf(f_timing_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e", scheiner_species_phylogenetic_functional_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_leinster_cobbold_diversity){fprintf(f_timing_ptr, "\tleinster_cobbold_diversity_alpha%.10e", leinster_cobbold_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_lexicographic){fprintf(f_timing_ptr, "\tlexicographic");}
			if(!enable_iterative_distance_computation && enable_functional_evenness){fprintf(f_timing_ptr, "\tfunctional_evenness");}
			if(!enable_iterative_distance_computation && enable_functional_dispersion){fprintf(f_timing_ptr, "\tfunctional_dispersion");}
			if(!enable_iterative_distance_computation && enable_functional_divergence_modified){fprintf(f_timing_ptr, "\tfunctional_divergence_modified");}
		}
		if(enable_non_disparity_functions){
			if(enable_shannon_weaver_entropy){fprintf(f_timing_ptr, "\tshannon_weaver_entropy");}
			if(enable_good_entropy){fprintf(f_timing_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", good_alpha, good_beta);}
			if(enable_renyi_entropy){fprintf(f_timing_ptr, "\trenyi_entropy_alpha%.4e", renyi_alpha);}
			if(enable_patil_taillie_entropy){fprintf(f_timing_ptr, "\tpatil_taillie_entropy_alpha%.4e", patil_taillie_alpha);}
			if(enable_q_logarithmic_entropy){fprintf(f_timing_ptr, "\tq_logarithmic_entropy_alpha%.4e", q_logarithmic_q);}
			if(enable_simpson_index){fprintf(f_timing_ptr, "\tsimpson_index");}
			if(enable_simpson_dominance_index){fprintf(f_timing_ptr, "\tsimpson_dominance_index");}
			if(enable_hill_number_standard){fprintf(f_timing_ptr, "\thill_number_standard_alpha%.4e", hill_number_standard_alpha);}
			if(enable_hill_evenness){fprintf(f_timing_ptr, "\thill_evenness_alpha%.4e_beta%.4e", hill_evenness_alpha, hill_evenness_beta);}
			if(enable_berger_parker_index){fprintf(f_timing_ptr, "\tberger_parker_index");}
			if(enable_junge1994_page22){fprintf(f_timing_ptr, "\tjunge1994_page22");}
			if(enable_brillouin_diversity){fprintf(f_timing_ptr, "\tbrillouin_diversity");}
			if(enable_mcintosh_index){fprintf(f_timing_ptr, "\tmcintosh_index");}
			if(enable_sw_entropy_over_log_n_species_pielou1975){fprintf(f_timing_ptr, "\tsw_entropy_over_log_n_species_pielou1975");}
			if(enable_sw_e_heip){fprintf(f_timing_ptr, "\tsw_e_heip");}
			if(enable_sw_e_one_minus_d){fprintf(f_timing_ptr, "\te_one_minus_D");}
			if(enable_sw_e_one_over_ln_d_williams1964){fprintf(f_timing_ptr, "\te_one_over_ln_D_williams1964");}
			if(enable_sw_e_minus_ln_d_pielou1977){fprintf(f_timing_ptr, "\te_minus_ln_D_pielou1977");}
			if(enable_sw_f_2_1_alatalo1981){fprintf(f_timing_ptr, "\tsw_f_2_1_alatalo1981");}
			if(enable_sw_g_2_1_molinari1989){fprintf(f_timing_ptr, "\tsw_g_2_1_molinari1989");}
			if(enable_sw_e_bulla1994){fprintf(f_timing_ptr, "\tsw_e_bulla1994");}
			if(enable_sw_o_bulla1994){fprintf(f_timing_ptr, "\tsw_o_bulla1994");}
			if(enable_sw_e_mci_pielou1969){fprintf(f_timing_ptr, "\tsw_e_mci_pielou1969");}
			if(enable_sw_e_prime_camargo1993){fprintf(f_timing_ptr, "\tsw_e_prime_camargo1993");}
			if(enable_sw_e_var_smith_and_wilson1996_original){fprintf(f_timing_ptr, "\tsw_e_var_smith_and_wilson1996_original");}
		}
		fprintf(f_timing_ptr, "\n");
	}

	FILE* f_memory_ptr = NULL;
	if(enable_output_memory){
		f_memory_ptr = fopen(output_path_memory, "w");
		if(f_memory_ptr == NULL){fprintf(stderr, "Failed to open file: %s\n", output_path_memory); return EXIT_FAILURE;}
		// fprintf(f_memory_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
		fprintf(f_memory_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn");

		if(enable_disparity_functions){
			// ----
			if(enable_functional_evenness){fprintf(f_memory_ptr, "\tm_mst_creation");}
			if(enable_distance_computation){fprintf(f_memory_ptr, "\tdist_matrix_computation");}
			if(enable_functional_evenness){fprintf(f_memory_ptr, "\tdist_heap_computation");}
			if(enable_distance_computation){fprintf(f_memory_ptr, "\tdist_stat_computation");}
			if(enable_functional_evenness){fprintf(f_memory_ptr, "\tmst_computation");}
			// ----
	
			if(enable_stirling){fprintf(f_memory_ptr, "\tstirling_alpha%.10e_beta%.10e", stirling_alpha, stirling_beta);}
			if(!enable_iterative_distance_computation && enable_ricotta_szeidl){fprintf(f_memory_ptr, "\tricotta_szeidl_alpha%.10e", ricotta_szeidl_alpha);}
			if(enable_pairwise){fprintf(f_memory_ptr, "\tpairwise");}
			if(!enable_iterative_distance_computation && enable_chao_et_al_functional_diversity){fprintf(f_memory_ptr, "\tchao_et_al_functional_diversity_alpha%.10e", chao_et_al_functional_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_scheiner_species_phylogenetic_functional_diversity){fprintf(f_memory_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e", scheiner_species_phylogenetic_functional_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_leinster_cobbold_diversity){fprintf(f_memory_ptr, "\tleinster_cobbold_diversity_alpha%.10e", leinster_cobbold_diversity_alpha);}
			if(!enable_iterative_distance_computation && enable_lexicographic){fprintf(f_memory_ptr, "\tlexicographic");}
			if(!enable_iterative_distance_computation && enable_functional_evenness){fprintf(f_memory_ptr, "\tfunctional_evenness");}
			if(!enable_iterative_distance_computation && enable_functional_dispersion){fprintf(f_memory_ptr, "\tfunctional_dispersion");}
			if(!enable_iterative_distance_computation && enable_functional_divergence_modified){fprintf(f_memory_ptr, "\tfunctional_divergence_modified");}
		}
		if(enable_non_disparity_functions){
			if(enable_shannon_weaver_entropy){fprintf(f_memory_ptr, "\tshannon_weaver_entropy");}
			if(enable_good_entropy){fprintf(f_memory_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", good_alpha, good_beta);}
			if(enable_renyi_entropy){fprintf(f_memory_ptr, "\trenyi_entropy_alpha%.4e", renyi_alpha);}
			if(enable_patil_taillie_entropy){fprintf(f_memory_ptr, "\tpatil_taillie_entropy_alpha%.4e", patil_taillie_alpha);}
			if(enable_q_logarithmic_entropy){fprintf(f_memory_ptr, "\tq_logarithmic_entropy_alpha%.4e", q_logarithmic_q);}
			if(enable_simpson_index){fprintf(f_memory_ptr, "\tsimpson_index");}
			if(enable_simpson_dominance_index){fprintf(f_memory_ptr, "\tsimpson_dominance_index");}
			if(enable_hill_number_standard){fprintf(f_memory_ptr, "\thill_number_standard_alpha%.4e", hill_number_standard_alpha);}
			if(enable_hill_evenness){fprintf(f_memory_ptr, "\thill_evenness_alpha%.4e_beta%.4e", hill_evenness_alpha, hill_evenness_beta);}
			if(enable_berger_parker_index){fprintf(f_memory_ptr, "\tberger_parker_index");}
			if(enable_junge1994_page22){fprintf(f_memory_ptr, "\tjunge1994_page22");}
			if(enable_brillouin_diversity){fprintf(f_memory_ptr, "\tbrillouin_diversity");}
			if(enable_mcintosh_index){fprintf(f_memory_ptr, "\tmcintosh_index");}
			if(enable_sw_entropy_over_log_n_species_pielou1975){fprintf(f_memory_ptr, "\tsw_entropy_over_log_n_species_pielou1975");}
			if(enable_sw_e_heip){fprintf(f_memory_ptr, "\tsw_e_heip");}
			if(enable_sw_e_one_minus_d){fprintf(f_memory_ptr, "\te_one_minus_D");}
			if(enable_sw_e_one_over_ln_d_williams1964){fprintf(f_memory_ptr, "\te_one_over_ln_D_williams1964");}
			if(enable_sw_e_minus_ln_d_pielou1977){fprintf(f_memory_ptr, "\te_minus_ln_D_pielou1977");}
			if(enable_sw_f_2_1_alatalo1981){fprintf(f_memory_ptr, "\tsw_f_2_1_alatalo1981");}
			if(enable_sw_g_2_1_molinari1989){fprintf(f_memory_ptr, "\tsw_g_2_1_molinari1989");}
			if(enable_sw_e_bulla1994){fprintf(f_memory_ptr, "\tsw_e_bulla1994");}
			if(enable_sw_o_bulla1994){fprintf(f_memory_ptr, "\tsw_o_bulla1994");}
			if(enable_sw_e_mci_pielou1969){fprintf(f_memory_ptr, "\tsw_e_mci_pielou1969");}
			if(enable_sw_e_prime_camargo1993){fprintf(f_memory_ptr, "\tsw_e_prime_camargo1993");}
			if(enable_sw_e_var_smith_and_wilson1996_original){fprintf(f_memory_ptr, "\tsw_e_var_smith_and_wilson1996_original");}
		}
		fprintf(f_memory_ptr, "\n");
	}


	int32_t num_nodes = 0;
	int64_t num_sentences = 0;
	int64_t max_num_sentences = 10000;
	int8_t use_max_num_sentences = 0;

	uint64_t num_documents = 0;
	uint64_t max_num_documents = 10000;
	uint8_t use_max_num_documents = 0;

	int64_t num_all_sentences = 0; // even ignored sentences
	int64_t max_num_all_sentences = 50000;
	int8_t use_max_num_all_sentences = 0;

	int64_t previous_g_num_nodes = -1;

	int8_t mst_initialised = 0;

	enum {
		CUPT,
		JSONL
	};

	for(int32_t i = 0 ; i < num_input_paths && ((!use_max_num_all_sentences) || num_all_sentences < max_num_all_sentences) && ((!use_max_num_sentences) || num_sentences < max_num_sentences) && ((!use_max_num_documents) || num_documents < max_num_documents); i++){
		num_documents++;
		/* // DO NOT REMOVE
		if(input_paths_true_positives == NULL){
			printf("processing %s\n", input_paths[i]);
		} else {
			printf("processing %s (true positives: %s)\n", input_paths[i], input_paths_true_positives[i]);
		}
		*/
		size_t input_path_len = strlen(input_paths[i]);
		struct cupt_sentence_iterator csi;
		struct cupt_sentence_iterator csi_tp;
		struct jsonl_document_iterator jdi;
		int32_t current_file_format = 0;
		if(strcmp(&(input_paths[i][input_path_len - 5]), ".cupt") == 0){
			current_file_format = CUPT;
			printf("CUPT: %s\n", input_paths[i]);
			int32_t err = create_cupt_sentence_iterator(&csi, input_paths[i]);
			if(err != 0){
				perror("failed to call create_cupt_sentence_iterator\n");
				return 1;
			}
			if(input_paths_true_positives != NULL){
				err = create_cupt_sentence_iterator(&csi_tp, input_paths_true_positives[i]);
				if(err != 0){
					perror("failed to call create_cupt_sentence_iterator\n");
				free_cupt_sentence_iterator(&csi);
					return 1;
				}
			}
			err = iterate_cupt_sentence_iterator(&csi);
			if(err != 0){
				perror("failed to call iterate_cupt_sentence_iterator\n");
				free_cupt_sentence_iterator(&csi);
				if(input_paths_true_positives != NULL){
					free_cupt_sentence_iterator(&csi_tp);
				}
				return 1;
			}
			if(input_paths_true_positives != NULL){
				err = iterate_cupt_sentence_iterator(&csi_tp);
				if(err != 0){
					perror("failed to call iterate_cupt_sentence_iterator\n");
					free_cupt_sentence_iterator(&csi);
					free_cupt_sentence_iterator(&csi_tp);
					return 1;
				}
			}
		} else if(strcmp(&(input_paths[i][input_path_len - 6]), ".jsonl") == 0){
			current_file_format = JSONL;
			printf("JSONL: %s\n", input_paths[i]);
			if(create_jsonl_document_iterator(&jdi, input_paths[i], jsonl_content_key) != 0){
				perror("failed to call create_jsonl_document_iterator\n");
				return 1;
			}
		} else {
			printf("unknown file type: %s\n", input_paths[i]);
			return 1;
		}

		int8_t found_at_least_one_mwe = 0;
		if(current_file_format == CUPT){
			int32_t err = 0;
			while((!csi.file_is_done) && ((!use_max_num_all_sentences) || num_all_sentences < max_num_all_sentences) && ((!use_max_num_sentences) || num_sentences < max_num_sentences)){
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
				for(int32_t j = 0 ; j < csi.current_sentence.num_tokens ; j++){
					if(target_column == UD_MWE){
						if(strcmp(csi.current_sentence.tokens[j].mwe, "") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "_") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "-") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "*") == 0){
							continue;
						}
	
						char* strtok_placeholder = NULL;
						strtok_placeholder = strtok(csi.current_sentence.tokens[j].mwe, ";");
						while(strtok_placeholder != NULL){
							int64_t mwe_num = strtol(strtok_placeholder, NULL, 10);
							if(mwe_num < max_mwe){
								size_t bytes_to_cpy = strlen(csi.current_sentence.tokens[j].lemma);
								if(bytes_to_cpy > size_token_mwe - 1){
									bytes_to_cpy = size_token_mwe - 1;
								}
								memcpy(&(mwe[mwe_num * max_tokens_per_mwe * size_token_mwe + mwe_lengths[mwe_num] * size_token_mwe]), csi.current_sentence.tokens[j].lemma, bytes_to_cpy);
								mwe_lengths[mwe_num]++;
							}
	
							if(input_paths_true_positives != NULL && (strcmp(csi_tp.current_sentence.tokens[j].mwe, "") == 0 || strcmp(csi_tp.current_sentence.tokens[j].mwe, "_") == 0 || strcmp(csi_tp.current_sentence.tokens[j].mwe, "-") == 0 || strcmp(csi_tp.current_sentence.tokens[j].mwe, "*") == 0)){
								mwe_correct_span[mwe_num] = 0;
							}
	
							strtok_placeholder = strtok(NULL, ";");
						}
	
						continue;
					}
	
	
					int32_t index;
					int32_t index_in_discarded = -1;
					const int32_t key_size = 256;
					char key[key_size];
					memset(key, '\0', key_size);
					size_t len = 0;
					if(NUM_CONLLU_COLUMNS_TO_ADD == 0){
						switch(target_column){
							case UD_FORM:
								index = word2vec_key_to_index(w2v, csi.current_sentence.tokens[j].form);
								index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, csi.current_sentence.tokens[j].form);
								len = strlen(csi.current_sentence.tokens[j].form);
								if(len > key_size - 1){
									len = key_size - 1;
								}
								memcpy(key, csi.current_sentence.tokens[j].form, len);
								break;
							case UD_LEMMA:
								index = word2vec_key_to_index(w2v, csi.current_sentence.tokens[j].lemma);
								index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, csi.current_sentence.tokens[j].lemma);
								len = strlen(csi.current_sentence.tokens[j].lemma);
								if(len > key_size - 1){
									len = key_size - 1;
								}
								memcpy(key, csi.current_sentence.tokens[j].lemma, len);
								break;
							default:
								index = -1;
								perror("target_column not properly defined\n");
								break;
						}
					} else {
						const int32_t max_size_alt_key = 256;
						char key_alt[max_size_alt_key];
						memset(key_alt, '\0', max_size_alt_key);
						int32_t key_alt_index = 0;
						for(int32_t k = 0 ; k < NUM_CONLLU_COLUMNS_TO_ADD ; k++){
							switch(CONLLU_COLUMNS_TO_ADD[k]){
								case CONLLU_COLUMN_UPOS:
									if(key_alt_index + 6 - 1 < max_size_alt_key - 1){
										memcpy(&(key_alt[key_alt_index]), "_UPOS-", 6);
										key_alt_index += 6;
										size_t size_upos = strlen(csi.current_sentence.tokens[j].upos);
										if(key_alt_index + size_upos < max_size_alt_key - 1){
											memcpy(&(key_alt[key_alt_index]), csi.current_sentence.tokens[j].upos, size_upos);
											key_alt_index += size_upos;
										}
									}
									break;
								case CONLLU_COLUMN_DEPREL:
									if(key_alt_index + 8 - 1 < max_size_alt_key - 1){
										memcpy(&(key_alt[key_alt_index]), "_DEPREL-", 8);
										key_alt_index += 8;
										size_t size_deprel = strlen(csi.current_sentence.tokens[j].deprel);
										if(key_alt_index + size_deprel < max_size_alt_key - 1){
											memcpy(&(key_alt[key_alt_index]), csi.current_sentence.tokens[j].deprel, size_deprel);
											key_alt_index += size_deprel;
										}
									}
									break;
								default:
									perror("unknown conllu column\n");
									break;
							}
						}
						if(CONLLU_ADD_FORM){
							if(key_alt_index < max_size_alt_key - 1 && NUM_CONLLU_COLUMNS_TO_ADD > 0){
								key_alt[key_alt_index] = '_';
								key_alt_index++;
							}
							size_t len_form = strlen(csi.current_sentence.tokens[j].form);
							if(key_alt_index + len_form < max_size_alt_key - 1){
								memcpy(&(key_alt[key_alt_index]), csi.current_sentence.tokens[j].form, len_form);
							}
						}
						memset(key, '\0', len);
						size_t bytes_to_cpy = strlen(key_alt);
						if(bytes_to_cpy > key_size - 1){
							bytes_to_cpy = key_size - 1;
						}
						memcpy(key, key_alt, bytes_to_cpy);
	
						index = word2vec_key_to_index(w2v, key);
						index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, key);
					}
	
					if(index != -1){
						if(w2v->keys[index].active_in_current_graph == 0){
							w2v->keys[index].active_in_current_graph = 1;
							num_nodes++;
	
							if(g.num_nodes == g.capacity){
								if(request_more_capacity_graph(&g) != 0){
									perror("failed to call request_more_capacity_graph\n");
									return 1;
								}
							}
							struct graph_node local_node;
							if(create_graph_node(&local_node, g.nodes[0].num_dimensions, FP32) != 0){
								perror("failed to call create_graph_node\n");
								return 1;
							}
							local_node.word2vec_entry_pointer = &(w2v->keys[index]);
							local_node.vector.fp32 = w2v->keys[index].vector;
	
							local_node.num_dimensions = w2v->num_dimensions;
							local_node.already_considered = 0;
							local_node.relative_proportion = 1.0;
							local_node.absolute_proportion = 1;
							g.nodes[g.num_nodes] = local_node;
							w2v->keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
							w2v->keys[index].graph_node_index = g.num_nodes;
							zipfian_n++;
							g.num_nodes++;
						} else {
							g.nodes[w2v->keys[index].graph_node_index].absolute_proportion++;
						}
						w2v->keys[index].num_occurrences++; // ? mutex ?
					} else {
						if(index_in_discarded == -1){
							struct sorted_array_str_int_element elem;
							size_t bytes_to_cpy = strlen(key);
							if(bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1){
								bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
							}
							memcpy(&elem.key, key, bytes_to_cpy);
							elem.value = 1;
							if(insert_sorted_array(&sorted_array_discarded_because_not_in_vector_database, &elem, 0) != 0){
								perror("failed to call insert_sorted_array\n");
								return 1;
							}
						} else {
							((struct sorted_array_str_int_element*) sorted_array_discarded_because_not_in_vector_database.bfr)[index_in_discarded].value++;
						}
					}
				}
	
				found_at_least_one_mwe = 0;
				if(target_column == UD_MWE){
					for(int32_t k = 0 ; k < max_mwe ; k++){
						if(mwe_lengths[k] == 0){
							continue;
						}
		
						if(input_paths_true_positives != NULL && !(mwe_correct_span[k])){
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
						for(int32_t m = 0 ; m < mwe_lengths[k] ; m++){
							bytes_to_cpy = strlen(&(mwe[k * max_tokens_per_mwe * size_token_mwe + m * size_token_mwe]));
							if(m < mwe_lengths[k] - 1){
								bytes_to_cpy++;
							}
							if(bytes_to_cpy > (size_t) (bfr_size - 1 - index_bfr)){
								bytes_to_cpy = bfr_size - 1 - index_bfr;
							}
							if(m < mwe_lengths[k] - 1){
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
						if(index != -1){
							found_at_least_one_mwe = 1;
							if(w2v->keys[index].active_in_current_graph == 0){
								w2v->keys[index].active_in_current_graph = 1;
								if(g.num_nodes == g.capacity){
									err = request_more_capacity_graph(&g);
									if(err != 0){
										perror("failed to call request_more_capacity_graph\n");
										return 1;
									}
								}
								struct graph_node local_graph_node;
								local_graph_node.num_dimensions = (int16_t) w2v->num_dimensions;
								local_graph_node.already_considered = 0;
								local_graph_node.relative_proportion = 1.0;
								local_graph_node.absolute_proportion = 1;
								local_graph_node.vector.fp32 = w2v->keys[index].vector;
								local_graph_node.word2vec_entry_pointer = &(w2v->keys[index]);
								g.nodes[g.num_nodes] = local_graph_node;
								w2v->keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
								w2v->keys[index].graph_node_index = g.num_nodes;
								zipfian_n++;
								g.num_nodes++;
								num_nodes++;
							} else {
								g.nodes[w2v->keys[index].graph_node_index].absolute_proportion++;
							}
							w2v->keys[index].num_occurrences++;
						} else {
							int32_t index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, bfr);
							
							if(index_in_discarded == -1){
								struct sorted_array_str_int_element elem;
								size_t bytes_to_cpy = strlen(bfr);
								if(bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1){
									bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
								}
								memcpy(&elem.key, bfr, bytes_to_cpy);
								elem.value = 1;
								if(insert_sorted_array(&sorted_array_discarded_because_not_in_vector_database, &elem, 0) != 0){
									perror("failed to call insert_sorted_array\n");
									return 1;
								}
							} else {
								((struct sorted_array_str_int_element*) sorted_array_discarded_because_not_in_vector_database.bfr)[index_in_discarded].value++;
							}
						} 
					}
		
					if(found_at_least_one_mwe){
						num_sentences++;
					}
				}

				num_all_sentences++;
				if(target_column != UD_MWE){
					num_sentences = num_all_sentences;
				}	

				// sentence level recomputation
				if((target_column != UD_MWE || found_at_least_one_mwe) && (enable_sentence_count_recompute_step && (((!sentence_recompute_step_use_log10) && num_sentences % sentence_count_recompute_step == 0) || (sentence_recompute_step_use_log10 && ((uint64_t) num_sentences) >= stacked_sentence_count_target))) && g.num_nodes > 1){
					printf("found_at_least_one_mwe: %i; g.num_nodes: %lu\n", found_at_least_one_mwe, g.num_nodes);
					double absolute_proportion_sum = 0.0;
					for(uint64_t p = 0 ; p < g.num_nodes ; p++){
						absolute_proportion_sum += (double) g.nodes[p].absolute_proportion;
					}
					for(uint64_t p = 0 ; p < g.num_nodes ; p++){
						g.nodes[p].relative_proportion = ((double) g.nodes[p].absolute_proportion) / absolute_proportion_sum;
					}
				
					double best_s = -1.0;
					int32_t err = zipfian_fit_from_graph(&g, &best_s);
					if(err != 0){
						perror("failed to call zipfian_fit_from_graph\n");
						goto free_bfr_exit_failure;
					}
				
					if(best_s != previous_best_s || g.num_nodes != ((uint64_t) previous_g_num_nodes)){
						printf("best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li\n", best_s, g.num_nodes, num_sentences, num_documents);
						err = apply_diversity_functions_to_graph(&g, &mst, &heap, f_ptr, f_timing_ptr, f_memory_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database,
	// target_column,
	num_row_threads,
	num_matrix_threads,
	enable_stirling,
	enable_ricotta_szeidl,
	enable_pairwise,
	enable_lexicographic,
	enable_chao_et_al_functional_diversity,
	enable_scheiner_species_phylogenetic_functional_diversity,
	enable_leinster_cobbold_diversity,
	enable_multithreaded_matrix_generation,
	enable_timings,
	enable_iterative_distance_computation,
	enable_multithreaded_row_generation,
	row_generation_batch_size,
// 	enable_sentence_count_recompute_step,
// 	sentence_recompute_step_use_log10,
// 	enable_document_count_recompute_step,
// 	document_recompute_step_use_log10,
 	enable_output_timing,
 	enable_output_memory,
// 	sentence_count_recompute_step,
// 	sentence_count_recompute_step_log10,
// 	document_count_recompute_step,
// 	document_count_recompute_step_log10,
	enable_functional_evenness,
	enable_functional_dispersion,
	enable_functional_divergence_modified,
	enable_non_disparity_functions,
	enable_disparity_functions,
	enable_shannon_weaver_entropy,
	enable_good_entropy,
	enable_renyi_entropy,
	enable_patil_taillie_entropy,
	enable_q_logarithmic_entropy,
	enable_simpson_index,
	enable_simpson_dominance_index,
	enable_hill_number_standard,
	enable_hill_evenness,
	enable_berger_parker_index,
	enable_junge1994_page22,
	enable_brillouin_diversity,
	enable_mcintosh_index,
	enable_sw_entropy_over_log_n_species_pielou1975,
	enable_sw_e_heip,
	enable_sw_e_one_minus_d,
	enable_sw_e_one_over_ln_d_williams1964,
	enable_sw_e_minus_ln_d_pielou1977,
	enable_sw_f_2_1_alatalo1981,
	enable_sw_g_2_1_molinari1989,
	enable_sw_e_bulla1994,
	enable_sw_o_bulla1994,
	enable_sw_e_mci_pielou1969,
	enable_sw_e_prime_camargo1993,
	enable_sw_e_var_smith_and_wilson1996_original,
	stirling_alpha,
	stirling_beta,
	ricotta_szeidl_alpha,
	chao_et_al_functional_diversity_alpha,
	scheiner_species_phylogenetic_functional_diversity_alpha,
	leinster_cobbold_diversity_alpha,
	good_alpha,
	good_beta,
	renyi_alpha,
	patil_taillie_alpha,
	q_logarithmic_q,
	hill_number_standard_alpha,
	hill_evenness_alpha,
	hill_evenness_beta
								);
						if(err != 0){
							perror("failed to call apply_diversity_functions_to_graph\n");
							return 1;
						}
						previous_best_s = best_s;
					} else { // end of comparison best_s?
						printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
					}

					if(sentence_recompute_step_use_log10){
						stacked_sentence_count_log10 += sentence_count_recompute_step_log10;
						stacked_sentence_count_target = (uint64_t) floor(pow(10.0, stacked_sentence_count_log10));
						printf("New sentence count target: %lu (10.0^%.3f)\n", stacked_sentence_count_target, stacked_sentence_count_log10);
					}
				}



				err = iterate_cupt_sentence_iterator(&csi);
				if(err != 0){
					perror("failed to call iterate_cupt_sentence_iterator\n");
					free_cupt_sentence_iterator(&csi);
					if(input_paths_true_positives != NULL){
						free_cupt_sentence_iterator(&csi_tp);
					}
					return 1;
				}
				if(input_paths_true_positives != NULL){
					err = iterate_cupt_sentence_iterator(&csi_tp);
					if(err != 0){
						perror("failed to call iterate_cupt_sentence_iterator\n");
						free_cupt_sentence_iterator(&csi);
						free_cupt_sentence_iterator(&csi_tp);
						return 1;
					}
				}
			}

			// document level recomputation

			if((target_column != UD_MWE || found_at_least_one_mwe) && (enable_document_count_recompute_step && (((!document_recompute_step_use_log10) && num_documents % document_count_recompute_step == 0) || (document_recompute_step_use_log10 && num_documents >= stacked_document_count_target)) && g.num_nodes > 1)){
				printf("found_at_least_one_mwe: %i; g.num_nodes: %lu\n", found_at_least_one_mwe, g.num_nodes);
				double absolute_proportion_sum = 0.0;
				for(uint64_t p = 0 ; p < g.num_nodes ; p++){
					absolute_proportion_sum += (double) g.nodes[p].absolute_proportion;
				}
				for(uint64_t p = 0 ; p < g.num_nodes ; p++){
					g.nodes[p].relative_proportion = ((double) g.nodes[p].absolute_proportion) / absolute_proportion_sum;
				}
			
				double best_s = -1.0;
				int32_t err = zipfian_fit_from_graph(&g, &best_s);
				if(err != 0){
					perror("failed to call zipfian_fit_from_graph\n");
					goto free_bfr_exit_failure;
				}
			
				if(best_s != previous_best_s || g.num_nodes != ((uint64_t) previous_g_num_nodes)){
					printf("best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li\n", best_s, g.num_nodes, num_sentences, num_documents);
					err = apply_diversity_functions_to_graph(&g, &mst, &heap, f_ptr, f_timing_ptr, f_memory_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database,
	// target_column,
	num_row_threads,
	num_matrix_threads,
	enable_stirling,
	enable_ricotta_szeidl,
	enable_pairwise,
	enable_lexicographic,
	enable_chao_et_al_functional_diversity,
	enable_scheiner_species_phylogenetic_functional_diversity,
	enable_leinster_cobbold_diversity,
	enable_multithreaded_matrix_generation,
	enable_timings,
	enable_iterative_distance_computation,
	enable_multithreaded_row_generation,
	row_generation_batch_size,
// 	enable_sentence_count_recompute_step,
// 	sentence_recompute_step_use_log10,
// 	enable_document_count_recompute_step,
// 	document_recompute_step_use_log10,
	enable_output_timing,
	enable_output_memory,
// 	sentence_count_recompute_step,
// 	sentence_count_recompute_step_log10,
// 	document_count_recompute_step,
// 	document_count_recompute_step_log10,
	enable_functional_evenness,
	enable_functional_dispersion,
	enable_functional_divergence_modified,
	enable_non_disparity_functions,
	enable_disparity_functions,
	enable_shannon_weaver_entropy,
	enable_good_entropy,
	enable_renyi_entropy,
	enable_patil_taillie_entropy,
	enable_q_logarithmic_entropy,
	enable_simpson_index,
	enable_simpson_dominance_index,
	enable_hill_number_standard,
	enable_hill_evenness,
	enable_berger_parker_index,
	enable_junge1994_page22,
	enable_brillouin_diversity,
	enable_mcintosh_index,
	enable_sw_entropy_over_log_n_species_pielou1975,
	enable_sw_e_heip,
	enable_sw_e_one_minus_d,
	enable_sw_e_one_over_ln_d_williams1964,
	enable_sw_e_minus_ln_d_pielou1977,
	enable_sw_f_2_1_alatalo1981,
	enable_sw_g_2_1_molinari1989,
	enable_sw_e_bulla1994,
	enable_sw_o_bulla1994,
	enable_sw_e_mci_pielou1969,
	enable_sw_e_prime_camargo1993,
	enable_sw_e_var_smith_and_wilson1996_original,
	stirling_alpha,
	stirling_beta,
	ricotta_szeidl_alpha,
	chao_et_al_functional_diversity_alpha,
	scheiner_species_phylogenetic_functional_diversity_alpha,
	leinster_cobbold_diversity_alpha,
	good_alpha,
	good_beta,
	renyi_alpha,
	patil_taillie_alpha,
	q_logarithmic_q,
	hill_number_standard_alpha,
	hill_evenness_alpha,
	hill_evenness_beta
							);
					if(err != 0){
						perror("failed to call apply_diversity_functions_to_graph\n");
						return 1;
					}
					previous_best_s = best_s;
				} else { // end of comparison best_s?
					printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
				}

				if(document_recompute_step_use_log10){
					stacked_document_count_log10 += document_count_recompute_step_log10;
					stacked_document_count_target = (uint64_t) floor(pow(10.0, stacked_document_count_log10));
					printf("New document count target: %lu (10.0^%.3f)\n", stacked_document_count_target, stacked_document_count_log10);
				}
			}
		} else if(current_file_format == JSONL){
			while(!(jdi.file_is_done)){
				memset(jdi.current_document.identifier, '\0', jdi.current_document.identifier_size); // ?
				jdi.current_document.identifier_size = 0;
				memset(jdi.current_document.text, '\0', jdi.current_document.text_size); // ?
				jdi.current_document.text_size = 0;
				if(iterate_jsonl_document_iterator(&jdi) != 0){
					perror("failed to call iterate_jsonl_document_iterator\n");
					return 1;
				}
				if(jdi.current_document.text_size == 0 || jdi.current_document.identifier_size == 0){
					continue;
				}
				while(!(jdi.current_document.reached_last_token)){
					if(iterate_document_current_token(&(jdi.current_document)) != 0){
						perror("failed to call iterate_document_current_token\n");
						return 1;
					}

					// add to graph
					int32_t index = word2vec_key_to_index(w2v, jdi.current_document.current_token);
					if(index != -1){
						if(w2v->keys[index].active_in_current_graph == 0){
							w2v->keys[index].active_in_current_graph = 1;
							num_nodes++;

							if(g.num_nodes == g.capacity){
								if(request_more_capacity_graph(&g) != 0){
									perror("failed to call request_more_capacity_graph\n");
									return 1;
								}
							}
							struct graph_node local_node;
							if(create_graph_node(&local_node, g.nodes[0].num_dimensions, FP32) != 0){
								perror("failed to call create_graph_node\n");
								return 1;
							}
							local_node.word2vec_entry_pointer = &(w2v->keys[index]);
							local_node.vector.fp32 = w2v->keys[index].vector;
							local_node.num_dimensions = w2v->num_dimensions;
							local_node.already_considered = 0;
							local_node.relative_proportion = 1.0;
							local_node.absolute_proportion = 1;
							g.nodes[g.num_nodes] = local_node;
							w2v->keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
							w2v->keys[index].graph_node_index = g.num_nodes;
							zipfian_n++;
							g.num_nodes++;
						} else {	
							g.nodes[w2v->keys[index].graph_node_index].absolute_proportion++;
						}
					}
				}
				if((target_column != UD_MWE || found_at_least_one_mwe) && (enable_document_count_recompute_step && (((!document_recompute_step_use_log10) && num_documents % document_count_recompute_step == 0) || (document_recompute_step_use_log10 && num_documents >= stacked_document_count_target)) && g.num_nodes > 1)){
					double absolute_proportion_sum = 0.0;
					for(uint64_t p = 0 ; p < g.num_nodes ; p++){
						absolute_proportion_sum += (double) g.nodes[p].absolute_proportion;
					}
					for(uint64_t p = 0 ; p < g.num_nodes ; p++){
						g.nodes[p].relative_proportion = ((double) g.nodes[p].absolute_proportion) / absolute_proportion_sum;
					}
				
					double best_s = -1.0;
					int32_t err = zipfian_fit_from_graph(&g, &best_s);
					if(err != 0){
						perror("failed to call zipfian_fit_from_graph\n");
						goto free_bfr_exit_failure;
					}
				
					if(best_s != previous_best_s || g.num_nodes != ((uint64_t) previous_g_num_nodes)){
						printf("best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li\n", best_s, g.num_nodes, num_sentences, num_documents);
						err = apply_diversity_functions_to_graph(&g, &mst, &heap, f_ptr, f_timing_ptr, f_memory_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database,
	// target_column,
	num_row_threads,
	num_matrix_threads,
	enable_stirling,
	enable_ricotta_szeidl,
	enable_pairwise,
	enable_lexicographic,
	enable_chao_et_al_functional_diversity,
	enable_scheiner_species_phylogenetic_functional_diversity,
	enable_leinster_cobbold_diversity,
	enable_multithreaded_matrix_generation,
	enable_timings,
	enable_iterative_distance_computation,
	enable_multithreaded_row_generation,
	row_generation_batch_size,
// 	enable_sentence_count_recompute_step,
// 	sentence_recompute_step_use_log10,
// 	enable_document_count_recompute_step,
// 	document_recompute_step_use_log10,
	enable_output_timing,
	enable_output_memory,
// 	sentence_count_recompute_step,
// 	sentence_count_recompute_step_log10,
// 	document_count_recompute_step,
// 	document_count_recompute_step_log10,
	enable_functional_evenness,
	enable_functional_dispersion,
	enable_functional_divergence_modified,
	enable_non_disparity_functions,
	enable_disparity_functions,
	enable_shannon_weaver_entropy,
	enable_good_entropy,
	enable_renyi_entropy,
	enable_patil_taillie_entropy,
	enable_q_logarithmic_entropy,
	enable_simpson_index,
	enable_simpson_dominance_index,
	enable_hill_number_standard,
	enable_hill_evenness,
	enable_berger_parker_index,
	enable_junge1994_page22,
	enable_brillouin_diversity,
	enable_mcintosh_index,
	enable_sw_entropy_over_log_n_species_pielou1975,
	enable_sw_e_heip,
	enable_sw_e_one_minus_d,
	enable_sw_e_one_over_ln_d_williams1964,
	enable_sw_e_minus_ln_d_pielou1977,
	enable_sw_f_2_1_alatalo1981,
	enable_sw_g_2_1_molinari1989,
	enable_sw_e_bulla1994,
	enable_sw_o_bulla1994,
	enable_sw_e_mci_pielou1969,
	enable_sw_e_prime_camargo1993,
	enable_sw_e_var_smith_and_wilson1996_original,
	stirling_alpha,
	stirling_beta,
	ricotta_szeidl_alpha,
	chao_et_al_functional_diversity_alpha,
	scheiner_species_phylogenetic_functional_diversity_alpha,
	leinster_cobbold_diversity_alpha,
	good_alpha,
	good_beta,
	renyi_alpha,
	patil_taillie_alpha,
	q_logarithmic_q,
	hill_number_standard_alpha,
	hill_evenness_alpha,
	hill_evenness_beta
								);
						if(err != 0){
							perror("failed to call apply_diversity_functions_to_graph\n");
							return 1;
						}
						previous_best_s = best_s;
					} else { // end of comparison best_s?
						printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
					}

					if(document_recompute_step_use_log10){
						stacked_document_count_log10 += document_count_recompute_step_log10;
						stacked_document_count_target = (uint64_t) floor(pow(10.0, stacked_document_count_log10));
						printf("New document count target: %lu (10.0^%.3f)\n", stacked_document_count_target, stacked_document_count_log10);
					}
				}

				num_documents++;
			}
		}

		if(current_file_format == CUPT){
			free_cupt_sentence_iterator(&csi);
			if(input_paths_true_positives != NULL){
				free_cupt_sentence_iterator(&csi_tp);
			}
		} else if(current_file_format == JSONL){
			free_jsonl_document_iterator(&jdi);
		}
	}

	free_graph(&g);
	free_graph_distance_heap(&heap);
	free_minimum_spanning_tree(&mst);

	free_sorted_array(&sorted_array_discarded_because_not_in_vector_database);

	for(int32_t i = 0 ; i < num_files ; i++){
		free(bfr[i]);
	}

	fclose(f_ptr);
	if(f_timing_ptr != NULL){fclose(f_timing_ptr);}
	if(f_memory_ptr != NULL){fclose(f_memory_ptr);}

	return 0;

	free_bfr_exit_failure:
	for(int32_t i = 0 ; i < num_files ; i++){
		free(bfr[i]);
	}
	return 1;
}

int32_t main(int32_t argc, char** argv){
	char* argv_w2v_path = NULL;
	char* argv_jsonl_content_key = NULL;
	char* argv_input_path = NULL;
	char* argv_output_path = NULL;
	char* argv_output_path_timing = NULL;
	char* argv_output_path_memory = NULL;
	uint32_t argv_target_column = TARGET_COLUMN;
	uint32_t argv_num_row_threads = NUM_ROW_THREADS;
	uint32_t argv_num_matrix_threads = NUM_MATRIX_THREADS;
	uint8_t argv_enable_stirling = ENABLE_STIRLING;
	uint8_t argv_enable_ricotta_szeidl = ENABLE_RICOTTA_SZEIDL;
	uint8_t argv_enable_pairwise = ENABLE_PAIRWISE;
	uint8_t argv_enable_lexicographic = ENABLE_LEXICOGRAPHIC;
	uint8_t argv_enable_chao_et_al_functional_diversity = ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY;
	uint8_t argv_enable_scheiner_species_phylogenetic_functional_diversity = ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY;
	uint8_t argv_enable_leinster_cobbold_diversity = ENABLE_LEINSTER_COBBOLD_DIVERSITY;
	uint8_t argv_enable_multithreaded_matrix_generation = ENABLE_MULTITHREADED_MATRIX_GENERATION;
	uint8_t argv_enable_timings = ENABLE_TIMINGS;
	uint8_t argv_enable_iterative_distance_computation = ENABLE_ITERATIVE_DISTANCE_COMPUTATION;
	uint8_t argv_enable_multithreaded_row_generation = ENABLE_MULTITHREADED_ROW_GENERATION;
	uint8_t argv_row_generation_batch_size = ROW_GENERATION_BATCH_SIZE;
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
	uint8_t argv_enable_sw_e_var_smith_and_wilson1996_original = ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL;
	double argv_stirling_alpha = STIRLING_ALPHA;
	double argv_stirling_beta = STIRLING_BETA;
	double argv_ricotta_szeidl_alpha = RICOTTA_SZEIDL_ALPHA;
	double argv_chao_et_al_functional_diversity_alpha = CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA;
	double argv_scheiner_species_phylogenetic_functional_diversity_alpha = SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA;
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

	for(int32_t i = 1 ; i < argc ; i++){
		if(strncmp(argv[i], "--w2v_path=", 11) == 0){argv_w2v_path = argv[i] + 11;}
		else if(strncmp(argv[i], "--target_column=", 16) == 0){
			// argv_target_column = argv[i] + 16;
			if(strcmp(argv[i] + 16, "UD_MWE") == 0){argv_target_column = UD_MWE;}
			else if(strcmp(argv[i] + 16, "UD_FORM") == 0){argv_target_column = UD_FORM;}
			else if(strcmp(argv[i] + 16, "UD_LEMMA") == 0){argv_target_column = UD_LEMMA;}
			else {fprintf(stderr, "Unknown target_column: %s\n", argv[i] + 16); return 1;}
		}
		else if(strncmp(argv[i], "--num_row_threads=", 18) == 0){argv_num_row_threads = (uint32_t) strtol(argv[i] + 18, NULL, 10);}
		else if(strncmp(argv[i], "--num_matrix_threads=", 21) == 0){argv_num_matrix_threads = (uint32_t) strtol(argv[i] + 21, NULL, 10);}
		else if(strncmp(argv[i], "--jsonl_content_key=", 20) == 0){argv_jsonl_content_key = argv[i] + 20;}
		else if(strncmp(argv[i], "--input_path=", 13) == 0){argv_input_path = argv[i] + 13;}
		else if(strncmp(argv[i], "--output_path=", 14) == 0){argv_output_path = argv[i] + 14;}
		else if(strncmp(argv[i], "--output_path_timing=", 21) == 0){argv_output_path_timing = argv[i] + 21;}
		else if(strncmp(argv[i], "--output_path_memory=", 21) == 0){argv_output_path_memory = argv[i] + 21;}
		else if(strncmp(argv[i], "--enable_multithreaded_matrix_generation=", 41) == 0){argv_enable_multithreaded_matrix_generation = (argv[i][41] == '1');}
		else if(strncmp(argv[i], "--enable_timings=", 17) == 0){argv_enable_timings = (argv[i][18] == '1');}
		else if(strncmp(argv[i], "--enable_iterative_distance_computation=", 40) == 0){argv_enable_iterative_distance_computation = (argv[i][40] == '1');}
		else if(strncmp(argv[i], "--enable_sentence_count_recompute_step=", 39) == 0){argv_enable_sentence_count_recompute_step = (argv[i][39] == '1');}
		else if(strncmp(argv[i], "--enable_document_count_recompute_step=", 39) == 0){argv_enable_document_count_recompute_step = (argv[i][39] == '1');}
		else if(strncmp(argv[i], "--sentence_recompute_step_use_log10=", 36) == 0){argv_sentence_recompute_step_use_log10 = (argv[i][36] == '1');}
		else if(strncmp(argv[i], "--document_recompute_step_use_log10=", 36) == 0){argv_document_recompute_step_use_log10 = (argv[i][36] == '1');}
		else if(strncmp(argv[i], "--enable_output_timing=", 23) == 0){argv_enable_output_timing = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_output_memory=", 23) == 0){argv_enable_output_memory = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_stirling=", 18) == 0){argv_enable_stirling = (argv[i][18] == '1');}
		else if(strncmp(argv[i], "--enable_ricotta_szeidl=", 24) == 0){argv_enable_ricotta_szeidl = (argv[i][24] == '1');}
		else if(strncmp(argv[i], "--enable_pairwise=", 18) == 0){argv_enable_pairwise = (argv[i][18] == '1');}
		else if(strncmp(argv[i], "--enable_lexicographic=", 23) == 0){argv_enable_lexicographic = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_chao_et_al_functional_diversity=", 41) == 0){argv_enable_chao_et_al_functional_diversity = (argv[i][41] == '1');}
		else if(strncmp(argv[i], "--enable_scheiner_species_phylogenetic_functional_diversity=", 60) == 0){argv_enable_scheiner_species_phylogenetic_functional_diversity = (argv[i][60] == '1');}
		else if(strncmp(argv[i], "--enable_leinster_cobbold_diversity=", 36) == 0){argv_enable_leinster_cobbold_diversity = (argv[i][36] == '1');}
		else if(strncmp(argv[i], "--enable_functional_evenness=", 29) == 0){argv_enable_functional_evenness = (argv[i][29] == '1');}
		else if(strncmp(argv[i], "--enable_functional_dispersion=", 31) == 0){argv_enable_functional_dispersion = (argv[i][31] == '1');}
		else if(strncmp(argv[i], "--enable_functional_divergence_modified=", 40) == 0){argv_enable_functional_divergence_modified = (argv[i][40] == '1');}
		else if(strncmp(argv[i], "--enable_non_disparity_functions=", 33) == 0){argv_enable_non_disparity_functions = (argv[i][33] == '1');}
		else if(strncmp(argv[i], "--enable_disparity_functions=", 29) == 0){argv_enable_disparity_functions = (argv[i][29] == '1');}
		else if(strncmp(argv[i], "--enable_shannon_weaver_entropy=", 32) == 0){argv_enable_shannon_weaver_entropy = (argv[i][32] == '1');}
		else if(strncmp(argv[i], "--enable_good_entropy=", 22) == 0){argv_enable_good_entropy = (argv[i][22] == '1');}
		else if(strncmp(argv[i], "--enable_renyi_entropy=", 23) == 0){argv_enable_renyi_entropy = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_patil_taillie_entropy=", 31) == 0){argv_enable_patil_taillie_entropy = (argv[i][31] == '1');}
		else if(strncmp(argv[i], "--enable_q_logarithmic_entropy=", 31) == 0){argv_enable_q_logarithmic_entropy = (argv[i][31] == '1');}
		else if(strncmp(argv[i], "--enable_simpson_index=", 23) == 0){argv_enable_simpson_index = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_simpson_dominance_index=", 33) == 0){argv_enable_simpson_dominance_index = (argv[i][33] == '1');}
		else if(strncmp(argv[i], "--enable_hill_number_standard=", 30) == 0){argv_enable_hill_number_standard = (argv[i][30] == '1');}
		else if(strncmp(argv[i], "--enable_hill_evenness=", 23) == 0){argv_enable_hill_evenness = (argv[i][23] == '1');}
		else if(strncmp(argv[i], "--enable_berger_parker_index=", 29) == 0){argv_enable_berger_parker_index = (argv[i][29] == '1');}
		else if(strncmp(argv[i], "--enable_junge1994_page22=", 26) == 0){argv_enable_junge1994_page22 = (argv[i][26] == '1');}
		else if(strncmp(argv[i], "--enable_brillouin_diversity=", 29) == 0){argv_enable_brillouin_diversity = (argv[i][29] == '1');}
		else if(strncmp(argv[i], "--enable_mcintosh_index=", 24) == 0){argv_enable_mcintosh_index = (argv[i][24] == '1');}
		else if(strncmp(argv[i], "--enable_sw_entropy_over_log_n_species_pielou1975=", 50) == 0){argv_enable_sw_entropy_over_log_n_species_pielou1975 = (argv[i][50] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_heip=", 19) == 0){argv_enable_sw_e_heip = (argv[i][19] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_one_over_ln_d_williams1964=", 41) == 0){argv_enable_sw_e_one_over_ln_d_williams1964 = (argv[i][41] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_minus_ln_d_pielou1977=", 36) == 0){argv_enable_sw_e_minus_ln_d_pielou1977 = (argv[i][36] == '1');}
		else if(strncmp(argv[i], "--enable_sw_f_2_1_alatalo1981=", 30) == 0){argv_enable_sw_f_2_1_alatalo1981 = (argv[i][30] == '1');}
		else if(strncmp(argv[i], "--enable_sw_g_2_1_molinari1989=", 31) == 0){argv_enable_sw_g_2_1_molinari1989 = (argv[i][31] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_bulla1994=", 24) == 0){argv_enable_sw_e_bulla1994 = (argv[i][24] == '1');}
		else if(strncmp(argv[i], "--enable_sw_o_bulla1994=", 24) == 0){argv_enable_sw_o_bulla1994 = (argv[i][24] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_mci_pielou1969=", 29) == 0){argv_enable_sw_e_mci_pielou1969 = (argv[i][29] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_prime_camargo1993=", 32) == 0){argv_enable_sw_e_prime_camargo1993 = (argv[i][32] == '1');}
		else if(strncmp(argv[i], "--enable_sw_e_var_smith_and_wilson1996_original=", 48) == 0){argv_enable_sw_e_var_smith_and_wilson1996_original = (argv[i][48] == '1');}
		else if(strncmp(argv[i], "--stirling_alpha=", 17) == 0){argv_stirling_alpha = strtod(argv[i] + 17, NULL);}
		else if(strncmp(argv[i], "--stirling_beta=", 16) == 0){argv_stirling_beta = strtod(argv[i] + 16, NULL);}
		else if(strncmp(argv[i], "--ricotta_szeidl_alpha=", 23) == 0){argv_ricotta_szeidl_alpha = strtod(argv[i] + 23, NULL);}
		else if(strncmp(argv[i], "--chao_et_al_functional_diversity_alpha=", 40) == 0){argv_chao_et_al_functional_diversity_alpha = strtod(argv[i] + 40, NULL);}
		else if(strncmp(argv[i], "--scheiner_species_phylogenetic_functional_diversity_alpha=", 59) == 0){argv_scheiner_species_phylogenetic_functional_diversity_alpha = strtod(argv[i] + 59, NULL);}
		else if(strncmp(argv[i], "--leinster_cobbold_diversity_alpha=", 35) == 0){argv_leinster_cobbold_diversity_alpha = strtod(argv[i] + 35, NULL);}
		else if(strncmp(argv[i], "--renyi_alpha=", 14) == 0){argv_renyi_alpha = strtod(argv[i] + 14, NULL);}
		else if(strncmp(argv[i], "--patil_taillie_alpha=", 22) == 0){argv_patil_taillie_alpha = strtod(argv[i] + 22, NULL);}
		else if(strncmp(argv[i], "--hill_number_standard_alpha=", 29) == 0){argv_hill_number_standard_alpha = strtod(argv[i] + 29, NULL);}
		else if(strncmp(argv[i], "--hill_evenness_alpha=", 22) == 0){argv_hill_evenness_alpha = strtod(argv[i] + 22, NULL);}
		else if(strncmp(argv[i], "--hill_evenness_beta=", 21) == 0){argv_hill_evenness_beta = strtod(argv[i] + 21, NULL);}
		else if(strncmp(argv[i], "--good_alpha=", 13) == 0){argv_good_alpha = strtod(argv[i] + 13, NULL);}
		else if(strncmp(argv[i], "--good_beta=", 12) == 0){argv_good_beta = strtod(argv[i] + 12, NULL);}
		else if(strncmp(argv[i], "--row_generation_batch_size=", 28) == 0){argv_row_generation_batch_size = strtol(argv[i] + 28, NULL, 10);}
		else if(strncmp(argv[i], "--sentence_count_recompute_step=", 32) == 0){argv_sentence_count_recompute_step = strtol(argv[i] + 32, NULL, 10);}
		else if(strncmp(argv[i], "--sentence_count_recompute_step_log10=", 38) == 0){argv_sentence_count_recompute_step_log10 = strtod(argv[i] + 38, NULL);}
		else if(strncmp(argv[i], "--document_count_recompute_step=", 32) == 0){argv_document_count_recompute_step = strtol(argv[i] + 32, NULL, 10);}
		else if(strncmp(argv[i], "--document_count_recompute_step_log10=", 38) == 0){argv_document_count_recompute_step_log10 = strtod(argv[i] + 38, NULL);}
		else if(strncmp(argv[i], "--force_timing_and_memory_to_output_path=", 41) == 0){argv_force_timing_and_memory_to_output_path = (argv[i][41] == '1');}
		else {fprintf(stderr, "Unknown argument: %s\n", argv[i]); return 1;}
	}

	if(argv_w2v_path == NULL){argv_w2v_path = W2V_PATH;}
	if(argv_jsonl_content_key == NULL){argv_jsonl_content_key = JSONL_CONTENT_KEY;}
	if(argv_input_path == NULL){argv_input_path = INPUT_PATH;}
	if(argv_output_path == NULL){argv_output_path = OUTPUT_PATH;}
	if(argv_output_path_timing == NULL){argv_output_path_timing = OUTPUT_PATH_TIMING;}
	if(argv_output_path_memory == NULL){argv_output_path_memory = OUTPUT_PATH_MEMORY;}

	if(argv_force_timing_and_memory_to_output_path){
		size_t path_len, delta, alloc_size, len_suffix;
		char* ptr_last_slash;

		ptr_last_slash = strrchr(argv_output_path, '/');
		if(ptr_last_slash == NULL){ptr_last_slash = argv_output_path - 1;}
		// ptr_last_slash++;

		delta = ((size_t) (ptr_last_slash - argv_output_path)) + 1;
		len_suffix = strlen("measurement_output_timing.tsv");
		path_len = delta + len_suffix;
		alloc_size = path_len + 1;
		argv_output_path_timing = malloc(alloc_size);
		if(argv_output_path_timing == NULL){goto malloc_fail;}
		memset(argv_output_path_timing, '\0', alloc_size);
		memcpy(argv_output_path_timing, argv_output_path, delta);
		memcpy(argv_output_path_timing + delta, "measurement_output_timing.tsv", len_suffix);
		
		len_suffix = strlen("measurement_output_memory.tsv");
		path_len = delta + len_suffix;
		alloc_size = path_len + 1;
		argv_output_path_memory = malloc(alloc_size);
		if(argv_output_path_memory == NULL){free(argv_output_path_timing); goto malloc_fail;}
		memset(argv_output_path_memory, '\0', alloc_size);
		memcpy(argv_output_path_memory, argv_output_path, delta);
		memcpy(argv_output_path_memory + delta, "measurement_output_memory.tsv", len_suffix);
	}

	printf("w2v_path: %s\n", argv_w2v_path);
	printf("jsonl_content_key: %s\n", argv_jsonl_content_key);
	printf("input_path: %s\n", argv_input_path);
	printf("output_path: %s\n", argv_output_path);
	printf("force_timing_and_memory_to_output_path: %u\n", argv_force_timing_and_memory_to_output_path);
	printf("output_path_timing: %s\n", argv_output_path_timing);
	printf("output_path_memory: %s\n", argv_output_path_memory);

	printf("target_column: %u\n", argv_target_column);
	printf("num_row_threads: %u\n", argv_num_row_threads);
	printf("num_matrix_threads: %u\n", argv_num_matrix_threads);

	printf("enable_multithreaded_matrix_generation: %u\n", argv_enable_multithreaded_matrix_generation);
	printf("enable_timings: %u\n", argv_enable_timings);
	printf("enable_iterative_distance_computation: %u\n", argv_enable_iterative_distance_computation);
	printf("enable_multithreaded_row_generation: %u\n", argv_enable_multithreaded_row_generation);
	printf("row_generation_batch_size: %u\n", argv_row_generation_batch_size);
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
	printf("enable_scheiner_species_phylogenetic_functional_diversity: %u\n", argv_enable_scheiner_species_phylogenetic_functional_diversity);
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
	printf("scheiner_species_phylogenetic_functional_diversity_alpha: %f\n", argv_scheiner_species_phylogenetic_functional_diversity_alpha);
	printf("leinster_cobbold_diversity_alpha: %f\n", argv_leinster_cobbold_diversity_alpha);
	printf("good_alpha: %f\n", argv_good_alpha);
	printf("good_beta: %f\n", argv_good_beta);
	printf("renyi_alpha: %f\n", argv_renyi_alpha);
	printf("patil_taillie_alpha: %f\n", argv_patil_taillie_alpha);
	printf("q_logarithmic_q: %f\n", argv_q_logarithmic_q);
	printf("hill_number_standard_alpha: %f\n", argv_hill_number_standard_alpha);
	printf("hill_evenness_alpha: %f\n", argv_hill_evenness_alpha);
	printf("hill_evenness_beta: %f\n", argv_hill_evenness_beta);

	jsonl_init_tokenization();

	int32_t err;

	struct word2vec w2v;
	err = load_word2vec_binary(&w2v, argv_w2v_path);
	if(err != 0){
		fprintf(stderr, "failed to call load_word2vec binary: %s\n", argv_w2v_path);
		if(argv_force_timing_and_memory_to_output_path){
			free(argv_output_path_timing);
			free(argv_output_path_memory);
		}
		return 1;
	}

	err = measurement(
		&w2v,
		argv_jsonl_content_key,
		argv_input_path,
		argv_output_path,
		argv_output_path_timing,
		argv_output_path_memory,
		argv_target_column,
		argv_num_row_threads,
		argv_num_matrix_threads,
		argv_enable_stirling,
		argv_enable_ricotta_szeidl,
		argv_enable_pairwise,
		argv_enable_lexicographic,
		argv_enable_chao_et_al_functional_diversity,
		argv_enable_scheiner_species_phylogenetic_functional_diversity,
		argv_enable_leinster_cobbold_diversity,
		argv_enable_multithreaded_matrix_generation,
		argv_enable_timings,
		argv_enable_iterative_distance_computation,
		argv_enable_multithreaded_row_generation,
		argv_row_generation_batch_size,
		argv_enable_sentence_count_recompute_step,
		argv_sentence_recompute_step_use_log10,
		argv_enable_document_count_recompute_step,
		argv_document_recompute_step_use_log10,
		argv_enable_output_timing,
		argv_enable_output_memory,
		argv_sentence_count_recompute_step,
		argv_sentence_count_recompute_step_log10,
		argv_document_count_recompute_step,
		argv_document_count_recompute_step_log10,
		argv_enable_functional_evenness,
		argv_enable_functional_dispersion,
		argv_enable_functional_divergence_modified,
		argv_enable_non_disparity_functions,
		argv_enable_disparity_functions,
		argv_enable_shannon_weaver_entropy,
		argv_enable_good_entropy,
		argv_enable_renyi_entropy,
		argv_enable_patil_taillie_entropy,
		argv_enable_q_logarithmic_entropy,
		argv_enable_simpson_index,
		argv_enable_simpson_dominance_index,
		argv_enable_hill_number_standard,
		argv_enable_hill_evenness,
		argv_enable_berger_parker_index,
		argv_enable_junge1994_page22,
		argv_enable_brillouin_diversity,
		argv_enable_mcintosh_index,
		argv_enable_sw_entropy_over_log_n_species_pielou1975,
		argv_enable_sw_e_heip,
		argv_enable_sw_e_one_minus_d,
		argv_enable_sw_e_one_over_ln_d_williams1964,
		argv_enable_sw_e_minus_ln_d_pielou1977,
		argv_enable_sw_f_2_1_alatalo1981,
		argv_enable_sw_g_2_1_molinari1989,
		argv_enable_sw_e_bulla1994,
		argv_enable_sw_o_bulla1994,
		argv_enable_sw_e_mci_pielou1969,
		argv_enable_sw_e_prime_camargo1993,
		argv_enable_sw_e_var_smith_and_wilson1996_original,
		argv_stirling_alpha,
		argv_stirling_beta,
		argv_ricotta_szeidl_alpha,
		argv_chao_et_al_functional_diversity_alpha,
		argv_scheiner_species_phylogenetic_functional_diversity_alpha,
		argv_leinster_cobbold_diversity_alpha,
		argv_good_alpha,
		argv_good_beta,
		argv_renyi_alpha,
		argv_patil_taillie_alpha,
		argv_q_logarithmic_q,
		argv_hill_number_standard_alpha,
		argv_hill_evenness_alpha,
		argv_hill_evenness_beta
	);
	if(err != 0){
		perror("failed to call measurement\n");
		if(argv_force_timing_and_memory_to_output_path){
			free(argv_output_path_timing);
			free(argv_output_path_memory);
		}
		return 1;
	}

	free_word2vec(&w2v);

	if(argv_force_timing_and_memory_to_output_path){
		free(argv_output_path_timing);
		free(argv_output_path_memory);
	}
	return 0;

	malloc_fail:
	perror("malloc failed\n");
	return 1;
}

