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

#include "cfgparser/parser.h"

#define BFR_SIZE 256
#define MAX_FILES 1024

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

enum {
	CONLLU_COLUMN_UPOS,
	CONLLU_COLUMN_DEPREL
};
const int32_t CONLLU_COLUMNS_TO_ADD[2] = {CONLLU_COLUMN_UPOS, CONLLU_COLUMN_DEPREL};
// const int32_t NUM_CONLLU_COLUMNS_TO_ADD = (int32_t) (sizeof(CONLLU_COLUMNS_TO_ADD) / sizeof(int32_t));
const int32_t NUM_CONLLU_COLUMNS_TO_ADD = 0;

const int32_t CONLLU_ADD_FORM = 0;

int32_t apply_diversity_functions_to_graph_no_macros(struct graph* g, struct minimum_spanning_tree* mst, struct graph_distance_heap* heap, FILE* f_ptr, int64_t* previous_g_num_nodes_p, int64_t* num_sentences_p, int64_t* num_all_sentences_p, double* best_s, int8_t* mst_initialised, uint64_t i, struct sorted_array* sorted_array_discarded_because_not_in_vector_database, \
	const int16_t NUM_ROW_THREADS, \
	const int16_t NUM_MATRIX_THREADS, \
	const int8_t ENABLE_STIRLING, \
	const int8_t ENABLE_RICOTTA_SZEIDL, \
	const int8_t ENABLE_PAIRWISE, \
	const int8_t ENABLE_LEXICOGRAPHIC, \
	const int8_t ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY, \
	const int8_t ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY, \
	const int8_t ENABLE_LEINSTER_COBBOLD_DIVERSITY, \
	const int8_t ENABLE_FUNCTIONAL_EVENNESS, \
	const int8_t ENABLE_FUNCTIONAL_DISPERSION, \
	const int8_t ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED, \
	const int8_t ENABLE_NON_DISPARITY_FUNCTIONS, \
	const int8_t ENABLE_DISPARITY_FUNCTIONS, \
	const int8_t ENABLE_SHANNON_WEAVER_ENTROPY, \
	const int8_t ENABLE_GOOD_ENTROPY, \
	const int8_t ENABLE_RENYI_ENTROPY, \
	const int8_t ENABLE_PATIL_TAILLIE_ENTROPY, \
	const int8_t ENABLE_Q_LOGARITHMIC_ENTROPY, \
	const int8_t ENABLE_SIMPSON_INDEX, \
	const int8_t ENABLE_SIMPSON_DOMINANCE_INDEX, \
	const int8_t ENABLE_HILL_NUMBER_STANDARD, \
	const int8_t ENABLE_HILL_EVENNESS, \
	const int8_t ENABLE_BERGER_PARKER_INDEX, \
	const int8_t ENABLE_JUNGE1994_PAGE20, \
	const int8_t ENABLE_JUNGE1994_PAGE22, \
	const int8_t ENABLE_BRILLOUIN_DIVERSITY, \
	const int8_t ENABLE_MCINTOSH_INDEX, \
	const int8_t ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975, \
	const int8_t ENABLE_SW_E_HEIP, \
	const int8_t ENABLE_SW_E_ONE_MINUS_D, \
	const int8_t ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964, \
	const int8_t ENABLE_SW_E_MINUS_LN_D_PIELOU1977, \
	const int8_t ENABLE_SW_F_2_1_ALATALO1981, \
	const int8_t ENABLE_SW_G_2_1_MOLINARI1989, \
	const int8_t ENABLE_SW_E_BULLA1994, \
	const int8_t ENABLE_SW_O_BULLA1994, \
	const int8_t ENABLE_SW_E_MCI_PIELOU1969, \
	const int8_t ENABLE_SW_E_PRIME_CAMARGO1993, \
	const int8_t ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL, \
	const int8_t ENABLE_NHC_E_Q, \
	const int8_t ENABLE_TIMINGS, \
	const int8_t ENABLE_ITERATIVE_DISTANCE_COMPUTATION, \
	const int8_t ENABLE_MULTITHREADED_ROW_GENERATION, \
	const int8_t ENABLE_MULTITHREADED_MATRIX_GENERATION, \
	const char* W2V_PATH \
){
	int32_t err;
	if(ENABLE_ITERATIVE_DISTANCE_COMPUTATION){
		size_t local_malloc_size = g->num_nodes * sizeof(float);
		float* vector = (float*) malloc(local_malloc_size);
		if(vector == NULL){goto malloc_fail;}
		memset(vector, '\0', local_malloc_size);

		local_malloc_size = g->num_nodes * sizeof(uint8_t);
		uint8_t* used = (uint8_t*) malloc(local_malloc_size);
		if(used == NULL){goto malloc_fail;}
		memset(used, '\0', local_malloc_size);
			
		struct iterative_state_stirling_from_graph iter_state_stirling;
		struct iterative_state_pairwise_from_graph iter_state_pairwise;
		if(create_iterative_state_stirling_from_graph(&iter_state_stirling, g, STIRLING_ALPHA, STIRLING_BETA) != 0){
			perror("failled to call create_iterateive_state_stirling_from_graph\n");
			return 1;
		}
		if(create_iterative_state_pairwise_from_graph(&iter_state_pairwise, g) != 0){
			perror("failled to call create_iterateive_state_pairwise_from_graph\n");
			return 1;
		}

		int32_t i_index = 0;
		double sum = 0.0;
		for(uint64_t h = 0 ; h < g->num_nodes ; h++){
			if(ENABLE_MULTITHREADED_ROW_GENERATION){
				distance_row_from_graph_multithread(g, i_index, vector, NUM_ROW_THREADS);
			} else {
				distance_row_from_graph(g, i_index, vector);
			}
			for(uint64_t p = i_index + 1 ; p < g->num_nodes ; p++){
				sum += (double) vector[p];
			}

			if(ENABLE_STIRLING){iterate_iterative_state_stirling_from_graph(&iter_state_stirling, vector);}
			if(ENABLE_PAIRWISE){iterate_iterative_state_pairwise_from_graph(&iter_state_pairwise, vector);}

			i_index = h + 1;
			iter_state_stirling.i = i_index;
			iter_state_pairwise.i = i_index;
		}
		free(vector);
		free(used);

		// finalise_iterative_state_stirling_from_graph(&iter_state_stirling);
		finalise_iterative_state_pairwise_from_graph(&iter_state_pairwise);

		double mu_dist = sum / ((double) (g->num_nodes * (g->num_nodes - 1) / 2));

		fprintf(f_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%c", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes, mu_dist, '?'); // recomputing sigma dist would be expensive

		if(ENABLE_STIRLING){printf("[log] [end iter] stirling: %f\n", iter_state_stirling.result);}
		fprintf(f_ptr, "\t%.10e", iter_state_stirling.result);
		if(ENABLE_PAIRWISE){printf("[log] [end iter] pairwise: %f\n", iter_state_pairwise.result);}
		fprintf(f_ptr, "\t%.10e", iter_state_pairwise.result);

		fprintf(f_ptr, "\n");


	} else {
		struct matrix m_mst;
		if(ENABLE_DISPARITY_FUNCTIONS && ENABLE_FUNCTIONAL_EVENNESS){
			if(create_matrix(&m_mst, g->num_nodes, g->num_nodes, FP64) != 0){
				perror("failed to call create_matrix for MST\n");
				return 1;
			}
			memset(m_mst.active, '\0', g->num_nodes * g->num_nodes * sizeof(uint8_t));
			memset(m_mst.active_final, '\0', g->num_nodes * g->num_nodes * sizeof(uint8_t));
		}

		struct matrix m;
		if(ENABLE_DISPARITY_FUNCTIONS){
			if(create_matrix(&m, (uint32_t) g->num_nodes, (uint32_t) g->num_nodes, FP32) != 0){
				perror("failed to call create_matrix\n");
				return 1;
			}
			time_t t = time(NULL);
			if(ENABLE_MULTITHREADED_MATRIX_GENERATION){
				if(distance_matrix_from_graph_multithread(g, &m, NUM_MATRIX_THREADS) != 0){
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
			time_t delta_t = time(NULL) - t;
			if(ENABLE_TIMINGS){
				printf("[log] [time] Computed matrix in %lis\n", delta_t);
			}

			if((*mst_initialised)){
				free_graph_distance_heap(heap);
				if(ENABLE_FUNCTIONAL_EVENNESS){
					free_minimum_spanning_tree(mst);
				}
			}
		}


		struct graph_distance_heap local_heap;

		double mu_dist;
		double sigma_dist;
		float mu_dist_fp32;
		float sigma_dist_fp32;
		if(ENABLE_DISPARITY_FUNCTIONS && ENABLE_FUNCTIONAL_EVENNESS){
			err = create_graph_distance_heap(&local_heap, g, GRAPH_NODE_FP32, &m);
			if(err != 0){
				perror("failed to call create_graph_distance_heap\n");
				return EXIT_FAILURE;
			}
	
			(*heap) = local_heap;
		}

		if(ENABLE_DISPARITY_FUNCTIONS){
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
	
			fprintf(f_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%.10e\t%.10e", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes, mu_dist, sigma_dist);
		} else {
			fprintf(f_ptr, "%lu\t%li\t%li\t%s\t%li\t%.10e\t%lu\t%s\t%s", i+1, (*num_sentences_p), (*num_all_sentences_p), W2V_PATH, sorted_array_discarded_because_not_in_vector_database->num_elements, (*best_s), g->num_nodes, "?", "?");
		}

		if(ENABLE_DISPARITY_FUNCTIONS && ENABLE_FUNCTIONAL_EVENNESS){
			time_t t = time(NULL);
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

			time_t delta_t = time(NULL) - t;
			if(ENABLE_TIMINGS){
				printf("[log] [time] Computed MST in %lis\n", delta_t);
			}
		}
		(*previous_g_num_nodes_p) = g->num_nodes;

		if(ENABLE_DISPARITY_FUNCTIONS){
			time_t t;
			time_t delta_t;
			if(ENABLE_STIRLING){
				double stirling;
				t = time(NULL);
				err = stirling_from_graph(g, &stirling, STIRLING_ALPHA, STIRLING_BETA, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call stirling_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed Stirling in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", stirling);
			}
	
			if(ENABLE_RICOTTA_SZEIDL){
				double ricotta_szeidl;
				t = time(NULL);
				err = ricotta_szeidl_from_graph(g, &ricotta_szeidl, RICOTTA_SZEIDL_ALPHA, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call ricotta_szeidl_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed Ricotta-Szeidl in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", ricotta_szeidl);
			}
	
			if(ENABLE_PAIRWISE){
				double pairwise;
				t = time(NULL);
				err = pairwise_from_graph(g, &pairwise, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call pairwise_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed pairwise in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", pairwise);
			}
	
			if(ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY){
				double chao_et_al_functional_diversity;
				double chao_et_al_functional_hill_number;
				t = time(NULL);
				err = chao_et_al_functional_diversity_from_graph(g, &chao_et_al_functional_diversity, &chao_et_al_functional_hill_number, CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call chao_et_al_functional_diversity_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed Chao et al. in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e\t%.10e", chao_et_al_functional_diversity, chao_et_al_functional_hill_number);
			}
	
			if(ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY){
				double scheiner_species_phylogenetic_functional_diversity;
				double scheiner_species_phylogenetic_functional_hill_number;
				t = time(NULL);
				err = scheiner_species_phylogenetic_functional_diversity_from_graph(g, &scheiner_species_phylogenetic_functional_diversity, &scheiner_species_phylogenetic_functional_hill_number, SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call scheiner_species_phylogenetic_functional_diversity_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed Scheiner in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e\t%.10e", scheiner_species_phylogenetic_functional_diversity, scheiner_species_phylogenetic_functional_hill_number);
			}
	
			if(ENABLE_LEINSTER_COBBOLD_DIVERSITY){
				double leinster_cobbold_diversity;
				double leinster_cobbold_hill_number;
				t = time(NULL);
				err = leinster_cobbold_diversity_from_graph(g, &leinster_cobbold_diversity, &leinster_cobbold_hill_number, LEINSTER_COBBOLD_DIVERSITY_ALPHA, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call leinster_cobbold_diversity_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed Leinster-Cobbold in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e\t%.10e", leinster_cobbold_diversity, leinster_cobbold_hill_number);
			}
	
			if(ENABLE_LEXICOGRAPHIC){
				double lexicographic;
				long double lexicographic_hybrid_scheiner;
				t = time(NULL);
				err = lexicographic_from_graph(g, &lexicographic, &lexicographic_hybrid_scheiner, GRAPH_NODE_FP32, &m);
				if(err != 0){
					perror("failed to call lexicographic_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed lexicographic in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e\t%.10Le", lexicographic, lexicographic_hybrid_scheiner);
			}
	
			if(ENABLE_FUNCTIONAL_EVENNESS){
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
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed FEve in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", functional_evenness);
			}
	
			if(ENABLE_FUNCTIONAL_DISPERSION){
				double functional_dispersion;
				t = time(NULL);
				err = functional_dispersion_from_graph(g, &functional_dispersion, GRAPH_NODE_FP32);
				if(err != 0){
					perror("failed to call functional_dispersion_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed FDisp in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", functional_dispersion);
			}
	
			if(ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED){
				double functional_divergence_modified;
				t = time(NULL);
				err = functional_divergence_modified_from_graph(g, &functional_divergence_modified, GRAPH_NODE_FP32);
				if(err != 0){
					perror("failed to call functional_divergence_modified_from_graph\n");
					return EXIT_FAILURE;
				}
				delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){
					printf("[log] [time] Computed FDiv in %lis\n", delta_t);
				}
				fprintf(f_ptr, "\t%.10e", functional_divergence_modified);
			}
		}

		if(ENABLE_NON_DISPARITY_FUNCTIONS){
			if(ENABLE_SHANNON_WEAVER_ENTROPY){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				shannon_weaver_entropy_from_graph(g, &res_entropy, &res_hill_number);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed SW entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
			}
			if(ENABLE_GOOD_ENTROPY){
				double res;
				time_t t = time(NULL);
				good_entropy_from_graph(g, &res, GOOD_ALPHA, GOOD_BETA);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Good entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_RENYI_ENTROPY){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				renyi_entropy_from_graph(g, &res_entropy, &res_hill_number, RENYI_ALPHA);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Renyi entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
			}
			if(ENABLE_PATIL_TAILLIE_ENTROPY){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				patil_taillie_entropy_from_graph(g, &res_entropy, &res_hill_number, PATIL_TAILLIE_ALPHA);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Patil-Taillie entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
			}
			if(ENABLE_Q_LOGARITHMIC_ENTROPY){
				double res_entropy;
				double res_hill_number;
				time_t t = time(NULL);
				q_logarithmic_entropy_from_graph(g, &res_entropy, &res_hill_number, Q_LOGARITHMIC_Q);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed q-logarithmic entropy in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_entropy, res_hill_number);
			}
			if(ENABLE_SIMPSON_INDEX){
				double res;
				time_t t = time(NULL);
				simpson_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Simpson index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SIMPSON_DOMINANCE_INDEX){
				double res;
				time_t t = time(NULL);
				simpson_dominance_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Simpson dominance index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_HILL_NUMBER_STANDARD){
				double res;
				time_t t = time(NULL);
				hill_number_standard_from_graph(g, &res, HILL_NUMBER_STANDARD_ALPHA);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Hill number (standard) in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_HILL_EVENNESS){
				double res;
				time_t t = time(NULL);
				hill_evenness_from_graph(g, &res, HILL_EVENNESS_ALPHA, HILL_EVENNESS_BETA);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Hill evenness in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_BERGER_PARKER_INDEX){
				double res;
				time_t t = time(NULL);
				berger_parker_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Berger Parker index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_JUNGE1994_PAGE22){
				double res;
				time_t t = time(NULL);
				junge1994_page22_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Junge 1994 p22 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_BRILLOUIN_DIVERSITY){
				double res;
				time_t t = time(NULL);
				brillouin_diversity_from_graph(g, &res);
				// printf("res: %f\n", res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed Brillouin diversity in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_MCINTOSH_INDEX){
				double res;
				time_t t = time(NULL);
				mcintosh_index_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed McIntosh index in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975){
				double res;
				time_t t = time(NULL);
				sw_entropy_over_log_n_species_pielou1975_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) entropy over log n species Pielou 1975 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_HEIP){
				double res;
				time_t t = time(NULL);
				sw_e_heip_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E Heip in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_ONE_MINUS_D){
				double res;
				time_t t = time(NULL);
				sw_e_one_minus_D_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E one minus D in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964){
				double res;
				time_t t = time(NULL);
				sw_e_one_over_D_williams1964_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E one over ln D Williams 1964 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_MINUS_LN_D_PIELOU1977){
				double res;
				time_t t = time(NULL);
				sw_e_minus_ln_D_pielou1977_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E minus ln D Pielou 1977 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_F_2_1_ALATALO1981){
				double res;
				time_t t = time(NULL);
				sw_f_2_1_alatalo1981_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) F_2_1 Alatalo 1981 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_G_2_1_MOLINARI1989){
				double res;
				time_t t = time(NULL);
				sw_g_2_1_molinari1989_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) G_2_1 Molinari 1989 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_BULLA1994){
				double res;
				time_t t = time(NULL);
				sw_e_bulla1994_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E Bulla 1994 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_O_BULLA1994){
				double res;
				time_t t = time(NULL);
				sw_o_bulla1994_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) O bulla 1994 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_MCI_PIELOU1969){
				double res;
				time_t t = time(NULL);
				sw_e_mci_pielou1969_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E MCI Pielou 1969 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_PRIME_CAMARGO1993){
				double res;
				time_t t = time(NULL);
				sw_e_prime_camargo1993_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E prime Camargo 1993 in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL){
				double res;
				time_t t = time(NULL);
				sw_e_var_smith_and_wilson1996_original_from_graph(g, &res);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed (SW) E var Smith and Wilson 1996 original in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e", res);
			}
			if(ENABLE_NHC_E_Q){
				double res_nhc;
				double res_e_q;
				time_t t = time(NULL);
				nhc_e_q_from_graph(g, &res_nhc, &res_e_q);
				time_t delta_t = time(NULL) - t;
				if(ENABLE_TIMINGS){printf("[log] [time] Computed NHC_E_Q in %lis\n", delta_t);}
				fprintf(f_ptr, "\t%.10e\t%.10e", res_nhc, res_e_q);
			}
		}

		fprintf(f_ptr, "\n");

		if(ENABLE_DISPARITY_FUNCTIONS){
			free_matrix(&m);
			if(ENABLE_FUNCTIONAL_EVENNESS){
				free_matrix(&m_mst);
			}
		}
	}

	return 0;

	malloc_fail:
	perror("malloc failed\n");
	return 1;
}

int32_t measurement_from_cfg(const struct cfg* const config){
	int32_t UD_COLUMN;
	char* value_ud_column = cfg_get_value(config, "TARGET_COLUMN");
	if(value_ud_column == NULL){UD_COLUMN = UD_FORM;}
	else if(strcmp(value_ud_column, "UD_FORM") == 0){UD_COLUMN = UD_FORM;}
	else if(strcmp(value_ud_column, "UD_MWE") == 0){UD_COLUMN = UD_MWE;}
	else {fprintf(stderr, "Unknown UD_COLUMN: %s\n", value_ud_column); return 1;}


	const int8_t ENABLE_STIRLING = (strcmp(cfg_get_value(config, "ENABLE_STIRLING"), "1") == 0);
	const int8_t ENABLE_RICOTTA_SZEIDL = (strcmp(cfg_get_value(config, "ENABLE_RICOTTA_SZEIDL"), "1") == 0);
	const int8_t ENABLE_PAIRWISE = (strcmp(cfg_get_value(config, "ENABLE_PAIRWISE"), "1") == 0);
	const int8_t ENABLE_LEXICOGRAPHIC = (strcmp(cfg_get_value(config, "ENABLE_LEXICOGRAPHIC"), "1") == 0);
	const int8_t ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY = (strcmp(cfg_get_value(config, "ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY"), "1") == 0);
	const int8_t ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY = (strcmp(cfg_get_value(config, "ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY"), "1") == 0);
	const int8_t ENABLE_LEINSTER_COBBOLD_DIVERSITY = (strcmp(cfg_get_value(config, "ENABLE_LEINSTER_COBBOLD_DIVERSITY"), "1") == 0);
	const int8_t ENABLE_FUNCTIONAL_EVENNESS = (strcmp(cfg_get_value(config, "ENABLE_FUNCTIONAL_EVENNESS"), "1") == 0);
	const int8_t ENABLE_FUNCTIONAL_DISPERSION = (strcmp(cfg_get_value(config, "ENABLE_FUNCTIONAL_DISPERSION"), "1") == 0);
	const int8_t ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED = (strcmp(cfg_get_value(config, "ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED"), "1") == 0);
	const int8_t ENABLE_NON_DISPARITY_FUNCTIONS = (strcmp(cfg_get_value(config, "ENABLE_NON_DISPARITY_FUNCTIONS"), "1") == 0);
	const int8_t ENABLE_DISPARITY_FUNCTIONS = (strcmp(cfg_get_value(config, "ENABLE_DISPARITY_FUNCTIONS"), "1") == 0);
	const int8_t ENABLE_SHANNON_WEAVER_ENTROPY = (strcmp(cfg_get_value(config, "ENABLE_SHANNON_WEAVER_ENTROPY"), "1") == 0);
	const int8_t ENABLE_GOOD_ENTROPY = (strcmp(cfg_get_value(config, "ENABLE_GOOD_ENTROPY"), "1") == 0);
	const int8_t ENABLE_RENYI_ENTROPY = (strcmp(cfg_get_value(config, "ENABLE_RENYI_ENTROPY"), "1") == 0);
	const int8_t ENABLE_PATIL_TAILLIE_ENTROPY = (strcmp(cfg_get_value(config, "ENABLE_PATIL_TAILLIE_ENTROPY"), "1") == 0);
	const int8_t ENABLE_Q_LOGARITHMIC_ENTROPY = (strcmp(cfg_get_value(config, "ENABLE_Q_LOGARITHMIC_ENTROPY"), "1") == 0);
	const int8_t ENABLE_SIMPSON_INDEX = (strcmp(cfg_get_value(config, "ENABLE_SIMPSON_INDEX"), "1") == 0);
	const int8_t ENABLE_SIMPSON_DOMINANCE_INDEX = (strcmp(cfg_get_value(config, "ENABLE_SIMPSON_DOMINANCE_INDEX"), "1") == 0);
	const int8_t ENABLE_HILL_NUMBER_STANDARD = (strcmp(cfg_get_value(config, "ENABLE_HILL_NUMBER_STANDARD"), "1") == 0);
	const int8_t ENABLE_HILL_EVENNESS = (strcmp(cfg_get_value(config, "ENABLE_HILL_EVENNESS"), "1") == 0);
	const int8_t ENABLE_BERGER_PARKER_INDEX = (strcmp(cfg_get_value(config, "ENABLE_BERGER_PARKER_INDEX"), "1") == 0);
	const int8_t ENABLE_JUNGE1994_PAGE20 = (strcmp(cfg_get_value(config, "ENABLE_JUNGE1994_PAGE20"), "1") == 0);
	const int8_t ENABLE_JUNGE1994_PAGE22 = (strcmp(cfg_get_value(config, "ENABLE_JUNGE1994_PAGE22"), "1") == 0);
	const int8_t ENABLE_BRILLOUIN_DIVERSITY = (strcmp(cfg_get_value(config, "ENABLE_BRILLOUIN_DIVERSITY"), "1") == 0);
	const int8_t ENABLE_MCINTOSH_INDEX = (strcmp(cfg_get_value(config, "ENABLE_MCINTOSH_INDEX"), "1") == 0);
	const int8_t ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975 = (strcmp(cfg_get_value(config, "ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975"), "1") == 0);
	const int8_t ENABLE_SW_E_HEIP = (strcmp(cfg_get_value(config, "ENABLE_SW_E_HEIP"), "1") == 0);
	const int8_t ENABLE_SW_E_ONE_MINUS_D = (strcmp(cfg_get_value(config, "ENABLE_SW_E_ONE_MINUS_D"), "1") == 0);
	const int8_t ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964 = (strcmp(cfg_get_value(config, "ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964"), "1") == 0);
	const int8_t ENABLE_SW_E_MINUS_LN_D_PIELOU1977 = (strcmp(cfg_get_value(config, "ENABLE_SW_E_MINUS_LN_D_PIELOU1977"), "1") == 0);
	const int8_t ENABLE_SW_F_2_1_ALATALO1981 = (strcmp(cfg_get_value(config, "ENABLE_SW_F_2_1_ALATALO1981"), "1") == 0);
	const int8_t ENABLE_SW_G_2_1_MOLINARI1989 = (strcmp(cfg_get_value(config, "ENABLE_SW_G_2_1_MOLINARI1989"), "1") == 0);
	const int8_t ENABLE_SW_E_BULLA1994 = (strcmp(cfg_get_value(config, "ENABLE_SW_E_BULLA1994"), "1") == 0);
	const int8_t ENABLE_SW_O_BULLA1994 = (strcmp(cfg_get_value(config, "ENABLE_SW_O_BULLA1994"), "1") == 0);
	const int8_t ENABLE_SW_E_MCI_PIELOU1969 = (strcmp(cfg_get_value(config, "ENABLE_SW_E_MCI_PIELOU1969"), "1") == 0);
	const int8_t ENABLE_SW_E_PRIME_CAMARGO1993 = (strcmp(cfg_get_value(config, "ENABLE_SW_E_PRIME_CAMARGO1993"), "1") == 0);
	const int8_t ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL = (strcmp(cfg_get_value(config, "ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL"), "1") == 0);
	const int8_t ENABLE_NHC_E_Q = (strcmp(cfg_get_value(config, "ENABLE_NHC_E_Q"), "1") == 0);

	const int8_t ENABLE_TIMINGS = (strcmp(cfg_get_value(config, "ENABLE_TIMINGS"), "1") == 0);
	const int8_t ENABLE_ITERATIVE_DISTANCE_COMPUTATION = (strcmp(cfg_get_value(config, "ENABLE_ITERATIVE_DISTANCE_COMPUTATION"), "1") == 0);
	const int8_t ENABLE_MULTITHREADED_ROW_GENERATION = (strcmp(cfg_get_value(config, "ENABLE_MULTITHREADED_ROW_GENERATION"), "1") == 0);
	const int8_t ENABLE_MULTITHREADED_MATRIX_GENERATION = (strcmp(cfg_get_value(config, "ENABLE_MULTITHREADED_MATRIX_GENERATION"), "1") == 0);

	const int8_t ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP = (strcmp(cfg_get_value(config, "ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP"), "1") == 0);
	const int8_t ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP = (strcmp(cfg_get_value(config, "ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP"), "1") == 0);

	const int32_t SENTENCE_COUNT_RECOMPUTE_STEP = (int32_t) strtol(cfg_get_value(config, "SENTENCE_COUNT_RECOMPUTE_STEP"), NULL, 10);
	const int32_t DOCUMENT_COUNT_RECOMPUTE_STEP = (int32_t) strtol(cfg_get_value(config, "DOCUMENT_COUNT_RECOMPUTE_STEP"), NULL, 10);
	const int16_t NUM_ROW_THREADS = (int16_t) strtol(cfg_get_value(config, "NUM_ROW_THREADS"), NULL, 10);
	const int16_t NUM_MATRIX_THREADS = (int16_t) strtol(cfg_get_value(config, "NUM_MATRIX_THREADS"), NULL, 10);

	const char* INPUT_PATH = cfg_get_value(config, "INPUT_PATH");
	const char* OUTPUT_PATH = cfg_get_value(config, "OUTPUT_PATH");
	const char* W2V_PATH = cfg_get_value(config, "W2V_PATH");

	const char* JSONL_CONTENT_KEY = cfg_get_value(config, "JSONL_CONTENT_KEY");




	struct word2vec w2v;
	int32_t err = load_word2vec_binary(&w2v, W2V_PATH);
	if(err != 0){
		perror("failed to call load_word2vec binary\n");
		return 1;
	}



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

	for(uint64_t i = 0 ; i < w2v.num_vectors ; i++){
		w2v.keys[i].active_in_current_graph = 0;
		w2v.keys[i].num_occurrences = 0;
		w2v.keys[i].graph_node_pointer = NULL; // !
	}

	FILE* f_paths_ptr = fopen(INPUT_PATH, "r");
	if(f_paths_ptr == NULL){
		fprintf(stderr, "cannot open %s\n", INPUT_PATH);
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
	fclose(f_paths_ptr);

	int32_t num_input_paths = num_files;
	char** input_paths = bfr;
	char** input_paths_true_positives = NULL;

	FILE* f_ptr = fopen(OUTPUT_PATH, "w");
	if(f_ptr == NULL){
		perror("failed to open file\n");
		return EXIT_FAILURE;
	}
	fprintf(f_ptr, "num_active_files\tnum_active_sentences\tnum_all_sentences\tw2v\tnum_discarded_types\ts\tn\tmu_dist\tsigma_dist");
	if(ENABLE_DISPARITY_FUNCTIONS){
		if(ENABLE_STIRLING){fprintf(f_ptr, "\tstirling_alpha%.10e_beta%.10e", STIRLING_ALPHA, STIRLING_BETA);}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_RICOTTA_SZEIDL){fprintf(f_ptr, "\tricotta_szeidl_alpha%.10e", RICOTTA_SZEIDL_ALPHA);}
		if(ENABLE_PAIRWISE){fprintf(f_ptr, "\tpairwise");}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY){fprintf(f_ptr, "\tchao_et_al_functional_diversity_alpha%.10e\tchao_et_al_functional_hill_number_alpha%.10e", CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA, CHAO_ET_AL_FUNCTIONAL_DIVERSITY_ALPHA);}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY){fprintf(f_ptr, "\tscheiner_species_phylogenetic_functional_diversity_alpha%.10e\tscheiner_species_phylogenetic_functional_hill_number_alpha%.10e", SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA, SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY_ALPHA);}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_LEINSTER_COBBOLD_DIVERSITY){fprintf(f_ptr, "\tleinster_cobbold_diversity_alpha%.10e\tleinster_cobbold_hill_number_alpha%.10e", LEINSTER_COBBOLD_DIVERSITY_ALPHA, LEINSTER_COBBOLD_DIVERSITY_ALPHA);}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_LEXICOGRAPHIC){fprintf(f_ptr, "\tlexicographic\tlexicographic_hybrid_scheiner");}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_FUNCTIONAL_EVENNESS){fprintf(f_ptr, "\tfunctional_evenness");}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_FUNCTIONAL_DISPERSION){fprintf(f_ptr, "\tfunctional_dispersion");}
		if(!ENABLE_ITERATIVE_DISTANCE_COMPUTATION && ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED){fprintf(f_ptr, "\tfunctional_divergence_modified");}
	}


	if(ENABLE_NON_DISPARITY_FUNCTIONS){
		if(ENABLE_SHANNON_WEAVER_ENTROPY){fprintf(f_ptr, "\tshannon_weaver_entropy\tshannon_weaver_hill_number");}
		if(ENABLE_GOOD_ENTROPY){fprintf(f_ptr, "\tgood_entropy_alpha%.4e_beta%.4e", GOOD_ALPHA, GOOD_BETA);}
		if(ENABLE_RENYI_ENTROPY){fprintf(f_ptr, "\trenyi_entropy_alpha%.4e\trenyi_hill_number_alpha%.4e", RENYI_ALPHA, RENYI_ALPHA);}
		if(ENABLE_PATIL_TAILLIE_ENTROPY){fprintf(f_ptr, "\tpatil_taillie_entropy_alpha%.4e\tpatil_taillie_hill_number_alpha%.4e", PATIL_TAILLIE_ALPHA, PATIL_TAILLIE_ALPHA);}
		if(ENABLE_Q_LOGARITHMIC_ENTROPY){fprintf(f_ptr, "\tq_logarithmic_entropy_alpha%.4e\tq_logarithmic_hill_number_alpha%.4e", Q_LOGARITHMIC_Q, Q_LOGARITHMIC_Q);}
		if(ENABLE_SIMPSON_INDEX){fprintf(f_ptr, "\tsimpson_index");}
		if(ENABLE_SIMPSON_DOMINANCE_INDEX){fprintf(f_ptr, "\tsimpson_dominance_index");}
		if(ENABLE_HILL_NUMBER_STANDARD){fprintf(f_ptr, "\thill_number_standard_alpha%.4e", HILL_NUMBER_STANDARD_ALPHA);}
		if(ENABLE_HILL_EVENNESS){fprintf(f_ptr, "\thill_evenness_alpha%.4e_beta%.4e", HILL_EVENNESS_ALPHA, HILL_EVENNESS_BETA);}
		if(ENABLE_BERGER_PARKER_INDEX){fprintf(f_ptr, "\tberger_parker_index");}
		if(ENABLE_JUNGE1994_PAGE22){fprintf(f_ptr, "\tjunge1994_page22");}
		if(ENABLE_BRILLOUIN_DIVERSITY){fprintf(f_ptr, "\tbrillouin_diversity");}
		if(ENABLE_MCINTOSH_INDEX){fprintf(f_ptr, "\tmcintosh_index");}
		if(ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975){fprintf(f_ptr, "\tsw_entropy_over_log_n_species_pielou1975");}
		if(ENABLE_SW_E_HEIP){fprintf(f_ptr, "\tsw_e_heip");}
		if(ENABLE_SW_E_ONE_MINUS_D){fprintf(f_ptr, "\te_one_minus_D");}
		if(ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964){fprintf(f_ptr, "\te_one_over_ln_D_williams1964");}
		if(ENABLE_SW_E_MINUS_LN_D_PIELOU1977){fprintf(f_ptr, "\te_minus_ln_D_pielou1977");}
		if(ENABLE_SW_F_2_1_ALATALO1981){fprintf(f_ptr, "\tsw_f_2_1_alatalo1981");}
		if(ENABLE_SW_G_2_1_MOLINARI1989){fprintf(f_ptr, "\tsw_g_2_1_molinari1989");}
		if(ENABLE_SW_E_BULLA1994){fprintf(f_ptr, "\tsw_e_bulla1994");}
		if(ENABLE_SW_O_BULLA1994){fprintf(f_ptr, "\tsw_o_bulla1994");}
		if(ENABLE_SW_E_MCI_PIELOU1969){fprintf(f_ptr, "\tsw_e_mci_pielou1969");}
		if(ENABLE_SW_E_PRIME_CAMARGO1993){fprintf(f_ptr, "\tsw_e_prime_camargo1993");}
		if(ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL){fprintf(f_ptr, "\tsw_e_var_smith_and_wilson1996_original");}
		if(ENABLE_NHC_E_Q){fprintf(f_ptr, "\tNHC\tE_Q");}
	}

	fprintf(f_ptr, "\n");


	int32_t num_nodes = 0;
	int64_t num_sentences = 0;
	int64_t max_num_sentences = 10000;
	int8_t use_max_num_sentences = 0;

	int64_t num_documents = 0;
	int64_t max_num_documents = 10000;
	int8_t use_max_num_documents = 0;

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
		if(input_paths_true_positives == NULL){
			printf("processing %s\n", input_paths[i]);
		} else {
			printf("processing %s (true positives: %s)\n", input_paths[i], input_paths_true_positives[i]);
		}
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
			if(create_jsonl_document_iterator(&jdi, input_paths[i], JSONL_CONTENT_KEY) != 0){
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
					if(UD_COLUMN == UD_MWE){
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
						switch(UD_COLUMN){
							case UD_FORM:
								index = word2vec_key_to_index(&w2v, csi.current_sentence.tokens[j].form);
								index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, csi.current_sentence.tokens[j].form);
								len = strlen(csi.current_sentence.tokens[j].form);
								if(len > key_size - 1){
									len = key_size - 1;
								}
								memcpy(key, csi.current_sentence.tokens[j].form, len);
								break;
							case UD_LEMMA:
								index = word2vec_key_to_index(&w2v, csi.current_sentence.tokens[j].lemma);
								index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, csi.current_sentence.tokens[j].lemma);
								len = strlen(csi.current_sentence.tokens[j].lemma);
								if(len > key_size - 1){
									len = key_size - 1;
								}
								memcpy(key, csi.current_sentence.tokens[j].lemma, len);
								break;
							default:
								index = -1;
								perror("UD_COLUMN not properly defined\n");
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
	
						index = word2vec_key_to_index(&w2v, key);
						index_in_discarded = key_to_index_sorted_array(&sorted_array_discarded_because_not_in_vector_database, key);
					}
	
					if(index != -1){
						if(w2v.keys[index].active_in_current_graph == 0){
							w2v.keys[index].active_in_current_graph = 1;
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
							local_node.word2vec_entry_pointer = &(w2v.keys[index]);
							local_node.vector.fp32 = w2v.keys[index].vector;
	
							local_node.num_dimensions = w2v.num_dimensions;
							local_node.already_considered = 0;
							local_node.relative_proportion = 1.0;
							local_node.absolute_proportion = 1;
							g.nodes[g.num_nodes] = local_node;
							w2v.keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
							w2v.keys[index].graph_node_index = g.num_nodes;
							zipfian_n++;
							g.num_nodes++;
						} else {
							g.nodes[w2v.keys[index].graph_node_index].absolute_proportion++;
						}
						w2v.keys[index].num_occurrences++; // ? mutex ?
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
				if(UD_COLUMN == UD_MWE){
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
		
						int32_t index = word2vec_key_to_index(&w2v, bfr);
						if(index != -1){
							found_at_least_one_mwe = 1;
							if(w2v.keys[index].active_in_current_graph == 0){
								w2v.keys[index].active_in_current_graph = 1;
								if(g.num_nodes == g.capacity){
									err = request_more_capacity_graph(&g);
									if(err != 0){
										perror("failed to call request_more_capacity_graph\n");
										return 1;
									}
								}
								struct graph_node local_graph_node;
								local_graph_node.num_dimensions = (int16_t) w2v.num_dimensions;
								local_graph_node.already_considered = 0;
								local_graph_node.relative_proportion = 1.0;
								local_graph_node.absolute_proportion = 1;
								local_graph_node.vector.fp32 = w2v.keys[index].vector;
								local_graph_node.word2vec_entry_pointer = &(w2v.keys[index]);
								g.nodes[g.num_nodes] = local_graph_node;
								w2v.keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
								w2v.keys[index].graph_node_index = g.num_nodes;
								zipfian_n++;
								g.num_nodes++;
								num_nodes++;
							} else {
								g.nodes[w2v.keys[index].graph_node_index].absolute_proportion++;
							}
							w2v.keys[index].num_occurrences++;
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
				if(UD_COLUMN != UD_MWE){
					num_sentences = num_all_sentences;
				}	

				// sentence level recomputation
				if((UD_COLUMN != UD_MWE || found_at_least_one_mwe) && ((!ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP) || num_sentences % SENTENCE_COUNT_RECOMPUTE_STEP == 0) && g.num_nodes > 1){
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
						// err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database);
						err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database, 
							NUM_ROW_THREADS,
							NUM_MATRIX_THREADS,
							ENABLE_STIRLING,
							ENABLE_RICOTTA_SZEIDL,
							ENABLE_PAIRWISE,
							ENABLE_LEXICOGRAPHIC,
							ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY,
							ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY,
							ENABLE_LEINSTER_COBBOLD_DIVERSITY,
							ENABLE_FUNCTIONAL_EVENNESS,
							ENABLE_FUNCTIONAL_DISPERSION,
							ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED,
							ENABLE_NON_DISPARITY_FUNCTIONS,
							ENABLE_DISPARITY_FUNCTIONS,
							ENABLE_SHANNON_WEAVER_ENTROPY,
							ENABLE_GOOD_ENTROPY,
							ENABLE_RENYI_ENTROPY,
							ENABLE_PATIL_TAILLIE_ENTROPY,
							ENABLE_Q_LOGARITHMIC_ENTROPY,
							ENABLE_SIMPSON_INDEX,
							ENABLE_SIMPSON_DOMINANCE_INDEX,
							ENABLE_HILL_NUMBER_STANDARD,
							ENABLE_HILL_EVENNESS,
							ENABLE_BERGER_PARKER_INDEX,
							ENABLE_JUNGE1994_PAGE20,
							ENABLE_JUNGE1994_PAGE22,
							ENABLE_BRILLOUIN_DIVERSITY,
							ENABLE_MCINTOSH_INDEX,
							ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975,
							ENABLE_SW_E_HEIP,
							ENABLE_SW_E_ONE_MINUS_D,
							ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964,
							ENABLE_SW_E_MINUS_LN_D_PIELOU1977,
							ENABLE_SW_F_2_1_ALATALO1981,
							ENABLE_SW_G_2_1_MOLINARI1989,
							ENABLE_SW_E_BULLA1994,
							ENABLE_SW_O_BULLA1994,
							ENABLE_SW_E_MCI_PIELOU1969,
							ENABLE_SW_E_PRIME_CAMARGO1993,
							ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL,
							ENABLE_NHC_E_Q,
							ENABLE_TIMINGS,
							ENABLE_ITERATIVE_DISTANCE_COMPUTATION,
							ENABLE_MULTITHREADED_ROW_GENERATION,
							ENABLE_MULTITHREADED_MATRIX_GENERATION,
							W2V_PATH
						);
						if(err != 0){
							perror("failed to call apply_diversity_functions_to_graph_no_macros\n");
							return 1;
						}
						previous_best_s = best_s;
					} else { // end of comparison best_s?
						printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
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

			if((UD_COLUMN != UD_MWE || found_at_least_one_mwe) && ((!ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP) || num_documents % DOCUMENT_COUNT_RECOMPUTE_STEP == 0) && g.num_nodes > 1){
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
					// err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database);
					err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database,
						NUM_ROW_THREADS,
						NUM_MATRIX_THREADS,
						ENABLE_STIRLING,
						ENABLE_RICOTTA_SZEIDL,
						ENABLE_PAIRWISE,
						ENABLE_LEXICOGRAPHIC,
						ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY,
						ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY,
						ENABLE_LEINSTER_COBBOLD_DIVERSITY,
						ENABLE_FUNCTIONAL_EVENNESS,
						ENABLE_FUNCTIONAL_DISPERSION,
						ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED,
						ENABLE_NON_DISPARITY_FUNCTIONS,
						ENABLE_DISPARITY_FUNCTIONS,
						ENABLE_SHANNON_WEAVER_ENTROPY,
						ENABLE_GOOD_ENTROPY,
						ENABLE_RENYI_ENTROPY,
						ENABLE_PATIL_TAILLIE_ENTROPY,
						ENABLE_Q_LOGARITHMIC_ENTROPY,
						ENABLE_SIMPSON_INDEX,
						ENABLE_SIMPSON_DOMINANCE_INDEX,
						ENABLE_HILL_NUMBER_STANDARD,
						ENABLE_HILL_EVENNESS,
						ENABLE_BERGER_PARKER_INDEX,
						ENABLE_JUNGE1994_PAGE20,
						ENABLE_JUNGE1994_PAGE22,
						ENABLE_BRILLOUIN_DIVERSITY,
						ENABLE_MCINTOSH_INDEX,
						ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975,
						ENABLE_SW_E_HEIP,
						ENABLE_SW_E_ONE_MINUS_D,
						ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964,
						ENABLE_SW_E_MINUS_LN_D_PIELOU1977,
						ENABLE_SW_F_2_1_ALATALO1981,
						ENABLE_SW_G_2_1_MOLINARI1989,
						ENABLE_SW_E_BULLA1994,
						ENABLE_SW_O_BULLA1994,
						ENABLE_SW_E_MCI_PIELOU1969,
						ENABLE_SW_E_PRIME_CAMARGO1993,
						ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL,
						ENABLE_NHC_E_Q,
						ENABLE_TIMINGS,
						ENABLE_ITERATIVE_DISTANCE_COMPUTATION,
						ENABLE_MULTITHREADED_ROW_GENERATION,
						ENABLE_MULTITHREADED_MATRIX_GENERATION,
						W2V_PATH
					);
					if(err != 0){
						perror("failed to call apply_diversity_functions_to_graph_no_macros\n");
						return 1;
					}
					previous_best_s = best_s;
				} else { // end of comparison best_s?
					printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
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
					int32_t index = word2vec_key_to_index(&w2v, jdi.current_document.current_token);
					if(index != -1){
						if(w2v.keys[index].active_in_current_graph == 0){
							w2v.keys[index].active_in_current_graph = 1;
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
							local_node.word2vec_entry_pointer = &(w2v.keys[index]);
							local_node.vector.fp32 = w2v.keys[index].vector;
							local_node.num_dimensions = w2v.num_dimensions;
							local_node.already_considered = 0;
							local_node.relative_proportion = 1.0;
							local_node.absolute_proportion = 1;
							g.nodes[g.num_nodes] = local_node;
							w2v.keys[index].graph_node_pointer = &(g.nodes[g.num_nodes]);
							w2v.keys[index].graph_node_index = g.num_nodes;
							zipfian_n++;
							g.num_nodes++;
						} else {	
							g.nodes[w2v.keys[index].graph_node_index].absolute_proportion++;
						}
					}
				}
				if((UD_COLUMN != UD_MWE || found_at_least_one_mwe) && (num_documents % DOCUMENT_COUNT_RECOMPUTE_STEP == 0) && g.num_nodes > 1){
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
						// err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database);
						err = apply_diversity_functions_to_graph_no_macros(&g, &mst, &heap, f_ptr, &previous_g_num_nodes, &num_sentences, &num_all_sentences, &best_s, &mst_initialised, i, &sorted_array_discarded_because_not_in_vector_database,
							NUM_ROW_THREADS,
							NUM_MATRIX_THREADS,
							ENABLE_STIRLING,
							ENABLE_RICOTTA_SZEIDL,
							ENABLE_PAIRWISE,
							ENABLE_LEXICOGRAPHIC,
							ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY,
							ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY,
							ENABLE_LEINSTER_COBBOLD_DIVERSITY,
							ENABLE_FUNCTIONAL_EVENNESS,
							ENABLE_FUNCTIONAL_DISPERSION,
							ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED,
							ENABLE_NON_DISPARITY_FUNCTIONS,
							ENABLE_DISPARITY_FUNCTIONS,
							ENABLE_SHANNON_WEAVER_ENTROPY,
							ENABLE_GOOD_ENTROPY,
							ENABLE_RENYI_ENTROPY,
							ENABLE_PATIL_TAILLIE_ENTROPY,
							ENABLE_Q_LOGARITHMIC_ENTROPY,
							ENABLE_SIMPSON_INDEX,
							ENABLE_SIMPSON_DOMINANCE_INDEX,
							ENABLE_HILL_NUMBER_STANDARD,
							ENABLE_HILL_EVENNESS,
							ENABLE_BERGER_PARKER_INDEX,
							ENABLE_JUNGE1994_PAGE20,
							ENABLE_JUNGE1994_PAGE22,
							ENABLE_BRILLOUIN_DIVERSITY,
							ENABLE_MCINTOSH_INDEX,
							ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975,
							ENABLE_SW_E_HEIP,
							ENABLE_SW_E_ONE_MINUS_D,
							ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964,
							ENABLE_SW_E_MINUS_LN_D_PIELOU1977,
							ENABLE_SW_F_2_1_ALATALO1981,
							ENABLE_SW_G_2_1_MOLINARI1989,
							ENABLE_SW_E_BULLA1994,
							ENABLE_SW_O_BULLA1994,
							ENABLE_SW_E_MCI_PIELOU1969,
							ENABLE_SW_E_PRIME_CAMARGO1993,
							ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL,
							ENABLE_NHC_E_Q,
							ENABLE_TIMINGS,
							ENABLE_ITERATIVE_DISTANCE_COMPUTATION,
							ENABLE_MULTITHREADED_ROW_GENERATION,
							ENABLE_MULTITHREADED_MATRIX_GENERATION,
							W2V_PATH
						);
						if(err != 0){
							perror("failed to call apply_diversity_functions_to_graph_no_macros\n");
							return 1;
						}
						previous_best_s = best_s;
					} else { // end of comparison best_s?
						printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", best_s, previous_best_s, g.num_nodes, previous_g_num_nodes);
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


	for(int32_t i = 0 ; i < num_files ; i++){
		free(bfr[i]);
	}

	fclose(f_ptr);

	return 0;

	free_bfr_exit_failure:
	for(int32_t i = 0 ; i < num_files ; i++){
		free(bfr[i]);
	}
	return 1;
}

/*
int32_t main(void){
	int32_t err;

	struct word2vec w2v;
	err = load_word2vec_binary(&w2v, W2V_PATH);
	if(err != 0){
		perror("failed to call load_word2vec binary\n");
		return 1;
	}

	err = measurement(&w2v, TARGET_COLUMN);
	if(err != 0){
		perror("failed to call measurement\n");
		return 1;
	}

	return 0;
}
*/

