#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "sorted_array/array.h"
#include "cupt/parser.h"
#include "distributions.h"
#include "measurement.h"
#include "logging.h"

int32_t cupt_to_graph(const uint64_t i, const char * const filename, struct measurement_configuration * const mcfg, struct measurement_structure_references * const sref, struct measurement_mutables * const mmut){
    struct cupt_sentence_iterator csi = {0};
    if(create_cupt_sentence_iterator(&csi, filename) != 0){
        perror("failed to call create_cupt_sentence_iterator\n");
        return 1;
    }

    const int32_t log_bfr_size = 256;
    char log_bfr[log_bfr_size];

    if(iterate_cupt_sentence_iterator(&csi) != 0){
        perror("failed to call iterate_cupt_sentence_iterator\n");
        return 1;
    }

    int8_t found_at_least_one_mwe = 0;
    while(!(csi.file_is_done)){
		for(int32_t j = 0 ; j < csi.current_sentence.num_tokens ; j++){
    		int32_t index;
    		int32_t index_in_discarded = -1;
    		const int32_t key_size = 256;
    		char key[key_size];
    		memset(key, '\0', key_size);
    		size_t len = 0;
    		switch(mcfg->target_column){
    			case UD_FORM:
    				len = strlen(csi.current_sentence.tokens[j].form);
    				if(len > key_size - 1){
    					len = key_size - 1;
    				}
    				memcpy(key, csi.current_sentence.tokens[j].form, len);
    				break;
    			case UD_LEMMA:
    				len = strlen(csi.current_sentence.tokens[j].lemma);
    				if(len > key_size - 1){
    					len = key_size - 1;
    				}
    				memcpy(key, csi.current_sentence.tokens[j].lemma, len);
    				break;
    			default:
    				index = -1;
    				perror("target_column not properly defined\n");
    				return 1;
    		}
    
					index = word2vec_key_to_index(sref->w2v, key);
	
					if(index != -1){
                        pthread_mutex_lock(&(sref->w2v->keys[index].mutex));
						if(sref->w2v->keys[index].active_in_current_graph == 0){
							sref->w2v->keys[index].active_in_current_graph = 1;
							// num_nodes++;
	
                            pthread_mutex_lock(&(sref->g->mutex_nodes));
							if(sref->g->num_nodes == sref->g->capacity){
								if(request_more_capacity_graph(sref->g) != 0){
									perror("failed to call request_more_capacity_graph\n");
									return 1;
								}
							}
							struct graph_node local_node;
							if(create_graph_node(&local_node, sref->g->nodes[0].num_dimensions, FP32) != 0){
								perror("failed to call create_graph_node\n");
								return 1;
							}
							local_node.word2vec_entry_pointer = &(sref->w2v->keys[index]);
							local_node.vector.fp32 = sref->w2v->keys[index].vector;
	
							local_node.num_dimensions = sref->w2v->num_dimensions;
							local_node.already_considered = 0;
							local_node.relative_proportion = 1.0;
							local_node.absolute_proportion = 1;
							sref->g->nodes[sref->g->num_nodes] = local_node;
							sref->w2v->keys[index].graph_node_pointer = &(sref->g->nodes[sref->g->num_nodes]);
							sref->w2v->keys[index].graph_node_index = sref->g->num_nodes;
							sref->g->num_nodes++;
                            pthread_mutex_unlock(&(sref->g->mutex_nodes));
						} else {
                            pthread_mutex_lock(&(sref->g->nodes[sref->w2v->keys[index].graph_node_index].mutex_local_node));
							sref->g->nodes[sref->w2v->keys[index].graph_node_index].absolute_proportion++;
                            pthread_mutex_unlock(&(sref->g->nodes[sref->w2v->keys[index].graph_node_index].mutex_local_node));
						}
						sref->w2v->keys[index].num_occurrences++; // ? mutex ?
                        pthread_mutex_unlock(&(sref->w2v->keys[index].mutex));
					} else {
                        pthread_mutex_lock(&(sref->sorted_array_discarded_because_not_in_vector_database->mutex));
					    index_in_discarded = key_to_index_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database, key); // before, was in upper level
						if(index_in_discarded == -1){
							struct sorted_array_str_int_element elem;
							size_t bytes_to_cpy = strlen(key);
							if(bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1){
								bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
							}
							memcpy(&elem.key, key, bytes_to_cpy);
							elem.value = 1;
							if(insert_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database, &elem, 0) != 0){
								perror("failed to call insert_sorted_array\n");
								return 1;
							}
						} else {
							((struct sorted_array_str_int_element*) sref->sorted_array_discarded_because_not_in_vector_database->bfr)[index_in_discarded].value++;
						}
                        pthread_mutex_unlock(&(sref->sorted_array_discarded_because_not_in_vector_database->mutex));
					}
        }

                pthread_mutex_lock(&mmut->mutex);
                pthread_mutex_lock(&sref->g->mutex_nodes);
				// sentence level recomputation
				if((mcfg->target_column != UD_MWE || found_at_least_one_mwe) && (mcfg->steps.sentence.enable_count_recompute_step && (((!mcfg->steps.sentence.use_log10) && mmut->sentence.num % mcfg->steps.sentence.recompute_step == 0) || (mcfg->steps.sentence.use_log10 && mmut->sentence.num >= mmut->sentence.count_target))) && sref->g->num_nodes > 1){
					// printf("found_at_least_one_mwe: %i; g->num_nodes: %lu\n", found_at_least_one_mwe, sref->g->num_nodes);
                    memset(log_bfr, '\0', log_bfr_size);
					snprintf(log_bfr, log_bfr_size, "found_at_least_one_mwe: %i; g->num_nodes: %lu", found_at_least_one_mwe, sref->g->num_nodes);
                    info_format(__FILE__, __func__, __LINE__, log_bfr);

                    compute_graph_relative_proportions(sref->g);
				
					int32_t err = zipfian_fit_from_graph(sref->g, &mmut->best_s);
					if(err != 0){
						perror("failed to call zipfian_fit_from_graph\n");
						goto panic_exit;
					}
				
					if(mmut->best_s != mmut->prev_best_s || sref->g->num_nodes != ((uint64_t) mmut->prev_num_nodes)){
						memset(log_bfr, '\0', log_bfr_size);
						snprintf(log_bfr, log_bfr_size, "best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li", mmut->best_s, sref->g->num_nodes, mmut->sentence.num, mmut->document.num);
						info_format(__FILE__, __func__, __LINE__, log_bfr);

                        err = apply_diversity_functions_to_graph(i, mcfg, sref, mmut);
						if(err != 0){
							perror("failed to call apply_diversity_functions_to_graph\n");
							return 1;
						}
						mmut->prev_best_s = mmut->best_s;
					} else { // end of comparison best_s?
						printf("ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n", mmut->best_s, mmut->prev_best_s, sref->g->num_nodes, mmut->prev_num_nodes);
					}

					if(mcfg->steps.sentence.use_log10){
						mmut->sentence.stacked_log += mcfg->steps.sentence.recompute_step_log10;
						mmut->sentence.count_target = (uint64_t) floor(pow(10.0, mmut->sentence.stacked_log));
                        memset(log_bfr, '\0', log_bfr_size);
						snprintf(log_bfr, log_bfr_size, "New sentence count target: %lu (10.0^%.3f)", mmut->sentence.count_target, mmut->sentence.stacked_log);
                        info_format(__FILE__, __func__, __LINE__, log_bfr);
					}
                }
                pthread_mutex_unlock(&sref->g->mutex_nodes);
        
                mmut->sentence.num++;

                pthread_mutex_unlock(&mmut->mutex);


        if(iterate_cupt_sentence_iterator(&csi) != 0){
            perror("failed to call iterate_cupt_sentence_iterator\n");
            return 1;
        } 
    }

    free_cupt_sentence_iterator(&csi);

    return 0;

    panic_exit:

    free_cupt_sentence_iterator(&csi);

    return 1;
}

void * cupt_to_graph_thread(void * args){
    if(cupt_to_graph(
        ((struct measurement_file_thread *) args)->i,
        ((struct measurement_file_thread *) args)->filename,
        ((struct measurement_file_thread *) args)->mcfg,
        ((struct measurement_file_thread *) args)->sref,
        ((struct measurement_file_thread *) args)->mmut
    ) != 0){
        perror("Failed to call cupt_to_graph in cupt_to_graph_thread\n");
        exit(1);
    }
    return NULL;
}

