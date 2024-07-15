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
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cupt/parser.h"
#include "distributions.h"
#include "graph.h"
#include "logging.h"
#include "measurement.h"
#include "sorted_array/array.h"
#include "unicode/utf8.h"

int32_t cupt_to_graph(const uint64_t i, const char *const filename, const char *const filename_tp,
                      struct measurement_configuration *const mcfg, struct measurement_structure_references *const sref,
                      struct measurement_mutables *const mmut) {
  const int32_t log_bfr_size = 256;
  char log_bfr[log_bfr_size];

  struct cupt_sentence_iterator csi = {0};
  struct cupt_sentence_iterator csi_tp = {0};
  if (create_cupt_sentence_iterator(&csi, filename) != 0) {
    goto failure_create_cupt_sentence_iterator;
  }
  if (filename_tp != NULL) {
    if (create_cupt_sentence_iterator(&csi_tp, filename_tp) != 0) {
      goto failure_create_cupt_sentence_iterator;
    }
  }

  if (iterate_cupt_sentence_iterator(&csi) != 0) {
    goto failure_iterate_cupt_sentence_iterator;
  }
  if (filename_tp != NULL) {
    if (iterate_cupt_sentence_iterator(&csi_tp) != 0) {
      goto failure_iterate_cupt_sentence_iterator;
    }
  }

  int8_t found_at_least_one_mwe = 0;
  while (!(csi.file_is_done)) {
    const size_t max_mwe = 32;
    const size_t max_tokens_per_mwe = 32;
    const size_t size_token_mwe = 32;
    const size_t size_mwe_class = 16;
    size_t mwe_bfr_size = max_mwe * max_tokens_per_mwe * size_token_mwe;
    size_t mwe_bfr_size_mwe_class = max_mwe * size_mwe_class;

    char mwe[mwe_bfr_size];
    memset(mwe, '\0', mwe_bfr_size);

    char mwe_tp[mwe_bfr_size];
    memset(mwe_tp, '\0', mwe_bfr_size);

    char mwe_class[mwe_bfr_size_mwe_class];
    memset(mwe_class, '\0', mwe_bfr_size_mwe_class);

    /*
    char mwe_class_tp[mwe_bfr_size_mwe_class];
    memset(mwe_class_tp, '\0', mwe_bfr_size_mwe_class);
    */

    int32_t mwe_lengths[max_mwe];
    memset(mwe_lengths, '\0', max_mwe * sizeof(int32_t));

    int32_t mwe_tp_lengths[max_mwe];
    memset(mwe_tp_lengths, '\0', max_mwe * sizeof(int32_t));

    int8_t mwe_correct_span[max_mwe];
    memset(mwe_correct_span, 1, max_mwe * sizeof(int8_t));

    for (int32_t j = 0; j < csi.current_sentence.num_tokens; j++) {
      if (mcfg->target_column == UD_MWE) {
        if (strcmp(csi.current_sentence.tokens[j].mwe, "") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "_") == 0 ||
            strcmp(csi.current_sentence.tokens[j].mwe, "-") == 0 || strcmp(csi.current_sentence.tokens[j].mwe, "*") == 0) {
          continue;
        }

        char *strtok_placeholder = NULL;
        strtok_placeholder = strtok(csi.current_sentence.tokens[j].mwe, ";");
        while (strtok_placeholder != NULL) {
          size_t mwe_num = (size_t)strtol(strtok_placeholder, NULL, 10);
          if (mwe_num < max_mwe) {
            size_t bytes_to_cpy = strlen(csi.current_sentence.tokens[j].lemma);
            if (bytes_to_cpy > size_token_mwe - 1) {
              bytes_to_cpy = size_token_mwe - 1;
            }
            memcpy(&(mwe[mwe_num * max_tokens_per_mwe * size_token_mwe + mwe_lengths[mwe_num] * size_token_mwe]),
                   csi.current_sentence.tokens[j].lemma, bytes_to_cpy);
            mwe_lengths[mwe_num]++;

            char *ptr_colon = strchr(strtok_placeholder, ':');
            if (ptr_colon != NULL) {
              bytes_to_cpy = strlen(ptr_colon + 1);
              if (bytes_to_cpy > size_mwe_class - 1) {
                bytes_to_cpy = size_mwe_class - 1;
              }
              memcpy(&(mwe_class[mwe_num * size_mwe_class]), ptr_colon + 1, bytes_to_cpy);
            }
          }

          /**/ // DO NOT REMOVE
          if (filename_tp != NULL && (strcmp(csi_tp.current_sentence.tokens[j].mwe, "") == 0 ||
                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "_") == 0 ||
                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "-") == 0 ||
                                      strcmp(csi_tp.current_sentence.tokens[j].mwe, "*") == 0)) {
            mwe_correct_span[mwe_num] = 0;
          }
          /**/

          strtok_placeholder = strtok(NULL, ";");
        }

        continue;
      }

      int32_t index;
      int32_t index_in_discarded = -1;
      const size_t key_size = 256;
      char key[key_size];
      memset(key, '\0', key_size);
      size_t len = 0;
      switch (mcfg->target_column) {
      case UD_FORM:
        len = strlen(csi.current_sentence.tokens[j].form);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].form, len);
        break;
      case UD_LEMMA:
        len = strlen(csi.current_sentence.tokens[j].lemma);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].lemma, len);
        break;
      case UD_UPOS:
        len = strlen(csi.current_sentence.tokens[j].upos);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].upos, len);
        break;
      case UD_XPOS:
        len = strlen(csi.current_sentence.tokens[j].xpos);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].xpos, len);
        break;
      case UD_FEATS:
        len = strlen(csi.current_sentence.tokens[j].feats);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].feats, len);
        break;
      case UD_HEAD:
        len = strlen(csi.current_sentence.tokens[j].head);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].head, len);
        break;
      case UD_DEPREL:
        len = strlen(csi.current_sentence.tokens[j].deprel);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].deprel, len);
        break;
      case UD_DEPS:
        len = strlen(csi.current_sentence.tokens[j].deps);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].deps, len);
        break;
      case UD_MISC:
        len = strlen(csi.current_sentence.tokens[j].misc);
        if (len > key_size - 1) {
          len = key_size - 1;
        }
        memcpy(key, csi.current_sentence.tokens[j].misc, len);
        break;
      default:
        index = -1;
        perror("target_column not properly defined\n");
        return 1;
      }

      // UTF-8 normalisation

      if (mcfg->enable_token_utf8_normalisation) {
        char *utf8_normalised_key = utf8_normalise(key);
        if (utf8_normalised_key != key) {
          // printf("Normalising \"%s\" to \"%s\"\n", key, utf8_normalised_key);
          size_t bytes_to_cpy = strlen(utf8_normalised_key);
          if (bytes_to_cpy > key_size - 1) {
            bytes_to_cpy = key_size - 1;
          }
          memcpy(key, utf8_normalised_key, bytes_to_cpy);
          key[bytes_to_cpy] = '\0';
        }
      }

      // retrieval in word2vec

      index = word2vec_key_to_index(sref->w2v, key);

      if (index != -1) {
        pthread_mutex_lock(&(sref->w2v->keys[index].mutex));
        if (sref->w2v->keys[index].active_in_current_graph == 0) {
          sref->w2v->keys[index].active_in_current_graph = 1;
          // num_nodes++;

          pthread_mutex_lock(&(sref->g->mutex_nodes));
          if (sref->g->num_nodes == sref->g->capacity) {
            if (request_more_capacity_graph(sref->g) != 0) {
              perror("failed to call request_more_capacity_graph\n");
              return 1;
            }
          }
          struct graph_node local_node = {0};
          if (create_graph_node(&local_node, sref->g->nodes[0].num_dimensions, FP32) != 0) {
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
        index_in_discarded = key_to_index_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database,
                                                       key); // before, was in upper level
        if (index_in_discarded == -1) {
          struct sorted_array_str_int_element elem;
          if (create_sorted_array_str_int_element(&elem) != 0) {
            perror("failed to call create_sorted_array_str_int_element\n");
            return 1;
          }
          size_t bytes_to_cpy = strlen(key);
          if (bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1) {
            bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
          }
          memcpy(&elem.key, key, bytes_to_cpy);
          elem.value = 1;
          if (insert_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database, &elem, 0) != 0) {
            perror("failed to call insert_sorted_array\n");
            return 1;
          }
        } else {
          ((struct sorted_array_str_int_element *)
               sref->sorted_array_discarded_because_not_in_vector_database->bfr)[index_in_discarded]
              .value++;
        }
        pthread_mutex_unlock(&(sref->sorted_array_discarded_because_not_in_vector_database->mutex));
      }
    }

    // MWE

    if (mcfg->target_column == UD_MWE) {
      uint8_t in_this_specific_sentence_found_at_least_one_correct_mwe = 0;
      uint8_t in_this_specific_sentence_found_at_least_one_mwe = 0;
      for (size_t k = 0; k < max_mwe; k++) {
        if (mwe_lengths[k] == 0) {
          continue;
        }
        in_this_specific_sentence_found_at_least_one_mwe = 1;

        /**/ // DO NOT REMOVE
        if (filename_tp != NULL && !(mwe_correct_span[k])) {
          continue;
        }
        in_this_specific_sentence_found_at_least_one_correct_mwe = 1;
        /**/

        qsort(&(mwe[k * max_tokens_per_mwe * size_token_mwe]), mwe_lengths[k], size_token_mwe, void_strcmp);

        const int32_t bfr_size = 512;
        int32_t index_bfr = 0;
        char bfr[bfr_size];
        memset(bfr, '\0', bfr_size);

        // memcpy(bfr + index_bfr, "_MWE_", 5);
        // index_bfr += 5;
        memcpy(bfr + index_bfr, "_MWE-", 5);
        index_bfr += 5;

        const size_t len_mwe_class = strlen(&(mwe_class[k * size_mwe_class]));
        memcpy(bfr + index_bfr, &(mwe_class[k * size_mwe_class]), len_mwe_class);
        index_bfr += len_mwe_class;

        bfr[index_bfr] = '_';
        index_bfr++;

        size_t bytes_to_cpy = 0;
        for (int32_t m = 0; m < mwe_lengths[k]; m++) {
          bytes_to_cpy = strlen(&(mwe[k * max_tokens_per_mwe * size_token_mwe + m * size_token_mwe]));
          if (m < mwe_lengths[k] - 1) {
            bytes_to_cpy++;
          }
          if (bytes_to_cpy > (size_t)(bfr_size - 1 - index_bfr)) {
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

        // printf("created key: %s\n", bfr);

        // pthread_mutex_lock(&sref->w2v->mutex);
        int32_t index = word2vec_key_to_index(sref->w2v, bfr);
        if (index != -1) {
          found_at_least_one_mwe = 1;

          pthread_mutex_lock(&sref->w2v->keys[index].mutex);
          pthread_mutex_lock(&sref->g->mutex_nodes);

          if (sref->w2v->keys[index].active_in_current_graph == 0) {
            sref->w2v->keys[index].active_in_current_graph = 1;

            if (sref->g->num_nodes == sref->g->capacity) {
              int32_t err = request_more_capacity_graph(sref->g);
              if (err != 0) {
                perror("failed to call request_more_capacity_graph\n");
                return 1;
              }
            }
            struct graph_node local_graph_node = {0};
            if (create_graph_node(&local_graph_node, sref->g->nodes[0].num_dimensions, FP32) !=
                0) { // copied and modified from above
              perror("failed to call create_graph_node\n");
              return 1;
            }
            local_graph_node.num_dimensions = (int16_t)sref->w2v->num_dimensions;
            local_graph_node.already_considered = 0;
            local_graph_node.relative_proportion = 1.0;
            local_graph_node.absolute_proportion = 1;
            local_graph_node.vector.fp32 = sref->w2v->keys[index].vector;
            local_graph_node.word2vec_entry_pointer = &(sref->w2v->keys[index]);
            sref->g->nodes[sref->g->num_nodes] = local_graph_node;
            sref->w2v->keys[index].graph_node_pointer = &(sref->g->nodes[sref->g->num_nodes]);
            sref->w2v->keys[index].graph_node_index = sref->g->num_nodes;
            sref->g->num_nodes++;
          } else {
            sref->g->nodes[sref->w2v->keys[index].graph_node_index].absolute_proportion++;
          }

          pthread_mutex_unlock(&sref->g->mutex_nodes);
          pthread_mutex_unlock(&sref->w2v->keys[index].mutex);

          sref->w2v->keys[index].num_occurrences++;
        } else {
          pthread_mutex_lock(&sref->sorted_array_discarded_because_not_in_vector_database->mutex);
          int32_t index_in_discarded =
              key_to_index_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database, bfr);

          if (index_in_discarded == -1) {
            struct sorted_array_str_int_element elem;
            if (create_sorted_array_str_int_element(&elem) != 0) {
              perror("failed to call create_sorted_array_str_int_element\n");
              return 1;
            }
            size_t bytes_to_cpy = strlen(bfr);
            if (bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1) {
              bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
            }
            memcpy(&elem.key, bfr, bytes_to_cpy);
            elem.value = 1;
            if (insert_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database, &elem, 0) != 0) {
              perror("failed to call insert_sorted_array\n");
              return 1;
            }
          } else {
            ((struct sorted_array_str_int_element *)
                 sref->sorted_array_discarded_because_not_in_vector_database->bfr)[index_in_discarded]
                .value++;
          }
          pthread_mutex_unlock(&sref->sorted_array_discarded_because_not_in_vector_database->mutex);
        }
        // pthread_mutex_unlock(&sref->w2v->mutex);
      }

      /*
                  if(found_at_least_one_mwe){ // for efficient recomputation
          pthread_mutex_lock(&mmut->mutex);
                          // mmut->sentence.num_all++;
          pthread_mutex_unlock(&mmut->mutex);
                  }
      */
      if (in_this_specific_sentence_found_at_least_one_mwe) {
        pthread_mutex_lock(&mmut->mutex);
        mmut->sentence.num_containing_mwe++;
        if (in_this_specific_sentence_found_at_least_one_correct_mwe) {
          mmut->sentence.num_containing_mwe_tp_only++;
        }
        pthread_mutex_unlock(&mmut->mutex);
      }
    }

    pthread_mutex_lock(&mmut->mutex);
    pthread_mutex_lock(&sref->g->mutex_nodes);
    // sentence level recomputation
    // if((mcfg->target_column != UD_MWE || found_at_least_one_mwe) && (mcfg->steps.sentence.enable_count_recompute_step &&
    // (((!mcfg->steps.sentence.use_log10) && mmut->sentence.num % mcfg->steps.sentence.recompute_step == 0) ||
    // (mcfg->steps.sentence.use_log10 && mmut->sentence.num >= mmut->sentence.count_target))) && sref->g->num_nodes > 1){ // DO
    // NOT REMOVE
    if ((mcfg->target_column != UD_MWE || found_at_least_one_mwe) &&
        (mcfg->steps.sentence.enable_count_recompute_step &&
         (((!mcfg->steps.sentence.use_log10) && mmut->sentence.num_all % mcfg->steps.sentence.recompute_step == 0) ||
          (mcfg->steps.sentence.use_log10 && mmut->sentence.num_all >= mmut->sentence.count_target))) &&
        sref->g->num_nodes > 1) {
      // printf("found_at_least_one_mwe: %i; g->num_nodes: %lu\n", found_at_least_one_mwe, sref->g->num_nodes);
      memset(log_bfr, '\0', log_bfr_size);
      snprintf(log_bfr, log_bfr_size, "found_at_least_one_mwe: %i; g->num_nodes: %lu", found_at_least_one_mwe,
               sref->g->num_nodes);
      info_format(__FILE__, __func__, __LINE__, log_bfr);

      compute_graph_relative_proportions(sref->g);

      int32_t err = zipfian_fit_from_graph(sref->g, &mmut->best_s);
      if (err != 0) {
        perror("failed to call zipfian_fit_from_graph\n");
        goto panic_exit;
      }

      if (mmut->best_s != mmut->prev_best_s || sref->g->num_nodes != ((uint64_t)mmut->prev_num_nodes)) {
        memset(log_bfr, '\0', log_bfr_size);
        // snprintf(log_bfr, log_bfr_size, "best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li", mmut->best_s,
        // sref->g->num_nodes, mmut->sentence.num, mmut->document.num); // DO NOT REMOVE
        snprintf(log_bfr, log_bfr_size, "best_s: %f; num_nodes: %lu; num_sentences: %li; num_documents: %li", mmut->best_s,
                 sref->g->num_nodes, mmut->sentence.num_all, mmut->document.num_all);
        info_format(__FILE__, __func__, __LINE__, log_bfr);

        err = apply_diversity_functions_to_graph(i, mcfg, sref, mmut);
        if (err != 0) {
          perror("failed to call apply_diversity_functions_to_graph\n");
          return 1;
        }
        mmut->prev_best_s = mmut->best_s;
      } else { // end of comparison best_s?
        printf(
            "ignoring because best_s (%.12f) == previous_best_s (%.12f) && g->num_nodes (%lu) == previous_g_num_nodes (%li)\n",
            mmut->best_s, mmut->prev_best_s, sref->g->num_nodes, mmut->prev_num_nodes);
      }

      if (mcfg->steps.sentence.use_log10) {
        mmut->sentence.stacked_log += mcfg->steps.sentence.recompute_step_log10;
        mmut->sentence.count_target = (uint64_t)floor(pow(10.0, mmut->sentence.stacked_log));
        memset(log_bfr, '\0', log_bfr_size);
        snprintf(log_bfr, log_bfr_size, "New sentence count target: %lu (10.0^%.3f)", mmut->sentence.count_target,
                 mmut->sentence.stacked_log);
        info_format(__FILE__, __func__, __LINE__, log_bfr);
      }
    }
    pthread_mutex_unlock(&sref->g->mutex_nodes);

    mmut->sentence.num_all++;

    pthread_mutex_unlock(&mmut->mutex);

    if (iterate_cupt_sentence_iterator(&csi) != 0) {
      goto failure_iterate_cupt_sentence_iterator;
    }
    if (filename_tp != NULL) {
      if (iterate_cupt_sentence_iterator(&csi_tp) != 0) {
        goto failure_iterate_cupt_sentence_iterator;
      }
    }
  }

  free_cupt_sentence_iterator(&csi);
  if (filename_tp != NULL) {
    free_cupt_sentence_iterator(&csi_tp);
  }

  return 0;

failure_iterate_cupt_sentence_iterator:
  perror("failed to call iterate_cupt_sentence_iterator\n");
  goto panic_exit;

failure_create_cupt_sentence_iterator:
  perror("failed to call create_cupt_sentence_iterator\n");

panic_exit:

  free_cupt_sentence_iterator(&csi);
  if (filename_tp != NULL) {
    free_cupt_sentence_iterator(&csi_tp);
  }

  return 1;
}

void *cupt_to_graph_thread(void *args) {
  if (cupt_to_graph(((struct measurement_file_thread *)args)->i, ((struct measurement_file_thread *)args)->filename,
                    ((struct measurement_file_thread *)args)->filename_tp, ((struct measurement_file_thread *)args)->mcfg,
                    ((struct measurement_file_thread *)args)->sref, ((struct measurement_file_thread *)args)->mmut) != 0) {
    perror("Failed to call cupt_to_graph in cupt_to_graph_thread\n");
    exit(1);
  }
  return NULL;
}
