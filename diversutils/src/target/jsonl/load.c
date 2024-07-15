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

#include "cupt/constants.h"
#include "distributions.h"
#include "filter.h"
#include "graph.h"
#include "jsonl/parser.h"
#include "logging.h"
#include "measurement.h"
#include "sorted_array/array.h"
#include "unicode/utf8.h"

int32_t jsonl_to_graph(const uint64_t i, const char *const filename, struct measurement_configuration *const mcfg,
                       struct measurement_structure_references *const sref, struct measurement_mutables *const mmut) {
  struct jsonl_document_iterator jdi = {0};
  if (create_jsonl_document_iterator(&jdi, filename, mcfg->jsonl_content_key) != 0) {
    perror("failed to call create_jsonl_document_iterator\n");
    return 1;
  }

  const int32_t log_bfr_size = 256;
  char log_bfr[log_bfr_size];

  int8_t found_at_least_one_mwe = 0;
  while (!(jdi.file_is_done)) {
    memset(jdi.current_document.identifier, '\0', jdi.current_document.identifier_size); // ?
    jdi.current_document.identifier_size = 0;
    memset(jdi.current_document.text, '\0', jdi.current_document.text_size); // ?
    jdi.current_document.text_size = 0;
    if (iterate_jsonl_document_iterator(&jdi) != 0) {
      perror("failed to call iterate_jsonl_document_iterator\n");
      return 1;
    }
    if (jdi.current_document.text_size == 0 || jdi.current_document.identifier_size == 0) {
      continue;
    }

#if ENABLE_FILTER_ON_JSONL_DOCUMENTS == 1
    if (filter_substitute_all(&jdi.current_document.text, &jdi.current_document.text_size)) {
      perror("failed to call filter_substitute_all\n");
      return 1;
    };
    jdi.current_document.text_size = strlen(jdi.current_document.text);
#endif

    while (!(jdi.current_document.reached_last_token)) {
      if (iterate_document_current_token(&(jdi.current_document)) != 0) {
        perror("failed to call iterate_document_current_token\n");
        return 1;
      }

      // UTF-8 normalisation

      if (mcfg->enable_token_utf8_normalisation) {
        char *utf8_normalised_key = utf8_normalise(jdi.current_document.current_token);
        if (utf8_normalised_key != jdi.current_document.current_token) {
          // printf("Normalising \"%s\" to \"%s\"\n", jdi.current_document.current_token, utf8_normalised_key);
          size_t bytes_to_cpy = strlen(utf8_normalised_key);
          if (bytes_to_cpy > JSONL_CURRENT_TOKEN_BUFFER_SIZE - 1) {
            bytes_to_cpy = JSONL_CURRENT_TOKEN_BUFFER_SIZE - 1;
          }
          memcpy(jdi.current_document.current_token, utf8_normalised_key, bytes_to_cpy);
          jdi.current_document.current_token[bytes_to_cpy] = '\0';
        }
      }

      // add to graph
      int32_t index = word2vec_key_to_index(sref->w2v, jdi.current_document.current_token);
      if (index != -1) {
        pthread_mutex_lock(&(sref->w2v->keys[index].mutex));
        if (sref->w2v->keys[index].active_in_current_graph == 0) {
          sref->w2v->keys[index].active_in_current_graph = 1;

          pthread_mutex_lock(&(sref->g->mutex_nodes));
          if (sref->g->num_nodes == sref->g->capacity) {
            if (request_more_capacity_graph(sref->g) != 0) {
              perror("failed to call request_more_capacity_graph\n");
              return 1;
            }
          }
          struct graph_node local_node;
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
        pthread_mutex_unlock(&(sref->w2v->keys[index].mutex));
      } else { // added from cupt
        pthread_mutex_lock(&(sref->sorted_array_discarded_because_not_in_vector_database->mutex));
        int32_t index_in_discarded =
            key_to_index_sorted_array(sref->sorted_array_discarded_because_not_in_vector_database,
                                      jdi.current_document.current_token); // before, was in upper level
        if (index_in_discarded == -1) {
          struct sorted_array_str_int_element elem;
          if (create_sorted_array_str_int_element(&elem) != 0) {
            perror("failed to call create_sorted_array_str_int_element\n");
            return 1;
          }
          size_t bytes_to_cpy = strlen(jdi.current_document.current_token);
          if (bytes_to_cpy > SORTED_ARRAY_DEFAULT_KEY_SIZE - 1) {
            bytes_to_cpy = SORTED_ARRAY_DEFAULT_KEY_SIZE - 1;
          }
          memcpy(&elem.key, jdi.current_document.current_token, bytes_to_cpy);
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

    pthread_mutex_lock(&mmut->mutex);
    pthread_mutex_lock(&sref->g->mutex_nodes);
    // if((mcfg->target_column != UD_MWE || found_at_least_one_mwe) && (mcfg->steps.document.enable_count_recompute_step &&
    // (((!mcfg->steps.document.use_log10) && mmut->document.num % mcfg->steps.document.recompute_step == 0) ||
    // (mcfg->steps.document.use_log10 && mmut->document.num >= mmut->document.count_target)) && sref->g->num_nodes > 1)){ // DO
    // NOT REMOVE
    if ((mcfg->target_column != UD_MWE || found_at_least_one_mwe) &&
        (mcfg->steps.document.enable_count_recompute_step &&
         (((!mcfg->steps.document.use_log10) && mmut->document.num_all % mcfg->steps.document.recompute_step == 0) ||
          (mcfg->steps.document.use_log10 && mmut->document.num_all >= mmut->document.count_target)) &&
         sref->g->num_nodes > 1)) {
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

      if (mcfg->steps.document.use_log10) {
        mmut->document.stacked_log += mcfg->steps.document.recompute_step_log10;
        mmut->document.count_target = (uint64_t)floor(pow(10.0, mmut->document.stacked_log));
        memset(log_bfr, '\0', log_bfr_size);
        snprintf(log_bfr, log_bfr_size, "New document count target: %lu (10.0^%.3f)", mmut->document.count_target,
                 mmut->document.stacked_log);
        info_format(__FILE__, __func__, __LINE__, log_bfr);
      }
    }
    pthread_mutex_unlock(&sref->g->mutex_nodes);

    // mmut->document.num++; // DO NOT REMOVE
    mmut->document.num_all++;
    pthread_mutex_unlock(&mmut->mutex);
  }

  free_jsonl_document_iterator(&jdi);

  return 0;

panic_exit:

  free_jsonl_document_iterator(&jdi);

  return 1;
}

void *jsonl_to_graph_thread(void *args) {
  if (jsonl_to_graph(((struct measurement_file_thread *)args)->i, ((struct measurement_file_thread *)args)->filename,
                     ((struct measurement_file_thread *)args)->mcfg, ((struct measurement_file_thread *)args)->sref,
                     ((struct measurement_file_thread *)args)->mmut) != 0) {
    perror("Failed to call jsonl_to_graph in jsonl_to_graph_thread\n");
    exit(1);
  }
  return NULL;
}
