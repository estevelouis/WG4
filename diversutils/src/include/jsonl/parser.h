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

#ifndef JSONL_PARSER_H
#define JSONL_PARSER_H

#include <stdint.h>
#include <stdio.h>

#include "jsonl/constants.h"
#if TOKENIZATION_METHOD == 0
#include <regex.h>
#endif

struct document {
  size_t identifier_size;
  size_t text_size;
  size_t identifier_capacity;
  size_t text_capacity;
#if TOKENIZATION_METHOD == 0
  regex_t reg;
  regoff_t latest_rm_eo;
#elif TOKENIZATION_METHOD == 1
  FILE *tmp_udpipe_output_file;
#elif TOKENIZATION_METHOD == 2
  FILE *tmp_udpipe_output_file;
  // int32_t output_pipefd[2];
  char *heap_char_output;
#endif
  int8_t reached_last_token;
  int8_t usable;
  char current_token[JSONL_CURRENT_TOKEN_BUFFER_SIZE];
  char *identifier;
  char *text;
};

struct jsonl_document_iterator {
  FILE *file_ptr;
  int8_t file_is_open;
  int8_t file_is_done;
  char bfr_read[JSONL_FILE_READ_BUFFER_SIZE];
  char content_key[JSONL_CONTENT_KEY_BUFFER_SIZE];
  struct document current_document;
  int32_t document_to_free;
  char *current_line;
  size_t current_line_cardinality;
  size_t current_line_capacity;
  // int64_t offset_start; // included
  // int64_t offset_end; // excluded
};

void jsonl_init_tokenization(void);

#if (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
int32_t launch_udpipe(struct document *const doc);
#endif

int32_t create_document(struct document *doc);
int32_t iterate_document_current_token(struct document *const doc);
void free_document(struct document *doc);
int32_t create_jsonl_document_iterator(struct jsonl_document_iterator *jdi, const char *file_name,
                                       const char *const content_key);
void free_jsonl_document_iterator(struct jsonl_document_iterator *jdi);
int32_t iterate_jsonl_document_iterator(struct jsonl_document_iterator *restrict const jdi);
int32_t jsonl_document_iterator_request_more_capacity_current_line(struct jsonl_document_iterator *const jdi);

#endif
