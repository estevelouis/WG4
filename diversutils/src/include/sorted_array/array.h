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

#ifndef SORTED_ARRAY_ARRAY_H
#define SORTED_ARRAY_ARRAY_H

#include <pthread.h>
#include <stdint.h>

#include "sorted_array/constants.h"

struct sorted_array_default_element {
  char key[SORTED_ARRAY_DEFAULT_KEY_SIZE];
  char value[SORTED_ARRAY_DEFAULT_VALUE_SIZE];
};

struct sorted_array_str_int_element {
  char key[SORTED_ARRAY_DEFAULT_KEY_SIZE];
  int64_t value;
};

struct sorted_array_int_int_element {
  int64_t key;
  int64_t value;
};

struct sorted_array_int_str_element {
  int64_t key;
  char value[SORTED_ARRAY_DEFAULT_VALUE_SIZE];
};

struct sorted_array {
  void *bfr;
  size_t element_size;
  int64_t num_elements;
  int64_t num_elements_sorted;
  int64_t capacity;
  int32_t (*comp)(const void *a, const void *b);
  int32_t to_free;
  pthread_mutex_t mutex;
};

int32_t create_sorted_array_default_element(struct sorted_array_default_element *elem);
int32_t create_sorted_array_str_int_element(struct sorted_array_str_int_element *elem);
int32_t create_sorted_array_int_int_element(struct sorted_array_int_int_element *elem);
int32_t create_sorted_array_int_str_element(struct sorted_array_int_str_element *elem);
int32_t create_sorted_array(struct sorted_array *const array, int64_t num_elements, size_t element_size,
                            int32_t (*comp)(const void *a, const void *b));
int32_t recreate_sorted_array(struct sorted_array *array, int64_t num_elements, size_t element_size,
                              int32_t (*comp)(const void *a, const void *b));
void free_sorted_array(struct sorted_array *array);
int32_t key_to_index_sorted_array(const struct sorted_array *const array, const void *const key);
void *key_to_elem_sorted_array(const struct sorted_array *const array, const void *const key);
int32_t bulk_insert_sorted_array(struct sorted_array *array, void *alt_array, size_t n);
int32_t insert_sorted_array(struct sorted_array *array, void *alt, int8_t insert_mode);
int32_t sorted_array_default_cmp(const void *a, const void *b);
int32_t sorted_array_str_int_cmp(const void *a, const void *b);
int32_t sorted_array_int_int_cmp(const void *a, const void *b);
int32_t sorted_array_int_str_cmp(const void *a, const void *b);

#endif
