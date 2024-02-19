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

#include<math.h>
#include<assert.h>

#include "sorted_array/constants.h"

struct sorted_array_default_element {
	char key[SORTED_ARRAY_DEFAULT_KEY_SIZE];
	char value[SORTED_ARRAY_DEFAULT_VALUE_SIZE];
};

int32_t create_sorted_array_default_element(struct sorted_array_default_element* elem){
	memset(elem->key, '\0', SORTED_ARRAY_DEFAULT_KEY_SIZE);
	memset(elem->value, '\0', SORTED_ARRAY_DEFAULT_VALUE_SIZE);
	return 0;
}

struct sorted_array_str_int_element {
	char key[SORTED_ARRAY_DEFAULT_KEY_SIZE];
	int64_t value;
};

int32_t create_sorted_array_str_int_element(struct sorted_array_str_int_element* elem){
	memset(elem->key, '\0', SORTED_ARRAY_DEFAULT_KEY_SIZE);
	elem->value = 0;
	return 0;
}

struct sorted_array_int_int_element {
	int64_t key;
	int64_t value;
};

int32_t create_sorted_array_int_int_element(struct sorted_array_int_int_element* elem){
	elem->key = 0;
	elem->value = 0;
	return 0;
}

struct sorted_array_int_str_element {
	int64_t key;
	char value[SORTED_ARRAY_DEFAULT_VALUE_SIZE];
};

int32_t create_sorted_array_int_str_element(struct sorted_array_int_str_element* elem){
	elem->key = 0;
	memset(elem->value, '\0', SORTED_ARRAY_DEFAULT_VALUE_SIZE);
	return 0;
}


struct sorted_array {
	void* bfr;
	void* bfr_cpy;
	size_t element_size;
	int64_t num_elements;
	int64_t capacity;
	int32_t (*comp)(const void* a, const void* b);
	int32_t to_free;
};

int32_t create_sorted_array(struct sorted_array* array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b)){
	size_t malloc_size;
	if(num_elements > 0){
		malloc_size = num_elements * element_size;
		array->capacity = num_elements;
		array->num_elements = num_elements;
	} else {
		malloc_size = SORTED_ARRAY_STEP * element_size;
		array->capacity = SORTED_ARRAY_STEP;
		array->num_elements = 0;
	}
	void* malloc_pointer = malloc(malloc_size);
	if(malloc_pointer == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(malloc_pointer, '\0', malloc_size);
	array->bfr = malloc_pointer;
	array->bfr_cpy = array->bfr;
	array->to_free = 1;
	array->element_size = element_size;
	array->comp = comp;
	if(array->bfr == NULL){
		printf("array->bfr == NULL\n");
		exit(1);
	}
	return 0;
}

int32_t recreate_sorted_array(struct sorted_array* array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b)){
	size_t malloc_size;
	if(num_elements > 0){
		malloc_size = num_elements * element_size;
		array->capacity = num_elements;
		array->num_elements = num_elements;
	} else {
		malloc_size = SORTED_ARRAY_STEP * element_size;
		array->capacity = SORTED_ARRAY_STEP;
		array->num_elements = 0;
	}
	assert(array->bfr == array->bfr_cpy);
	if(array->bfr != array->bfr_cpy){
		perror("!\n");
		exit(1);
	}
	if(malloc_size == 0){
		printf("new_size == 0\n");
		exit(1);
	}
	array->bfr = realloc(array->bfr, malloc_size);
	if(array->bfr == NULL){
		perror("failed to realloc\n");
		/* perror("failed to malloc\n"); */
		return 1;
	}
	array->bfr_cpy = array->bfr;
	memset(array->bfr, '\0', malloc_size);
	array->to_free = 1;
	array->capacity = num_elements;
	array->num_elements = num_elements;
	array->element_size = element_size;
	array->comp = comp;
	if(array->bfr == NULL){
		printf("array->bfr == NULL\n");
		exit(1);
	}
	return 0;
}

void free_sorted_array(struct sorted_array* array){
	assert(array->bfr == array->bfr_cpy);
	if(array->bfr != array->bfr_cpy){
		perror("!\n");
		exit(1);
	}
	if(array->to_free == 1){
		free(array->bfr);
		array->to_free = 0;
	}
}

int32_t key_to_index_sorted_array(struct sorted_array* array, void* key){
	if(array->num_elements == 0){
		return -1;
	}
	int32_t lower_bound = 0;
	int32_t upper_bound = array->num_elements;

	while(lower_bound <= upper_bound){
		int32_t middle_index = lower_bound + floor((upper_bound - lower_bound) / 2);
		int32_t cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (middle_index * array->element_size)), key); // compliance
		if(cmp_result == 0){
			return middle_index;
		} else if(cmp_result < 0){
			lower_bound = middle_index + 1;
		} else if(cmp_result > 0){
			upper_bound = middle_index - 1;
		}
	}

	return -1;
}

void* key_to_elem_sorted_array(struct sorted_array* array, void* key){
	int32_t index = key_to_index_sorted_array(array, key);
	return ((void*) (((uint64_t) array->bfr) + (index * array->element_size))); // compliance
}

int32_t bulk_insert_sorted_array(struct sorted_array* array, void* alt_array, size_t n){
	size_t available_space = (array->capacity - array->num_elements) * array->element_size;
	size_t size_delta = n - available_space;
	if(size_delta > 0){
		size_t new_size = array->capacity * array->element_size + (size_delta + (array->element_size - (size_delta % array->element_size)) % array->element_size);
	    	if(new_size == 0){
	    		printf("new_size == 0\n");
	    		exit(1);
	    	}
		array->bfr = realloc(array->bfr, new_size);
		if(array->bfr == NULL){
			perror("failed to realloc\n");
			return 1;
		}
		memset((void*) (((uint64_t) array->bfr) + (array->capacity * array->element_size)), '\0', new_size - (array->capacity * array->element_size)); // compliance
		array->capacity = new_size / array->element_size;
	}
	memcpy((void*) (((uint64_t) array->bfr) + (array->num_elements * array->element_size)), alt_array, n); // compliance
	array->num_elements += n / array->element_size;

	qsort(array->bfr, array->num_elements, n, array->comp);

	return 0;
}

int32_t insert_sorted_array(struct sorted_array* array, void* alt, int8_t insert_mode){
	/*
	 * insert_mode == 0 -> insert whatever happens
	 * insert_mode == 1 -> insert only if not present
	 * insert_mode == 2 -> overwrite if present, otherwise insert
	 * */
	if(insert_mode > 0){
		int32_t index = key_to_index_sorted_array(array, alt);
		if(insert_mode == 1){
			if(index != -1){
				return 0;
			}
		} else if(insert_mode == 2){
			if(index != -1){
				memcpy((void*) (((uint64_t) array->bfr) + (index * array->element_size)), alt, array->element_size); // compliance
				return 0;
			}
		}
	}
	if(array->num_elements == array->capacity){
		size_t new_size = (array->capacity + SORTED_ARRAY_STEP) * array->element_size;
		array->bfr = realloc(array->bfr, new_size);
		if(array->bfr == NULL){
			perror("failed to realloc\n");
			return 1;
		}
		memset((void*) (((uint64_t) array->bfr) + (array->capacity * array->element_size)), '\0', SORTED_ARRAY_STEP * array->element_size); // compliance
		array->capacity += SORTED_ARRAY_STEP;
	}
	memcpy((void*) (((uint64_t) array->bfr) + (array->num_elements * array->element_size)), alt, array->element_size); // compliance
	array->num_elements++;

	qsort(array->bfr, array->num_elements, array->element_size, array->comp);

	return 0;
}

int32_t sorted_array_default_cmp(const void* a, const void* b){
	size_t size_a = sizeof(a);
	size_t size_b = sizeof(b);
	size_t min_size = size_a;
	if(size_b < size_a){
		min_size = size_b;
	}
	for(size_t i = 0 ; i < min_size ; i++){
		uint8_t a_;
		uint8_t b_;
		memcpy(&(a_), (void*) (((uint64_t) a) + i), 1); // compliance
		memcpy(&(b_), (void*) (((uint64_t) b) + i), 1); // compliance
		if(a_ != b_){
			return (int32_t) (a_ - b_);
		}
	}
	return 0;
}

int32_t sorted_array_str_int_cmp(const void* a, const void* b){
	return sorted_array_default_cmp(a, b);
}

int32_t sorted_array_int_int_cmp(const void* a, const void* b){
	return (int32_t) (((struct sorted_array_int_int_element*) a)->key - ((struct sorted_array_int_int_element*) b)->key);
}

int32_t sorted_array_int_str_cmp(const void* a, const void* b){
	return sorted_array_int_int_cmp(a, b);
}

#endif
