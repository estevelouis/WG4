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

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sorted_array/array.h"
#include "sorted_array/constants.h"

int32_t create_sorted_array_default_element(struct sorted_array_default_element* elem){
	memset(elem->key, '\0', SORTED_ARRAY_DEFAULT_KEY_SIZE);
	memset(elem->value, '\0', SORTED_ARRAY_DEFAULT_VALUE_SIZE);
	return 0;
}

int32_t create_sorted_array_str_int_element(struct sorted_array_str_int_element* elem){
	memset(elem->key, '\0', SORTED_ARRAY_DEFAULT_KEY_SIZE);
	elem->value = 0;
	return 0;
}

int32_t create_sorted_array_int_int_element(struct sorted_array_int_int_element* elem){
	elem->key = 0;
	elem->value = 0;
	return 0;
}

int32_t create_sorted_array_int_str_element(struct sorted_array_int_str_element* elem){
	elem->key = 0;
	memset(elem->value, '\0', SORTED_ARRAY_DEFAULT_VALUE_SIZE);
	return 0;
}

int32_t create_sorted_array(struct sorted_array * const array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b), const uint8_t method){
	size_t malloc_size;
	if(num_elements > 0){
		malloc_size = num_elements * element_size;
		array->capacity = num_elements;
	} else {
		malloc_size = SORTED_ARRAY_STEP * element_size;
		array->capacity = SORTED_ARRAY_STEP;
	}
	array->num_elements = 0;
	void* malloc_pointer = malloc(malloc_size);
	if(malloc_pointer == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(malloc_pointer, '\0', malloc_size);
	array->bfr = malloc_pointer;
	array->to_free = 1;
	array->element_size = element_size;
	array->comp = comp;
    array->num_elements_sorted = 0;
	if(array->bfr == NULL){
		printf("array->bfr == NULL\n");
		exit(1);
	}
    array->method = method;

    if(pthread_mutex_init(&(array->mutex), NULL) != 0){
        perror("Failed to call pthread_mutex_init in create_sorted_array\n");
        return 1;
    }
    // array->mutex = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;

	return 0;
}

int32_t recreate_sorted_array(struct sorted_array* array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b), const uint8_t method){
	size_t malloc_size;
	if(num_elements > 0){
		malloc_size = num_elements * element_size;
		array->capacity = num_elements;
	} else {
		malloc_size = SORTED_ARRAY_STEP * element_size;
		array->capacity = SORTED_ARRAY_STEP;
	}
	array->num_elements = 0;
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
	memset(array->bfr, '\0', malloc_size);
	array->to_free = 1;
	array->capacity = num_elements;
	array->num_elements = num_elements;
	array->element_size = element_size;
	array->comp = comp;
    array->num_elements_sorted = 0;
	if(array->bfr == NULL){
		printf("array->bfr == NULL\n");
		exit(1);
	}
    array->method = method;
    pthread_mutex_destroy(&(array->mutex));
    pthread_mutex_init(&(array->mutex), NULL);
	return 0;
}

void free_sorted_array(struct sorted_array* array){
	if(array->to_free == 1){
		free(array->bfr);
		array->to_free = 0;
	}
    pthread_mutex_destroy(&(array->mutex));
}

int32_t key_to_index_sorted_array_tree(const struct sorted_array * const array, const void * const key){
    int64_t current_index = 0;

    while(1){
        if(current_index >= array->num_elements){break;}
        int32_t cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (current_index * array->element_size)), key);
        if(cmp_result == 0){
            return current_index;
        } else if(cmp_result > 0){
            current_index = current_index * 2 + 1;
        } else if(cmp_result < 0){
            current_index = current_index * 2 + 2;
        }
    }

    return -1;
}

int32_t key_to_index_sorted_array_linear(const struct sorted_array * const array, const void * const key){
	if(array->num_elements == 0){return -1;}
	int32_t lower_bound = 0;
	int32_t upper_bound = array->num_elements_sorted - 1;

	while(lower_bound <= upper_bound){
		int32_t middle_index = lower_bound + ((upper_bound - lower_bound) >> 1);
		int32_t cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (middle_index * array->element_size)), key); // compliance
		if(cmp_result == 0){
			return middle_index;
		} else if(cmp_result < 0){
			lower_bound = middle_index + 1;
		} else if(cmp_result > 0){
			upper_bound = middle_index - 1;
		}
	}

    for(int64_t i = array->num_elements_sorted ; i < array->num_elements ; i++){
		int32_t cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (i * array->element_size)), key);
        
        if(cmp_result == 0){return i;}
    }

	return -1;
}

int32_t key_to_index_sorted_array(const struct sorted_array * const array, const void * const key){
    switch (array->method) {
        case SORTED_ARRAY_METHOD_LINEAR:
            return key_to_index_sorted_array_linear(array, key);
        case SORTED_ARRAY_METHOD_TREE:
            return key_to_index_sorted_array_tree(array, key);
        default:
            fprintf(stderr, "incorrect sorted_array method: %u\n", array->method);
            exit(1);
    }
}

void* key_to_elem_sorted_array(const struct sorted_array* const array, const void* const key){
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

int32_t sorted_array_request_more_capacity(struct sorted_array * const array){
    size_t new_size = (array->capacity + SORTED_ARRAY_STEP) * array->element_size;
    array->bfr = realloc(array->bfr, new_size);
    if(array->bfr == NULL){
    	perror("failed to realloc\n");
    	return 1;
    }
    memset((void*) (((uint64_t) array->bfr) + (array->capacity * array->element_size)), '\0', SORTED_ARRAY_STEP * array->element_size); // compliance
    array->capacity += SORTED_ARRAY_STEP;
    return 0;
}

int32_t insert_sorted_array_tree(struct sorted_array* array, void* alt, int8_t insert_mode){
	/*
	 * insert_mode == 0 -> insert whatever happens
	 * insert_mode == 1 -> insert only if not present
	 * insert_mode == 2 -> overwrite if present, otherwise insert
	 * */

    int32_t index = key_to_index_sorted_array_tree(array, alt);
    int32_t index_parent, is_left_child, is_right_child, index_child_left, index_child_right, index_current, cmp_result, direction;

    enum {
        DIRECTION_UP = 0,
        DIRECTION_DOWN = 1,
        DIRECTION_STOP = 2,
    };

    direction = DIRECTION_UP;

    if(insert_mode == 2){
        if(index != -1){
            memcpy((void*) (((uint64_t) array->bfr) + (index * array->element_size)), alt, array->element_size);
            return 0;
        }
    }

    assert(array->element_size > 0);

    void * swap_buffer = malloc(array->element_size);
    if(swap_buffer == NULL){
        perror("alloc failed\n");
        return 1;
    }

    if(!(insert_mode == 0 || insert_mode == 1)){
        fprintf(stderr, "incorrect insert_mode: %i\n", insert_mode);
	free(swap_buffer);
        return 1;
    }
    if(index == -1 || insert_mode == 0){
    	if(array->num_elements == array->capacity){
            if(sorted_array_request_more_capacity(array) != 0){
                perror("Failed to call sorted_array_request_more_capacity\n");
		free(swap_buffer);
                return 1;
            }
    	}
        
        memcpy((void*) (((uint64_t) array->bfr) + (array->num_elements * array->element_size)), alt, array->element_size);

        index_current = array->num_elements;
        if(index_current == 0){
            direction = DIRECTION_DOWN;
        }

        array->num_elements++;

        while(1){

            if(direction == DIRECTION_UP){
                index_parent = ((index_current - 1) % 2) >> 1;
                is_left_child = index_current % 2;
                is_right_child = (!is_left_child) & (index_parent != 0);
                cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (index_parent * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)));
                if((is_left_child && cmp_result < 0) || (is_right_child && cmp_result > 0)){
                    memcpy(swap_buffer, (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), array->element_size);
                    memcpy((void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_parent * array->element_size)), array->element_size);
                    memcpy((void*) (((uint64_t) array->bfr) + (index_parent * array->element_size)), swap_buffer, array->element_size);
                    index_current = index_parent;
                    if(index_current == 0){
                        direction = DIRECTION_DOWN;
                    }
                } else {
                    direction = DIRECTION_DOWN;
                }
            } else if(direction == DIRECTION_DOWN){
                index_child_left = index_current * 2 + 1;
                index_child_right = index_current * 2 + 2;
                
                // left

                if(index_child_left < array->num_elements){
                    cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (index_child_left * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)));
                    if(cmp_result > 0){
                        memcpy(swap_buffer, (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), array->element_size);
                        memcpy((void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_child_left * array->element_size)), array->element_size);
                        memcpy((void*) (((uint64_t) array->bfr) + (index_child_left * array->element_size)), swap_buffer, array->element_size);
                        index_current = index_child_left;
                        continue;
                    }
                }

                // right

                if(index_child_right < array->num_elements){
                    cmp_result = array->comp((void*) (((uint64_t) array->bfr) + (index_child_right * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)));
                    if(cmp_result < 0){
                        memcpy(swap_buffer, (void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), array->element_size);
                        memcpy((void*) (((uint64_t) array->bfr) + (index_current * array->element_size)), (void*) (((uint64_t) array->bfr) + (index_child_right * array->element_size)), array->element_size);
                        memcpy((void*) (((uint64_t) array->bfr) + (index_child_right * array->element_size)), swap_buffer, array->element_size);
                        index_current = index_child_right;
                        continue;
                    }
                }

                direction = DIRECTION_STOP;
            } else {
                break;
            }
        }
    }

    free(swap_buffer);

    return 0;
}

int32_t insert_sorted_array_linear(struct sorted_array* array, void* alt, int8_t insert_mode){
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
		} else {
            fprintf(stderr, "Unknown insert mode: %i\n", insert_mode);
            return 1;
        }
	}
	if(array->num_elements == array->capacity){
        if(sorted_array_request_more_capacity(array) != 0){
            perror("Failed to call sorted_array_request_more_capacity\n");
            return 1;
        }
	}
	memcpy((void*) (((uint64_t) array->bfr) + (array->num_elements * array->element_size)), alt, array->element_size); // compliance
	array->num_elements++;

    if(array->num_elements - array->num_elements_sorted >= SORTED_ARRAY_UNSORTED_PART_MAX_SIZE){
	    qsort(array->bfr, array->num_elements, array->element_size, array->comp);
        array->num_elements_sorted = array->num_elements;
    }

	return 0;
}

int32_t insert_sorted_array(struct sorted_array* array, void* alt, int8_t insert_mode){
    switch (array->method) {
        case SORTED_ARRAY_METHOD_LINEAR:
            return insert_sorted_array_linear(array, alt, insert_mode);
        case SORTED_ARRAY_METHOD_TREE:
            return insert_sorted_array_tree(array, alt, insert_mode);
        default:
            fprintf(stderr, "incorrect sorted_array method: %u\n", array->method);
            exit(1);
    }
}

int32_t sorted_array_default_cmp(const void* a, const void* b){
    return strcmp(((struct sorted_array_default_element*) a)->key, ((struct sorted_array_default_element*) b)->key);
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

