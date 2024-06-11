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
	void* bfr;
	size_t element_size;
	int64_t num_elements;
    int64_t num_elements_sorted;
	int64_t capacity;
	int32_t (*comp)(const void* a, const void* b);
	int32_t to_free;
    pthread_mutex_t mutex;
};

int32_t create_sorted_array_default_element(struct sorted_array_default_element* elem);
int32_t create_sorted_array_str_int_element(struct sorted_array_str_int_element* elem);
int32_t create_sorted_array_int_int_element(struct sorted_array_int_int_element* elem);
int32_t create_sorted_array_int_str_element(struct sorted_array_int_str_element* elem);
int32_t create_sorted_array(struct sorted_array * const array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b));
int32_t recreate_sorted_array(struct sorted_array* array, int64_t num_elements, size_t element_size, int32_t (*comp)(const void* a, const void* b));
void free_sorted_array(struct sorted_array* array);
int32_t key_to_index_sorted_array(const struct sorted_array * const array, const void* const key);
void* key_to_elem_sorted_array(const struct sorted_array* const array, const void* const key);
int32_t bulk_insert_sorted_array(struct sorted_array* array, void* alt_array, size_t n);
int32_t insert_sorted_array(struct sorted_array* array, void* alt, int8_t insert_mode);
int32_t sorted_array_default_cmp(const void* a, const void* b);
int32_t sorted_array_str_int_cmp(const void* a, const void* b);
int32_t sorted_array_int_int_cmp(const void* a, const void* b);
int32_t sorted_array_int_str_cmp(const void* a, const void* b);

#endif
