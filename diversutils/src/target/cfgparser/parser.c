#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "cfgparser/parser.h"

int32_t create_cfg_from_file(struct cfg* const c, const char* const path){
	size_t alloc_size = CFGPARSER_BFR_SIZE_READ * sizeof(char);

	char* bfr = (char*) malloc(alloc_size);
	if(bfr == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(bfr, '\0', alloc_size);

	FILE* f = fopen(path, "r");
	if(f == NULL){
		fprintf(stderr, "failed to open %s\n", path);
		free(bfr);
		return 1;
	}

	c->keys = NULL;
	c->values = NULL;
	c->num_entries = 0;

	while(fgets(bfr, CFGPARSER_BFR_SIZE_READ, f)){
		if(bfr[0] == '#' || bfr[0] == '\n' || bfr[0] == '\r'){
			continue;
		}
		char* equal_pointer = strchr(bfr, '=');
		if(equal_pointer == NULL){
			fprintf(stderr, "failed to parse a line: %s\n", bfr);
			continue;
		}
		char* end_of_key_pointer = equal_pointer;
		while(end_of_key_pointer != bfr && (end_of_key_pointer - 1)[0] == ' '){
			end_of_key_pointer--;
		}
		size_t key_size = end_of_key_pointer - bfr;
		char* value_start_pointer = equal_pointer + 1;
		while(value_start_pointer[0] == ' '){
			value_start_pointer++;
		}
		size_t value_size = strlen(value_start_pointer);
		while(value_size > 0 && value_start_pointer[value_size - 1] == '\n'){
			value_size--;
		}


		size_t local_alloc_size = (c->num_entries + 1) * sizeof(char*);
		c->keys = realloc(c->keys, local_alloc_size);
		if(c->keys == NULL){
			perror("failed to realloc\n");
			free(bfr);
			return 1;
		}
		// memset((void*) &(c->keys[c->num_entries]), '\0', sizeof(char*));

		// size_t local_alloc_size = (c->num_entries + 1) * sizeof(char*);
		c->values = realloc(c->values, local_alloc_size);
		if(c->values == NULL){
			perror("failed to realloc\n");
			free(bfr);
			return 1;
		}

		local_alloc_size = (key_size + 1) * sizeof(char);
		c->keys[c->num_entries] = malloc(local_alloc_size);
		if(c->keys[c->num_entries] == NULL){
			perror("failed to malloc for key\n");
			free(bfr);
			return 1;
		}
		strncpy(c->keys[c->num_entries], bfr, key_size);
		c->keys[c->num_entries][key_size] = '\0';

		local_alloc_size = (value_size + 1) * sizeof(char);
		c->values[c->num_entries] = malloc(local_alloc_size);
		if(c->values[c->num_entries] == NULL){
			perror("failed to malloc for value\n");
			free(bfr);
			return 1;
		}
		strncpy(c->values[c->num_entries], value_start_pointer, value_size);
		c->values[c->num_entries][value_size] = '\0';

		c->num_entries++;		
	}

	fclose(f);
	free(bfr);

	for(uint32_t i = 0 ; i < c->num_entries - 1 ; i++){
		uint32_t j = i + 1;
		while(j >= 1){
			int32_t str_cmp_result = strcmp(c->keys[j-1], c->keys[j]);
			if(str_cmp_result > 0){
				char* placeholder = c->keys[j-1];
				c->keys[j-1] = c->keys[j];
				c->keys[j] = placeholder;
				placeholder = c->values[j-1];
				c->values[j-1] = c->values[j];
				c->values[j] = placeholder;
				j--;
			} else {
				break;
			}
		}
	}

	for(uint32_t i = 0 ; i < c->num_entries ; i++){
		printf("[cfg] %s: %s\n", c->keys[i], c->values[i]);
	}

	return 0;
}

void free_cfg(struct cfg* const c){
	for(uint32_t i = 0 ; i < c->num_entries ; i++){
		free(c->keys[i]);
		free(c->values[i]);
	}
	free(c->keys);
	free(c->values);
}

char* cfg_get_value(const struct cfg* const c, const char* const key){
	int32_t lower_bound = 0;
	int32_t higher_bound = (int32_t) c->num_entries;
	while(lower_bound <= higher_bound){
		int32_t index = lower_bound + floor(((float) (higher_bound - lower_bound)) / 2.0f);
		int32_t str_cmp_result = strcmp(key, c->keys[index]);
		if(str_cmp_result > 0){
			lower_bound = index + 1;
		} else if(str_cmp_result < 0){
			higher_bound = index - 1;
		} else {
			return c->values[index];
		}
	}
	return NULL;
}

