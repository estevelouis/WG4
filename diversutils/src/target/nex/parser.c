#ifndef NEX_PARSER_H
#define NEX_PARSER_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "nex/constants.h"
#include "nex/parser.h"
#include "graph.h"

struct phyl_node {
	unsigned char key[NEX_PHYL_NODE_KEY_BUFFER_SIZE];
	uint16_t forth_degree;
	uint16_t forth_capacity;
	struct phyl_node* neighbours;
	struct phyl_node* parent;
};

int32_t create_phyl_node(struct phyl_node* const pn){
	memset(pn->key, '\0', NEX_PHYL_NODE_KEY_BUFFER_SIZE * sizeof(unsigned char));
	pn->forth_degree = 0;
	pn->forth_capacity = NEX_PHYL_NODE_FORTH_CAPACITY_STEP;
	size_t alloc_size = pn->forth_capacity * sizeof(struct phyl_node);
	pn->neighbours = malloc(alloc_size);
	if(pn->neighbours == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(pn->neighbours, '\0', alloc_size);

	pn->parent = NULL;

	return 0;
}

int32_t request_more_capacity_phyl_node(struct phyl_node* const pn){
	uint16_t new_capacity = pn->forth_capacity + NEX_PHYL_NODE_FORTH_CAPACITY_STEP;
	size_t alloc_size = ((size_t) new_capacity) * sizeof(struct phyl_node);
	pn->neighbours = realloc(pn->neighbours, alloc_size);
	if(pn->neighbours == NULL){
		perror("realloc failed\n");
		return 1;
	}
	memset(&(pn->neighbours[pn->forth_capacity]), '\0', NEX_PHYL_NODE_FORTH_CAPACITY_STEP * sizeof(struct phyl_node));
	pn->forth_capacity = new_capacity;

	return 0;
}

void free_recursively_phyl_node(struct phyl_node* const pn){
	for(uint16_t i = 0 ; i < pn->forth_degree ; i++){
		free_recursively_phyl_node(&(pn->neighbours[i]));
	}
	free(pn->neighbours);
}

struct tree_path {
	char key[NEX_PHYL_NODE_KEY_BUFFER_SIZE];
	unsigned char keys[NEX_TREE_PATH_MAX_DEPTH][NEX_PHYL_NODE_KEY_BUFFER_SIZE];
	uint16_t length;
};

void compute_path_recursively_phyl_node(struct phyl_node* const pn, struct tree_path* all_paths, uint16_t current_depth, unsigned char* current_path_bfr, uint8_t reset_path_counter){
	static uint32_t path_counter = 0;

	if(reset_path_counter){
		path_counter = 0;
	}

	memcpy(&(current_path_bfr[current_depth * NEX_PHYL_NODE_KEY_BUFFER_SIZE]), pn->key, NEX_PHYL_NODE_KEY_BUFFER_SIZE);

	memcpy(all_paths[path_counter].key, pn->key, NEX_PHYL_NODE_KEY_BUFFER_SIZE);
	memcpy(all_paths[path_counter].keys, current_path_bfr, (current_depth + 1) * NEX_PHYL_NODE_KEY_BUFFER_SIZE);

	all_paths[path_counter].length = current_depth + 1;

	path_counter++;
	for(uint16_t i = 0 ; i < pn->forth_degree ; i++){
		compute_path_recursively_phyl_node(&(pn->neighbours[i]), all_paths, current_depth + 1, current_path_bfr, 0);
	}

	memset(&(current_path_bfr[current_depth * NEX_PHYL_NODE_KEY_BUFFER_SIZE]), '\0', NEX_PHYL_NODE_KEY_BUFFER_SIZE);
}

int32_t cmp_tree_path_keys(const void* a, const void* b){
	return strcmp(((struct tree_path*) a)->key, ((struct tree_path*) b)->key);
}

int32_t parse_nex_file(const char* const path, struct matrix** ptr_matrix, unsigned char** ptr_keys_ordered){
	uint32_t num_created_nodes = 0;

	struct phyl_node pn_root;
	if(create_phyl_node(&pn_root) != 0){
		perror("failed to call create_phyl_node\n");
		return 1;
	}
	num_created_nodes++;
	memcpy(pn_root.key, "root", 4);

	FILE* fp = fopen(path, "r");
	if(fp == NULL){
		perror("failed to call fopen\n");
		return 1;
	}

	unsigned char* bfr;
	size_t alloc_size = NEX_FILE_READ_BUFFER_SIZE * sizeof(unsigned char);
	bfr = malloc(alloc_size);
	if(bfr == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(bfr, '\0', alloc_size);

	uint8_t reached_begin_tree;
	while(fgets(bfr, NEX_FILE_READ_BUFFER_SIZE, fp)){
		if(strncmp(bfr, "BEGIN TREES;", 12) == 0){
			reached_begin_tree = 1;
			break;
		}
	}

	if(!reached_begin_tree){
		perror("could not find 'BEGIN TREES;'");
		free(bfr);
		fclose(fp);
		return 1;
	}

	uint16_t depth = 0;
	uint16_t max_depth = 0;
	
	while(1){
		if(fgets(bfr, NEX_FILE_READ_BUFFER_SIZE, fp) == NULL){
			perror("failed to call fgets\n");
			free(bfr);
			fclose(fp);
			return 1;
		}
		if(strncmp(bfr, "END;", 4) == 0){
			break;
		}

		unsigned char* ptr_tree = strstr(bfr, "tree ");
		if(ptr_tree == NULL){
			continue;
		}

		unsigned char* ptr_tree_name = ptr_tree + 5;
		unsigned char* ptr_equal = strchr(ptr_tree_name, '=');
		unsigned char* ptr_tree_name_end = ptr_equal;

		/* // DO NOT REMOVE

		while(ptr_tree_name_end >= ptr_tree_name && (ptr_tree_name_end - 1)[0] == ' '){
			ptr_tree_name_end--;
		}

		size_t tree_name_length = ptr_tree_name_end - ptr_tree_name;
		if(tree_name_length > NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1){
			tree_name_length = NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1;
		}
		
		struct phyl_node new_node;
		if(create_phyl_node(&new_node) != 0){
			perror("failed to call create_node\n");
			free(bfr);
			fclose(fp);
			return 1;
		}
		memcpy(new_node.key, ptr_tree_name, tree_name_length);

		printf("new tree: %s\n", new_node.key);

		if(pn_root.forth_degree == pn_root.forth_capacity){
			if(request_more_capacity_phyl_node(&pn_root) != 0){
				perror("failed to call request_more_capacity_phyl_node\n");
				free(bfr);
				fclose(fp);
				return 1;
			}
		}

		new_node.parent = &pn_root;

		pn_root.neighbours[pn_root.forth_degree] = new_node;
		pn_root.forth_degree++;

		*/

		// start filling the local tree

		struct phyl_node* pn_current = &pn_root;

		unsigned char* ptr_char = ptr_equal;
		while(!((*ptr_char) == '\0' || (*ptr_char) == '(')){
			ptr_char++;
		}

		while(1){
			while((*ptr_char) != ';'){
				struct phyl_node local_node;
				size_t bytes_to_cpy;
				size_t key_len;
				unsigned char* end = NULL;
				unsigned char* next_parenthesis_closing = NULL;
				unsigned char* next_parenthesis_opening = NULL;
				unsigned char* next_comma = NULL;
				unsigned char* next_semi_colon = NULL;
				switch((*ptr_char)){
					case '(':
						struct phyl_node local_node;
						if(create_phyl_node(&local_node) != 0){
							perror("failed to call create_phyl_node\n");
							free(bfr);
							fclose(fp);
							return 1;
						}
						num_created_nodes++;
						if(pn_current->forth_degree == pn_current->forth_capacity){
							if(request_more_capacity_phyl_node(pn_current) != 0){
								perror("failed to call request_mode_capacity_phyl_node\n");
								free(bfr);
								fclose(fp);
								return 1;
							}
						}

						local_node.parent = pn_current;
						pn_current->neighbours[pn_current->forth_degree] = local_node;
						pn_current->forth_degree++;

						pn_current = &(pn_current->neighbours[pn_current->forth_degree - 1]);
						depth++;
						if(depth > max_depth){
							max_depth = depth;
						}

						ptr_char++;

						break;
					case ')':
						ptr_char++;

						end = NULL;
						next_parenthesis_closing = strchr(ptr_char, ')');
						next_parenthesis_opening = strchr(ptr_char, '(');
						next_comma = strchr(ptr_char, ',');
						next_semi_colon = strchr(ptr_char, ';');

						if(next_parenthesis_closing != NULL && (end == NULL || next_parenthesis_closing < end)){end = next_parenthesis_closing;}
						if(next_parenthesis_opening != NULL && (end == NULL || next_parenthesis_opening < end)){end = next_parenthesis_opening;}
						if(next_comma != NULL && (end == NULL || next_comma < end)){end = next_comma;}
						if(next_semi_colon != NULL && (end == NULL || next_semi_colon < end)){end = next_semi_colon;}

						key_len;
						if(end == NULL){
							key_len = strlen(ptr_char);
						} else {
							key_len = end - ptr_char;
						}
						bytes_to_cpy = key_len;
						if(bytes_to_cpy > NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1){
							bytes_to_cpy = NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1;
						}
						memcpy(pn_current->key, ptr_char, bytes_to_cpy);

						if(end == NULL){
							ptr_char += key_len;
						} else {
							ptr_char = end;
						}
						
						pn_current = pn_current->parent;
						depth--;

						break;
					case ',':
						ptr_char++;

						break;
					case ';':
						break;
					default:
						end = NULL;
						next_parenthesis_closing = strchr(ptr_char, ')');
						next_parenthesis_opening = strchr(ptr_char, '(');
						next_comma = strchr(ptr_char, ',');
						next_semi_colon = strchr(ptr_char, ';');

						if(next_parenthesis_closing != NULL && (end == NULL || next_parenthesis_closing < end)){end = next_parenthesis_closing;}
						if(next_parenthesis_opening != NULL && (end == NULL || next_parenthesis_opening < end)){end = next_parenthesis_opening;}
						if(next_comma != NULL && (end == NULL || next_comma < end)){end = next_comma;}
						if(next_semi_colon != NULL && (end == NULL || next_semi_colon < end)){end = next_semi_colon;}

						key_len;
						if(end == NULL){
							key_len = strlen(ptr_char);
						} else {
							key_len = end - ptr_char;
						}
						bytes_to_cpy = key_len;
						if(bytes_to_cpy > NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1){
							bytes_to_cpy = NEX_PHYL_NODE_KEY_BUFFER_SIZE - 1;
						}

						if(key_len == 0){
							perror("empty key\n");
							free(bfr);
							fclose(fp);
							return 1;
						}

						struct phyl_node local_node_;
						if(create_phyl_node(&local_node_) != 0){
							perror("failed to call create_phyl_node\n");
							free(bfr);
							fclose(fp);
							return 1;
						}
						num_created_nodes++;

						memcpy(local_node_.key, ptr_char, bytes_to_cpy);

						if(pn_current->forth_degree == pn_current->forth_capacity){
							if(request_more_capacity_phyl_node(pn_current) != 0){
								perror("failed to call request_more_capacity_phyl_node\n");
								free(bfr);
								fclose(fp);
								return 1;
							}
						}

						pn_current->neighbours[pn_current->forth_degree] = local_node_;
						pn_current->forth_degree++;

						if(end == NULL){
							ptr_char += key_len;
						} else {
							ptr_char = end;
						}

						break;
				}
			}
			if((*ptr_char) == '\0'){
				break;
			}
			if((*ptr_char) != ';'){
				ptr_char = fgets(bfr, NEX_FILE_READ_BUFFER_SIZE, fp);
				if(ptr_char == NULL){
					perror("failed to call fgets\n");
					free(bfr);
					fclose(fp);
					return 1;
				}
			} else {
				break;
			}
		}
	}

	free(bfr);
	fclose(fp);

	alloc_size = NEX_TREE_PATH_MAX_DEPTH * NEX_PHYL_NODE_KEY_BUFFER_SIZE;
	unsigned char* current_path_bfr = malloc(alloc_size);
	if(current_path_bfr == NULL){
		perror("failed to malloc\n");
		free_recursively_phyl_node(&pn_root);
		return 1;
	}
	memset(current_path_bfr, '\0', alloc_size);

	alloc_size = (num_created_nodes) * sizeof(struct tree_path);
	struct tree_path* all_paths = malloc(alloc_size);
	if(all_paths == NULL){
		perror("failed to malloc\n");
		free_recursively_phyl_node(&pn_root);
		return 1;
	}
	memset(all_paths, '\0', alloc_size);

	compute_path_recursively_phyl_node(&pn_root, all_paths, 0, current_path_bfr, 1);
	
	free_recursively_phyl_node(&pn_root);

	struct matrix m;
	if(create_matrix(&m, num_created_nodes, num_created_nodes, FP32) != 0){
		perror("failed to call create_matrix\n");
		return 1;
	}

	qsort(all_paths, num_created_nodes, sizeof(struct tree_path), cmp_tree_path_keys);

	unsigned char* keys_ordered;
	alloc_size = num_created_nodes * NEX_PHYL_NODE_KEY_BUFFER_SIZE;
	keys_ordered = (unsigned char*) malloc(alloc_size);
	if(keys_ordered == NULL){
		perror("malloc failed\n");
		return 1;
	}

	for(uint32_t i = 0 ; i < num_created_nodes ; i++){
		/* // DO NOT REMOVE
		printf("key=%s;\t", all_paths[i].key);
		for(uint16_t j = 0 ; j < all_paths[i].length ; j++){
			printf("%s ", all_paths[i].keys[j]);
		}
		printf("\n");
		*/

		memcpy(&(keys_ordered[i * NEX_PHYL_NODE_KEY_BUFFER_SIZE]), all_paths[i].key, NEX_PHYL_NODE_KEY_BUFFER_SIZE);

		for(uint32_t j = 0 ; j < num_created_nodes ; j++){
			uint16_t k_limit = all_paths[i].length;
			uint16_t k = 0;
			if(all_paths[j].length < k_limit){
				k_limit = all_paths[j].length;
			}
			for(k = 0 ; k < k_limit ; k++){
				int32_t strcmp_result = strcmp(all_paths[i].keys[k], all_paths[j].keys[k]);
				if(strcmp_result != 0){
					break;
				}
			}

			uint16_t delta;
			if(all_paths[j].length > all_paths[i].length){
				delta = all_paths[j].length - k;
			} else {
				delta = all_paths[i].length - k;
			}

			m.bfr.fp32[i * m.b + j] = (float) delta;
		}
	}

	(*ptr_matrix) = &m;
	(*ptr_keys_ordered) = keys_ordered;

	free(current_path_bfr);
	free(all_paths);

	return 0;
}

#endif
