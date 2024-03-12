#ifndef TEST_ENTROPY_H
#define TEST_ENTROPY_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_general.h"
#include "graph.h"
#include "dfunctions.h"
#include "distances.h"

int32_t test_shannon_weaver_entropy(void){
	int32_t result = 0;
	
	for(uint64_t n = 1 ; n < 1e4 ; n *= 10){
		struct graph g;
		if(create_graph(&g, n, 0, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
		if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}

		for(uint64_t i = 0 ; i < n ; i++){
			g.nodes[i] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
		}

		double res_entropy;
		double res_hill_number;
		double target;
		const size_t log_bfr_size = 256;
		char log_bfr[log_bfr_size];
	
		compute_graph_relative_proportions(&g);
		shannon_weaver_entropy_from_graph(&g, &res_entropy, &res_hill_number);
		res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
		target = log((double) n);
		target = round(target * 1000000.0) / 1000000.0;
		memset(log_bfr, '\0', log_bfr_size);
		if(res_entropy == target){
			snprintf(log_bfr, log_bfr_size, "SW entropy test (even distribution / %lu elements): OK (%f == %f)", n, res_entropy, target);
			info_format(__FILE__, __func__, __LINE__, log_bfr);
		} else {
			snprintf(log_bfr, log_bfr_size, "SW entropy test (even distribution / %lu elements): FAIL (!(%f == %f))", n, res_entropy, target);
			error_format(__FILE__, __func__, __LINE__, log_bfr);
			result = 1;
		}
	
		if(n > 1){
			g.nodes[0].absolute_proportion = 2;
			compute_graph_relative_proportions(&g);
			shannon_weaver_entropy_from_graph(&g, &res_entropy, &res_hill_number);
			res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
			memset(log_bfr, '\0', log_bfr_size);
			if(res_entropy < target){
				snprintf(log_bfr, log_bfr_size, "SW entropy test (uneven distribution / %lu elements): OK (%f < %f)", n, res_entropy, target);
				info_format(__FILE__, __func__, __LINE__, log_bfr);
			} else {
				snprintf(log_bfr, log_bfr_size, "SW entropy test (uneven distribution / %lu elements): FAIL (!(%f < %f))", n, res_entropy, target);
				error_format(__FILE__, __func__, __LINE__, log_bfr);
				result = 1;
			}
		}
	
		free_graph(&g);
	}
	return result;
}

int32_t test_renyi_entropy(void){
	int32_t result = 0;
	
	for(double alpha = 0.0 ; alpha < 3.0 ; alpha += 0.25){
		for(uint64_t n = 1 ; n < 1e4 ; n *= 10){
			struct graph g;
			if(create_graph(&g, n, 0, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
			if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}
	
			for(uint64_t i = 0 ; i < n ; i++){
				g.nodes[i] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
			}
	
			double res_entropy;
			double res_hill_number;
			double target;
			const size_t log_bfr_size = 256;
			char log_bfr[log_bfr_size];
		
			compute_graph_relative_proportions(&g);
			renyi_entropy_from_graph(&g, &res_entropy, &res_hill_number, alpha);
			res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
			target = log((double) n);
			target = round(target * 1000000.0) / 1000000.0;
			memset(log_bfr, '\0', log_bfr_size);
			if(res_entropy == target){
				snprintf(log_bfr, log_bfr_size, "R entropy test (even distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_entropy, target);
				info_format(__FILE__, __func__, __LINE__, log_bfr);
			} else {
				snprintf(log_bfr, log_bfr_size, "R entropy test (even distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_entropy, target);
				error_format(__FILE__, __func__, __LINE__, log_bfr);
				result = 1;
			}
		
			if(n > 1 && alpha > 0.0){
				g.nodes[0].absolute_proportion = 2;
				compute_graph_relative_proportions(&g);
				renyi_entropy_from_graph(&g, &res_entropy, &res_hill_number, alpha);
				res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
				memset(log_bfr, '\0', log_bfr_size);
				if(res_entropy < target){
					snprintf(log_bfr, log_bfr_size, "R entropy test (uneven distribution / alpha = %f / %lu elements): OK (%f < %f)", alpha, n, res_entropy, target);
					info_format(__FILE__, __func__, __LINE__, log_bfr);
				} else {
					snprintf(log_bfr, log_bfr_size, "R entropy test (uneven distribution / alpha = %f / %lu elements): FAIL (!(%f < %f))", alpha, n, res_entropy, target);
					error_format(__FILE__, __func__, __LINE__, log_bfr);
					result = 1;
				}
			}
		
			free_graph(&g);
		}
	}
	return result;
}

int32_t test_patil_taillie_entropy(void){
	int32_t result = 0;
	
	for(double alpha = -1.0 ; alpha <= 2.0 ; alpha += 0.25){
		for(uint64_t n = 1 ; n < 1e4 ; n *= 10){
			struct graph g;
			if(create_graph(&g, n, 0, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
			if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}
	
			for(uint64_t i = 0 ; i < n ; i++){
				g.nodes[i] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
			}
	
			double res_entropy;
			double res_hill_number;
			double target;
			const size_t log_bfr_size = 256;
			char log_bfr[log_bfr_size];
		
			compute_graph_relative_proportions(&g);
			patil_taillie_entropy_from_graph(&g, &res_entropy, &res_hill_number, alpha);
			res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
			if(alpha == 0.0){
				target = log((double) n);
			} else {
				target = (1.0 - pow(1.0 / ((double) n), alpha)) / alpha;
			}
			target = round(target * 1000000.0) / 1000000.0;
			memset(log_bfr, '\0', log_bfr_size);
			if(res_entropy == target){
				snprintf(log_bfr, log_bfr_size, "PT entropy test (even distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_entropy, target);
				info_format(__FILE__, __func__, __LINE__, log_bfr);
			} else {
				snprintf(log_bfr, log_bfr_size, "PT entropy test (even distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_entropy, target);
				error_format(__FILE__, __func__, __LINE__, log_bfr);
				result = 1;
			}
		
			if(n > 1 && alpha > 0.0){
				g.nodes[0].absolute_proportion = 10;
				compute_graph_relative_proportions(&g);
				patil_taillie_entropy_from_graph(&g, &res_entropy, &res_hill_number, alpha);
				res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
				memset(log_bfr, '\0', log_bfr_size);
				if(res_entropy < target){
					snprintf(log_bfr, log_bfr_size, "PT entropy test (uneven distribution / alpha = %f / %lu elements): OK (%f < %f)", alpha, n, res_entropy, target);
					info_format(__FILE__, __func__, __LINE__, log_bfr);
				} else {
					snprintf(log_bfr, log_bfr_size, "PT entropy test (uneven distribution / alpha = %f / %lu elements): FAIL (!(%f < %f))", alpha, n, res_entropy, target);
					error_format(__FILE__, __func__, __LINE__, log_bfr);
					result = 1;
				}
			}
		
			free_graph(&g);
		}
	}
	return result;
}

int32_t test_q_logarithmic_entropy(void){
	int32_t result = 0;
	
	for(double q = 0.0 ; q <= 3.0 ; q += 0.25){
		for(uint64_t n = 1 ; n < 1e4 ; n *= 10){
			struct graph g;
			if(create_graph(&g, n, 0, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
			if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}
	
			for(uint64_t i = 0 ; i < n ; i++){
				g.nodes[i] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
			}
	
			double res_entropy;
			double res_hill_number;
			double target;
			const size_t log_bfr_size = 256;
			char log_bfr[log_bfr_size];
		
			compute_graph_relative_proportions(&g);
			q_logarithmic_entropy_from_graph(&g, &res_entropy, &res_hill_number, q);
			res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
			target = q_logarithm((double) n, q);
			target = round(target * 1000000.0) / 1000000.0;
			memset(log_bfr, '\0', log_bfr_size);
			if(res_entropy == target){
				snprintf(log_bfr, log_bfr_size, "q-log entropy test (even distribution / q = %f / %lu elements): OK (%f == %f)", q, n, res_entropy, target);
				info_format(__FILE__, __func__, __LINE__, log_bfr);
			} else {
				snprintf(log_bfr, log_bfr_size, "q-log entropy test (even distribution / q = %f / %lu elements): FAIL (!(%f == %f))", q, n, res_entropy, target);
				error_format(__FILE__, __func__, __LINE__, log_bfr);
				result = 1;
			}
		
			if(n > 1 && q > 0.0){
				g.nodes[0].absolute_proportion = 10;
				compute_graph_relative_proportions(&g);
				q_logarithmic_entropy_from_graph(&g, &res_entropy, &res_hill_number, q);
				res_entropy = round(res_entropy * 1000000.0) / 1000000.0;
				memset(log_bfr, '\0', log_bfr_size);
				if(res_entropy < target){
					snprintf(log_bfr, log_bfr_size, "q-log entropy test (uneven distribution / q = %f / %lu elements): OK (%f < %f)", q, n, res_entropy, target);
					info_format(__FILE__, __func__, __LINE__, log_bfr);
				} else {
					snprintf(log_bfr, log_bfr_size, "q-log entropy test (uneven distribution / q = %f / %lu elements): FAIL (!(%f < %f))", q, n, res_entropy, target);
					error_format(__FILE__, __func__, __LINE__, log_bfr);
					result = 1;
				}
			}
		
			free_graph(&g);
		}
	}
	return result;
}

#endif
