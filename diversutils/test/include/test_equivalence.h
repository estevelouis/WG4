#ifndef TEST_EQUIVALENCE_H
#define TEST_EQUIVALENCE_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_general.h"
#include "graph.h"
#include "dfunctions.h"
#include "distances.h"

int32_t test_equivalence_entropy(void){
	int32_t result = 0;
	
	for(double alpha = 0.0 ; alpha < 3.0 ; alpha += 0.25){
		for(uint64_t n = 1 ; n < 1e4 ; n *= 10){
			struct graph g;
			if(create_graph(&g, n, 0, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
			if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}
	
			for(uint64_t i = 0 ; i < n ; i++){
				g.nodes[i] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
			}
	
			double res_shannon_weaver_entropy, res_shannon_weaver_hill_number;
			double res_renyi_entropy, res_renyi_hill_number;
			double res_patil_taillie_entropy, res_patil_taillie_hill_number;
			double res_q_logarithmic_entropy, res_q_logarithmic_hill_number;
			const size_t log_bfr_size = 256;
			char log_bfr[log_bfr_size];
		
			compute_graph_relative_proportions(&g);
			shannon_weaver_entropy_from_graph(&g, &res_shannon_weaver_entropy, &res_shannon_weaver_hill_number);
			renyi_entropy_from_graph(&g, &res_renyi_entropy, &res_renyi_hill_number, alpha);
			patil_taillie_entropy_from_graph(&g, &res_patil_taillie_entropy, &res_patil_taillie_hill_number, alpha - 1.0);
			q_logarithmic_entropy_from_graph(&g, &res_q_logarithmic_entropy, &res_q_logarithmic_hill_number, alpha);
			res_shannon_weaver_entropy = round(res_shannon_weaver_entropy * 1000000.0) / 1000000.0;
			res_renyi_entropy = round(res_renyi_entropy * 1000000.0) / 1000000.0;
			res_patil_taillie_entropy = round(res_patil_taillie_entropy * 1000000.0) / 1000000.0;
			res_q_logarithmic_entropy = round(res_q_logarithmic_entropy * 1000000.0) / 1000000.0;
			memset(log_bfr, '\0', log_bfr_size);
			if(res_patil_taillie_entropy == res_q_logarithmic_entropy){
				snprintf(log_bfr, log_bfr_size, "PT = q-log equivalence test (even distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_patil_taillie_entropy, res_q_logarithmic_entropy);
				info_format(__FILE__, __func__, __LINE__, log_bfr);
			} else {
				snprintf(log_bfr, log_bfr_size, "PT = q-log equivalence test (even distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_patil_taillie_entropy, res_q_logarithmic_entropy);
				error_format(__FILE__, __func__, __LINE__, log_bfr);
				result = 1;
			}
			if(n > 1){
				if(alpha == 1.0){
					if(res_shannon_weaver_entropy == res_q_logarithmic_entropy){
						snprintf(log_bfr, log_bfr_size, "SW = q-log equivalence test (even distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						info_format(__FILE__, __func__, __LINE__, log_bfr);
					} else {
						snprintf(log_bfr, log_bfr_size, "SW = q-log equivalence test (even distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						error_format(__FILE__, __func__, __LINE__, log_bfr);
						result = 1;
					}
				} else {
					if(res_shannon_weaver_entropy != res_q_logarithmic_entropy){
						snprintf(log_bfr, log_bfr_size, "SW != q-log equivalence test (even distribution / alpha = %f / %lu elements): OK (%f != %f)", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						info_format(__FILE__, __func__, __LINE__, log_bfr);
					} else {
						snprintf(log_bfr, log_bfr_size, "SW != q-log equivalence test (even distribution / alpha = %f / %lu elements): FAIL (!(%f != %f))", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						error_format(__FILE__, __func__, __LINE__, log_bfr);
						result = 1;
					}
				}
			}
		
			// uneven
			if(n > 1 && alpha > 0.0){
				g.nodes[0].absolute_proportion = 2;
				shannon_weaver_entropy_from_graph(&g, &res_shannon_weaver_entropy, &res_shannon_weaver_hill_number);
				renyi_entropy_from_graph(&g, &res_renyi_entropy, &res_renyi_hill_number, alpha);
				patil_taillie_entropy_from_graph(&g, &res_patil_taillie_entropy, &res_patil_taillie_hill_number, alpha - 1.0);
				q_logarithmic_entropy_from_graph(&g, &res_q_logarithmic_entropy, &res_q_logarithmic_hill_number, alpha);
				res_shannon_weaver_entropy = round(res_shannon_weaver_entropy * 1000000.0) / 1000000.0;
				res_renyi_entropy = round(res_renyi_entropy * 1000000.0) / 1000000.0;
				res_patil_taillie_entropy = round(res_patil_taillie_entropy * 1000000.0) / 1000000.0;
				res_q_logarithmic_entropy = round(res_q_logarithmic_entropy * 1000000.0) / 1000000.0;
				memset(log_bfr, '\0', log_bfr_size);
				if(res_patil_taillie_entropy == res_q_logarithmic_entropy){
					snprintf(log_bfr, log_bfr_size, "PT = q-log equivalence test (uneven distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_patil_taillie_entropy, res_q_logarithmic_entropy);
					info_format(__FILE__, __func__, __LINE__, log_bfr);
				} else {
					snprintf(log_bfr, log_bfr_size, "PT = q-log equivalence test (uneven distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_patil_taillie_entropy, res_q_logarithmic_entropy);
					error_format(__FILE__, __func__, __LINE__, log_bfr);
					result = 1;
				}
				if(alpha == 1.0){
					if(res_shannon_weaver_entropy == res_q_logarithmic_entropy){
						snprintf(log_bfr, log_bfr_size, "SW = q-log equivalence test (uneven distribution / alpha = %f / %lu elements): OK (%f == %f)", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						info_format(__FILE__, __func__, __LINE__, log_bfr);
					} else {
						snprintf(log_bfr, log_bfr_size, "SW = q-log equivalence test (uneven distribution / alpha = %f / %lu elements): FAIL (!(%f == %f))", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						error_format(__FILE__, __func__, __LINE__, log_bfr);
						result = 1;
					}
				} else {
					if(res_shannon_weaver_entropy != res_q_logarithmic_entropy){
						snprintf(log_bfr, log_bfr_size, "SW != q-log equivalence test (uneven distribution / alpha = %f / %lu elements): OK (%f != %f)", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						info_format(__FILE__, __func__, __LINE__, log_bfr);
					} else {
						snprintf(log_bfr, log_bfr_size, "SW != q-log equivalence test (uneven distribution / alpha = %f / %lu elements): FAIL (!(%f != %f))", alpha, n, res_shannon_weaver_entropy, res_q_logarithmic_entropy);
						error_format(__FILE__, __func__, __LINE__, log_bfr);
						result = 1;
					}
				}
			}
		
			free_graph(&g);
		}
	}
	return result;
}

#endif
