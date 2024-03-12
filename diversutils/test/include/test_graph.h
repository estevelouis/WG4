#ifndef TEST_GRAPH_H
#define TEST_GRAPH_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_general.h"
#include "graph.h"
#include "dfunctions.h"
#include "distances.h"

int32_t test_compute_graph_relative_proportions(void){
	struct graph g;
	if(create_graph(&g, 4, 100, 0) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call create_graph"); return 1;}
	if(request_more_capacity_graph(&g) != 0){error_format(__FILE__, __func__, __LINE__, "failed to call request_more_capacity_graph"); free_graph(&g); return 1;}

	g.nodes[0] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
	g.nodes[1] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
	g.nodes[2] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };
	g.nodes[3] = (struct graph_node) { .num_dimensions = 0, .already_considered = 0, .relative_proportion = 0.0, .absolute_proportion = 1, .word2vec_entry_pointer = NULL, .neighbours = NULL, .capacity_neighbours = 0, .num_neighbours = 0, .mutex_local_node = PTHREAD_MUTEX_INITIALIZER, .vector = NULL };

	int32_t result = 0;

	compute_graph_relative_proportions(&g);
	if(g.nodes[0].relative_proportion != 0.25 || g.nodes[1].relative_proportion != 0.25 || g.nodes[2].relative_proportion != 0.25 || g.nodes[3].relative_proportion != 0.25){
		error_format(__FILE__, __func__, __LINE__, "Compute graph relative proportions (even distribution): FAIL");
		result = 1;
	} else {
		info_format(__FILE__, __func__, __LINE__, "Compute graph relative proportions (even distribution): OK");
	}

	g.nodes[0].absolute_proportion = 2;
	compute_graph_relative_proportions(&g);
	if(g.nodes[0].relative_proportion != 0.4 || g.nodes[1].relative_proportion != 0.2 || g.nodes[2].relative_proportion != 0.2 || g.nodes[3].relative_proportion != 0.2){
		error_format(__FILE__, __func__, __LINE__, "Compute graph relative proportions (uneven distribution): FAIL");
		result = 1;
	} else {
		info_format(__FILE__, __func__, __LINE__, "Compute graph relative proportions (uneven distribution): OK");
	}

	free_graph(&g);
	return result;
}

#endif
