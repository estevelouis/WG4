#ifndef CUPT_LOAD_H
#define CUPT_LOAD_H

#include <pthread.h>
#include <stdint.h>

#include "graph.h"
#include "sorted_array/array.h"

int32_t cupt_to_graph(const uint64_t i, const char * const filename, struct measurement_configuration * const mcfg, struct measurement_structure_references * const sref, struct measurement_mutables * const mmut);

void * cupt_to_graph_thread(void * args);

#endif
