#ifndef JSONL_LOAD_H
#define JSONL_LOAD_H

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>

#include "graph.h"
#include "sorted_array/array.h"

int32_t jsonl_to_graph(const uint64_t i, const char * const filename, struct measurement_configuration * const mcfg, struct measurement_structure_references * const sref, struct measurement_mutables * const mmut);

void * jsonl_to_graph_thread(void * args);

#endif
