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

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

#include "graph.h"
#include "cfgparser/parser.h"

#include "measurement.h"

enum {
	ID_ENTROPY_SHANNON_WEAVER,
	ID_ENTROPY_Q_LOGARITHMIC,
	ID_ENTROPY_PATIL_TAILLIE,
	ID_ENTROPY_RENYI,
	ID_ENTROPY_GOOD,
};

static uint32_t num_graphs = 0;
static struct graph* global_graphs = NULL;
static uint8_t* global_graphs_freed = NULL;
static struct cfg* configurations = NULL;

static PyObject* interface_measurement_from_cfg(PyObject* self, PyObject* args){
	char* config_path;

	if(!PyArg_ParseTuple(args, "s", &config_path)){
		return NULL;
	}

	struct cfg c;
	if(create_cfg_from_file(&c, config_path) != 0){
		perror("failed to call create_cfg_from_file\n");
		return NULL;
	}

	if(measurement_from_cfg(&c) != 0){
		perror("failed to call measurement_from_cfg\n");
		return NULL;
	}

	PyObject* res;
	res = Py_BuildValue("i", 0);
	return res;
}

static PyObject* interface_create_empty_graph(PyObject* self, PyObject* args){
	PyObject* res;
	size_t alloc_size;
	int32_t num_nodes;
	int32_t num_dimensions;

	if(!PyArg_ParseTuple(args, "ii", &num_nodes, &num_dimensions)){
		return NULL;
	}

	alloc_size = (num_graphs + 1) * sizeof(struct graph);
	global_graphs = realloc(global_graphs, alloc_size);
	if(global_graphs == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_graphs[num_graphs]), '\0', sizeof(struct graph));

	alloc_size = (num_graphs + 1) * sizeof(uint8_t);
	global_graphs_freed = realloc(global_graphs_freed, alloc_size);
	if(global_graphs_freed == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_graphs_freed[num_graphs]), '\0', sizeof(uint8_t));

	alloc_size = (num_graphs + 1) * sizeof(struct cfg);
	configurations = realloc(configurations, alloc_size);
	if(configurations == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(configurations[num_graphs]), '\0', sizeof(struct cfg));

	if(create_graph(&(global_graphs[num_graphs]), num_nodes, (int16_t) num_dimensions, FP32) != 0){
		perror("failed to call create_graph\n");
		goto exit_failure;
	}

	res = Py_BuildValue("i", num_graphs);

	num_graphs++;

	return res;

	exit_failure:
	res = Py_BuildValue("i", -1);
	return res;
}

static PyObject* interface_create_graph(PyObject* self, PyObject* args){
	PyObject* res;
	size_t alloc_size;
	int32_t num_nodes;
	int32_t num_dimensions;
	char* config_path;

	if(!PyArg_ParseTuple(args, "iis", &num_nodes, &num_dimensions, &config_path)){
		return NULL;
	}

	alloc_size = (num_graphs + 1) * sizeof(struct graph);
	global_graphs = realloc(global_graphs, alloc_size);
	if(global_graphs == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_graphs[num_graphs]), '\0', sizeof(struct graph));

	alloc_size = (num_graphs + 1) * sizeof(uint8_t);
	global_graphs_freed = realloc(global_graphs_freed, alloc_size);
	if(global_graphs_freed == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_graphs_freed[num_graphs]), '\0', sizeof(uint8_t));

	alloc_size = (num_graphs + 1) * sizeof(struct cfg);
	configurations = realloc(configurations, alloc_size);
	if(configurations == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(configurations[num_graphs]), '\0', sizeof(struct cfg));

	if(create_graph(&(global_graphs[num_graphs]), num_nodes, (int16_t) num_dimensions, FP32) != 0){
		perror("failed to call create_graph\n");
		goto exit_failure;
	}

	if(create_cfg_from_file(&(configurations[num_graphs]), config_path) != 0){
		perror("failed to call create_cfg_from_file\n");
		goto exit_failure;
	}

	res = Py_BuildValue("i", num_graphs);

	num_graphs++;

	return res;

	exit_failure:
	res = Py_BuildValue("i", -1);
	return res;
}

static PyObject* interface_free_graph(PyObject* self, PyObject* args){
	int32_t index;

	if(!PyArg_ParseTuple(args, "i", &index)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_graphs){
		perror("index too high\n");
		return NULL;
	}

	if(global_graphs_freed[index] == 0){
		free_graph(&(global_graphs[index]));
		free_cfg(&(configurations[index]));
		global_graphs_freed[index] = 1;
	} else {
		printf("graph already freed\n");
	}

	PyObject* res = Py_BuildValue("i", 0);

	return res;
}

static PyObject* interface_add_node(PyObject* self, PyObject* args){
	int32_t index;
	int32_t absolute_proportion;

	if(!PyArg_ParseTuple(args, "ii", &index, &absolute_proportion)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_graphs){
		perror("index too high\n");
		return NULL;
	}

	if(absolute_proportion < 0){
		perror("absolute proportion cannot be negative\n");
		return NULL;
	}

	if(global_graphs[index].num_nodes == global_graphs[index].capacity){
		if(request_more_capacity_graph(&(global_graphs[index])) != 0){
			perror("failed to call request_more_capacity_graph\n");
			return NULL;
		}
	}

	if(create_graph_node(&(global_graphs[index].nodes[global_graphs[index].num_nodes]), global_graphs[index].num_dimensions, FP32) != 0){
		perror("failed to call create_graph_node\n");
		return NULL;
	}
	global_graphs[index].nodes[global_graphs[index].num_nodes].absolute_proportion = (uint32_t) absolute_proportion;
	global_graphs[index].num_nodes++;

	PyObject* res = Py_BuildValue("i", 0);

	return res;
}

static PyObject* interface_compute_relative_proportions(PyObject* self, PyObject* args){
	int32_t index;

	if(!PyArg_ParseTuple(args, "i", &index)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_graphs){
		perror("index too high\n");
		return NULL;
	}

	compute_graph_relative_proportions(&(global_graphs[index]));

	PyObject* res = Py_BuildValue("i", 0);

	return res;
}

static PyObject* interface_individual_measure(PyObject* self, PyObject* args){
	int32_t index;
	int32_t id_function;
	double alpha = 1.0, beta = 1.0;

	if(!PyArg_ParseTuple(args, "ii|dd", &index, &id_function, &alpha, &beta)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_graphs){
		perror("index too high\n");
		return NULL;
	}

	double res1, res2;
	switch(id_function){
		case ID_ENTROPY_SHANNON_WEAVER:
			shannon_weaver_entropy_from_graph(&(global_graphs[index]), &res1, &res2);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_Q_LOGARITHMIC:
			q_logarithmic_entropy_from_graph(&(global_graphs[index]), &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_PATIL_TAILLIE:
			patil_taillie_entropy_from_graph(&(global_graphs[index]), &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_RENYI:
			renyi_entropy_from_graph(&(global_graphs[index]), &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_GOOD:
			good_entropy_from_graph(&(global_graphs[index]), &res1, alpha, beta);
			return Py_BuildValue("(d)", res1);
		default:
			perror("unknown diversity function\n");
			return NULL;
	}
}

static PyObject* interface_cfg_get_value(PyObject* self, PyObject* args){
	int32_t index;
	char* key;

	if(!PyArg_ParseTuple(args, "is", &index, &key)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_graphs){
		perror("index too high\n");
		return NULL;
	}

	PyObject* res;

	char* value = cfg_get_value(&(configurations[index]), key);
	if(value == NULL){
		res = Py_BuildValue("s", "UNKNOWN_KEY");
	} else {
		res = Py_BuildValue("s", value);
	}
	return res;
}

static PyMethodDef diversutilsmethods[] = {
	{"create_graph", interface_create_graph, METH_VARARGS, "Create a graph. This returns the graph index. ARGS: num_nodes, num_dimensions, config_path."},
	{"create_empty_graph", interface_create_empty_graph, METH_VARARGS, "Create a an graph. This returns the graph index. ARGS: num_nodes, num_dimensions."},
	{"free_graph", interface_free_graph, METH_VARARGS, "Free a graph. This requires the graph index."},
	{"cfg_get_value", interface_cfg_get_value, METH_VARARGS, "Provide a graph index, and a key to fetch."},
	{"measurement_from_cfg", interface_measurement_from_cfg, METH_VARARGS, "Call measurement function."},
	{"add_node", interface_add_node, METH_VARARGS, "Add a node (args: graph index, number of dimensions, absolute proportion)."},
	{"individual_measure", interface_individual_measure, METH_VARARGS, "Compute an individual measure. ARGS: graph index, measure index, order 1 (optional), order 2 (optional)."},
	{"compute_relative_proportion", interface_compute_relative_proportions, METH_VARARGS, "Compute relative proportions. ARGS: graph index."},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef diversutilsmodule = {
	PyModuleDef_HEAD_INIT,
	"diversutils",
	NULL, // docs
	-1,
	diversutilsmethods,
	NULL, // PyModuleDef_Slot* m_slots, but NULL if single phase initialisation
	NULL, // traverseproc m_traverse, "A traversal function to call during GC traversal of the module object, or NULL if not needed"
	NULL, // inquiry m_clear, "A clear function to call during GC clearing of the module object, or NULL if not needed"
	NULL, // freefunc m_free, "A function to call during deallocation of the module object, or NULL if not needed"
};

PyMODINIT_FUNC PyInit_diversutils(void){
	PyObject* mod = PyModule_Create(&diversutilsmodule);
	PyModule_AddIntConstant(mod, "ENTROPY_SHANNON_WEAVER", ID_ENTROPY_SHANNON_WEAVER);
	PyModule_AddIntConstant(mod, "ENTROPY_Q_LOGARITHMIC", ID_ENTROPY_Q_LOGARITHMIC);
	PyModule_AddIntConstant(mod, "ENTROPY_PATIL_TAILLIE", ID_ENTROPY_PATIL_TAILLIE);
	PyModule_AddIntConstant(mod, "ENTROPY_RENYI", ID_ENTROPY_RENYI);
	PyModule_AddIntConstant(mod, "ENTROPY_GOOD", ID_ENTROPY_GOOD);
	return mod;
}
