#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

#include "graph.h"
#include "cfgparser/parser.h"

#include "measurement.h"

static uint32_t num_graphs = 0;
static struct graph* global_graphs = NULL;
static uint8_t* global_graphs_freed = NULL;
static struct cfg* configurations = NULL;

static PyObject* diversutils_test(PyObject* self, PyObject* args){
	const char* s;

	if(!PyArg_ParseTuple(args, "s", &s)){
		return NULL;
	}

	printf("input string: %s\n", s);

	// return 0;
	// PyObject* res = Py_BuildValue("d", 0);
	PyObject* res = Py_BuildValue("i", 0);
	return res;
}

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
	{"test", diversutils_test, METH_VARARGS, "test (description)"},
	{"create_graph", interface_create_graph, METH_VARARGS, "Create a graph. This returns the graph index."},
	{"free_graph", interface_free_graph, METH_VARARGS, "Free a graph. This requires the graph index."},
	{"cfg_get_value", interface_cfg_get_value, METH_VARARGS, "Provide a graph index, and a key to fetch."},
	{"measurement_from_cfg", interface_measurement_from_cfg, METH_VARARGS, "Call measurement function."},
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
	/*
	PyObject* m;

	m = PyModule_Create(&
	*/
	return PyModule_Create(&diversutilsmodule);
}
