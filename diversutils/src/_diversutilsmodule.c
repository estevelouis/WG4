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
#include "dfunctions.h"

#include "cupt/parser.h"
#include "cupt/load.h"
#include "jsonl/parser.h"
#include "jsonl/load.h"

// for Python bundle
#ifndef DEFINE_UNICODE_CONSTANTS
#define DEFINE_UNICODE_CONSTANTS
#endif
#include "unicode/unicode.h"

enum {
	ID_ENTROPY_SHANNON_WEAVER,
	ID_ENTROPY_Q_LOGARITHMIC,
	ID_ENTROPY_PATIL_TAILLIE,
	ID_ENTROPY_RENYI,
	ID_ENTROPY_GOOD,
	ID_INDEX_SIMPSON_DOMINANCE,
	ID_INDEX_SIMPSON,
	ID_INDEX_RICHNESS,
	ID_INDEX_SPECIES_COUNT,
	ID_INDEX_HILL_EVENNESS,
	ID_INDEX_SHANNON_EVENNESS,
	ID_INDEX_BERGER_PARKER,
	ID_INDEX_JUNGE1994_PAGE22,
	ID_INDEX_BRILLOUIN,
	ID_INDEX_MCINTOSH,
	ID_INDEX_E_HEIP,
	ID_INDEX_ONE_MINUS_D,
	ID_INDEX_ONE_OVER_D_WILLIAMS1964,
	ID_INDEX_E_MINUS_LN_D_PIELOU1977,
	ID_INDEX_F_2_1_ALATALO1981,
	ID_INDEX_G_2_1_MOLINARI1989,
	ID_INDEX_O_BULLA1994,
	ID_INDEX_E_BULLA1994,
	ID_INDEX_E_MCI_PIELOU1969,
	ID_INDEX_E_PRIME_CAMARGO1993,
	ID_INDEX_E_VAR_SMITH_AND_WILSON1996,
	ID_DISPARITY_PAIRWISE,
	ID_DISPARITY_CHAO_ET_AL_FUNCTIONAL,
	ID_DISPARITY_LEINSTER_COBBOLD,
	ID_DISPARITY_SCHEINER,
	ID_DISPARITY_STIRLING,
	ID_DISPARITY_RICOTTA_SZEIDL,
};

static uint32_t num_graphs = 0;
static struct graph* global_graphs = NULL;
static uint8_t* global_graphs_freed = NULL;
static int32_t* global_graphs_word2vec_bindings = NULL;

static uint32_t num_w2v = 0;
static uint8_t* global_w2v_freed = NULL;
static struct word2vec* global_word2vecs = NULL;

static struct cfg* configurations = NULL;

/* // DO NOT REMOVE
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
*/

static PyObject* interface_create_empty_graph(PyObject* self, PyObject* args){
    (void) self;

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

	alloc_size = (num_graphs + 1) * sizeof(int32_t);
	global_graphs_word2vec_bindings = realloc(global_graphs_word2vec_bindings, alloc_size);
	if(global_graphs_word2vec_bindings == NULL){perror("failed to realloc\n"); goto exit_failure;}
	global_graphs_word2vec_bindings[num_graphs] = -1;

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
    (void) self;

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

	alloc_size = (num_graphs + 1) * sizeof(int32_t);
	global_graphs_word2vec_bindings = realloc(global_graphs_word2vec_bindings, alloc_size);
	if(global_graphs_word2vec_bindings == NULL){perror("failed to realloc\n"); goto exit_failure;}
	global_graphs_word2vec_bindings[num_graphs] = -1;

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
    (void) self;

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

static PyObject* interface_free_w2v(PyObject* self, PyObject* args){
    (void) self;

	int32_t index;

	if(!PyArg_ParseTuple(args, "i", &index)){
		return NULL;
	}

	if(index < 0){
		perror("index must be >= 0\n");
		return NULL;
	}
	if(((uint32_t) index) >= num_w2v){
		perror("index too high\n");
		return NULL;
	}

	if(global_w2v_freed[index] == 0){
		printf("trying to free\n");
		free_word2vec(&(global_word2vecs[index]));
		global_w2v_freed[index] = 1;
	} else {
		printf("w2v already freed\n");
	}

	PyObject* res = Py_BuildValue("i", 0);

	return res;
}

static PyObject* interface_load_w2v(PyObject* self, PyObject* args){
    (void) self;

	char* w2v_path;
	size_t alloc_size;
	PyObject* res;

	if(!PyArg_ParseTuple(args, "s", &w2v_path)){
		return NULL;
	}

	alloc_size = (num_w2v + 1) * sizeof(struct word2vec);
	global_word2vecs = realloc(global_word2vecs, alloc_size);
	if(global_word2vecs == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_word2vecs[num_w2v]), '\0', sizeof(struct word2vec));

	alloc_size = (num_w2v + 1) * sizeof(uint8_t);
	global_w2v_freed = realloc(global_w2v_freed, alloc_size);
	if(global_w2v_freed == NULL){perror("failed to realloc\n"); goto exit_failure;}
	memset((void*) &(global_w2v_freed[num_w2v]), '\0', sizeof(uint8_t));

	if(load_word2vec_binary(&(global_word2vecs[num_w2v]), w2v_path) != 0){
		perror("failed to call load_word2vec_binary\n");
		return NULL;
	}

	res = Py_BuildValue("i", num_w2v);

	num_w2v++;

	return res;

	exit_failure:
	res = Py_BuildValue("i", -1);
	return res;
}

static PyObject* interface_bind_w2v(PyObject* self, PyObject* args){
    (void) self;

	int32_t index_g, index_w2v;

	if(!PyArg_ParseTuple(args, "ii", &index_g, &index_w2v)){
		return NULL;
	}

	if(index_g < 0){perror("index_g must be >= 0\n"); return NULL;}
	if(((uint32_t) index_g) >= num_graphs){perror("index_g too high\n"); return NULL;}
	if(index_w2v < 0){perror("index_w2v must be >= 0\n"); return NULL;}
	if(((uint32_t) index_w2v) >= num_w2v){perror("index_w2v too high\n"); return NULL;}

	global_graphs_word2vec_bindings[index_g] = index_w2v;
	
	printf("connected graph %i and w2v %i\n", index_g, index_w2v);

	return Py_BuildValue("i", 0);
}

static PyObject* interface_add_node(PyObject* self, PyObject* args){
    (void) self;

	int32_t index;
	int32_t absolute_proportion;
	char* key;
	int32_t w2v_index;

	struct word2vec_entry* word2vec_entry_pointer;
	float* vector_pointer;

	key = NULL;

	if(!PyArg_ParseTuple(args, "ii|s", &index, &absolute_proportion, &key)){return NULL;}
	if(index < 0){perror("index must be >= 0\n"); return NULL;}
	if(((uint32_t) index) >= num_graphs){perror("index too high\n"); return NULL;}
	if(absolute_proportion < 0){perror("absolute proportion cannot be negative\n"); return NULL;}

	word2vec_entry_pointer = NULL;
	vector_pointer = NULL;

	if(key != NULL){
		if(global_graphs_word2vec_bindings[index] < 0){perror("No Word2Vec bound to graph\n"); return NULL;}
		w2v_index = word2vec_key_to_index(&(global_word2vecs[global_graphs_word2vec_bindings[index]]), key);
		if(w2v_index == -1){
			fprintf(stdout, "unknown key: %s\n", key);
			return NULL;
		}
		word2vec_entry_pointer = &(global_word2vecs[global_graphs_word2vec_bindings[index]].keys[w2v_index]);
		vector_pointer = global_word2vecs[global_graphs_word2vec_bindings[index]].keys[w2v_index].vector;
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
	global_graphs[index].nodes[global_graphs[index].num_nodes].word2vec_entry_pointer = word2vec_entry_pointer;
	global_graphs[index].nodes[global_graphs[index].num_nodes].vector.fp32 = vector_pointer;
	global_graphs[index].nodes[global_graphs[index].num_nodes].num_dimensions = global_graphs[index].num_dimensions;

	global_graphs[index].num_nodes++;

	PyObject* res = Py_BuildValue("i", 0);

	return res;
}

static PyObject* interface_compute_relative_proportions(PyObject* self, PyObject* args){
    (void) self;

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

static PyObject* individual_measure(struct graph * const g, const int32_t id_function, const double alpha, const double beta){
	double res1, res2;
	switch(id_function){
		case ID_ENTROPY_SHANNON_WEAVER:
			shannon_weaver_entropy_from_graph(g, &res1, &res2);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_Q_LOGARITHMIC:
			q_logarithmic_entropy_from_graph(g, &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_PATIL_TAILLIE:
			patil_taillie_entropy_from_graph(g, &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_RENYI:
			renyi_entropy_from_graph(g, &res1, &res2, alpha);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_ENTROPY_GOOD:
			good_entropy_from_graph(g, &res1, alpha, beta);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SIMPSON_DOMINANCE:
			simpson_dominance_index_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SIMPSON:
			simpson_index_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_RICHNESS:
			richness_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SPECIES_COUNT:
			species_count_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_HILL_EVENNESS:
			hill_evenness_from_graph(g, &res1, alpha, beta);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SHANNON_EVENNESS:
			shannon_evenness_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_BERGER_PARKER:
			berger_parker_index_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_JUNGE1994_PAGE22:
			junge1994_page22_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_BRILLOUIN:
			brillouin_diversity_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_MCINTOSH:
			mcintosh_index_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_HEIP:
			sw_e_heip_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_ONE_MINUS_D:
			sw_e_one_minus_D_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_ONE_OVER_D_WILLIAMS1964:
			sw_e_one_over_D_williams1964_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_MINUS_LN_D_PIELOU1977:
			sw_e_minus_ln_D_pielou1977_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_F_2_1_ALATALO1981:
			sw_f_2_1_alatalo1981_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_G_2_1_MOLINARI1989:
			sw_g_2_1_molinari1989_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_O_BULLA1994:
			sw_o_bulla1994_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_BULLA1994:
			sw_e_bulla1994_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_MCI_PIELOU1969:
			sw_e_mci_pielou1969_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_PRIME_CAMARGO1993:
			sw_e_prime_camargo1993_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_VAR_SMITH_AND_WILSON1996:
			sw_e_var_smith_and_wilson1996_original_from_graph(g, &res1);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_PAIRWISE:
			pairwise_from_graph(g, &res1, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_CHAO_ET_AL_FUNCTIONAL:
			if(chao_et_al_functional_diversity_from_graph(g, &res1, &res2, alpha, FP32, NULL) != 0){
				perror("failed to call chao_et_al_functional_diversity_from_graph\n");
				return NULL;
			}
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_LEINSTER_COBBOLD:
			leinster_cobbold_diversity_from_graph(g, &res1, &res2, alpha, FP32, NULL);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_SCHEINER:
			if(scheiner_species_phylogenetic_functional_diversity_from_graph(g, &res1, &res2, alpha, FP32, NULL) != 0){
				perror("failed to call scheiner_species_phylogenetic_functional_diversity_from_graph\n");
				return NULL;
			}
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_STIRLING:
			stirling_from_graph(g, &res1, alpha, beta, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_RICOTTA_SZEIDL:
			ricotta_szeidl_from_graph(g, &res1, alpha, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		default:
			perror("unknown diversity function\n");
			return NULL;
	}
}

static PyObject* interface_individual_measure(PyObject* self, PyObject* args){
    (void) self;

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
	
    /*
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
		case ID_INDEX_SIMPSON_DOMINANCE:
			simpson_dominance_index_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SIMPSON:
			simpson_index_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_RICHNESS:
			richness_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SPECIES_COUNT:
			species_count_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_HILL_EVENNESS:
			hill_evenness_from_graph(&(global_graphs[index]), &res1, alpha, beta);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_SHANNON_EVENNESS:
			shannon_evenness_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_BERGER_PARKER:
			berger_parker_index_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_JUNGE1994_PAGE22:
			junge1994_page22_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_BRILLOUIN:
			brillouin_diversity_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_MCINTOSH:
			mcintosh_index_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_HEIP:
			sw_e_heip_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_ONE_MINUS_D:
			sw_e_one_minus_D_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_ONE_OVER_D_WILLIAMS1964:
			sw_e_one_over_D_williams1964_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_MINUS_LN_D_PIELOU1977:
			sw_e_minus_ln_D_pielou1977_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_F_2_1_ALATALO1981:
			sw_f_2_1_alatalo1981_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_G_2_1_MOLINARI1989:
			sw_g_2_1_molinari1989_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_O_BULLA1994:
			sw_o_bulla1994_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_BULLA1994:
			sw_e_bulla1994_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_MCI_PIELOU1969:
			sw_e_mci_pielou1969_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_PRIME_CAMARGO1993:
			sw_e_prime_camargo1993_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_INDEX_E_VAR_SMITH_AND_WILSON1996:
			sw_e_var_smith_and_wilson1996_original_from_graph(&(global_graphs[index]), &res1);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_PAIRWISE:
			pairwise_from_graph(&(global_graphs[index]), &res1, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_CHAO_ET_AL_FUNCTIONAL:
			if(chao_et_al_functional_diversity_from_graph(&(global_graphs[index]), &res1, &res2, alpha, FP32, NULL) != 0){
				perror("failed to call chao_et_al_functional_diversity_from_graph\n");
				return NULL;
			}
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_LEINSTER_COBBOLD:
			leinster_cobbold_diversity_from_graph(&(global_graphs[index]), &res1, &res2, alpha, FP32, NULL);
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_SCHEINER:
			if(scheiner_species_phylogenetic_functional_diversity_from_graph(&(global_graphs[index]), &res1, &res2, alpha, FP32, NULL) != 0){
				perror("failed to call scheiner_species_phylogenetic_functional_diversity_from_graph\n");
				return NULL;
			}
			return Py_BuildValue("(dd)", res1, res2);
		case ID_DISPARITY_STIRLING:
			stirling_from_graph(&(global_graphs[index]), &res1, alpha, beta, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		case ID_DISPARITY_RICOTTA_SZEIDL:
			ricotta_szeidl_from_graph(&(global_graphs[index]), &res1, alpha, FP32, NULL);
			return Py_BuildValue("(d)", res1);
		default:
			perror("unknown diversity function\n");
			return NULL;
	}
    */

    return individual_measure(&(global_graphs[index]), id_function, alpha, beta);
}

/* // DO NOT REMOVE
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
*/

static PyObject* interface_score_file(PyObject* self, PyObject* args){
    (void) self;

    /*
    PyListObject * listFiles;
    PyListObject * listFunctions;
    PyListObject * listResults;
    */
    PyObject * listFiles;
    PyObject * listFunctions;
    PyObject * listResults;
    int32_t w2v_index = -1;
    int32_t cardinality_files = 0;
    int32_t cardinality_functions = 0;

    // if(!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &listFiles, &PyList_Type, &listFunctions)){ // works
    if(!PyArg_ParseTuple(args, "O!O!i", &PyList_Type, &listFiles, &PyList_Type, &listFunctions, &w2v_index)){
        fprintf(stderr, "Failed to parse arguments!\n");
        return NULL;
    }

    if(PyList_Check(listFiles) != 1 || PyList_Check(listFunctions) != 1){
        fprintf(stderr, "Both arguments must be lists.\n");
        return NULL;
    }

    cardinality_files = PyList_Size(listFiles);
    cardinality_functions = PyList_Size(listFunctions);

    reset_word2vec_active_in_current_graph(&global_word2vecs[w2v_index]);

    struct graph g = {0};

    if(create_graph_empty(&g) != 0){
        fprintf(stderr, "Failed to call create_graph_empty.\n");
        return NULL;
    }

    struct sorted_array sorted_array_discarded_because_not_in_vector_database = {0};
    if(create_sorted_array(&sorted_array_discarded_because_not_in_vector_database, 0, sizeof(struct sorted_array_str_int_element), sorted_array_str_int_cmp) != 0){
        fprintf(stderr, "Failed to call create_sorted_array.\n");
        return NULL;
    }

    for(int32_t i = 0 ; i < cardinality_files ; i++){
        char * s;
        // if(!PyArg_ParseTuple(PyList_GetItem(listFiles, i), "s", &s)){
        if(!PyArg_Parse(PyList_GetItem(listFiles, i), "s", &s)){
            fprintf(stderr, "Failed to get file name at index %i.\n", i);
            return NULL;
        }
        printf("Processing file %i/%i: %s\n", i, cardinality_files, s);

        struct measurement_configuration mcfg = {
            .target_column = UD_FORM,
            .enable_token_utf8_normalisation = 0,
            .jsonl_content_key = "text",
            .io = (struct measurement_io) {
                .input_path = s,
                .jsonl_content_key = "text",
            },
		.threading = (struct measurement_threading) {0},
		.steps = (struct measurement_step_parameters) { .sentence = (struct measurement_step) {0}, .document = (struct measurement_step) {0} },
        };

        struct measurement_structure_references sref = {
            .g = &g,
            .sorted_array_discarded_because_not_in_vector_database = &sorted_array_discarded_because_not_in_vector_database,
            .w2v = &(global_word2vecs[w2v_index]),
        };

        // struct measurement_mutables mmut = { .best_s = 0.0, .prev_best_s = 0.0, .prev_num_nodes = 0, .sentence = (struct measurement_mutable_counters) {0}, .document = (struct measurement_mutable_counters) {0}, .mst_initialised = 0, .mutex = (pthread_mutex_t) {0}, };
        struct measurement_mutables mmut = { .best_s = 0.0, .prev_best_s = 0.0, .prev_num_nodes = 0, .sentence = (struct measurement_mutable_counters) {0}, .document = (struct measurement_mutable_counters) {0}, .mst_initialised = 0, };
        if(pthread_mutex_init(&mmut.mutex, NULL) != 0){
            fprintf(stderr, "Failed to call pthread_mutex_init.\n");
            free_graph(&g);
            return NULL;
        }

        size_t len_s = strlen(s);

        if(strcmp(s + len_s - 5, ".cupt") == 0 || strcmp(s + len_s - 7, ".conllu") == 0){
            if(cupt_to_graph(i, s, &mcfg, &sref, &mmut) != 0){
                fprintf(stderr, "Failed to call cupt_to_graph for %s.\n", s);
                pthread_mutex_destroy(&mmut.mutex);
                free_graph(&g);
                return NULL;
            }
        } else if(strcmp(s + len_s - 6, ".jsonl") == 0){
            if(jsonl_to_graph(i, s, &mcfg, &sref, &mmut) != 0){
                fprintf(stderr, "Failed to call jsonl_to_graph for %s.\n", s);
                pthread_mutex_destroy(&mmut.mutex);
                free_graph(&g);
                return NULL;
            }
        }

        pthread_mutex_destroy(&mmut.mutex);
    }

    listResults = PyList_New(cardinality_functions);
    if(listResults == NULL){
        fprintf(stderr, "Failed to create a new list.\n");
        free_graph(&g);
        return NULL;
    }

    for(int32_t j = 0 ; j < cardinality_functions ; j++){
        PyObject * id_function_obj = PyList_GetItem(listFunctions, j);
        int32_t id_function = -1;
        if(!PyArg_Parse(id_function_obj, "i", &id_function)){
            fprintf(stderr, "Failed to transform Python object to int32_t.\n");
            return NULL;
        }
        PyObject * diversity_score = individual_measure(&g, id_function, 1.0, 1.0);
        if(PyList_SetItem(listResults, j, diversity_score) != 0){
            fprintf(stderr, "Failed to set diversity score at index %i.\n", j);
            return NULL;
        }
    }

    free_sorted_array(&sorted_array_discarded_because_not_in_vector_database);
    free_graph(&g);

    return listResults;
}

void interface_free_globals(void * args){
    (void) args;
	if(global_graphs != NULL){
        for(uint32_t i = 0 ; i < num_graphs ; i++){
            if(!global_graphs_freed[i]){
                free_graph(&global_graphs[i]);
            }
        }
        free(global_graphs);
    }
	if(global_graphs_freed != NULL){free(global_graphs_freed);}
	if(global_graphs_word2vec_bindings != NULL){free(global_graphs_word2vec_bindings);}
	if(global_word2vecs != NULL){
        for(uint32_t i = 0 ; i < num_w2v ; i++){
            if(!global_w2v_freed[i]){
                free_word2vec(&global_word2vecs[i]);
            }
        }
        free(global_word2vecs);
    }
	if(global_w2v_freed != NULL){free(global_w2v_freed);}
	if(configurations != NULL){free(configurations);}

	// return Py_BuildValue("i", 0);
}

static PyMethodDef diversutilsmethods[] = {
	// {"__del__", interface_free_globals, METH_VARARGS, "__del__"},
	{"create_graph", interface_create_graph, METH_VARARGS, "Create a graph. This returns the graph index. ARGS: num_nodes, num_dimensions, config_path."},
	{"create_empty_graph", interface_create_empty_graph, METH_VARARGS, "Create a an graph. This returns the graph index. ARGS: num_nodes, num_dimensions."},
	{"free_graph", interface_free_graph, METH_VARARGS, "Free a graph. This requires the graph index."},
	// {"cfg_get_value", interface_cfg_get_value, METH_VARARGS, "Provide a graph index, and a key to fetch."}, // DO NOT REMOVE
	// {"measurement_from_cfg", interface_measurement_from_cfg, METH_VARARGS, "Call measurement function."}, // DO NOT REMOVE
	{"add_node", interface_add_node, METH_VARARGS, "Add a node (args: graph index, number of dimensions, absolute proportion)."},
	{"individual_measure", interface_individual_measure, METH_VARARGS, "Compute an individual measure. ARGS: graph index, measure index, order 1 (optional), order 2 (optional)."},
	{"compute_relative_proportion", interface_compute_relative_proportions, METH_VARARGS, "Compute relative proportions. ARGS: graph index."},
	{"bind_w2v", interface_bind_w2v, METH_VARARGS, "Bind a Word2Vec binary to a graph. ARGS: graph index, Word2Vec index."},
	{"load_w2v", interface_load_w2v, METH_VARARGS, "Load a Word2Vec binary. ARGS: Word2Vec index."},
	{"free_w2v", interface_free_w2v, METH_VARARGS, "Free a Word2Vec binary. ARGS: Word2Vec index."},
    {"score_file", interface_score_file, METH_VARARGS, "Compute scores for one or more files. ARGS: List of files, list of measures."},
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
	// NULL, // freefunc m_free, "A function to call during deallocation of the module object, or NULL if not needed"
	interface_free_globals, // freefunc m_free, "A function to call during deallocation of the module object, or NULL if not needed"
};

PyMODINIT_FUNC PyInit__diversutils(void){
	PyObject* mod = PyModule_Create(&diversutilsmodule);
	PyModule_AddIntConstant(mod, "DF_ENTROPY_SHANNON_WEAVER", ID_ENTROPY_SHANNON_WEAVER);
	PyModule_AddIntConstant(mod, "DF_ENTROPY_Q_LOGARITHMIC", ID_ENTROPY_Q_LOGARITHMIC);
	PyModule_AddIntConstant(mod, "DF_ENTROPY_PATIL_TAILLIE", ID_ENTROPY_PATIL_TAILLIE);
	PyModule_AddIntConstant(mod, "DF_ENTROPY_RENYI", ID_ENTROPY_RENYI);
	PyModule_AddIntConstant(mod, "DF_ENTROPY_GOOD", ID_ENTROPY_GOOD);
	PyModule_AddIntConstant(mod, "DF_INDEX_SIMPSON_DOMINANCE", ID_INDEX_SIMPSON_DOMINANCE);
	PyModule_AddIntConstant(mod, "DF_INDEX_SIMPSON", ID_INDEX_SIMPSON);
	PyModule_AddIntConstant(mod, "DF_INDEX_RICHNESS", ID_INDEX_RICHNESS);
	PyModule_AddIntConstant(mod, "DF_INDEX_SPECIES_COUNT", ID_INDEX_SPECIES_COUNT);
	PyModule_AddIntConstant(mod, "DF_INDEX_HILL_EVENNESS", ID_INDEX_HILL_EVENNESS);
	PyModule_AddIntConstant(mod, "DF_INDEX_SHANNON_EVENNESS", ID_INDEX_SHANNON_EVENNESS);
	PyModule_AddIntConstant(mod, "DF_INDEX_BERGER_PARKER", ID_INDEX_BERGER_PARKER);
	PyModule_AddIntConstant(mod, "DF_INDEX_JUNGE1994_PAGE22", ID_INDEX_JUNGE1994_PAGE22);
	PyModule_AddIntConstant(mod, "DF_INDEX_BRILLOUIN", ID_INDEX_BRILLOUIN);
	PyModule_AddIntConstant(mod, "DF_INDEX_MCINTOSH", ID_INDEX_MCINTOSH);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_HEIP", ID_INDEX_E_HEIP);
	PyModule_AddIntConstant(mod, "DF_INDEX_ONE_MINUS_D", ID_INDEX_ONE_MINUS_D);
	PyModule_AddIntConstant(mod, "DF_INDEX_ONE_OVER_D_WILLIAMS1964", ID_INDEX_ONE_OVER_D_WILLIAMS1964);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_MINUS_LN_D_PIELOU1977", ID_INDEX_E_MINUS_LN_D_PIELOU1977);
	PyModule_AddIntConstant(mod, "DF_INDEX_F_2_1_ALATALO1981", ID_INDEX_F_2_1_ALATALO1981);
	PyModule_AddIntConstant(mod, "DF_INDEX_G_2_1_MOLINARI1989", ID_INDEX_G_2_1_MOLINARI1989);
	PyModule_AddIntConstant(mod, "DF_INDEX_O_BULLA1994", ID_INDEX_O_BULLA1994);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_BULLA1994", ID_INDEX_E_BULLA1994);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_MCI_PIELOU1969", ID_INDEX_E_MCI_PIELOU1969);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_PRIME_CAMARGO1993", ID_INDEX_E_PRIME_CAMARGO1993);
	PyModule_AddIntConstant(mod, "DF_INDEX_E_VAR_SMITH_AND_WILSON1996", ID_INDEX_E_VAR_SMITH_AND_WILSON1996);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_PAIRWISE", ID_DISPARITY_PAIRWISE);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_CHAO_ET_AL_FUNCTIONAL", ID_DISPARITY_CHAO_ET_AL_FUNCTIONAL);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_LEINSTER_COBBOLD", ID_DISPARITY_LEINSTER_COBBOLD);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_SCHEINER", ID_DISPARITY_SCHEINER);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_STIRLING", ID_DISPARITY_STIRLING);
	PyModule_AddIntConstant(mod, "DF_DISPARITY_RICOTTA_SZEIDL", ID_DISPARITY_RICOTTA_SZEIDL);
	return mod;
}
