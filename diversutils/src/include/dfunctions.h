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

#ifndef DFUNCTIONS_H
#define DFUNCTIONS_H

#include "general_constants.h"
#include "graph.h"

double q_logarithm(const double, const double);
void shannon_weaver_entropy_from_graph(const struct graph *const, double *const, double *const);
void good_entropy_from_graph(const struct graph *const, double *const, double, double);
void renyi_entropy_from_graph(const struct graph *const, double *const, double *const, double);
void patil_taillie_entropy_from_graph(const struct graph *const, double *const, double *const, double);
void q_logarithmic_entropy_from_graph(const struct graph *const, double *const, double *const, double);
void simpson_dominance_index_from_graph(const struct graph *const, double *const);
void simpson_index_from_graph(const struct graph *const, double *const);
void richness_from_graph(const struct graph *const, double *const);
void species_count_from_graph(const struct graph *const, double *const);
void hill_number_standard_from_graph(const struct graph *const, double *const, double);
void hill_evenness_from_graph(const struct graph *const, double *const, double, double);
void berger_parker_index_from_graph(const struct graph *const, double *const);
void shannon_evenness_from_graph(const struct graph *const, double *const);
void junge1994_page22_from_graph(const struct graph *const, double *const);
void brillouin_diversity_from_graph(const struct graph *const, double *const);
void mcintosh_index_from_graph(const struct graph *const, double *const);
void sw_entropy_over_log_n_species_pielou1975_from_graph(const struct graph *const, double *const);
void sw_e_heip_from_graph(const struct graph *const, double *const);
void sw_e_one_minus_D_from_graph(const struct graph *const, double *const);
void sw_e_one_over_D_williams1964_from_graph(const struct graph *const, double *const);
void sw_e_minus_ln_D_pielou1977_from_graph(const struct graph *const, double *const);
void sw_f_2_1_alatalo1981_from_graph(const struct graph *const, double *const);
void sw_g_2_1_molinari1989_from_graph(const struct graph *const, double *const);
void sw_o_bulla1994_from_graph(const struct graph *const, double *const);
void sw_e_bulla1994_from_graph(const struct graph *const, double *const);
void sw_e_mci_pielou1969_from_graph(const struct graph *const, double *const);
void sw_e_var_smith_and_wilson1996_original_from_graph(const struct graph *const, double *const);

void sw_e_prime_camargo1993_from_graph(const struct graph *const, double *const);
struct sw_e_prime_camargo1993_thread_args {
  double sum_local;
  int16_t thread_number;
  int16_t num_threads;
  const struct graph *g;
};
void *sw_e_prime_camargo1993_thread(void *const args);
int32_t sw_e_prime_camargo1993_from_graph_multithread(const struct graph *g, double *const res, const int16_t num_threads);

/* ======== MULTITHREAD ======== */

struct non_disparity_multithread_args {
  double (*function_transform_proportion)(const double, const double, const double);
  double (*function_agregate_local)(const double, const double);
  double thread_result;
  struct graph *g;
  int32_t start_index;
  int32_t end_index;
  double agregation_identity;
  double order0;
  double order1;
};

int32_t non_disparity_multithread(double (*function_transform_proportion)(const double, const double, const double),
                                  double (*function_agregate_local)(const double, const double),
                                  double (*function_finalise)(const double, const double, const double),
                                  double agregation_identity, struct graph *g, double order0, double order1, double *res0,
                                  double *res1, int32_t num_threads);
double agregate_local_add(const double a, const double b);
double agregate_local_multiply(const double a, const double b);
double entropy_shannon_weaver_transform_proportion(const double p, const double order0, const double order1);
double entropy_good_transform_proportion(const double p, const double order0, const double order1);
double entropy_renyi_transform_proportion(const double p, const double order0, const double order1);
double entropy_renyi_finalise(const double p, const double order0, const double order1);
double entropy_shannon_weaver_to_hill_number(const double x);
double entropy_renyi_to_hill_number(const double x);

#endif
