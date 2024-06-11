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

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "dfunctions.h"
#include "general_constants.h"
#include "graph.h"

#define LOGARITHMIC_BASE E

// static const double NORMALISATION_BASE = log(LOGARITHMIC_BASE);

#define NORMALISATION_BASE log(LOGARITHMIC_BASE)

double q_logarithm(const double x, const double q) {
  if (q == 1.0) {
    return log(x);
  } else {
    return (pow(x, 1.0 - q) - 1) / (1.0 - q);
  }
}

void shannon_weaver_entropy_from_graph(const struct graph *const g, double *const res_entropy, double *const res_hill_number) {
  double loc_res = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    if (g->nodes[i].relative_proportion <= 0.0) {
      continue;
    }
    loc_res += g->nodes[i].relative_proportion * (log(g->nodes[i].relative_proportion) / log(LOGARITHMIC_BASE));
  }
  loc_res *= -1.0;
  (*res_entropy) = loc_res;
  (*res_hill_number) = pow(LOGARITHMIC_BASE, loc_res);
}

void good_entropy_from_graph(const struct graph *const g, double *const res, double alpha, double beta) {
  double loc_res = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    loc_res +=
        pow(g->nodes[i].relative_proportion, alpha) * pow(-log(g->nodes[i].relative_proportion) / log(LOGARITHMIC_BASE), beta);
  }
  (*res) = loc_res;
}

void renyi_entropy_from_graph(const struct graph *const g, double *const res_entropy, double *const res_hill_number,
                              double alpha) {
  if (alpha == 1.0) {
    shannon_weaver_entropy_from_graph(g, res_entropy, res_hill_number);
  } else {
    double loc_res = 0.0;
    for (uint64_t i = 0; i < g->num_nodes; i++) {
      loc_res += pow(g->nodes[i].relative_proportion, alpha);
    }
    loc_res = (1.0 / (1.0 - alpha)) * (log(loc_res) / log(LOGARITHMIC_BASE));
    (*res_entropy) = loc_res;
    (*res_hill_number) = pow(LOGARITHMIC_BASE, loc_res);
  }
}

void patil_taillie_entropy_from_graph(const struct graph *const g, double *const res_entropy, double *const res_hill_number,
                                      double alpha) {
  if (alpha == 0.0) {
    shannon_weaver_entropy_from_graph(g, res_entropy, res_hill_number);
  } else {
    double loc_res = 1.0;
    for (uint64_t i = 0; i < g->num_nodes; i++) {
      loc_res -= pow(g->nodes[i].relative_proportion, alpha + 1.0);
    }
    loc_res /= alpha;
    (*res_entropy) = loc_res;
    (*res_hill_number) = 1.0 / pow(1.0 - (alpha * loc_res), 1.0 / alpha);
  }
}

void q_logarithmic_entropy_from_graph(const struct graph *const g, double *const res_entropy, double *const res_hill_number,
                                      double q) {
  double loc_res = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    loc_res += g->nodes[i].relative_proportion * q_logarithm(1.0 / g->nodes[i].relative_proportion, q);
  }
  (*res_entropy) = loc_res;
  if (q == 1.0) {
    (*res_hill_number) = pow(LOGARITHMIC_BASE, loc_res); // ?
  } else {
    (*res_hill_number) = pow(1.0 - (q - 1.0) * loc_res, 1.0 / (1.0 - q));
  }
}

void simpson_dominance_index_from_graph(const struct graph *const g, double *const res) {
  double loc_res = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    loc_res += pow(g->nodes[i].relative_proportion, 2.0);
  }
  (*res) = loc_res;
}

void simpson_index_from_graph(const struct graph *const g, double *const res) {
  double loc_res;
  simpson_dominance_index_from_graph(g, &loc_res);
  (*res) = 1.0 - loc_res;
}

void richness_from_graph(const struct graph *const g, double *const res) {
  (*res) = (double)g->num_nodes;
}

void species_count_from_graph(const struct graph *const g, double *const res) {
  (*res) = (double)(g->num_nodes - 1);
}

void hill_number_standard_from_graph(const struct graph *const g, double *const res, double alpha) {
  double renyi_entropy;
  double hill_number;
  renyi_entropy_from_graph(g, &renyi_entropy, &hill_number, alpha);
  (*res) = hill_number;
}

void hill_evenness_from_graph(const struct graph *const g, double *const res, double alpha, double beta) {
  double loc_res_upper;
  double loc_res_lower;
  hill_number_standard_from_graph(g, &loc_res_upper, alpha);
  hill_number_standard_from_graph(g, &loc_res_lower, beta);
  (*res) = loc_res_upper / loc_res_lower;
}

void berger_parker_index_from_graph(const struct graph *const g, double *const res) {
  double loc_res = g->nodes[0].relative_proportion;
  for (uint64_t i = 1; i < g->num_nodes; i++) {
    if (g->nodes[i].relative_proportion > loc_res) {
      loc_res = g->nodes[i].relative_proportion;
    }
  }
  (*res) = loc_res;
}

void shannon_evenness_from_graph(const struct graph *const g, double *const res) {
  double sw_entropy;
  double hill_number;
  shannon_weaver_entropy_from_graph(g, &sw_entropy, &hill_number);
  double max_entropy = log((double)g->num_nodes) / log(LOGARITHMIC_BASE);
  (*res) = sw_entropy / max_entropy;
}

void junge1994_page22_from_graph(const struct graph *const g, double *res) {
  double sum = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum += pow(g->nodes[i].relative_proportion, 2.0);
  }
  (*res) = 1.0 - pow(sum, 0.5);
}

void brillouin_diversity_from_graph(const struct graph *const g, double *const res) {
  /*
  double sum_left = 0.0;
  double sum_right = 0.0;
  for(uint64_t i = 0 ; i < g->num_nodes ; i++){
          sum_left += log((double) (i+1)) / log(LOGARITHMIC_BASE);
          int64_t fact = 1;
          for(uint64_t j = 1 ; j <= g->nodes[i].absolute_proportion ; j++){
                  fact *= (int64_t) j;
          }
          printf("fact(%i): %li\n", g->nodes[i].absolute_proportion, fact);
          sum_right += log((double) fact) / log(LOGARITHMIC_BASE);
  }
  (*res) = sum_left - sum_right;
  */

  double sum_left = 0.0;
  double sum_right = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum_left += log((double)(i + 1)) / log(LOGARITHMIC_BASE);
    for (uint64_t j = 0; j < (uint64_t)g->nodes[i].absolute_proportion; j++) {
      sum_right += log((double)(j + 1)) / log(LOGARITHMIC_BASE);
    }
  }
  (*res) = sum_left - sum_right;
}

void mcintosh_index_from_graph(const struct graph *const g, double *const res) {
  double sum = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum += pow(g->nodes[i].relative_proportion, 2.0);
  }
  (*res) = 1.0 - pow(sum, 0.5);
}

void sw_entropy_over_log_n_species_pielou1975_from_graph(const struct graph *const g, double *const res) {
  // should it necessarily be logarithmic of base e?
  shannon_evenness_from_graph(g, res);
}

void sw_e_heip_from_graph(const struct graph *const g, double *const res) {
  double sw_entropy;
  double hill_number;
  shannon_weaver_entropy_from_graph(g, &sw_entropy, &hill_number);
  (*res) = (pow(E, sw_entropy) - 1.0) / ((double)(g->num_nodes - 1));
}

void sw_e_one_minus_D_from_graph(const struct graph *const g, double *const res) {
  double dom_index;
  simpson_dominance_index_from_graph(g, &dom_index);
  (*res) = (1.0 - dom_index) / (1.0 - (1.0 / ((double)g->num_nodes)));
}

void sw_e_one_over_D_williams1964_from_graph(const struct graph *const g, double *const res) {
  double dom_index;
  simpson_dominance_index_from_graph(g, &dom_index);
  (*res) = (1.0 / dom_index) / ((double)g->num_nodes);
}

void sw_e_minus_ln_D_pielou1977_from_graph(const struct graph *const g, double *const res) {
  double dom_index;
  simpson_dominance_index_from_graph(g, &dom_index);
  (*res) = (-(log(dom_index) / log(LOGARITHMIC_BASE))) / (log((double)g->num_nodes) / log(LOGARITHMIC_BASE));
}

void sw_f_2_1_alatalo1981_from_graph(const struct graph *const g, double *const res) {
  double dom_index;
  simpson_dominance_index_from_graph(g, &dom_index);
  double sw_entropy;
  double hill_number;
  shannon_weaver_entropy_from_graph(g, &sw_entropy, &hill_number);
  (*res) = ((1.0 / dom_index) - 1.0) / (pow(LOGARITHMIC_BASE, sw_entropy) - 1.0);
}

void sw_g_2_1_molinari1989_from_graph(const struct graph *const g, double *const res) {
  double f_2_1;
  sw_f_2_1_alatalo1981_from_graph(g, &f_2_1);
  if (f_2_1 > pow(0.5, 0.5)) {
    (*res) = f_2_1 * 0.636611 * asin(f_2_1);
  } else {
    (*res) = pow(f_2_1, 3.0);
  }
}

void sw_o_bulla1994_from_graph(const struct graph *const g, double *const res) {
  double sum = 0.0;
  double division = 1.0 / ((double)g->num_nodes);
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    if (g->nodes[i].relative_proportion < division) {
      sum += g->nodes[i].relative_proportion;
    } else {
      sum += division;
    }
  }
  (*res) = sum;
}

void sw_e_bulla1994_from_graph(const struct graph *const g, double *const res) {
  double o_bulla;
  sw_o_bulla1994_from_graph(g, &o_bulla);
  (*res) = (o_bulla - (1.0 / ((double)g->num_nodes))) / (1.0 - (1.0 / ((double)g->num_nodes)));
}

void sw_e_mci_pielou1969_from_graph(const struct graph *const g, double *const res) {
  double sum_x = 0.0;
  double sum_x_square = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    sum_x += g->nodes[i].relative_proportion;
    sum_x_square += pow(g->nodes[i].relative_proportion, 2.0);
  }
  (*res) = (sum_x - pow(sum_x_square, 0.5)) / (sum_x - (sum_x / pow((double)g->num_nodes, 0.5)));
}

void sw_e_prime_camargo1993_from_graph(const struct graph *const g, double *const res) {
  double sum = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    for (uint64_t j = i + 1; j < g->num_nodes; j++) {
      double val = g->nodes[i].relative_proportion - g->nodes[j].relative_proportion;
      if (val < 0.0) {
        val *= -1.0;
      }
      sum += val / ((double)g->num_nodes);
    }
  }
  (*res) = 1.0 - sum;
}

void sw_e_var_smith_and_wilson1996_original_from_graph(const struct graph *const g, double *const res) {
  double inner_sum = 0.0;
  double outer_sum = 0.0;
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    inner_sum += (log(g->nodes[i].relative_proportion) / log(LOGARITHMIC_BASE)) / ((double)g->num_nodes);
  }
  for (uint64_t i = 0; i < g->num_nodes; i++) {
    outer_sum += pow((log(g->nodes[i].relative_proportion) / log(LOGARITHMIC_BASE)) - inner_sum, 2.0) / ((double)g->num_nodes);
  }
  (*res) = 1.0 - ((2.0 / PI) * atan(outer_sum));
}

/* ======== MULTITHREAD ======== */

void *non_disparity_thread(void *args) {
  double (*function_transform_proportion)(const double, const double, const double) =
      ((struct non_disparity_multithread_args *)args)->function_transform_proportion;
  double (*function_agregate_local)(const double, const double) =
      ((struct non_disparity_multithread_args *)args)->function_agregate_local;
  // double (* function_finalise)(double, double, double),
  // double thread_result = ((struct non_disparity_multithread_args *) args)->thread_result;
  struct graph *g = ((struct non_disparity_multithread_args *)args)->g;
  const int32_t start_index = ((struct non_disparity_multithread_args *)args)->start_index;
  const int32_t end_index = ((struct non_disparity_multithread_args *)args)->end_index;
  double local_result = ((struct non_disparity_multithread_args *)args)->agregation_identity;
  double order0 = ((struct non_disparity_multithread_args *)args)->order0;
  double order1 = ((struct non_disparity_multithread_args *)args)->order1;

  for (int32_t i = start_index; i < end_index; i++) {
    double local_transformation = function_transform_proportion(g->nodes[i].relative_proportion, order0, order1);
    if (isnan(local_transformation)) {
      continue;
    }
    local_result = function_agregate_local(local_result, local_transformation);
  }

  ((struct non_disparity_multithread_args *)args)->thread_result = local_result;

  return NULL;
}

int32_t non_disparity_multithread(double (*function_transform_proportion)(const double, const double, const double),
                                  double (*function_agregate_local)(const double, const double),
                                  double (*function_finalise)(const double, const double, const double),
                                  double agregation_identity, struct graph *g, double order0, double order1, double *res0,
                                  double *res1, int32_t num_threads) {
  pthread_t threads[num_threads];
  struct non_disparity_multithread_args args[num_threads];

  int32_t start_index = 0;
  int32_t base_step = (int32_t)floor(((double)g->num_nodes) / ((double)num_threads));
  for (int32_t i = 0; i < num_threads; i++) {
    // int32_t num_elements_for_current_thread = base_step + ((g->num_nodes % num_threads) < num_threads);
    int32_t num_elements_for_current_thread = base_step + (i < (int32_t)(g->num_nodes % num_threads));
    int32_t end_index = start_index + num_elements_for_current_thread;
    // printf("%i to %i\n", start_index, end_index);
    args[i] = (struct non_disparity_multithread_args){
        .function_transform_proportion = function_transform_proportion,
        .function_agregate_local = function_agregate_local,
        .thread_result = agregation_identity,
        .g = g,
        .start_index = start_index,
        .end_index = end_index,
        .agregation_identity = agregation_identity,
        .order0 = order0,
        .order1 = order1,
    };
    if (pthread_create(&(threads[i]), NULL, non_disparity_thread, &(args[i])) != 0) {
      perror("Failed to call pthread_create\n");
      return 1;
    }
    start_index = end_index;
  }

  double global_result = agregation_identity;

  for (int32_t i = 0; i < num_threads; i++) {
    if (pthread_join(threads[i], NULL) != 0) {
      perror("Failed to call pthread_join\n");
      return 1;
    }
    // printf("thread %i finished with result %f\n", i, args[i].thread_result);
    global_result = function_agregate_local(global_result, args[i].thread_result);
  }

  if (function_finalise != NULL) {
    global_result = function_finalise(global_result, order0, order1);
  }

  (*res0) = global_result;
  (void)res1; // ...

  return 0;
}

double agregate_local_add(const double a, const double b) {
  return a + b;
}

double agregate_local_multiply(const double a, const double b) {
  return a * b;
}

double entropy_shannon_weaver_transform_proportion(const double p, const double order0, const double order1) {
  (void)order0;
  (void)order1;

  return p * (-(log(p) / NORMALISATION_BASE));
}
double entropy_shannon_weaver_to_hill_number(const double h) {
  return pow(LOGARITHMIC_BASE, h);
}

double entropy_good_transform_proportion(const double p, const double order0, const double order1) {
  return pow(p, order0) * pow(-(log(p) / NORMALISATION_BASE), order1);
}

double entropy_renyi_transform_proportion(const double p, const double order0, const double order1) {
  (void)order1;

  if (order0 == 1.0) {
    return entropy_shannon_weaver_transform_proportion(p, order0, order1);
  } else {
    return pow(p, order0);
  }
}

double entropy_renyi_finalise(const double x, const double order0, const double order1) {
  (void)order1;

  if (order0 == 1.0) {
    return x;
  } else {
    return (log(x) / NORMALISATION_BASE) / (1.0 - order0);
  }
}

double entropy_renyi_to_hill_number(const double h) {
  return pow(LOGARITHMIC_BASE, h);
}
