#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>

#include "graph.h"

struct zipfian_distribution {
  union {
    float *fp32;
    double *fp64;
  } vector;
  double s;
  uint32_t n;
};

int32_t create_zipfian_distribution(struct zipfian_distribution *z, double s, uint32_t n, const int8_t fp_mode);
void free_zipfian_distribution(struct zipfian_distribution *z, const int8_t fp_mode);
double mean_squared_error(double *v, double *w, uint32_t n);
int32_t double_cmp_reverse(const void *a, const void *b);
int32_t zipfian_fit(double *v, uint32_t n, double *result);
int32_t zipfian_fit_from_graph(struct graph *g, double *result);

#endif
