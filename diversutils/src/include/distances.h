#ifndef DISTANCES_H
#define DISTANCES_H

#include <stdint.h>

#include "graph.h" // matrix

double minkowski_distance(double *restrict a, double *restrict b, int n, double order);
float minkowski_distance_fp32(float *restrict a, float *restrict b, int n, float order);
float cosine_distance_fp32(const float *restrict const a, const float *restrict const b, const int32_t n);
float cosine_distance_norm_fp32(float *restrict a, float *restrict b, int32_t n);
#if ENABLE_AVX256 == 1
float cosine_distance_fp32_avx(const float *restrict const a, const float *restrict const b, const int32_t n);
#endif
double cosine_distance(const double *restrict const a, const double *restrict const b, const int n);
double cosine_distance_norm(double *restrict a, double *restrict b, int32_t n);
double chebyshev_distance(double *restrict a, double *restrict b, int n);
double canberra_distance(double *restrict a, double *restrict b, int n);
double bray_curtis_distance(double *restrict a, double *restrict b, int n);
double angular_minkowski_distance(double *restrict a, double *restrict b, int32_t n, double order);
double rootless_angular_minkowski_distance(double *restrict a, double *restrict b, int32_t n, double order);
float angular_minkowski_distance_fp32(float *restrict a, float *restrict b, int32_t n, float order);
float rootless_angular_minkowski_distance_fp32(float *restrict a, float *restrict b, int32_t n, float order);

#endif
