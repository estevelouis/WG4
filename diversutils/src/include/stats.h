#ifndef STATS_H
#define STATS_H

#include <stdint.h>

float avg_fp32(const float* p, const int32_t n);
float std_fp32(const float* p, const int32_t n);
void avg_and_std_fp32(const float* p, const int32_t n, float* avg_p, float* std_p);
double avg_fp64(const double* p, const int32_t n);
double std_fp64(const double* p, const int32_t n);
void avg_and_std_fp64(const double* p, const int32_t n, double* avg_p, double* std_p);

#endif

