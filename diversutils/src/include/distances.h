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

#ifndef DISTANCES_H
#define DISTANCES_H

#include <stdint.h>
#if (ENABLE_AVX256 == 1 || ENABLE_AVX512 == 1)
#include <immintrin.h>
#endif

#include "graph.h" // matrix

double minkowski_distance(double *restrict a, double *restrict b, int n, double order);
float minkowski_distance_fp32(float *restrict a, float *restrict b, int n, float order);
float cosine_distance_fp32(const float *restrict const a, const float *restrict const b, const int32_t n);
float cosine_distance_norm_fp32(float *restrict a, float *restrict b, int32_t n);
#if ENABLE_AVX256 == 1
float cosine_distance_fp32_avx256(const float *restrict const a, const float *restrict const b, const int32_t n);
#endif
#if ENABLE_AVX512 == 1
float cosine_distance_fp32_avx512(const float *restrict const a, const float *restrict const b, const int32_t n);
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
