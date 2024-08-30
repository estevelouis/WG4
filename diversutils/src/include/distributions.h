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

#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>

#include "graph.h"

struct zipfian_distribution {
	union {
		float* fp32;
		double* fp64;
	} vector;
	double s;
	uint32_t n;
};

int32_t create_zipfian_distribution(struct zipfian_distribution* z, double s, uint32_t n, const int8_t fp_mode);
void free_zipfian_distribution(struct zipfian_distribution* z, const int8_t fp_mode);
double mean_squared_error(double* v, double* w, uint32_t n);
int32_t double_cmp_reverse(const void* a, const void* b);
int32_t zipfian_fit(double* v, uint32_t n, double* result);
int32_t zipfian_fit_from_graph(struct graph* g, double* result);

#endif
