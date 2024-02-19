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

#ifndef STATS_H
#define STATS_H

#include<stdint.h>
#include<math.h>

float avg_fp32(const float* p, const int32_t n){
	float sum = 0.0f;
	for(int32_t i = 0 ; i < n ; i++){
		sum += p[i];
	}
	sum /= (float) n;
	return sum;
}

float std_fp32(const float* p, const int32_t n){
	float avg = avg_fp32(p, n);
	float sum = 0.0f;
	for(int32_t i = 0 ; i < n ; i++){
		sum += powf(p[i] - avg, 2.0f);
	}
	sum /= (float) n;
	float result = powf(sum, 0.5f);
	return result;
}

void avg_and_std_fp32(const float* p, const int32_t n, float* avg_p, float* std_p){
	float avg = avg_fp32(p, n);
	float sum = 0.0f;
	for(int32_t i = 0 ; i < n ; i++){
		sum += powf(p[i] - avg, 2.0f);
	}
	sum /= (float) n;
	float result = powf(sum, 0.5f);
	(*avg_p) = avg;
	(*std_p) = result;
}

double avg_fp64(const double* p, const int32_t n){
	double sum = 0.0;
	for(int32_t i = 0 ; i < n ; i++){
		sum += p[i];
	}
	sum /= (double) n;
	return sum;
}

double std_fp64(const double* p, const int32_t n){
	double avg = avg_fp64(p, n);
	double sum = 0.0;
	for(int32_t i = 0 ; i < n ; i++){
		sum += pow(p[i] - avg, 2.0);
	}
	sum /= (double) n;
	double result = pow(sum, 0.5);
	return result;
}

void avg_and_std_fp64(const double* p, const int32_t n, double* avg_p, double* std_p){
	double avg = avg_fp64(p, n);
	double sum = 0.0;
	for(int32_t i = 0 ; i < n ; i++){
		sum += pow(p[i] - avg, 2.0);
	}
	sum /= (double) n;
	double result = pow(sum, 0.5);
	(*avg_p) = avg;
	(*std_p) = result;
}


#endif
