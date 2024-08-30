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
#if (ENABLE_AVX256 == 1 || ENABLE_AVX512 == 1)
#include <immintrin.h>
#endif

#include "distances.h"

// double minkowski_distance(double* restrict a, double* restrict b, int n, double order){
double minkowski_distance(double* restrict a, double* restrict b, int n, double order){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double sum_ = 0.0;
	for(int i = 0 ; i < n ; i++){
		double diff = a[i] - b[i];
		if(diff < 0.0){
			diff *= -1.0;
		}
		sum_ += pow(diff, order);
	}
	return pow(sum_, 1.0 / order);
}

float minkowski_distance_fp32(float* restrict a, float* restrict b, int n, float order){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	// printf("n: %i\n", n);
	float sum_ = 0.0f;
	for(int i = 0 ; i < n ; i++){
		float diff = a[i] - b[i];
		if(diff < 0.0f){
			diff *= -1.0f;
		}
		sum_ += (float) pow(diff, order);
	}
	return (float) pow(sum_, 1.0f / order);
}

float cosine_distance_fp32(const float* restrict const a, const float* restrict const b, int32_t n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	// printf("n: %i\n", n);
	float upper_sum = 0.0;
	float lower_sum_a = 0.0;
	float lower_sum_b = 0.0;
	for(int32_t i = 0 ; i < n ; i++){
		upper_sum += a[i] * b[i];
		lower_sum_a += powf(a[i], 2.0f);
		lower_sum_b += powf(b[i], 2.0f);
	}
	float cosine_similarity = (upper_sum / (sqrtf(lower_sum_a) * sqrtf(lower_sum_b)));
	return (1.0f - cosine_similarity);
}

float cosine_distance_norm_fp32(float* restrict a, float* restrict b, int32_t n){
	return cosine_distance_fp32(a, b, n) / 2.0f;
}

#if ENABLE_AVX256 == 1
float cosine_distance_fp32_avx256(const float* restrict const a, const float* restrict const b, int32_t n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	float upper_sum = 0.0;
	float lower_sum_a = 0.0;
	float lower_sum_b = 0.0;

	__m256 avx256_upper_sum = _mm256_setzero_ps();
	__m256 avx256_lower_sum_a = _mm256_setzero_ps();
	__m256 avx256_lower_sum_b = _mm256_setzero_ps();

	float vec_upper_sum[8];
	float vec_lower_sum_a[8];
	float vec_lower_sum_b[8];

	int32_t i = 0;
	while(i < n){
		__m256 avx256_a = _mm256_loadu_ps(&(a[i])); // had to load unaligned
		__m256 avx256_b = _mm256_loadu_ps(&(b[i]));

		__m256 avx256_local_upper_sum = _mm256_mul_ps(avx256_a, avx256_b);
		avx256_upper_sum = _mm256_add_ps(avx256_upper_sum, avx256_local_upper_sum);

		// ----

		__m256 avx256_local_lower_sum_a = _mm256_mul_ps(avx256_a, avx256_a);
		avx256_lower_sum_a = _mm256_add_ps(avx256_lower_sum_a, avx256_local_lower_sum_a);

		__m256 avx256_local_lower_sum_b = _mm256_mul_ps(avx256_b, avx256_b);
		avx256_lower_sum_b = _mm256_add_ps(avx256_lower_sum_b, avx256_local_lower_sum_b);

		i += 8;
	}

	_mm256_storeu_ps(vec_upper_sum, avx256_upper_sum);
	_mm256_storeu_ps(vec_lower_sum_a, avx256_lower_sum_a);
	_mm256_storeu_ps(vec_lower_sum_b, avx256_lower_sum_b);

	for(int32_t j = 0 ; j < 8 ; j++){
		upper_sum += vec_upper_sum[j];
		lower_sum_a += vec_lower_sum_a[j];
		lower_sum_b += vec_lower_sum_b[j];
	}
	while(i < n){
		upper_sum += a[i] * b[i];
		lower_sum_a += powf(a[i], 2.0f);
		lower_sum_b += powf(b[i], 2.0f);
		i++;
	}
	// printf("upper_sum: %f; lower_sum_a: %f; lower_sum_b: %f\n", upper_sum, lower_sum_a, lower_sum_b);
	float cosine_similarity = (upper_sum / (sqrtf(lower_sum_a) * sqrtf(lower_sum_b)));
	return (1.0f - cosine_similarity);
}
#endif

#if ENABLE_AVX512 == 1
float cosine_distance_fp32_avx512(const float* restrict const a, const float* restrict const b, int32_t n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	float upper_sum = 0.0;
	float lower_sum_a = 0.0;
	float lower_sum_b = 0.0;

	__m512 avx512_upper_sum = _mm512_setzero_ps();
	__m512 avx512_lower_sum_a = _mm512_setzero_ps();
	__m512 avx512_lower_sum_b = _mm512_setzero_ps();

	float vec_upper_sum[16];
	float vec_lower_sum_a[16];
	float vec_lower_sum_b[16];

	int32_t i = 0;
	while(i < n){
		__m512 avx512_a = _mm512_loadu_ps(&(a[i])); // had to load unaligned
		__m512 avx512_b = _mm512_loadu_ps(&(b[i]));

		__m512 avx512_local_upper_sum = _mm512_mul_ps(avx512_a, avx512_b);
		avx512_upper_sum = _mm512_add_ps(avx512_upper_sum, avx512_local_upper_sum);

		// ----

		__m512 avx512_local_lower_sum_a = _mm512_mul_ps(avx512_a, avx512_a);
		avx512_lower_sum_a = _mm512_add_ps(avx512_lower_sum_a, avx512_local_lower_sum_a);

		__m512 avx512_local_lower_sum_b = _mm512_mul_ps(avx512_b, avx512_b);
		avx512_lower_sum_b = _mm512_add_ps(avx512_lower_sum_b, avx512_local_lower_sum_b);

		i += 16;
	}

	_mm512_storeu_ps(vec_upper_sum, avx512_upper_sum);
	_mm512_storeu_ps(vec_lower_sum_a, avx512_lower_sum_a);
	_mm512_storeu_ps(vec_lower_sum_b, avx512_lower_sum_b);

	for(int32_t j = 0 ; j < 16 ; j++){
		upper_sum += vec_upper_sum[j];
		lower_sum_a += vec_lower_sum_a[j];
		lower_sum_b += vec_lower_sum_b[j];
	}
	while(i < n){
		upper_sum += a[i] * b[i];
		lower_sum_a += powf(a[i], 2.0f);
		lower_sum_b += powf(b[i], 2.0f);
		i++;
	}
	// printf("upper_sum: %f; lower_sum_a: %f; lower_sum_b: %f\n", upper_sum, lower_sum_a, lower_sum_b);
	float cosine_similarity = (upper_sum / (sqrtf(lower_sum_a) * sqrtf(lower_sum_b)));
	return (1.0f - cosine_similarity);
}
#endif

double cosine_distance(const double* restrict const a, const double* restrict const b, int n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double upper_sum = 0.0;
	double lower_sum_a = 0.0;
	double lower_sum_b = 0.0;
	for(int i = 0 ; i < n ; i++){
		upper_sum += a[i] * b[i];
		lower_sum_a += pow(a[i], 2.0);
		lower_sum_b += pow(b[i], 2.0);
	}
	double cosine_similarity = (upper_sum / (pow(lower_sum_a, 0.5) * pow(lower_sum_b, 0.5)));
	return (1.0 - cosine_similarity);
}

double cosine_distance_norm(double* restrict a, double* restrict b, int32_t n){
	return cosine_distance(a, b, n) / 2.0;
}

double chebyshev_distance(double* restrict a, double* restrict b, int n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double result = 0.0;
	for(int i = 0 ; i < n ; i++){
		double diff = a[i] - b[i];
		if(diff < 0.0){
			diff *= -1.0;
		}
		if(diff > result){
			result = diff;
		}
	}
	return result;
}

double canberra_distance(double* restrict a, double* restrict b, int n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double result = 0.0;
	for(int i = 0 ; i < n ; i++){
		double diff = a[i] - b[i];
		if(diff < 0.0){
			diff *= -1.0;
		}
		double a_ = a[i];
		if(a_ < 0.0){
			a_ *= -1.0;
		}
		double b_ = b[i];
		if(b_ < 0.0){
			b_ *= -1.0;
		}
		if(a_ + b_ != 0.0){
			result += diff / (a_ + b_);
		}
	}
	return result;
}

double bray_curtis_distance(double* restrict a, double* restrict b, int n){
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double upper_sum = 0.0;
	double lower_sum = 0.0;
	for(int i = 0 ; i < n ; i++){
		if(a[i] < b[i]){
			upper_sum += a[i];
		}
		else{
			upper_sum += b[i];
		}
		lower_sum += a[i] + b[i];
	}
	return (1.0 - ((2.0 * upper_sum) / lower_sum));
}

/*
double rootless_angular_minkowski_distance(double* restrict a, double* restrict b, int n, double order){
	// Lenz and Cornelis (2023)
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double a_p_norm = 0.0;
	double b_p_norm = 0.0;
	for(int i = 0 ; i < n ; i++){
		double a_ = a[i];
		double b_ = b[i];
		a_ = pow(a_, order);
		b_ = pow(b_, order);
		if(a_ < 0.0){a_ *= -1.0;}
		if(b_ < 0.0){b_ *= -1.0;}
		a_p_norm += a_;
		b_p_norm += b_;
	}
	a_p_norm = pow(a_p_norm, 1.0 / order);
	b_p_norm = pow(b_p_norm, 1.0 / order);

	
	double global_norm = 0.0;
	for(int i = 0 ; i < n ; i++){
		double value = 0.0;
		if(b_p_norm != 0.0){
			value = (b[i] / b_p_norm);
		}
		if(a_p_norm != 0.0){
			value -= (a[i] / a_p_norm);
		}

		value = pow(value, order);
		if(value < 0.0){value *= -1.0;}
		global_norm += value;
	}
	// return global_norm;
	return (1.0 / order) * global_norm;
	
}

double angular_minkowski_distance(double* restrict a, double* restrict b, int n, double order){
	// Lenz and Cornelis (2023)
	return pow(rootless_angular_minkowski_distance(a, b, n, order), 1.0 / order);
}
*/

double angular_minkowski_distance(double* restrict a, double* restrict b, int32_t n, double order){
	// Lenz and Cornelis (2023)
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	double a_p_norm = 0.0;
	double b_p_norm = 0.0;
	for(int i = 0 ; i < n ; i++){
		double a_ = a[i];
		double b_ = b[i];
		a_ = pow(a_, order);
		b_ = pow(b_, order);
		if(a_ < 0.0){a_ *= -1.0;}
		if(b_ < 0.0){b_ *= -1.0;}
		a_p_norm += a_;
		b_p_norm += b_;
	}
	a_p_norm = pow(a_p_norm, 1.0 / order);
	b_p_norm = pow(b_p_norm, 1.0 / order);

	
	double global_norm = 0.0;
	for(int i = 0 ; i < n ; i++){
		double value = 0.0;
		if(b_p_norm != 0.0){
			value = (b[i] / b_p_norm);
		}
		if(a_p_norm != 0.0){
			value -= (a[i] / a_p_norm);
		}
		if(value < 0.0){
			value *= -1.0;
		}

		value = pow(value, order);
		if(value < 0.0){value *= -1.0;}
		global_norm += value;
	}
	/*
	// return global_norm;
	return (1.0 / order) * global_norm;
	*/
	double result = pow(global_norm, 1.0 / order);
	return result;
}

double rootless_angular_minkowski_distance(double* restrict a, double* restrict b, int32_t n, double order){
	return pow(angular_minkowski_distance(a, b, n, order), order);
}

float angular_minkowski_distance_fp32(float* restrict a, float* restrict b, int32_t n, float order){
	// Lenz and Cornelis (2023)
	#if MST_SANITY_TESTING == 1
	if(n > 2){n = 2;}
	#endif
	float a_p_norm = 0.0f;
	float b_p_norm = 0.0f;
	for(int i = 0 ; i < n ; i++){
		float a_ = a[i];
		float b_ = b[i];
		a_ = powf(a_, order);
		b_ = powf(b_, order);
		if(a_ < 0.0f){a_ *= -1.0f;}
		if(b_ < 0.0f){b_ *= -1.0f;}
		a_p_norm += a_;
		b_p_norm += b_;
	}
	a_p_norm = powf(a_p_norm, 1.0f / order);
	b_p_norm = powf(b_p_norm, 1.0f / order);

	
	float global_norm = 0.0f;
	for(int i = 0 ; i < n ; i++){
		float value = 0.0f;
		if(b_p_norm != 0.0f){
			value = (b[i] / b_p_norm);
		}
		if(a_p_norm != 0.0f){
			value -= (a[i] / a_p_norm);
		}
		if(value < 0.0f){
			value *= -1.0f;
		}

		value = powf(value, order);
		if(value < 0.0f){value *= -1.0f;}
		global_norm += value;
	}
	/*
	// return global_norm;
	return (1.0 / order) * global_norm;
	*/
	float result = powf(global_norm, 1.0f / order);
	return result;
}

float rootless_angular_minkowski_distance_fp32(float* restrict a, float* restrict b, int32_t n, float order){
	return powf(angular_minkowski_distance_fp32(a, b, n, order), order);
}

