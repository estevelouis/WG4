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

#include <stddef.h> // size_t
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sanitize.h"

int32_t sanitize_for_shell(const char* const s, const size_t n, char** const z){
	size_t alloc_size = n + 1024;
	size_t i, j;

	(*z) = malloc(alloc_size * sizeof(char));
	if((*z) == NULL){
		perror("malloc failed\n");
		return 1;
	}
	memset((*z), '\0', alloc_size * sizeof(char));

	j = 0;
	for(i = 0 ; i < n ; i++){
		if(j + 2 >= alloc_size){alloc_size += 1024; (*z) = realloc((*z), alloc_size * sizeof(char)); if((*z) == NULL){perror("failed to realloc\n"); return 1;}}
		switch(s[i]){
			case '\\':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = '\\';
				j++;
				break;
			case '\n':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = 'n';
				j++;
				break;
			case '\r':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = 'r';
				j++;
				break;
			case '"':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = '"';
				j++;
				break;
			default:
				(*z)[j] = s[i];
				j++;
				break;
		}
	}
	return 0;
}

