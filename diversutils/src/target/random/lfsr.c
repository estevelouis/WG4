/*
 *      DiversUtils - Functions to measure diversity
 *
 * Copyright (c) 2025  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
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

#include <stdint.h>
#include <stdio.h>

#include "random/lfsr.h"

int32_t lfsr_init(struct lfsr * restrict const rng, const uint64_t seed){
    if(seed == 0){
        perror("LFSR seed cannot be 0\n");
        return 1;
    }
    rng->status = seed;
    return 0;
}

void lfsr_rand(struct lfsr * restrict const rng){
    uint64_t output = 0;
    for(uint64_t i = 0 ; i < 64 ; i++){
        uint64_t new_output_bit = rng->status % 2;
        output |= new_output_bit << i;

        // https://datacipy.cz/lfsr_table.pdf
        uint64_t new_status_bit = (rng->status ^ (rng->status >> 1) ^ (rng->status >> 3) ^ (rng->status >> 4) ^ 1);
        rng->status = (rng->status >> 1) | (new_status_bit << 63);
    }

    rng->output = output;
}
