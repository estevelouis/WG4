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

#ifndef CUPT_CONSTANTS_H
#define CUPT_CONSTANTS_H

#include<stdint.h>

const int32_t TOKEN_BUFFER_SIZE = 32;
const int32_t SENTENCE_BUFFER_SIZE = 128;

const int32_t TOKEN_ID_RAW_SIZE = 16;
const int32_t TOKEN_FORM_SIZE = 32;
const int32_t TOKEN_LEMMA_SIZE = 32;
const int32_t TOKEN_UPOS_SIZE = 16;
const int32_t TOKEN_XPOS_SIZE = 16;
const int32_t TOKEN_FEATS_SIZE = 64;
const int32_t TOKEN_DEPREL_SIZE = 32;
const int32_t TOKEN_DEPS_SIZE = 64;
const int32_t TOKEN_MWE_SIZE = 32;

#endif
