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

#define TOKEN_ID_RAW_SIZE 8
#define TOKEN_FORM_SIZE 16
#define TOKEN_LEMMA_SIZE 16
#define TOKEN_UPOS_SIZE 8
#define TOKEN_XPOS_SIZE 16
#define TOKEN_FEATS_SIZE 1024
#define TOKEN_HEAD_SIZE 16
#define TOKEN_DEPREL_SIZE 16
#define TOKEN_DEPS_SIZE 16
#define TOKEN_MISC_SIZE 16
#define TOKEN_MWE_SIZE 16

#define TOKEN_SERIALIZE_SIZE (TOKEN_ID_RAW_SIZE + TOKEN_FORM_SIZE + TOKEN_LEMMA_SIZE + TOKEN_UPOS_SIZE + TOKEN_XPOS_SIZE + TOKEN_FEATS_SIZE + TOKEN_HEAD_SIZE + TOKEN_DEPREL_SIZE + TOKEN_DEPS_SIZE + TOKEN_MISC_SIZE + TOKEN_MWE_SIZE + 10 + 1)



#define SENTENCE_TOKEN_CAPACITY_STEP 16
/* #define SENTENCE_TOKEN_CAPACITY_STEP 32 */
/* #define SENTENCE_TOKEN_CAPACITY_STEP 64 */
/* #define SENTENCE_TOKEN_CAPACITY_STEP 512 */
#define SENTENCE_SENT_ID_SIZE 4096

// #define FILE_READ_BUFFER_SIZE 256
#define FILE_READ_BUFFER_SIZE 2048

#endif
