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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "udpipe/download.h"
#include "udpipe/sources/source.h"
#include "udpipe/sources/udpipe1/all.h"

int32_t ensure_udpipe_model_is_available(const char * const isolang){
    size_t n = sizeof(all_udpipe_sources) / sizeof(all_udpipe_sources[0]);
    for(size_t i = 0 ; i < n ; i++){
        if(strcmp(isolang, all_udpipe_sources[i].isolang) == 0){
            const int32_t bfr_size = 512;
            char bfr[bfr_size];
            memset(bfr, '\0', bfr_size);

            snprintf(bfr, bfr_size, "curl -o data/udpipe/udpipe1/%s %s%s", all_udpipe_sources[i].filename, all_udpipe_sources[i].url_download, all_udpipe_sources[i].filename);

            puts(bfr);

            if(system(bfr) != 0){
                fprintf(stderr, "Failed to launch \"%s\"\n", bfr);
                return 1;
            }

            return 0;
        }
    }
    
    fprintf(stderr, "Unknown lang: \"%s\"\n", isolang);
    return 1;
}
