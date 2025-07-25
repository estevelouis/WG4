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

#ifndef UDPIPE_INTERFACE_HPP
#define UDPIPE_INTERFACE_HPP

#include <iostream>
#include <string>
#include <sstream> // std::istringstream

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

extern "C" void custom_format(const char* const log_type, const int32_t log_num, const char* const filename, const char* const func, const int32_t line, const char* const msg);
extern "C" void error_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);
extern "C" void warning_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);
extern "C" void info_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);

#include "udpipe.h"
#include "udpipe/interface/size.h"

#include "udpipe/sources/source.h"
#include "udpipe/sources/udpipe1/all.h"

#ifndef UDPIPE_STORAGE_DIR
#define UDPIPE_STORAGE_DIR "$HOME/.local/udpipe/udpipe1/models/"
#endif

static ufal::udpipe::pipeline* global_pipeline;

extern "C" void ensure_proper_udpipe_pipeline_size(){
	if(UDPIPE_PIPELINE_SIZE != sizeof(ufal::udpipe::pipeline)){
		fprintf(stderr, "UDPIPE_PIPELINE_SIZE != sizeof(ufal::udpipe::pipeline); %i != %li\n", UDPIPE_PIPELINE_SIZE, sizeof(ufal::udpipe::pipeline));
		exit(1);
	}
}

extern "C" int32_t ensure_udpipe_model_is_available(const char * const isolang, char * const bfr_p, const size_t bfr_p_size){
    size_t n = sizeof(all_udpipe_sources) / sizeof(all_udpipe_sources[0]);
    for(size_t i = 0 ; i < n ; i++){
        if(strcmp(isolang, all_udpipe_sources[i].isolang) == 0){
            const int32_t bfr_size = 512;
            char bfr[bfr_size];

            memset(bfr, '\0', bfr_size);
            snprintf(bfr, bfr_size, "set -eu ; mkdir -p %s", UDPIPE_STORAGE_DIR);
            info_format(__FILE__, __func__, __LINE__, bfr);
            if(system(bfr) != 0){
                fprintf(stderr, "Failed to launch \"%s\"\n", bfr);
                return 1;
            }

            memset(bfr, '\0', bfr_size);
            snprintf(bfr, bfr_size, "set -eu ; if ! [ -f \"%s%s\" ] ; then wget --no-clobber -O %s%s %s%s ; fi", UDPIPE_STORAGE_DIR, all_udpipe_sources[i].filename, UDPIPE_STORAGE_DIR, all_udpipe_sources[i].filename, all_udpipe_sources[i].url_download, all_udpipe_sources[i].filename);
            info_format(__FILE__, __func__, __LINE__, bfr);
            if(system(bfr) != 0){
                fprintf(stderr, "Failed to launch \"%s\"\n", bfr);
                return 1;
            }

            memset(bfr, '\0', bfr_size);
            snprintf(bfr, bfr_size, "set -eu ; cd %s ; md5sum --ignore-missing -c %s/src/include/udpipe/sources/udpipe1/md5_checksums.txt", UDPIPE_STORAGE_DIR, COMPILATION_DIR);
            info_format(__FILE__, __func__, __LINE__, bfr);
            if(system(bfr) != 0){
                fprintf(stderr, "Failed to launch \"%s\"\n", bfr);
                return 1;
            }

            memset(bfr_p, '\0', bfr_p_size);
            snprintf(bfr_p, bfr_p_size, "%s%s", UDPIPE_STORAGE_DIR, all_udpipe_sources[i].filename);

            return 0;
        }
    }
    
    fprintf(stderr, "Unknown lang: \"%s\"\n", isolang);
    return 1;
}

extern "C" void udpipe_pipeline_create_global(const char* model_name, const char* input, const char* tagger, const char* parser, const char* output){
    const size_t bfr_size = 512;
    char bfr[bfr_size];

    memset(bfr, '\0', bfr_size);
    snprintf(bfr, bfr_size, "Creating pipeline with \"%s\"", model_name);
    info_format(__FILE__, __func__, __LINE__, bfr);

    if(strlen(model_name) == 3){
        if(ensure_udpipe_model_is_available(model_name, bfr, bfr_size) != 0){
            fprintf(stderr, "Failed to call ensure_udpipe_model_is_available for \"%s\"\n", model_name);
            exit(1);
        }
    } else {
        memset(bfr, '\0', bfr_size);
        snprintf(bfr, bfr_size, "%s", model_name);
    }

    info_format(__FILE__, __func__, __LINE__, bfr);

	static ufal::udpipe::pipeline local_pipeline = ufal::udpipe::pipeline(
		ufal::udpipe::model::load(bfr),
		std::string(input),
		std::string(tagger),
		std::string(parser),
		std::string(output)
	);
	global_pipeline = std::addressof(local_pipeline);
}

extern "C" void udpipe_pipeline_process(const char* raw_txt, FILE** const pointer_file, void** const pointer_heap_char){
	std::string raw_txt_to_string(raw_txt);
	std::istringstream iss (raw_txt_to_string);
	std::ostringstream oss("");
	std::string error_string;

	global_pipeline->process(iss, oss, error_string);

	if(error_string != ""){
		std::cerr << "error: " << error_string << std::endl;
		exit(1);
	}

	std::string str_from_oss = oss.str();
	size_t str_from_oss_length = str_from_oss.length();
	
	(*pointer_heap_char) = realloc((*pointer_heap_char), (str_from_oss_length + 1) * sizeof(char));
	if((*pointer_heap_char) == NULL){
		perror("realloc failed\n");
		exit(1);
	}
	memcpy((*pointer_heap_char), str_from_oss.c_str(), str_from_oss_length);
	((char*) (*pointer_heap_char))[str_from_oss_length] = '\0';
	if(pointer_file != NULL){
	  (*pointer_file) = fmemopen((void*) (*pointer_heap_char), str_from_oss_length, "r");
	}
}

#endif
