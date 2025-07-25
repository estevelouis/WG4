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

#include <stdio.h>
extern int pclose (FILE *__stream); // somehow not detected

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h> // ?
#include <sys/stat.h>
#include <unistd.h> // ?
#if TOKENIZATION_METHOD == 0
#include <regex.h>
#endif

#include "jsonl/constants.h"
#include "sorted_array/array.h"
#include "sanitize.h"
#if (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
#include "udpipe/interface/cinterface.h" // external declarations
#include "udpipe/interface/conversion.h" // transformation to pipes
#endif

#include "jsonl/parser.h"

#if TOKENIZATION_METHOD == 0
// const char* const REGEX_STD_PATTERN = "[a-zéèêàç\\-]+'?|[\\.\\-,;:!\\?]";
// const char* const REGEX_STD_PATTERN = "([a-zéèêàç\\-]+'?)|[\\.\\-,;:!\\?]";
// const char* const REGEX_STD_PATTERN = "([a-zéèêàç\\-]+'?)|[,;:!\\?]"; // works
// const char* const REGEX_STD_PATTERN = "[a-zéèêàç\\-]+'?|[,;:!\\?]"; // works but not for accents?
// const char* const REGEX_STD_PATTERN = "[^\\s]+"; // approximate but faster
// const char* const REGEX_STD_PATTERN = "\\S(\\B\\S)*"; // slower
// const char* const REGEX_STD_PATTERN = "[a-z]+";
// const char* const REGEX_STD_PATTERN = "[^ \t\r\n]+'?"; // working?
const char* const REGEX_STD_PATTERN = "[^ \t\r\n,;:!\\?\\./']+'?";
const size_t NMATCH = 1;
// const size_t NMATCH = 2;
const int32_t REGEX_FLAGS = REG_EXTENDED | REG_ICASE;
#elif TOKENIZATION_METHOD == 1
const char* udpipe_repo_directory = "./udpipe";
const char* udpipe_model_directory = "./udpipe/sandbox_models";
#elif TOKENIZATION_METHOD == 2
#else
	#error "Unknown TOKENIZATION_METHOD"
#endif

void jsonl_init_tokenization(void){
	#if TOKENIZATION_METHOD == 1
	const int32_t bfr_size = 256;
	char bfr[bfr_size];
	memset(bfr, '\0', bfr_size);
	printf("UDPipe initialization ...\n");
	printf("Checking repository ...\n");
	snprintf(bfr, bfr_size, "if [ -d \"%s\" ]; then\ncurrentdirectory=$(pwd)\ncd %s\ngit pull\ncd src\nmake\ncd $currentdirectory\nelse\ngit clone https://github.com/ufal/udpipe %s\ncurrentdirectory=$(pwd)\ncd udpipe/src\nmake\ncd $currentdirectory\nfi\n", udpipe_repo_directory, udpipe_repo_directory, udpipe_repo_directory);
	// printf("trying to run: %s\n", bfr);
	if(system(bfr) != 0){goto system_failure;}
	printf("Creating model directory ...\n");

	struct stat local_st = (struct stat) {};
	if(stat(udpipe_model_directory, &local_st) == -1){
		if(mkdir(udpipe_model_directory, 0664) != 0){
			fprintf(stderr, "failed to create %s\n", udpipe_model_directory);
			exit(1);
		}
	}

	return;

	system_failure:
		fprintf(stderr, "Failed to launch process: %s\n", bfr);
		exit(1);
	#endif
}

#if (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
int32_t launch_udpipe(struct document* const doc){
	#if TOKENIZATION_METHOD == 1
	const char* model_name = "english-ewt-ud-2.5-191206.udpipe";
	// const char* model_name = "french-sequoia-ud-2.5-191206.udpipe";

	const int32_t path_model_bfr_size = 1024;
	char path_model_bfr[path_model_bfr_size];
	memset(path_model_bfr, '\0', path_model_bfr_size);
	snprintf(path_model_bfr, path_model_bfr_size, "%s/%s", udpipe_model_directory, model_name);

	char* sanitized_text;
	if(sanitize_for_shell(doc->text, doc->text_size, &sanitized_text) != 0){
		perror("failed to sanitize text\n");
		return 1;
	}
	const int32_t udpipe_call_bfr_size = 1024 + strlen(sanitized_text);
	size_t alloc_size = udpipe_call_bfr_size * sizeof(char);
	char* udpipe_call_bfr = malloc(alloc_size);
	if(udpipe_call_bfr == NULL){perror("malloc failed\n"); return 1;}
	memset(udpipe_call_bfr, '\0', udpipe_call_bfr_size);
	// snprintf(udpipe_call_bfr, udpipe_call_bfr_size, "echo \"%s\" | %s/src/udpipe --tokenize %s --output=vertical --input=horizontal\n", sanitized_text, udpipe_repo_directory, path_model_bfr); // works
	snprintf(udpipe_call_bfr, udpipe_call_bfr_size, "echo \"%s\" | %s/src/udpipe --tokenize %s --output=vertical --input=horizontal --immediate\n", sanitized_text, udpipe_repo_directory, path_model_bfr); // works
	doc->tmp_udpipe_output_file = popen(udpipe_call_bfr, "r");
	if(doc->tmp_udpipe_output_file == NULL){
		fprintf(stderr, "Failed to call UDPipe: %s\n", udpipe_call_bfr);
		free(udpipe_call_bfr);
		return 1;
	}
	free(sanitized_text); // malloc by sanitize_for_shell
	free(udpipe_call_bfr);
	#elif TOKENIZATION_METHOD == 2
	// if(transform_to_pipes(local_pipeline, doc->text, &(doc->tmp_udpipe_output_file), &(doc->heap_char_output), &(doc->output_pipefd)) != 0){
	if(transform_to_pipes(doc->text, &(doc->tmp_udpipe_output_file), &(doc->heap_char_output)) != 0){
		perror("Failed to call transform_to_pipes\n");
		return 1;
	}
	#endif
	return 0;
}
#endif

int32_t create_document(struct document* doc){
	#if TOKENIZATION_METHOD == 1
	/*
	static uint32_t document_counter = 0;
	if(document_counter < UINT32_MAX){
		document_counter++;
	} else {
		perror("document_counter reached UINT32_MAX\n");
		return 1;
	}
	*/
	#endif

	size_t malloc_size = JSONL_DOCUMENT_IDENTIFIER_DEFAULT_CAPACITY * sizeof(char);
	doc->identifier = (char*) malloc(malloc_size);
	if(doc->identifier == NULL){goto malloc_fail;}
	memset(doc->identifier, '\0', malloc_size);
	doc->identifier_capacity = JSONL_DOCUMENT_IDENTIFIER_DEFAULT_CAPACITY;
	doc->identifier_size = 0;

	malloc_size = JSONL_DOCUMENT_TEXT_DEFAULT_CAPACITY * sizeof(char);
	doc->text = (char*) malloc(malloc_size);
	if(doc->text == NULL){goto malloc_fail;}
	memset(doc->text, '\0', malloc_size);
	doc->text_capacity = JSONL_DOCUMENT_TEXT_DEFAULT_CAPACITY;
	doc->text_size = 0;

	#if TOKENIZATION_METHOD == 0
	if(regcomp(&(doc->reg), REGEX_STD_PATTERN, REGEX_FLAGS) != 0){
		perror("failed to call regcomp\n");
		return 1;
	}
	doc->latest_rm_eo = 0;
	#endif
	doc->reached_last_token = 0;
	doc->usable = 1;

	memset(doc->current_token, '\0', JSONL_CURRENT_TOKEN_BUFFER_SIZE);

	return 0;

	malloc_fail:
	if(doc->identifier != NULL){
		free(doc->identifier);
	}
	perror("[err] failed to malloc\n");
	return 1;
}

int32_t iterate_document_current_token(struct document* const doc){
	#if TOKENIZATION_METHOD == 0
	regmatch_t pmatch[NMATCH];
	int32_t err = regexec(&(doc->reg), doc->text + doc->latest_rm_eo, NMATCH, pmatch, 0);
	if(err == REG_NOMATCH){
		doc->reached_last_token = 1;
		return 0;
	} else if(err != 0){
		perror("failed to call regexec\n");
		return 1;
	}
	size_t len = pmatch[0].rm_eo - pmatch[0].rm_so;
	size_t objective_eo = doc->latest_rm_eo + pmatch[0].rm_eo;
	if(len > JSONL_CURRENT_TOKEN_BUFFER_SIZE - 1){
		len = JSONL_CURRENT_TOKEN_BUFFER_SIZE - 1;
	}
	memcpy(doc->current_token, doc->text + doc->latest_rm_eo + pmatch[0].rm_so, len);
	doc->current_token[len] = '\0';
	doc->latest_rm_eo = objective_eo;
	return 0;
	#elif (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
	if(doc->tmp_udpipe_output_file != NULL && feof(doc->tmp_udpipe_output_file)){
		doc->reached_last_token = 1;
		doc->usable = 0;
		pclose(doc->tmp_udpipe_output_file);
		doc->tmp_udpipe_output_file = NULL;
		return 0;
	}
	if(doc->tmp_udpipe_output_file == NULL){
		launch_udpipe(doc);
		doc->reached_last_token = 0;
		doc->usable = 0;
	}

	while(fgets(doc->current_token, JSONL_CURRENT_TOKEN_BUFFER_SIZE, doc->tmp_udpipe_output_file) && doc->current_token[0] == '\n'){}
	for(int32_t i = 0 ; i < JSONL_CURRENT_TOKEN_BUFFER_SIZE ; i++){
		if(doc->current_token[i] == '\n'){
			doc->current_token[i] = '\0';
			break;
		} else if(doc->current_token[i] == '\0'){
			break;
		}
	}
	// printf("newly found token: %s\n", doc->current_token);
	return 0;
	#endif
}


void free_document(struct document* doc){
	free(doc->identifier);
	free(doc->text);
	#if TOKENIZATION_METHOD == 0
	regfree(&(doc->reg));
	#elif TOKENIZATION_METHOD == 1
	pclose(doc->tmp_udpipe_output_file);
	#elif TOKENIZATION_METHOD == 2
	if(doc->tmp_udpipe_output_file != NULL){fclose(doc->tmp_udpipe_output_file);}
	if(doc->heap_char_output != NULL){free(doc->heap_char_output);}

	#endif
}

int32_t jsonl_document_iterator_request_more_capacity_current_line(struct jsonl_document_iterator * const jdi){
    const size_t new_capacity = jdi->current_line_capacity + JSONL_FILE_READ_BUFFER_SIZE;
    const size_t alloc_size = new_capacity * sizeof(char);
    void * realloc_ptr = realloc(jdi->current_line, alloc_size);
    if(realloc_ptr == NULL){
        perror("Failed to realloc\n");
        return 1;
    }
    jdi->current_line = realloc_ptr;
    memset(&(jdi->current_line[jdi->current_line_capacity]), '\0', JSONL_FILE_READ_BUFFER_SIZE * sizeof(char));

    jdi->current_line_capacity = new_capacity;

    return 0;
}

int32_t create_jsonl_document_iterator(struct jsonl_document_iterator* jdi, const char* file_name, const char* const content_key){
    if(strcmp(file_name, "-") == 0){
        jdi->file_ptr = stdin;
    } else {
	    jdi->file_ptr = fopen(file_name, "r");
    }
	if(jdi->file_ptr == NULL){
		fprintf(stderr, "[err] failed to open file (%s); errno: %i\n", file_name, errno);
		return 1;
	}
	jdi->file_is_open = 1;
	jdi->file_is_done = 0;
	memset(jdi->bfr_read, '\0', JSONL_FILE_READ_BUFFER_SIZE);
	memset(jdi->content_key, '\0', JSONL_CONTENT_KEY_BUFFER_SIZE);
	size_t bytes_to_cpy = strlen(content_key);
	if(bytes_to_cpy > JSONL_CONTENT_KEY_BUFFER_SIZE - 1){
		bytes_to_cpy = JSONL_CONTENT_KEY_BUFFER_SIZE - 1;
	}
	memcpy(jdi->content_key, content_key, bytes_to_cpy);

	memset(&(jdi->current_document), '\0', sizeof(struct document));

	if(create_document(&(jdi->current_document)) != 0){
		perror("[err] failed to call create_document\n");
		return 1;
	}

	jdi->document_to_free = 1;

    const size_t alloc_size = JSONL_FILE_READ_BUFFER_SIZE * sizeof(char);
    jdi->current_line = malloc(alloc_size);
    if(jdi->current_line == NULL){
        perror("Failed to malloc\n");
        return 1;
    }
    memset(jdi->current_line, '\0', alloc_size);
    jdi->current_line_cardinality = 0;
    jdi->current_line_capacity = JSONL_FILE_READ_BUFFER_SIZE;

	return 0;
}

void free_jsonl_document_iterator(struct jsonl_document_iterator* jdi){
	free_document(&(jdi->current_document));
	if(jdi->file_ptr != stdin){fclose(jdi->file_ptr);}
    free(jdi->current_line);
}

int32_t iterate_jsonl_document_iterator(struct jsonl_document_iterator* restrict const jdi){
	memset(jdi->bfr_read, '\0', JSONL_FILE_READ_BUFFER_SIZE);
    memset(jdi->current_line, '\0', jdi->current_line_capacity); // cardinality instead?
    jdi->current_line_cardinality = 0;

	int8_t in_key = 0;
	int8_t in_value = 0;
	int8_t escape_char = 0;

	char* key;
	char* value;
	int32_t key_size = 0;
	int32_t value_size = 0;

	const int32_t key_capacity_step = 16;
	const int32_t value_capacity_step = 64;

	int32_t key_capacity = key_capacity_step;
	int32_t value_capacity = value_capacity_step;

	int8_t colon_passed = 0;

	int8_t in_brackets = 0;

	#if TOKENIZATION_METHOD == 0
	jdi->current_document.latest_rm_eo = 0;
	#endif
	jdi->current_document.reached_last_token = 0;
	jdi->current_document.usable = 1;

	size_t malloc_size = key_capacity * sizeof(char);
	key = (char*) malloc(malloc_size);
	if(key == NULL){goto malloc_fail;}
	memset(key, '\0', malloc_size);

	malloc_size = value_capacity * sizeof(char);
	value = (char*) malloc(malloc_size);
	if(value == NULL){goto malloc_fail;}
	memset(value, '\0', malloc_size);

	if(feof(jdi->file_ptr)){
		jdi->file_is_done = 1;
		jdi->current_document.usable = 0;
		free(key);
		free(value);
		return 0;
	}

	while(fgets(jdi->bfr_read, JSONL_FILE_READ_BUFFER_SIZE, jdi->file_ptr) != NULL && strlen(jdi->bfr_read) != 0){
        // update jdi current line
        const size_t len_bfr_read = strlen(jdi->bfr_read);
        if(jdi->current_line_cardinality + len_bfr_read + 1 >= jdi->current_line_capacity){
            if(jsonl_document_iterator_request_more_capacity_current_line(jdi) != 0){
                perror("Failed to call jsonl_document_iterator_request_more_capacity_current_line\n");
                return 1;
            }
        }
        memcpy(&(jdi->current_line[jdi->current_line_cardinality]), jdi->bfr_read, len_bfr_read);
        jdi->current_line_cardinality += len_bfr_read;
        

		int32_t i = 0;
		int32_t bytes_read = 0;
		while(i < JSONL_FILE_READ_BUFFER_SIZE && jdi->bfr_read[i] != '\0'){

            // update jdi current line
            /*
            if(jdi->current_line_cardinality >= jdi->current_line_capacity - 1){ // -1 to ensure null byte at the end
                if(jsonl_document_iterator_request_more_capacity_current_line(jdi) != 0){
                    perror("Failed to call jsonl_document_iterator_request_more_capacity_current_line\n");
                    return 1;
                }
            }
            jdi->current_line[jdi->current_line_cardinality] = jdi->bfr_read[i];
            jdi->current_line_cardinality++;
            */

			bytes_read++;
			if(jdi->bfr_read[i] == '\n'){
				free(key);
				free(value);
				/* // DO NOT REMOVE
				if(strcmp(key, "id") == 0){
					printf("setting id: %s\n", value);
					if(jdi->current_document.identifier != NULL){free(jdi->current_document.identifier);}
					jdi->current_document.identifier = value;
					jdi->current_document.identifier_size = value_size;
					jdi->current_document.identifier_capacity = value_capacity;
				} else if(strcmp(key, "text") == 0){
					printf("setting text: %s\n", value);
					if(jdi->current_document.text != NULL){free(jdi->current_document.text);}
					jdi->current_document.text = value;
					jdi->current_document.text_size = value_size;
					jdi->current_document.text_capacity = value_capacity;
				} else {
					free(value);
				}
				free(key);
				*/
				if(feof(jdi->file_ptr)){
					jdi->file_is_done = 1;
				}

                // memset(jdi->current_line, '\0', jdi->current_line_cardinality * sizeof(char)); // should not be here

				return 0; // ?
			}

			if(in_brackets){
				while(i < JSONL_FILE_READ_BUFFER_SIZE - 1 && (escape_char || jdi->bfr_read[i] != ']') && jdi->bfr_read[i] != '\0'){
					if(escape_char){
						escape_char = 0;
						i++;
					} else {
						if(jdi->bfr_read[i] == '\\'){
							i++;
						}
						i++;
					}
				}	
				if(i < JSONL_FILE_READ_BUFFER_SIZE - 1 && jdi->bfr_read[i] == ']' && escape_char == 0){
					in_brackets = 0;
					in_value = 0;
				}
			} else if(in_key == 1){
				if(escape_char){
					if(key_size >= key_capacity - 1){
						int32_t new_capacity = key_capacity + key_capacity_step;
						size_t realloc_size = new_capacity * sizeof(char);
						void * realloc_ptr = (char*) realloc(key, realloc_size);
						if(realloc_ptr == NULL){goto realloc_fail;}
						key = realloc_ptr;
						memset(&(key[key_capacity]), '\0', key_capacity_step);
						key_capacity = new_capacity;
					}
					key[key_size] = jdi->bfr_read[i];
					key_size++;
					escape_char = 0;
				} else {
					switch(jdi->bfr_read[i]){
						case '\\':
							escape_char = 1;
							break;
						case '"':
							in_key = 0;
							break;
						default:
							if(key_size >= key_capacity - 1){
								int32_t new_capacity = key_capacity + key_capacity_step;
								size_t realloc_size = new_capacity * sizeof(char);
								void * realloc_ptr = (char*) realloc(key, realloc_size);
								if(realloc_ptr == NULL){goto realloc_fail;}
								key = realloc_ptr;
								memset(&(key[key_capacity]), '\0', key_capacity_step);
								key_capacity = new_capacity;
							}
							key[key_size] = jdi->bfr_read[i];
							key_size++;
							break;

					}
				}
			} else if(in_value == 1 || ('0' <= jdi->bfr_read[i] && jdi->bfr_read[i] <= '9')){
				if(escape_char){
					if(value_size >= value_capacity - 1){
						int32_t new_capacity = value_capacity + value_capacity_step;
						size_t realloc_size = new_capacity * sizeof(char);
						void * realloc_ptr = (char*) realloc(value, realloc_size);
						if(realloc_ptr == NULL){goto realloc_fail;}
						value = realloc_ptr;
						memset(&(value[value_capacity]), '\0', value_capacity_step);
						value_capacity = new_capacity;
					}
					switch(jdi->bfr_read[i]){
						case 'n':
							value[value_size] = '\n';
							break;
						case 'r':
							value[value_size] = '\r';
							break;
						case 't':
							value[value_size] = '\t';
							break;
						default:
							value[value_size] = jdi->bfr_read[i];
							break;
					}
					value_size++;
					escape_char = 0;
				} else {
					switch(jdi->bfr_read[i]){
						case '\\':
							escape_char = 1;
							break;
						case '"':
							if(strcmp(key, "id") == 0){
								jdi->current_document.identifier = realloc(jdi->current_document.identifier, value_size + 1);
								memcpy(jdi->current_document.identifier, value, value_size);
								jdi->current_document.identifier[value_size] = '\0';
								jdi->current_document.identifier_size = value_size;
								jdi->current_document.identifier_capacity = value_size + 1;
							} else if(strcmp(key, jdi->content_key) == 0){
								jdi->current_document.text = realloc(jdi->current_document.text, value_size + 1);
								memcpy(jdi->current_document.text, value, value_size);
								jdi->current_document.text[value_size] = '\0';
								jdi->current_document.text_size = value_size;
								jdi->current_document.text_capacity = value_size + 1;
							}
							in_value = 0;
							colon_passed = 0;
							break;
						default:
							if(value_size >= value_capacity - 1){
								int32_t new_capacity = value_capacity + value_capacity_step;
								size_t realloc_size = new_capacity * sizeof(char);
								void * realloc_ptr = (char*) realloc(value, realloc_size);
								if(realloc_ptr == NULL){goto realloc_fail;}
								value = realloc_ptr;
								memset(&(value[value_capacity]), '\0', value_capacity_step);
								value_capacity = new_capacity;
							}
							value[value_size] = jdi->bfr_read[i];
							value_size++;
							break;

					}
				}
			} else if(jdi->bfr_read[i] == ':'){
				colon_passed = 1;
			} else if(jdi->bfr_read[i] == '"'){
				if(colon_passed){
					in_value = 1;
				} else {
					memset(key, '\0', key_size);
					key_size = 0;
					memset(value, '\0', value_size);
					value_size = 0;
					in_key = 1;
				}
			} else if(jdi->bfr_read[i] == '['){
				i++;
				while(i < JSONL_FILE_READ_BUFFER_SIZE - 1 && (escape_char || jdi->bfr_read[i] != ']') && jdi->bfr_read[i] != '\0'){
					if(escape_char){
						escape_char = 0;
						i++;
					} else {
						if(jdi->bfr_read[i] == '\\'){
							// i++;
							escape_char = 1;
						}
						i++;
					}
				}
				if(i < JSONL_FILE_READ_BUFFER_SIZE - 1 && jdi->bfr_read[i] != ']' && escape_char == 0){
					in_brackets = 1;
				}
			} else if(jdi->bfr_read[i] == ',' || jdi->bfr_read[i] == '}'){
				if(strcmp(key, "id") == 0){
					jdi->current_document.identifier = realloc(jdi->current_document.identifier, value_size + 1);
					memcpy(jdi->current_document.identifier, value, value_size);
					jdi->current_document.identifier[value_size] = '\0';
					jdi->current_document.identifier_size = value_size;
					jdi->current_document.identifier_capacity = value_size + 1;
				} else if(strcmp(key, jdi->content_key) == 0){
					jdi->current_document.text = realloc(jdi->current_document.text, value_size + 1);
					memcpy(jdi->current_document.text, value, value_size);
					jdi->current_document.text[value_size] = '\0';
					jdi->current_document.text_size = value_size;
					jdi->current_document.text_capacity = value_size + 1;
				}
				in_value = 0;
				colon_passed = 0;
			}
			i++;
		}
		if(feof(jdi->file_ptr)){
			jdi->file_is_done = 1;
		} else {
			memset(jdi->bfr_read, '\0', bytes_read);
		}
	} 
	free(key);
	free(value);

	if(feof(jdi->file_ptr)){
		jdi->file_is_done = 1;
	}

	#if TOKENIZATION_METHOD == 1
	if(launch_udpipe(&(jdi->current_document)) != 0){
		perror("Failed to call launch_udpipe\n");
		return 1;
	}
	#endif

	return 0;

	realloc_fail:
	perror("[err] failed to realloc\n");
	return 1;
	malloc_fail:
	perror("[err] failed to malloc\n");
	return 1;
}

