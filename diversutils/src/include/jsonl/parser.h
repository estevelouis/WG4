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

#ifndef JSONL_PARSER_H
#define JSONL_PARSER_H

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<errno.h>
#include<regex.h>

#include "jsonl/constants.h"

#include "sorted_array/array.h"

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

struct document {
	int32_t identifier_size;
	int32_t text_size;
	int32_t identifier_capacity;
	int32_t text_capacity;
	regex_t reg;
	regoff_t latest_rm_eo;
	int8_t reached_last_token;
	int8_t usable;
	char current_token[JSONL_CURRENT_TOKEN_BUFFER_SIZE];
	char* identifier;
	char* text;
};

int32_t create_document(struct document* doc){
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

	if(regcomp(&(doc->reg), REGEX_STD_PATTERN, REGEX_FLAGS) != 0){
		perror("failed to call regcomp\n");
		return 1;
	}
	doc->latest_rm_eo = 0;
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
}


void free_document(struct document* doc){
	free(doc->identifier);
	free(doc->text);
	regfree(&(doc->reg));
}

struct jsonl_document_iterator {
	FILE* file_ptr;
	int8_t file_is_open;
	int8_t file_is_done;
	char bfr_read[JSONL_FILE_READ_BUFFER_SIZE];
	char content_key[JSONL_CONTENT_KEY_BUFFER_SIZE];
	struct document current_document;
	int32_t document_to_free;
};

int32_t create_jsonl_document_iterator(struct jsonl_document_iterator* jdi, const char* file_name, const char* const content_key){
	jdi->file_ptr = fopen(file_name, "r");
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

	return 0;
}

void free_jsonl_document_iterator(struct jsonl_document_iterator* jdi){
	free_document(&(jdi->current_document));
	fclose(jdi->file_ptr);
}

int32_t iterate_jsonl_document_iterator(struct jsonl_document_iterator* restrict const jdi){
	memset(jdi->bfr_read, '\0', JSONL_FILE_READ_BUFFER_SIZE);

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

	jdi->current_document.latest_rm_eo = 0;
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
		int32_t i = 0;
		int32_t bytes_read = 0;
		while(jdi->bfr_read[i] != '\0'){
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
						key = (char*) realloc(key, realloc_size);
						if(key == NULL){goto realloc_fail;}
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
								key = (char*) realloc(key, realloc_size);
								if(key == NULL){goto realloc_fail;}
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
						value = (char*) realloc(value, realloc_size);
						if(value == NULL){goto realloc_fail;}
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
								value = (char*) realloc(value, realloc_size);
								if(value == NULL){goto realloc_fail;}
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

	return 0;

	realloc_fail:
	perror("[err] failed to realloc\n");
	return 1;
	malloc_fail:
	perror("[err] failed to malloc\n");
	return 1;
}

#endif
				
