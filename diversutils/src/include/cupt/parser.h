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

#ifndef CUPT_PARSER_H
#define CUPT_PARSER_H

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>

#include "cupt/constants.h"

/* BASED ON UD FORMAT */
struct token {
	char id_raw[TOKEN_ID_RAW_SIZE];
	char form[TOKEN_FORM_SIZE];
	char lemma[TOKEN_LEMMA_SIZE];
	char upos[TOKEN_UPOS_SIZE];
	char xpos[TOKEN_XPOS_SIZE];
	char feats[TOKEN_FEATS_SIZE];
	char head[TOKEN_HEAD_SIZE];
	char deprel[TOKEN_DEPREL_SIZE];
	char deps[TOKEN_DEPS_SIZE];
	char misc[TOKEN_MISC_SIZE];
	char mwe[TOKEN_MWE_SIZE];
};


int32_t create_token(struct token* t, char* id_raw, char* form, char* lemma, char* upos, char* xpos, char* feats, char* head, char* deprel, char* deps, char* misc, char* mwe){
	memset(t->id_raw, '\0', TOKEN_ID_RAW_SIZE);
	memset(t->form, '\0', TOKEN_FORM_SIZE);
	memset(t->lemma, '\0', TOKEN_LEMMA_SIZE);
	memset(t->upos, '\0', TOKEN_UPOS_SIZE);
	memset(t->xpos, '\0', TOKEN_XPOS_SIZE);
	memset(t->feats, '\0', TOKEN_FEATS_SIZE);
	memset(t->head, '\0', TOKEN_HEAD_SIZE);
	memset(t->deprel, '\0', TOKEN_DEPREL_SIZE);
	memset(t->deps, '\0', TOKEN_DEPS_SIZE);
	memset(t->misc, '\0', TOKEN_MISC_SIZE);
	memset(t->mwe, '\0', TOKEN_MWE_SIZE);

	strncpy(t->id_raw, id_raw, TOKEN_ID_RAW_SIZE - 1);
	strncpy(t->form, form, TOKEN_FORM_SIZE - 1);
	strncpy(t->lemma, lemma, TOKEN_LEMMA_SIZE - 1);
	strncpy(t->upos, upos, TOKEN_UPOS_SIZE - 1);
	strncpy(t->xpos, xpos, TOKEN_XPOS_SIZE - 1);
	strncpy(t->feats, feats, TOKEN_FEATS_SIZE - 1);
	strncpy(t->head, head, TOKEN_HEAD_SIZE - 1);
	strncpy(t->deprel, deprel, TOKEN_DEPREL_SIZE - 1);
	strncpy(t->deps, deps, TOKEN_DEPS_SIZE - 1);
	strncpy(t->misc, misc, TOKEN_MISC_SIZE - 1);
	strncpy(t->mwe, mwe, TOKEN_MWE_SIZE - 1);

	return 0;
}

int32_t serialize_token(struct token* t, char* bfr, int32_t* bytes_written, int8_t null_terminate){
	size_t offset = 0;
	size_t length = 0;

	length = strlen(t->id_raw);
	strcpy(&(bfr[offset]), t->id_raw);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->form);
	strcpy(&(bfr[offset]), t->form);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->lemma);
	strcpy(&(bfr[offset]), t->lemma);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->upos);
	strcpy(&(bfr[offset]), t->upos);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->xpos);
	strcpy(&(bfr[offset]), t->xpos);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->feats);
	strcpy(&(bfr[offset]), t->feats);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->head);
	strcpy(&(bfr[offset]), t->head);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->deprel);
	strcpy(&(bfr[offset]), t->deprel);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->deps);
	strcpy(&(bfr[offset]), t->deps);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->misc);
	strcpy(&(bfr[offset]), t->misc);
	offset += length;
	bfr[offset] = '\t';
	offset++;

	length = strlen(t->mwe);
	strcpy(&(bfr[offset]), t->mwe);
	offset += length;
	if(null_terminate){
		bfr[offset] = '\0';
		offset++;
	}

	(*bytes_written) = (int32_t) offset;

	return 0;
}


/* ======== SENTENCE ======== */

struct sentence {
	struct token* tokens;
	int32_t capacity;
	int32_t num_tokens;
	char sentence_id[SENTENCE_SENT_ID_SIZE];
};

int32_t create_sentence(struct sentence* s, char* sentence_id){
	void* malloc_pointer;

	memset(s->sentence_id, '\0', SENTENCE_SENT_ID_SIZE);

	strncpy(s->sentence_id, sentence_id, SENTENCE_SENT_ID_SIZE - 1);

	size_t malloc_size = SENTENCE_TOKEN_CAPACITY_STEP * sizeof(struct token);
	malloc_pointer = malloc(malloc_size);
	if(malloc_pointer == NULL){
		perror("malloc failed\n");
		return 1;
	}
	memset(malloc_pointer, '\0', malloc_size);
	s->tokens = (struct token*) malloc_pointer;
	s->capacity = SENTENCE_TOKEN_CAPACITY_STEP;
	s->num_tokens = 0;

	if(s->tokens == NULL){
		printf("(create_sentence) s->tokens == NULL\n");
		exit(1);
	}

	return 0;
}

int32_t recreate_sentence(struct sentence* s, char* sentence_id){
	memset(s->sentence_id, '\0', SENTENCE_SENT_ID_SIZE);

	strncpy(s->sentence_id, sentence_id, SENTENCE_SENT_ID_SIZE - 1);

	size_t malloc_size = SENTENCE_TOKEN_CAPACITY_STEP * sizeof(struct token);
	if(s->tokens == NULL){
		printf("s->tokens == NULL\n");
		printf("s->sentence_id: %s\n", s->sentence_id);
		exit(1);
		/* // do not remove
		void* malloc_pointer = malloc(malloc_size);
		if(malloc_pointer == NULL){
			perror("malloc failed\n");
			return 1;
		}
		printf("malloc: %p\n", malloc_pointer);
		s->tokens = (struct token*) malloc_pointer;
		*/
	} else {
		s->tokens = (struct token*) realloc((void*) s->tokens, malloc_size);
		if(s->tokens == NULL){
			perror("malloc failed\n");
			return 1;
		}
	}
	memset(s->tokens, '\0', malloc_size);
	s->capacity = SENTENCE_TOKEN_CAPACITY_STEP;
	s->num_tokens = 0;

	if(s->tokens == NULL){
		printf("(recreate_sentence) s->tokens == NULL\n");
		exit(1);
	}

	return 0;
}

void free_sentence(struct sentence* s){
	free(s->tokens);
}

int32_t add_token_to_sentence(struct sentence* s, struct token* t){
	if(s->num_tokens >= s->capacity){
		size_t new_size = (size_t) ((((size_t) s->capacity) + SENTENCE_TOKEN_CAPACITY_STEP) * sizeof(struct token)); // compliance
		s->tokens = realloc(s->tokens, new_size);
		if(s->tokens == NULL){
			perror("realloc failed\n");
			return 1;
		}
		memset(&(s->tokens[s->capacity]), '\0', (size_t) (new_size - (((size_t) s->capacity) * sizeof(struct token))));
		s->capacity += SENTENCE_TOKEN_CAPACITY_STEP;
	}
	s->tokens[s->num_tokens] = *t;
	s->num_tokens++;
	return 0;
}

size_t buffer_size_for_serialization_of_sentence(struct sentence* s){
	size_t output = 0;
	output += 10; /* # sent_id= */
	output += SENTENCE_SENT_ID_SIZE;
	output += sizeof(struct token) * ((size_t) s->num_tokens);
	return output;
}

int32_t serialize_sentence(struct sentence* s, char* bfr, int32_t* bytes_written){
	size_t offset = 0;
	size_t length = 0;
	char* sent_id_init = "# sent_id=";
	int32_t bytes_written_by_token = 0;
	int32_t err = 0;

	length = strlen(sent_id_init);
	strcpy(&(bfr[offset]), sent_id_init);
	offset += length;

	length = strlen(s->sentence_id);
	strcpy(&(bfr[offset]), s->sentence_id);
	offset += length;

	bfr[offset] = '\n';
	offset++;
	for(int32_t i = 0 ; i < s->num_tokens ; i++){
		err = serialize_token(&(s->tokens[i]), &(bfr[offset]), &bytes_written_by_token, i == s->num_tokens - 1);
		if(err != 0){
			perror("failed to serialize token during sentence serialization\n");
			return 1;
		}
		offset += (size_t) bytes_written_by_token;
		bfr[offset] = '\n';
		offset++;
	}

	(*bytes_written) = (int32_t) offset;

	return 0;
}

/* ======== FILE READER ======= */

struct cupt_sentence_iterator {
	FILE* file_ptr;
	int8_t file_is_open;
	int8_t file_is_done;
	char bfr_read[FILE_READ_BUFFER_SIZE];
	struct sentence current_sentence;
	int32_t sentence_to_free;
};

int32_t create_cupt_sentence_iterator(struct cupt_sentence_iterator* const csi, const char* const file_name){
	int32_t err = 0;

	FILE* fp = fopen(file_name, "r");
	if(fp == NULL){
		perror("failed to open file\n");
		printf("errno: %i\n", errno);
		return 1;
	}
	csi->file_ptr = fp;
	csi->file_is_open = 1;
	csi->file_is_done = 0;
	memset(csi->bfr_read, '\0', FILE_READ_BUFFER_SIZE);

	memset(&(csi->current_sentence), '\0', sizeof(struct sentence));

	csi->current_sentence.tokens = NULL;
	
	err = create_sentence(&(csi->current_sentence), "default_sentence");
	if(err != 0){
		perror("failed to create default_sentence in create_cupt_sentence_iterator\n");
		return 1;
	}

	csi->sentence_to_free = 1;

	return 0;
}

void free_cupt_sentence_iterator(struct cupt_sentence_iterator* csi){
	fclose(csi->file_ptr);
}

int32_t iterate_cupt_sentence_iterator(struct cupt_sentence_iterator* restrict csi){
	if(csi->current_sentence.tokens == NULL){
		printf("csi->current_sentence.tokens == NULL\n");
		exit(1);
	}
	static int32_t num_calls = 0;
	num_calls++;
	const char* sent_id_pattern = "# sent_id = ";
	const size_t sent_id_pattern_len = strlen(sent_id_pattern);
	const char* source_sent_id_pattern = "# source_sent_id = ";
	const size_t source_sent_id_pattern_len = strlen(source_sent_id_pattern);

	int32_t err = 0;
	void* strstr_ptr;
	char new_sent_id[SENTENCE_SENT_ID_SIZE];

	memset(new_sent_id, '\0', SENTENCE_SENT_ID_SIZE);

	while(fgets(csi->bfr_read, FILE_READ_BUFFER_SIZE, csi->file_ptr)){
		strstr_ptr = strstr(csi->bfr_read, sent_id_pattern);
		int64_t len_to_consider;
		if(strstr_ptr == NULL){
			strstr_ptr = strstr(csi->bfr_read, source_sent_id_pattern);
			if(strstr_ptr == NULL){
				continue;
			}
			len_to_consider = (int64_t) source_sent_id_pattern_len; // compliance
		} else {
			len_to_consider = (int64_t) sent_id_pattern_len; // compliance
		}
		char* p_null = strchr(&(csi->bfr_read[len_to_consider]), '\0');
		char* p_linefeed = strchr(&(csi->bfr_read[len_to_consider]), '\n');
		size_t bytes_to_cpy = (size_t) (p_null - &(csi->bfr_read[len_to_consider]));
		if(p_linefeed < p_null){
			bytes_to_cpy = (size_t) (p_linefeed - &(csi->bfr_read[len_to_consider]));
		}
		if(SENTENCE_SENT_ID_SIZE - 1 < bytes_to_cpy){
			bytes_to_cpy = SENTENCE_SENT_ID_SIZE - 1;
		}
		strncpy(new_sent_id, &(csi->bfr_read[len_to_consider]), bytes_to_cpy);
		break;
	}
	err = recreate_sentence(&(csi->current_sentence), new_sent_id);
	if(err != 0){
		perror("failed to create sentence in file parsing\n");
		return 1;
	}

	char id_raw[TOKEN_ID_RAW_SIZE];
	char form[TOKEN_FORM_SIZE];
	char lemma[TOKEN_LEMMA_SIZE];
	char upos[TOKEN_UPOS_SIZE];
	char xpos[TOKEN_XPOS_SIZE];
	char feats[TOKEN_FEATS_SIZE];
	char head[TOKEN_HEAD_SIZE];
	char deprel[TOKEN_DEPREL_SIZE];
	char deps[TOKEN_DEPS_SIZE];
	char misc[TOKEN_MISC_SIZE];
	char mwe[TOKEN_MWE_SIZE];
	memset(id_raw, '\0', TOKEN_ID_RAW_SIZE);
	memset(form, '\0', TOKEN_FORM_SIZE);
	memset(lemma, '\0', TOKEN_LEMMA_SIZE);
	memset(upos, '\0', TOKEN_UPOS_SIZE);
	memset(xpos, '\0', TOKEN_XPOS_SIZE);
	memset(feats, '\0', TOKEN_FEATS_SIZE);
	memset(head, '\0', TOKEN_HEAD_SIZE);
	memset(deprel, '\0', TOKEN_DEPREL_SIZE);
	memset(deps, '\0', TOKEN_DEPS_SIZE);
	memset(misc, '\0', TOKEN_MISC_SIZE);
	memset(mwe, '\0', TOKEN_MWE_SIZE);

	int32_t current_column = 0;
	int32_t current_index_in_column = 0;
	int32_t column_sizes[11] = {TOKEN_ID_RAW_SIZE, TOKEN_FORM_SIZE, TOKEN_LEMMA_SIZE, TOKEN_UPOS_SIZE, TOKEN_XPOS_SIZE, TOKEN_FEATS_SIZE, TOKEN_HEAD_SIZE, TOKEN_DEPREL_SIZE, TOKEN_DEPS_SIZE, TOKEN_MISC_SIZE, TOKEN_MWE_SIZE};
	char* columns[11] = {id_raw, form, lemma, upos, xpos, feats, head, deprel, deps, misc, mwe};

	// printf("number of columns: %li\n", sizeof(column_sizes) / sizeof(column_sizes[0]));

	int8_t previous_line_finished = 1;
	if(feof(csi->file_ptr)){
		csi->file_is_done = 1;
		return 0;
	}
	int8_t currently_in_text_line = 0;
	while(fgets(csi->bfr_read, FILE_READ_BUFFER_SIZE, csi->file_ptr)){
		if(previous_line_finished && (csi->bfr_read[0] == '\n' || csi->bfr_read[0] == '\r')){
			break;
		}
		if(currently_in_text_line == 1){
			for(int32_t i = 0 ; i < FILE_READ_BUFFER_SIZE ; i++){
				if(csi->bfr_read[i] == '\0'){
					break;
				}
				if(csi->bfr_read[i] == '\r' || csi->bfr_read[i] == '\n'){
					currently_in_text_line = 0;
					break;
				}
			}
			memset(csi->bfr_read, '\0', FILE_READ_BUFFER_SIZE);
			continue;
		}
		if(strncmp(csi->bfr_read, "# text =", strlen("# text =")) == 0 || strncmp(csi->bfr_read, "# text=", strlen("# text=")) == 0){
			if(strchr(csi->bfr_read, '\n') != NULL || strchr(csi->bfr_read, '\r') != NULL){
				currently_in_text_line = 0;
			} else {
				currently_in_text_line = 1;
			}
			memset(csi->bfr_read, '\0', FILE_READ_BUFFER_SIZE);
			continue;
		}

		previous_line_finished = 0;
		int32_t i = 0;
		while(i < FILE_READ_BUFFER_SIZE){
			if(csi->bfr_read[i] == '\0'){
				// if(i < FILE_READ_BUFFER_SIZE - 1 && csi->bfr_read[i + 1] != '\0'){i++; continue;}
				break;
			}

			int32_t unicode_length = 1;
			if(((unsigned char) csi->bfr_read[i]) >= 0b00000000 && ((unsigned char) csi->bfr_read[i]) < 0b11000000){unicode_length = 1;}
			else if(((unsigned char) csi->bfr_read[i]) >= 0b11000000 && ((unsigned char) csi->bfr_read[i]) < 0b11100000){unicode_length = 2;}
			else if(((unsigned char) csi->bfr_read[i]) >= 0b11100000 && ((unsigned char) csi->bfr_read[i]) < 0b11110000){unicode_length = 3;}
			else {unicode_length = 4;}

			if(unicode_length == 1){
				if(csi->bfr_read[i] == '\t'){
					current_column++;
					if(current_column >= sizeof(column_sizes) / sizeof(column_sizes[0])){
						printf("Failed to parse a token correctly (sent_id=%s)\n", csi->current_sentence.sentence_id);
						current_column = 0;
						current_index_in_column = 0;
						previous_line_finished = 1; // ?
						break;
					}
					current_index_in_column = 0;
				} else if(csi->bfr_read[i] == '\r' || csi->bfr_read[i] == '\n'){
					struct token t;

					/*
					if(current_column != sizeof(column_sizes) / sizeof(column_sizes[0])){
						perror("Failed to parse a token\n");
						i += unicode_length;
					}
					*/

					if(current_column >= 9){
						err = create_token(&t, id_raw, form, lemma, upos, xpos, feats, head, deprel, deps, misc, mwe);
						if(err != 0){
							perror("failed to create token in parsing\n");
							return 1;
						}
						err = add_token_to_sentence(&(csi->current_sentence), &t);
						if(err != 0){
							perror("failed to add token to sentence in parsing\n");
							return 1;
						}
					} else {
						perror("Failed to parse a token\n");
					}					

					memset(id_raw, '\0', TOKEN_ID_RAW_SIZE);
					memset(form, '\0', TOKEN_FORM_SIZE);
					memset(lemma, '\0', TOKEN_LEMMA_SIZE);
					memset(upos, '\0', TOKEN_UPOS_SIZE);
					memset(xpos, '\0', TOKEN_XPOS_SIZE);
					memset(feats, '\0', TOKEN_FEATS_SIZE);
					memset(head, '\0', TOKEN_HEAD_SIZE);
					memset(deprel, '\0', TOKEN_DEPREL_SIZE);
					memset(deps, '\0', TOKEN_DEPS_SIZE);
					memset(misc, '\0', TOKEN_MISC_SIZE);
					memset(mwe, '\0', TOKEN_MWE_SIZE);
	
					current_column = 0;
					current_index_in_column = 0;
					previous_line_finished = 1;
				} else {
					if(current_index_in_column < column_sizes[current_column] - 1){
						columns[current_column][current_index_in_column] = csi->bfr_read[i];
					}
					current_index_in_column++;
				}
			} else {
				if(current_index_in_column + (unicode_length - 1) < column_sizes[current_column] - 1){
					for(int32_t j = 0 ; j < unicode_length ; j++){
						columns[current_column][current_index_in_column + j] = csi->bfr_read[i + j];
					}
				}
				current_index_in_column += unicode_length;
				
			}

			i += unicode_length;
		}
		memset(csi->bfr_read, '\0', FILE_READ_BUFFER_SIZE);
	}

	csi->sentence_to_free = 1;

	if(feof(csi->file_ptr)){
		csi->file_is_done = 1;
	}

	return 0;	
}

#endif
