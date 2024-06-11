#ifndef JSONL_PARSER_H
#define JSONL_PARSER_H

#include <stdint.h>
#include <stdio.h>

#include "jsonl/constants.h"
#if TOKENIZATION_METHOD == 0
#include <regex.h>
#endif

struct document {
	int32_t identifier_size;
	int32_t text_size;
	int32_t identifier_capacity;
	int32_t text_capacity;
	#if TOKENIZATION_METHOD == 0
	regex_t reg;
	regoff_t latest_rm_eo;
	#elif TOKENIZATION_METHOD == 1
	FILE* tmp_udpipe_output_file;
	#elif TOKENIZATION_METHOD == 2
	FILE* tmp_udpipe_output_file;
	// int32_t output_pipefd[2];
	char* heap_char_output;
	#endif
	int8_t reached_last_token;
	int8_t usable;
	char current_token[JSONL_CURRENT_TOKEN_BUFFER_SIZE];
	char* identifier;
	char* text;
};

struct jsonl_document_iterator {
	FILE* file_ptr;
	int8_t file_is_open;
	int8_t file_is_done;
	char bfr_read[JSONL_FILE_READ_BUFFER_SIZE];
	char content_key[JSONL_CONTENT_KEY_BUFFER_SIZE];
	struct document current_document;
	int32_t document_to_free;
    // int64_t offset_start; // included
    // int64_t offset_end; // excluded
};

void jsonl_init_tokenization();

#if (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
int32_t launch_udpipe(struct document* const doc);
#endif

int32_t create_document(struct document* doc);
int32_t iterate_document_current_token(struct document* const doc);
void free_document(struct document* doc);
int32_t create_jsonl_document_iterator(struct jsonl_document_iterator* jdi, const char* file_name, const char* const content_key);
void free_jsonl_document_iterator(struct jsonl_document_iterator* jdi);
int32_t iterate_jsonl_document_iterator(struct jsonl_document_iterator* restrict const jdi);

#endif
