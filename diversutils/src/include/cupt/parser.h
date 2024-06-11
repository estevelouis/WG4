#ifndef CUPT_PARSER_H
#define CUPT_PARSER_H

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

struct sentence {
  struct token *tokens;
  int32_t capacity;
  int32_t num_tokens;
  char sentence_id[SENTENCE_SENT_ID_SIZE];
};

struct cupt_sentence_iterator {
  FILE *file_ptr;
  int8_t file_is_open;
  int8_t file_is_done;
  char bfr_read[FILE_READ_BUFFER_SIZE];
  struct sentence current_sentence;
  int32_t sentence_to_free;
};

int32_t create_token(struct token *const, char *const, char *const, char *const, char *const, char *const, char *const,
                     char *const, char *const, char *const, char *const, char *const);

int32_t serialize_token(struct token *const, char *const, int32_t *const, const int8_t);

int32_t create_sentence(struct sentence *const, const char *const);

void free_sentence(struct sentence *const);

int32_t add_token_to_sentence(struct sentence *const, struct token *const);

size_t buffer_size_for_serialization_of_sentence(struct sentence *const);

int32_t serialize_sentence(struct sentence *const, char *const, int32_t *const);

int32_t create_cupt_sentence_iterator(struct cupt_sentence_iterator *const, const char *const);

void free_cupt_sentence_iterator(struct cupt_sentence_iterator *const);

int32_t iterate_cupt_sentence_iterator(struct cupt_sentence_iterator *const restrict);

#endif
