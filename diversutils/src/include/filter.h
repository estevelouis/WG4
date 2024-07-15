#ifndef FILTER_H
#define FILTER_H

/**/
#ifndef PCRE2_CODE_UNIT_WIDTH
#define PCRE2_CODE_UNIT_WIDTH 8
#endif
#include <pcre2.h>
/**/

/*
#include <regex.h>
*/

#include <stdint.h>

#ifndef ENABLE_FILTER
#define ENABLE_FILTER 0
#endif

struct filter {
  /*
  const char * const pattern_match;
  const char * const pattern_substitution;
  const int32_t flags;
  regex_t regex;
  */
  pcre2_code *regex;
  PCRE2_SPTR pattern_match;
  PCRE2_SIZE length;
  const uint32_t options_compile;
  PCRE2_SPTR replacement;
  PCRE2_SIZE rlength;
  const uint32_t options_replace;
  /*
  int * errorcode;
  PCRE2_SIZE * erroroffset;
  pcre2_compile_context * context;
  */
};

int32_t filter_compile(struct filter *const f);
void filter_free(struct filter *const f);
int32_t filter_ready(void);
void filter_release(void);
int32_t filter_substitute_all(char **const s, size_t *const s_len);

#endif
