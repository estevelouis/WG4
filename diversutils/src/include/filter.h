#ifndef FILTER_H
#define FILTER_H

#ifndef ENABLE_FILTER
#define ENABLE_FILTER 0
#endif

#if ENABLE_FILTER == 1

#ifndef PCRE2_CODE_UNIT_WIDTH
#define PCRE2_CODE_UNIT_WIDTH 8
#endif
#include <pcre2.h>

#include <stdint.h>

struct filter {
  pcre2_code *regex;
  PCRE2_SPTR pattern_match;
  PCRE2_SIZE length;
  const uint32_t options_compile;
  PCRE2_SPTR replacement;
  PCRE2_SIZE rlength;
  const uint32_t options_replace;
};

int32_t filter_compile(struct filter *const f);
void filter_free(struct filter *const f);
int32_t filter_ready(void);
void filter_release(void);
int32_t filter_substitute_all(char **const s, size_t *const s_len);

#endif

#endif
