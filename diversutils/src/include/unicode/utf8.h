#ifndef UTF8_H
#define UTF8_H

#include <stdint.h>

#include "unicode/unicode.h"

int8_t utf8_get_length(const char *const s);

int32_t utf8_is_well_formed(const char *const s, int8_t *const len_p);

unicode_t utf8_to_unicode(const char *const s, int8_t *const len_p);

int32_t utf8_to_unicode_subset(const char *const s, int8_t *const len_p);

char *utf8_normalise(char *const s);

#endif
