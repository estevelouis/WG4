#ifndef CFGPARSER_PARSER_H
#define CFGPARSER_PARSER_H

#include <stdint.h>

#define CFGPARSER_BFR_SIZE_READ 1024

struct cfg {
  char **keys;
  char **values;
  uint32_t num_entries;
};

int32_t create_cfg_from_file(struct cfg *const, const char *const);

void free_cfg(struct cfg *const);

char *cfg_get_value(const struct cfg *const, const char *const);

#endif
