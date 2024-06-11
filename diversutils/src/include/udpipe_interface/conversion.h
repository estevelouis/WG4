#ifndef UDPIPE_INTERFACE_CONVERSION_H
#define UDPIPE_INTERFACE_CONVERSION_H

#include <stdint.h>
#include <stdio.h>

int32_t transform_to_pipes(const char *const raw_text, FILE **pointer_to_output_pipe, char **pointer_to_heap_char);

#endif
