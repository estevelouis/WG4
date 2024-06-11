#ifndef TEST_GENERAL_H
#define TEST_GENERAL_H

#include <stdint.h>

#include "logging.h"

struct test_case {
	// const int32_t function(void);
	int32_t (*function)(void);
	int8_t can_fail;
};

#endif
