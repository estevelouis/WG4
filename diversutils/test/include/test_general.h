#ifndef TEST_GENERAL_H
#define TEST_GENERAL_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

static int32_t num_calls_info = 0;
static int32_t num_calls_warning = 0;
static int32_t num_calls_error = 0;

struct test_case {
	// const int32_t function(void);
	int32_t (*function)(void);
	int8_t can_fail;
};

void custom_format(const char* const log_type, const int32_t log_num, const char* const filename, const char* const func, const int32_t line, const char* const msg){
	const int32_t isodate_capacity = 64;
	char isodate[isodate_capacity];
	memset(isodate, '\0', isodate_capacity);
	time_t now;
	struct tm* timeinfo;
	time(&now);
	timeinfo = localtime(&now);

	strftime(isodate, isodate_capacity, "%Y-%m-%dT%H:%M:%S%z", timeinfo);

	fprintf(
		stderr,
		"%s / [%s %i\x1b[0m] in function <\x1b[1m\x1b[34m%s\x1b[0m@\x1b[1m\x1b[33m%s:%i\x1b[0m>: %s\n",
		isodate,
		log_type,
		log_num,
		// num_calls,
		func,
		filename,
		line,
		msg
	);
}

void error_format(const char* const filename, const char* const func, const int32_t line, const char* const msg){
	num_calls_error++;

	custom_format("\x1b[31m\x1b[1mERROR  ", num_calls_error, filename, func, line, msg);
}

void warning_format(const char* const filename, const char* const func, const int32_t line, const char* const msg){
	num_calls_warning++;

	custom_format("\x1b[33m\x1b[1mWARNING", num_calls_warning, filename, func, line, msg);
}

void info_format(const char* const filename, const char* const func, const int32_t line, const char* const msg){
	num_calls_info++;

	custom_format("\x1b[32m\x1b[1mINFO   ", num_calls_info, filename, func, line, msg);
}

#endif
