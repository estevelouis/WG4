#ifndef LOGGING_H
#define LOGGING_H

#include <stdint.h>

void custom_format(const char* const log_type, const int32_t log_num, const char* const filename, const char* const func, const int32_t line, const char* const msg);
void error_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);
void warning_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);
void info_format(const char* const filename, const char* const func, const int32_t line, const char* const msg);

#endif
