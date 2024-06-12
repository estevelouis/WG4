/*
 *      DiversUtils - Functions to measure diversity
 *
 * Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int32_t num_calls_info = 0;
static int32_t num_calls_warning = 0;
static int32_t num_calls_error = 0;

void custom_format(const char *const log_type, const int32_t log_num, const char *const filename, const char *const func,
                   const int32_t line, const char *const msg) {
  const int32_t isodate_capacity = 64;
  char isodate[isodate_capacity];
  memset(isodate, '\0', isodate_capacity);
  time_t now;
  struct tm *timeinfo;
  time(&now);
  timeinfo = localtime(&now);

  strftime(isodate, isodate_capacity, "%Y-%m-%dT%H:%M:%S%z", timeinfo);

  fprintf(stderr, "%s [%s %i\x1b[0m] <\x1b[1m\x1b[34m%s\x1b[0m@\x1b[1m\x1b[33m%s:%i\x1b[0m> %s\n", isodate, log_type, log_num,
          // num_calls,
          func, filename, line, msg);
}

void error_format(const char *const filename, const char *const func, const int32_t line, const char *const msg) {
  num_calls_error++;

  custom_format("\x1b[31m\x1b[1mERROR  ", num_calls_error, filename, func, line, msg);
}

void warning_format(const char *const filename, const char *const func, const int32_t line, const char *const msg) {
  num_calls_warning++;

  custom_format("\x1b[33m\x1b[1mWARNING", num_calls_warning, filename, func, line, msg);
}

void info_format(const char *const filename, const char *const func, const int32_t line, const char *const msg) {
  num_calls_info++;

  custom_format("\x1b[32m\x1b[1mINFO   ", num_calls_info, filename, func, line, msg);
}
