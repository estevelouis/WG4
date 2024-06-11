#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cpu.h"

int32_t get_cpu_info(struct cpu_info *const result) {
  FILE *f;
  const int32_t bfr_size = 2048;
  char bfr[bfr_size];
  const char processor_key[] = "processor\t: ";
  const size_t processor_key_length = strlen(processor_key);
  const char flag_key[] = "flags\t\t: ";
  const size_t flag_key_length = strlen(flag_key);
  const char avx256_key[] = "avx2";
  const size_t avx256_key_length = strlen(avx256_key);
  const char avx512_key[] = "avx512";
  const size_t avx512_key_length = strlen(avx512_key);
  int64_t max_proc, current_proc;
  uint8_t avx256_capable, avx512_capable;
  char *strtok_placeholder;

  f = fopen("/proc/cpuinfo", "r");
  if (f == NULL) {
    perror("Failed to open \"/proc/cpuinfo\"\n");
    return 1;
  }

  memset(bfr, '\0', bfr_size);

  max_proc = -1;
  avx256_capable = 0;
  avx512_capable = 0;
  while (fgets(bfr, bfr_size, f)) {
    if (strncmp(processor_key, bfr, processor_key_length) == 0) {
      current_proc = strtol(&(bfr[processor_key_length]), NULL, 10);
      if (current_proc > max_proc) {
        max_proc = current_proc;
      }
    } else if (strncmp(flag_key, bfr, flag_key_length) == 0) {
      strtok_placeholder = strtok(&(bfr[flag_key_length]), " ");
      while ((avx256_capable == 0 || avx512_capable == 0) && strtok_placeholder != NULL) {
        // printf("current_token: %s\n", strtok_placeholder);
        if (strncmp(avx256_key, strtok_placeholder, avx256_key_length) == 0) {
          avx256_capable = 1;
        } else if (strncmp(avx512_key, strtok_placeholder, avx512_key_length) == 0) {
          avx512_capable = 1;
        }
        strtok_placeholder = strtok(NULL, " ");
      }
    }
  }

  fclose(f);

  if (max_proc == -1) {
    perror("Could not find number of virtual cores in \"/proc/cpuinfo\"\n");
    return 1;
  }

  (*result) = (struct cpu_info){
      .cardinality_virtual_cores = (uint16_t)(max_proc + 1),
      .avx256_capable = avx256_capable,
      .avx512_capable = avx512_capable,
  };

  return 0;
}
