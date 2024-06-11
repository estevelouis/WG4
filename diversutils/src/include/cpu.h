#ifndef CPU_H
#define CPU_H

#include <stdint.h>

struct cpu_info {
  uint16_t cardinality_virtual_cores;
  uint8_t avx256_capable;
  uint8_t avx512_capable;
};

int32_t get_cpu_info(struct cpu_info *const);

#endif
