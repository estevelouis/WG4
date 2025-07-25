#ifndef LFSR_H
#define LFSR_H

#include<stdint.h>

struct lfsr {
    uint64_t status;
    // __uint128 status;
    uint64_t output;
};

int32_t lfsr_init(struct lfsr * restrict const rng, const uint64_t seed);
void lfsr_rand(struct lfsr * restrict const rng);

#endif
