#ifndef UDPIPE_INTERFACE_CINTERFACE_H
#define UDPIPE_INTERFACE_CINTERFACE_H

#include "udpipe_interface/size.h"
#include <stdint.h>

struct udpipe_pipeline {
  uint8_t placeholder[UDPIPE_PIPELINE_SIZE];
};

extern void ensure_proper_udpipe_pipeline_size();

/*
extern void udpipe_pipeline_create(const char*, const char*, const char*, const char*, const char*, struct udpipe_pipeline*);
*/

extern void udpipe_pipeline_create_global(const char *, const char *, const char *, const char *, const char *);

/*
extern void udpipe_pipeline_print_global_info();
*/

extern int32_t udpipe_pipeline_process(const char *, FILE **const, void **const); // in

#endif
