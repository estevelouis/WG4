#ifndef UDPIPE_INTERFACE_CONVERSION_H
#define UDPIPE_INTERFACE_CONVERSION_H

#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "udpipe_interface/cinterface.h"

int32_t transform_to_pipes(
	// struct udpipe_pipeline* const local_pipeline,
	const char* const raw_text,
	FILE** pointer_to_output_pipe,
	char** pointer_to_heap_char //,
	// int32_t output_pipefd[1][2]
){
	udpipe_pipeline_process(raw_text, pointer_to_output_pipe, (void**) pointer_to_heap_char);

	return 0;

	/*
	pipe_call_failed:
	fprintf(stderr, "Call to create pipe failed\n");
	return 1;

	closing_fd_failed:
	perror("Closing output pipefd failed\n");
	return 1;
	*/
}

#endif
