#include <stdint.h>
#include <stdio.h>
// #include <string.h>
// #include <unistd.h>

#include "udpipe_interface/cinterface.h"

int32_t transform_to_pipes(const char *const raw_text, FILE **pointer_to_output_pipe, char **pointer_to_heap_char) {
  udpipe_pipeline_process(raw_text, pointer_to_output_pipe, (void **)pointer_to_heap_char);

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
