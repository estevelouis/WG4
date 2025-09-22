#include <fcntl.h>
#include <unistd.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "cupt/parser.h"

#if (TOKENIZATION_METHOD == 1 || TOKENIZATION_METHOD == 2)
#include "udpipe/interface/cinterface.h" // external declarations
#include "udpipe/interface/conversion.h" // transformation to pipes
#endif

void write_check(int fd, char * str, size_t n){
  ssize_t res = write(fd, str, n);
  if(res == -1){
    perror("Failed to write: return status is -1");
    exit(1);
  }
}

int32_t cupt_ec_process_sentence(const csi_t* csi, char ** const tsv, const size_t bfr_size_tsv){
  static uint64_t call_count = 0;

  printf("starting cupt_ec_process_sentence\n");

  sent_t* const sent = &(csi->current_sentence);

  // generate random identifier

  const size_t bfr_size_rand_iter = 64;
  char rand_iter[bfr_size_rand_iter + 1];
  srand(time(NULL));
  for(size_t i = 0 ; i < bfr_size_rand_iter ; i++){
    int rand_gen = rand();
    rand_iter[i] = '0' + (rand_gen % 10);
  }
  rand_iter[bfr_size_rand_iter] = '\0';

  // prepare paths for io

  const size_t bfr_size_path = 512;
  char path_in[bfr_size_path];
  memset(path_in, '\0', bfr_size_path);
  snprintf(path_in, bfr_size_path, "/lustre/fswork/projects/rech/lpe/uki74yt/diverscontrol-tmp/cupt_ec_process_sentence_%s-input.cupt", rand_iter);
  char path_out[bfr_size_path];
  memset(path_out, '\0', bfr_size_path);
  snprintf(path_out, bfr_size_path, "/lustre/fswork/projects/rech/lpe/uki74yt/diverscontrol-tmp/cupt_ec_process_sentence_%s-output.cupt", rand_iter);

  snprintf(*tsv, bfr_size_tsv, "%s.tsv", path_out);

  printf("%s\n%s\n%s\n", path_in, path_out, *tsv);

  /**/
  // write sent_id
  
//  /*
//  const size_t bfr_size_sent_id = 8192;
//  char sent_id[bfr_size_sent_id];
//  memset(sent_id, '\0', bfr_size_sent_id);
//  snprintf(sent_id, bfr_size_sent_id, "# sent_id = %s\n", sent->sentence_id);
//  */
//
//  int fd_in = open(path_in, O_CREAT | O_WRONLY | O_TRUNC, 0600);
//  // write_check(fd_in, sent_id, strlen(sent_id));
//
//  const int32_t n = sent->num_tokens;
//
//  for(int32_t i = 0 ; i < n ; i++){
//    write_check(fd_in, sent->tokens[i].id_raw, strlen(sent->tokens[i].id_raw));
//    write_check(fd_in, " ", 1);
//    /*
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].form, strlen(sent->tokens[i].form));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].lemma, strlen(sent->tokens[i].lemma));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].upos, strlen(sent->tokens[i].upos));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].xpos, strlen(sent->tokens[i].xpos));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].feats, strlen(sent->tokens[i].feats));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].head, strlen(sent->tokens[i].head));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].deprel, strlen(sent->tokens[i].deprel));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].deps, strlen(sent->tokens[i].deps));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].misc, strlen(sent->tokens[i].misc));
//    write_check(fd_in, "\n", 1);
//    write_check(fd_in, sent->tokens[i].mwe, strlen(sent->tokens[i].mwe));
//    write_check(fd_in, "\n\n", 2);
//    */
//  }
//  /**/

  FILE * f_in = fopen(path_in, "w");
  if(f_in == NULL){perror("failed to call fopen for f_in\n"); exit(1);}

  printf("regenerating readable text\n");

  const size_t bfr_write_size = 65536;
  char * bfr_write = malloc(bfr_write_size);
  if(bfr_write == NULL){perror("malloc failed\n"); return 1;}
  memset(bfr_write, '\0', bfr_write_size);
  printf("malloced bfr_write\n");

  int32_t j = 0;
  int32_t n = csi->current_sentence.num_tokens;
  for(int32_t i = 0 ; i < n ; i++){
    size_t available_bfr_size = (bfr_write_size - 1) - j;
    if(available_bfr_size < 0){break;}
    snprintf(bfr_write + j, available_bfr_size, "%s ", sent->tokens[i].form);
    j += strlen(sent->tokens[i].form) + 1;
  }
  // write_check(fd_in, bfr_write, j - 1); // for verification purposes, not required
  fputs(bfr_write, f_in);

  char * bfr_out = malloc(bfr_write_size);
  if(bfr_out == NULL){perror("malloc failed\n"); return 1;}
  memset(bfr_out, '\0', bfr_write_size);

  FILE * f_out = fopen(path_out, "w");
  printf("calling transform_to_pipes\n");
  // if(transform_to_pipes(bfr_write, &f_out, NULL) != 0){
  // if(transform_to_pipes(bfr_write, &f_out, &bfr_out) != 0){
  if(transform_to_pipes(bfr_write, NULL, &bfr_out) != 0){
    perror("failed to call transform_to_pipes\n");
    exit(1);
  }
  printf("udpipe output: %s\n", bfr_out);
  fputs(bfr_out, f_out);
  fclose(f_out);

  free(bfr_write);
  free(bfr_out);


  /*
  const size_t bfr_write_size = 65536;
  char * bfr_write = malloc(bfr_write_size);
  if(bfr_write == NULL){perror("malloc failed\n"); return 1;}
  memset(bfr_write, '\0', bfr_write_size);
  int32_t bytes_written = 0;
  if(serialize_sentence(&(csi->current_sentence), bfr_write, &bytes_written) != 0){
    perror("failed to call serialize_sentence\n");
    close(fd_in);
    return 1;
  }
  write_check(fd_in, bfr_write, bytes_written);
  free(bfr_write);
  */

  // close(fd_in);
  fclose(f_in);







  const size_t bfr_size_cmd = 512;
  char cmd[bfr_size_cmd];
  memset(cmd, '\0', bfr_size_cmd);
  snprintf(cmd, bfr_size_cmd, "python3 ${THESIS}/other_repos/STARK/stark.py --config=%s --input=%s --output=%s", csi->ec_cfg, path_out, *tsv);

  printf("command: %s\n", cmd);

  system(cmd);

  call_count++;

  return 0;
}
