// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/text_buffer.h"
#include "core/unordered_map.h"
#include "polyglot/tough2_file.h"

struct tough2_file_t 
{
  string_string_unordered_map_t* blocks;

  // Parsing machinery.
  bool valid;
  char error[1024];
};

static void parse_rocks(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_multi(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_param(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_indom(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_incon(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_foft(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_coft(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_goft(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_diffu(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_selec(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_rpcap(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_times(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_gener(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_eleme(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_conne(text_buffer_t* buffer, tough2_file_t* file)
{
}

static void parse_meshm(text_buffer_t* buffer, tough2_file_t* file)
{
}

tough2_file_t* tough2_file_new(const char* filename)
{
  tough2_file_t* file = polymec_malloc(sizeof(tough2_file_t));
  file->blocks = string_string_unordered_map_new();
  file->valid = true;
  file->error[0] = '\0';

  // Open the file and parse the blocks.
  text_buffer_t* buffer = text_buffer_from_file(filename);
  int pos = 0, line_length;
  char* line;
  while (text_buffer_next_nonempty(buffer, &pos, &line, &line_length))
  {
    char block[line_length+1];
    string_copy_from_raw(line, line_length, block);
    if (strcasecmp(block, "rocks"))
      parse_rocks(buffer, file);
    else if (strcasecmp(block, "multi"))
      parse_multi(buffer, file);
    else if (strcasecmp(block, "param"))
      parse_param(buffer, file);
    else if (strcasecmp(block, "indom"))
      parse_indom(buffer, file);
    else if (strcasecmp(block, "incon"))
      parse_incon(buffer, file);
    else if (strcasecmp(block, "foft"))
      parse_foft(buffer, file);
    else if (strcasecmp(block, "coft"))
      parse_coft(buffer, file);
    else if (strcasecmp(block, "goft"))
      parse_goft(buffer, file);
    else if (strcasecmp(block, "diffu"))
      parse_diffu(buffer, file);
    else if (strcasecmp(block, "selec"))
      parse_selec(buffer, file);
    else if (strcasecmp(block, "rpcap"))
      parse_rpcap(buffer, file);
    else if (strcasecmp(block, "times"))
      parse_times(buffer, file);
    else if (strcasecmp(block, "gener"))
      parse_gener(buffer, file);
    else if (strcasecmp(block, "eleme"))
      parse_eleme(buffer, file);
    else if (strcasecmp(block, "conne"))
      parse_conne(buffer, file);
    else if (strcasecmp(block, "meshm"))
      parse_meshm(buffer, file);
    else 
      string_string_unordered_map_insert_with_k_dtor(file->blocks, string_dup(block), NULL, string_free);
    if (!file->valid)
      break;
  }

  text_buffer_free(buffer);
  if (!file->valid)
  {
    tough2_file_close(file);
    file = NULL;
  }
  return file;
}

bool tough2_file_is_valid(tough2_file_t* file)
{
  return file->valid;
}

char* tough2_file_error(tough2_file_t* file)
{
  return &file->error[0];
}

void tough2_file_close(tough2_file_t* file)
{
  string_string_unordered_map_free(file->blocks);
  polymec_free(file);
}

bool tough2_file_contains_block(tough2_file_t* file, 
                                const char* block_name)
{
  return string_string_unordered_map_contains(file->blocks, (char*)block_name);
}

point_cloud_t* tough2_file_read_mesh(tough2_file_t* file)
{
}

bool tough2_file_contains_mesh(tough2_file_t* file)
{
  return (tough2_file_contains_block(file, "eleme") && 
          tough2_file_contains_block(file, "conne")) || 
         tough2_file_contains_block(file, "meshm");
}

tough2_file_rock_t* tough2_file_read_rocks(tough2_file_t* file, int* num_entries)
{
}

bool tough2_file_contains_rocks(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "rocks");
}

void tough2_file_read_multi(tough2_file_t* file,
                            int* nk, int* neq, int* nph, int* nb, int* nkin)
{
}

bool tough2_file_contains_multi(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "multi");
}

void tough2_file_read_param(tough2_file_t* file,
                            int* noite, int* kdata, int* mcyc, 
                            real_t* msec, int* mcypr, int* mop,
                            real_t* texp, real_t* be,
                            real_t* tstart, real_t* timax,
                            real_t* delten, real_t* deltmx,
                            char* elst, real_t* gf, real_t* redlt,
                            real_t* scale, real_t* dlt, 
                            real_t* re1, real_t* re2, real_t* u,
                            real_t* wup, real_t* wnr, real_t* dfac,
                            real_t* dep)
{
}

bool tough2_file_contains_param(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "param");
}

void tough2_file_read_indom(tough2_file_t* file,
                            char* mat,
                            real_t* X)
{
}

bool tough2_file_contains_indom(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "indom");
}

real_t* tough2_file_read_incon(tough2_file_t* file)
{
}

bool tough2_file_contains_incon(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "incon");
}

char** tough2_file_read_foft(tough2_file_t* file)
{
}

bool tough2_file_contains_foft(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "foft");
}

char** tough2_file_read_coft(tough2_file_t* file)
{
}

bool tough2_file_contains_coft(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "coft");
}

char** tough2_file_read_goft(tough2_file_t* file)
{
}

bool tough2_file_contains_goft(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "goft");
}

tough2_file_diffusion_t* tough2_file_read_diffu(tough2_file_t* file, int* num_entries)
{
}

bool tough2_file_contains_diffu(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "diffu");
}

void tough2_file_read_selec(tough2_file_t* file, real_t* selec, int* size)
{
}

bool tough2_file_contains_selec(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "selec");
}

void tough2_file_read_rpcap(tough2_file_t* file, 
                            tough2_file_rel_perm_t* rel_perm_type,
                            real_t* rel_perm_params,
                            tough2_file_cap_pres_t* cap_pres_type,
                            real_t* cap_pres_params)
{
}

bool tough2_file_contains_rpcap(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "rpcap");
}

real_t* tough2_file_read_times(tough2_file_t* file, int* num_times)
{
}

bool tough2_file_contains_times(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "times");
}

tough2_file_source_t* tough2_file_read_gener(tough2_file_t* file, int* num_sources)
{
}

bool tough2_file_contains_gener(tough2_file_t* file)
{
  return tough2_file_contains_block(file, "gener");
}

