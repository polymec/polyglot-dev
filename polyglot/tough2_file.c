// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/text_buffer.h"
#include "core/unordered_map.h"
#include "polyglot/tough2_file.h"

struct tough2_file_t 
{
  // The "title" of the problem.
  char title[81];

  // A mapping of names of blocks to their text content.
  string_ptr_unordered_map_t* blocks;

  // Some useful metadata.
  int nk, nph;

  // Parsing machinery.
  bool valid;
  char error[1024];
};

static void parse_fixed_length_block(const char* block_name, 
                                     text_buffer_t* buffer, 
                                     int* pos, 
                                     int num_lines, 
                                     tough2_file_t* file)
{
  string_array_t* lines = string_array_new();
  for (int i = 0; i < num_lines; ++i)
  {
    char* line;
    int line_length;
    bool found = text_buffer_next(buffer, pos, &line, &line_length);
    if (!found)
    {
      file->valid = false;
      snprintf(file->error, 1023, 
               "File terminated unexpectedly in %s block.", block_name);
      break;
    }
    string_array_append_with_dtor(lines, string_ndup(line, line_length), string_free);
  }
  if (file->valid)
  {
    string_array_t** text_p = (string_array_t**)string_ptr_unordered_map_get(file->blocks, (char*)block_name);
    if (text_p == NULL)
      string_ptr_unordered_map_insert_with_kv_dtors(file->blocks, (char*)block_name, lines, string_free, DTOR(string_array_free));
    else
    {
      // Transfer our lines to the existing entry.
      for (int i = 0; i < lines->size; ++i)
        string_array_append_with_dtor(*text_p, lines->data[i], lines->dtors[i]);
      string_array_release_data_and_free(lines);
    }
  }
  else
    string_array_free(lines);
}

static void parse_variable_length_block(const char* block_name, 
                                        text_buffer_t* buffer, 
                                        int* pos, 
                                        tough2_file_t* file)
{
  string_array_t* lines = string_array_new();
  char* line;
  int line_length;
  bool file_ended = true;
  while (text_buffer_next(buffer, pos, &line, &line_length))
  {
    file_ended = false;
    if ((line_length == 0) || (strncmp(line, "+++", 3) == 0))
      break;
    string_array_append_with_dtor(lines, string_ndup(line, line_length), string_free);
  }
  if (file_ended)
  {
    file->valid = false;
    snprintf(file->error, 1023, "File terminated unexpectedly in %s block.", block_name);
    string_array_free(lines);
  }
  else if (string_array_empty(lines))
  {
    file->valid = false;
    snprintf(file->error, 1023, "%s block is empty.", block_name);
    string_array_free(lines);
  }
  else 
  {
    string_array_t** text_p = (string_array_t**)string_ptr_unordered_map_get(file->blocks, (char*)block_name);
    if (text_p == NULL)
      string_ptr_unordered_map_insert_with_kv_dtors(file->blocks, string_dup(block_name), lines, string_free, DTOR(string_array_free));
    else
    {
      // Transfer our lines to the existing entry.
      for (int i = 0; i < lines->size; ++i)
        string_array_append_with_dtor(*text_p, lines->data[i], lines->dtors[i]);
      string_array_release_data_and_free(lines);
    }
  }
}

static void parse_rocks(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("ROCKS", buffer, pos, file);
}

static void parse_multi(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_fixed_length_block("MULTI", buffer, pos, 1, file);
}

static void parse_param(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_fixed_length_block("PARAM", buffer, pos, 4, file);
}

static void parse_indom(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("INDOM", buffer, pos, file);
}

static void parse_incon(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("INCON", buffer, pos, file);
}

static void parse_foft(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("FOFT", buffer, pos, file);
}

static void parse_coft(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("COFT", buffer, pos, file);
}

static void parse_goft(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("GOFT", buffer, pos, file);
}

static void parse_diffu(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  if (file->nk == -1)
  {
    file->valid = false;
    snprintf(file->error, 1023, "DIFFU block precedes MULTI block.");
  }
  else
    parse_fixed_length_block("DIFFU", buffer, pos, file->nk, file);
}

static void parse_selec(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  char* line;
  int line_length;
  bool got_header = text_buffer_next(buffer, pos, &line, &line_length);
  if (!got_header)
  {
    file->valid = false;
    snprintf(file->error, 1023, "File ended unexpectedly in SELEC block.");
  }
  else
  {
    int ie;
    int num_matches = sscanf(line, "%5d", &ie);
    if (num_matches < 1)
    {
      file->valid = false;
      snprintf(file->error, 1023, "Invalid SELEC.1 record.");
    }
    else
    {
      string_array_t* lines = string_array_new();
      string_array_append_with_dtor(lines, string_ndup(line, line_length), string_free);
      string_ptr_unordered_map_insert_with_kv_dtors(file->blocks, "SELEC", lines, string_free, DTOR(string_array_free));
      parse_fixed_length_block("SELEC", buffer, pos, ie, file);
    }
  }
}

static void parse_rpcap(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_fixed_length_block("RPCAP", buffer, pos, 2, file);
}

static void parse_times(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  char* line;
  int line_length;
  bool got_header = text_buffer_next(buffer, pos, &line, &line_length);
  if (!got_header)
  {
    file->valid = false;
    snprintf(file->error, 1023, "File ended unexpectedly in TIMES block.");
  }
  else
  {
    int iti, ite;
    double delaf, tinter;
    int num_matches = sscanf(line, "%5d%5d%10lg%10lg", &iti, &ite, &delaf, &tinter);
    if (num_matches < 4)
    {
      file->valid = false;
      snprintf(file->error, 1023, "Invalid TIMES.1 record.");
    }
    else
    {
      string_array_t* lines = string_array_new();
      string_array_append_with_dtor(lines, string_ndup(line, line_length), string_free);
      string_ptr_unordered_map_insert_with_kv_dtors(file->blocks, "TIMES", lines, string_free, DTOR(string_array_free));
      parse_fixed_length_block("TIMES", buffer, pos, iti, file);
    }
  }
}

static void parse_gener(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("GENER", buffer, pos, file);
}

static void parse_eleme(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("ELEME", buffer, pos, file);
}

static void parse_conne(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("CONNE", buffer, pos, file);
}

static void parse_meshm(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  parse_variable_length_block("MESHM", buffer, pos, file);
}

tough2_file_t* tough2_file_new(const char* filename)
{
  tough2_file_t* file = polymec_malloc(sizeof(tough2_file_t));
  file->blocks = string_ptr_unordered_map_new();
  file->title[0] = '\0';
  file->valid = true;
  file->error[0] = '\0';

  // Open the file and parse the blocks.
  text_buffer_t* buffer = text_buffer_from_file(filename);
  int pos = 0, line_length;
  char* line;
  while (text_buffer_next_nonempty(buffer, &pos, &line, &line_length))
  {
    // Set the title if it's not set yet.
    if (strlen(file->title) == 0)
    {
      strncpy(file->title, line, 80);
      continue;
    }

    // Interpret the next block.
    char block[line_length+1];
    string_copy_from_raw(line, line_length, block);
    if (strcasecmp(block, "rocks"))
      parse_rocks(buffer, &pos, file);
    else if (strcasecmp(block, "multi"))
      parse_multi(buffer, &pos, file);
    else if (strcasecmp(block, "param"))
      parse_param(buffer, &pos, file);
    else if (strcasecmp(block, "indom"))
      parse_indom(buffer, &pos, file);
    else if (strcasecmp(block, "incon"))
      parse_incon(buffer, &pos, file);
    else if (strcasecmp(block, "foft"))
      parse_foft(buffer, &pos, file);
    else if (strcasecmp(block, "coft"))
      parse_coft(buffer, &pos, file);
    else if (strcasecmp(block, "goft"))
      parse_goft(buffer, &pos, file);
    else if (strcasecmp(block, "diffu"))
      parse_diffu(buffer, &pos, file);
    else if (strcasecmp(block, "selec"))
      parse_selec(buffer, &pos, file);
    else if (strcasecmp(block, "rpcap"))
      parse_rpcap(buffer, &pos, file);
    else if (strcasecmp(block, "times"))
      parse_times(buffer, &pos, file);
    else if (strcasecmp(block, "gener"))
      parse_gener(buffer, &pos, file);
    else if (strcasecmp(block, "eleme"))
      parse_eleme(buffer, &pos, file);
    else if (strcasecmp(block, "conne"))
      parse_conne(buffer, &pos, file);
    else if (strcasecmp(block, "meshm"))
      parse_meshm(buffer, &pos, file);
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
  string_ptr_unordered_map_free(file->blocks);
  polymec_free(file);
}

bool tough2_file_contains_block(tough2_file_t* file, 
                                const char* block_name)
{
  return string_ptr_unordered_map_contains(file->blocks, (char*)block_name);
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

