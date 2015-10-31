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
  // Parsing machinery.
  bool valid;
  char error[1024];
  char filename[FILENAME_MAX];
  int line;

  // Validated text (block name -> arrays of text).
  string_ptr_unordered_map_t* text;
};

static string_array_t* fixed_length_block(const char* block_name, 
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
  if (!file->valid)
  {
    string_array_free(lines);
    lines = NULL;
  }
  return lines;
}

static string_array_t* variable_length_block(const char* block_name, 
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
    return NULL;
  }
  else 
    return lines;
}

static void add_text_if_valid(bool valid,
                              const char* block_name, 
                              string_array_t* lines, 
                              tough2_file_t* file)
{
  if (valid)
  {
    string_ptr_unordered_map_insert_with_kv_dtors(file->text, string_dup(block_name),
                                                  lines, string_free, DTOR(string_array_free));
  }
  else
    string_array_free(lines);
}

//------------------------------------------------------------------------
// These functions parse/validate the text in each block.
//------------------------------------------------------------------------

static void parse_incon(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  string_array_t* lines = variable_length_block("INCON", buffer, pos, file);
}

static void parse_eleme(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  string_array_t* lines = variable_length_block("ELEME", buffer, pos, file);
}

static void parse_conne(text_buffer_t* buffer, int* pos, tough2_file_t* file)
{
  string_array_t* lines = variable_length_block("CONNE", buffer, pos, file);
}

static void parse_tough2_file(const char* filename, tough2_file_t* file)
{
  // Open the file and parse the blocks.
  text_buffer_t* buffer = text_buffer_from_file(filename);
  if (buffer == NULL)
  {
    file->valid = false;
    snprintf(file->error, FILENAME_MAX, "Input file '%s' does not exist!", filename);
    return;
  }

  // Step through each line of the input file and parse it.
  int pos = 0, line_length;
  char* line;
  while (text_buffer_next_nonempty(buffer, &pos, &line, &line_length))
  {
//    // Tokenize the next block into a single text string.
//    char block_name[line_length+1];
//    string_copy_from_raw(line, line_length, block_name);
//    char* block_text;
//    int block_type = tokenize_block(block_name, buffer, &pos, &block_text);

    // If there was a problem, bug out.
    if (!file->valid)
      break;

    ++(file->line);
  }
  
  text_buffer_free(buffer);
}

tough2_file_t* tough2_file_new(const char* filename)
{
  tough2_file_t* file = polymec_malloc(sizeof(tough2_file_t));
  snprintf(file->filename, FILENAME_MAX, filename);
  file->line = 1;
  file->valid = true;
  file->error[0] = '\0';
  file->text = string_ptr_unordered_map_new();

  parse_tough2_file(filename, file);

  if (!file->valid)
  {
    tough2_file_close(file);
    polymec_error("%s is not a valid TOUGH2 input file.", filename);
  }

  return file;
}

void tough2_file_close(tough2_file_t* file)
{
  string_ptr_unordered_map_free(file->text);
  polymec_free(file);
}

point_cloud_t* tough2_file_read_mesh(tough2_file_t* file)
{
}

bool tough2_file_contains_mesh(tough2_file_t* file)
{
  return (string_ptr_unordered_map_contains(file->text, "eleme") && 
          string_ptr_unordered_map_contains(file->text, "conne"));
}

real_t* tough2_file_read_incon(tough2_file_t* file, point_cloud_t* mesh)
{
  ASSERT(mesh != NULL);
}

bool tough2_file_contains_incon(tough2_file_t* file)
{
  return string_ptr_unordered_map_contains(file->text, "incon");
}

