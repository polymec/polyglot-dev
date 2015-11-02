// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/text_buffer.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "model/neighbor_pairing.h"
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
    return NULL;
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

static void parse_tough2_file(const char* filename, tough2_file_t* file)
{
  // Open the file and parse the blocks.
  text_buffer_t* buffer = text_buffer_from_file(filename);
  if (buffer == NULL)
  {
    file->valid = false;
    snprintf(file->error, 1023, "Input file '%s' does not exist!", filename);
    return;
  }

  // Step through each line of the input file and parse it.
  int pos = 0, line_length;
  char* line;
  while (text_buffer_next_nonempty(buffer, &pos, &line, &line_length))
  {
    // Find the block name.
    char block_name[line_length+1];
    string_copy_from_raw(line, line_length, block_name);

    // Parse the blocks of interest.
    if ((strcasecmp(block_name, "incon") == 0) ||
        (strcasecmp(block_name, "eleme") == 0) ||
        (strcasecmp(block_name, "conne") == 0))
    {
      string_array_t* lines = variable_length_block(block_name, buffer, &pos, file);
      if (file->valid)
      {
        string_ptr_unordered_map_insert_with_kv_dtors(file->text, string_dup(block_name),
                                                      lines, string_free, DTOR(string_array_free));
      }
      else
        string_array_free(lines);
    }

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
  strncpy(file->filename, filename, FILENAME_MAX);
  file->line = 1;
  file->valid = true;
  file->error[0] = '\0';
  file->text = string_ptr_unordered_map_new();

  parse_tough2_file(filename, file);

  if (!file->valid)
  {
    tough2_file_close(file);
    polymec_error("%s is not a valid TOUGH2 input file:\n%s", filename, file->error);
  }

  return file;
}

void tough2_file_close(tough2_file_t* file)
{
  string_ptr_unordered_map_free(file->text);
  polymec_free(file);
}

static serializer_t* t2_field_serializer(point_cloud_t* mesh)
{
  // FIXME
  return NULL;
}

static bool string_is_blank(char* str)
{
  bool blank = true;
  for (int j = 0; j < strlen(str); ++j)
  {
    if (str[j] != ' ')
    {
      blank = false;
      break;
    }
  }
  return blank;
}

static void pair_free(int* pair)
{
  polymec_free(pair);
}

point_cloud_t* tough2_file_read_mesh(tough2_file_t* file)
{
  if (!tough2_file_contains_mesh(file))
    return NULL;

  string_array_t* eleme_lines = *string_ptr_unordered_map_get(file->text, "eleme");
  int num_elements = eleme_lines->size;
  point_cloud_t* mesh = point_cloud_new(MPI_COMM_WORLD, num_elements);

  // Add element properties.
  real_t* volume = polymec_malloc(sizeof(real_t) * num_elements);
  real_t* ahtx = polymec_malloc(sizeof(real_t) * num_elements);
  real_t* pmx = polymec_malloc(sizeof(real_t) * num_elements);
  point_cloud_set_property(mesh, "volume", volume, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "ahtx", ahtx, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "pmx", pmx, t2_field_serializer(mesh));

  // Add a mapping of element names to indices.
  string_int_unordered_map_t* elem_index_map = string_int_unordered_map_new();
  point_cloud_set_property(mesh, "elem_index", elem_index_map, NULL);

  // This keeps track of miscellaneous tags we want to add.
  string_ptr_unordered_map_t* misc_tags = string_ptr_unordered_map_new();

  // Parse element data.
  log_debug("tough2_file_read_mesh: Parsing %d elements...", num_elements);
  bool have_ahtx = false, have_pmx = false;
  for (int i = 0; i < eleme_lines->size; ++i)
  {
    char* line = eleme_lines->data[i];
    int offset = 0;

    // Parse the element name.
    char el[4], ne[3];
    strncpy(el, &line[offset], 3); offset += 3;
    strncpy(ne, &line[offset], 2); offset += 2;
    if (!string_is_number((const char*)ne))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: invalid element name.", i);
      break;
    }
    char elem_name[6];
    snprintf(elem_name, 6, "%s%s", el, ne);

    // Parse the NSEQ and NADD values, and reject any non-blank entries!
    char nseq[6], nadd[6];
    strncpy(nseq, &line[offset], 5); offset += 5;
    strncpy(nadd, &line[offset], 5); offset += 5;
    if (!string_is_blank(nseq) || !string_is_blank(nadd))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: NSEQ and NADD are not allowed.", i);
      break;
    }

    // Parse the material name and tag the element accordingly.
    char ma1[4], ma2[3];
    strncpy(ma1, &line[offset], 3); offset += 3;
    strncpy(ma2, &line[offset], 2); offset += 2;
    int_array_t* tag = NULL;
    if (string_is_blank(ma1))
    {
      int_array_t** tag_p;
      if (string_is_blank(ma2))
        tag_p = (int_array_t**)string_ptr_unordered_map_get(misc_tags, "ROCK0");
      else
      {
        char ma[6];
        snprintf(ma, 5, "ROCK%s", ma2);
        tag_p = (int_array_t**)string_ptr_unordered_map_get(misc_tags, ma);
      }

      if (tag_p == NULL)
      {
        tag = int_array_new();
        string_ptr_unordered_map_insert_with_kv_dtors(misc_tags, string_dup("ROCK0"), tag, string_free, DTOR(int_array_free));
      }
      else
        tag = *tag_p;
    }
    else
    {
      int_array_t** tag_p;
      char ma[6];
      snprintf(ma, 5, "%s%s", ma1, ma2);
      tag_p = (int_array_t**)string_ptr_unordered_map_get(misc_tags, ma);
      if (tag_p == NULL)
      {
        tag = int_array_new();
        string_ptr_unordered_map_insert_with_kv_dtors(misc_tags, string_dup(ma), tag, string_free, DTOR(int_array_free));
      }
      else
        tag = *tag_p;
    }
    int_array_append(tag, i);

    // Parse the volume, interface area, permeability modifier.
    char volx_str[11], ahtx_str[11], pmx_str[11];
    strncpy(volx_str, &line[offset], 10); offset += 10;
    strncpy(ahtx_str, &line[offset], 10); offset += 10;
    strncpy(pmx_str, &line[offset], 10); offset += 10;
    if (!string_is_number(volx_str))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: invalid volume.", i);
      break;
    }
    else
      volume[i] = (real_t)atof(volx_str);

    if (string_is_number(ahtx_str))
    {
      have_ahtx = true;
      ahtx[i] = (real_t)atof(ahtx_str);
    }
    if (string_is_number(pmx_str))
    {
      have_pmx = true;
      pmx[i] = (real_t)atof(pmx_str);
    }

    // Parse the coordinates.
    char x[11], y[11], z[11];
    strncpy(x, &line[offset], 10); offset += 10;
    strncpy(y, &line[offset], 10); offset += 10;
    strncpy(z, &line[offset], 10); offset += 10;
    if (!string_is_number(x))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: invalid x coordinate.", i);
      break;
    }
    if (!string_is_number(y))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: invalid y coordinate.", i);
      break;
    }
    if (!string_is_number(z))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing element %d: invalid z coordinate.", i);
      break;
    }
    mesh->points[i].x = (real_t)atof(x);
    mesh->points[i].y = (real_t)atof(y);
    mesh->points[i].z = (real_t)atof(z);
  }
  if (!file->valid)
  {
    point_cloud_free(mesh);
    return NULL;
  }

  // Prune properties if needed.
  if (!have_ahtx)
    point_cloud_delete_property(mesh, "ahtx");
  if (!have_pmx)
    point_cloud_delete_property(mesh, "pmx");

  // Add miscellaneous element tags.
  int pos = 0;
  char* tag;
  void* data;
  while (string_ptr_unordered_map_next(misc_tags, &pos, &tag, &data))
  {
    int_array_t* elements = data;
    point_cloud_create_tag(mesh, tag, elements->size);
    for (int i = 0; i < elements->size; ++i)
      tag[i] = elements->data[i];
  }
  string_ptr_unordered_map_free(misc_tags);

  // Now parse connection data.
  string_array_t* conne_lines = *string_ptr_unordered_map_get(file->text, "conne");
  int num_connections = conne_lines->size;

  // Do a first pass through the connections and construct the index space, tossing out bad connections.
  int_unordered_set_t* bad_connections = int_unordered_set_new();
  int_pair_unordered_set_t* recorded_pairs = int_pair_unordered_set_new();
  int* pair_pool = polymec_malloc(sizeof(int) * 2 * num_connections); // Memory "pool" for pairs.
  for (int i = 0; i < num_connections; ++i)
  {
    char* line = eleme_lines->data[i];
    int offset = 0;

    // Parse the element names.
    char el1[6], el2[6];
    strncpy(el1, &line[offset], 5); offset += 5;
    strncpy(el2, &line[offset], 5); offset += 5;

    if (!string_int_unordered_map_contains(elem_index_map, el1))
    {
      int_unordered_set_insert(bad_connections, i);
      log_debug("tough2_file_read_mesh: Pruning connection %s -> %s (element %s does not exist).", el1, el2, el1);
    }
    else if (!string_int_unordered_map_contains(elem_index_map, el2))
    {
      int_unordered_set_insert(bad_connections, i);
      log_debug("tough2_file_read_mesh: Pruning connection %s -> %s (element %s does not exist).", el1, el2, el2);
    }

    // Figure out the element indices and create the pair.
    int e1 = *string_int_unordered_map_get(elem_index_map, el1);
    int e2 = *string_int_unordered_map_get(elem_index_map, el2);
    int pair1[2] = {e1, e2}, pair2[2] = {e2, e1};
    if (int_pair_unordered_set_contains(recorded_pairs, pair1) || 
        int_pair_unordered_set_contains(recorded_pairs, pair2))
    {
      int_unordered_set_insert(bad_connections, i);
      log_debug("tough2_file_read_mesh: Pruning duplicate connection %s -> %s.");
    }
    else
    {
      int* pair = &pair_pool[2*i];
      pair[0] = e1;
      pair[1] = e2;
      int_pair_unordered_set_insert(recorded_pairs, pair);
    }
  }
  if (bad_connections->size > 0)
    log_debug("tough2_file_read_mesh: Pruned %d invalid/duplicate connections.", bad_connections->size);
  num_connections -= bad_connections->size;
  int_pair_unordered_set_free(recorded_pairs);
  polymec_free(pair_pool);

  // Construct properties for connections.
  real_t* d1 = polymec_malloc(sizeof(real_t) * num_connections);
  real_t* d2 = polymec_malloc(sizeof(real_t) * num_connections);
  real_t* area = polymec_malloc(sizeof(real_t) * num_connections);
  real_t* beta = polymec_malloc(sizeof(real_t) * num_connections);
  real_t* sigma = polymec_malloc(sizeof(real_t) * num_connections);
  point_cloud_set_property(mesh, "d1", d1, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "d2", d2, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "area", area, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "beta", beta, t2_field_serializer(mesh));
  point_cloud_set_property(mesh, "sigma", sigma, t2_field_serializer(mesh));

  // Now construct the neighbor pairing.
  log_debug("tough2_file_read_mesh: Parsing %d connections...", num_connections);
  int* pairs = polymec_malloc(sizeof(int) * 2 * num_connections);
  int k = 0;
  bool have_distances = false, have_beta = false, have_sigma = false;
  for (int i = 0; i < num_connections; ++i)
  {
    if (int_unordered_set_contains(bad_connections, i)) 
      continue;

    char* line = eleme_lines->data[i];
    int offset = 0;

    // Record the pair of elements.
    char el1[6], el2[6];
    strncpy(el1, &line[offset], 5); offset += 5;
    strncpy(el2, &line[offset], 5); offset += 5;
    int e1 = *string_int_unordered_map_get(elem_index_map, el1);
    int e2 = *string_int_unordered_map_get(elem_index_map, el2);
    pairs[2*k]   = e1;
    pairs[2*k+1] = e2;

    // Parse the NSEQ, NAD1, and NAD2 values, and reject any non-blank entries!
    char nseq[6], nad1[6], nad2[6];
    strncpy(nseq, &line[offset], 5); offset += 5;
    strncpy(nad1, &line[offset], 5); offset += 5;
    strncpy(nad2, &line[offset], 5); offset += 5;
    if (!string_is_blank(nseq) || !string_is_blank(nad1) || !string_is_blank(nad2))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: NSEQ, NAD1, and NAD2 are not allowed.", k);
      break;
    }

    // Parse ISOT.
    char isot[6];
    strncpy(isot, &line[offset], 5); offset += 5;
    if (!string_is_blank(isot) && !string_is_number(isot))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: Invalid ISOT entry (must be 1, 2, or 3).", k);
      break;
    }

    // Parse D1, D2.
    char d1_str[11], d2_str[11];
    strncpy(d1_str, &line[offset], 10); offset += 10;
    strncpy(d2_str, &line[offset], 10); offset += 10;
    if (string_is_number(d1_str))
    {
      if (!string_is_number(d2_str))
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing connection %d: Invalid D2.", k);
        break;
      }
      have_distances = true;
      d1[k] = (real_t)atof(d1_str);
      d2[k] = (real_t)atof(d2_str);
    }
    else if (string_is_number(d2_str))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: Invalid D1.", k);
      break;
    }

    // Parse AREAX.
    char areax[11];
    strncpy(areax, &line[offset], 10); offset += 10;
    if (!string_is_number(areax))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: Invalid AREAX.", k);
      break;
    }
    area[k] = (real_t)atof(areax);

    // Parse BETAX.
    char betax[11];
    strncpy(betax, &line[offset], 10); offset += 10;
    if (string_is_number(betax))
    {
      have_beta = true;
      beta[k] = (real_t)atof(betax);
    }
    else if (!string_is_blank(betax))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: Invalid BETAX.", k);
      break;
    }

    // Parse SIGX.
    char sigx[11];
    strncpy(sigx, &line[offset], 10); offset += 10;
    if (string_is_number(sigx))
    {
      real_t sig = (real_t)atof(sigx);
      if (sig < 0.0)
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing connection %d: SIGX must be non-negative.", k);
        break;
      }
      have_sigma = true;
      sigma[k] = sig;
    }
    else if (!string_is_blank(sigx))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing connection %d: Invalid SIGX.", k);
      break;
    }

    ++k;
  }

  // Set up the neighbor pairing.
  neighbor_pairing_t* neighbors = unweighted_neighbor_pairing_new("connections", num_connections, pairs, 
                                                                  exchanger_new(mesh->comm));
  point_cloud_set_property(mesh, "connections", neighbors, neighbor_pairing_serializer());

  // Prune properties if needed.
  if (!have_distances)
  {
    point_cloud_delete_property(mesh, "d1");
    point_cloud_delete_property(mesh, "d2");
  }
  if (!have_beta)
    point_cloud_delete_property(mesh, "beta");
  if (!have_sigma)
    point_cloud_delete_property(mesh, "sigma");

  // Clean up.
  int_unordered_set_free(bad_connections);

  return mesh;
}

bool tough2_file_contains_mesh(tough2_file_t* file)
{
  return (string_ptr_unordered_map_contains(file->text, "eleme") && 
          string_ptr_unordered_map_contains(file->text, "conne"));
}

int_ptr_unordered_map_t* tough2_file_read_incon(tough2_file_t* file, 
                                                point_cloud_t* mesh,
                                                int num_primaries)
{
  ASSERT(mesh != NULL);
  ASSERT(num_primaries > 0);
  string_array_t* incon_lines = *string_ptr_unordered_map_get(file->text, "incon");

  // Do we have an even number of records?
  if ((incon_lines->size % 2) != 0)
  {
    file->valid = false;
    snprintf(file->error, 1023, "Each record in INCON must be two lines.");
    return NULL;
  }

  // Try to find the mesh's element name -> index map.
  string_int_unordered_map_t* elem_index_map = point_cloud_property(mesh, "elem_index");
  if (elem_index_map == NULL)
  {
    file->valid = false;
    snprintf(file->error, 1023, "Could not find 'elem_index' property in mesh.");
    return NULL;
  }

  // Count up the element entries (taking into account NSEQ and NADD).
  int num_ics = 0;
  for (int i = 0; i < incon_lines->size/2; ++i)
  {
    char* line = incon_lines->data[i];
    int offset = 0;

    // Parse the element name.
    char el[4], ne_str[3];
    strncpy(el, &line[offset], 3); offset += 3;
    strncpy(ne_str, &line[offset], 2); offset += 2;
    if (!string_is_number((const char*)ne_str))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing initial condition %d: invalid element name.", num_ics);
      break;
    }
    int ne = atoi(ne_str);
    if (ne < 0)
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing initial condition %d: invalid element code.", num_ics);
      break;
    }

    // Cross-reference it with the mesh.
    char elem_name[6];
    snprintf(elem_name, 6, "%s%s", el, ne_str);
    
    // Parse the NSEQ, NAD1, and NAD2 values, and reject any non-blank entries!
    char nseq_str[6], nadd_str[6];
    strncpy(nseq_str, &line[offset], 5); offset += 5;
    strncpy(nadd_str, &line[offset], 5); offset += 5;
    int nseq = 0, nadd = 1;
    if (string_is_number(nseq_str))
    {
      if (nadd <= 0)
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing initial condition %d: NSEQ must be a positive number.", num_ics);
        break;
      }
      if (string_is_number(nadd_str))
      {
        nadd = atoi(nadd_str);
        if (nadd <= 0)
        {
          file->valid = false;
          snprintf(file->error, 1023, "Error parsing initial condition %d: NADD must be a positive number.", num_ics);
          break;
        }
      }
      else if (!string_is_blank(nadd_str))
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing initial condition %d: NADD must be a number.", num_ics);
        break;
      }
    }
    else if (!string_is_blank(nseq_str))
    {
      file->valid = false;
      snprintf(file->error, 1023, "Error parsing initial condition %d: NSEQ must be a number.", num_ics);
      break;
    }

    for (int j = 0; j <= nseq; ++j, ++num_ics)
    {
      char elem_name[6];
      snprintf(elem_name, 6, "%s%d", el, ne + j*nadd);
      int* index_p = string_int_unordered_map_get(elem_index_map, elem_name);
      if (index_p == NULL)
      {
        file->valid = false;
        snprintf(file->error, 1023, "Initial condition %d refers to non-existent element %s.", num_ics, elem_name);
        break;
      }
    }
    if (!file->valid)
      break;
  }

  if (!file->valid)
    return NULL;

  int_ptr_unordered_map_t* incon = int_ptr_unordered_map_new();
  real_t* ic_pool = polymec_malloc(sizeof(real_t) * (num_primaries + 1) * num_ics);
  bool first_var = true;
  log_debug("tough2_file_read_mesh: Parsing %d initial conditions...", num_ics);
  int ic = 0;
  for (int i = 0; i < incon_lines->size/2; ++i)
  {
    // Parse the first line.
    char* line = incon_lines->data[2*i];
    int line_len = strlen(line);
    int offset = 0;

    // Parse the element name(s) and find indices.
    char el[4], ne_str[3];
    strncpy(el, &line[offset], 3); offset += 3;
    strncpy(ne_str, &line[offset], 2); offset += 2;
    int ne = atoi(ne_str);
    char elem_name[6];
    snprintf(elem_name, 6, "%s%s", el, ne_str);

    char nseq_str[6], nadd_str[6];
    strncpy(nseq_str, &line[offset], 5); offset += 5;
    strncpy(nadd_str, &line[offset], 5); offset += 5;
    int nseq = 0, nadd = 1;
    if (string_is_number(nseq_str))
      nseq = atoi(nseq_str);
    if (string_is_number(nadd_str))
      nadd = atoi(nadd_str);

    // Parse the porosity.
    real_t porosity = 0.0;
    if (offset < line_len)
    {
      char porx[16];
      strncpy(porx, &line[offset], 15); offset += 15;
      if (string_is_number(porx))
      {
        porosity = (real_t)atof(porx);
        if ((porosity < 0.0) || (porosity > 1.0))
        {
          file->valid = false;
          snprintf(file->error, 1023, "Error parsing initial condition %d: invalid porosity.", ic);
          break;
        }
      }
      else if (!string_is_blank(porx))
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing initial condition %d: invalid porosity.", ic);
        break;
      }
    }

    // Parse the second line.
    line = incon_lines->data[2*i+1];
    line_len = strlen(line);
    offset = 0;

    // Primary variables.
    real_t primaries[num_primaries];
    for (int p = 0; p < num_primaries; ++p)
    {
      char var[21];
      strncpy(var, &line[offset], 20); offset += 20;
      if (string_is_number(var))
        primaries[p] = atof(var);
      else
      {
        file->valid = false;
        snprintf(file->error, 1023, "Error parsing initial condition %d: invalid primary variable %d.", ic, p);
        break;
      }
    }
    if (!file->valid)
      break;

    // Now set things up.
    for (int j = 0; j <= nseq; ++j, ++ic)
    {
      char elem_name[6];
      snprintf(elem_name, 6, "%s%d", el, ne + j*nadd);
      int elem = *string_int_unordered_map_get(elem_index_map, elem_name);

      // Copy the data into place.
      real_t* ic_data = &ic_pool[(num_primaries+1) * ic];
      for (int p = 0; p < num_primaries; ++p)
        ic_data[p] = primaries[p];
      ic_data[num_primaries] = porosity;

      // Insert the data. The first entry into the map gets the privilege of managing memory.
      if (first_var)
      {
        ASSERT(ic_data == ic_pool);
        int_ptr_unordered_map_insert_with_v_dtor(incon, elem, ic_data, polymec_free);
        first_var = false;
      }
      else
        int_ptr_unordered_map_insert(incon, elem, ic_data);
    }
  }

  if (!file->valid)
  {
    int_ptr_unordered_map_free(incon);
    if (first_var)
      polymec_free(ic_pool);
    incon = NULL;
  }

  return incon;
}

bool tough2_file_contains_incon(tough2_file_t* file)
{
  return string_ptr_unordered_map_contains(file->text, "incon");
}

char* tough2_file_error(tough2_file_t* file)
{
  return &file->error[0];
}
