// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "polyglot/exodus_file.h"

#if POLYMEC_HAVE_MPI
#include "mpi.h"
#define PARALLEL_AWARE_EXODUS 1
#include "exodusII_par.h"
#endif 

#include "exodusII.h"

// This helper function converts the given element identifier string and number of nodes to 
// our own element enumerated type.
static fe_mesh_element_t get_element_type(const char* elem_type_id, int num_nodes_per_elem)
{
  if (string_ncasecmp(elem_type_id, "nfaced", 6) == 0)
  {
    ASSERT(num_nodes_per_elem == 0);
    return FE_POLYHEDRON;
  }
  else if (string_ncasecmp(elem_type_id, "tetra", 5) == 0)
  {
    if (num_nodes_per_elem == 4)
      return FE_TETRAHEDRON_4;
    else if (num_nodes_per_elem == 8)
      return FE_TETRAHEDRON_8;
    else if (num_nodes_per_elem == 10)
      return FE_TETRAHEDRON_10;
    else 
    {
      ASSERT(num_nodes_per_elem == 14);
      return FE_TETRAHEDRON_14;
    }
  }
  else if (string_ncasecmp(elem_type_id, "pyramid", 7) == 0)
  {
    if (num_nodes_per_elem == 5)
      return FE_PYRAMID_5;
    else 
    {
      ASSERT(num_nodes_per_elem == 13);
      return FE_PYRAMID_13;
    }
  }
  else if (string_ncasecmp(elem_type_id, "wedge", 5) == 0)
  {
    if (num_nodes_per_elem == 6)
      return FE_WEDGE_6;
    else if (num_nodes_per_elem == 15)
      return FE_WEDGE_15;
    else 
    {
      ASSERT(num_nodes_per_elem == 16);
      return FE_WEDGE_15;
    }
  }
  else if (string_ncasecmp(elem_type_id, "hex", 3) == 0)
  {
    if (num_nodes_per_elem == 8)
      return FE_HEXAHEDRON_8;
    else if (num_nodes_per_elem == 9)
      return FE_HEXAHEDRON_9;
    else if (num_nodes_per_elem == 20)
      return FE_HEXAHEDRON_20;
    else 
    {
      ASSERT(num_nodes_per_elem == 27);
      return FE_HEXAHEDRON_27;
    }
  }
  else
    return FE_INVALID;
}

struct exodus_file_t 
{
#if POLYMEC_HAVE_MPI
  MPI_Comm comm;        // Parallel communicator.
  MPI_Info mpi_info;    // Parallel info returned by NetCDF.
#endif
  int ex_id;            // Exodus file descriptor/identifier.
  float ex_version;     // Exodus database version number.
  int ex_real_size;     // Word size of data in the present Exodus file.
  int last_time_index;  // Index of most-recently added time written to file.

  int num_nodes, num_edges, num_faces, num_elem, 
      num_elem_blocks, num_face_blocks, num_edge_blocks;

  // Variable names.
  int num_node_vars, num_node_set_vars,
      num_edge_vars, num_edge_set_vars,
      num_face_vars, num_face_set_vars,
      num_elem_vars, num_elem_set_vars, num_side_set_vars;
  char **node_var_names, **node_set_var_names,
       **edge_var_names, **edge_set_var_names,
       **face_var_names, **face_set_var_names,
       **elem_var_names, **elem_set_var_names, **side_set_var_names;
};

bool exodus_file_query(const char* filename,
                       int* real_size,
                       float* version,
                       int* num_mpi_processes,
                       real_array_t* times)
{
  bool valid = true;
  int my_real_size = (int)sizeof(real_t);
  *real_size = 0;
#if POLYMEC_HAVE_MPI
  MPI_Info info;
  MPI_Info_create(&info);
  int id = ex_open_par(filename, EX_READ, &my_real_size,
                       real_size, version, 
                       MPI_COMM_WORLD, info);
#else
  int id = ex_open(filename, EX_READ, &my_real_size,
                   &real_size, &version);
#endif
  if (id < 0)
    valid = false;
  else
  {
    // Make sure that the file has 3D data.
    ex_init_params mesh_info;
    int status = ex_get_init_ext(id, &mesh_info);
    if ((status < 0) || (mesh_info.num_dim != 3))
      valid = false;
    else
    {
      // Query the number of processes for which this file has data.
      int num_proc_in_file;
      char file_type[2];
      ex_get_init_info(id, num_mpi_processes, &num_proc_in_file, file_type);
      ASSERT(*num_mpi_processes == num_proc_in_file);

      if (times != NULL)
      {
        // Ask for the times within the file.
        int num_times = ex_inquire_int(id, EX_INQ_TIME);
        real_array_resize(times, num_times);
        ex_get_all_times(id, times->data);
      }
    }

    ex_close(id);
  }

#if POLYMEC_HAVE_MPI
  MPI_Info_free(&info);
#endif

  return valid;
}

static void fetch_variable_names(int ex_id, ex_entity_type obj_type, int* num_vars, char*** var_names)
{
  ex_get_variable_param(ex_id, obj_type, num_vars);
  *var_names = polymec_malloc(sizeof(char*) * (*num_vars));
  for (int i = 0; i < (*num_vars); ++i)
    (*var_names)[i] = polymec_malloc(sizeof(char) * (MAX_NAME_LENGTH+1));
  ex_get_variable_names(ex_id, obj_type, *num_vars, *var_names);
}

static void fetch_all_variable_names(exodus_file_t* file)
{
  fetch_variable_names(file->ex_id, EX_NODAL, &file->num_node_vars, &file->node_var_names);
  fetch_variable_names(file->ex_id, EX_NODE_SET, &file->num_node_set_vars, &file->node_set_var_names);
  fetch_variable_names(file->ex_id, EX_EDGE_BLOCK, &file->num_edge_vars, &file->edge_var_names);
  fetch_variable_names(file->ex_id, EX_EDGE_SET, &file->num_edge_set_vars, &file->edge_set_var_names);
  fetch_variable_names(file->ex_id, EX_FACE_BLOCK, &file->num_face_vars, &file->face_var_names);
  fetch_variable_names(file->ex_id, EX_FACE_SET, &file->num_face_set_vars, &file->face_set_var_names);
  fetch_variable_names(file->ex_id, EX_ELEM_BLOCK, &file->num_elem_vars, &file->elem_var_names);
  fetch_variable_names(file->ex_id, EX_ELEM_SET, &file->num_elem_set_vars, &file->elem_set_var_names);
  fetch_variable_names(file->ex_id, EX_SIDE_SET, &file->num_side_set_vars, &file->side_set_var_names);
}

static void free_variable_names(int ex_id, int num_vars, char** var_names)
{
  for (int i = 0; i < num_vars; ++i)
    string_free(var_names[i]);
  polymec_free(var_names);
}

static void free_all_variable_names(exodus_file_t* file)
{
  free_variable_names(file->ex_id, file->num_node_vars, file->node_var_names);
  free_variable_names(file->ex_id, file->num_node_set_vars, file->node_set_var_names);
  free_variable_names(file->ex_id, file->num_edge_vars, file->edge_var_names);
  free_variable_names(file->ex_id, file->num_edge_set_vars, file->edge_set_var_names);
  free_variable_names(file->ex_id, file->num_face_vars, file->face_var_names);
  free_variable_names(file->ex_id, file->num_face_set_vars, file->face_set_var_names);
  free_variable_names(file->ex_id, file->num_elem_vars, file->elem_var_names);
  free_variable_names(file->ex_id, file->num_elem_set_vars, file->elem_set_var_names);
  free_variable_names(file->ex_id, file->num_side_set_vars, file->side_set_var_names);
}

static exodus_file_t* open_exodus_file(MPI_Comm comm,
                                       const char* filename,
                                       int mode)
{
  exodus_file_t* file = polymec_malloc(sizeof(exodus_file_t));
  file->last_time_index = 0;
  file->comm = comm;
  int real_size = (int)sizeof(real_t);
  file->ex_real_size = 0;
#if POLYMEC_HAVE_MPI
  MPI_Info_create(&file->mpi_info);
  file->ex_id = ex_open_par(filename, mode, &real_size,
                            &file->ex_real_size, &file->ex_version, 
                            file->comm, file->mpi_info);
#else
  file->ex_id = ex_open(filename, mode, &real_size,
                        &file->ex_real_size, &file->ex_version);
#endif
  if (file->ex_id >= 0)
  {
    // Read all the available variable names.
    fetch_all_variable_names(file);

    // Get information from the file.
    ex_init_params mesh_info;
    int status = ex_get_init_ext(file->ex_id, &mesh_info);
    if ((status >= 0) && (mesh_info.num_dim == 3))
    {
      file->num_nodes = mesh_info.num_nodes;
      file->num_elem = mesh_info.num_elem;
      file->num_faces = mesh_info.num_face;
      file->num_edges = mesh_info.num_edge;
      file->num_elem_blocks = mesh_info.num_elem_blk;
      file->num_face_blocks = mesh_info.num_face_blk;
      file->num_edge_blocks = mesh_info.num_edge_blk;
    }
  }
  else
  {
    polymec_free(file);
    file = NULL;
  }

  return file;
}

exodus_file_t* exodus_file_new(MPI_Comm comm,
                               const char* filename)
{
  return open_exodus_file(comm, filename, EX_WRITE);
}

exodus_file_t* exodus_file_open(MPI_Comm comm,
                                const char* filename)
{
  return open_exodus_file(comm, filename, EX_READ);
}

void exodus_file_close(exodus_file_t* file)
{
  // Clean up.
  free_all_variable_names(file);
#if POLYMEC_HAVE_MPI
  MPI_Info_free(&file->mpi_info);
#endif

  ex_close(file->ex_id);
}

void exodus_file_write_fe_mesh(exodus_file_t* file,
                               fe_mesh_t* mesh)
{
}

fe_mesh_t* exodus_file_read_fe_mesh(exodus_file_t* file)
{
  // Create the "host" FE mesh.
  fe_mesh_t* mesh = fe_mesh_new(file->comm, file->num_nodes);

  // Go over the element blocks and feel out the data.
  for (int elem_block = 1; elem_block <= file->num_elem_blocks; ++elem_block)
  {
    char elem_type_name[MAX_NAME_LENGTH+1];
    int num_elem, num_nodes_per_elem, num_faces_per_elem;
    ex_get_block(file->ex_id, EX_ELEM_BLOCK, elem_block, 
                 elem_type_name, &num_elem,
                 &num_nodes_per_elem, NULL,
                 &num_faces_per_elem, NULL);

    // Get the type of element for this block.
    fe_mesh_element_t elem_type = get_element_type(elem_type_name, num_nodes_per_elem);
    fe_block_t* block = NULL;
    char block_name[MAX_NAME_LENGTH+1];
    if (elem_type == FE_POLYHEDRON)
    {
      // Dig up the face block corresponding to this element block.
      // FIXME: Do we know that it shares the same ID?
      int face_block = elem_block;
      char face_type[MAX_NAME_LENGTH+1];
      int num_faces, num_nodes;
      ex_get_block(file->ex_id, EX_FACE_BLOCK, face_block, face_type, &num_faces,
                   &num_nodes, NULL, NULL, NULL);
      if (string_ncasecmp(face_type, "nsided", 6) != 0)
      {
        fe_mesh_free(mesh);
        ex_close(file->ex_id);
        polymec_error("Invalid face type for polyhedral element block.");
      }

      // Find the number of faces for each element in the block.
      int* num_elem_faces = polymec_malloc(sizeof(int) * num_elem);
      ex_get_entity_count_per_polyhedra(file->ex_id, EX_ELEM_BLOCK, elem_block, 
                                          num_elem_faces);

      // Get the element->face connectivity.
      int elem_face_size = 0;
      for (int i = 0; i < num_elem; ++i)
        elem_face_size += num_elem_faces[i];
      int* elem_faces = polymec_malloc(sizeof(int) * elem_face_size);
      ex_get_conn(file->ex_id, EX_ELEM_BLOCK, elem_block, NULL, NULL, elem_faces);

      // Find the number of nodes for each face in the block.
      int* num_face_nodes = polymec_malloc(sizeof(int) * num_faces);
      ex_get_entity_count_per_polyhedra(file->ex_id, EX_FACE_BLOCK, face_block, 
                                        num_face_nodes);

      // Get the face->node connectivity.
      int face_node_size = 0;
      for (int i = 0; i < num_faces; ++i)
        face_node_size += num_face_nodes[i];
      int* face_nodes = polymec_malloc(sizeof(int) * face_node_size);
      ex_get_conn(file->ex_id, EX_FACE_BLOCK, face_block, face_nodes, NULL, NULL);

      // Create the element block.
      block = fe_polyhedral_block_new(num_elem, num_elem_faces, elem_faces);
    }
    else if (elem_type != FE_INVALID)
    {
      // Get the element's nodal mapping.
      int* node_conn = polymec_malloc(sizeof(int) * num_elem * num_nodes_per_elem);
      ex_get_conn(file->ex_id, EX_ELEM_BLOCK, elem_block, node_conn, NULL, NULL);
      
      // Build the element block.
      block = fe_block_new(num_elem, elem_type, node_conn);
    }
    else
    {
      fe_mesh_free(mesh);
      ex_close(file->ex_id);
      polymec_error("Block %d contains an invalid (3D) element type.", elem_block);
    }

    // Fish out the element block name if it has one, or make a default.
    ex_get_name(file->ex_id, EX_ELEM_BLOCK, elem_block, block_name);
    if (strlen(block_name) == 0)
      sprintf(block_name, "block_%d", elem_block);

    // Add the element block to the mesh.
    fe_mesh_add_block(mesh, block_name, block);
  }

  // Fetch node positions and compute geometry.
  real_t x[file->num_nodes], y[file->num_nodes], z[file->num_nodes];
  ex_get_coord(file->ex_id, x, y, z);
  point_t* X = fe_mesh_node_coordinates(mesh);
  for (int n = 0; n < file->num_nodes; ++n)
  {
    X[n].x = x[n];
    X[n].y = y[n];
    X[n].z = z[n];
  }

  return mesh;
}

int exodus_file_write_time(exodus_file_t* file, real_t time)
{
  int next_index = file->last_time_index + 1;
  int status = ex_put_time(file->ex_id, next_index, &time);
  if (status >= 0)
    file->last_time_index = next_index;
  else 
    next_index = status;
  return next_index;
}

bool exodus_file_next_time(exodus_file_t* file,
                           int* pos,
                           int* time_index,
                           real_t* time)
{
  if (*pos >= file->last_time_index)
    return false;

  int next_index = *pos + 1;
  int status = ex_get_time(file->ex_id, next_index, time);
  if (status >= 0)
  {
    *time_index = *pos = next_index;
    return true;
  }
  else
    return false;
}

void exodus_file_write_element_field(exodus_file_t* file,
                                     int time_index,
                                     const char* field_name,
                                     real_t* field_data)
{
}

real_t* exodus_file_read_element_field(exodus_file_t* file,
                                       int time_index,
                                       const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->num_elem_vars)
  {
    if (strcmp(field_name, file->elem_var_names[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->num_elem_vars)
  {
    int offset = 0;
    real_t* field = polymec_malloc(sizeof(real_t) * file->num_elem);
    memset(field, 0, sizeof(real_t) * file->num_elem);
    for (int i = 1; i <= file->num_elem_blocks; ++i)
    {
      int N;
      ex_get_block(file->ex_id, EX_ELEM_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
      ex_get_var(file->ex_id, time_index, EX_ELEM_BLOCK, index+1, i, N, &field[offset]);
      offset += N;
    }
    return field;
  }
  else
    return NULL;
}

bool exodus_file_contains_element_field(exodus_file_t* file, 
                                        int time_index,
                                        const char* field_name)
{
  int index = 0;
  while (index < file->num_elem_vars)
  {
    if (strcmp(field_name, file->elem_var_names[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_face_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
}

real_t* exodus_file_read_face_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->num_face_vars)
  {
    if (strcmp(field_name, file->face_var_names[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->num_face_vars)
  {
    int offset = 0;
    real_t* field = polymec_malloc(sizeof(real_t) * file->num_faces);
    memset(field, 0, sizeof(real_t) * file->num_faces);
    for (int i = 1; i <= file->num_face_blocks; ++i)
    {
      int N;
      ex_get_block(file->ex_id, EX_FACE_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
      ex_get_var(file->ex_id, time_index, EX_FACE_BLOCK, index+1, i, N, &field[offset]);
      offset += N;
    }
    return field;
  }
  else
    return NULL;
}

bool exodus_file_contains_face_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
  int index = 0;
  while (index < file->num_face_vars)
  {
    if (strcmp(field_name, file->face_var_names[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_edge_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
}

real_t* exodus_file_read_edge_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->num_edge_vars)
  {
    if (strcmp(field_name, file->edge_var_names[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->num_edge_vars)
  {
    int offset = 0;
    real_t* field = polymec_malloc(sizeof(real_t) * file->num_edges);
    memset(field, 0, sizeof(real_t) * file->num_edges);
    for (int i = 1; i <= file->num_edge_blocks; ++i)
    {
      int N;
      ex_get_block(file->ex_id, EX_EDGE_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
      ex_get_var(file->ex_id, time_index, EX_EDGE_BLOCK, index+1, i, N, &field[offset]);
      offset += N;
    }
    return field;
  }
  else
    return NULL;
}

bool exodus_file_contains_edge_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
  int index = 0;
  while (index < file->num_edge_vars)
  {
    if (strcmp(field_name, file->edge_var_names[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_node_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
}

real_t* exodus_file_read_node_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->num_node_vars)
  {
    if (strcmp(field_name, file->node_var_names[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->num_node_vars)
  {
    real_t* field = polymec_malloc(sizeof(real_t) * file->num_nodes);
    memset(field, 0, sizeof(real_t) * file->num_nodes);
    ex_get_var(file->ex_id, time_index, EX_NODAL, index+1, 1, file->num_nodes, field);
    return field;
  }
  else
    return NULL;
}

bool exodus_file_contains_node_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
  int index = 0;
  while (index < file->num_node_vars)
  {
    if (strcmp(field_name, file->node_var_names[index]) == 0)
      return true;
  }
  return false;
}

