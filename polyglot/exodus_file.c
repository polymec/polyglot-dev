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

static void get_elem_info(fe_mesh_element_t elem_type,
                          char* elem_type_name,
                          int* num_elem_nodes)
{
  // Figure out the element type name.
  switch(elem_type)
  {
    case FE_TETRAHEDRON_4: 
    case FE_TETRAHEDRON_8:
    case FE_TETRAHEDRON_10: 
    case FE_TETRAHEDRON_14: strcpy(elem_type_name, "TETRA"); break;
    case FE_PYRAMID_5: 
    case FE_PYRAMID_13: strcpy(elem_type_name, "PYRAMID"); break;
    case FE_WEDGE_6: 
    case FE_WEDGE_15: strcpy(elem_type_name, "WEDGE"); break;
    case FE_HEXAHEDRON_8: 
    case FE_HEXAHEDRON_9: 
    case FE_HEXAHEDRON_20:
    case FE_HEXAHEDRON_27: strcpy(elem_type_name, "HEX"); break;
    default: strcpy(elem_type_name, "UNKNOWN");
  }

  // Figure out the number of nodes.
  switch(elem_type)
  {
    case FE_TETRAHEDRON_4: *num_elem_nodes = 4; break;
    case FE_PYRAMID_5: *num_elem_nodes = 5; break;
    case FE_WEDGE_6: *num_elem_nodes = 6; break;
    case FE_TETRAHEDRON_8:
    case FE_HEXAHEDRON_8: *num_elem_nodes = 8; break;
    case FE_HEXAHEDRON_9: *num_elem_nodes = 9; break;
    case FE_TETRAHEDRON_10: *num_elem_nodes = 10; break;
    case FE_PYRAMID_13: *num_elem_nodes = 13; break;
    case FE_TETRAHEDRON_14: *num_elem_nodes = 14; break;
    case FE_WEDGE_15: *num_elem_nodes = 15; break;
    case FE_HEXAHEDRON_20: *num_elem_nodes = 20; break;
    case FE_HEXAHEDRON_27: *num_elem_nodes = 27; break;
    default: *num_elem_nodes = -1;
  }
}

struct exodus_file_t 
{
  char title[MAX_NAME_LENGTH+1];

#if POLYMEC_HAVE_MPI
  MPI_Comm comm;        // Parallel communicator.
  MPI_Info mpi_info;    // Parallel info returned by NetCDF.
#endif

  int ex_id;            // Exodus file descriptor/identifier.
  float ex_version;     // Exodus database version number.
  int ex_real_size;     // Word size of data in the present Exodus file.
  int last_time_index;  // Index of most-recently added time written to file.

  // Set to true if we're writing to an Exodus file, false if not.
  bool writing;

  int num_nodes, num_edges, num_faces, num_elem, 
      num_elem_blocks, num_face_blocks, num_edge_blocks,
      num_elem_sets, num_face_sets, num_edge_sets, num_node_sets, num_side_sets;

  // Variable names.
  string_array_t *node_var_names, *node_set_var_names,
                 *edge_var_names, *edge_set_var_names,
                 *face_var_names, *face_set_var_names,
                 *elem_var_names, *elem_set_var_names, *side_set_var_names;
};

bool exodus_file_query(const char* filename,
                       size_t* real_size,
                       float* version,
                       int* num_mpi_processes,
                       real_array_t* times)
{
  bool valid = true;
  bool is_parallel = false;
  int my_real_size = (int)sizeof(real_t);
  int io_real_size = 0;
#if POLYMEC_HAVE_MPI
  MPI_Info info;
  MPI_Info_create(&info);
  int id = ex_open_par(filename, EX_READ, &my_real_size,
                       &io_real_size, version, 
                       MPI_COMM_WORLD, info);

  // Did that work? If not, try the serial opener.
  if (id < 0)
  {
    MPI_Info_free(&info);
    id = ex_open(filename, EX_READ, &my_real_size,
                 &io_real_size, version);
  }
  else
    is_parallel = true;
#else
  int id = ex_open(filename, EX_READ, &my_real_size,
                   &io_real_size, version);
#endif

  if (id < 0)
    valid = false;
  else
  {
    *real_size = (size_t)io_real_size;

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
      if (is_parallel)
      {
        ASSERT(*num_mpi_processes == num_proc_in_file);
      }

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
  if (is_parallel)
    MPI_Info_free(&info);
#endif

  return valid;
}

static void fetch_variable_names(int ex_id, ex_entity_type obj_type, string_array_t* var_names)
{
  int num_vars;
  ex_get_variable_param(ex_id, obj_type, &num_vars);
  for (int i = 0; i < num_vars; ++i)
    string_array_append_with_dtor(var_names, (char*)polymec_malloc(sizeof(char) * (MAX_NAME_LENGTH+1)), string_free);
  ex_get_variable_names(ex_id, obj_type, num_vars, var_names->data);
}

static void fetch_all_variable_names(exodus_file_t* file)
{
  fetch_variable_names(file->ex_id, EX_NODAL, file->node_var_names);
  fetch_variable_names(file->ex_id, EX_NODE_SET, file->node_set_var_names);
  fetch_variable_names(file->ex_id, EX_EDGE_BLOCK, file->edge_var_names);
  fetch_variable_names(file->ex_id, EX_EDGE_SET, file->edge_set_var_names);
  fetch_variable_names(file->ex_id, EX_FACE_BLOCK, file->face_var_names);
  fetch_variable_names(file->ex_id, EX_FACE_SET, file->face_set_var_names);
  fetch_variable_names(file->ex_id, EX_ELEM_BLOCK, file->elem_var_names);
  fetch_variable_names(file->ex_id, EX_ELEM_SET, file->elem_set_var_names);
  fetch_variable_names(file->ex_id, EX_SIDE_SET, file->side_set_var_names);
}

static void free_all_variable_names(exodus_file_t* file)
{
  string_array_free(file->node_var_names);
  string_array_free(file->node_set_var_names);
  string_array_free(file->edge_var_names);
  string_array_free(file->edge_set_var_names);
  string_array_free(file->face_var_names);
  string_array_free(file->face_set_var_names);
  string_array_free(file->elem_var_names);
  string_array_free(file->elem_set_var_names);
  string_array_free(file->side_set_var_names);
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

  // Did that work? If not, try the serial opener.
  if (file->ex_id < 0)
  {
    file->ex_id = ex_open(filename, mode, &real_size,
                          &file->ex_real_size, &file->ex_version);
  }
#else
  file->ex_id = ex_open(filename, mode, &real_size,
                        &file->ex_real_size, &file->ex_version);
#endif
  if (file->ex_id >= 0)
  {
    file->writing = (mode == EX_WRITE);
    file->node_var_names = string_array_new();
    file->node_set_var_names = string_array_new();
    file->edge_var_names = string_array_new();
    file->edge_set_var_names = string_array_new();
    file->face_var_names = string_array_new();
    file->face_set_var_names = string_array_new();
    file->elem_var_names = string_array_new();
    file->elem_set_var_names = string_array_new();
    file->side_set_var_names = string_array_new();

    if (!file->writing)
    {
      // Read all the available variable names.
      fetch_all_variable_names(file);

      // Get information from the file.
      ex_init_params mesh_info;
      int status = ex_get_init_ext(file->ex_id, &mesh_info);
      if ((status >= 0) && (mesh_info.num_dim == 3))
      {
        strncpy(file->title, mesh_info.title, MAX_NAME_LENGTH);
        file->num_nodes = mesh_info.num_nodes;
        file->num_elem = mesh_info.num_elem;
        file->num_faces = mesh_info.num_face;
        file->num_edges = mesh_info.num_edge;
        file->num_elem_blocks = mesh_info.num_elem_blk;
        file->num_face_blocks = mesh_info.num_face_blk;
        file->num_edge_blocks = mesh_info.num_edge_blk;
        file->num_elem_sets = mesh_info.num_elem_sets;
        file->num_face_sets = mesh_info.num_face_sets;
        file->num_edge_sets = mesh_info.num_edge_sets;
        file->num_side_sets = mesh_info.num_side_sets;
      }
    }
    else
    {
      // By default, the title of the database is its filename.
      strncpy(file->title, filename, MAX_NAME_LENGTH);
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
  if (file->writing)
  {
    // Write a QA record.
    char* qa_record[1][4];
    qa_record[0][0] = string_dup(polymec_executable_name());
    qa_record[0][1] = string_dup(polymec_executable_name());
    time_t invocation_time = polymec_invocation_time();
    struct tm* time_data = localtime(&invocation_time);
    char date[20], instant[20];
    snprintf(date, 19, "%2d/%2d/%2d", time_data->tm_mon, time_data->tm_mday, 
             time_data->tm_year % 100);
    qa_record[0][1] = string_dup(date);
    snprintf(instant, 19, "%2d:%2d:%2d", time_data->tm_hour, time_data->tm_min, 
             time_data->tm_sec % 60);
    qa_record[0][1] = string_dup(instant);
    ex_put_qa(file->ex_id, 1, qa_record);
    for (int i = 0; i < 4; ++i)
      string_free(qa_record[0][i]);
  }

  // Clean up.
  free_all_variable_names(file);
#if POLYMEC_HAVE_MPI
  MPI_Info_free(&file->mpi_info);
#endif

  ex_close(file->ex_id);
}

char* exodus_file_title(exodus_file_t* file)
{
  return file->title;
}

void exodus_file_set_title(exodus_file_t* file, const char* title)
{
  strncpy(file->title, title, MAX_NAME_LENGTH);
}

static void write_set(exodus_file_t* file, 
                      ex_entity_type set_type,
                      int set_id,
                      char* set_name,
                      int* set,
                      int set_size)
{
  ex_put_name(file->ex_id, EX_ELEM_SET, (ex_entity_id)set_id, set_name);
  int num_dist_factors = 0;
  ex_put_set_param(file->ex_id, EX_ELEM_SET, (ex_entity_id)set_id, set_size, num_dist_factors);
  if (set_type != EX_SIDE_SET)
    ex_put_set(file->ex_id, EX_ELEM_SET, (ex_entity_id)set_id, set, NULL);
  else
  {
    int elems[set_size], faces[set_size];
    for (int i = 0; i < set_size; ++i)
    {
      elems[i] = set[2*i];
      faces[i] = set[2*i+1];
    }
    ex_put_set(file->ex_id, EX_ELEM_SET, (ex_entity_id)set_id, elems, faces);
  }
}

void exodus_file_write_fe_mesh(exodus_file_t* file,
                               fe_mesh_t* mesh)
{
  ASSERT(file->writing);

  // Count up the number of polyhedral blocks.
  int num_blocks = fe_mesh_num_blocks(mesh);
  int num_poly_blocks = 0;
  int pos = 0;
  char* block_name;
  fe_block_t* block;
  while (fe_mesh_next_block(mesh, &pos, &block_name, &block))
  {
    fe_mesh_element_t elem_type = fe_block_element_type(block);
    if (elem_type == FE_POLYHEDRON)
      ++num_poly_blocks;
  }

  // Write out information about elements, faces, edges, nodes.
  file->num_nodes = fe_mesh_num_nodes(mesh);
  ex_init_params params;
  strcpy(params.title, file->title);
  params.num_dim = 3;
  params.num_nodes = file->num_nodes;
  int num_edges = fe_mesh_num_edges(mesh);
  params.num_edge = num_edges;
  params.num_edge_blk = 0;
  int num_faces = fe_mesh_num_faces(mesh);
  params.num_face = num_faces;
  params.num_face_blk = num_poly_blocks;
  int num_elem = fe_mesh_num_elements(mesh);
  params.num_elem = num_elem;
  params.num_elem_blk = num_blocks;
  params.num_elem_sets = file->num_elem_sets = fe_mesh_num_element_sets(mesh);
  params.num_face_sets = file->num_face_sets = fe_mesh_num_face_sets(mesh);
  params.num_edge_sets = file->num_edge_sets = fe_mesh_num_edge_sets(mesh);
  params.num_node_sets = file->num_node_sets = fe_mesh_num_node_sets(mesh);
  params.num_side_sets = file->num_side_sets = fe_mesh_num_side_sets(mesh);
  params.num_elem_maps = 0;
  params.num_face_maps = 0;
  params.num_edge_maps = 0;
  params.num_node_maps = 0;
  ex_put_init_ext(file->ex_id, &params);

  // If we have any polyhedral element blocks, we write out a single face 
  // block that incorporates all of the polyhedral elements.
  if (num_poly_blocks > 0)
  {
    int num_poly_faces = 0, tot_poly_face_nodes = 0;
    // FIXME

    // Write an "nsided" face block.
    ex_put_block(file->ex_id, EX_FACE_BLOCK, 1, "nsided",
                 num_poly_faces, tot_poly_face_nodes, 0, 0, 0);
    ex_put_name(file->ex_id, EX_FACE_BLOCK, 1, "face_block");

    // Write face->node connectivity information.
  }

  // Go over the element blocks and write out the data.
  pos = 0;
  while (fe_mesh_next_block(mesh, &pos, &block_name, &block))
  {
    int elem_block = pos;
    int num_elem = fe_block_num_elements(block);
    fe_mesh_element_t elem_type = fe_block_element_type(block);
    if (elem_type == FE_POLYHEDRON)
    {
      // Count up the faces in the block and write the block information.
      int tot_num_elem_faces = 0;
      int faces_per_elem[num_elem];
      for (int i = 0; i < num_elem; ++i)
      {
        faces_per_elem[i] = fe_block_num_element_faces(block, i);
        tot_num_elem_faces += faces_per_elem[i];
      }
      ex_put_block(file->ex_id, EX_ELEM_BLOCK, elem_block, "nfaced", 
                   num_elem, 0, 0, tot_num_elem_faces, 0);

      // Write elem->face connectivity information.
      int elem_faces[tot_num_elem_faces], offset = 0;
      for (int i = 0; i < num_elem; ++i)
      {
        fe_block_get_element_faces(block, i, &elem_faces[offset]);
        offset += faces_per_elem[i];
      }
      ex_put_conn(file->ex_id, EX_ELEM_BLOCK, elem_block, NULL, NULL, elem_faces);
      ex_put_entity_count_per_polyhedra(file->ex_id, EX_ELEM_BLOCK, elem_block, faces_per_elem); 
    }
    else if (elem_type != FE_INVALID)
    {
      // Get element information.
      char elem_type_name[MAX_NAME_LENGTH+1];
      int num_nodes_per_elem;
      get_elem_info(elem_type, elem_type_name, &num_nodes_per_elem);

      // Write the block.
      ex_put_block(file->ex_id, EX_ELEM_BLOCK, elem_block, elem_type_name, 
                   num_elem, num_nodes_per_elem, 0, 0, 0);

      // Write the elem->node connectivity.
      int elem_nodes[num_elem * num_nodes_per_elem], offset = 0;
      for (int i = 0; i < num_elem; ++i)
      {
        fe_block_get_element_nodes(block, i, &elem_nodes[offset]);
        offset += num_nodes_per_elem;
      }
      ex_put_conn(file->ex_id, EX_ELEM_BLOCK, elem_block, elem_nodes, NULL, NULL);
    }

    // Set the element block name.
    ex_put_name(file->ex_id, EX_ELEM_BLOCK, elem_block, block_name);
  }

  // Set node positions.
  real_t x[file->num_nodes], y[file->num_nodes], z[file->num_nodes];
  point_t* X = fe_mesh_node_coordinates(mesh);
  for (int n = 0; n < file->num_nodes; ++n)
  {
    x[n] = X[n].x;
    y[n] = X[n].y;
    z[n] = X[n].z;
  }
  ex_put_coord(file->ex_id, x, y, z);

  // Write sets of entities.
  int *set, set_size;
  char* set_name;
  pos = 0;
  while (fe_mesh_next_element_set(mesh, &pos, &set_name, &set, &set_size))
    write_set(file, EX_ELEM_SET, pos, set_name, set, set_size);
  pos = 0;
  while (fe_mesh_next_face_set(mesh, &pos, &set_name, &set, &set_size))
    write_set(file, EX_FACE_SET, pos, set_name, set, set_size);
  pos = 0;
  while (fe_mesh_next_edge_set(mesh, &pos, &set_name, &set, &set_size))
    write_set(file, EX_EDGE_SET, pos, set_name, set, set_size);
  pos = 0;
  while (fe_mesh_next_node_set(mesh, &pos, &set_name, &set, &set_size))
    write_set(file, EX_NODE_SET, pos, set_name, set, set_size);
  pos = 0;
  while (fe_mesh_next_side_set(mesh, &pos, &set_name, &set, &set_size))
    write_set(file, EX_SIDE_SET, pos, set_name, set, set_size);
}

static void fetch_set(exodus_file_t* file, 
                      ex_entity_type set_type,
                      ex_entity_id set_id,
                      fe_mesh_t* mesh,
                      int* (*create_set)(fe_mesh_t* mesh, const char* name, int))
{
  char set_name[MAX_NAME_LENGTH+1];
  ex_get_name(file->ex_id, EX_ELEM_SET, set_id, set_name);
  int set_size, num_dist_factors;
  ex_get_set_param(file->ex_id, EX_ELEM_SET, set_id, &set_size, &num_dist_factors);
  int* set = create_set(mesh, set_name, set_size);
  ex_get_set(file->ex_id, EX_ELEM_SET, set_id, set, NULL);
}

fe_mesh_t* exodus_file_read_fe_mesh(exodus_file_t* file)
{
  // Create the "host" FE mesh.
  fe_mesh_t* mesh = fe_mesh_new(file->comm, file->num_nodes);

  // Count up the number of polyhedral blocks.
  int num_poly_blocks = 0;
  for (int elem_block = 1; elem_block <= file->num_elem_blocks; ++elem_block)
  {
    char elem_type_name[MAX_NAME_LENGTH+1];
    int num_elem, num_nodes_per_elem, num_faces_per_elem;
    ex_get_block(file->ex_id, EX_ELEM_BLOCK, elem_block, 
                 elem_type_name, &num_elem,
                 &num_nodes_per_elem, NULL,
                 &num_faces_per_elem, NULL);
    fe_mesh_element_t elem_type = get_element_type(elem_type_name, num_nodes_per_elem);
    if (elem_type == FE_POLYHEDRON)
      ++num_poly_blocks;
  }

  // If we have any polyhedral element blocks, we read a single face 
  // block that incorporates all of the polyhedral elements.
  if (num_poly_blocks > 0)
  {
    // Dig up the face block corresponding to this element block.
    char face_type[MAX_NAME_LENGTH+1];
    int num_faces, num_nodes;
    ex_get_block(file->ex_id, EX_FACE_BLOCK, 1, face_type, &num_faces,
                 &num_nodes, NULL, NULL, NULL);
    if (string_ncasecmp(face_type, "nsided", 6) != 0)
    {
      fe_mesh_free(mesh);
      ex_close(file->ex_id);
      polymec_error("Invalid face type for polyhedral element block.");
    }

    // Find the number of nodes for each face in the block.
    int* num_face_nodes = polymec_malloc(sizeof(int) * num_faces);
    ex_get_entity_count_per_polyhedra(file->ex_id, EX_FACE_BLOCK, 1, 
                                      num_face_nodes);

    // Read face->node connectivity information.
    int face_node_size = 0;
    for (int i = 0; i < num_faces; ++i)
      face_node_size += num_face_nodes[i];
    int* face_nodes = polymec_malloc(sizeof(int) * face_node_size);
    ex_get_conn(file->ex_id, EX_FACE_BLOCK, 1, face_nodes, NULL, NULL);
    fe_mesh_set_face_nodes(mesh, num_faces, num_face_nodes, face_nodes);

    // Clean up.
    polymec_free(num_face_nodes);
  }

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

  // Fetch sets of entities.
  for (int i = 1; i <= file->num_elem_sets; ++i)
    fetch_set(file, EX_ELEM_SET, i, mesh, fe_mesh_create_element_set);
  for (int i = 1; i <= file->num_face_sets; ++i)
    fetch_set(file, EX_FACE_SET, i, mesh, fe_mesh_create_face_set);
  for (int i = 1; i <= file->num_edge_sets; ++i)
    fetch_set(file, EX_EDGE_SET, i, mesh, fe_mesh_create_edge_set);
  for (int i = 1; i <= file->num_node_sets; ++i)
    fetch_set(file, EX_NODE_SET, i, mesh, fe_mesh_create_edge_set);
  for (int i = 1; i <= file->num_side_sets; ++i)
    fetch_set(file, EX_SIDE_SET, i, mesh, fe_mesh_create_side_set);

  return mesh;
}

int exodus_file_write_time(exodus_file_t* file, real_t time)
{
  ASSERT(file->writing);
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
  ASSERT(file->writing);

  // Find the variable index if it already exists.
  int index = 0;
  while (index < file->elem_var_names->size)
  {
    if (strcmp(field_name, file->elem_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Append the variable to our list if we didn't find it.
  if (index >= file->elem_var_names->size)
    string_array_append_with_dtor(file->elem_var_names, string_dup(field_name), string_free);

  // Insert the data.
  int offset = 0;
  for (int i = 1; i <= file->num_elem_blocks; ++i)
  {
    int N;
    ex_get_block(file->ex_id, EX_ELEM_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
    ex_put_var(file->ex_id, time_index, EX_ELEM_BLOCK, index+1, i, N, &field_data[offset]);
    offset += N;
  }
}

real_t* exodus_file_read_element_field(exodus_file_t* file,
                                       int time_index,
                                       const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->elem_var_names->size)
  {
    if (strcmp(field_name, file->elem_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->elem_var_names->size)
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
  while (index < file->elem_var_names->size)
  {
    if (strcmp(field_name, file->elem_var_names->data[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_face_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
  ASSERT(file->writing);

  // Find the variable index if it already exists.
  int index = 0;
  while (index < file->face_var_names->size)
  {
    if (strcmp(field_name, file->face_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Append the variable to our list if we didn't find it.
  if (index >= file->face_var_names->size)
    string_array_append_with_dtor(file->face_var_names, string_dup(field_name), string_free);

  // Insert the data.
  int offset = 0;
  for (int i = 1; i <= file->num_face_blocks; ++i)
  {
    int N;
    ex_get_block(file->ex_id, EX_FACE_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
    ex_put_var(file->ex_id, time_index, EX_FACE_BLOCK, index+1, i, N, &field_data[offset]);
    offset += N;
  }
}

real_t* exodus_file_read_face_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->face_var_names->size)
  {
    if (strcmp(field_name, file->face_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->face_var_names->size)
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
  while (index < file->face_var_names->size)
  {
    if (strcmp(field_name, file->face_var_names->data[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_edge_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
  ASSERT(file->writing);

  // Find the variable index if it already exists.
  int index = 0;
  while (index < file->edge_var_names->size)
  {
    if (strcmp(field_name, file->edge_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Append the variable to our list if we didn't find it.
  if (index >= file->edge_var_names->size)
    string_array_append_with_dtor(file->edge_var_names, string_dup(field_name), string_free);

  // Insert the data.
  int offset = 0;
  for (int i = 1; i <= file->num_edge_blocks; ++i)
  {
    int N;
    ex_get_block(file->ex_id, EX_EDGE_BLOCK, i, NULL, &N, NULL, NULL, NULL, NULL);
    ex_put_var(file->ex_id, time_index, EX_EDGE_BLOCK, index+1, i, N, &field_data[offset]);
    offset += N;
  }
}

real_t* exodus_file_read_edge_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->edge_var_names->size)
  {
    if (strcmp(field_name, file->edge_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->edge_var_names->size)
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
  while (index < file->edge_var_names->size)
  {
    if (strcmp(field_name, file->edge_var_names->data[index]) == 0)
      return true;
  }
  return false;
}

void exodus_file_write_node_field(exodus_file_t* file,
                                  int time_index,
                                  const char* field_name,
                                  real_t* field_data)
{
  ASSERT(file->writing);

  // Find the variable index if it already exists.
  int index = 0;
  while (index < file->node_var_names->size)
  {
    if (strcmp(field_name, file->node_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Append the variable to our list if we didn't find it.
  if (index >= file->node_var_names->size)
    string_array_append_with_dtor(file->node_var_names, string_dup(field_name), string_free);

  // Insert the data.
  ex_put_var(file->ex_id, time_index, EX_NODE_BLOCK, index+1, 1, file->num_nodes, field_data);
}

real_t* exodus_file_read_node_field(exodus_file_t* file,
                                    int time_index,
                                    const char* field_name)
{
  // Find the variable index.
  int index = 0;
  while (index < file->node_var_names->size)
  {
    if (strcmp(field_name, file->node_var_names->data[index]) == 0)
      break;
    ++index;
  }

  // Fetch the field data.
  if (index < file->node_var_names->size)
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
  while (index < file->node_var_names->size)
  {
    if (strcmp(field_name, file->node_var_names->data[index]) == 0)
      return true;
  }
  return false;
}

