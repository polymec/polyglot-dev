// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyglot/exodus_file.h"

#if POLYMEC_HAVE_MPI
#include "mpi.h"
#define PARALLEL_AWARE_EXODUS 1
#include "exodusII_par.h"
#endif 

#include "exodusII.h"

struct exodus_file_t 
{
#if POLYMEC_HAVE_MPI
  MPI_Comm comm;      // Parallel communicator.
  MPI_Info mpi_info;  // ???
#endif
  int ex_id;          // Exodus file descriptor/identifier.
  float ex_version;   // Exodus database version number.
  int ex_real_size;   // Word size of data in the present Exodus file.
};

exodus_file_t* exodus_file_new(MPI_Comm comm,
                               const char* filename)
{
  exodus_file_t* file = polymec_malloc(sizeof(exodus_file_t));
#if POLYMEC_HAVE_MPI
  MPI_Info info;
  file->comm = comm;
  int real_size = (int)sizeof(real_t);
  file->ex_id = ex_open_par(filename, EX_WRITE, &real_size,
                            &file->ex_real_size, &file->ex_version, 
                            file->comm, file->mpi_info);
#else
  file->ex_id = ex_open(filename, EX_WRITE, &real_size,
                        &file->ex_real_size, &file->ex_version);
#endif
  if (file->ex_id < 0)
  {
    polymec_free(file);
    file = NULL;
  }
  return file;
}

exodus_file_t* exodus_file_open(MPI_Comm comm,
                                const char* filename)
{
  exodus_file_t* file = polymec_malloc(sizeof(exodus_file_t));
#if POLYMEC_HAVE_MPI
  MPI_Info info;
  file->comm = comm;
  int real_size = (int)sizeof(real_t);
  file->ex_id = ex_open_par(filename, EX_READ, &real_size,
                            &file->ex_real_size, &file->ex_version, 
                            file->comm, file->mpi_info);
#else
  file->ex_id = ex_open(filename, EX_READ, &real_size,
                        &file->ex_real_size, &file->ex_version);
#endif
  if (file->ex_id < 0)
  {
    polymec_free(file);
    file = NULL;
  }
  return file;
}

void exodus_file_close(exodus_file_t* file)
{
  ex_close(file->ex_id);
}

void exodus_file_write_mesh(exodus_file_t* file,
                            mesh_t* mesh)
{
}

mesh_t* exodus_file_read_mesh(exodus_file_t* file)
{
}

int exodus_file_write_time(exodus_file_t* file, real_t time)
{
}

bool exodus_file_next_time(exodus_file_t* file,
                           int* pos,
                           int* time_index,
                           real_t* time)
{
}

void exodus_file_write_scalar_cell_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data)
{
}

real_t* exodus_file_read_scalar_cell_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name)
{
}

void exodus_file_write_cell_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components)
{
}

real_t* exodus_file_read_cell_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components)
{
}

bool exodus_file_contains_cell_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
}

void exodus_file_write_scalar_face_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data)
{
}

real_t* exodus_file_read_scalar_face_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name)
{
}

void exodus_file_write_face_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components)
{
}

real_t* exodus_file_read_face_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components)
{
}

bool exodus_file_contains_face_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
}

void exodus_file_write_scalar_node_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data)
{
}

real_t* exodus_file_read_scalar_node_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name)
{
}

void exodus_file_write_node_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components)
{
}

real_t* exodus_file_read_node_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components)
{
}

bool exodus_file_contains_node_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name)
{
}

