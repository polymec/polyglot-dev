// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYGLOT_EXODUS_FILE_H
#define POLYGLOT_EXODUS_FILE_H

#include "polyglot/fe_mesh.h"

// The Exodus file class provides an interface for reading and writing Exodus II
// files which are NetCDF files that follow the Exodus II finite element 
// conventions.

// This type provides the interface for Exodus II files.
typedef struct exodus_file_t exodus_file_t;

// Queries the Exodus file with the given name, fetching the size of real numbers in the file, 
// (floating point) version number for the specification, the number of MPI processes for which 
// it has data, and (if times is non-NULL) an array of times for which the file contains data.
// Returns true if the given file is a valid Exodus file, false if it is not.
bool exodus_file_query(const char* filename,
                       int* real_size,
                       float* version,
                       int* num_mpi_processes,
                       real_array_t* times);

// Creates and opens a new Exodus file for writing simulation data, 
// returning the Exodus file object. 
exodus_file_t* exodus_file_new(MPI_Comm comm, const char* filename);

// Opens an existing Exodus file for reading simulation data,
// returning the Exodus file object. 
exodus_file_t* exodus_file_open(MPI_Comm comm, const char* filename);

// Closes and destroys the given Exodus file.
void exodus_file_close(exodus_file_t* file);

// Writes a finite element mesh to the given Exodus file, overwriting 
// any existing mesh there. All cells (or "elements") are written to a single 
// element block within the Exodus mesh.
void exodus_file_write_fe_mesh(exodus_file_t* file,
                               fe_mesh_t* mesh);

// Reads a finite element mesh from the given Exodus file, returning 
// a newly-allocated object.
fe_mesh_t* exodus_file_read_mesh(exodus_file_t* file);

// Writes a time value to the mesh, returning a newly-created time index 
// that can associate field data to this time.
int exodus_file_write_time(exodus_file_t* file, real_t time);

// Traverse the avaiable times in the Exodus mesh, returning true if the 
// traversal should continue and false if not, and storing the time index 
// and time for the next available set of data in their respective pointers.
// Set *pos to 0 to reset the iteration.
bool exodus_file_next_time(exodus_file_t* file,
                           int* pos,
                           int* time_index,
                           real_t* time);

// Writes a named scalar cell-centered field to the given Exodus file, 
// associated it the time identified by the given time index.
void exodus_file_write_scalar_cell_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data);

// Reads a named scalar cell-centered field from the Exodus file, 
// returning a newly-allocated array of field data associated with the time 
// for the given time index.
real_t* exodus_file_read_scalar_cell_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name);

// Writes a named multicomponent cell-centered field, to the given Exodus file. 
// The field data is interpreted to be in component-minor order and is associated
// with the time identified by the given time index..
void exodus_file_write_cell_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components);

// Reads a named multicomponent cell-centered field from the Exodus file, returning a 
// newly-allocated array of field data for the given time index.
real_t* exodus_file_read_cell_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components);

// Returns true if the given Exodus file contains a (scalar) cell-centered field 
// with the given name and associated with the given time index, false otherwise.
bool exodus_file_contains_cell_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name);

// Writes a named scalar face-centered field to the given Exodus file, 
// associated it the time identified by the given time index.
void exodus_file_write_scalar_face_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data);

// Reads a named scalar face-centered field from the Exodus file, 
// returning a newly-allocated array of field data associated with the time 
// for the given time index.
real_t* exodus_file_read_scalar_face_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name);

// Writes a named multicomponent face-centered field, to the given Exodus file. 
// The field data is interpreted to be in component-minor order and is associated
// with the time identified by the given time index..
void exodus_file_write_face_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components);

// Reads a named multicomponent face-centered field from the Exodus file, returning a 
// newly-allocated array of field data for the given time index.
real_t* exodus_file_read_face_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components);

// Returns true if the given Exodus file contains a (scalar) face-centered field 
// with the given name and associated with the given time index, false otherwise.
bool exodus_file_contains_face_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name);

// Writes a named scalar node-centered field to the given Exodus file, 
// associated it the time identified by the given time index.
void exodus_file_write_scalar_node_field(exodus_file_t* file,
                                         int time_index,
                                         const char* field_name,
                                         real_t* field_data);

// Reads a named scalar node-centered field from the Exodus file, 
// returning a newly-allocated array of field data associated with the time 
// for the given time index.
real_t* exodus_file_read_scalar_node_field(exodus_file_t* file,
                                           int time_index,
                                           const char* field_name);

// Writes a named multicomponent node-centered field, to the given Exodus file. 
// The field data is interpreted to be in component-minor order and is associated
// with the time identified by the given time index..
void exodus_file_write_node_field(exodus_file_t* file,
                                  int time_index,
                                  const char** field_component_names,
                                  real_t* field_data,
                                  int num_components);

// Reads a named multicomponent node-centered field from the Exodus file, returning a 
// newly-allocated array of field data for the given time index.
real_t* exodus_file_read_node_field(exodus_file_t* file,
                                    int time_index,
                                    const char** field_component_names,
                                    int num_components);

// Returns true if the given Exodus file contains a (scalar) node-centered field 
// with the given name and associated with the given time index, false otherwise.
bool exodus_file_contains_node_field(exodus_file_t* file, 
                                     int time_index,
                                     const char* field_name);

#endif
