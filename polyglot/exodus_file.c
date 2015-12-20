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
  MPI_Comm comm;        // Parallel communicator.
  MPI_Info mpi_info;    // ???
#endif
  int ex_id;            // Exodus file descriptor/identifier.
  float ex_version;     // Exodus database version number.
  int ex_real_size;     // Word size of data in the present Exodus file.
  int last_time_index;  // Index of most-recently added time written to file.
};

exodus_file_t* exodus_file_new(MPI_Comm comm,
                               const char* filename)
{
  exodus_file_t* file = polymec_malloc(sizeof(exodus_file_t));
  file->last_time_index = 0;
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
  file->last_time_index = 0;
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
  mesh_t* mesh = NULL;

  // Get information from the file.
  char title[MAX_LINE_LENGTH];
  ex_init_params mesh_info;
  int status = ex_get_init_ext(file->ex_id, &mesh_info);
  if ((status >= 0) && (mesh_info.num_dim == 3))
  {
    int num_nodes = mesh_info.num_nodes;
    int num_elem_blocks = mesh_info.num_elem_blk;
    int num_face_blocks = mesh_info.num_face_blk;
    int num_edge_blocks = mesh_info.num_edge_blk;

    // Count the number of cells in the mesh.
    int num_cells = 0, num_ghost_cells = 0;
    // FIXME: Need to get ghost cells somehow?
    int elem_counts[num_elem_blocks];
    char elem_types[num_elem_blocks][MAX_NAME_LENGTH];
    int num_nodes_per_elem[num_elem_blocks], 
        num_faces_per_elem[num_elem_blocks], 
        num_attr_per_elem[num_elem_blocks];
    for (int elem_block = 1; elem_block <= num_elem_blocks; ++elem_block)
    {
      int block_index = elem_block-1;
      int status = ex_get_block(file->ex_id, EX_ELEM_BLOCK, elem_block, 
                                elem_types[block_index], 
                                &elem_counts[block_index], 
                                &num_nodes_per_elem[block_index], NULL,
                                &num_faces_per_elem[block_index],
                                &num_attr_per_elem[block_index]);
      if (status < 0)
        return NULL;
      num_cells += elem_counts[block_index];
    }

    // Count the number of faces in the mesh.
    int num_faces = 0;
    int face_counts[num_face_blocks];
    char face_types[num_face_blocks][MAX_NAME_LENGTH];
    int num_nodes_per_face[num_face_blocks], 
        num_attr_per_face[num_face_blocks];
    for (int face_block = 1; face_block <= num_face_blocks; ++face_block)
    {
      int block_index = face_block-1;
      int status = ex_get_block(file->ex_id, EX_FACE_BLOCK, face_block, 
                                face_types[block_index], 
                                &face_counts[block_index], 
                                &num_nodes_per_face[block_index], 
                                NULL, NULL, &num_attr_per_face[block_index]);
      if (status < 0)
        return NULL;
      num_faces += face_counts[face_block-1];
    }

    // Now create the mesh and allocate space for its connectivity.
    mesh = mesh_new(file->comm, num_cells, num_ghost_cells, num_faces, num_nodes);
    int cell_index = 0;
    mesh->cell_face_offsets[0] = 0;
    for (int elem_block = 1; elem_block <= num_elem_blocks; ++elem_block)
    {
      int block_index = elem_block - 1;
      int num_elem = elem_counts[block_index];
      if (string_casecmp(elem_types[block_index], "nfaced") == 0)
      {
        // The cells in this block are polyhedral.
        int elem_face_count[num_elem];
        int status = ex_get_entity_count_per_polyhedra(file->ex_id, EX_ELEM_BLOCK,
                                                       elem_block, elem_face_count);
        if (status < 0)
          return NULL;

        // Fill in the cell-face offsets.
        for (int i = cell_index; i < cell_index + num_elem; ++i)
          mesh->cell_face_offsets[i+1] = mesh->cell_face_offsets[i] + elem_face_count[i];
      }
      else
      {
        // Fill in the cell-face offsets.
        for (int i = cell_index; i < cell_index + num_elem; ++i)
          mesh->cell_face_offsets[i+1] = mesh->cell_face_offsets[i] + num_faces_per_elem[block_index];
      }
      cell_index += num_elem;
    }

    int face_index = 0;
    mesh->face_node_offsets[0] = 0;
    for (int face_block = 1; face_block <= num_face_blocks; ++face_block)
    {
      int block_index = face_block - 1;
      if (string_casecmp(face_types[block_index], "nsided") == 0)
      {
        // The faces in this block are polygonal.
        int num_faces = face_counts[block_index];
        int face_node_count[num_faces];
        int status = ex_get_entity_count_per_polyhedra(file->ex_id, EX_FACE_BLOCK,
                                                       face_block, face_node_count);
        if (status < 0)
          return NULL;

        // Fill in the face-node offsets.
        for (int i = face_index; i < face_index + num_faces; ++i)
          mesh->face_node_offsets[i+1] = mesh->face_node_offsets[i] + face_node_count[i];
      }
      else
      {
        // Fill in the face-node offsets.
        for (int i = face_index; i < face_index + num_faces; ++i)
          mesh->face_node_offsets[i+1] = mesh->face_node_offsets[i] + num_nodes_per_face[block_index];
      }
      face_index += num_faces;
    }
    mesh_reserve_connectivity_storage(mesh);

    // Now fill in the connectivity.
    for (int elem_block = 1; elem_block <= num_elem_blocks; ++elem_block)
    {
      // Fetch the connectivity information. Search first for elem->face 
      // connectivity.
      int block_index = elem_block - 1;
      int offset = mesh->cell_face_offsets[block_index];
      int tot_num_cell_faces = mesh->cell_face_offsets[block_index+1] - offset;
      int face_conn[tot_num_cell_faces];
      int status = ex_get_conn(file->ex_id, EX_ELEM_BLOCK, elem_block, NULL, NULL, face_conn);
      if (status >= 0)
      {
        // Great! We have faces for cells. Copy them into place.
        memcpy(&mesh->cell_faces[offset], face_conn, sizeof(int) * tot_num_cell_faces);
      }
      else
      {
        // We don't have cell->face connectivity, but we do have cell->node 
        // connectivity. I guess we have to put on our finite element hat
        // and construct the faces manually.
        polymec_error("exodus_file_read_mesh: Nodal finite element meshes are not yet supported.");
      }
    }
    
    // Face->node connectivity.
    for (int face_block = 1; face_block <= num_face_blocks; ++face_block)
    {
      int block_index = face_block - 1;
      int offset = mesh->face_node_offsets[block_index];
      int tot_num_face_nodes = mesh->face_node_offsets[block_index+1] - offset;
      int node_conn[tot_num_face_nodes];
      int status = ex_get_conn(file->ex_id, EX_FACE_BLOCK, face_block, node_conn, NULL, NULL);
      memcpy(&mesh->face_nodes[offset], node_conn, sizeof(int) * tot_num_face_nodes);
    }

    // Edge->node connectivity.
    if (num_edge_blocks > 0)
    {
      for (int edge_block = 1; edge_block <= num_edge_blocks; ++edge_block)
      {
        int block_index = edge_block - 1;
        int offset = 0; // FIXME
        int tot_num_edge_nodes = 0; // FIXME
        int node_conn[tot_num_edge_nodes];
        int status = ex_get_conn(file->ex_id, EX_EDGE_BLOCK, edge_block, node_conn, NULL, NULL);
        memcpy(&mesh->edge_nodes[offset], node_conn, sizeof(int) * tot_num_edge_nodes);
      }
    }
    else
    {
      // No edge information was found in the mesh, so we create our own.
      mesh_construct_edges(mesh);
    }

    // Fetch node positions and compute geometry.
    real_t x[num_nodes], y[num_nodes], z[num_nodes];
    ex_get_coord(file->ex_id, x, y, z);
    for (int n = 0; n < mesh->num_nodes; ++n)
    {
      mesh->nodes[n].x = x[n];
      mesh->nodes[n].y = y[n];
      mesh->nodes[n].z = z[n];
    }
    mesh_compute_geometry(mesh);

    return mesh;
  }
  else
    return NULL;
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

