// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "polyglot/fe_mesh.h"

struct fe_block_t 
{
  int num_elem;
  fe_mesh_element_t elem_type;

  fe_mesh_face_t* elem_face_types;
  int num_elem_faces;
  int* elem_face_offsets;
  int* elem_faces;

  int num_elem_nodes;
  int* elem_node_offsets;
  int* elem_nodes;

  int num_elem_edges;
  int* elem_edge_offsets;
  int* elem_edges;
};

fe_block_t* fe_block_new(int num_elements,
                         fe_mesh_element_t type,
                         int* elem_node_indices)
{
  ASSERT(num_elements > 0);
  ASSERT(elem_node_indices != NULL);
  fe_block_t* block = polymec_malloc(sizeof(fe_block_t));
  block->num_elem = num_elements;
  block->elem_type = type;
  return block;
}

fe_block_t* fe_polyhedral_block_new(int num_elements,
                                    int* num_elem_faces,
                                    fe_mesh_face_t* elem_face_types,
                                    int* face_node_indices)
{
  ASSERT(num_elements > 0);
  ASSERT(num_elem_faces != NULL);
  ASSERT(elem_face_types != NULL);
  ASSERT(face_node_indices != NULL);
  fe_block_t* block = polymec_malloc(sizeof(fe_block_t));
  block->num_elem = num_elements;
  block->elem_type = FE_POLYHEDRON;
  return block;
}

void fe_block_free(fe_block_t* block)
{
  polymec_free(block);
}

fe_block_t* fe_block_clone(fe_block_t* block)
{
  fe_block_t* copy = polymec_malloc(sizeof(fe_block_t));
  copy->num_elem = block->num_elem;
  copy->elem_type = block->elem_type;
  return copy;
}

fe_mesh_element_t fe_block_element_type(fe_block_t* block)
{
  return block->elem_type;
}

int fe_block_num_elements(fe_block_t* block)
{
  return block->num_elem;
}

int fe_block_num_element_nodes(fe_block_t* block, int elem_index)
{
  int offset = block->elem_node_offsets[elem_index];
  return block->elem_node_offsets[elem_index+1] - offset;
}

void fe_block_get_element_nodes(fe_block_t* block, 
                                int elem_index, 
                                int* elem_nodes)
{
  int offset = block->elem_node_offsets[elem_index];
  int num_nodes = block->elem_node_offsets[elem_index+1] - offset;
  memcpy(elem_nodes, &block->elem_nodes[offset], sizeof(int) * num_nodes);
}

int fe_block_num_element_faces(fe_block_t* block, int elem_index)
{
  int offset = block->elem_face_offsets[elem_index];
  return block->elem_face_offsets[elem_index+1] - offset;
}

void fe_block_get_element_faces(fe_block_t* block, 
                                int elem_index, 
                                int* elem_faces)
{
  int offset = block->elem_face_offsets[elem_index];
  int num_faces = block->elem_face_offsets[elem_index+1] - offset;
  memcpy(elem_faces, &block->elem_faces[offset], sizeof(int) * num_faces);
}

int fe_block_num_element_edges(fe_block_t* block, int elem_index)
{
  int offset = block->elem_edge_offsets[elem_index];
  return block->elem_edge_offsets[elem_index+1] - offset;
}

void fe_block_get_element_edges(fe_block_t* block, 
                                int elem_index, 
                                int* elem_edges)
{
  int offset = block->elem_edge_offsets[elem_index];
  int num_faces = block->elem_edge_offsets[elem_index+1] - offset;
  memcpy(elem_edges, &block->elem_edges[offset], sizeof(int) * num_faces);
}

serializer_t* fe_block_serializer()
{
  return serializer_new("fe_block", NULL, NULL, NULL, NULL);
}

struct fe_mesh_t 
{
  MPI_Comm comm;
  ptr_array_t* blocks;
  string_array_t* block_names;

  // mesh -> block element index mapping.
  int_array_t* block_elem_offsets;

  // Nodal coordinates.
  int num_nodes;
  point_t* node_coords;

  // Face-related connectivity.
  int num_faces;
  int* face_edge_offsets;
  int* face_edges;
  int* face_node_offsets;
  int* face_nodes;

  // Edge-related connectivity.
  int num_edges;
  int* edge_node_offsets;
  int* edge_nodes;
};

fe_mesh_t* fe_mesh_new(MPI_Comm comm, int num_nodes)
{
  ASSERT(num_nodes >= 4);
  fe_mesh_t* mesh = polymec_malloc(sizeof(fe_mesh_t));
  mesh->comm = comm;
  mesh->num_nodes = num_nodes;
  mesh->blocks = ptr_array_new();
  mesh->block_names = string_array_new();
  mesh->block_elem_offsets = int_array_new();
  int_array_append(mesh->block_elem_offsets, 0);
  mesh->node_coords = polymec_malloc(sizeof(point_t) * mesh->num_nodes);
  return mesh;
}

void fe_mesh_free(fe_mesh_t* mesh)
{
  ptr_array_free(mesh->blocks);
  string_array_free(mesh->block_names);
  int_array_free(mesh->block_elem_offsets);
  polymec_free(mesh->node_coords);
  polymec_free(mesh);
}

fe_mesh_t* fe_mesh_clone(fe_mesh_t* mesh)
{
  fe_mesh_t* copy = polymec_malloc(sizeof(fe_mesh_t));
  copy->comm = mesh->comm;
  copy->num_nodes = mesh->num_nodes;
  copy->blocks = ptr_array_new();
  for (int i = 0; i < mesh->blocks->size; ++i)
    copy->blocks->data[i] = fe_block_clone(mesh->blocks->data[i]);
  copy->block_names = string_array_new();
  for (int i = 0; i < mesh->block_names->size; ++i)
    copy->block_names->data[i] = string_dup(mesh->block_names->data[i]);
  copy->block_elem_offsets = int_array_new();
  for (int i = 0; i < mesh->block_elem_offsets->size; ++i)
    copy->block_elem_offsets->data[i] = mesh->block_elem_offsets->data[i];
  copy->node_coords = polymec_malloc(sizeof(point_t) * copy->num_nodes);
  memcpy(copy->node_coords, mesh->node_coords, sizeof(point_t) * copy->num_nodes);
  return copy;
}

void fe_mesh_add_block(fe_mesh_t* mesh, 
                       const char* name,
                       fe_block_t* block)
{
  ptr_array_append_with_dtor(mesh->blocks, block, DTOR(fe_block_free));
  string_array_append_with_dtor(mesh->block_names, string_dup(name), string_free);
  int num_elements = mesh->block_elem_offsets->data[mesh->block_elem_offsets->size] + fe_block_num_elements(block);
  int_array_append(mesh->block_elem_offsets, num_elements);
}

int fe_mesh_num_blocks(fe_mesh_t* mesh)
{
  return mesh->blocks->size;
}

bool fe_mesh_next_block(fe_mesh_t* mesh, 
                        int* pos, 
                        char** block_name, 
                        fe_block_t** block)
{
  if (*pos >= mesh->blocks->size)
    return false;

  *block = mesh->blocks->data[*pos];
  *block_name = mesh->block_names->data[*pos];
  ++(*pos);
  return true;
}

int fe_mesh_num_elements(fe_mesh_t* mesh)
{
  return mesh->block_elem_offsets->data[mesh->block_elem_offsets->size-1];
}

int fe_mesh_num_element_nodes(fe_mesh_t* mesh, int elem_index)
{
  ASSERT(elem_index >= 0);
  ASSERT(elem_index < fe_mesh_num_elements(mesh));

  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  return fe_block_num_element_nodes(block, e);
}

void fe_mesh_get_element_nodes(fe_mesh_t* mesh, 
                               int elem_index, 
                               int* elem_nodes)
{
  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  fe_block_get_element_nodes(block, e, elem_nodes);
}

int fe_mesh_num_element_faces(fe_mesh_t* mesh, int elem_index)
{
  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  return fe_block_num_element_faces(block, e);
}

void fe_mesh_get_element_faces(fe_mesh_t* mesh, 
                               int elem_index, 
                               int* elem_faces)
{
  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  fe_block_get_element_faces(block, e, elem_faces);
}

int fe_mesh_num_element_edges(fe_mesh_t* mesh, int elem_index)
{
  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  return fe_block_num_element_edges(block, e);
}

void fe_mesh_get_element_edges(fe_mesh_t* mesh, 
                               int elem_index, 
                               int* elem_edges)
{
  // Find the block that houses this element.
  int b = 0;
  while (mesh->block_elem_offsets->data[b+1] < elem_index)
    ++b;

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  fe_block_get_element_edges(block, e, elem_edges);
}

int fe_mesh_num_faces(fe_mesh_t* mesh)
{
  return mesh->num_faces;
}

int fe_mesh_num_face_edges(fe_mesh_t* mesh,
                           int face_index)
{
  int offset = mesh->face_edge_offsets[face_index];
  return mesh->face_edge_offsets[face_index+1] - offset;
}

void fe_mesh_get_face_edges(fe_mesh_t* mesh, 
                            int face_index, 
                            int* face_edges)
{
  int offset = mesh->face_edge_offsets[face_index];
  int num_edges = mesh->face_edge_offsets[face_index+1] - offset;
  memcpy(face_edges, &mesh->face_edges[offset], sizeof(int) * num_edges);
}

int fe_mesh_num_face_nodes(fe_mesh_t* mesh,
                           int face_index)
{
  int offset = mesh->face_node_offsets[face_index];
  return mesh->face_node_offsets[face_index+1] - offset;
}

void fe_mesh_get_face_nodes(fe_mesh_t* mesh, 
                            int face_index, 
                            int* face_nodes)
{
  int offset = mesh->face_node_offsets[face_index];
  int num_nodes = mesh->face_node_offsets[face_index+1] - offset;
  memcpy(face_nodes, &mesh->face_nodes[offset], sizeof(int) * num_nodes);
}

int fe_mesh_num_edges(fe_mesh_t* mesh)
{
  return mesh->num_edges;
}

int fe_mesh_num_edge_nodes(fe_mesh_t* mesh,
                           int edge_index)
{
  int offset = mesh->edge_node_offsets[edge_index];
  return mesh->edge_node_offsets[edge_index+1] - offset;
}

void fe_mesh_get_edge_nodes(fe_mesh_t* mesh, 
                            int edge_index, 
                            int* edge_nodes)
{
  int offset = mesh->edge_node_offsets[edge_index];
  int num_nodes = mesh->edge_node_offsets[edge_index+1] - offset;
  memcpy(edge_nodes, &mesh->edge_nodes[offset], sizeof(int) * num_nodes);
}

int fe_mesh_num_nodes(fe_mesh_t* mesh)
{
  return mesh->num_nodes;
}

point_t* fe_mesh_node_coordinates(fe_mesh_t* mesh)
{
  return mesh->node_coords;
}

serializer_t* fe_mesh_serializer()
{
  return serializer_new("fe_mesh", NULL, NULL, NULL, NULL);
}

