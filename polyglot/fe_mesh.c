// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/array_utils.h"
#include "core/tagger.h"
#include "polyglot/fe_mesh.h"

struct fe_block_t 
{
  int num_elem;
  fe_mesh_element_t elem_type;

  int* elem_face_offsets;
  int* elem_faces;

  int* elem_node_offsets;
  int* elem_nodes;

  int* elem_edge_offsets;
  int* elem_edges;
};

static int get_num_elem_nodes(fe_mesh_element_t elem_type)
{
  switch(elem_type)
  {
    case FE_TETRAHEDRON_4: return 4;
    case FE_PYRAMID_5: return 5;
    case FE_WEDGE_6: return 6;
    case FE_TETRAHEDRON_8:
    case FE_HEXAHEDRON_8: return 8;
    case FE_HEXAHEDRON_9: return 9;
    case FE_TETRAHEDRON_10: return 10;
    case FE_PYRAMID_13: return 13;
    case FE_TETRAHEDRON_14: return 14;
    case FE_WEDGE_15: return 15;
    case FE_HEXAHEDRON_20: return 20;
    case FE_HEXAHEDRON_27: return 27;
    default: return -1;
  }
}

fe_block_t* fe_block_new(int num_elem,
                         fe_mesh_element_t type,
                         int* elem_node_indices)
{
  ASSERT(num_elem > 0);
  ASSERT(elem_node_indices != NULL);
  fe_block_t* block = polymec_malloc(sizeof(fe_block_t));
  block->num_elem = num_elem;
  block->elem_type = type;

  // Element nodes.
  int num_elem_nodes = get_num_elem_nodes(type);
  block->elem_node_offsets = polymec_malloc(sizeof(int) * num_elem_nodes * num_elem);
  block->elem_node_offsets[0] = 0;
  for (int i = 0; i < num_elem; ++i)
    block->elem_node_offsets[i+1] = block->elem_node_offsets[i] + num_elem_nodes;
  block->elem_nodes = polymec_malloc(sizeof(int) * block->elem_node_offsets[num_elem]);
  memcpy(block->elem_nodes, elem_node_indices, sizeof(int) * block->elem_node_offsets[num_elem]);

  // Elements don't understand their faces.
  block->elem_face_offsets = NULL;
  block->elem_faces = NULL;

  block->elem_edge_offsets = NULL;
  block->elem_edges = NULL;

  return block;
}

fe_block_t* polyhedral_fe_block_new(int num_elem,
                                    int* num_elem_faces,
                                    int* elem_face_indices)
{
  ASSERT(num_elem > 0);
  ASSERT(num_elem_faces != NULL);
  ASSERT(elem_face_indices != NULL);
  fe_block_t* block = polymec_malloc(sizeof(fe_block_t));
  block->num_elem = num_elem;
  block->elem_type = FE_POLYHEDRON;

  // Element faces.
  int tot_elem_faces = 0;
  for (int i = 0; i < num_elem; ++i)
    tot_elem_faces += num_elem_faces[i];
  block->elem_face_offsets = polymec_malloc(sizeof(int) * tot_elem_faces);
  block->elem_face_offsets[0] = 0;
  for (int i = 0; i < num_elem; ++i)
    block->elem_face_offsets[i+1] = block->elem_face_offsets[i] + num_elem_faces[i];
  block->elem_faces = polymec_malloc(sizeof(int) * block->elem_face_offsets[num_elem]);
  memcpy(block->elem_faces, elem_face_indices, sizeof(int) * block->elem_face_offsets[num_elem]);

  // Element nodes are not determined until the block is added to the mesh.
  block->elem_node_offsets = NULL;
  block->elem_nodes = NULL;

  block->elem_edge_offsets = NULL;
  block->elem_edges = NULL;

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

  // Nodal positions.
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

  // Entity sets.
  tagger_t* elem_sets;
  tagger_t* face_sets;
  tagger_t* edge_sets;
  tagger_t* node_sets;
  tagger_t* side_sets;
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
  memset(mesh->node_coords, 0, sizeof(point_t) * mesh->num_nodes);

  mesh->num_faces = 0;
  mesh->face_node_offsets = NULL;
  mesh->face_nodes = NULL;
  mesh->face_edge_offsets = NULL;
  mesh->face_edges = NULL;

  mesh->num_edges = 0;
  mesh->edge_node_offsets = NULL;
  mesh->edge_nodes = NULL;

  mesh->elem_sets = tagger_new();
  mesh->face_sets = tagger_new();
  mesh->edge_sets = tagger_new();
  mesh->node_sets = tagger_new();
  mesh->side_sets = tagger_new();

  return mesh;
}

void fe_mesh_free(fe_mesh_t* mesh)
{
  tagger_free(mesh->elem_sets);
  tagger_free(mesh->face_sets);
  tagger_free(mesh->edge_sets);
  tagger_free(mesh->node_sets);
  tagger_free(mesh->side_sets);

  if (mesh->face_nodes != NULL)
  {
    polymec_free(mesh->face_nodes);
    polymec_free(mesh->face_node_offsets);
  }

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

MPI_Comm fe_mesh_comm(fe_mesh_t* mesh)
{
  return mesh->comm;
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
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return -1;
  }

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
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return;
  }

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  fe_block_get_element_nodes(block, e, elem_nodes);
}

int fe_mesh_num_element_faces(fe_mesh_t* mesh, int elem_index)
{
  // Find the block that houses this element.
  int b = 0;
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return -1;
  }

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
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return;
  }

  // Now ask the block about the element.
  fe_block_t* block = mesh->blocks->data[b];
  int e = elem_index - mesh->block_elem_offsets->data[b];
  fe_block_get_element_faces(block, e, elem_faces);
}

int fe_mesh_num_element_edges(fe_mesh_t* mesh, int elem_index)
{
  // Find the block that houses this element.
  int b = 0;
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return -1;
  }

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
  if (mesh->block_elem_offsets->size == 2)
  {
    while (((b+1) < mesh->block_elem_offsets->size) && 
           (mesh->block_elem_offsets->data[b+1] < elem_index)) ++b;
    if (b == mesh->block_elem_offsets->size-1)
      return;
  }

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

void fe_mesh_set_face_nodes(fe_mesh_t* mesh, 
                            int num_faces,
                            int* num_face_nodes, 
                            int* face_nodes)
{
  ASSERT(num_faces > 0);
  mesh->num_faces = num_faces;
  mesh->face_node_offsets = polymec_malloc(sizeof(int) * (mesh->num_faces+1));
  mesh->face_node_offsets[0] = 0;
  for (int i = 0; i < num_faces; ++i)
    mesh->face_node_offsets[i+1] = mesh->face_node_offsets[i] + num_face_nodes[i];
  mesh->face_nodes = polymec_malloc(sizeof(int) * mesh->face_node_offsets[num_faces]);
  memcpy(mesh->face_nodes, face_nodes, sizeof(int) * mesh->face_node_offsets[num_faces]);
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

point_t* fe_mesh_node_positions(fe_mesh_t* mesh)
{
  return mesh->node_coords;
}

int fe_mesh_num_element_sets(fe_mesh_t* mesh)
{
  // Count up the tags in the appropriate tagger.
  int pos = 0, *tag, size, num_tags = 0;
  char* tag_name;
  while (tagger_next_tag(mesh->elem_sets, &pos, &tag_name, &tag, &size))
    ++num_tags;
  return num_tags;
}

int* fe_mesh_create_element_set(fe_mesh_t* mesh, const char* name, int size)
{
  return tagger_create_tag(mesh->elem_sets, name, size);
}

bool fe_mesh_next_element_set(fe_mesh_t* mesh, int* pos, char** name, int** set, int* size)
{
  return tagger_next_tag(mesh->elem_sets, pos, name, set, size);
}

int fe_mesh_num_face_sets(fe_mesh_t* mesh)
{
  // Count up the tags in the appropriate tagger.
  int pos = 0, *tag, size, num_tags = 0;
  char* tag_name;
  while (tagger_next_tag(mesh->face_sets, &pos, &tag_name, &tag, &size))
    ++num_tags;
  return num_tags;
}

int* fe_mesh_create_face_set(fe_mesh_t* mesh, const char* name, int size)
{
  return tagger_create_tag(mesh->face_sets, name, size);
}

bool fe_mesh_next_face_set(fe_mesh_t* mesh, int* pos, char** name, int** set, int* size)
{
  return tagger_next_tag(mesh->face_sets, pos, name, set, size);
}

int fe_mesh_num_edge_sets(fe_mesh_t* mesh)
{
  // Count up the tags in the appropriate tagger.
  int pos = 0, *tag, size, num_tags = 0;
  char* tag_name;
  while (tagger_next_tag(mesh->edge_sets, &pos, &tag_name, &tag, &size))
    ++num_tags;
  return num_tags;
}

int* fe_mesh_create_edge_set(fe_mesh_t* mesh, const char* name, int size)
{
  return tagger_create_tag(mesh->edge_sets, name, size);
}

bool fe_mesh_next_edge_set(fe_mesh_t* mesh, int* pos, char** name, int** set, int* size)
{
  return tagger_next_tag(mesh->edge_sets, pos, name, set, size);
}

int fe_mesh_num_node_sets(fe_mesh_t* mesh)
{
  // Count up the tags in the appropriate tagger.
  int pos = 0, *tag, size, num_tags = 0;
  char* tag_name;
  while (tagger_next_tag(mesh->node_sets, &pos, &tag_name, &tag, &size))
    ++num_tags;
  return num_tags;
}

int* fe_mesh_create_node_set(fe_mesh_t* mesh, const char* name, int size)
{
  return tagger_create_tag(mesh->node_sets, name, size);
}

bool fe_mesh_next_node_set(fe_mesh_t* mesh, int* pos, char** name, int** set, int* size)
{
  return tagger_next_tag(mesh->node_sets, pos, name, set, size);
}

int fe_mesh_num_side_sets(fe_mesh_t* mesh)
{
  // Count up the tags in the appropriate tagger.
  int pos = 0, *tag, size, num_tags = 0;
  char* tag_name;
  while (tagger_next_tag(mesh->side_sets, &pos, &tag_name, &tag, &size))
    ++num_tags;
  return num_tags;
}

int* fe_mesh_create_side_set(fe_mesh_t* mesh, const char* name, int size)
{
  return tagger_create_tag(mesh->side_sets, name, 2*size);
}

bool fe_mesh_next_side_set(fe_mesh_t* mesh, int* pos, char** name, int** set, int* size)
{
  return tagger_next_tag(mesh->side_sets, pos, name, set, size);
}

serializer_t* fe_mesh_serializer()
{
  return serializer_new("fe_mesh", NULL, NULL, NULL, NULL);
}

mesh_t* mesh_from_fe_mesh(fe_mesh_t* fe_mesh)
{
  // Feel out the faces for the finite element mesh. Do we need to create 
  // them ourselves, or are they already all there?
  int num_cells = fe_mesh_num_elements(fe_mesh);
  int num_faces = fe_mesh_num_faces(fe_mesh);
  int* cell_face_offsets = NULL;
  int* cell_faces = NULL;
  int* face_node_offsets = NULL;
  int* face_nodes = NULL;
  if (num_faces == 0)
  {
    // This is a nodal finite element mesh, so we need to assemble the faces ourselves.
    polymec_not_implemented("nodal finite element mesh conversion");

    // Traverse the element blocks and create faces.
    int pos = 0, elem_offset = 0;
    char* block_name;
    fe_block_t* block;
    while (fe_mesh_next_block(fe_mesh, &pos, &block_name, &block))
    {
      int num_block_elem = fe_block_num_elements(block);
      fe_mesh_element_t elem_type = fe_block_element_type(block);

      elem_offset += num_block_elem;
    }

    // Now identify the nodes for each face.
  }
  else
  {
    // Fill in these arrays block by block.
    cell_face_offsets = polymec_malloc(sizeof(int) * (num_cells + 1));
    cell_face_offsets[0] = 0;
    int pos = 0;
    char* block_name;
    fe_block_t* block;
    int block_cell_offset = 0;
    while (fe_mesh_next_block(fe_mesh, &pos, &block_name, &block))
    {
      int num_block_elem = fe_block_num_elements(block);
      for (int i = 0; i < num_block_elem; ++i)
        cell_face_offsets[block_cell_offset+i] = cell_face_offsets[block_cell_offset+i-1] + block->elem_face_offsets[i];
      block_cell_offset += num_block_elem;
    }

    cell_faces = polymec_malloc(sizeof(int) * cell_face_offsets[num_cells]);
    pos = 0, block_cell_offset = 0;
    while (fe_mesh_next_block(fe_mesh, &pos, &block_name, &block))
    {
      int num_block_elem = fe_block_num_elements(block);
      memcpy(&cell_faces[cell_face_offsets[block_cell_offset]], block->elem_faces, sizeof(int) * block->elem_face_offsets[num_block_elem]);
      block_cell_offset += num_block_elem;
    }

    // We just borrow these from the mesh. Theeenks!
    face_node_offsets = fe_mesh->face_node_offsets;
    face_nodes = fe_mesh->face_nodes;
  }
  ASSERT(cell_face_offsets != NULL);
  ASSERT(cell_faces != NULL);
  ASSERT(face_node_offsets != NULL);
  ASSERT(face_nodes != NULL);

  int num_ghost_cells = 0; // FIXME!
  mesh_t* mesh = mesh_new(fe_mesh_comm(fe_mesh), 
                          num_cells, num_ghost_cells, 
                          fe_mesh_num_faces(fe_mesh),
                          fe_mesh_num_nodes(fe_mesh));
  memcpy(mesh->cell_face_offsets, cell_face_offsets, sizeof(int) * (mesh->num_cells+1));
  memcpy(mesh->face_node_offsets, face_node_offsets, sizeof(int) * (mesh->num_faces+1));
  mesh_reserve_connectivity_storage(mesh);
  memcpy(mesh->cell_faces, cell_faces, sizeof(int) * (mesh->cell_face_offsets[mesh->num_cells]));
  memcpy(mesh->face_nodes, face_nodes, sizeof(int) * (mesh->face_node_offsets[mesh->num_faces]));

  // Set up face->cell connectivity.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    for (int f = mesh->cell_face_offsets[c]; f < mesh->cell_face_offsets[c+1]; ++f)
    {
      int face = mesh->cell_faces[f];
      if (mesh->face_cells[2*face] == -1)
        mesh->face_cells[2*face] = c;
      else
        mesh->face_cells[2*face+1] = c;
    }
  }

  // Set up face->edge connectivity and edge->node connectivity (if provided).
  if (fe_mesh->face_edges != NULL)
  {
    memcpy(mesh->face_edge_offsets, fe_mesh->face_edge_offsets, sizeof(int) * (mesh->num_faces+1));
    mesh->face_edges = polymec_malloc(sizeof(int) * fe_mesh->face_edge_offsets[mesh->num_faces]);
    memcpy(mesh->face_edges, fe_mesh->face_edges, sizeof(int) * fe_mesh->face_edge_offsets[mesh->num_faces]);
  }
  else
  {
    // Construct edges if we didn't find them.
    mesh_construct_edges(mesh);
  }

  // Copy the node positions into place.
  memcpy(mesh->nodes, fe_mesh_node_positions(fe_mesh), sizeof(point_t) * mesh->num_nodes);

  // Calculate geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  polymec_free(cell_face_offsets);
  polymec_free(cell_faces);
  if (fe_mesh_num_faces(fe_mesh) == 0)
  {
    polymec_free(face_node_offsets);
    polymec_free(face_nodes);
  }

  return mesh;
}

fe_mesh_t* fe_mesh_from_mesh(mesh_t* fv_mesh,
                             string_array_t* element_block_tags)
{
  fe_mesh_t* fe_mesh = fe_mesh_new(fv_mesh->comm, fv_mesh->num_nodes);

  if ((element_block_tags != NULL) && (element_block_tags->size > 1))
  {
    // Block-by-block construction.
    for (int b = 0; b < element_block_tags->size; ++b)
    {
      char* tag_name = element_block_tags->data[b];
      int num_elem;
      int* block_tag = mesh_tag(fv_mesh->cell_tags, tag_name, &num_elem);
      int* num_elem_faces = polymec_malloc(sizeof(int) * num_elem);
      int elem_faces_size = 0;
      for (int i = 0; i < num_elem; ++i)
      {
        int nef = fv_mesh->cell_face_offsets[block_tag[i]+1] - fv_mesh->cell_face_offsets[block_tag[i]];
        num_elem_faces[i] = nef;
        elem_faces_size += nef;
      }
      int* elem_faces = polymec_malloc(sizeof(int) * elem_faces_size);
      int offset = 0;
      for (int i = 0; i < num_elem; ++i)
      {
        for (int f = 0; f < num_elem_faces[i]; ++f, ++offset)
          elem_faces[offset+f] = fv_mesh->cell_faces[fv_mesh->cell_face_offsets[block_tag[i]+f]];
      }
      fe_block_t* block = polyhedral_fe_block_new(num_elem, num_elem_faces, fv_mesh->cell_faces);
      fe_mesh_add_block(fe_mesh, tag_name, block);
      polymec_free(num_elem_faces);
      polymec_free(elem_faces);
    }
  }
  else
  {
    // One honking block.
    int num_elem = fv_mesh->num_cells;
    int* num_elem_faces = polymec_malloc(sizeof(int) * num_elem);
    for (int i = 0; i < num_elem; ++i)
      num_elem_faces[i] = fv_mesh->cell_face_offsets[i+1] - fv_mesh->cell_face_offsets[i];
    fe_block_t* block = polyhedral_fe_block_new(num_elem, num_elem_faces, fv_mesh->cell_faces);
    fe_mesh_add_block(fe_mesh, "block_1", block);
    polymec_free(num_elem_faces);
  }

  // Copy coordinates.
  memcpy(fe_mesh_node_positions(fe_mesh), fv_mesh->nodes, sizeof(point_t) * fv_mesh->num_nodes);

  return fe_mesh;
}

