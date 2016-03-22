// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyglot/import_tetgen_mesh.h"
#include "polyglot/interpreter_register_polyglot_functions.h"
#include "polyglot/exodus_file.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This file contains definitions for functions supported by Polyglot.

// Tetgen mesh factory method.
static int mesh_factory_tetgen(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isstring(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.tetgen(mesh_prefix).");
  }

  // Use the mesh prefix to generate filenames.
  const char* mesh_prefix = lua_tostring(lua, 1);
  char node_file[512], ele_file[512], face_file[512], neigh_file[512];
  snprintf(node_file, 512, "%s.node", mesh_prefix);
  snprintf(ele_file, 512, "%s.ele", mesh_prefix);
  snprintf(face_file, 512, "%s.face", mesh_prefix);
  snprintf(neigh_file, 512, "%s.neigh", mesh_prefix);
  mesh_t* mesh = import_tetgen_mesh(MPI_COMM_WORLD, node_file, ele_file, 
                                    face_file, neigh_file);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

// read_exodus_mesh(args) -- This function reads a mesh from a the given 
// Exodus file on disk. 
static int lua_read_exodus_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isstring(lua, 1))
  {
    return luaL_error(lua, "read_exodus_mesh: invalid arguments. Usage:\n"
                      "mesh = read_exodus_mesh(filename)");
  }

  // Get the argument(s).
  char* filename = string_dup(lua_tostring(lua, 1));

  // Do some checking.
  if (!file_exists(filename))
    return luaL_error(lua, "read_exodus_mesh: file does not exist.");
  size_t real_size;
  float version;
  int num_procs;
  bool valid = exodus_file_query(filename, &real_size, &version,
                                 &num_procs, NULL);
  if (!valid)
    return luaL_error(lua, "read_exodus_mesh: file contains an invalid mesh.");

  // Do our business.
  exodus_file_t* file = exodus_file_open(MPI_COMM_WORLD, filename);
  fe_mesh_t* fe_mesh = exodus_file_read_mesh(file);
  exodus_file_close(file);
  mesh_t* mesh = mesh_from_fe_mesh(fe_mesh);
  fe_mesh_free(fe_mesh);

  lua_pushmesh(lua, mesh);
  return 1;
}

#if 0
int mesh_factory_pebi(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 4))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.pebi(cell_centers, cell_volumes, faces) OR\n"
                      "mesh = mesh_factory.pebi(cell_centers, cell_volumes, faces, face_centers)");
  }
  if (!lua_ispointlist(lua, 1))
    return luaL_error(lua, "cell_centers must be a list of points.");
  if (!lua_issequence(lua, 2))
    return luaL_error(lua, "cell_volumes must be a sequence of cell volumes.");
  if (!lua_istable(lua, 3))
    return luaL_error(lua, "faces must be a table of 3-tuples containing (cell1, cell2, area).");
  if ((num_args == 4) && !lua_ispointlist(lua, 4))
    return luaL_error(lua, "face_centers must be a list of points.");

  // Get the arguments.

  // Cell center point list.
  int num_cells;
  point_t* cell_centers = lua_topointlist(lua, 1, &num_cells);

  // Cell volume list.
  int num_cell_volumes;
  real_t* cell_volumes = lua_tosequence(lua, 2, &num_cell_volumes);
  if (num_cell_volumes != num_cells)
    return luaL_error(lua, "Number of cell volumes (%d) does not match number of cells (%d).", num_cell_volumes, num_cells);

  // Mine the faces table for all its 3-tuples.
  int num_faces = luaL_len(lua, 3);
  ASSERT(num_faces >= 4);
  real_t** faces_table_entries = polymec_malloc(sizeof(real_t*)*num_faces);
  lua_pushnil(lua);
  int face = 0;
  while (lua_next(lua, 3))
  {
    // Key is at index -2, value is at -1.
    static const int key_index = -2;
    static const int val_index = -1;
    bool key_is_number = lua_isnumber(lua, key_index);
    bool val_is_sequence = lua_issequence(lua, val_index);
    if (!key_is_number || !val_is_sequence)
    {
      lua_pop(lua, 1);
      polymec_free(faces_table_entries);
      return luaL_error(lua, "Found non-numeric entries in faces table.");
    }
    int tuple_len;
    faces_table_entries[face] = lua_tosequence(lua, val_index, &tuple_len);
    if (tuple_len != 3)
    {
      lua_pop(lua, 1);
      polymec_free(faces_table_entries);
      return luaL_error(lua, "Tuple at index %d of faces table has %d values (should be 3).", key_index, tuple_len);
    }
    ++face;
    lua_pop(lua, 1);
  }
  ASSERT(face == num_faces);

  // Check the faces data.
  for (int f = 0; f < num_faces; ++f)
  {
    real_t* face_tuple = faces_table_entries[f];
    if (face_tuple[0] < 0.0)
      return luaL_error(lua, "Invalid first cell for face %d: %d (must be non-negative).", f, (int)face_tuple[0]);
    if ((face_tuple[1] < 0.0) && (face_tuple[1] != -1.0))
      return luaL_error(lua, "Invalid second cell for face %d: %d (must be non-negative or -1).", f, (int)face_tuple[1]);
    if (face_tuple[2] < 0.0)
      return luaL_error(lua, "Invalid area for face %d: %g.", f, face_tuple[2]);
  }

  // Face centers?
  int num_face_centers = 0;
  point_t* face_centers = NULL;
  if (num_args == 4)
  {
    face_centers = lua_topointlist(lua, 1, &num_face_centers);
    if (num_face_centers != num_faces)
      return luaL_error(lua, "Number of face centers (%d) does not match number of faces (%d).", num_face_centers, num_faces);
  }

  // Shuffle the faces data into canonical form.
  int* faces = polymec_malloc(2 * sizeof(int) * num_faces);
  real_t* face_areas = polymec_malloc(sizeof(real_t) * num_faces);
  for (int f = 0; f < num_faces; ++f)
  {
    faces[2*f] = faces_table_entries[f][0];
    faces[2*f+1] = faces_table_entries[f][1];
    face_areas[f] = faces_table_entries[f][2];
    polymec_free(faces_table_entries[f]);
  }
  polymec_free(faces_table_entries);

  // Create the mesh.
  mesh_t* mesh = create_pebi_mesh(MPI_COMM_WORLD, cell_centers, cell_volumes, num_cells, 
                                  faces, face_areas, face_centers, num_faces);
  polymec_free(face_areas);
  polymec_free(faces);
  polymec_free(cell_centers);

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

int mesh_factory_dual(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 4 ) || (num_args != 5) || 
      !lua_ismesh(lua, 1) || !lua_isstringlist(lua, 2) || 
      !lua_isstringlist(lua, 3) || !lua_isstringlist(lua, 4) || 
      ((num_args == 5) && !lua_isstringlist(lua, 5)))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.dual(original_mesh, external_model_face_tags, model_edge_tags, model_vertex_tags) OR"
                      "mesh = mesh_factory.dual(original_mesh, external_model_face_tags, internal_model_face_tags, model_edge_tags, model_vertex_tags) OR");
  }

  mesh_t* orig_mesh = lua_tomesh(lua, 1);
  int num_external_model_face_tags, num_internal_model_face_tags = 0,
      num_model_edge_tags, num_model_vertex_tags;
  char** external_model_face_tags = lua_tostringlist(lua, 2, &num_external_model_face_tags);
  char** internal_model_face_tags = NULL;
  char** model_edge_tags; 
  char** model_vertex_tags;
  if (num_args == 4)
  {
    model_edge_tags = lua_tostringlist(lua, 3, &num_model_edge_tags);
    model_vertex_tags = lua_tostringlist(lua, 4, &num_model_vertex_tags);
  }
  else
  {
    internal_model_face_tags = lua_tostringlist(lua, 3, &num_external_model_face_tags);
    model_edge_tags = lua_tostringlist(lua, 4, &num_model_edge_tags);
    model_vertex_tags = lua_tostringlist(lua, 5, &num_model_vertex_tags);
  }

  // Make sure the mesh contains the given tags.
  for (int i = 0; i < num_external_model_face_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->face_tags, external_model_face_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain face tag '%s'.", external_model_face_tags[i]);
  }
  for (int i = 0; i < num_internal_model_face_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->face_tags, internal_model_face_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain face tag '%s'.", internal_model_face_tags[i]);
  }
  for (int i = 0; i < num_model_edge_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->edge_tags, model_edge_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain edge tag '%s'.", model_edge_tags[i]);
  }
  for (int i = 0; i < num_model_vertex_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->node_tags, model_vertex_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain node tag '%s'.", model_vertex_tags[i]);
  }

  // For now, we only support duals of tet meshes.
  if (!mesh_has_feature(orig_mesh, MESH_IS_TETRAHEDRAL))
    return luaL_error(lua, "mesh_factory.dual: A dual mesh can only be created from a tetrahedral mesh.");

  mesh_t* mesh = create_dual_mesh(MPI_COMM_WORLD, orig_mesh, 
                                  external_model_face_tags, num_external_model_face_tags,
                                  internal_model_face_tags, num_internal_model_face_tags,
                                  model_edge_tags, num_model_edge_tags,
                                  model_vertex_tags, num_model_vertex_tags);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

#endif

void interpreter_register_polyglot_functions(interpreter_t* interp)
{
  if (!interpreter_has_global_table(interp, "mesh_factory"))
    interpreter_register_global_table(interp, "mesh_factory", NULL);
  interpreter_register_global_method(interp, "mesh_factory", "tetgen", mesh_factory_tetgen, NULL);
//  interpreter_register_global_method(interp, "mesh_factory", "pebi", mesh_factory_pebi, NULL);
//  interpreter_register_global_method(interp, "mesh_factory", "dual", mesh_factory_dual, NULL);
  interpreter_register_function(interp, "read_exodus_mesh", lua_read_exodus_mesh, NULL);
}

