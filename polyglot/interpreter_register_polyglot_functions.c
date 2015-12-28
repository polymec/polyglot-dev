// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/interpreter.h"
#include "polyglot/interpreter_register_polyglot_functions.h"
#include "polyglot/exodus_file.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This file contains definitions for functions supported by Polyglot.

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

void interpreter_register_polyglot_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "read_exodus_mesh", lua_read_exodus_mesh, NULL);
}

