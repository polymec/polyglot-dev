// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "polyglot/exodus_file.h"

void test_mesh_from_fe_mesh(void** state)
{
  exodus_file_t* file = exodus_file_open(MPI_COMM_WORLD, "test-3d.exo");
  fe_mesh_t* fe_mesh = exodus_file_read_mesh(file);
  exodus_file_close(file);
  mesh_t* fv_mesh = mesh_from_fe_mesh(fe_mesh);
  fe_mesh_free(fe_mesh);
  assert_int_equal(4, fv_mesh->num_cells);
  assert_int_equal(19, fv_mesh->num_faces);
  assert_int_equal(22, fv_mesh->num_nodes);
  assert_true(mesh_has_tag(fv_mesh->node_tags, "nset_1"));
  assert_true(mesh_has_tag(fv_mesh->node_tags, "nset_2"));
  mesh_free(fv_mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_mesh_from_fe_mesh)
  };
  return run_tests(tests);
}
