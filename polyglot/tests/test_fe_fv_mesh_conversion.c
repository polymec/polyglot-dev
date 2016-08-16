// Copyright (c) 2015-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "geometry/create_uniform_mesh.h"
#include "polyglot/exodus_file.h"

static void test_mesh_from_fe_mesh(void** state)
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

static void test_fe_mesh_from_mesh(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};
  mesh_t* fv_mesh = create_uniform_mesh(MPI_COMM_SELF, 10, 10, 10, &bbox);
  fe_mesh_t* fe_mesh = fe_mesh_from_mesh(fv_mesh, NULL);
  mesh_free(fv_mesh);
  assert_int_equal(1000, fe_mesh_num_elements(fe_mesh));
  assert_int_equal(1, fe_mesh_num_blocks(fe_mesh));
  assert_int_equal(1331, fe_mesh_num_nodes(fe_mesh));
  assert_int_equal(3300, fe_mesh_num_faces(fe_mesh));
  fe_mesh_free(fe_mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_mesh_from_fe_mesh),
    cmocka_unit_test(test_fe_mesh_from_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
