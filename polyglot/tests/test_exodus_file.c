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

void test_exodus_file_query(void** state)
{
  size_t real_size;
  float version;
  int num_mpi_processes;
  real_array_t* times = real_array_new();
  bool result = exodus_file_query("test-nfaced.exo", &real_size, &version,
                                  &num_mpi_processes, times);
  assert_true(result);
  assert_true(real_size == sizeof(float));
  assert_true(version >= 6.09);
  assert_true(times->size == 0);
  real_array_free(times);
}

void test_read_exodus_file(void** state)
{
  exodus_file_t* file = exodus_file_open(MPI_COMM_WORLD, "test-nfaced.exo");
  assert_true(file != NULL);
  assert_true(strcmp(exodus_file_title(file), "This is a test") == 0);
  fe_mesh_t* mesh = exodus_file_read_mesh(file);
  assert_int_equal(14, fe_mesh_num_nodes(mesh));
  assert_int_equal(1, fe_mesh_num_blocks(mesh));
  assert_int_equal(3, fe_mesh_num_elements(mesh));
  assert_int_equal(0, fe_mesh_num_element_sets(mesh));
  assert_int_equal(0, fe_mesh_num_face_sets(mesh));
  assert_int_equal(0, fe_mesh_num_edge_sets(mesh));
  assert_int_equal(0, fe_mesh_num_node_sets(mesh));
  assert_int_equal(0, fe_mesh_num_side_sets(mesh));

  int elem_faces[10];

  assert_int_equal(5, fe_mesh_num_element_faces(mesh, 0));
  fe_mesh_get_element_faces(mesh, 0, elem_faces);
  assert_int_equal(0, elem_faces[0]);
  assert_int_equal(1, elem_faces[1]);
  assert_int_equal(2, elem_faces[2]);
  assert_int_equal(3, elem_faces[3]);
  assert_int_equal(4, elem_faces[4]);

  assert_int_equal(5, fe_mesh_num_element_faces(mesh, 1));
  fe_mesh_get_element_faces(mesh, 1, elem_faces);
  assert_int_equal(3, elem_faces[0]);
  assert_int_equal(5, elem_faces[1]);
  assert_int_equal(6, elem_faces[2]);
  assert_int_equal(7, elem_faces[3]);
  assert_int_equal(8, elem_faces[4]);

  assert_int_equal(7, fe_mesh_num_element_faces(mesh, 2));
  fe_mesh_get_element_faces(mesh, 2, elem_faces);
  assert_int_equal(7, elem_faces[0]);
  assert_int_equal(9, elem_faces[1]);
  assert_int_equal(10, elem_faces[2]);
  assert_int_equal(11, elem_faces[3]);
  assert_int_equal(12, elem_faces[4]);
  assert_int_equal(13, elem_faces[5]);
  assert_int_equal(14, elem_faces[6]);

  int pos = 0;
  char* block_name;
  fe_block_t* block;
  assert_true(fe_mesh_next_block(mesh, &pos, &block_name, &block));
  assert_true(strcmp(block_name, "nfaced_1") == 0);
  assert_int_equal(5, fe_block_num_element_faces(block, 0));

  fe_mesh_free(mesh);
  exodus_file_close(file);
}

void test_write_exodus_file(void** state)
{
  exodus_file_t* file = exodus_file_new(MPI_COMM_WORLD, "test-nfaced-2.exo");
  assert_true(file != NULL);
  exodus_file_close(file);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exodus_file_query),
    unit_test(test_read_exodus_file),
    unit_test(test_write_exodus_file)
  };
  return run_tests(tests);
}
