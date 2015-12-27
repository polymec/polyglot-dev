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

void test_read_poly_exodus_file(void** state)
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

void test_write_poly_exodus_file(void** state)
{
  fe_mesh_t* mesh = fe_mesh_new(MPI_COMM_WORLD, 14);
  int num_elem_faces[3] = {5, 5, 7};
  int elem_face_indices[] = {0, 1, 2, 3, 4, 3, 5, 6, 7, 8, 7, 9, 10, 11, 12, 13, 14};
  fe_block_t* block = polyhedral_fe_block_new(3, num_elem_faces, elem_face_indices);
  fe_mesh_add_block(mesh, "nfaced_1", block);
  
  int num_face_nodes[15] = {2, 2, 3, 3, 3, 2, 2, 3, 3, 4, 4, 3, 3, 3, 3};
  int face_nodes[] = {4, 5, 7, 1, 0, 3, 5, 1, 3, 7, 7, 3, 0, 4, 0, 1, 5, 4, 4, 7, 6, 0,
                      2, 3, 6, 7, 3, 2, 6, 2, 0, 4, 7, 3, 13, 9, 11, 6, 10, 8, 12, 2, 6, 7,
                      11, 10, 10, 11, 9, 8, 8, 9, 13, 12, 12, 13, 3, 2};
  fe_mesh_set_face_nodes(mesh, 15, num_face_nodes, face_nodes);

  point_t* X = fe_mesh_node_coordinates(mesh);
  X[ 0].x = 0.0; X[ 0].y = 0.0; X[ 0].z = 0.0;
  X[ 1].x = 2.0; X[ 1].y = 0.0; X[ 1].z = 0.0;
  X[ 2].x = 0.0; X[ 2].y = 2.0; X[ 2].z = 0.0;
  X[ 3].x = 2.0; X[ 3].y = 2.0; X[ 3].z = 0.0;
  X[ 4].x = 0.0; X[ 4].y = 0.0; X[ 4].z = 2.0;
  X[ 5].x = 2.0; X[ 5].y = 0.0; X[ 5].z = 2.0;
  X[ 6].x = 0.0; X[ 6].y = 2.0; X[ 6].z = 2.0;
  X[ 7].x = 2.0; X[ 7].y = 2.0; X[ 7].z = 2.0;
  X[ 8].x = 0.0; X[ 8].y = 3.5; X[ 8].z = 1.0;
  X[ 9].x = 2.0; X[ 9].y = 3.5; X[ 9].z = 1.0;
  X[10].x = 0.0; X[10].y = 3.0; X[10].z = 1.0;
  X[11].x = 2.0; X[11].y = 3.0; X[11].z = 1.0;
  X[12].x = 0.0; X[12].y = 3.0; X[12].z = 0.5;
  X[13].x = 2.0; X[13].y = 3.0; X[13].z = 0.5;

  exodus_file_t* file = exodus_file_new(MPI_COMM_WORLD, "test-nfaced-2.exo");
  assert_true(file != NULL);
  exodus_file_set_title(file, "This is a test");
  exodus_file_write_mesh(file, mesh);
  exodus_file_close(file);

  fe_mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exodus_file_query),
    unit_test(test_read_poly_exodus_file),
    unit_test(test_write_poly_exodus_file)
  };
  return run_tests(tests);
}
