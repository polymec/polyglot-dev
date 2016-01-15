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
#include "polyglot/exodus_file.h"

void test_exodus_file_query(void** state)
{
  size_t real_size;
  float version;
  int num_mpi_processes;
  real_array_t* times = real_array_new();
  bool result = exodus_file_query("test.exo", &real_size, &version,
                                  &num_mpi_processes, times);
  assert_false(result); // has 2D elements!

  result = exodus_file_query("test-nfaced.exo", &real_size, &version,
                             &num_mpi_processes, times);
  assert_true(result);
  assert_true(real_size == sizeof(float));
  assert_true(version >= 6.09);
  assert_true(times->size == 0);
  real_array_free(times);
}

void test_write_exodus_file(void** state)
{
  fe_mesh_t* mesh = fe_mesh_new(MPI_COMM_WORLD, 22);
  int elem1_node_indices[] = {0, 1, 2, 3, 4, 5, 6, 7};
  fe_block_t* block = fe_block_new(1, FE_HEXAHEDRON, 8, elem1_node_indices);
  fe_mesh_add_block(mesh, "block_1", block);
  int elem2_node_indices[] = {8, 9, 10, 11};
  block = fe_block_new(1, FE_TETRAHEDRON, 4, elem2_node_indices);
  fe_mesh_add_block(mesh, "block_2", block);
  int elem3_node_indices[] = {12, 13, 14, 15, 16, 17};
  block = fe_block_new(1, FE_WEDGE, 6, elem3_node_indices);
  fe_mesh_add_block(mesh, "block_3", block);
  int elem4_node_indices[] = {18, 19, 20, 21};
  block = fe_block_new(1, FE_TETRAHEDRON, 4, elem4_node_indices);
  fe_mesh_add_block(mesh, "block_4", block);
  
  point_t* X = fe_mesh_node_positions(mesh);
  X[ 0].x =  0.0; X[ 0].y =  0.0; X[ 0].z =  0.0;
  X[ 1].x = 10.0; X[ 1].y =  0.0; X[ 1].z =  0.0;
  X[ 2].x = 10.0; X[ 2].y =  0.0; X[ 2].z = -10.0;
  X[ 3].x =  1.0; X[ 3].y =  0.0; X[ 3].z = -10.0;
  X[ 4].x =  0.0; X[ 4].y = 10.0; X[ 4].z =  0.0;
  X[ 5].x = 10.0; X[ 5].y = 10.0; X[ 5].z =  0.0;
  X[ 6].x = 10.0; X[ 6].y = 10.0; X[ 6].z = -10.0;
  X[ 7].x =  1.0; X[ 7].y = 10.0; X[ 7].z = -10.0;

  X[ 8].x =  0.0; X[ 8].y =  0.0; X[ 8].z =  0.0;
  X[ 9].x =  1.0; X[ 9].y =  0.0; X[ 9].z =  5.0;
  X[10].x = 10.0; X[10].y =  0.0; X[10].z =  2.0;
  X[11].x =  7.0; X[11].y =  5.0; X[11].z =  3.0;

  X[12].x =  3.0; X[12].y =  0.0; X[12].z =  6.0;
  X[13].x =  6.0; X[13].y =  0.0; X[13].z =  0.0;
  X[14].x =  0.0; X[14].y =  0.0; X[14].z =  0.0;
  X[15].x =  3.0; X[15].y =  2.0; X[15].z =  6.0;
  X[16].x =  6.0; X[16].y =  2.0; X[16].z =  2.5;
  X[17].x =  2.0; X[17].y =  2.0; X[17].z =  0.0;

  X[18].x =  2.7; X[18].y =  1.7; X[18].z =  2.7;
  X[19].x =  6.0; X[19].y =  1.7; X[19].z =  3.3;
  X[20].x =  5.7; X[20].y =  1.7; X[20].z =  1.7;
  X[21].x =  3.7; X[21].y =  0.0; X[21].z =  2.3;

  // Node sets.
  int* ns1 = fe_mesh_create_node_set(mesh, "nset_1", 5);
  ns1[0] = 1; ns1[1] = 2; ns1[2] = 3; ns1[3] = 4; ns1[4] = 5;
  int* ns2 = fe_mesh_create_node_set(mesh, "nset_2", 3);
  ns2[0] = 11; ns2[1] = 12; ns2[2] = 13; 

  // Side set.
  int* ss1 = fe_mesh_create_side_set(mesh, "sset_1", 1);
  ss1[0] = 0; ss1[1] = 5;

  exodus_file_t* file = exodus_file_new(MPI_COMM_WORLD, "test-3d.exo");
  assert_true(file != NULL);
  exodus_file_set_title(file, "This is a test");
  exodus_file_write_mesh(file, mesh);
  exodus_file_close(file);

  fe_mesh_free(mesh);
}

void test_read_exodus_file(void** state)
{
  exodus_file_t* file = exodus_file_open(MPI_COMM_WORLD, "test-3d.exo");
  assert_true(file != NULL);
  assert_true(strcmp(exodus_file_title(file), "This is a test") == 0);
  fe_mesh_t* mesh = exodus_file_read_mesh(file);
  assert_int_equal(22, fe_mesh_num_nodes(mesh));
  assert_int_equal(4, fe_mesh_num_blocks(mesh));
  assert_int_equal(4, fe_mesh_num_elements(mesh));
  assert_int_equal(0, fe_mesh_num_element_sets(mesh));
  assert_int_equal(0, fe_mesh_num_face_sets(mesh));
  assert_int_equal(0, fe_mesh_num_edge_sets(mesh));
  assert_int_equal(2, fe_mesh_num_node_sets(mesh));
  assert_int_equal(1, fe_mesh_num_side_sets(mesh));

  int elem_nodes[10];

  assert_int_equal(8, fe_mesh_num_element_nodes(mesh, 0));
  fe_mesh_get_element_nodes(mesh, 0, elem_nodes);
  assert_int_equal(0, elem_nodes[0]);
  assert_int_equal(1, elem_nodes[1]);
  assert_int_equal(2, elem_nodes[2]);
  assert_int_equal(3, elem_nodes[3]);
  assert_int_equal(4, elem_nodes[4]);
  assert_int_equal(5, elem_nodes[5]);
  assert_int_equal(6, elem_nodes[6]);
  assert_int_equal(7, elem_nodes[7]);

  assert_int_equal(4, fe_mesh_num_element_nodes(mesh, 1));
  fe_mesh_get_element_nodes(mesh, 1, elem_nodes);
  assert_int_equal(8, elem_nodes[0]);
  assert_int_equal(9, elem_nodes[1]);
  assert_int_equal(10, elem_nodes[2]);
  assert_int_equal(11, elem_nodes[3]);

  assert_int_equal(6, fe_mesh_num_element_nodes(mesh, 2));
  fe_mesh_get_element_nodes(mesh, 2, elem_nodes);
  assert_int_equal(12, elem_nodes[0]);
  assert_int_equal(13, elem_nodes[1]);
  assert_int_equal(14, elem_nodes[2]);
  assert_int_equal(15, elem_nodes[3]);
  assert_int_equal(16, elem_nodes[4]);
  assert_int_equal(17, elem_nodes[5]);

  assert_int_equal(4, fe_mesh_num_element_nodes(mesh, 3));
  fe_mesh_get_element_nodes(mesh, 3, elem_nodes);
  assert_int_equal(18, elem_nodes[0]);
  assert_int_equal(19, elem_nodes[1]);
  assert_int_equal(20, elem_nodes[2]);
  assert_int_equal(21, elem_nodes[3]);

  int pos = 0;
  char* block_name;
  fe_block_t* block;
  assert_true(fe_mesh_next_block(mesh, &pos, &block_name, &block));
  assert_true(strcmp(block_name, "block_1") == 0);
  assert_true(fe_block_element_type(block) == FE_HEXAHEDRON);
  assert_true(fe_mesh_next_block(mesh, &pos, &block_name, &block));
  assert_true(strcmp(block_name, "block_2") == 0);
  assert_true(fe_block_element_type(block) == FE_TETRAHEDRON);
  assert_true(fe_mesh_next_block(mesh, &pos, &block_name, &block));
  assert_true(strcmp(block_name, "block_3") == 0);
  assert_true(fe_block_element_type(block) == FE_WEDGE);
  assert_true(fe_mesh_next_block(mesh, &pos, &block_name, &block));
  assert_true(strcmp(block_name, "block_4") == 0);
  assert_true(fe_block_element_type(block) == FE_TETRAHEDRON);

  // Node sets.
  int* set, set_size;
  char* set_name;
  pos = 0;
  int num_node_sets = 0;
  while (fe_mesh_next_node_set(mesh, &pos, &set_name, &set, &set_size))
  {
    ++num_node_sets;
    assert_true((strcmp(set_name, "nset_1") == 0) || 
                (strcmp(set_name, "nset_2") == 0));
    if (strcmp(set_name, "nset_1") == 0)
    {
      assert_int_equal(5, set_size);
      assert_int_equal(1, set[0]);
      assert_int_equal(2, set[1]);
      assert_int_equal(3, set[2]);
      assert_int_equal(4, set[3]);
      assert_int_equal(5, set[4]);
    }
    else
    {
      assert_int_equal(3, set_size);
      assert_int_equal(11, set[0]);
      assert_int_equal(12, set[1]);
      assert_int_equal(13, set[2]);
    }
  }
  assert_int_equal(2, num_node_sets);

  fe_mesh_free(mesh);
  exodus_file_close(file);
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

  point_t* X = fe_mesh_node_positions(mesh);
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
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_exodus_file_query),
    cmocka_unit_test(test_write_exodus_file),
    cmocka_unit_test(test_read_exodus_file),
    cmocka_unit_test(test_read_poly_exodus_file),
    cmocka_unit_test(test_write_poly_exodus_file)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
