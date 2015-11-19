include(add_polyglot_test)

add_mpi_polyglot_test(test_import_tetgen_mesh test_import_tetgen_mesh.c 1 2 4)
