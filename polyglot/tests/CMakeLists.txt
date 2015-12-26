include(add_polyglot_executable)
include(add_polyglot_test)

# Tetgen mesh import.
add_mpi_polyglot_test(test_import_tetgen_mesh test_import_tetgen_mesh.c 1 2 4)

# Exodus test data.
add_polyglot_executable(testwt-nfaced ${PROJECT_SOURCE_DIR}/3rdparty/exodus/exodus/cbind/test/testwt-nfaced.c)
add_test(NAME generate_exodus_data
         COMMAND testwt-nfaced)

# Exodus tests.
add_polyglot_test(test_exodus_file test_exodus_file.c)
set_tests_properties(test_exodus_file PROPERTIES DEPENDS generate_exodus_data)
