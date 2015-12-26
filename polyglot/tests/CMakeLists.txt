include(add_polyglot_executable)
include(add_polyglot_test)

# Tetgen mesh import.
add_mpi_polyglot_test(test_import_tetgen_mesh test_import_tetgen_mesh.c 1 2 4)

# Exodus test data.
set(exo_test_exes testwt-nfaced)
foreach(exo_test_exe ${exo_test_exes})
  add_polyglot_executable(${exo_test_exe} ${PROJECT_SOURCE_DIR}/3rdparty/exodus/exodus/cbind/test/${exo_test_exe}.c)
  set(exo_test_script "${exo_test_script}${CMAKE_CURRENT_BINARY_DIR}/${exo_test_exe}\n")
endforeach()
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/generate_exodus_data.sh ${exo_test_script})
add_test(NAME generate_exodus_data
         COMMAND sh generate_exodus_data.sh)

# Exodus tests.
add_polyglot_test(test_exodus_file test_exodus_file.c)
set_tests_properties(test_exodus_file PROPERTIES DEPENDS generate_exodus_data)
