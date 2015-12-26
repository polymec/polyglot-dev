# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polyglot_test exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} cmockery ${POLYGLOT_LIBRARIES})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  add_test(${exe} ${exe})
  set_tests_properties(${exe} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

# This function adds a parallel unit test executable to be built using gtest.
# The procs argument is a list of numbers of processes to be run.
# 1 test run will be generated for each processor number value.
function(add_mpi_polyglot_test exe)
  foreach (arg ${ARGN})
    if (arg MATCHES ".c")
      list(APPEND sources ${arg})
    else()
      list(APPEND procs ${arg})
    endif()
  endforeach()
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} cmockery ${POLYGLOT_LIBRARIES})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  if (POLYMEC_HAVE_MPI EQUAL 1)
    foreach (proc ${procs})
      add_test(${exe}_${proc}_proc ${POLYMEC_MPIEXEC} ${POLYMEC_MPIEXEC_NUMPROC_FLAG} ${proc} ${POLYMEC_MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${POLYMEC_MPIEXEC_POSTFLAGS})
      set_tests_properties(${exe}_${proc}_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endforeach()
  else()
    # We only add a single-process test case when MPI is not present.
    add_test(${exe}_1_proc ${CMAKE_CURRENT_BINARY_DIR}/${exe})
    set_tests_properties(${exe}_1_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endfunction()

