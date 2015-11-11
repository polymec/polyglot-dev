include_directories(${PROJECT_BINARY_DIR}/include)

# Build the netcdf library.
set(NETCDF_LIBRARY "${PROJECT_BINARY_DIR}/lib/libnetcdf.a")
if (NOT EXISTS ${NETCDF_LIBRARY})
  message("Preparing netcdf...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/netcdf ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring netcdf library...")
  execute_process(COMMAND "${CMAKE_COMMAND}" -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DUSE_HDF5=ON -DHDF5_LIB=${POLYMEC_HDF5_LIBRARY} -DHDF5_HL_LIB=${POLYMEC_HDF5_HL_LIBRARY} -DHDF5_INCLUDE_DIR=${POLYMEC_HDF5_INCLUDE_DIR} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DENABLE_DAP=OFF -DBUILD_SHARED_LIBS=OFF -DFIND_SHARED_LIBS=OFF -DENABLE_DYNAMIC_LOADING=ON -DENABLE_TESTS=OFF -DBUILD_UTILITIES=OFF -DBUILD_EXAMPLES=OFF -DNC_HAVE_PARALLEL_HDF5=${POLYMEC_HAVE_MPI} -DHDF5_IS_PARALLEL=${POLYMEC_HAVE_MPI}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  OUTPUT_VARIABLE netcdf_config_log ERROR_VARIABLE netcdf_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config.log ${netcdf_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config_errors.log ${netcdf_config_err})
    message(FATAL_ERROR "Error configuring netcdf library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config*.log for details.")
  endif()

  message("Building and installing netcdf library...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  OUTPUT_VARIABLE netcdf_build_log ERROR_VARIABLE netcdf_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build.log ${netcdf_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build_errors.log ${netcdf_build_err})
    message(FATAL_ERROR "Error building netcdf library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build*.log for details.")
  endif()
endif()

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  if (${NETCDF_LIBRARY} STREQUAL "${PROJECT_BINARY_DIR}/lib/libnetcdf.a")
    install(FILES ${PROJECT_BINARY_DIR}/lib/libnetcdf.a DESTINATION lib)
  endif()
endif()
set(POLYGLOT_TP_LIBS ${NETCDF_LIBRARY};${POLYGLOT_TP_LIBS})

# Build the Exodus library.
if (NOT EXISTS "${PROJECT_BINARY_DIR}/lib/libexodus.a")
  message("Preparing exodus...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/exodus ${CMAKE_CURRENT_BINARY_DIR}/exodus
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Adjusting exodus CMakeLists.txt...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt exodus_CMakeLists)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt.old ${exodus_CMakeLists})
  # Remove the Fortran bindings from the build.
  string(REPLACE "forbind" "" exodus_CMakeLists ${exodus_CMakeLists})
  # Remove the Python bindings stuff.
  string(REPLACE "if (PYTHON)" "if (0)" exodus_CMakeLists ${exodus_CMakeLists})
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt ${exodus_CMakeLists})
 
  message("Configuring exodus library...")
  execute_process(COMMAND env NETCDF_DIR=${PROJECT_BINARY_DIR} ${CMAKE_COMMAND} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED=OFF 
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus
                  OUTPUT_VARIABLE exodus_config_log ERROR_VARIABLE exodus_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_config.log ${exodus_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_config_errors.log ${exodus_config_err})
    message(FATAL_ERROR "Error configuring exodus library. See ${CMAKE_CURRENT_BINARY_DIR}/exodus_config*.log for details.")
  endif()

  message("Building and installing exodus library...")
  execute_process(COMMAND env NETCDF_DIR=${PROJECT_BINARY_DIR} make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus
                  OUTPUT_VARIABLE exodus_build_log ERROR_VARIABLE exodus_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_build.log ${exodus_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_build_errors.log ${exodus_build_err})
    message(FATAL_ERROR "Error building exodus library. See ${CMAKE_CURRENT_BINARY_DIR}/exodus_build*.log for details.")
  endif()
endif()
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(FILES ${PROJECT_BINARY_DIR}/lib/libexodus.a DESTINATION lib)
endif()
set(POLYGLOT_TP_LIBS exodus;${POLYGLOT_TP_LIBS})

# Build the Lemon parser generator, which we use to generate parsers for 
# certain troublesome input formats. See the documentation at 
# http://www.hwaci.com/sw/lemon/lemon.html
# and the tutorial at 
# http://souptonuts.sourceforge.net/readme_lemon_tutorial.html
# for details on usage.
add_executable(lemon lemon.c)

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYGLOT_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYGLOT_LIBRARIES ${POLYGLOT_TP_LIBS};${POLYGLOT_LIBRARIES} PARENT_SCOPE)
set(POLYGLOT_INCDIRS ${POLYGLOT_INCDIRS};${POLYGLOT_TP_INCDIRS} PARENT_SCOPE)

