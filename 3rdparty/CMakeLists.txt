include_directories(${PROJECT_BINARY_DIR}/include)

if (BUILD_SHARED_LIBS)
  set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(LIB_SUFFIX .a)
endif()

# Build the netcdf library.
set(NETCDF_LIBRARY "${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX}")
if (NOT EXISTS ${NETCDF_LIBRARY})
  message(STATUS "Preparing netcdf...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/netcdf ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Configuring netcdf library...")
  # The amazing amount of crap required in this CMake command suggests that NetCDF/HDF5/CMake 
  # coverage is not particularly well tested. :-/
  string(REPLACE "-pedantic-errors" "" NETCDF_C_FLAGS ${CMAKE_C_FLAGS})
  string(REPLACE "-Werror-implicit-function-declaration" "" NETCDF_C_FLAGS ${NETCDF_C_FLAGS})
  execute_process(COMMAND "${CMAKE_COMMAND}" -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_C_FLAGS=${NETCDF_C_FLAGS}
                  -DUSE_HDF5=ON -DHDF5_C_LIBRARY=${POLYMEC_HDF5_LIBRARY} -DHDF5_HL_LIBRARY=${POLYMEC_HDF5_HL_LIBRARY} 
                  -DHDF5_C_LIBRARIES=${POLYMEC_HDF5_LIBRARY} -DHDF5_HL_LIBRARIES=${POLYMEC_HDF5_HL_LIBRARY} 
                  -DHDF5_INCLUDE_DIR=${POLYMEC_HDF5_INCLUDE_DIR} -DENABLE_PARALLEL=${POLYMEC_HAVE_MPI} -DHDF5_IS_PARALLEL_MPIO=${POLYMEC_HAVE_MPI} 
                  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DENABLE_DAP=OFF -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} 
                  -DFIND_SHARED_LIBS=${BUILD_SHARED_LIBS} -DENABLE_DYNAMIC_LOADING=${BUILD_SHARED_LIBS} -DENABLE_TESTS=OFF -DBUILD_UTILITIES=ON -DBUILD_EXAMPLES=OFF 
                  -DNC_HAVE_PARALLEL_HDF5=${POLYMEC_HAVE_MPI} -DHDF5_IS_PARALLEL=${POLYMEC_HAVE_MPI} -DCMAKE_FIND_LIBRARY_SUFFIXES=${LIB_SUFFIX}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  OUTPUT_VARIABLE netcdf_config_log ERROR_VARIABLE netcdf_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config.log ${netcdf_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config_errors.log ${netcdf_config_err})
    message(FATAL_ERROR "Error configuring netcdf library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_config*.log for details.")
  endif()

  message(STATUS "Building and installing netcdf library...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf
                  OUTPUT_VARIABLE netcdf_build_log ERROR_VARIABLE netcdf_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build.log ${netcdf_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build_errors.log ${netcdf_build_err})
    message(FATAL_ERROR "Error building netcdf library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_build*.log for details.")
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libnetcdf${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX})
  endif()

endif()

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  if (${NETCDF_LIBRARY} STREQUAL "${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX}")
    file(GLOB nc_libs "${PROJECT_BINARY_DIR}/lib/libnetcdf*${LIB_SUFFIX}*")
    install(FILES ${nc_libs} DESTINATION lib)
    foreach(nc_util nc-config;nccopy;ncdump;ncgen;ncgen3)
      list(APPEND nc_utils ${PROJECT_BINARY_DIR}/bin/${nc_util})
    endforeach()
    install(FILES ${nc_utils} DESTINATION bin
            PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;
                        GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE)
  endif()
endif()
set(POLYGLOT_TP_LIBS netcdf;${POLYGLOT_TP_LIBS})

# Build the Exodus library.
if (NOT EXISTS "${PROJECT_BINARY_DIR}/lib/libexodus${LIB_SUFFIX}")
  message(STATUS "Preparing exodus...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/exodus ${CMAKE_CURRENT_BINARY_DIR}/exodus
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message(STATUS "Adjusting exodus CMakeLists.txt...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt exodus_CMakeLists)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt.old ${exodus_CMakeLists})
  # Remove the Fortran bindings from the build.
  string(REPLACE "add_subdirectory(forbind)" "#add_subdirectory(forbind)" exodus_CMakeLists ${exodus_CMakeLists})
  # Remove the Python bindings stuff.
  string(REPLACE "if (PYTHON)" "if (0)" exodus_CMakeLists ${exodus_CMakeLists})
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/CMakeLists.txt ${exodus_CMakeLists})
 
  message(STATUS "Configuring exodus library...")
  string(REPLACE "-pedantic-errors" "" EXODUS_C_FLAGS ${CMAKE_C_FLAGS})
  execute_process(COMMAND ${CMAKE_COMMAND} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${EXODUS_C_FLAGS}
                  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} 
                  -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED=${BUILD_SHARED_LIBS} 
                  -DNETCDF_LIBRARY=${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX} 
                  -DNETCDF_INCLUDE_DIR=${PROJECT_BINARY_DIR}/include 
                  -DHDF5_LIBRARY=${POLYMEC_HDF5_LIBRARY} -DHDF5HL_LIBRARY=${POLYMEC_HDF5_HL_LIBRARY}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus
                  OUTPUT_VARIABLE exodus_config_log ERROR_VARIABLE exodus_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_config.log ${exodus_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_config_errors.log ${exodus_config_err})
    message(FATAL_ERROR "Error configuring exodus library. See ${CMAKE_CURRENT_BINARY_DIR}/exodus_config*.log for details.")
  endif()

  message(STATUS "Building and installing exodus library...")
  execute_process(COMMAND env NETCDF_DIR=${PROJECT_BINARY_DIR} make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus
                  OUTPUT_VARIABLE exodus_build_log ERROR_VARIABLE exodus_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_build.log ${exodus_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/exodus_build_errors.log ${exodus_build_err})
    message(FATAL_ERROR "Error building exodus library. See ${CMAKE_CURRENT_BINARY_DIR}/exodus_build*.log for details.")
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libexodus${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libexodus${LIB_SUFFIX})
  endif()

endif()
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  foreach (header exodusII.h;exodusII_cfg.h;exodusII_int.h)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/cbind/include/${header} DESTINATION include)
  endforeach()
  file(GLOB exo_libs "${PROJECT_BINARY_DIR}/lib/libexodus*${LIB_SUFFIX}*")
  install(FILES ${exo_libs} DESTINATION lib)
  if (POLYMEC_HAVE_MPI)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/cbind/include/exodusII_par.h DESTINATION include)
  endif()
endif()
set(POLYGLOT_TP_LIBS exodus;${POLYGLOT_TP_LIBS})

# Build the Exodus parallel library if we need it.
if (POLYMEC_HAVE_MPI) 
  # We need to manually copy the parallel header.
  execute_process(COMMAND cmake -E copy ${CMAKE_CURRENT_SOURCE_DIR}/exodus/exodus/cbind/include/exodusII_par.h 
                                        ${PROJECT_BINARY_DIR}/include
                  RESULT_VARIABLE stat)
  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/cbind/include/exodusII_par.h DESTINATION include)
  endif()

  # Build an exodus_par library.
  foreach (src ex_open_par.c;ex_create_par.c)
    list(APPEND exodus_par_sources ${CMAKE_CURRENT_BINARY_DIR}/exodus/exodus/cbind/src/${src})
  endforeach()
  add_library(exodus_par ${exodus_par_sources})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(exodus_par exodus;netcdf)
  endif()
  set(POLYGLOT_TP_LIBS exodus_par;${POLYGLOT_TP_LIBS})
endif()

# Add all the libraries to the build system at large.
set(POLYGLOT_LIBRARIES ${POLYGLOT_TP_LIBS};${POLYGLOT_LIBRARIES} PARENT_SCOPE)
set(POLYGLOT_INCDIRS ${POLYGLOT_INCDIRS};${POLYGLOT_TP_INCDIRS} PARENT_SCOPE)

