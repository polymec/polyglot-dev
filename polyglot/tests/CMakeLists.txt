include(add_polyglot_executable)
include(add_polyglot_test)

# Tetgen mesh import.
add_mpi_polyglot_test(test_import_tetgen_mesh test_import_tetgen_mesh.c 1 2 4)

# Exodus test data.
set(exo_test_exes testwt-nfaced;testwt)
foreach(exo_test_exe ${exo_test_exes})
  # Use these Exodus test programs to generate some Exodus data. We needed to 
  # modify their source to force them to create NetCDF4/HDF5 data by changing the 
  # file creation mode from NC_CLOBBER to NC_CLOBBER | NC_NETCDF4.
  add_polyglot_executable(${exo_test_exe} ${CMAKE_CURRENT_SOURCE_DIR}/${exo_test_exe}.c)

  # Here's the test script.
  set(exo_test_script "${exo_test_script}${CMAKE_CURRENT_BINARY_DIR}/${exo_test_exe}\n")
endforeach()
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/generate_exodus_data.sh ${exo_test_script})
add_test(NAME generate_exodus_data
         COMMAND sh generate_exodus_data.sh)

# Exodus tests.
add_polyglot_test(test_exodus_file test_exodus_file.c)
set_tests_properties(test_exodus_file PROPERTIES DEPENDS generate_exodus_data)

# CF format tests.
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cf_test_data.nc)
  set(cf_test_data_url https://www.unidata.ucar.edu/software/netcdf/examples/sresa1b_ncar_ccsm3-example.nc)
  message("Downloading NetCDF CF test data: ${cf_test_data_url} -> ${CMAKE_CURRENT_SOURCE_DIR}/cf_test_data.nc")
  file(DOWNLOAD ${cf_test_data_url} ${CMAKE_CURRENT_SOURCE_DIR}/cf_test_data.nc)
endif()
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cf_test_data.nc)
  add_polyglot_test(test_cf_file test_cf_file.c)
endif()

# FE <--> FV mesh conversion.
add_polyglot_test(test_fe_fv_mesh_conversion test_fe_fv_mesh_conversion.c)
set_tests_properties(test_fe_fv_mesh_conversion PROPERTIES DEPENDS test_exodus_file)
