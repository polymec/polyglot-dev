include(add_polyglot_executable)

# Executable.
add_polyglot_executable(polymesher polymesher.c write_gnuplot_points.c)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/polymesher DESTINATION bin)
endif()
