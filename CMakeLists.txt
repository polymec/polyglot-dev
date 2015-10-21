# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

# Project and version numbers.
project (polyglot)
set (POLYGLOT_MAJOR_VERSION 0)
set (POLYGLOT_MINOR_VERSION 6)
set (POLYGLOT_PATCH_VERSION 0)
set (POLYGLOT_VERSION "${POLYGLOT_MAJOR_VERSION}.${POLYGLOT_MINOR_VERSION}.${POLYGLOT_PATCH_VERSION}")
message("-- Building polyglot (v${POLYGLOT_VERSION})")

# Use our own CMake stuff.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Look for Polymec and use its CMake settings.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${POLYMEC_PREFIX}/share/polymec/")
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

# Do we have polyamri?
include(polyamri)
if (DEFINED POLYAMRI_VERSION)
  message("-- Found polyamri library in ${POLYMEC_PREFIX} (v${POLYAMRI_VERSION})")
  set(HAVE_POLYAMRI 1)
endif()

# Generate polyglot_version.h with the right version numbers.
add_custom_target(update_version_h ALL
                  python ${POLYMEC_PREFIX}/share/polymec/update_version_h.py polyglot ${POLYGLOT_VERSION} ${PROJECT_BINARY_DIR}/polyglot/polyglot_version.h)

# Build 3rd-party libraries.
add_subdirectory(3rdparty)

# Look for headers in the top directory.
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

# Look for libraries in these directories.
link_directories("${POLYMEC_PREFIX}/lib")

enable_testing()

add_subdirectory(polyglot)

# Generate a polyglot.cmake file that contains installation information.
set(CMAKE_INSTALL_PREFIX ${POLYMEC_PREFIX})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polyglot.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polyglot.cmake"
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polyglot.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyglot_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyglot_test.cmake DESTINATION share/polymec)

