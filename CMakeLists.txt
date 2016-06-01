# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Project and version numbers.
project (polyglot)
set (POLYGLOT_MAJOR_VERSION 1)
set (POLYGLOT_MINOR_VERSION 0)
set (POLYGLOT_PATCH_VERSION 0)
set (POLYGLOT_VERSION "${POLYGLOT_MAJOR_VERSION}.${POLYGLOT_MINOR_VERSION}.${POLYGLOT_PATCH_VERSION}")
message("-- Building polyglot (v${POLYGLOT_VERSION})")

# Use our own CMake stuff.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Look for headers in the top directory.
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/include)

# Look for Polymec and use its CMake settings.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${POLYMEC_PREFIX}/share/polymec/")
include(polymec)
if (NOT DEFINED POLYMEC_VERSION)
  message(FATAL_ERROR "Could not find polymec library. Please use the polymec=/path/to/polymec option.")
endif()
message("-- Found polymec library in ${POLYMEC_PREFIX} (v${POLYMEC_VERSION})")

# Build everything the way Polymec does.
set(CMAKE_C_COMPILER ${POLYMEC_C_COMPILER})
set(CMAKE_C_FLAGS ${POLYMEC_C_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${POLYMEC_EXE_LINKER_FLAGS})
set(BUILD_SHARED_LIBS ${POLYMEC_BUILD_SHARED_LIBS})
set(CMAKE_SHARED_LINKER_FLAGS ${POLYMEC_SHARED_LINKER_FLAGS})
set(CMAKE_BUILD_TYPE ${POLYMEC_BUILD_TYPE})
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${POLYMEC_PREFIX}/lib")

# Do we have polyamri?
if (EXISTS ${POLYMEC_PREFIX}/share/polymec/polyamri.cmake)
  include(polyamri)
  if (DEFINED POLYAMRI_VERSION)
    message("-- Found polyamri library in ${POLYMEC_PREFIX} (v${POLYAMRI_VERSION})")
    set(HAVE_POLYAMRI 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYGLOT_HAVE_POLYAMRI")
  endif()
endif()

# Generate polyglot_version.h with the right version numbers.
add_custom_target(update_version_h ALL
                  python ${POLYMEC_PREFIX}/share/polymec/update_version_h.py polyglot ${POLYGLOT_VERSION} ${PROJECT_BINARY_DIR}/polyglot/polyglot_version.h)

# Generate polyamri.h at the top level.
add_custom_target(generate_top_header ALL
                  python ${PROJECT_SOURCE_DIR}/tools/generate_top_header.py ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
                  
# Build 3rd-party libraries.
link_directories("${PROJECT_BINARY_DIR}/lib")
add_subdirectory(3rdparty)

enable_testing()

# Library proper.
add_subdirectory(polyglot)

# We bundle a mesh generator/editor/importer/exporter.
add_subdirectory(polymesher)

# Generate a polyglot.cmake file that contains installation information.
set(CMAKE_INSTALL_PREFIX ${POLYMEC_PREFIX})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polyglot.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polyglot.cmake"
)
install(FILES ${PROJECT_BINARY_DIR}/include/polyglot.h DESTINATION include/polymec)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polyglot.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyglot_executable.cmake DESTINATION share/polymec)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polyglot_test.cmake DESTINATION share/polymec)

