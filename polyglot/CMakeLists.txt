
# Library.
set(POLYGLOT_SOURCES polyglot.c import_tetgen_mesh.c 
                     fe_mesh.c exodus_file.c cf_file.c 
                     interpreter_register_polyglot_functions.c)
if (HAVE_POLYAMRI)
  include(add_polyamri_library)
  add_polyamri_library(polyglot ${POLYGLOT_SOURCES})
else()
  include(add_polymec_library)
  add_polymec_library(polyglot ${POLYGLOT_SOURCES})
endif()
add_dependencies(polyglot update_version_h) # <-- needed on Mac(?!)

if (HAVE_POLYAMRI)
  set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYAMRI_LIBRARIES} PARENT_SCOPE)
  set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYAMRI_LIBRARIES})
else()
  set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYMEC_LIBRARIES} PARENT_SCOPE)
  set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYMEC_LIBRARIES})
endif()

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyglot DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/polyglot)

