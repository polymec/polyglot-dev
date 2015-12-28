# Library.
add_library(polyglot polyglot.c fe_mesh.c exodus_file.c cf_file.c 
            import_tetgen_mesh.c tough2_file.c 
            interpreter_register_polyglot_functions.c)

set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_LIBRARIES};${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyglot DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polyglot)

