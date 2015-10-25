
# We use the Lemon parser generator to generate parsers for certain input 
# formats.
add_custom_target(tough2_file_parser.c ALL
                  cp ${CMAKE_CURRENT_SOURCE_DIR}/tough2_file_parser.y ${CMAKE_CURRENT_BINARY_DIR}
                  COMMAND lemon -T${CMAKE_CURRENT_SOURCE_DIR}/lempar.c tough2_file_parser.y 
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Library.
add_library(polyglot polyglot.c tough2_file.c cf_file.c)
add_dependencies(polyglot tough2_file_parser.c)

set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_TP_LIBRARIES};${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_TP_LIBRARIES};${POLYMEC_LIBRARIES})

#add_subdirectory(tests)

# Installation targets.
install(TARGETS polyglot DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polyglot)

