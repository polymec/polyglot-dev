
# Library.
add_library(polyglot polyglot.c tough2_file.c nc_cf_file.c)

set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_TP_LIBRARIES};${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYGLOT_LIBRARIES polyglot;${POLYGLOT_TP_LIBRARIES};${POLYMEC_LIBRARIES})

#add_subdirectory(tests)

# Installation targets.
install(TARGETS polyglot DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polyglot)

