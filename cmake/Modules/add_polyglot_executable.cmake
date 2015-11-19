function(add_polyglot_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYGLOT_LIBRARIES})
endfunction(add_polyglot_executable)

