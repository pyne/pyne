macro( add_lib_to_pyne _name _source )
  # add the library
  add_library(${_name}             ${_source})
  # add it to the list of pyne libraries
  set(PYNE_LIBRARIES ${PYNE_LIBRARIES} ${_name})
endmacro()

macro( install_lib _name )
  # install it
  install(TARGETS ${_name}
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
endmacro()

macro( print_logo )
  set(cat_prog cat)
  if(WIN32)
    set(cat_prog type)
  endif(WIN32)
  exec_program(${cat_prog} ARGS ${PROJECT_SOURCE_DIR}/cmake/logo.txt 
    OUTPUT_VARIABLE variable)
  message("${variable}")
endmacro()
