set(HDF5_SO hdf5)
if(WIN32)
    set(HDF5_SO hdf5)
endif(WIN32)

macro( add_lib_to_pyne _name _source )
  # add the library
  add_library(${_name}             ${_source})
  # add it to the list of pyne libraries
  set(PYNE_LIBRARIES ${PYNE_LIBRARIES} ${_name})
endmacro()

macro( install_lib _name )
  # install it
  set(lib_type LIBRARY)
  if(WIN32)
    set(lib_type RUNTIME)
  endif(WIN32)
  install(TARGETS ${_name} ${lib_type} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endmacro()

macro( print_logo )
  set(cat_prog cat)
  if(WIN32)
    set(cat_prog type)
  endif(WIN32)
  message("PRINT CAT!!!")
  #execute_process(COMMAND ${cat_prog} ${PROJECT_SOURCE_DIR}\\cmake\\logo.txt)
  execute_process(COMMAND ${cat_prog} ${PROJECT_SOURCE_DIR}/cmake/logo.txt OUTPUT_VARIABLE variable)
  message("var = ${variable}")
endmacro()
