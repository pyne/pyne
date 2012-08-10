macro add_lib_to_pyne( _name _source )
# add the library
add_library(${_name} ${_source})

# add it to the list of pyne libraries
set(PYNE_LIBRARIES ${PYNE_LIBRARIES} ${_name})
endmacro()

macro install_lib( _name )
# install it and add it to cpack
install(TARGETS ${_name}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  component libs
  )
endmacro()
