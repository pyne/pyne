set(PYNE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_IS_AMALGAMATION")

add_lib_to_pyne(pyne jsoncpp.cpp pyne.cpp)
install_lib(pyne)

add_lib_to_pyne(pyne_nucname nucname.cpp)
target_link_libraries(pyne_nucname ${PYNE_LIBRARIES})
install_lib(pyne_nucname)

add_lib_to_pyne(pyne_data data.cpp)
target_link_libraries(pyne_data ${PYNE_LIBRARIES})
install_lib(pyne_data)

add_lib_to_pyne(pyne_material material.cpp)
target_link_libraries(pyne_material ${PYNE_LIBRARIES})
install_lib(pyne_material)

add_lib_to_pyne(pyne_enrichment enrichment_cascade.cpp enrichment_symbolic.cpp enrichment.cpp)
target_link_libraries(pyne_enrichment material)
install_lib(pyne_enrichment)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
