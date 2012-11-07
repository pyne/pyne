set(PYNE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# pyne
add_library(pyne jsoncpp.cpp pyne.cpp)
install_lib(pyne)

# pyne_nucname
add_library(pyne_nucname nucname.cpp)
target_link_libraries(pyne_nucname pyne)
install_lib(pyne_nucname)

# pyne_data
add_library(pyne_data data.cpp)
target_link_libraries(pyne_data pyne pyne_nucname hdf5)
install_lib(pyne_data)

# pyne_material
add_library(pyne_material material.cpp)
target_link_libraries(pyne_material pyne pyne_nucname pyne_data hdf5)
install_lib(pyne_material)

# pyne_enrichment
add_library(pyne_enrichment enrichment_cascade.cpp 
                            enrichment_symbolic.cpp enrichment.cpp)
target_link_libraries(pyne_enrichment pyne_material)
install_lib(pyne_enrichment)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION 
#    ${CMAKE_BINARY_DIR}/pyne/includes/ FILES_MATCHING PATTERN "*.h")
