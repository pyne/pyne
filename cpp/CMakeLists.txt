set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# pyne
add_library(pyne jsoncpp.cpp pyne.cpp)
set_target_properties(pyne PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
install_lib(pyne)

# pyne_nucname
add_library(pyne_nucname nucname.cpp)
set_target_properties(pyne_nucname PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
target_link_libraries(pyne_nucname pyne)
install_lib(pyne_nucname)

# pyne_rxname
add_library(pyne_rxname rxname.cpp)
set_target_properties(pyne_rxname PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
target_link_libraries(pyne_rxname pyne pyne_nucname)
install_lib(pyne_rxname)

# pyne_data
add_library(pyne_data data.cpp)
set_target_properties(pyne_data PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
target_link_libraries(pyne_data pyne_nucname pyne_rxname pyne hdf5)
install_lib(pyne_data)

# pyne_material
add_library(pyne_material material.cpp)
set_target_properties(pyne_material PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
target_link_libraries(pyne_material pyne pyne_nucname pyne_data hdf5)
install_lib(pyne_material)

# pyne_enrichment
add_library(pyne_enrichment enrichment_cascade.cpp 
                            enrichment_symbolic.cpp enrichment.cpp)
set_target_properties(pyne_enrichment PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
target_link_libraries(pyne_enrichment pyne_material)
install_lib(pyne_enrichment)

# pyne_tally
add_library(pyne_tally tally.cpp)
set_target_properties(pyne_tally PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyne/lib")
install_lib(pyne_tally)


# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

#file(COPY "${CMAKE_CURRENT_SOURCE_DIR}" DESTINATION 
#    "${CMAKE_BINARY_DIR}/pyne/includes/" FILES_MATCHING PATTERN "*.h")
