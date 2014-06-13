set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_library(pyne 
    "dagmc_bridge.cpp"
    "data.cpp"
    "enrichment.cpp"
    "enrichment_cascade.cpp"
    "enrichment_symbolic.cpp"
    "jsoncpp.cpp"
    "material.cpp"
    "nucname.cpp"
    "rxname.cpp"
    "tally.cpp"
    "utils.cpp"
    )
target_link_libraries(pyne hdf5)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
