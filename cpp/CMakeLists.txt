set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(PYNE_SRCS
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
if (MOAB_FOUND)
    set(PYNE_SRCS "${PYNE_SRCS}" "dagmc_bridge.cpp")    
endif (MOAB_FOUND)

# compile and link library
add_library(pyne ${PYNE_SRCS})
target_link_libraries(pyne hdf5)
if (MOAB_FOUND)
    target_link_libraries(pyne dagmc MOAB)
endif (MOAB_FOUND)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
