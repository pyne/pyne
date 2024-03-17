set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# Build normal cython modules
set(PYNE_CYTHON_MODULES
    ace
    bins
    data
    dtypes
    endf
    enrichment
    extra_types
    jsoncpp
    material
    material_library
    nucname
    particle
    pyne_config
    rxname
    source
    stlcontainers
    tally
    transmuters
    _utils
    )
if(MOAB_FOUND)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} source_sampling)
endif(MOAB_FOUND)

message(STATUS "Modules to Cythonize: ${PYNE_CYTHON_MODULES}")

foreach(cymod ${PYNE_CYTHON_MODULES})
  set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/${cymod}.pyx"
                              PROPERTIES CYTHON_IS_CXX TRUE)
  add_cython_target(${cymod} CXX)
  add_library(${cymod} MODULE ${cymod})
  python_extension_module(${cymod})
  target_link_libraries(${cymod} pyne)
  install(TARGETS ${cymod} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")
endforeach()

#
# Cram is actually a C, not C++ file
#
add_cython_target(cram C)
add_library(cram MODULE ${cram})
python_extension_module(cram)
target_link_libraries(cram pyne)
install(TARGETS cram LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

#
# dagmc requires a rename and is only built if we have MOAB.
#
if(DAGMC_FOUND)
  set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/dagmc.pyx"
                              PROPERTIES CYTHON_IS_CXX TRUE)
  add_cython_target(_dagmc CXX "${PROJECT_SOURCE_DIR}/src/dagmc_bridge.cpp")
  add_library(_dagmc MODULE ${_dagmc})
  python_extension_module(_dagmc)
  target_link_libraries(_dagmc dagmc MOAB pyne)
  set_target_properties(_dagmc PROPERTIES OUTPUT_NAME dagmc)
  install(TARGETS _dagmc LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")
endif(DAGMC_FOUND)

#
# Spatial Solver
#
message("-- F2PY F90 FLAGS: ${F2PY_F90FLAGS}")

if(BUILD_SPATIAL_SOLVER)
  add_custom_target(transport_spatial_methods ALL
                    DEPENDS transport_spatial_methods_so pyne)
  add_custom_command(OUTPUT transport_spatial_methods_so
    COMMAND ${F2PY_EXECUTABLE} -c -m transport_spatial_methods
            -I${PROJECT_BINARY_DIR}/src -L${PROJECT_BINARY_DIR}/src
            ${F2PY_F90FLAGS}
            --f90flags="-fdefault-real-8"
            --no-wrap-functions
            -lpyne
            ${PROJECT_SOURCE_DIR}/src/transport_spatial_methods/3d/main.f90
    DEPENDS pyne
    )
endif()
