# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# extra_types
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/extra_types.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(extra_types extra_types.pyx)

# dtypes
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/dtypes.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(dtypes dtypes.pyx)

# STL containers
# If the pyx file is a C++ file, we should specify that here.
# then, add the module
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/stlcontainers.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(stlcontainers stlcontainers.pyx)

# jsoncpp
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/jsoncpp.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(jsoncpp jsoncpp.pyx)
target_link_libraries(jsoncpp pyne)

# pyne_config
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/pyne_config.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(pyne_config pyne_config.pyx)
target_link_libraries(pyne_config pyne)

# _utils
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/_utils.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(_utils _utils.pyx)
target_link_libraries(_utils pyne)

# dagmc
if (MOAB_FOUND)
    set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/dagmc.pyx"
                                PROPERTIES CYTHON_IS_CXX TRUE)
    cython_add_module(_dagmc dagmc.pyx "${PROJECT_SOURCE_DIR}/src/dagmc_bridge.cpp")
    target_link_libraries(_dagmc dagmc MOAB pyne)
    set_target_properties(_dagmc PROPERTIES OUTPUT_NAME dagmc)
endif (MOAB_FOUND)

# endf
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/endf.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(endf endf.pyx)
target_link_libraries(endf pyne)

# nucname
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/nucname.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(nucname nucname.pyx)
target_link_libraries(nucname pyne)

# rxname
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/rxname.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(rxname rxname.pyx)
target_link_libraries(rxname pyne)

# particle
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/particle.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(particle particle.pyx)
target_link_libraries(particle pyne)


# data
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/data.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(data data.pyx)
target_link_libraries(data pyne)

# material
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/material.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(material material.pyx)
target_link_libraries(material pyne hdf5)

# ace
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/ace.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(ace ace.pyx)
target_link_libraries(ace pyne)

# enrichment
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/enrichment.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(enrichment enrichment.pyx)
target_link_libraries(enrichment pyne)

# tally
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/tally.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(tally tally.pyx)
target_link_libraries(tally pyne)


# source_sampling
if (MOAB_FOUND)
    set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/source_sampling.pyx"
                                PROPERTIES CYTHON_IS_CXX TRUE)
    cython_add_module(source_sampling source_sampling.pyx)
    target_link_libraries(source_sampling pyne)
endif (MOAB_FOUND)

# bins
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/bins.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(bins bins.pyx)
target_link_libraries(bins pyne)

# figure out which f2py to use
find_program(F2PY_C f2py)
if(${PYTHON_VERSION_MAJOR} GREATER 2)
  find_program(F2PY_C2 f2py3)
  if(F2PY_C2)
    set(F2PY_C ${F2PY_C2})
  endif(F2PY_C2)
elseif(${PYTHON_VERSION_MAJOR} LESS 3)
  #because arch is dumb
  find_program(F2PY_C2 f2py2)
  if(F2PY_C2)
    set(F2PY_C ${F2PY_C2})
  endif(F2PY_C2)
endif(${PYTHON_VERSION_MAJOR} GREATER 2)



if(BUILD_SPATIAL_SOLVER)
    # TRANSPORT_SPATIAL_METHODS wrapper
    # I don't know which method is better...
    #add_custom_target(transport_spatial_methods ALL DEPENDS ${PROJECT_BINARY_DIR}/pyne/transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX})
    #add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/pyne/transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX}
    add_custom_target(transport_spatial_methods ALL DEPENDS transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_custom_command(OUTPUT transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX}
        COMMAND ${F2PY_C} -m transport_spatial_methods -I${PROJECT_BINARY_DIR}/src -L${PROJECT_BINARY_DIR}/src
                    --f90flags="-fdefault-real-8"
                    -c ${PROJECT_SOURCE_DIR}/src/transport_spatial_methods/3d/main.f90
                    -lpyne
                    #--debug --debug-capi
                    #-lpyne --debug-capi
        DEPENDS pyne
        )
    target_link_libraries(transport_spatial_methods pyne)
endif()
