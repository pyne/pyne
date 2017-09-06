# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# extra_types
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/extra_types.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(extra_types extra_types.pyx)
install(TARGETS extra_types LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# dtypes
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/dtypes.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(dtypes dtypes.pyx)
install(TARGETS dtypes LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# STL containers
# If the pyx file is a C++ file, we should specify that here.
# then, add the module
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/stlcontainers.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(stlcontainers stlcontainers.pyx)
install(TARGETS stlcontainers LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# jsoncpp
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/jsoncpp.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(jsoncpp jsoncpp.pyx)
target_link_libraries(jsoncpp pyne)
install(TARGETS jsoncpp LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# pyne_config
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/pyne_config.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(pyne_config pyne_config.pyx)
target_link_libraries(pyne_config pyne)
install(TARGETS pyne_config LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# _utils
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/_utils.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(_utils _utils.pyx)
target_link_libraries(_utils pyne)
install(TARGETS _utils LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# dagmc
if (MOAB_FOUND)
  set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/dagmc.pyx"
                              PROPERTIES CYTHON_IS_CXX TRUE)
  cython_add_module(_dagmc dagmc.pyx "${PROJECT_SOURCE_DIR}/src/dagmc_bridge.cpp")
  target_link_libraries(_dagmc dagmc MOAB pyne)
  set_target_properties(_dagmc PROPERTIES OUTPUT_NAME dagmc)
  install(TARGETS _dagmc LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")
endif (MOAB_FOUND)

# endf
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/endf.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(endf endf.pyx)
target_link_libraries(endf pyne)
install(TARGETS endf LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# nucname
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/nucname.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(nucname nucname.pyx)
target_link_libraries(nucname pyne)
install(TARGETS nucname LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# rxname
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/rxname.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(rxname rxname.pyx)
target_link_libraries(rxname pyne)
install(TARGETS rxname LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# particle
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/particle.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(particle particle.pyx)
target_link_libraries(particle pyne)
install(TARGETS particle LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")


# data
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/data.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(data data.pyx)
target_link_libraries(data pyne)
install(TARGETS data LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# material
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/material.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(material material.pyx)
target_link_libraries(material pyne hdf5)
install(TARGETS material LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# ace
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/ace.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(ace ace.pyx)
target_link_libraries(ace pyne)
install(TARGETS ace LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# enrichment
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/enrichment.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(enrichment enrichment.pyx)
target_link_libraries(enrichment pyne)
install(TARGETS enrichment LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# tally
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/tally.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(tally tally.pyx)
target_link_libraries(tally pyne)
install(TARGETS tally LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

# source_sampling
if (MOAB_FOUND)
  set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/source_sampling.pyx"
                              PROPERTIES CYTHON_IS_CXX TRUE)
  cython_add_module(source_sampling source_sampling.pyx)
  target_link_libraries(source_sampling pyne)
  install(TARGETS source_sampling LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")
endif (MOAB_FOUND)

# bins
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/bins.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(bins bins.pyx)
target_link_libraries(bins pyne)
install(TARGETS bins LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne")

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
