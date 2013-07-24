# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# extra_types
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/extra_types.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(extra_types extra_types.pyx)

# STL containers
# If the pyx file is a C++ file, we should specify that here.
# then, add the module
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/stlcontainers.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(stlcontainers stlcontainers.pyx)
                            
# jsoncpp
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/jsoncpp.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(jsoncpp jsoncpp.pyx)
target_link_libraries(jsoncpp pyne)

# pyne_config
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/pyne_config.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(pyne_config  pyne_config.pyx)
target_link_libraries(pyne_config pyne)

# _utils
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/_utils.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(_utils _utils.pyx)
target_link_libraries(_utils pyne)

# endf
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/endf.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(endf endf.pyx)
target_link_libraries(endf pyne_nucname)

# nucname
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/nucname.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(nucname nucname.pyx)
target_link_libraries(nucname pyne pyne_nucname)

# rxname
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/rxname.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(rxname rxname.pyx)
target_link_libraries(rxname pyne pyne_rxname)

# data
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/data.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(data data.pyx)
target_link_libraries(data pyne pyne_nucname pyne_data hdf5)

# material
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/material.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(material  material.pyx)
target_link_libraries(material pyne pyne_nucname pyne_data pyne_material hdf5)

# ace
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/ace.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(ace ace.pyx)
target_link_libraries(ace pyne pyne_nucname)

# enrichment
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/enrichment.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(enrichment enrichment.pyx)
target_link_libraries(enrichment pyne_enrichment pyne_material)
