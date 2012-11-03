# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})

# STL converters
# If the pyx file is a C++ file, we should specify that here.
# then, add the module
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/stlconverters.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(stlconverters stlconverters.pyx)
                            
# jsoncpp
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/jsoncpp.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(jsoncpp jsoncpp.pyx)
target_link_libraries(jsoncpp dl pyne)

# pyne_config
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/pyne_config.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(pyne_config  pyne_config.pyx)
target_link_libraries(pyne_config dl pyne)

# _utils
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/_utils.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(_utils _utils.pyx)

# nucname
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/nucname.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(nucname  nucname.pyx)
target_link_libraries(nucname dl  pyne)

# data
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/data.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(data data.pyx)
target_link_libraries(data dl pyne pyne_nucname pyne_data hdf5)

# material
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/material.pyx 
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(material  material.pyx)
target_link_libraries(material dl pyne pyne_nucname pyne_data pyne_material hdf5)

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
