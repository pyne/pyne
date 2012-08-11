# xs
set_source_files_properties(${PROJECT_SOURCE_DIR}/pyne/xs/Models.pyx PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module( models models.pyx models.cpp )
target_link_libraries( models  pyne pyne_nucname )
