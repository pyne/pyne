# Python extension modules
include_directories(${PYNE_INCLUDE_DIRS})
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# xs
set_source_files_properties("${PROJECT_SOURCE_DIR}/pyne/xs/models.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(models models.pyx)
target_link_libraries(models pyne)
install(TARGETS models LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pyne/xs")
