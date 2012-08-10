

add_lib_to_pyne(pyne pyne.cpp)
install_lib(pyne)

add_lib_to_pyne(nucname nucname.cpp)
target_link_libraries(nucname pyne)
install_lib(nucname)

