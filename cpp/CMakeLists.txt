

add_lib_to_pyne(pyne pyne.cpp)
install_lib(pyne)

add_lib_to_pyne(nucname nucname.cpp)
target_link_libraries(nucname ${PYNE_LIBRARIES})
install_lib(nucname)

add_lib_to_pyne(ace ace.cpp)
target_link_libraries(nucname ${PYNE_LIBRARIES})
install_lib(ace)


