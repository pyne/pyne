
add_lib_to_pyne(pyne pyne.cpp)
install_lib(pyne)

add_lib_to_pyne(pyne_nucname nucname.cpp)
target_link_libraries(pyne_nucname ${PYNE_LIBRARIES})
install_lib(pyne_nucname)

add_lib_to_pyne(pyne_data data.cpp)
target_link_libraries(pyne_data ${PYNE_LIBRARIES})
install_lib(pyne_data)

add_lib_to_pyne(pyne_material material.cpp)
target_link_libraries(pyne_material ${PYNE_LIBRARIES})
install_lib(pyne_material)

add_lib_to_pyne(pyne_enrichment enrichment.cpp)
target_link_libraries(pyne_enrichment material)
install_lib(pyne_enrichment)
