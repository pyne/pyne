
add_lib_to_pyne(pyne pyne.cpp)
install_lib(pyne)

add_lib_to_pyne(nucname nucname.cpp)
target_link_libraries(nucname ${PYNE_LIBRARIES})
install_lib(nucname)

add_lib_to_pyne(data data.cpp)
target_link_libraries(data ${PYNE_LIBRARIES})
install_lib(data)

add_lib_to_pyne(material material.cpp)
target_link_libraries(material ${PYNE_LIBRARIES})
install_lib(material)

add_lib_to_pyne(enrichment enrichment.cpp)
target_link_libraries(enrichment material)
install_lib(enrichment)
