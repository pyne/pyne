# Set PYNE version
set(PYNE_MAJOR_VERSION 0)
set(PYNE_MINOR_VERSION 7)
set(PYNE_PATCH_VERSION 6)
set(PYNE_VERSION ${PYNE_MAJOR_VERSION}.${PYNE_MINOR_VERSION}.${PYNE_PATCH_VERSION})

if( NOT (PYNE_SUBMODULE_PATH) )
   MESSAGE(STATUS "PyNE not used as submodule, using current directory")
   set(PYNE_SUBMODULE_PATH "./")
endif(NOT (PYNE_SUBMODULE_PATH))

# Configure Pyne and cpp version headers
configure_file(${PYNE_SUBMODULE_PATH}/pyne/pyne_version.py.in ${CMAKE_CURRENT_SOURCE_DIR}/${PYNE_SUBMODULE_PATH}/pyne/pyne_version.py)
configure_file(${PYNE_SUBMODULE_PATH}/src/pyne_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/${PYNE_SUBMODULE_PATH}/src/pyne_version.h)

MESSAGE(STATUS "Setting PyNE Version to ${PYNE_VERSION} by changing files in ${CMAKE_CURRENT_SOURCE_DIR}/${PYNE_SUBMODULE_PATH}")
