# Set PYNE version
set(PYNE_MAJOR_VERSION 0)
set(PYNE_MINOR_VERSION 7)
set(PYNE_PATCH_VERSION 6)
set(PYNE_VERSION ${PYNE_MAJOR_VERSION}.${PYNE_MINOR_VERSION}.${PYNE_PATCH_VERSION})

MESSAGE(STATUS "Setting PyNE Version to ${PYNE_VERSION} by changing files in ${PYNE_SUBMODULE_PATH}")

# Configure Pyne and cpp version headers
configure_file(${PYNE_SUBMODULE_PATH}pyne/pyne_version.py.in ${PYNE_SUBMODULE_PATH}pyne/pyne_version.py)
configure_file(${PYNE_SUBMODULE_PATH}src/pyne_version.h.in ${PYNE_SUBMODULE_PATH}src/pyne_version.h)

