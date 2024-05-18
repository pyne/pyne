# This code sets the following variables:
#
#  MOAB_CMAKE_CONFIG
#  MOAB_FOUND
#  MOAB_INCLUDE_DIRS
#  MOAB_LIBRARY_DIRS
#  MOAB_LIBRARIES

message(STATUS "MOAB_ROOT: ${MOAB_ROOT}")

# Find MOAB cmake config file
set(MOAB_SEARCH_DIRS)
file(GLOB MOAB_SEARCH_DIRS ${MOAB_SEARCH_DIRS} "${MOAB_ROOT}/lib*/cmake/MOAB")
string(REPLACE "\n" ";" MOAB_SEARCH_DIRS "${MOAB_SEARCH_DIRS}")
find_path(MOAB_CMAKE_CONFIG
  NAMES MOABConfig.cmake
  PATHS ${MOAB_SEARCH_DIRS}
  NO_DEFAULT_PATH
)

# Include MOAB config file
if (MOAB_CMAKE_CONFIG)
  set(MOAB_CMAKE_CONFIG ${MOAB_CMAKE_CONFIG}/MOABConfig.cmake)
  include(${MOAB_CMAKE_CONFIG})
  message(STATUS "MOAB_CMAKE_CONFIG: ${MOAB_CMAKE_CONFIG}")
  message(STATUS "MOAB Include directory: ${MOAB_INCLUDE_DIRS}")
  message(STATUS "MOAB Library directories: ${MOAB_LIBRARY_DIRS}")
  message(STATUS "MOAB Libraries: ${MOAB_LIBRARIES}")
  message(STATUS "MOAB library version: ${MOAB_VERSION}")
else ()
  message(WARNING "Could not find MOAB. Set -DMOAB_ROOT=/path/to/MOAB when running cmake or use the $MOAB_ROOT environment variable.")
endif ()
