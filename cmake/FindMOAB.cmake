# This code sets the following variables:
#
#  MOAB_FOUND
#  MOAB_INCLUDE_DIRS
#  MOAB_LIBRARIES

# Find MOAB cmake config file
set(MOAB_SEARCH_DIRS)
file(GLOB MOAB_SEARCH_DIRS ${MOAB_SEARCH_DIRS} "${MOAB_ROOT}/lib*/cmake/MOAB")
string(REPLACE "\n" ";" MOAB_SEARCH_DIRS "${MOAB_SEARCH_DIRS}")
find_path(MOAB_CMAKE_CONFIG
  NAMES MOABConfig.cmake
  PATHS ${MOAB_SEARCH_DIRS}
  NO_DEFAULT_PATH
)
if (MOAB_CMAKE_CONFIG)
  set(MOAB_CMAKE_CONFIG ${MOAB_CMAKE_CONFIG}/MOABConfig.cmake)
  message(STATUS "MOAB_CMAKE_CONFIG: ${MOAB_CMAKE_CONFIG}")
  # Include MOAB config
  include(${MOAB_CMAKE_CONFIG})
else ()
  message(STATUS "Could not find MOABConfig.cmake in ${MOAB_SEARCH_DIRS}.")
  set(MOAB_FOUND FALSE)
endif ()


