# Try to find MOAB
#
# Once done this will define
#
#  MOAB_FOUND - system has MOAB
#  MOAB_INCLUDE_DIRS - the MOAB include directory
#  MOAB_LIBRARIES - Link these to use MOAB
#  MOAB_DEFINITIONS - Compiler switches required for using MOAB

find_path(MOAB_CMAKE_CONFIG NAMES MOABConfig.cmake
          HINTS ${MOAB_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib cmake cmake/MOAB
          NO_DEFAULT_PATH)

message(STATUS "Found MOAB in ${MOAB_CMAKE_CONFIG}")

include(${MOAB_CMAKE_CONFIG}/MOABConfig.cmake)

string(STRIP "${DAGMC_LIBRARIES}" DAGMC_LIBRARIES)
string(REGEX REPLACE ";" " " DAGMC_LIBRARIES "${DAGMC_LIBRARIES}")
string(REGEX REPLACE "[ ]+" " " DAGMC_LIBRARIES "${DAGMC_LIBRARIES}")
separate_arguments(DAGMC_LIBRARIES)
list(REMOVE_DUPLICATES DAGMC_LIBRARIES)

message(STATUS "  MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
message(STATUS "  MOAB_LIBRARY_DIRS: ${MOAB_LIBRARY_DIRS}")
message(STATUS "  DAGMC_LIBRARIES: ${DAGMC_LIBRARIES}")
