# Try to find MOAB
#
# Once done this will define
#
# MOAB_FOUND        - boolean indicating that MOAB is found
# MOAB_INCLUDE_DIRS - include directories from which to pick up MOAB includes
# MOAB_LIBRARY_DIRS - location of MOAB libraries
# MOAB_LIBRARIES    - libraries need to link to MOAB
# DAGMC_LIBRARIES   - MOAB_LIBRARIES plus DAGMC library
# MOAB_CXX, MOAB_CC, MOAB_F77, MOAB_FC - compilers used to compile MOAB
# MOAB_CXXFLAGS, MOAB_CFLAGS, MOAB_FORTRAN_FLAGS - compiler flags used to compile MOAB

find_path(MOAB_CMAKE_CONFIG NAMES MOABConfig.cmake
          HINTS ${MOAB_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib cmake cmake/MOAB
          NO_DEFAULT_PATH)

if (MOAB_CMAKE_CONFIG)
  message(STATUS "Found MOAB in ${MOAB_CMAKE_CONFIG}")

  include(${MOAB_CMAKE_CONFIG}/MOABConfig.cmake)

  if (NOT MOAB_LIBRARY_DIRS)
    get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_INCLUDE_DIRS}/../lib ABSOLUTE)
  endif ()

  message(STATUS "  MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
  message(STATUS "  MOAB_LIBRARY_DIRS: ${MOAB_LIBRARY_DIRS}")
else ()
  set(MOAB_FOUND 0)
  message(STATUS "MOAB not found. Not building pyne with MOAB.")
endif ()
