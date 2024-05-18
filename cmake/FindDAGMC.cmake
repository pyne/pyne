# - Try to find DAGMC
# Once done this will define
#
#  DAGMC_CMAKE_CONFIG - The path to the DAGMCConfig.cmake file
#  DAGMC_FOUND - system has DAGMC
#  DAGMC_VERSION - the version of DAGMC found
#  DAGMC_INCLUDE_DIRS - the DAGMC include directory
#  DAGMC_LIBRARY_DIRS - the DAGMC library directory
#  DAGMC_LIBRARIES - Link these to use DAGMC
#  DAGMC_DEFINITIONS - Compiler switches required for using DAGMC
#
#  Copyright (c) 2010 Roman Putanowicz <putanowr@l5.pk.edu.pl>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

message(STATUS "DAGMC_ROOT: ${DAGMC_ROOT}")

# Find DAGMC cmake config
find_path(DAGMC_CMAKE_CONFIG NAMES DAGMCConfig.cmake
          HINTS ${DAGMC_ROOT}
          PATH_SUFFIXES lib Lib cmake dagmc lib/cmake/dagmc
          NO_DEFAULT_PATH
          )

# Include DAGMC config
if(DAGMC_CMAKE_CONFIG)
  set(DAGMC_CMAKE_CONFIG ${DAGMC_CMAKE_CONFIG}/DAGMCConfig.cmake)
  include(${DAGMC_CMAKE_CONFIG})
  message(STATUS "DAGMC_CMAKE_CONFIG: ${DAGMC_CMAKE_CONFIG}")
  message(STATUS "DAGMC Include directory: ${DAGMC_INCLUDE_DIRS}")
  message(STATUS "DAGMC Library directories: ${DAGMC_LIBRARY_DIRS}")
  message(STATUS "DAGMC Libraries: ${DAGMC_LIBRARIES}")
  message(STATUS "DAGMC library version: ${DAGMC_VERSION}")
else()
  message(WARNING "Could not find DAGMC config file. "
    "Set -DDAGMC_ROOT=/path/to/DAGMC when running cmake or use the $DAGMC_ROOT environment variable."
    )
endif()