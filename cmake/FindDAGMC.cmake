# - Try to find DAGMC
# Once done this will define
#
#  DAGMC_FOUND - system has DAGMC
#  DAGMC_INCLUDE_DIRS - the DAGMC include directory
#  DAGMC_LIBRARIES - Link these to use DAGMC
#  DAGMC_DEFINITIONS - Compiler switches required for using DAGMC
#
#  Copyright (c) 2010 Roman Putanowicz <putanowr@l5.pk.edu.pl>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#


message(STATUS ${DAGMC_ROOT})

find_path(DAGMC_CMAKE_CONFIG NAMES DAGMCConfig.cmake
          HINTS ${DAGMC_ROOT}
          PATH_SUFFIXES lib Lib cmake lib/cmake/
          NO_DEFAULT_PATH)

message(STATUS "Found DAGMC in ${DAGMC_CMAKE_CONFIG}")

if( "DAGMC_CMAKE_CONFIG-NOTFOUND" STREQUAL "${DAGMC_CMAKE_CONFIG}")
  set(DAGMC_FOUND FALSE)
else()
  include(${DAGMC_CMAKE_CONFIG}/DAGMCConfig.cmake)
  set(DAGMC_FOUND TRUE)
endif()
