# - Try to find MOAB
# Once done this will define
#
#  MOAB_FOUND - system has MOAB
#  MOAB_INCLUDE_DIRS - the MOAB include directory
#  MOAB_LIBRARIES - Link these to use MOAB
#  MOAB_DEFINITIONS - Compiler switches required for using MOAB
#
#  Copyright (c) 2010 Roman Putanowicz <putanowr@l5.pk.edu.pl>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)
  # in cache already
  set(MOAB_FOUND TRUE)
else (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)
  find_path(MOAB_INCLUDE_DIR NAMES MBiMesh.hpp
    HINTS ${MOAB_ROOT}/include ${DEPS_INCLUDE_HINTS}
    PATHS $ENV{HOME}/.local/include
    PATH_SUFFIXES include Include
    PATHS "${BASE_DIR}/include" "${BASE_DIR}/../install/include"
    ENV MOAB_ROOT
    NO_DEFAULT_PATH
    )
  find_path(MOAB_INCLUDE_DIR NAMES MBiMesh.hpp
            HINTS ${MOAB_ROOT}/include ${DEPS_INCLUDE_HINTS})

  find_library(MOAB_LIBRARY NAMES MOAB
    HINTS ${MOAB_ROOT}/lib ${DEPS_LIB_HINTS}
    PATHS $ENV{HOME}/.local/lib
        PATH_SUFFIXES lib Lib
    PATHS "${BASE_DIR_LIB}" "${BASE_DIR_LIB}/../../install/lib"
    ENV MOAB_ROOT
    NO_DEFAULT_PATH
    )
  find_library(MOAB_LIBRARY NAMES MOAB
               HINTS ${MOAB_ROOT}/lib ${DEPS_LIB_HINTS})

  set(MOAB_INCLUDE_DIRS
      ${MOAB_INCLUDE_DIR} CACHE PATH "Path to MOAB headers")

  set(MOAB_LIBRARIES
      ${MOAB_LIBRARY} CACHE STRING "Directories to be linked to use MOAB")

  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set MOAB_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(MOAB  DEFAULT_MSG
      MOAB_LIBRARY MOAB_INCLUDE_DIRS)
  if (MOAB_FOUND)
    message(STATUS "MOAB header files: ${MOAB_INCLUDE_DIRS}")
    message(STATUS "MOAB library: ${MOAB_LIBRARY}")
  endif (MOAB_FOUND)
  mark_as_advanced(MOAB_INCLUDE_DIRS MOAB_LIBRARIES)
endif (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)

