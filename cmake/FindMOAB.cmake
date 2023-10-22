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

if (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)
  # in cache already
  set(MOAB_FOUND TRUE)
else (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)
  find_path(MOAB_INCLUDE_DIR NAMES MBiMesh.hpp moab_export.h
            HINTS ${MOAB_ROOT}/include ${DEPS_INCLUDE_HINTS})

  find_library(MOAB_LIBRARY NAMES MOAB
              HINTS ${MOAB_ROOT}/lib ${DEPS_LIB_HINTS})

  set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIR} CACHE PATH "Path to MOAB headers")
  set(MOAB_LIBRARIES ${MOAB_LIBRARY} CACHE STRING "Directories to be linked to use MOAB")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MOAB  DEFAULT_MSG MOAB_LIBRARY MOAB_INCLUDE_DIRS)

  if (MOAB_FOUND)
    message(STATUS "MOAB header files: ${MOAB_INCLUDE_DIRS}")
    message(STATUS "MOAB library: ${MOAB_LIBRARY}")
  endif (MOAB_FOUND)

  mark_as_advanced(MOAB_INCLUDE_DIRS MOAB_LIBRARIES)

endif (MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)

