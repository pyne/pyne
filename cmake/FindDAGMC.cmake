# - Try to find DAGMC
# Once done this will define
#
#  DAGMC_FOUND - system has DAGMC
#  DAGMC_INCLUDE_DIRS - the DAGMC include directory
#  DAGMC_LIBRARIES - Link these to use DAGMC
#

if (DAGMC_LIBRARIES AND DAGMC_INCLUDE_DIRS)
  # in cache already
  set(DAGMC_FOUND TRUE)
else (DAGMC_LIBRARIES AND DAGMC_INCLUDE_DIRS)
  find_path(DAGMC_INCLUDE_DIR NAMES DagMC.hpp
    PATHS $ENV{HOME}/.local/include ${HOME}/opt/dagmc
        PATH_SUFFIXES include Include
    PATHS "${BASE_DIR}/include" 
    ENV DAGMC_ROOT
    NO_DEFAULT_PATH
    )
  find_path(DAGMC_INCLUDE_DIR NAMES DagMC.hpp)

  find_library(DAGMC_LIBRARY NAMES DAGMC 
    PATHS $ENV{HOME}/.local/lib 
        PATH_SUFFIXES lib Lib
    PATHS "${BASE_DIR_LIB}" 
    ENV DAGMC_ROOT
    NO_DEFAULT_PATH
    )
  find_library(DAGMC_LIBRARY NAMES DagMC)

  set(DAGMC_INCLUDE_DIRS
    ${DAGMC_INCLUDE_DIR} CACHE PATH "Path to DAGMC headers")

  set(DAGMC_LIBRARIES
      ${DAGMC_LIBRARY} CACHE STRING "Directories to be linked to use DAGMC")

  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set DAGMC_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(DAGMC  DEFAULT_MSG
                                  DAGMC_LIBRARY DAGMC_INCLUDE_DIR)
  if (DAGMC_FOUND)
    message(STATUS "DAGMC header files: ${DAGMC_INCLUDE_DIR}")
    message(STATUS "DAGMC library: ${DAGMC_LIBRARY}")
  endif (DAGMC_FOUND)
  mark_as_advanced(DAGMC_INCLUDE_DIRS DAGMC_LIBRARIES)
endif (DAGMC_LIBRARIES AND DAGMC_INCLUDE_DIRS)

