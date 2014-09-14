##---------------------------------------------------------------------------##
## File  : pyne/cmake/FindDagMC.cmake
## Author: Seth R Johnson
## Date  : Sun Sep 14 09:34:03 2014
##
## Find the location of the dagmc libraries.
#
# This module finds the DagMC library and header file.
# It creates the following CMAKE variables:
#
#  DAGMC_FOUND             - Whether the DagMC installation was found
#  DAGMC_LIBRARY           - Path to the DagMC library
#  IMESH_LIBRARY           - Path to the iMesh library
#  DAGMC_INCLUDE_DIR       - Path to the DagMC include directory
#
#  DAGMC_INCLUDE_DIRS      - Path to the DagMC include directory
#  DAGMC_LIBRARIES         - DagMC library target
#
##---------------------------------------------------------------------------##

find_path(DAGMC_INCLUDE_DIR DagMC.hpp)

# handle the QUIETLY and REQUIRED arguments and set DAGMC_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)

find_library(DAGMC_LIBRARY NAMES dagmc
  PATHS $ENV{HOME}/.local/lib
  PATHS "${BASE_DIR_LIB}" "${BASE_DIR_LIB}/../../install/lib"
  )

find_library(IMESH_LIBRARY NAMES iMesh
  PATHS $ENV{HOME}/.local/lib
  PATHS "${BASE_DIR_LIB}" "${BASE_DIR_LIB}/../../install/lib"
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(DagMC
  DAGMC_INCLUDE_DIR DAGMC_LIBRARY IMESH_LIBRARY)

if(DAGMC_FOUND)
  mark_as_advanced(DAGMC_INCLUDE_DIR DAGMC_LIBRARY IMESH_LIBRARY)
  set(DAGMC_INCLUDE_DIRS ${DAGMC_INCLUDE_DIR})
  set(DAGMC_LIBRARIES ${DAGMC_LIBRARY} ${IMESH_LIBRARY})
endif()

##---------------------------------------------------------------------------##
## end of pyne/cmake/FindDagMC.cmake
##---------------------------------------------------------------------------##
