######################################################################
#
# CMakeLists.txt for Geant4 of dagmc
#
# $Id: CMakeLists.txt 6426 2014-02-21 17:56:17Z cary $
#
# Copyright &copy; 2002-2012, University of Colorado and Tech-X Corporation
# See LICENSE file for conditions of use.
#
######################################################################

# Project information
project(dagmcg4)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8)

######################################################################
#
# Grab the scimake repo
#
######################################################################

if (EXISTS ${PROJECT_SOURCE_DIR}/scimake)
  message(STATUS "${PROJECT_SOURCE_DIR}/scimake found.  Update needed?")
else ()
  execute_process(COMMAND svn co https://svn.code.sf.net/p/scimake/code/trunk scimake
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
endif ()
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/scimake)
  message(FATAL_ERROR "scimake not found.")
endif ()

######################################################################
#
# Usual initialization stuff, sciInit includes txfindpackage and
# takes care of MPI
#
######################################################################

set(NOFORTRAN TRUE)
set(NO_CONFIG_H TRUE)
include(${PROJECT_SOURCE_DIR}/scimake/SciInit.cmake)
include_directories(${SCIMAKE_DIR}/include)
set(DAGMC_VERSION "${VERSION}-${PROJECT_REV}")
install(CODE "message(STATUS \"Installing into \${CMAKE_INSTALL_PREFIX}.\")")

######################################################################
#
# Always use rpath to greatest extent.
# See: http://www.itk.org/Wiki/CMake_RPATH_handling
# Add -DCMAKE_SKIP_RPATH:BOOL=TRUE to prevent any rpath handling
#
######################################################################

# SET(CMAKE_SKIP_BUILD_RPATH FALSE)
# Use, i.e. don't skip the full RPATH for the build tree
# Not needed since we build static libs

# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# When building, don't link with the install RPATH at build,
# but add later on when installing.
# Not needed since we build static libs

# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# The RPATH to be used when installing
# Not needed since we build static libs

if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()
# Add the automatically determined parts of the RPATH that
# point to directories outside the build tree to the install RPATH

###################################################
#
# Common packages
#
###################################################

find_package(SciHdf5 REQUIRED)
find_package(SciBoost REQUIRED thread chrono system)

###################################################
#
# Packages specific to geometry and monte carlo
#
###################################################

find_package(SciMoab REQUIRED)
find_package(SciGeant4 REQUIRED)

###################################################
#
# Add subdirectories
#
###################################################

message("")
set(SCI_INDENT "")
message(STATUS "Adding Geant4 subdir.")
add_subdirectory(Geant4)

###################################################
#
# Documentation
#
###################################################

find_package(SciDoxygen)
if (${DOXYGEN_FOUND})
  message(STATUS "Doxygen found.")
else ()
  message(WARNING "Doxygen not found.")
endif ()

###################################################
#
# Testing
#
###################################################

# Enable tests
# include(CTestConfig.cmake)
include(CTest)
enable_testing()

######################################################################
#
# Configure files
#
######################################################################

if (FALSE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DagMcConfig.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/DagMcConfig.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DagMcVersion.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/DagMcVersion.h)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/DagMcVersion.h
  ${CMAKE_CURRENT_BINARY_DIR}/DagMcConfig.h
  DESTINATION include
)
endif ()

######################################################################
#
# Packaging
#
######################################################################

# CPack version numbers for release tarball name.
set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH}-r${PROJECT_REV})

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DagMc Library")
set(CPACK_PACKAGE_VENDOR "No Vendor")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.rst)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "dagmc-${VERSION}-r${PROJECT_REV}"
    CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES
"/CVS/;/.svn/;.swp$;.#;/#;/build/;/serial/;/parallel/;~$;/autom4te.cache/;/.config;preconfig.out")
set(CPACK_GENERATOR TGZ)
include(CPack)

