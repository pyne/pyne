######################################################################
#
# CMakeLists.txt for dagmc/Geant4/dagsolid/src
#
# $Id: CMakeLists.txt 6346 2014-01-27 03:33:21Z cary $
#
######################################################################

if (NOT DEFINED CMAKE_PROJECT_NAME)
# Legacy file
  set(IS_TOPDIR TRUE)
  cmake_minimum_required(VERSION 2.6.4)
# project name
  project(dagsolid)

# set moab paths
  set(MOAB_INCLUDE ${MOAB_DIR}/include)
  set(MOAB_LIBS ${MOAB_DIR}/lib)
# set geant paths
  set(GEANT_INCLUDE ${GEANT_DIR}/include/Geant4/ )
  include_directories(../include ${GEANT_INCLUDE} ${MOAB_INCLUDE} ${HDF5_INCLUDE})
# Why the next line?
  set(moab_libs ${MOAB_DIR}/lib/libdagmc.so ${MOAB_DIR}/lib/libMOAB.so)
  add_library(dagsolid SHARED DagSolid.cc)
# Why the next line?
  set (CMAKE_INSTALL_PREFIX .)
  install(TARGETS dagsolid DESTINATION ${DAGSOLID_DIR}/lib)
else ()
  set(IS_TOPDIR FALSE)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${Geant4_INCLUDE_DIRS} ${Moab_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}
  )
  add_library(dagsolid DagSolid.cc)
  if (BUILD_SHARED_LIBS)
    set_target_properties(dagsolid PROPERTIES
      SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}
      VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
      COMPILE_DEFINITIONS dagsolid_EXPORTS
    )
    target_link_libraries(dagsolid
      ${Moab_MOAB_LIBRARY}
      ${Moab_dagmc_LIBRARY}
      ${Geant4_G4geometry_LIBRARY}
      ${Geant4_G4clhep_LIBRARY}
      ${Geant4_G4global_LIBRARY}
    )
  endif ()
  install(TARGETS dagsolid DESTINATION Geant4/lib)
endif ()

