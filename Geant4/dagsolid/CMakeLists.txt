######################################################################
#
# CMakeLists.txt for dagmc/Geant4/dagsolid
#
# $Id: CMakeLists.txt 6346 2014-01-27 03:33:21Z cary $
#
######################################################################

if (NOT DEFINED CMAKE_PROJECT_NAME)
  PROJECT(DagSolid)
  CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
  set(TOP_IS_DAGSOLID TRUE)
  if (DAGSOLID_DIR)
    message(STATUS "Setting CMAKE_INSTALL_PREFIX.")
    get_filename_component(CMAKE_INSTALL_PREFIX ${DAGSOLID_DIR}/../..
      REALPATH)
  endif ()
else ()
  set(TOP_IS_DAGSOLID FALSE)
endif ()
message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}.")
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}.")

# Fill in includes that are found in two different ways
if (MOAB_DIR AND NOT Moab_INCLUDE_DIRS)
  set(Moab_INCLUDE_DIRS ${MOAB_DIR}/include)
endif ()
if (GEANT_DIR AND NOT Geant4_INCLUDE_DIRS)
  set(Geant4_INCLUDE_DIRS ${GEANT_DIR}/include/Geant4)
endif ()
# hdf5 last as most problematic
if (NOT Hdf5_INCLUDE_DIRS)
  if (HDF5_DIR AND NOT Hdf5_INCLUDE_DIRS)
    set(Hdf5_INCLUDE_DIRS ${HDF5_DIR}/include)
  elseif (EXISTS ${CMAKE_INSTALL_PREFIX}/../hdf5)
    get_filename_component(Hdf5_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/../../../hdf5/include REALPATH)
  elseif (EXISTS ${Geant4_INCLUDE_DIRS}/../../../hdf5-sersh)
    get_filename_component(Hdf5_INCLUDE_DIRS
        ${Geant4_INCLUDE_DIRS}/../../../hdf5-sersh/include REALPATH)
  endif ()
endif ()
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIRS} ${Moab_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}
)

# Find libraries that can be found in two different ways
if (MOAB_DIR AND NOT Moab_MOAB_LIBRARY)
  find_library(Moab_MOAB_LIBRARY MOAB PATH ${MOAB_DIR}/lib NO_DEFAULT_PATH)
  find_library(Moab_dagmc_LIBRARY dagmc PATH ${MOAB_DIR}/lib NO_DEFAULT_PATH)
endif ()
if (GEANT_DIR AND NOT Geant4_G4geometry_LIBRARY)
  message(STATUS "Looking for library, G4geometry.")
  find_library(Geant4_G4geometry_LIBRARY G4geometry
    PATH ${GEANT_DIR}/lib64 ${GEANT_DIR}/lib NO_DEFAULT_PATH)
  find_library(Geant4_G4clhep_LIBRARY G4clhep
    PATH ${GEANT_DIR}/lib64 ${GEANT_DIR}/lib NO_DEFAULT_PATH)
  find_library(Geant4_G4global_LIBRARY G4global
    PATH ${GEANT_DIR}/lib64 ${GEANT_DIR}/lib NO_DEFAULT_PATH)
endif ()
foreach (var Hdf5_INCLUDE_DIRS Moab_INCLUDE_DIRS Geant4_INCLUDE_DIRS
  Moab_MOAB_LIBRARY Moab_dagmc_LIBRARY
  Geant4_G4geometry_LIBRARY Geant4_G4clhep_LIBRARY Geant4_G4global_LIBRARY)
  message (STATUS "${var} = ${${var}}.")
endforeach ()

# So we can see the levels of addition
set(SCI_INDENT "${SCI_INDENT}  ")
# disable threading
add_definitions(-DGTEST_HAS_PTHREAD=0)
message(STATUS "${SCI_INDENT}Adding src subdir.")
ADD_SUBDIRECTORY(src)
if (TOP_IS_DAGSOLID)
  message(STATUS "${SCI_INDENT}Adding test subdir.")
  ADD_SUBDIRECTORY(test)
endif ()

