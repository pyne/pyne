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
  ADD_SUBDIRECTORY(src)
  ADD_SUBDIRECTORY(test)
else ()
  set(SCI_INDENT "${SCI_INDENT}  ")
  message(STATUS "${SCI_INDENT}Adding src subdir.")
  ADD_SUBDIRECTORY(src)
endif ()

# disable threading
add_definitions(-DGTEST_HAS_PTHREAD=0)

