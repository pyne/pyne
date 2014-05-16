######################################################################
#
# CMakeLists.txt for dagmc/Geant4/dagsolid/src
#
# $Id: CMakeLists.txt 6346 2014-01-27 03:33:21Z cary $
#
######################################################################

add_library(dagsolid SHARED DagSolid.cc)
if (VERSION_MAJOR)
  set_target_properties(dagsolid PROPERTIES
    SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    COMPILE_DEFINITIONS dagsolid_EXPORTS
  )
endif ()
# Mac needs to link shared libs to shared libs
if (APPLE)
  target_link_libraries(dagsolid
    ${Moab_MOAB_LIBRARY}
    ${Moab_dagmc_LIBRARY}
    ${Geant4_G4geometry_LIBRARY}
    ${Geant4_G4clhep_LIBRARY}
    ${Geant4_G4global_LIBRARY}
  )
endif ()

install(TARGETS dagsolid DESTINATION Geant4/dagsolid/lib)

