cmake_minimum_required(VERSION 2.6.4)

# project name
project (dagsolid)


# make the dagsolid lib
add_library (dagsolid SHARED /data/opt/andy_dag/DAGMC/Geant4/dagsolid/src/DagSolid.cc)

# set moab paths
set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIBS ${MOAB_HOME}/lib)

set(GEANT_INCLUDE ${GEANT_HOME}/include/Geant4/)

# include files
include_directories(/data/opt/andy_dag/DAGMC/Geant4/dagsolid/include ${GEANT_INCLUDE} ${MOAB_INCLUDE} ${HDF5_INCLUDE})

# set libs
set(moab_libs ${MOAB_HOME}/lib/libdagmc.so ${MOAB_HOME}/lib/libMOAB.so)

# targets
install(TARGETS dagsolid DESTINATION lib)
install(FILES dagsolid.hh DESTINATION include) 