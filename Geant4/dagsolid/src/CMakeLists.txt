cmake_minimum_required(VERSION 2.6.4)

# project name
project (dagsolid)

# make the dagsolid lib
add_library (dagsolid SHARED ../src/DagSolid.cc)

# set moab paths
set(MOAB_INCLUDE ${MOAB_DIR}/include)
set(MOAB_LIBS ${MOAB_DIR}/lib)

# set geant paths
set(GEANT_INCLUDE ${GEANT_DIR}/include/Geant4/ )

# include files
include_directories(../include ${GEANT_INCLUDE} ${MOAB_INCLUDE} ${HDF5_INCLUDE})

# set libs
set(moab_libs ${MOAB_DIR}/lib/libdagmc.so ${MOAB_DIR}/lib/libMOAB.so)

# targets
set (CMAKE_INSTALL_PREFIX .)
install(TARGETS dagsolid DESTINATION ${DAGSOLID_DIR}/lib)

