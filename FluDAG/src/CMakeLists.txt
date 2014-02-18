# CMAKE build script for FluDag 

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
PROJECT (MyFluDAG)

IF (NOT DEFINED HOME) 
	SET (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

IF (NOT DEFINED FLUPRO)
	SET (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include", "${FLUPRO}/flukapro")

# make fludag source files into a library
ADD_LIBRARY(fludag SHARED
	    ${DAGMC_FLUDAG_SOURCE}/fluka_funcs.cpp
	    ${DAGMC_FLUDAG_SOURCE}/dagmc_utils.cpp
	    ${DAGMC_FLUDAG_SOURCE}/WrapMag.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapLookFX.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapReg.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapFlgfwr.cc

)

# set libraries needed to link to test cases
SET(FLUDAG_LIBRARIES
    pthread
    gfortran
    fludag
    ${FLUPRO}/libflukahp.a
    ${MOAB_HOME}/lib/libMOAB.so
    ${MOAB_HOME}/lib/libdagmc.so
    ${MOAB_HOME}/lib/libMOAB.so
)

ENABLE_LANGUAGE(Fortran)

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag mainFluDAG.cpp) 
TARGET_LINK_LIBRARIES(mainfludag ${FLUDAG_LIBRARIES})


