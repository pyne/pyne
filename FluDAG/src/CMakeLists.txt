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
	    ${FLUDAG_SOURCE}/fluka_funcs.cpp
	    ${FLUDAG_SOURCE}/dagmc_utils.cpp
	    ${FLUDAG_SOURCE}/WrapMag.cc
	    ${FLUDAG_SOURCE}/WrapLookFX.cc
	    ${FLUDAG_SOURCE}/WrapReg.cc
	    ${FLUDAG_SOURCE}/WrapFlgfwr.cc
)
# ADD_LIBRARY(parseflu SHARED
#	    ${DAGMC_FLUDAG_SOURCE}/dagmc_define.cpp
#	    ${DAGMC_FLUDAG_SOURCE}/parse_flu.cpp
#)

# set libraries needed to link to test cases
SET(FLUDAG_LIBRARIES
    pthread
    ${MOAB_HOME}/lib/libdagmc.so
    ${MOAB_HOME}/lib/libMOAB.so
    ${FLUPRO}/libflukahp.a
    fludag
    gfortran
)

ENABLE_LANGUAGE(Fortran)

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag mainFluDAG.cpp) 
TARGET_LINK_LIBRARIES(mainfludag ${FLUDAG_LIBRARIES})

# ADD_EXECUTABLE(mainparse dagmc_define.cpp parse_flu.cpp moab_utils.cpp UnitNumberManager.cpp) 
# TARGET_LINK_LIBRARIES(mainparse ${FLUDAG_LIBRARIES})
