PROJECT (mainFluDAG)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

IF (NOT DEFINED FLUPRO)
	set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

SET(MOAB_INCLUDE ${MOAB_HOME}/include)
SET(MOAB_LIB_DIR ${MOAB_HOME}/lib)
SET(moablibs ${MOAB_LIB_DIR}/libdagmc.so ${MOAB_LIB_DIR}/libMOAB.so)

# Where to look for includes
INCLUDE_DIRECTORIES(${MOAB_INCLUDE})
ENABLE_LANGUAGE(Fortran)

#SET (fluka_file ${FLUPRO}/usermvax/mgdraw.f)
SET (fluka_includes ${FLUPRO}/flukapro)
SET (fluka_libdir ${FLUPRO})
#ADD_CUSTOM_COMMAND(OUTPUT mgdraw.o
#       COMMAND gfortran -c ${fluka_file} -I${fluka_includes} 
#       DEPENDS ${FLUPRO}/usermvax/mgdraw.f
#       COMMENT "Custom command to compile mgdraw object with parenthetically named includes"
#       )
#ADD_CUSTOM_TARGET (mgdraw
#	        DEPENDS mgdraw.o)

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag mainFluDAG.cpp fluka_funcs.cpp)
TARGET_LINK_LIBRARIES(mainfludag ${moablibs} ${FLUPRO}/libflukahp.a gfortran)
