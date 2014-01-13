CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
PROJECT (MyFluDAG)

IF (NOT DEFINED HOME) 
	SET (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

IF (NOT DEFINED FLUPRO)
	SET (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

SET(MOAB_INCLUDE ${MOAB_HOME}/include)
SET(MOAB_LIB_DIR ${MOAB_HOME}/lib)
SET (moablibs ${MOAB_LIB_DIR}/libMOAB.so)
SET (moabdaglibs ${MOAB_LIB_DIR}/libdagmc.so ${moablibs})
SET (flukalib ${FLUPRO}/libflukahp.a)

# Where to look for includes
INCLUDE_DIRECTORIES(${MOAB_INCLUDE})
ENABLE_LANGUAGE(Fortran)

SET (fluka_includes ${FLUPRO}/flukapro)
SET (fluka_libdir ${FLUPRO})

# Keep these lines as a template for custom commands
#SET (fluka_file ${FLUPRO}/usermvax/mgdraw.f)
#ADD_CUSTOM_COMMAND(OUTPUT mgdraw.o
#       COMMAND gfortran -c ${fluka_file} -I${fluka_includes} 
#       DEPENDS ${FLUPRO}/usermvax/mgdraw.f
#       COMMENT "Custom command to compile mgdraw object with parenthetically named includes"
#       )
#ADD_CUSTOM_TARGET (mgdraw
#	        DEPENDS mgdraw.o)

SET (fludag_files fluka_funcs.cpp dagmc_utils.cpp )
SET (dagmc_define_files dagmc_define.cpp parse_flu.cpp chkerr.cpp moab_utils.cpp UnitNumberManager.cpp )
SET (wrapperfiles WrapMag.cc WrapLookFX.cc WrapReg.cc WrapFlgfwr.cc)

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag mainFluDAG.cpp ${fludag_files} ${wrapperfiles})
TARGET_LINK_LIBRARIES(mainfludag ${moabdaglibs} ${flukalib} gfortran)

# Create the executable that writes records
ADD_EXECUTABLE(mainparse ${dagmc_define_files} )
TARGET_LINK_LIBRARIES(mainparse ${moablibs} gfortran)

