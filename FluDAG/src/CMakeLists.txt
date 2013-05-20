CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT (MyFluDAG)

IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 
IF (NOT DEFINED FLUPRO)
	set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

IF (NOT DEFINED HOME OR NOT DEFINED FLUPRO)
	message ("Including local.cmake")
	INCLUDE(local.cmake)
ENDIF (NOT DEFINED HOME OR NOT DEFINED FLUPRO)

IF (NOT DEFINED MOAB_HOME)
	message ("Including local.cmake")
	INCLUDE(local.cmake)
ENDIF (NOT DEFINED MOAB_HOME)

SET(MOAB_INCLUDE ${MOAB_HOME}/include)
SET(MOAB_LIB_DIR ${MOAB_HOME}/lib)
SET (moablibs ${MOAB_LIB_DIR}/libMOAB.so ${MOAB_LIB_DIR}/libdagmc.so)


# Where to look for includes
include_directories(include ${MOAB_INCLUDE})

SET (commonfiles cpp/fluka_funcs.cpp cpp/dagmc_utils.cpp cpp/Input.cpp)
SET (wrapperfiles cpp/WrapReg2Name.cc cpp/WrapLookMG.cc cpp/WrapMag.cc 
		  cpp/WrapLookFX.cc cpp/WrapReg.cc  
		  cpp/WrapFlgfwr.cc 
		  )

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag cpp/mainFluDAG.cpp ${commonfiles} ${wrapperfiles})
SET (CMAKE_Fortran_Compiler gfortran)
ENABLE_LANGUAGE(Fortran)

TARGET_LINK_LIBRARIES(mainfludag ${moablibs} ${FLUPRO}/libflukahp.a gfortran)
