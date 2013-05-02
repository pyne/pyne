cmake_minimum_required(VERSION 2.8)
project (MyFluDAG)

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

set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIB_DIR ${MOAB_HOME}/lib)
set (moablibs ${MOAB_LIB_DIR}/libMOAB.so ${MOAB_LIB_DIR}/libdagmc.so)


# Where to look for includes
include_directories(include ${MOAB_INCLUDE})

set (commonfiles cpp/readVol.cpp cpp/fluka_funcs.cpp cpp/dagmc_utils.cpp cpp/Input.cpp)

# add_executable(readvol cpp/readVol.cpp cpp/fluka_funcs.cpp)
add_executable(readvol cpp/mainReadVol.cpp ${commonfiles})
target_link_libraries(readvol ${moablibs})

# This is not used yet, but may come into play when cmake can be used to link the fluka library
set (wrapperfiles cpp/WrapReg2Name.cc cpp/WrapLookMG.cc cpp/WrapMag.cc 
		  cpp/WrapLookFX.cc cpp/WrapReg.cc  
		  cpp/WrapFlgfwr.cc 
		  )

# Create the executable that can call flukam_
add_executable(mainfludag cpp/mainFluDAG.cpp ${commonfiles} ${wrapperfiles})
set (CMAKE_Fortran_Compiler gfortran)
ENABLE_LANGUAGE(Fortran)
target_link_libraries(mainfludag ${moablibs} ${FLUPRO}/libflukahp.a gfortran)
