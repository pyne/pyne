cmake_minimum_required(VERSION 2.8)
project (MyFluDAG)

IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 
IF (NOT DEFINED FLUPRO)
	set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIB_DIR ${MOAB_HOME}/lib)
set (moablibs ${MOAB_LIB_DIR}/libdagmc.so ${MOAB_LIB_DIR}/libMOAB.so)


# Where to look for includes
include_directories(include ${MOAB_INCLUDE})

set (commonfiles cpp/fluka_funcs.cpp cpp/dagmc_utils.cpp cpp/Input.cpp)
set (wrapperfiles cpp/WrapReg2Name.cc cpp/WrapLookMG.cc cpp/WrapMag.cc 
		  cpp/WrapLookFX.cc cpp/WrapReg.cc  
		  cpp/WrapFlgfwr.cc cpp/WrapDN.cc
		  )

# Create the executable that can call flukam_
add_executable(mainfludag cpp/mainFluDAG.cpp ${commonfiles} ${wrapperfiles})
set (CMAKE_Fortran_Compiler gfortran)
ENABLE_LANGUAGE(Fortran)
target_link_libraries(mainfludag ${moablibs} ${FLUPRO}/libflukahp.a gfortran)
