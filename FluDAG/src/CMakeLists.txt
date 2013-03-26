cmake_minimum_required(VERSION 2.8)
project (MyFluDAG)
set(HOME /filespace/people/z/zachman)
set(FLUPRO ${HOME}/flugg_bld/FLUKA)
set(MOAB_HOME ${HOME}/dagmc_bld/MOAB)

set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIB_DIR ${MOAB_HOME}/lib)
set (moablibs ${MOAB_LIB_DIR}/libMOAB.so ${MOAB_LIB_DIR}/libdagmc.so)


# Where to look for includes
include_directories(include ${MOAB_INCLUDE})

set (commonfiles cpp/readVol.cpp cpp/fluka_funcs.cpp cpp/dagmc_utils.cpp)

# add_executable(readvol cpp/readVol.cpp cpp/fluka_funcs.cpp)
add_executable(readvol cpp/mainReadVol.cpp ${commonfiles})
target_link_libraries(readvol ${moablibs})

# This is not used yet, but may come into play when cmake can be used to link the fluka library
set (wrapperfiles cpp/WrapInit.cc cpp/WrapLookDB.cc cpp/WrapSavHist.cc cpp/WrapLookZ.cc 
                  cpp/WrapDN.cc cpp/WrapReg2Name.cc cpp/WrapLookMG.cc cpp/WrapMag.cc 
		  cpp/WrapLookFX.cc cpp/WrapReg.cc cpp/WrapIniHist.cc cpp/WrapIncrHist.cc 
		  cpp/WrapFlgfwr.cc cpp/WrapG1RT.cc)

# Create the executable that can call flukam_
add_executable(mainfludag cpp/mainFluDAG.cpp ${commonfiles} ${wrapperfiles})
set (CMAKE_Fortran_Compiler gfortran)
ENABLE_LANGUAGE(Fortran)
target_link_libraries(mainfludag ${moablibs} ${FLUPRO}/libflukahp.a /usr/lib/libgfortran.so.3)
