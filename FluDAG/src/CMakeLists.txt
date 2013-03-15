cmake_minimum_required(VERSION 2.8)
project (MyFluDAG)
set(HOME /filespace/people/z/zachman)
set(MOAB_HOME ${HOME}/dagmc_bld/MOAB)

set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIBS ${MOAB_HOME}/lib)

# Where to look for includes
include_directories(include ${MOAB_INCLUDE})

# add_executable(readvol cpp/readVol.cpp cpp/fluka_funcs.cpp)
add_executable(readvol cpp/readVol.cpp cpp/fluka_funcs.cpp cpp/dagmc_utils.cpp)
target_link_libraries(readvol ${MOAB_LIBS}/libMOAB.so ${MOAB_LIBS}/libdagmc.so)

# This is not used yet, but may come into play when cmake can be used to link the fluka library
set (wrapperfiles cpp/WrapInit.cc cpp/WrapLookDB.cc cpp/WrapSavHist.cc cpp/WrapLookZ.cc cpp/WrapG1.cc 
                  cpp/WrapDN.cc cpp/WrapReg2Name.cc cpp/WrapNorml.cc cpp/WrapLookMG.cc cpp/WrapMag.cc 
		  cpp/WrapLookFX.cc cpp/WrapReg.cc cpp/WrapIniHist.cc cpp/WrapIncrHist.cc 
		  cpp/WrapFlgfwr.cc cpp/WrapG1RT.cc)
