# CMAKE build script for FluDag testing

PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${GTEST_HOME}/gtest-1.7.0/include")

# make fludag source files into a library
ADD_LIBRARY(fludag SHARED
	    ${DAGMC_FLUDAG_SOURCE}/fluka_funcs.cpp
	    ${DAGMC_FLUDAG_SOURCE}/dagmc_utils.cpp
	    ${DAGMC_FLUDAG_SOURCE}/Input.cpp
	    ${DAGMC_FLUDAG_SOURCE}/WrapMag.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapLookFX.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapReg.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapFlgfwr.cc
	    ${DAGMC_FLUDAG_SOURCE}/UnitNumberManager.cpp
 )


# set libraries needed to link to test cases
SET(LIBRARIES
    ${GTEST_HOME}/lib/libgtest.a
    ${GTEST_HOME}/lib/libgtest_main.a
    pthread
    fludag
    ${MOAB_HOME}/lib/libMOAB.so
    ${MOAB_HOME}/lib/libdagmc.so
    ${MOAB_HOME}/lib/libMOAB.so
)

# set up test cases
ADD_EXECUTABLE(test_UnitNumber test_UnitNumber.cpp)
TARGET_LINK_LIBRARIES(test_UnitNumber ${LIBRARIES})

ADD_EXECUTABLE(test_FlukaFuncs test_FlukaFuncs.cpp)
TARGET_LINK_LIBRARIES(test_FlukaFuncs ${LIBRARIES})

# This adds dashboard testing and also calls enable_testing()
# include( CTest )

# enable test cases
enable_testing()

# list test cases
add_test(test_UnitNumber test_UnitNumber)
add_test(test_FlukaFuncs test_FlukaFuncs)
