# CMAKE build script for FluDag testing

PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# FIND_PACKAGE(GTest REQUIRED PATHS ${GTEST_DIR})
# ADD_SUBDIRECTORY(gtest/gtest-1.7.0)

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${GTEST_HOME}/gtest-1.7.0/include")

# make mesh tally source files into a library
ADD_LIBRARY(fludag
	    ${DAGMC_FLUDAG_SOURCE}/fluka_funcs.cpp
	    ${DAGMC_FLUDAG_SOURCE}/dagmc_utils.cpp
	    ${DAGMC_FLUDAG_SOURCE}/Input.cpp
	    ${DAGMC_FLUDAG_SOURCE}/UnitNumberManager.cpp
	    ${DAGMC_FLUDAG_SOURCE}/WrapReg2Name.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapMag.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapLookFX.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapReg.cc
	    ${DAGMC_FLUDAG_SOURCE}/WrapFlgfwr.cc
 )


# set libraries
SET(LIBRARIES
    ${GTEST_HOME}/build/libgtest.a
    ${GTEST_HOME}/build/libgtest_main.a
    pthread
    fludag
    ${MOAB_HOME}/lib/libMOAB.so
)

# set up FluDAG test cases
ADD_EXECUTABLE(test_UnitNumber test_UnitNumber.cpp)
TARGET_LINK_LIBRARIES(test_UnitNumber ${LIBRARIES})


# enable FluDAG test cases
ENABLE_TESTING()

# list FluDAG test cases
ADD_TEST(test_UnitNumber test_UnitNumber)
