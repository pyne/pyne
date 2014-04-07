# CMAKE build script for FluDag testing

PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

SET (gtest_home "${CMAKE_CURRENT_SOURCE_DIR}/../../../gtest")
SET (gtest_dir "${gtest_home}/gtest-1.7.0")
SET (fludag_source "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${gtest_dir}/include")


# make fludag source files into a library
ADD_LIBRARY(fludag SHARED
	    ${fludag_source}/fluka_funcs.cpp
	    ${fludag_source}/dagmc_utils.cpp
	    ${fludag_source}/Input.cpp
	    ${fludag_source}/WrapMag.cc
	    ${fludag_source}/WrapLookFX.cc
	    ${fludag_source}/WrapReg.cc
	    ${fludag_source}/WrapFlgfwr.cc
	    ${fludag_source}/UnitNumberManager.cpp
 )


# set libraries needed to link to test cases
SET(LIBRARIES
    ${gtest_home}/lib/libgtest.a
    ${gtest_home}/lib/libgtest_main.a
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
