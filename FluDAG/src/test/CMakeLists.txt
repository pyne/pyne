# CMAKE build script for FluDag testing

PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

SET (gtest_home    "${CMAKE_CURRENT_SOURCE_DIR}/../../../gtest")
SET (gtest_dir     "${gtest_home}/gtest-1.7.0")
SET (fludag_source "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${gtest_dir}" "${gtest_dir}/include")

# This makes libgtest.a in the binary directory for an out-of-source gtest
# SET ( GTest "${gtest_dir}/src/gtest-all.cc")
# ADD_LIBRARY (gtest ${GTest})

# make fludag source files into a library
set(FLUDAG_SOURCE ${fludag_source})
ADD_LIBRARY(fludag SHARED
	    ${FLUDAG_SOURCE}/fluka_funcs.cpp
	    ${FLUDAG_SOURCE}/dagmc_utils.cpp
	    ${FLUDAG_SOURCE}/Input.cpp
	    ${FLUDAG_SOURCE}/WrapMag.cc
	    ${FLUDAG_SOURCE}/WrapLookFX.cc
	    ${FLUDAG_SOURCE}/WrapReg.cc
	    ${FLUDAG_SOURCE}/WrapFlgfwr.cc
	    ${FLUDAG_SOURCE}/UnitNumberManager.cpp
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

# To add a new file, add it to this lis.  
# SET (FLUDAG_TEST_SOURCE
#	"${CMAKE_CURRENT_SOURCE_DIR}/test_FlukaFuncs.cpp"
#	"${CMAKE_CURRENT_SOURCE_DIR}/test_UnitNubmer.cpp"
#	)

# ADD_EXECUTABLE (fludag_unit_tests fludag_unit_test_driver.cc ${FLUDAG_TEST_SOURCE})

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
