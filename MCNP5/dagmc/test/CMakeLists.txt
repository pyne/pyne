# CMAKE build script for MeshTally testing

PROJECT(MeshTallyUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# FIND_PACKAGE(GTest REQUIRED PATHS ${GTEST_DIR})

ADD_SUBDIRECTORY(gtest/gtest-1.7.0)

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" gtest/gtest-1.7.0/include)

# Locally defined vars
# set(HOME /filespace/people/z/zachman)
IF (NOT DEFINED HOME) 
	SET (HOME "$ENV{HOME}")
	message("Home is ${HOME}")
ENDIF (NOT DEFINED HOME)

SET(MOAB_HOME ${HOME}/opt/moab)

# Where to look for includes
INCLUDE_DIRECTORIES("${HOME}/opt/moab/include")

# make mesh tally source files into a library
ADD_LIBRARY(meshtally
    ${DAGMC_MESHTALLY_SOURCE}/Tally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/TallyData.cpp
    ${DAGMC_MESHTALLY_SOURCE}/TallyManager.cpp
    ${DAGMC_MESHTALLY_SOURCE}/MeshTally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/TrackLengthMeshTally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/KDEMeshTally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/KDEKernel.cpp
    ${DAGMC_MESHTALLY_SOURCE}/KDEMeshTally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/KDENeighborhood.cpp
    ${DAGMC_MESHTALLY_SOURCE}/MeshTally.cpp
    ${DAGMC_MESHTALLY_SOURCE}/PolynomialKernel.cpp
    ${DAGMC_MESHTALLY_SOURCE}/Quadrature.cpp
    ${DAGMC_MESHTALLY_SOURCE}/TallyEvent.cpp
)

# set libraries
SET(LIBRARIES
	    gtest
	gtest_main
	 pthread
    meshtally
    ${MOAB_HOME}/lib/libMOAB.so
)

# set up MeshTally test cases
ADD_EXECUTABLE(test_KDEKernel test_KDEKernel.cpp)
TARGET_LINK_LIBRARIES(test_KDEKernel ${LIBRARIES})

ADD_EXECUTABLE(test_KDEMeshTally test_KDEMeshTally.cpp)
TARGET_LINK_LIBRARIES(test_KDEMeshTally ${LIBRARIES})

ADD_EXECUTABLE(test_KDENeighborhood test_KDENeighborhood.cpp)
TARGET_LINK_LIBRARIES(test_KDENeighborhood ${LIBRARIES})

ADD_EXECUTABLE(test_PolynomialKernel test_PolynomialKernel.cpp)
TARGET_LINK_LIBRARIES(test_PolynomialKernel ${LIBRARIES})

ADD_EXECUTABLE(test_Quadrature test_Quadrature.cpp)
TARGET_LINK_LIBRARIES(test_Quadrature ${LIBRARIES})

# enable MeshTally test cases
ENABLE_TESTING()

# list MeshTally test cases
ADD_TEST(test_KDEKernel test_KDEKernel)
ADD_TEST(test_KDEMeshTally test_KDEMeshTally)
ADD_TEST(test_KDENeighborhood test_KDENeighborhood)
ADD_TEST(test_PolynomialKernel test_PolynomialKernel)
ADD_TEST(test_Quadrature test_Quadrature)
