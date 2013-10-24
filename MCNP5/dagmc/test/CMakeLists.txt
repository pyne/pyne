# CMAKE build script for DAGMCTally testing

PROJECT(DAGMCTallyUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# FIND_PACKAGE(GTest REQUIRED PATHS ${GTEST_DIR})
ADD_SUBDIRECTORY(gtest/gtest-1.7.0)

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" gtest/gtest-1.7.0/include)

# make DAGMC Tally source files into a library
ADD_LIBRARY(dagtally SHARED
    ${DAGMC_TALLY_SOURCE}/Tally.cpp
    ${DAGMC_TALLY_SOURCE}/TallyData.cpp
    ${DAGMC_TALLY_SOURCE}/TallyManager.cpp
    ${DAGMC_TALLY_SOURCE}/MeshTally.cpp
    ${DAGMC_TALLY_SOURCE}/TrackLengthMeshTally.cpp
    ${DAGMC_TALLY_SOURCE}/KDEMeshTally.cpp
    ${DAGMC_TALLY_SOURCE}/CellTally.cpp
    ${DAGMC_TALLY_SOURCE}/KDEKernel.cpp
    ${DAGMC_TALLY_SOURCE}/KDENeighborhood.cpp
    ${DAGMC_TALLY_SOURCE}/PolynomialKernel.cpp
    ${DAGMC_TALLY_SOURCE}/Quadrature.cpp
)

# set libraries
SET(LIBRARIES
    gtest
    gtest_main
    pthread
    dagtally
    ${MOAB_HOME}/lib/libMOAB.so
)

# set up DAGMC Tally test cases
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

# enable DAGMC Tally test cases
ENABLE_TESTING()

# list DAGMC Tally test cases
ADD_TEST(test_KDEKernel test_KDEKernel)
ADD_TEST(test_KDEMeshTally test_KDEMeshTally)
ADD_TEST(test_KDENeighborhood test_KDENeighborhood)
ADD_TEST(test_PolynomialKernel test_PolynomialKernel)
ADD_TEST(test_Quadrature test_Quadrature)
