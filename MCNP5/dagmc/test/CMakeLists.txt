# CMAKE build script for MeshTally testing

PROJECT(MeshTallyUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(GTest REQUIRED)

# make mesh tally source files into a library
ADD_LIBRARY(meshtally SHARED
    ${DAGMC_MESHTALLY_SOURCE}/KDEKernel.cpp
    ${DAGMC_MESHTALLY_SOURCE}/PolynomialKernel.cpp
    ${DAGMC_MESHTALLY_SOURCE}/Quadrature.cpp
)

# set libraries
SET(LIBRARIES
    gtest
    gtest_main
    pthread
    meshtally
)

# set up MeshTally test cases
ADD_EXECUTABLE(test_KDEKernel test_KDEKernel.cpp)
TARGET_LINK_LIBRARIES(test_KDEKernel ${LIBRARIES})

ADD_EXECUTABLE(test_PolynomialKernel test_PolynomialKernel.cpp)
TARGET_LINK_LIBRARIES(test_PolynomialKernel ${LIBRARIES})

ADD_EXECUTABLE(test_Quadrature test_Quadrature.cpp)
TARGET_LINK_LIBRARIES(test_Quadrature ${LIBRARIES})

# enable MeshTally test cases
ENABLE_TESTING()

# list MeshTally test cases
ADD_TEST(test_KDEKernel test_KDEKernel)
ADD_TEST(test_PolynomialKernel test_PolynomialKernel)
ADD_TEST(test_Quadrature test_Quadrature)
