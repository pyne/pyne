# CMAKE build script for MeshTally testing

PROJECT(MeshTallyUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(GTest REQUIRED)
# Locally defined vars
# set(HOME /filespace/people/z/zachman)
IF (NOT DEFINED HOME) 
	SET (HOME "$ENV{HOME}")
	message("Home is ${HOME}")
ENDIF (NOT DEFINED HOME)

SET(MOAB_HOME ${HOME}/dagmc_bld/MOAB)

# Where to look for includes
INCLUDE_DIRECTORIES("${HOME}/dagmc_bld/MOAB/include")

ADD_LIBRARY(tally 
	${DAGMC_MESHTALLY_SOURCE}/TallyEvent.cpp
	${DAGMC_MESHTALLY_SOURCE}/Tally.cpp
)

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
