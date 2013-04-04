cmake_minimum_required(VERSION 2.8)
project (TestMyFluDAG)
find_package(GTest REQUIRED)
IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

IF (NOT DEFINED FLUDAG)
	set (FLUDAG  ${HOME}/DAGMC/FluDAG)
	set (FLUDAG_SRC ${FLUDAG}/src)
	set (FLUDAG_TST	${FLUDAG}/testing)
ENDIF (NOT DEFINED FLUDAG)

set (FLUPRO "$ENV{FLUPRO}")



# Where to look for includes
#include_directories(../src/include ${MOAB_INCLUDE} ${GTEST_DIR})
IF(NOT DEFINED GTEST_DIR)
	set(GTEST_DIR "$ENV{GTEST_DIR}")
ENDIF(NOT DEFINED GTEST_DIR)

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${HOME}/DAGMC/FluDAG/src/include)

set (ifacefiles ${FLUDAG_TST}/iface/issue22.cpp)
add_executable(issuetest ${ifacefiles} ${FLUDAG_TST}/iface/gtest_main.cc ${FLUDAG_SRC}/cpp/Input.cpp)
target_link_libraries(issuetest ${PROJECT_NAME_STR} ${GTEST_BOTH_LIBRARIES} pthread)
ENABLE_TESTING()
#add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})

