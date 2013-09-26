# CMake build script for FluDAG testing

project (FluDAGUnitTests)
cmake_minimum_required(VERSION 2.8)

find_package(GTest REQUIRED)

# If FLUKA is installed, $FLUPRO is probably an environment varaiable
IF (NOT DEFINED FLUPRO)
	set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

# IF (NOT DEFINED MOAB_HOME OR NOT DEFINED TOPDIR)
#	INCLUDE(local.cmake)
#ENDIF (NOT DEFINED MOAB_HOME OR NOT DEFINED TOPDIR)

# IF (NOT DEFINED TOPDIR)
#	message("Please define the directory containing DAGMC/FluDAG.")
#	message("Definition can be done with the -D switch or in local.cmake.")
# ENDIF (NOT DEFINED TOPDIR)


IF (NOT DEFINED FLUDAG)
	set (FLUDAG  ${TOPDIR}/DAGMC/FluDAG)
	set (FLUDAG_SRC ${FLUDAG}/src)
	set (FLUDAG_TST	${FLUDAG}/testing)
ENDIF (NOT DEFINED FLUDAG)

set(MOAB_INCLUDE ${MOAB_HOME}/include)
set(MOAB_LIB_DIR ${MOAB_HOME}/lib)
set (moablibs ${MOAB_LIB_DIR}/libMOAB.so ${MOAB_LIB_DIR}/libdagmc.so)

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${FLUDAG}/src/cpp)
include_directories(${MOAB_INCLUDE})

set (gmain      ${FLUDAG_TST}/iface/gtest_main.cc)
set (testfile ${FLUDAG_TST}/iface/tests/issue22.cpp)
set (issue22_src ${FLUDAG_SRC}/cpp/Input.cpp ${FLUDAG_SRC}/cpp/fluka_funcs.cpp ${FLUDAG_SRC}/cpp/dagmc_utils.cpp)
add_executable(issuetest ${gmain} ${testfile} ${issue22_src} ${gmain})
target_link_libraries(issuetest ${PROJECT_NAME_STR} ${moablibs} ${GTEST_BOTH_LIBRARIES} pthread)

# Test call to nrml with concentric spheres geometry
set (wrappers ${FLUDAG_SRC}/cpp/fluka_funcs.cpp ${FLUDAG_SRC}/cpp/dagmc_utils.cpp ${gmain})
set (nrml_spheres_src ${FLUDAG_TST}/iface/tests/nrml_spheres_test.cpp)
add_executable(nrml_spheres ${nrml_spheres_src} ${wrappers})
target_link_libraries(nrml_spheres ${PROJECT_NAME_STR} ${moablibs} ${GTEST_BOTH_LIBRARIES} pthread)

# Test call to nrml with nested cube geometry
set (nrml_cubes_src ${FLUDAG_TST}/iface/tests/nrml_cubes_test.cpp)
add_executable(nrml_cubes ${nrml_cubes_src} ${wrappers})
target_link_libraries(nrml_cubes ${PROJECT_NAME_STR} ${moablibs} ${GTEST_BOTH_LIBRARIES} pthread)

# Test looking for material strings in a list of material strings
set(findMat_src ${FLUDAG_TST}/iface/tests/findMat.cpp)
# add_executable(findMat ${findMat_src} ${wrappers} ${gmain})
add_executable(findMat ${findMat_src} ${gmain})
target_link_libraries(findMat ${PROJECT_NAME_STR} ${GTEST_BOTH_LIBRARIES} pthread)

set (all_test_src ${nrml_cubes_src} ${nrml_spheres_src})
add_executable(all_tests ${all_test_src} ${wrappers})
target_link_libraries(all_tests ${PROJECT_NAME_STR} ${moablibs} ${GTEST_BOTH_LIBRARIES} pthread)

ENABLE_TESTING()


#add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})

