PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

SET (gtest_home    "${CMAKE_CURRENT_SOURCE_DIR}/../../gtest")
SET (gtest_dir     "${gtest_home}/gtest-1.7.0")
SET (fludag_source "${CMAKE_CURRENT_SOURCE_DIR}/../src")

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${gtest_dir}" "${gtest_dir}/include" ${fludag_source})

# Setup for gtest libraries
SET(cxx_base_flags "-Wall -Wshadow")
SET(cxx_exception_flags "-fexceptions")
SET(cxx_strict_flags "-Wextra")

IF (CMAKE_USE_PTHREADS_INIT)  # The pthreads library is available.
    SET(cxx_base_flags "${cxx_base_flags} -DGTEST_HAS_PTHREAD=1")
ENDIF()

# Make the gtest libraries
SET(cxx_exception "${CMAKE_CXX_FLAGS} ${cxx_base_flags} ${cxx_exception_flags}")
SET(cxx_strict "${cxx_exception} ${cxx_strict_flags}")
########################################################################
#
# Defines the gtest & gtest_main libraries.  User tests should link
# with one of them.
FUNCTION(cxx_library_with_type name type cxx_flags)
# type can be either STATIC or SHARED to denote a static or shared library.
# ARGN refers to additional arguments after 'cxx_flags'.
      ADD_LIBRARY(${name} ${type} ${ARGN})
      SET_TARGET_PROPERTIES(${name}
	    PROPERTIES
            COMPILE_FLAGS "${cxx_flags}")
    IF (CMAKE_USE_PTHREADS_INIT)
         TARGET_LINK_LIBRARIES(${name} ${CMAKE_THREAD_LIBS_INIT})
    ENDIF()
ENDFUNCTION()

FUNCTION(cxx_static_library name cxx_flags)
	cxx_library_with_type(${name} STATIC "${cxx_flags}" ${ARGN})
ENDFUNCTION()

cxx_static_library(gtest      ${cxx_strict} "${gtest_dir}/src/gtest-all.cc")
cxx_static_library(gtest_main ${cxx_strict} "${gtest_dir}/src/gtest_main.cc")
TARGET_LINK_LIBRARIES(gtest_main gtest)
SET (tgt "gtest")
ADD_CUSTOM_COMMAND (TARGET ${tgt}
	COMMAND cp 
		"${CMAKE_SOURCE_DIR}/slabs.h5m" 
		"${CMAKE_BINARY_DIR}"
	COMMENT "...copying test file from test source ${CMAKE_SOURCE_DIR} TO ${CMAKE_BINARY_DIR}"
	DEPENDS
	VERBATIM
 	)

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
    ${CMAKE_BINARY_DIR}/tests/libgtest.a
    ${CMAKE_BINARY_DIR}/tests/libgtest_main.a
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
