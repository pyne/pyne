PROJECT(FluDAGUnitTests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

# Note that CMAKE_SOURCE_DIR (without 'CURRNET') is the top-level dir containing CMakeLists.txt
SET (gtest_home    "${CMAKE_SOURCE_DIR}/../gtest")
SET (gtest_dir     "${gtest_home}/gtest-1.7.0")
SET (fludag_source "${CMAKE_SOURCE_DIR}/src")

# Where to look for includes
INCLUDE_DIRECTORIES("${MOAB_HOME}/include" "${gtest_dir}" "${gtest_dir}/include" ${fludag_source})

# Setup for gtest libraries
SET(cxx_base_flags "-Wall -Wshadow")
SET(cxx_exception_flags "-fexceptions")
# Ref: https://groups.google.com/forum/#!topic/googletestframework/8S3QSezcoRA
SET(cxx_strict_flags "-Wextra -Wno-missing-field-initializers")

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

cxx_static_library(gtest ${cxx_strict} "${gtest_dir}/src/gtest-all.cc")

SET (tgt "gtest")
ADD_CUSTOM_COMMAND (TARGET ${tgt}
	COMMAND ln "-s"
		"${CMAKE_CURRENT_SOURCE_DIR}/slabs.h5m" 
		"${CMAKE_CURRENT_BINARY_DIR}/slabs.h5m"
		COMMENT "...linking test file in ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_CURRENT_BINARY_DIR}/slabs.h5m"
	DEPENDS
	VERBATIM
 	)

# make fludag source files into a library
SET(fludag_source_files 
	    ${fludag_source}/fluka_funcs.cpp
 )

SET (fludag_test_source ${fludag_source_files}
	${CMAKE_CURRENT_SOURCE_DIR}/test_FlukaFuncs.cpp
	)

# set libraries needed to link to test cases
SET(LIBS
    pthread
    ${MOAB_HOME}/lib/libMOAB.so
    ${MOAB_HOME}/lib/libdagmc.so
    ${MOAB_HOME}/lib/libMOAB.so
)

ADD_EXECUTABLE (fludag_unit_tests 
	fludag_unit_test_driver.cc
	${fludag_test_source}
	)

TARGET_LINK_LIBRARIES (fludag_unit_tests ${LIBS} gtest)
	
# enable test cases
enable_testing()

# list test cases for compatibility with BaTLab
add_test(test_FlukaFuncs test_FlukaFuncs)
