# Defines the CMake commands/policies
cmake_minimum_required( VERSION 2.8.5 )

# Set the project name
project( pyne Fortran CXX)


#determine if spatial solver module should be built
IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  IF(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6" AND NOT APPLE)
    SET(BUILD_SPATIAL_SOLVER true)
  ELSE()
    SET(BUILD_SPATIAL_SOLVER false)
  ENDIF()
ENDIF()

IF(BUILD_SPATIAL_SOLVER)
# languages
enable_language(Fortran)

# Augment the Fortran implicit link libraries
message(STATUS
  "CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES = ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
if (APPLE)
  message(STATUS
    "Before fix: CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES = ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  set(LIBGCC_S)
# The previous method found the gcc_s library by version,
    # find_library(LIBGCC_S_PATH gcc_s.${gcc_s_ver}
# but this always found the library regardless of version, and then
# the name as set was different from what was found.  This new way
# ensures that the name is consistent with what was found.
  set(LIBGCC_S_PATH)
  foreach (gcc_s_ver 3 2 1)
    find_library(LIBGCC_S_PATH gcc_s.${gcc_s_ver}
      PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
      NODEFAULTPATH
    )
    if (LIBGCC_S_PATH)
      break()
    endif ()
  endforeach ()
  message(STATUS "LIBGCC_S_PATH = ${LIBGCC_S_PATH}.")
  if (LIBGCC_S_PATH)
    get_filename_component(LIBGCC_S ${LIBGCC_S_PATH} NAME)
# Pull off leading lib and trailing .dylib.
    string(REGEX REPLACE "^lib" "" LIBGCC_S ${LIBGCC_S})
    string(REGEX REPLACE "\\.dylib$" "" LIBGCC_S ${LIBGCC_S})
  endif ()
  message(STATUS "LIBGCC_S = ${LIBGCC_S}.")
  if (LIBGCC_S)
    set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${LIBGCC_S})
  endif ()
endif ()
message(STATUS
  "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES = ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
ENDIF(BUILD_SPATIAL_SOLVER)

# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Default to release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "The build type" FORCE)
endif(NOT CMAKE_BUILD_TYPE)


# quiets fortify_source warnings when not compiling with optimizations
# in linux distros where compilers were compiled with fortify_source enabled by
# default (e.g. Arch linux).
MESSAGE("-- Build type: ${CMAKE_BUILD_TYPE}")
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
IF(NOT ${BUILD_TYPE} STREQUAL "release")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
ENDIF()

# Include the HDF5 library and c++ headers
find_package( HDF5 COMPONENTS C REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
if(WIN32)
    # FindHDF5 finds the includes but not the libraries on Windows (MSYS).  Annoying!
    get_filename_component(_hdf5libdir ${HDF5_INCLUDE_DIRS} PATH)
    list(APPEND HDF5_LIBRARY_DIRS "${_hdf5libdir}/bin")
	list(APPEND HDF5_LIBRARY_DIRS "${_hdf5libdir}/lib")
elseif(NOT HDF5_LIBRARY_DIRS)
  list(GET HDF5_LIBRARIES 0 _hdf5lib)
  if (${_hdf5lib} STREQUAL debug)
    list(GET HDF5_LIBRARIES 1 _hdf5lib)
  endif ()
  get_filename_component(HDF5_LIBRARY_DIRS ${_hdf5lib} PATH)
endif(WIN32)
link_directories(${HDF5_LIBRARY_DIRS})
add_definitions(${HDF5_DEFINITIONS})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
message("-- HDF5_LIBRARY_DIRS = ${HDF5_LIBRARY_DIRS}")

# Use new Python library finder
find_package(PythonInterp REQUIRED)
find_package(PythonLibsNew REQUIRED)
message("-- PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")

# Check for MOAB library
find_package(MOAB)
message("-- MOAB Found: ${MOAB_FOUND}")
if(${MOAB_FOUND})
    include_directories(${MOAB_INCLUDE_DIR})
    get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
    link_directories(${MOAB_LIBRARY_DIRS})
endif(${MOAB_FOUND})

# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of cmake/UseCython.cmake.
include( UseCython )

execute_process(COMMAND "${CYTHON_EXECUTABLE}" "-V"
    ERROR_VARIABLE CYTHON_VERSION
    ERROR_STRIP_TRAILING_WHITESPACE)
message( STATUS "Cython Version: " ${CYTHON_VERSION} )


# Also, there are some custom pyne macros to make for less verbosity
include( PyneMacros )

# Beware of dragons
print_logo()

# This makes all the libraries build as SHARED
set(BUILD_SHARED_LIBS true)
message("-- CMake Install Prefix: ${CMAKE_INSTALL_PREFIX}")

# RPATH Settings
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
    if(HAVE_SETUPTOOLS)
      set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
    else(HAVE_SETUPTOOLS)
      set(CMAKE_INSTALL_RPATH "@loader_path/../../..")
    endif(HAVE_SETUPTOOLS)
else(APPLE)
    # For linux
    if(HAVE_SETUPTOOLS)
      set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib")
    else(HAVE_SETUPTOOLS)
      set(CMAKE_INSTALL_RPATH "\$ORIGIN/lib")
    endif(HAVE_SETUPTOOLS)
endif(APPLE)
message("-- RPATH: ${CMAKE_INSTALL_RPATH}")

# find numpy and include the numpy headers
find_package(Numpy REQUIRED)
include_directories("${NUMPY_INCLUDE_DIR}")

# Add JsonCpp Flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_IS_AMALGAMATION")

# With CMake, a clean separation can be made between the source tree and the
# build tree.  When all source is compiled, as with pure C/C++, the source is
# no-longer needed in the build tree.  However, with pure *.py source, the
# source is processed directly.  To handle this, we reproduce the availability
# of the source files in the build tree.
add_custom_target(ReplicatePythonSourceTree ALL "${CMAKE_COMMAND}" -P
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(src)
add_subdirectory(pyne)
add_subdirectory(pyne/xs)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

message("-- Copying C/C++ header files.")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/cpp/" DESTINATION
    "${CMAKE_BINARY_DIR}/pyne/include/" FILES_MATCHING PATTERN "*.h")

message("-- Copying Cython header files.")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/pyne/"
    DESTINATION "${CMAKE_BINARY_DIR}/pyne/"
    FILES_MATCHING PATTERN "*.pxd"
                   PATTERN "lib" EXCLUDE
                   PATTERN "include" EXCLUDE)

message("-- Copying scripts.")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/scripts" DESTINATION "${CMAKE_BINARY_DIR}")
