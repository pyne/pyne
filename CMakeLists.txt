# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Setup the project
project(pyne)

option ( ENABLE_SPATIAL_SOLVERS     "Should build AHOT spatial solvers?"                 ON )
option ( ENABLE_ENSDF_PROCESSING   "Should build ENSDF processing tools?"               ON )

include(PyneMacros)
pyne_print_logo()  # Beware of dragons
pyne_set_platform()
pyne_set_asm_platform()
pyne_set_build_spatial_solver()
pyne_set_build_type()
pyne_set_fast_compile()
pyne_configure_rpath()

# enable assembly
enable_language(ASM)

# This makes all the libraries build as SHARED
set(BUILD_SHARED_LIBS true)
message(STATUS "CMake Install Prefix: ${CMAKE_INSTALL_PREFIX}")


#
# Allow for dependnecies to exist in non-install location
#
if (NOT DEPS_ROOT_DIR)
  get_filename_component(compdir ${CMAKE_CXX_COMPILER} DIRECTORY)
  get_filename_component(DEPS_ROOT_DIR ${compdir} DIRECTORY)
endif (NOT DEPS_ROOT_DIR)
set(DEPS_HINTS HINTS "${DEPS_ROOT_DIR}")
set(DEPS_BIN_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/bin")
set(DEPS_LIB_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/lib")
set(DEPS_INCLUDE_HINTS HINTS "${DEPS_ROOT_DIR}/include")
message(STATUS "Dependency Root Dir: ${DEPS_ROOT_DIR}")
message(STATUS "Dependency Hints: ${DEPS_HINTS}")
message(STATUS "Dependency Binary Hints: ${DEPS_BIN_HINTS}")
message(STATUS "Dependency Library Hints: ${DEPS_LIB_HINTS}")
message(STATUS "Dependency Include Hints: ${DEPS_INCLUDE_HINTS}")

# Find HDF5
include(ExternalProject)

if(HDF5_ROOT)
  find_package(HDF5 REQUIRED)
else()
find_package(HDF5)
endif()

if(NOT HDF5_FOUND AND NOT HDF5_ROOT)
  MESSAGE(STATUS "HDF5 not found. Attempting to download and install.")
  set(HDF5_ROOT "${CMAKE_BINARY_DIR}/hdf5")
  ExternalProject_Add(hdf5
    PREFIX ${HDF5_ROOT}
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.11/src/hdf5-1.10.11.tar.gz
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    DOWNLOAD_EXTRACT_TIMESTAMP true
    BUILD_BYPRODUCTS "${HDF5_ROOT}/lib/libhdf5.a"
  )
  set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
  set(HDF5_LIBRARY_DIRS "${HDF5_ROOT}/lib")
  set(HDF5_LIBRARIES "${HDF5_ROOT}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(HDF5_C_LIBRARIES "${HDF5_ROOT}/lib/libhdf5.a")
  set(HDF5_H_LIBRARIES "${HDF5_ROOT}/lib/libhdf5_hl.a")
endif()

include_directories("${HDF5_INCLUDE_DIRS}")
if(NOT HDF5_LIBRARY_DIRS STREQUAL "")
  link_directories(${HDF5_LIBRARY_DIRS})
endif()
set(LIBS_HDF5 ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 library dirs: ${HDF5_LIBRARY_DIRS}")
message(STATUS "HDF5 Libraries: ${HDF5_C_LIBRARIES}")

# Set up shared and static versions of the HDF5 libraries
set(HDF5_LIBRARIES_SHARED ${HDF5_LIBRARIES})
if(BUILD_STATIC_LIBS)
  string(REPLACE ${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_STATIC_LIBRARY_SUFFIX}
          HDF5_LIBRARIES_STATIC "${HDF5_LIBRARIES_SHARED}"
  )
endif()
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_LIBRARIES_SHARED)
endif()

# Print out some more information about the HDF5 setup
MESSAGE("-- HDF5 libraries shared: ${HDF5_LIBRARIES_SHARED}")
MESSAGE("-- HDF5 libraries static: ${HDF5_LIBRARIES_STATIC}")

# Find MOAB
MESSAGE("-- MOAB root: ${MOAB_ROOT}")
find_package(MOAB)

if(${MOAB_FOUND})
  include_directories(${MOAB_INCLUDE_DIRS})
  get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
  link_directories(${MOAB_LIBRARY_DIRS})
endif(${MOAB_FOUND})


# Look for DAGMC if requested
if(${MOAB_FOUND})
  # Always look for DAGMC library
  if(NOT DAGMC_ROOT)
    set(DAGMC_ROOT "${DEPS_ROOT_DIR}")
  endif()
  MESSAGE("-- DAGMC Root: ${DAGMC_ROOT}")
  find_package(DAGMC)
  message("-- DAGMC Found: ${DAGMC_FOUND}")
  if(${DAGMC_FOUND})
    include_directories(${DAGMC_INCLUDE_DIRS})
    link_directories(${DAGMC_LIBRARY_DIRS})
  endif(${DAGMC_FOUND})
else(${MOAB_FOUND})
  message("-- DAGMC: Not searching because MOAB was not found.")
endif(${MOAB_FOUND})


# Find Python
find_package(
  Python 
  COMPONENTS Interpreter Development.Module NumPy 
  REQUIRED
)

pyne_setup_fortran()

# Find f2py, if building spatial solver
if (BUILD_SPATIAL_SOLVER)
  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_package(LAPACK)
  if(NOT BLAS_FOUND)
    MESSAGE(STATUS "BLAS not found. Attempting to download and install.")
    set(BLAS_ROOT "${CMAKE_BINARY_DIR}/blas")
    ExternalProject_Add(blas
      PREFIX ${BLAS_ROOT}
      URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.26/OpenBLAS-0.3.26.tar.gz
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      DOWNLOAD_EXTRACT_TIMESTAMP true
      BUILD_BYPRODUCTS "${BLAS_ROOT}/lib/libopenblas.a"
    )
    set(BLAS_INCLUDE_DIRS "${BLAS_ROOT}/include")
    set(BLAS_LIBRARY_DIRS "${BLAS_ROOT}/lib")
    set(BLAS_LIBRARIES "${BLAS_ROOT}/lib/libopenblas.a")
  endif()
  if(NOT LAPACK_FOUND)
    MESSAGE(STATUS "LAPACK not found. Attempting to download and install.")
    set(LAPACK_ROOT "${CMAKE_BINARY_DIR}/lapack")
    ExternalProject_Add(lapack
      PREFIX ${LAPACK_ROOT}
      URL https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.0.tar.gz
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      DOWNLOAD_EXTRACT_TIMESTAMP true
      BUILD_BYPRODUCTS "${LAPACK_ROOT}/lib/liblapack.a"
    )
    set(LAPACK_INCLUDE_DIRS "${LAPACK_ROOT}/include")
    set(LAPACK_LIBRARY_DIRS "${LAPACK_ROOT}/lib")
    set(LAPACK_LIBRARIES "${LAPACK_ROOT}/lib/liblapack.a")
  endif()
endif()


# Find Cython
find_program(CYTHON "cython")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/pyne")

# Add JsonCpp Flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_IS_AMALGAMATION")

#
# Build and install!
#
message(STATUS "LIBS: ${LIBS}")

# add sub dirs
add_subdirectory(src)
add_subdirectory(pyne)
add_subdirectory(pyne/xs)

# Collect scripts
file(GLOB SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/*")

# Filter scripts based on the OS
if(WIN32)
    list(FILTER SCRIPTS INCLUDE REGEX "\\.bat$")
else()
    list(FILTER SCRIPTS EXCLUDE REGEX "\\.bat$")
endif()

# Install scripts to SKBUILD_SCRIPTS_DIR
if (SKBUILD_SCRIPTS_DIR)
  install(FILES ${SCRIPTS} DESTINATION ${SKBUILD_SCRIPTS_DIR})
endif()
