# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.15)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Setup the project
project(pyne)

# Include the GNUInstallDirs module.
# This module provides installation directories as per the GNU Coding Standards.
include(GNUInstallDirs)

# Include the FeatureSummary module.
# This module handles the recording, reporting and display of enabled/disabled features.
include(FeatureSummary)

# Include the ExternalProject module
# This module handles downloading and building external projects
include(ExternalProject)

# Set build options that can be toggled by the user
option ( ENABLE_MOAB                "Should PyNE use MOAB?"                        ON )
option ( ENABLE_DAGMC               "Should PyNE use DAGMC?"                       ON )
option ( ENABLE_SPATIAL_SOLVERS     "Should build AHOT spatial solvers?"           ON )
option ( ENABLE_ENSDF_PROCESSING    "Should build ENSDF processing tools?"         ON )
option ( PYNE_FAST_COMPILE          "Should PyNE fast compile?"                    ON )

# Add feature info
add_feature_info(SpatialSolvers ENABLE_SPATIAL_SOLVERS "AHOT spatial solvers")
add_feature_info(ENSDFProcessing ENABLE_ENSDF_PROCESSING "ENSDF processing tools")
add_feature_info(PyneFastCompile PYNE_FAST_COMPILE "PyNE fast compile")

# Include PyneMacros file
include(PyneMacros)
print_pyne_logo()  # Beware of dragons
pyne_set_platform()
pyne_set_asm_platform()
pyne_set_build_spatial_solver()
pyne_set_build_type()

# This makes all the libraries build as SHARED
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

############################################################
# Setup HDF5                                               #
############################################################

# If HDF5 is installed in a non-standard location, set HDF5_ROOT
# to the root directory of the HDF5 installation.
# For example: export HDF5_ROOT=/path/to/hdf5
# If HDF5 is not installed, no warry, we will download and build it.
if(HDF5_ROOT)
  message(STATUS "HDF5 root: ${HDF5_ROOT}")
  find_package(HDF5 REQUIRED)
else()
  find_package(HDF5)
endif()

if(NOT HDF5_FOUND AND NOT HDF5_ROOT)
  message(STATUS "HDF5 not found. Will attempt to download and build.")
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

# Include the HDF5 libraries
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
message(STATUS "HDF5 libraries shared: ${HDF5_LIBRARIES_SHARED}")
message(STATUS "HDF5 libraries static: ${HDF5_LIBRARIES_STATIC}")

# Set package properties
set_package_properties(HDF5 PROPERTIES
  DESCRIPTION "HDF5 library"
  TYPE REQUIRED
  PURPOSE "Enables HDF5 support"
)

############################################################
# Setup MOAB                                               #
############################################################

# Check if MOAB is requested by user
if (ENABLE_MOAB)
  # If MOAB is installed in a non-standard location, set MOAB_ROOT
  # to the root directory of the MOAB installation.
  # For example: export MOAB_ROOT=/path/to/moab
  if (MOAB_ROOT)
    find_package(MOAB REQUIRED)
  else()
    find_package(MOAB)
  endif()

  if(MOAB_FOUND)
    include_directories(${MOAB_INCLUDE_DIRS})
    get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
    link_directories(${MOAB_LIBRARY_DIRS})
  else()
    message(WARNING "MOAB was not found. Will not use MOAB.")
    set(ENABLE_MOAB OFF)
  endif()
endif()

# Set package properties
set_package_properties(MOAB PROPERTIES
  DESCRIPTION "Mesh-Oriented datABase (MOAB)"
  TYPE RECOMMENDED
  PURPOSE "Enables MOAB support"
)

############################################################
# Setup DAGMC                                              #
############################################################

# Check if DAGMC is requested by user
if (ENABLE_DAGMC AND MOAB_FOUND)
  # If DAGMC is installed in a non-standard location, set DAGMC_ROOT
  # to the root directory of the DAGMC installation.
  # For example: export DAGMC_ROOT=/path/to/dagmc
  if (DAGMC_ROOT)
    find_package(DAGMC REQUIRED)
  else()
    find_package(DAGMC)
  endif()

  if(DAGMC_FOUND)
    include_directories(${DAGMC_INCLUDE_DIRS})
    link_directories(${DAGMC_LIBRARY_DIRS})
  else()
    message(WARNING "DAGMC was not found. Will not use DAGMC.")
    set(ENABLE_DAGMC OFF)
  endif()
elseif(ENABLE_DAGMC)
  message(WARNING "MOAB was not found. Will not use DAGMC.")
  set(ENABLE_DAGMC OFF)
endif()

# Set package properties
set_package_properties(DAGMC PROPERTIES
  DESCRIPTION "Direct Accelerated Geometry Monte Carlo (DAGMC)"
  TYPE RECOMMENDED
  PURPOSE "Enables DAGMC support"
)

############################################################
# Setup Python, NumPy and Cython                           #
############################################################

# Find Python and NumPy
find_package(
  Python 
  COMPONENTS Interpreter Development.Module NumPy 
  REQUIRED
)

# Set package properties
set_package_properties(Python PROPERTIES
  DESCRIPTION "Python interpreter"
  TYPE REQUIRED
  PURPOSE "Enables Python support"
)

# Set package properties
set_package_properties(NumPy PROPERTIES
  DESCRIPTION "NumPy"
  TYPE REQUIRED
  PURPOSE "Enables NumPy support"
)

find_program(CYTHON "cython")

# Set package properties
set_package_properties(Cython PROPERTIES
  DESCRIPTION "Cython"
  TYPE REQUIRED
  PURPOSE "Enables Cython support"
)

############################################################
# Setup Fortran, f2py and LAPACK                           #
############################################################

# Call macro to setup Fortran
pyne_setup_fortran()

# Check if spatial solver is requested by user
if (BUILD_SPATIAL_SOLVER)
  # Find f2py, if building spatial solver
  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Try to find LAPACK and BLAS
  # If not found, no warry, we'll try to download and build them
  find_package(LAPACK)
  if(NOT BLAS_FOUND AND NOT LAPACK_FOUND)
    message(STATUS "BLAS not found. Will attempt to download and build.")
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
    MESSAGE(STATUS "LAPACK not found. Will attempt to download and build.")
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

# Set package properties
set_package_properties(F2PY PROPERTIES
  DESCRIPTION "Utility that connects Python and Fortran code"
  TYPE REQUIRED
  PURPOSE "Enables f2py support"
)
set_package_properties(LAPACK PROPERTIES
  DESCRIPTION "LAPACK library"
  TYPE RECOMMENDED
  PURPOSE "Enables LAPACK support"
)
set_package_properties(BLAS PROPERTIES
  DESCRIPTION "BLAS library"
  TYPE OPTIONAL
  PURPOSE "Required by LAPACK build"
)

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

feature_summary(
  WHAT ALL #ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND
  DESCRIPTION "PyNE Configuration Summary"
  INCLUDE_QUIET_PACKAGES
  FILENAME ${CMAKE_BINARY_DIR}/summary.log
)
message(STATUS "Summary log written to ${CMAKE_BINARY_DIR}/summary.log")