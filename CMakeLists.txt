# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.15)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add JsonCpp Flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_IS_AMALGAMATION")

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
option ( DOWNLOAD_HDF5              "Should PyNE download HDF5?"                  ON )
option ( ENABLE_EIGEN3              "Should PyNE use Eigen3?"                     ON )
option ( DOWNLOAD_EIGEN3            "Should PyNE download Eigen3?"                ON )
option ( ENABLE_LAPACK              "Should PyNE use LAPACK?"                     ON )
option ( DOWNLOAD_LAPACK            "Should PyNE download LAPACK?"                ON )
option ( ENABLE_MOAB                "Should PyNE use MOAB?"                       ON )
option ( DOWNLOAD_MOAB              "Should PyNE download MOAB?"                  ON )
option ( ENABLE_DAGMC               "Should PyNE use DAGMC?"                      ON )
option ( DOWNLOAD_DAGMC             "Should PyNE download DAGMC?"                 ON )
option ( ENABLE_SPATIAL_SOLVERS     "Should build AHOT spatial solvers?"          ON )
option ( ENABLE_ENSDF_PROCESSING    "Should build ENSDF processing tools?"        ON )
option ( ENABLE_FAST_COMPILE        "Should PyNE fast compile?"                   ON )

# Add feature info
add_feature_info("Download HDF5" DOWNLOAD_HDF5 "Download HDF5")
add_feature_info("Enable Eigen3" ENABLE_EIGEN3 "Use Eigen3")
add_feature_info("Download Eigen3" DOWNLOAD_EIGEN3 "Download Eigen3")
add_feature_info("Enable LAPACK" ENABLE_LAPACK "Use LAPACK")
add_feature_info("Download LAPACK" DOWNLOAD_LAPACK "Download LAPACK")
add_feature_info("Enable MOAB" ENABLE_MOAB "Use MOAB")
add_feature_info("Download MOAB" DOWNLOAD_MOAB "Download MOAB")
add_feature_info("Enable DAGMC" ENABLE_DAGMC "Use DAGMC")
add_feature_info("Download DAGMC" DOWNLOAD_DAGMC "Download DAGMC")
add_feature_info("AHOT" ENABLE_SPATIAL_SOLVERS "AHOT spatial solvers")
add_feature_info("ENSDF" ENABLE_ENSDF_PROCESSING "ENSDF processing tools")
add_feature_info("Enable Fast Compile" ENABLE_FAST_COMPILE "Use pre-compiled binaries")

# Include PyneMacros file
include(PyneMacros)
print_pyne_logo()
pyne_configure_skbuild()
pyne_configure_rpath()
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
# If HDF5 is not installed, set DOWNLOAD_HDF5 to ON.
# PyNE will download and build HDF5 if DOWNLOAD_HDF5 is ON.
if(NOT DOWNLOAD_HDF5)
  message(STATUS "HDF5 root: ${HDF5_ROOT}")
  find_package(HDF5 REQUIRED)
else()
  message(STATUS "HDF5 will be downloaded and installed")

  # Configure HDF5 
  set(HDF5_VERSION "1.10.11")
  
  # Replace the dots with semicolons to create a list
  string(REPLACE "." ";" HDF5_VERSION_LIST ${HDF5_VERSION})

  # Convert the list into individual variables
  list(GET HDF5_VERSION_LIST 0 HDF5_VERSION_MAJOR)
  list(GET HDF5_VERSION_LIST 1 HDF5_VERSION_MINOR)
  list(GET HDF5_VERSION_LIST 2 HDF5_VERSION_PATCH)
  set(HDF5_ROOT "${CMAKE_BINARY_DIR}/hdf5")
  ExternalProject_Add(hdf5-project
    PREFIX ${HDF5_ROOT}
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DBUILD_STATIC_LIBS:BOOL=OFF
    DOWNLOAD_EXTRACT_TIMESTAMP true
    BUILD_BYPRODUCTS "${HDF5_ROOT}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
  set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
  set(HDF5_LIBRARY_DIRS "${HDF5_ROOT}/lib")

  # HDF5 is a shared library, so we need to install it
  install(DIRECTORY ${HDF5_LIBRARY_DIRS}
    DESTINATION .
    FILES_MATCHING PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
    )
endif()

# Include the HDF5 libraries
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")

link_directories(${HDF5_LIBRARY_DIRS})
message(STATUS "HDF5 library dirs: ${HDF5_LIBRARY_DIRS}")

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
# Setup Eigen3                                             #
############################################################

# Check if Eigen3 is requested by user
if(ENABLE_EIGEN3)
  # If Eigen3 is installed in a non-standard location, set EIGEN3_ROOT
  # to the root directory of the Eigen3 installation.
  # For example: export EIGEN3_ROOT=/path/to/eigen3
  # If Eigen is not installed, set DOWNLOAD_EIGEN3 to ON.
  # PyNE will download Eigen3 if DOWNLOAD_EIGEN3 is ON. (Recommended)
  if(NOT DOWNLOAD_EIGEN3)
    find_package(Eigen3 REQUIRED)
  else()
    message(STATUS "Eigen3 will be downloaded")

    # Configure Eigen3
    set(EIGEN3_VERSION "3.4.0")
    set(EIGEN3_ROOT "${CMAKE_BINARY_DIR}/eigen3")
    ExternalProject_Add(eigen3-project
      PREFIX ${EIGEN3_ROOT}
      URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN3_VERSION}/eigen-${EIGEN3_VERSION}.tar.bz2
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      DOWNLOAD_EXTRACT_TIMESTAMP true
    )
    set(EIGEN3_INCLUDE_DIRS "${EIGEN3_ROOT}/include/eigen3")  
  endif()

    # Include the Eigen3 libraries
    include_directories("${EIGEN3_INCLUDE_DIRS}")
    message(STATUS "Eigen3 include dirs: ${EIGEN3_INCLUDE_DIRS}")
endif()

# Set package properties
set_package_properties(EIGEN3 PROPERTIES
  DESCRIPTION "Eigen3 library"
  TYPE REQUIRED
  PURPOSE "Enables Eigen3 support"
  )


############################################################
# Setup MOAB                                               #
############################################################

# Check if MOAB is requested by user
if(ENABLE_MOAB)
  if(NOT ENABLE_EIGEN3)
    message(FATAL_ERROR "Eigen3 is required to build MOAB.")
  endif()
  # If MOAB is installed in a non-standard location, set MOAB_ROOT
  # to the root directory of the MOAB installation.
  # For example: export MOAB_ROOT=/path/to/moab
  # If MOAB is not installed, set DOWNLOAD_MOAB to ON.
  # PyNE will download and build MOAB if DOWNLOAD_MOAB is ON. (Recommended)
  if(NOT DOWNLOAD_MOAB)
    find_package(MOAB REQUIRED)
  else()
    message(STATUS "MOAB will be downloaded and built")

    # Configure MOAB version
    set(MOAB_VERSION "5.5.1")
    set(MOAB_ROOT "${CMAKE_BINARY_DIR}/moab")
    ExternalProject_Add(moab-project
      PREFIX ${MOAB_ROOT}
      GIT_REPOSITORY https://bitbucket.org/fathomteam/moab.git
      GIT_TAG ${MOAB_VERSION}
      CMAKE_ARGS 
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS:BOOL=ON
        -DENABLE_HDF5:BOOL=ON
        -DHDF5_ROOT:PATH=${HDF5_ROOT}
        -DEIGEN3_DIR:PATH=${EIGEN3_INCLUDE_DIRS}
        -DENABLE_BLASLAPACK:BOOL=OFF
        -DENABLE_FORTRAN:BOOL=OFF
      DOWNLOAD_EXTRACT_TIMESTAMP true
      BUILD_BYPRODUCTS "${MOAB_ROOT}/${CMAKE_INSTALL_LIBDIR}/libMOAB${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )
    set(MOAB_INCLUDE_DIRS "${MOAB_ROOT}/include")
    set(MOAB_LIBRARY_DIRS "${MOAB_ROOT}/${CMAKE_INSTALL_LIBDIR}")

    # HDF5 is needed to be build before MOAB
    if(DOWNLOAD_HDF5)
      add_dependencies(moab-project hdf5-project)
    endif()

    # EIGEN3 is needed to be build before MOAB
    if(DOWNLOAD_EIGEN3)
      add_dependencies(moab-project eigen3-project)
    endif()

    # MOAB is a shared library, so we need to install it
    install(DIRECTORY ${MOAB_LIBRARY_DIRS}
      DESTINATION .
      FILES_MATCHING PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
      )
  endif()
  include_directories(${MOAB_INCLUDE_DIRS})
  link_directories(${MOAB_LIBRARY_DIRS})
  message(STATUS "MOAB include dirs: ${MOAB_INCLUDE_DIRS}")
  message(STATUS "MOAB library dirs: ${MOAB_LIBRARY_DIRS}")
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
if(ENABLE_DAGMC)
  if(NOT ENABLE_EIGEN3)
    message(FATAL_ERROR "Eigen3 is required to build DAGMC python module.")
  endif()
  if(NOT ENABLE_MOAB)
    message(FATAL_ERROR "MOAB is required to build DAGMC.")
  endif()
  # If DAGMC is installed in a non-standard location, set DAGMC_ROOT
  # to the root directory of the DAGMC installation.
  # For example: export DAGMC_ROOT=/path/to/dagmc
  # If DAGMC is not installed, set DOWNLOAD_DAGMC to ON.
  # PyNE will download and build DAGMC if DOWNLOAD_DAGMC is ON. (Recommended)
  if(NOT DOWNLOAD_DAGMC)
    find_package(DAGMC REQUIRED)
  else()
    message(STATUS "DAGMC will be downloaded and built")

    # Configure DAGMC
    set(DAGMC_VERSION "3.2.3")
    set(DAGMC_ROOT "${CMAKE_BINARY_DIR}/DAGMC")
    ExternalProject_Add(dagmc-project
      PREFIX ${DAGMC_ROOT}
      GIT_REPOSITORY https://github.com/svalinn/DAGMC.git
      GIT_TAG v${DAGMC_VERSION}
      CMAKE_ARGS 
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DMOAB_DIR:PATH=${MOAB_ROOT}
        -DBUILD_STATIC_LIBS:BOOL=OFF
        -DBUILD_UWUW:BOOL=OFF
        -DBUILD_TALLY:BOOL=OFF
        -DBUILD_MAKE_WATERTIGHT:BOOL=OFF
        -DBUILD_OVERLAP_CHECK:BOOL=OFF
        -DBUILD_TESTS:BOOL=OFF
      BUILD_BYPRODUCTS "${DAGMC_ROOT}/lib/libdagmc${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )
    set(DAGMC_INCLUDE_DIRS "${DAGMC_ROOT}/include")
    set(DAGMC_LIBRARY_DIRS "${DAGMC_ROOT}/lib")
    
    # MOAB is needed to be build before DAGMC
    if(DOWNLOAD_MOAB)
      add_dependencies(dagmc-project moab-project)
    endif()

    # DAGMC is a shared library, so we need to install it
    install(DIRECTORY ${DAGMC_LIBRARY_DIRS}
      DESTINATION .
      FILES_MATCHING PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
      )
  endif()

  include_directories(${DAGMC_INCLUDE_DIRS})
  link_directories(${DAGMC_LIBRARY_DIRS})
  message(STATUS "DAGMC include dirs: ${DAGMC_INCLUDE_DIRS}")
  message(STATUS "DAGMC library dirs: ${DAGMC_LIBRARY_DIRS}")
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

find_program(CYTHON "cython")


############################################################
# Setup Fortran, LAPACK and f2py                           #
############################################################

# Call macro to setup Fortran
pyne_setup_fortran()

if(ENABLE_LAPACK)
  # If LAPACK is installed in a non-standard location, set LAPACK_ROOT
  # to the root directory of the LAPACK installation.
  # For example: export LAPACK_ROOT=/path/to/lapack
  # If LAPACK is not installed, set LAPACK_ROOT to OFF.
  # PyNE will download and build LAPACK if LAPACK_ROOT is ON. (Recommended)
  if(NOT DOWNLOAD_LAPACK)
    find_package(LAPACK REQUIRED)
  else()
    message(STATUS "LAPACK will be downloaded and built")

    # Configure LAPACK
    set(LAPACK_VERSION "3.12.0")
    set(LAPACK_ROOT "${CMAKE_BINARY_DIR}/lapack")
    ExternalProject_Add(lapack-project
      PREFIX ${LAPACK_ROOT}
      URL https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${LAPACK_VERSION}.tar.gz
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> 
        -DBUILD_SHARED_LIBS:BOOL=ON
        -DUSE_OPTIMIZED_BLAS:BOOL=ON
      BUILD_BYPRODUCTS "${LAPACK_ROOT}/lib/liblapack${CMAKE_SHARED_LIBRARY_SUFFIX}"
      )
    set(LAPACK_INCLUDE_DIRS "${LAPACK_ROOT}/include")
    set(LAPACK_LIBRARY_DIRS "${LAPACK_ROOT}/lib")

    # LAPACK is a shared library, so we need to install it
    install(DIRECTORY ${LAPACK_LIBRARY_DIRS}
      DESTINATION .
      FILES_MATCHING PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
      )
  endif()
  include_directories(${LAPACK_INCLUDE_DIRS})
  link_directories(${LAPACK_LIBRARY_DIRS})
  message(STATUS "LAPACK include dirs: ${LAPACK_INCLUDE_DIRS}")
  message(STATUS "LAPACK library dirs: ${LAPACK_LIBRARY_DIRS}")
endif()

set_package_properties(LAPACK PROPERTIES
  DESCRIPTION "LAPACK library"
  TYPE RECOMMENDED
  PURPOSE "Enables LAPACK support"
)

# Check if spatial solver is requested by user
if (BUILD_SPATIAL_SOLVER)
  if( NOT ENABLE_LAPACK)
    message(FATAL_ERROR "LAPACK is required to build spatial solvers.")
  endif()

  # Find f2py, if building spatial solver
  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import numpy.f2py; print(numpy.f2py.get_include())"
            OUTPUT_VARIABLE F2PY_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()


############################################################
# Setup Scripts                                            #
############################################################

# Collect scripts
file(GLOB SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/*")

# Filter scripts based on the OS
if(WIN32)
    list(FILTER SCRIPTS INCLUDE REGEX "\\.bat$")
else()
    list(FILTER SCRIPTS EXCLUDE REGEX "\\.bat$")
endif()

# Install scripts to bin directory
install(PROGRAMS ${SCRIPTS} DESTINATION ${CMAKE_INSTALL_BINDIR})


############################################################
# Add sub dirs                                             #
############################################################

# add sub dirs
add_subdirectory(src)
add_subdirectory(pyne)
add_subdirectory(pyne/xs)


############################################################
# Summary                                                  #
############################################################

# Write summary
feature_summary(
  WHAT ALL
  DESCRIPTION "PyNE Configuration Summary"
  INCLUDE_QUIET_PACKAGES
  FILENAME ${CMAKE_SOURCE_DIR}/summary.log
)
message(STATUS "Summary log written to ${CMAKE_SOURCE_DIR}/summary.log")