# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.0.0)

# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Setup the project
project(pyne)

# check for and enable c++11 support
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
ENDIF()

include(PyneMacros)
pyne_print_logo()  # Beware of dragons
pyne_set_platform()
pyne_setup_cxx()
pyne_set_asm_platform()
pyne_set_build_spatial_solver()
pyne_set_build_type()
pyne_set_fast_compile()
pyne_configure_rpath()

# enable assembly
enable_language(ASM)

# This makes all the libraries build as SHARED
set(BUILD_SHARED_LIBS true)
message("-- CMake Install Prefix: ${CMAKE_INSTALL_PREFIX}")


#
# Allow for dependnecies to exist in non-install location
#
if (NOT DEPS_ROOT_DIR)
  get_filename_component(compdir ${CMAKE_CXX_COMPILER} DIRECTORY)
  get_filename_component(DEPS_ROOT_DIR ${compdir} DIRECTORY)
endif (NOT DEPS_ROOT_DIR)
SET(DEPS_HINTS HINTS "${DEPS_ROOT_DIR}")
SET(DEPS_BIN_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/bin")
SET(DEPS_LIB_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/lib")
SET(DEPS_INCLUDE_HINTS HINTS "${DEPS_ROOT_DIR}/include")
MESSAGE("-- Dependency Root Dir (DEPS_ROOT_DIR): ${DEPS_ROOT_DIR}")
MESSAGE("-- Dependency Hints (DEPS_HINTS): ${DEPS_HINTS}")
MESSAGE("-- Dependency Binary Hints (DEPS_BIN_HINTS): ${DEPS_BIN_HINTS}")
MESSAGE("-- Dependency Library Hints (DEPS_LIB_HINTS): ${DEPS_LIB_HINTS}")
MESSAGE("-- Dependency Include Hints (DEPS_INCLUDE_HINTS): ${DEPS_INCLUDE_HINTS}")


#
# Find Libraries
#

# Find HDF5
MESSAGE("-- HDF5 Root: ${HDF5_ROOT}")
if(HDF5_ROOT)
  FIND_PACKAGE(HDF5 REQUIRED)
else()
  FIND_PACKAGE(HDF5)
endif()
if(NOT HDF5_FOUND AND NOT HDF5_ROOT)
  MESSAGE("--    HDF5 not found on empty HDF5_ROOT, "
          "setting to DEPS_ROOT_DIR: ${DEPS_ROOT_DIR}")
  set(HDF5_ROOT "${DEPS_ROOT_DIR}")
  FIND_PACKAGE(HDF5 REQUIRED)
endif()
ADD_DEFINITIONS(${HDF5_DEFINITIONS})
include_directories("${HDF5_INCLUDE_DIRS}")
if(NOT HDF5_LIBRARY_DIRS STREQUAL "")
  link_directories(${HDF5_LIBRARY_DIRS})
endif()
set(LIBS_HDF5 ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
MESSAGE("--    HDF5 Include directory: ${HDF5_INCLUDE_DIRS}")
MESSAGE("--    HDF5 Library directories: ${HDF5_LIBRARY_DIRS}")
MESSAGE("--    HDF5 Libraries: ${HDF5_C_LIBRARIES}")


# Look for MOAB if requested
if(WITH_MOAB)
  # user may have set a search path
  if(NOT MOAB_ROOT)
    set(MOAB_ROOT "${DEPS_ROOT_DIR}")
  endif()
  MESSAGE("-- MOAB Root: ${MOAB_ROOT}")
  find_package(MOAB)
  message("-- MOAB Found: ${MOAB_FOUND}")
  if(${MOAB_FOUND})
    include_directories(${MOAB_INCLUDE_DIRS})
    get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
    if(APPLE)
      link_directories(${MOAB_LIBRARY_DIRS}/lib)
    else(APPLE)
      link_directories(${MOAB_LIBRARY_DIRS})
    endif(APPLE)
  endif(${MOAB_FOUND})
endif()

# Look for DAGMC if requested
if(WITH_DAGMC)
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
endif()

# Find Python
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
if(PYTHON_VERSION_MAJOR EQUAL 2)
  set(SITE_PACKAGES_CODE
      "import site; print site.getusersitepackages().replace(site.USER_BASE, '${CMAKE_INSTALL_PREFIX}', 1)")
else(PYTHON_VERSION_MAJOR EQUAL 2)
  set(SITE_PACKAGES_CODE
      "import site; print(site.getsitepackages(['${CMAKE_INSTALL_PREFIX}'])[0])")
endif(PYTHON_VERSION_MAJOR EQUAL 2)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "${SITE_PACKAGES_CODE}"
                OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
message("-- PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")
message("-- PYTHON_SITE_PACKAGES: ${PYTHON_SITE_PACKAGES}")


# Find NumPy
find_package(Numpy REQUIRED)
include_directories("${NUMPY_INCLUDE_DIR}")


pyne_setup_fortran()

# Find f2py, if building spatial solver
IF(BUILD_SPATIAL_SOLVER)
  find_package(F2py REQUIRED)
  message("-- F2PY Executable: ${F2PY_EXECUTABLE}")
  message("-- F2PY Version: ${F2PY_VERSION}")

  FIND_PACKAGE(LAPACK REQUIRED)
  MESSAGE("-- Found LAPACK Linker Flags: ${LAPACK_LINKER_FLAGS}")
  MESSAGE("-- Found LAPACK Libraries: ${LAPACK_LIBRARIES}")
  MESSAGE("-- Found BLAS Libraries: ${BLAS_LIBRARIES}")
ENDIF(BUILD_SPATIAL_SOLVER)


# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of cmake/UseCython.cmake.
include(UseCython)
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/pyne")

# Add JsonCpp Flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_IS_AMALGAMATION")

#
# Build and install!
#
message(STATUS "LIBS: ${LIBS}")
# With CMake, a clean separation can be made between the source tree and the
# build tree.  When all source is compiled, as with pure C/C++, the source is
# no-longer needed in the build tree.  However, with pure *.py source, the
# source is processed directly.  To handle this, we reproduce the availability
# of the source files in the build tree.
add_custom_target(ReplicatePythonSourceTree ALL "${CMAKE_COMMAND}" -P
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# add sub dirs
add_subdirectory(src)
add_subdirectory(pyne)
add_subdirectory(pyne/xs)

# install python package
set(pyinstallscript "${CMAKE_CURRENT_BINARY_DIR}/cmake/SetupSubPyInstall.cmake")
configure_file(cmake/SetupSubPyInstall.cmake.in ${pyinstallscript} @ONLY)
install(SCRIPT ${pyinstallscript})
