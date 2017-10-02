# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.0.0)

# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Setup the project
project(pyne)
include(PyneMacros)
pyne_print_logo()  # Beware of dragons
pyne_set_platform()
pyne_set_build_spatial_solver()
pyne_set_build_type()
pyne_configure_rpath()

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
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS HL)
ADD_DEFINITIONS(${HDF5_DEFINITIONS})
include_directories("${HDF5_INCLUDE_DIRS}")
if(NOT HDF5_LIBRARY_DIRS STREQUAL "")
  link_directories(${HDF5_LIBRARY_DIRS})
endif()
set(LIBS_HDF5 ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
MESSAGE("--    HDF5 Root: ${HDF5_ROOT}")
MESSAGE("--    HDF5 Include directory: ${HDF5_INCLUDE_DIRS}")
MESSAGE("--    HDF5 Library directories: ${HDF5_LIBRARY_DIRS}")
MESSAGE("--    HDF5 Libraries: ${HDF5_C_LIBRARIES}")
MESSAGE("--    HDF5 High Level Libraries: ${HDF5_C_HL_LIBRARIES}")


# Find MOAB library
find_package(MOAB)
message("-- MOAB Found: ${MOAB_FOUND}")
if(${MOAB_FOUND})
  include_directories(${MOAB_INCLUDE_DIR})
  get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
  if(APPLE)
    link_directories(${MOAB_LIBRARY_DIRS}/lib)
  else(APPLE)
    link_directories(${MOAB_LIBRARY_DIRS})
  endif(APPLE)
endif(${MOAB_FOUND})


# Find Python
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                        "import site; print(site.getsitepackages(['${CMAKE_INSTALL_PREFIX}'])[0])"
                OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
message("-- PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")
message("-- PYTHON_SITE_PACKAGES: ${PYTHON_SITE_PACKAGES}")


# Find NumPy
find_package(Numpy REQUIRED)
include_directories("${NUMPY_INCLUDE_DIR}")


# Find f2py, if building spatial solver
IF(BUILD_SPATIAL_SOLVER)
  pyne_setup_fortran()
  find_package(F2py REQUIRED)
  message("-- F2PY Executable: ${F2PY_EXECUTABLE}")
  message("-- F2PY Version: ${F2PY_VERSION}")
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
