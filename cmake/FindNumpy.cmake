# - Find Numpy
# NumPy is the fundamental package needed for scientific computing with Python
# www.numpy.scipy.org
#
# The module defines the following variables:
#  NUMPY_FOUND - the system has numpy
#  NUMPY_INCLUDE_DIR - where to find numpy/arrayobject.h
#  NUMPY_INCLUDE_DIRS - numpy include directories
#  NUMPY_VERSION_STRING - version (ex. 1.2.3)
#  NUMPY_MAJOR_VERSION - major version (ex. 1)
#  NUMPY_MINOR_VERSION - minor version (ex. 2)
#  NUMPY_PATCH_VERSION - patch version (ex. 3)

#=============================================================================
# Copyright 2005-2012 EDF-EADS-Phimeca
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

# set NUMPY_INCLUDE_DIR
find_package ( PythonInterp )

if ( PYTHONINTERP_FOUND )
  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
                    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
endif ()

# set NUMPY_INCLUDE_DIRS
set ( NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR} )

# version
if ( PYTHONINTERP_FOUND )
  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.__version__)"
                    OUTPUT_VARIABLE NUMPY_VERSION_STRING
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

  if ( NUMPY_VERSION_STRING )
    string ( REGEX REPLACE "([0-9]+)\\..*" "\\1" NUMPY_MAJOR_VERSION ${NUMPY_VERSION_STRING} )
    string ( REGEX REPLACE "[0-9]+\\.([0-9]+).*" "\\1" NUMPY_MINOR_VERSION ${NUMPY_VERSION_STRING} )
    string ( REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" NUMPY_PATCH_VERSION ${NUMPY_VERSION_STRING} )
  endif ()

endif ()

# check version
set ( _NUMPY_VERSION_MATCH TRUE )
if ( Numpy_FIND_VERSION AND NUMPY_VERSION )
  if ( Numpy_FIND_VERSION_EXACT )
    if ( Numpy_FIND_VERSION VERSION_EQUAL NUMPY_VERSION_STRING )
    else()
      set ( _NUMPY_VERSION_MATCH FALSE)
    endif ()
  else ()
    if ( Numpy_FIND_VERSION VERSION_GREATER NUMPY_VERSION_STRING )
      set ( _NUMPY_VERSION_MATCH FALSE )
    endif ()
  endif ()
endif ()

message("-- NUMPY_VERSION_STRING = ${NUMPY_VERSION_STRING}")

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )
find_package_handle_standard_args ( Numpy DEFAULT_MSG
  NUMPY_VERSION_STRING
  _NUMPY_VERSION_MATCH
  NUMPY_INCLUDE_DIR
  NUMPY_INCLUDE_DIRS
)

mark_as_advanced (
  NUMPY_VERSION_STRING
  NUMPY_MAJOR_VERSION
  NUMPY_MINOR_VERSION
  NUMPY_PATCH_VERSION
  NUMPY_INCLUDE_DIR
  NUMPY_INCLUDE_DIRS
)
