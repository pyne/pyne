# Find the Cython compiler.
#
# This code sets the following variables:
#
#  CYTHON_EXECUTABLE
#  CYTHON_VERSION
#  CYTHON_VERSION_MAJOR
#  CYTHON_VERSION_MINOR
#  CYTHON_VERSION_MICRO
#
# See also UseCython.cmake

#=============================================================================
# Copyright 2011 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================

# Use the Cython executable that lives next to the Python executable
# if it is a local installation.
find_package( PythonInterp )
if( PYTHONINTERP_FOUND )
  get_filename_component( _python_path ${PYTHON_EXECUTABLE} PATH )
  find_program(CYTHON_EXECUTABLE
    NAMES cython cython.bat cython3 cython-2.7 cython-3.3 cython-3.4 cython-3.5 cython-3.6 cython-3.7
    HINTS ENV PATH ${_python_path}
    ${DEPS_BIN_HINTS}
    )
else()
  find_program(CYTHON_EXECUTABLE
    NAMES cython cython.bat cython3 cython-2.7 cython-3.3 cython-3.4 cython-3.5 cython-3.6 cython-3.7
    ${DEPS_BIN_HINTS}
    )
endif()

if( ${CYTHON_EXECUTABLE} STREQUAL "" OR ${CYTHON_EXECUTABLE} STREQUAL "CYTHON_EXECUTABLE-NOTFOUND" )
  SET( Cython_FOUND FALSE )
else()
  SET( Cython_FOUND TRUE )

  # get the version strings
  execute_process(COMMAND "${CYTHON_EXECUTABLE}" "-V"
                  ERROR_VARIABLE CYTHON_VERSION_RTN
                  ERROR_STRIP_TRAILING_WHITESPACE)
  string(REPLACE " " ";" CYTHON_VERSION_RTN_LIST ${CYTHON_VERSION_RTN})
  list(GET CYTHON_VERSION_RTN_LIST -1 CYTHON_VERSION)
  string(REPLACE "." ";" CYTHON_VERSION_LIST ${CYTHON_VERSION})
  list(LENGTH CYTHON_VERSION_LIST CYTHON_VERSION_N)
  list(GET CYTHON_VERSION_LIST 0 CYTHON_VERSION_MAJOR)
  list(GET CYTHON_VERSION_LIST 1 CYTHON_VERSION_MINOR)
  if(CYTHON_VERSION_N GREATER 2)
    list(GET CYTHON_VERSION_LIST 2 CYTHON_VERSION_MICRO)
  else(CYTHON_VERSION_N GREATER 2)
    set(CYTHON_VERSION_MICRO 0)
  endif(CYTHON_VERSION_N GREATER 2)
endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cython REQUIRED_VARS CYTHON_EXECUTABLE)

mark_as_advanced( CYTHON_EXECUTABLE )
