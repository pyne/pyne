# - Find PyNE libraries
# This module finds the libraries corresponding to the PyNE library, using the Python
# interpreter.
# This code sets the following variables:
#
#  PYNE_LIBS_FOUND            - have the PyNE libs been found
#  PYNE_PREFIX                - path to the PyNE installation
#  PYNE_LIBS_DIR              - path to the PyNE libs dir
#  PYNE_INCLUDE_DIR           - path to where PyNE header files are
#
# To use PyNE add the following lines to your main CMakeLists.txt file:
# 
#   # Find the Pyne installation
#   find_package(Pyne REQUIRED)
#   link_directories(${PYNE_LIBS_DIR})
#   include_directories(${PYNE_INCLUDE_DIR})
#
# Enjoy!

# Use the Python interpreter to find the libs.
if(Pyne_FIND_REQUIRED)
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

if(NOT PYTHONINTERP_FOUND)
    set(PYNE_LIBS_FOUND FALSE)
    return()
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "from pyne import pyne_config as pc; print(pc.prefix); print(pc.lib); print(pc.includes)"
    RESULT_VARIABLE _PYNE_SUCCESS
    OUTPUT_VARIABLE _PYNE_VALUES
    ERROR_VARIABLE _PYNE_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _PYNE_SUCCESS MATCHES 0)
    if(Pyne_FIND_REQUIRED)
        message(FATAL_ERROR
            "Pyne config failure:\n${_PYNE_ERROR_VALUE}")
    endif()
    set(PYNE_LIBS_FOUND FALSE)
    return()
endif()

# Convert the process output into a list
string(REGEX REPLACE ";" "\\\\;" _PYNE_VALUES ${_PYNE_VALUES})
string(REGEX REPLACE "\n" ";" _PYNE_VALUES ${_PYNE_VALUES})
list(GET _PYNE_VALUES 0 PYNE_PREFIX)
list(GET _PYNE_VALUES 1 PYNE_LIBS_DIR)
list(GET _PYNE_VALUES 2 PYNE_INCLUDE_DIR)

# Make sure all directory separators are '/'
string(REGEX REPLACE "\\\\" "/" PYNE_PREFIX ${PYNE_PREFIX})
string(REGEX REPLACE "\\\\" "/" PYNE_LIBS_DIR ${PYNE_LIBS_DIR})
string(REGEX REPLACE "\\\\" "/" PYNE_INCLUDE_DIR ${PYNE_INCLUDE_DIR})


MARK_AS_ADVANCED(
  PYNE_LIBS_DIR
  PYNE_INCLUDE_DIR
)

# We use PYNE_INCLUDE_DIR, PYNE_LIBS_DIR and PYNE_PREFIX for the
# cache entries because they are meant to specify the location of a single
# library. We now set the variables listed by the documentation for this
# module.
SET(PYNE_PREFIX      "${PYNE_PREFIX}")
SET(PYNE_LIBS_DIR    "${PYNE_LIBS_DIR}")
SET(PYNE_INCLUDE_DIR "${PYNE_INCLUDE_DIR}")


find_package_message(PYNE "Found PyNE: ${PYNE_PREFIX}" "${PYNE_PREFIX}")
