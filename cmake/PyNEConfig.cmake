# PyNEConfig.cmake
# 
# Configuration file for integrating PyNE with CMake projects.
# This script finds the Python interpreter, ensures its availability,
# and retrieves necessary PyNE details such as version, include paths, 
# and library paths. It also verifies that the Python environment is 
# properly set up for PyNE and checks for potential compatibility issues.

# This CMake script extracts the following PyNE details:
# - PyNE version (PyNE_VERSION)
# - PyNE include paths (PyNE_INCLUDE_DIRS)
# - PyNE library paths (PyNE_LIBRARY_DIRS)
# - PyNE libraries (PyNE_LIBRARIES)
# - PyNE extra libraries (PyNE_EXTRA_LIBRARIES)

# Find the Python interpreter and ensure it's available.
find_package(Python COMPONENTS Interpreter REQUIRED)

# Function to run Python commands and validate their execution.
function(run_python_command output_var command)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "${command}"
    OUTPUT_VARIABLE ${output_var}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result
  )
  # Check if the command was successful
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Failed to run Python command: ${command}")
  else()
    # Add the output variable to the parent scope
    set(${output_var} "${${output_var}}" PARENT_SCOPE)
  endif()
endfunction()

# Extract PyNE version, include paths, library paths, and extra libraries
run_python_command(PyNE_VERSION "import pyne; print(pyne.__version__)")
run_python_command(PyNE_INCLUDE_DIRS "import pyne; print(pyne.include_path)")
run_python_command(PyNE_LIBRARY_DIRS "import pyne; print(pyne.lib_path)")
run_python_command(PyNE_EXTRA_LIBRARIES "import pyne; print(' '.join(pyne.extra_lib))")

# Check if the wheel was repaired using auditwheel or delocate
if(PyNE_EXTRA_LIBRARIES)
  message(FATAL_ERROR
      "This build of PyNE is not supported. "
      "It appears that the wheel was repaired using tools like auditwheel or delocate, "
      "that modifies the shared libraries, which may cause problems.\n"
      "PYNE_EXTRA_LIBRARIES is not empty: ${PyNE_EXTRA_LIBRARIES}.\n"
      "To resolve this, please build PyNE from scratch. "
      "For more information, visit: https://pyne.io\n"
  )
endif()

# Find the library
find_library(PyNE_LIBRARIES 
  NAMES pyne 
  PATHS ${PyNE_LIBRARY_DIRS} 
  NO_DEFAULT_PATH
  )

# Include standard argument handling for finding packages
include(FindPackageHandleStandardArgs)

# Validates that the necessary variables (PyNE include paths and library) are set
find_package_handle_standard_args(PyNE
  REQUIRED_VARS PyNE_LIBRARIES PyNE_INCLUDE_DIRS
  VERSION_VAR PyNE_VERSION
  )