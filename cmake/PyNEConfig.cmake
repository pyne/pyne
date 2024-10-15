find_package(Python COMPONENTS Interpreter REQUIRED)

# Function to run Python commands and validate their execution.
function(run_python_command output_var command)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "${command}"
    OUTPUT_VARIABLE ${output_var}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Failed to run Python command: ${command}")
  endif()
endfunction()

# Extract PyNE paths and validate.
run_python_command(PYNE_VERSION "import pyne; print(pyne.__version__)")
run_python_command(PYNE_INCLUDE_DIRS "import pyne; print(pyne.include_path)")
run_python_command(PYNE_LIBRARY_DIRS "import pyne; print(pyne.lib_path)")
run_python_command(PYNE_EXTRA_LIBRARIES "import pyne; print(' '.join(pyne.extra_lib))")

if(PYNE_EXTRA_LIBRARIES)
  message(FATAL_ERROR
      "This build of PyNE is not supported. "
      "It appears that the wheel was repaired using tools like auditwheel or delocate, "
      "that modifies the shared libraries, which may cause problems.\n"
      "PYNE_EXTRA_LIBRARIES is not empty: ${PyNE_EXTRA_LIBRARIES}.\n"
      "To resolve this, please build PyNE from scratch. "
      "For more information, visit: https://pyne.io\n"
  )
endif()