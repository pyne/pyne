# Once done this will define
#  Numpy_FOUND        - system has NumPy
#  NUMPY_INCLUDE_DIR  - the NumPy include directory
#  NUMPY_VERSION      - the version of NumPy found

find_package(Python COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
    if(NOT NUMPY_VERSION_STRING)
        # Get numpy include directory and version
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import numpy; print(numpy.get_include()); print(numpy.__version__)"
                        RESULT_VARIABLE _NUMPY_SEARCH_SUCCESS
                        OUTPUT_VARIABLE _NUMPY_VALUES_OUTPUT
                        ERROR_VARIABLE _NUMPY_ERROR_VALUE
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

        # If numpy is found, set the include dir and version variables
        if(_NUMPY_SEARCH_SUCCESS EQUAL 0)
            string(REGEX REPLACE "\n" ";" _NUMPY_VALUES_LIST ${_NUMPY_VALUES_OUTPUT})
            list(GET _NUMPY_VALUES_LIST 0 NUMPY_INCLUDE_DIR)
            list(GET _NUMPY_VALUES_LIST 1 NUMPY_VERSION_STRING)
        endif()
    endif()

    find_package_handle_standard_args(Numpy REQUIRED_VARS NUMPY_INCLUDE_DIR VERSION_VAR NUMPY_VERSION_STRING)

    if(Numpy_FOUND AND NOT TARGET Numpy::Numpy)
        add_library(Numpy::Numpy UNKNOWN IMPORTED)
        set_target_properties(Numpy::Numpy PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NUMPY_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(NUMPY_INCLUDE_DIR NUMPY_VERSION_STRING)