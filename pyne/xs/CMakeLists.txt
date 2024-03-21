# Include the source and build directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)

# Build C++ file
add_custom_command(
    OUTPUT models.cxx
    DEPENDS models.pyx
    VERBATIM
    COMMAND "${CYTHON}" --cplus 
            -I ${CMAKE_SOURCE_DIR}/src 
            -I ${CMAKE_BINARY_DIR}/src 
            -I ${CMAKE_CURRENT_SOURCE_DIR} 
            -I ${CMAKE_CURRENT_BINARY_DIR} 
            -I ${Python_NumPy_INCLUDE_DIRS}
            "${CMAKE_CURRENT_SOURCE_DIR}/models.pyx" 
            --output-file "${CMAKE_CURRENT_BINARY_DIR}/models.cxx"
)

# Creating a Python module
python_add_library(
    models MODULE "${CMAKE_CURRENT_BINARY_DIR}/models.cxx"
    WITH_SOABI
)

# Add the pyne library as a dependency
add_dependencies(models ${CMAKE_PROJECT_NAME})

# Link the target with the pyne library
target_link_libraries(models PRIVATE ${CMAKE_PROJECT_NAME})

# Install the target in the pyne directory
install(TARGETS models DESTINATION pyne/xs)