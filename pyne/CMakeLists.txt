set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# Include the source and build directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)

# Build normal cython modules
set(PYNE_CYTHON_MODULES
    ace
    bins
    data
    dtypes
    endf
    enrichment
    extra_types
    jsoncpp
    material
    material_library
    nucname
    particle
    pyne_config
    rxname
    source
    stlcontainers
    tally
    transmuters
    _utils
    )
if(MOAB_FOUND)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} source_sampling)
endif()

if(DAGMC_FOUND)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} dagmc)
endif()

if(NOT WIN32)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} cram)
endif()

message(STATUS "Modules to Cythonize: ${PYNE_CYTHON_MODULES}")

foreach(PYNE_CYTHON_MODULE ${PYNE_CYTHON_MODULES})
    # Cram is actually a C, not C++ file
    if(NOT "${PYNE_CYTHON_MODULE}" STREQUAL "cram")
        set(PYNE_CYTHON_ARGS --cplus)
        set(PYNE_CYTHON_EXT cxx)
    else()
        set(PYNE_CYTHON_ARGS)
        set(PYNE_CYTHON_EXT c)
    endif()

    # Add different libraries
    if("${PYNE_CYTHON_MODULE}" STREQUAL "dagmc")
      set(PYNE_CYTHON_LIBS ${CMAKE_PROJECT_NAME} MOAB dagmc)
    elseif("${PYNE_CYTHON_MODULE}" STREQUAL "source_sampling")
      set(PYNE_CYTHON_LIBS ${CMAKE_PROJECT_NAME} MOAB)
    else()
      set(PYNE_CYTHON_LIBS ${CMAKE_PROJECT_NAME})
    endif()

    # Build C++/C file
    add_custom_command(
        OUTPUT ${PYNE_CYTHON_MODULE}.${PYNE_CYTHON_EXT}
        DEPENDS ${PYNE_CYTHON_MODULE}.pyx
        VERBATIM
        COMMAND "${CYTHON}" ${PYNE_CYTHON_ARGS}
                -I ${CMAKE_SOURCE_DIR}/src 
                -I ${CMAKE_BINARY_DIR}/src 
                -I ${CMAKE_CURRENT_SOURCE_DIR} 
                -I ${CMAKE_CURRENT_BINARY_DIR} 
                -I ${Python_NumPy_INCLUDE_DIRS}
                "${CMAKE_CURRENT_SOURCE_DIR}/${PYNE_CYTHON_MODULE}.pyx" 
                --output-file "${CMAKE_CURRENT_BINARY_DIR}/${PYNE_CYTHON_MODULE}.${PYNE_CYTHON_EXT}"
    )

    # Creating a Python module
    python_add_library(
        ${PYNE_CYTHON_MODULE} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${PYNE_CYTHON_MODULE}.${PYNE_CYTHON_EXT}"
        WITH_SOABI
    )

    # Add the pyne library as a dependency
    add_dependencies(${PYNE_CYTHON_MODULE} ${CMAKE_PROJECT_NAME})

    # Link the target with the pyne library
    target_link_libraries(${PYNE_CYTHON_MODULE} PRIVATE ${PYNE_CYTHON_LIBS})

    # Configure the RPATH
    set_target_properties(${PYNE_CYTHON_MODULE} PROPERTIES
      INSTALL_RPATH "$ORIGIN/core/lib"
    )

    # Install the target in the pyne directory
    install(TARGETS ${PYNE_CYTHON_MODULE} DESTINATION ../.)
endforeach()


if(BUILD_SPATIAL_SOLVER)
  include_directories(${F2PY_INCLUDE_DIR})
  # Common variables
  set(f2py_module_name "transport_spatial_methods")
  set(fortran_src_file "${PROJECT_SOURCE_DIR}/src/transport_spatial_methods/3d/main.f90")
  set(f2py_module_c "${f2py_module_name}module.c")

  # Generate sources
  add_custom_target(
    genpyf
    DEPENDS "${fortran_src_file}"
  )
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
    COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
            -I${PROJECT_BINARY_DIR}/src
            "${fortran_src_file}"
            -m "${f2py_module_name}"
            --lower # Important
    DEPENDS ${CMAKE_PROJECT_NAME}
  )

  # Set up target
  python_add_library( ${f2py_module_name} MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # Generated
    "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
    "${fortran_src_file}" # Fortran source(s)
  )

  # Depend on sources
  target_include_directories(${f2py_module_name} PUBLIC ${F2PY_INCLUDE_DIR})
  target_link_libraries(${f2py_module_name} PUBLIC Python::NumPy ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${CMAKE_PROJECT_NAME})
  
  # Configure the RPATH
  set_target_properties(${f2py_module_name} PROPERTIES
    INSTALL_RPATH "$ORIGIN/core/lib"
  )
  add_dependencies(${f2py_module_name} genpyf)
  install(TARGETS ${f2py_module_name} DESTINATION ../.)
endif()
