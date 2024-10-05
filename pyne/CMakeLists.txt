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
if(ENABLE_MOAB)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} source_sampling)
endif()

if(ENABLE_DAGMC)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} _dagmc)
endif()

if(NOT WIN32)
  set(PYNE_CYTHON_MODULES ${PYNE_CYTHON_MODULES} cram)
endif()

message(STATUS "Modules to Cythonize: ${PYNE_CYTHON_MODULES}")

foreach(PYNE_CYTHON_MODULE ${PYNE_CYTHON_MODULES})

  # _dagmc module needs Eigen
  if(PYNE_CYTHON_MODULE STREQUAL "_dagmc")
    set(PYNE_DAGMC_INCLUDES "-I${EIGEN3_INCLUDE_DIRS}")
  else()
    set(PYNE_DAGMC_INCLUDES)
  endif()

  # Cram is actually a C, not C++ file
  if(NOT "${PYNE_CYTHON_MODULE}" STREQUAL "cram")
      set(PYNE_CYTHON_ARGS --cplus)
      set(PYNE_CYTHON_EXT cxx)
  else()
      set(PYNE_CYTHON_ARGS)
      set(PYNE_CYTHON_EXT c)
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
            ${PYNE_DAGMC_INCLUDES}
            "${CMAKE_CURRENT_SOURCE_DIR}/${PYNE_CYTHON_MODULE}.pyx" 
            --output-file "${CMAKE_CURRENT_BINARY_DIR}/${PYNE_CYTHON_MODULE}.${PYNE_CYTHON_EXT}"
    )

  # Creating a Python module
  python_add_library(
    ${PYNE_CYTHON_MODULE} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${PYNE_CYTHON_MODULE}.${PYNE_CYTHON_EXT}"
    WITH_SOABI
    )

  # pyne needs to be built before the cython module
  add_dependencies(${PYNE_CYTHON_MODULE} ${CMAKE_PROJECT_NAME})

  # Link the target with the pyne library
  target_link_libraries(${PYNE_CYTHON_MODULE} PRIVATE ${CMAKE_PROJECT_NAME})

  # Configure the RPATH
  if(PYNE_PYTHON_MODULE_RPATH)
    set_target_properties(${PYNE_CYTHON_MODULE} PROPERTIES
      INSTALL_RPATH ${PYNE_PYTHON_MODULE_RPATH}
      )
  endif()

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
            "${fortran_src_file}"
            -m ${f2py_module_name}
            --lower # Important
    DEPENDS ${fortran_src_file} ${CMAKE_PROJECT_NAME}
    )

  # Set up target
  python_add_library( ${f2py_module_name} MODULE WITH_SOABI
    "${f2py_module_c}" # Generated
    "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
    "${fortran_src_file}" # Fortran source(s)
    )

  # Depend on sources
  target_include_directories(${f2py_module_name} PUBLIC ${F2PY_INCLUDE_DIR})
  target_link_libraries(${f2py_module_name} PUBLIC Python::NumPy ${CMAKE_PROJECT_NAME})
  
  # Configure the RPATH
  if(PYNE_PYTHON_MODULE_RPATH)
    set_target_properties(${f2py_module_name} PROPERTIES
      INSTALL_RPATH ${PYNE_PYTHON_MODULE_RPATH}
      )
  endif()
  add_dependencies(${f2py_module_name} genpyf ${CMAKE_PROJECT_NAME})
  install(TARGETS ${f2py_module_name} DESTINATION ../.)
endif()
