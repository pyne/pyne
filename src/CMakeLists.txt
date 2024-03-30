set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS}
    PARENT_SCOPE)

pyne_download_files()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(PYNE_SRCS
  "atomic_data.cpp"
  "data.cpp"
  "enrichment.cpp"
  "enrichment_cascade.cpp"
  "jsoncpp.cpp"
  "jsoncustomwriter.cpp"
  "material.cpp"
  "material_library.cpp"
  "nucname.cpp"
  "particle.cpp"
  "source.cpp"
  "rxname.cpp"
  "tally.cpp"
  "transmuters.cpp"
  "endftod.f90"
  "utils.cpp"
  )

# FIx: enrichment_symbolic.cpp doesn't compile on Windows
if(NOT WIN32)
  list(APPEND PYNE_SRCS "enrichment_symbolic.cpp")
endif()

set(TRANSPORT_SPATIAL_METHODS_SRCS
  "transport_spatial_methods/3d/ahotn_kernel_module.f90"
  "transport_spatial_methods/3d/angle.f90"
  "transport_spatial_methods/3d/check.f90"
  "transport_spatial_methods/3d/dgfem_kernel.f90"
  "transport_spatial_methods/3d/echo.f90"
  "transport_spatial_methods/3d/geompack.f90"
  "transport_spatial_methods/3d/igeompack.f90"
  "transport_spatial_methods/3d/inner.f90"
  "transport_spatial_methods/3d/invar.f90"
  "transport_spatial_methods/3d/output.f90"
  "transport_spatial_methods/3d/output_phi.f90"
  "transport_spatial_methods/3d/p.f90"
  "transport_spatial_methods/3d/precision_module.f90"
  "transport_spatial_methods/3d/read_inflow_ahotn.f90"
  "transport_spatial_methods/3d/read_inflow_dgfem.f90"
  "transport_spatial_methods/3d/read_inflow_sct_step.f90"
  "transport_spatial_methods/3d/readsrc.f90"
  "transport_spatial_methods/3d/readxs.f90"
  "transport_spatial_methods/3d/sct_module.f90"
  "transport_spatial_methods/3d/sct_step_kernel_module.f90"
  "transport_spatial_methods/3d/solvar.f90"
  "transport_spatial_methods/3d/solve.f90"
  "transport_spatial_methods/3d/sweep_ahotn_l.f90"
  "transport_spatial_methods/3d/sweep_ahotn_nefd.f90"
  "transport_spatial_methods/3d/sweep_dgfem.f90"
  "transport_spatial_methods/3d/sweep_sct_step.f90"
  "transport_spatial_methods/3d/timevar.f90"
  "transport_spatial_methods/3d/trackroutines.f90"
  "transport_spatial_methods/3d/trackstruct.f90"
  "transport_spatial_methods/3d/version.f90"
  #"transport_spatial_methods/3d/main.f90"
  )
set(MOAB_SRCS
  "measure.cpp"
  "source_sampling.cpp"
  )

# setup source files
if(BUILD_SPATIAL_SOLVER)
  set(PYNE_SRCS "${PYNE_SRCS}" "${TRANSPORT_SPATIAL_METHODS_SRCS}")
endif()

if(MOAB_FOUND)
  set(PYNE_SRCS "${PYNE_SRCS}" "${MOAB_SRCS}")
endif(MOAB_FOUND)

# Fast compile decay and cram
fast_compile(decay.cpp "-O0" "-O0" "-O0")
fast_compile(cram.c
  "-O0 -fcx-fortran-rules -fcx-limited-range -ftree-sra -ftree-ter -fexpensive-optimizations"
  "-O0 -ffast-math"
  "-O0"
  )

# compile and link library
add_library(${CMAKE_PROJECT_NAME} ${PYNE_SRCS})
if(NOT HDF5_FOUND)
  add_dependencies(${CMAKE_PROJECT_NAME} hdf5)
endif()
if("${LIBS_HDF5}" STREQUAL "")
  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE hdf5)
else()
  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${LIBS_HDF5})
endif()
if(BUILD_SPATIAL_SOLVER)
    if (NOT LAPACK_FOUND)
    add_dependencies(${CMAKE_PROJECT_NAME} lapack)
    endif()
    target_link_libraries(${CMAKE_PROJECT_NAME} 
      PRIVATE ${LAPACK_LIBRARIES}
      )
endif()
if(MOAB_FOUND)
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE MOAB)
endif(MOAB_FOUND)
if(DAGMC_FOUND)
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE dagmc MOAB)
endif(DAGMC_FOUND)

if(ENABLE_ENSDF_PROCESSING)
  add_executable(alphad 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/ALPHAD/alphad.f
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS alphad RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(delta 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/DELTA/delta.for
    )
  install(TARGETS delta RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(gtol 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/GTOL/gtol.f
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS gtol RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(bldhst 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/bldhst.for
    )
  install(TARGETS bldhst RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(hsicc 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/hsicc.for
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS hsicc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(hsmrg 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/hsmrg.for
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS hsmrg RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(seqhst 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/seqhst.for
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS seqhst RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(logft 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/LOGFT/logft.for
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS logft RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(radd 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/RADD/RadD.for
    )
  install(TARGETS radd RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(ruler 
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/RULER/ruler.f
    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f
    )
  install(TARGETS ruler RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")


#
# Install header files
#
file(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB hpp_files "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set(h_files ${h_files} ${hpp_files})
install(FILES ${h_files} 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pyne 
  COMPONENT headers
  )


#
# Install libpyne
#
install(TARGETS ${CMAKE_PROJECT_NAME} LIBRARY 
  DESTINATION ${CMAKE_INSTALL_LIBDIR} 
  COMPONENT ${CMAKE_PROJECT_NAME}
)

if(RPATH)
  # Configure the RPATH
  set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES
    INSTALL_RPATH ${RPATH}
  )
endif()
