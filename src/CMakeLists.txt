set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS}
    PARENT_SCOPE)

pyne_download_platform()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(PYNE_SRCS
  "atomic_data.cpp"
  "data.cpp"
  "enrichment.cpp"
  "enrichment_cascade.cpp"
  "enrichment_symbolic.cpp"
  "jsoncpp.cpp"
  "jsoncustomwriter.cpp"
  "material.cpp"
  "material_library.cpp"
  "nucname.cpp"
  "particle.cpp"
  "rxname.cpp"
  "tally.cpp"
  "transmuters.cpp"
  "utils.cpp"
  "endftod.f90"
  )
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


#check operating system
IF(NOT ${PYNE_PLATFORM} STREQUAL "__LINUX__")
  set(BUILD_SPATIAL_SOLVER false)
ENDIF()

# setup source files
IF(BUILD_SPATIAL_SOLVER)
  set(PYNE_SRCS "${PYNE_SRCS}" "${TRANSPORT_SPATIAL_METHODS_SRCS}")
ENDIF()

if(MOAB_FOUND)
  set(PYNE_SRCS "${PYNE_SRCS}" "${MOAB_SRCS}")
endif(MOAB_FOUND)

# Fast compile decay and cram
fast_compile(decay.cpp "-O0" "-O0" "-O0")
fast_compile(cram.c
  "-O0 -fcx-fortran-rules -fcx-limited-range -ftree-sra -ftree-ter -fexpensive-optimizations"
  "-O0 -ffast-math"
  "-O0")

# compile and link library
add_library(pyne ${PYNE_SRCS})
if("${LIBS_HDF5}" STREQUAL "")
  target_link_libraries(pyne hdf5)
else()
  target_link_libraries(pyne ${LIBS_HDF5})
endif()
IF(BUILD_SPATIAL_SOLVER)
    target_link_libraries(pyne ${LIBS_HDF5} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
ENDIF(BUILD_SPATIAL_SOLVER)
if(MOAB_FOUND)
    target_link_libraries(pyne MOAB)
endif(MOAB_FOUND)
if(DAGMC_FOUND)
    target_link_libraries(pyne dagmc MOAB)
endif(DAGMC_FOUND)

add_executable(alphad ${PROJECT_SOURCE_DIR}/src/ensdf_processing/ALPHAD/alphad.f
                      ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(delta ${PROJECT_SOURCE_DIR}/src/ensdf_processing/DELTA/delta.for)
add_executable(gtol ${PROJECT_SOURCE_DIR}/src/ensdf_processing/GTOL/gtol.f
                    ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(bldhst ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/bldhst.for)
add_executable(hsicc ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/hsicc.for
                     ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(hsmrg ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/hsmrg.for
                     ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(seqhst ${PROJECT_SOURCE_DIR}/src/ensdf_processing/HSICC/seqhst.for
                      ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(logft ${PROJECT_SOURCE_DIR}/src/ensdf_processing/LOGFT/logft.for
                     ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)
add_executable(radd ${PROJECT_SOURCE_DIR}/src/ensdf_processing/RADD/RadD.for)
add_executable(ruler ${PROJECT_SOURCE_DIR}/src/ensdf_processing/RULER/ruler.f
                     ${PROJECT_SOURCE_DIR}/src/ensdf_processing/nsdflib95.f)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")


#
# Install header files
#
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB hpp_files "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
SET(h_files ${h_files} ${hpp_files})
INSTALL(FILES ${h_files} DESTINATION include/pyne COMPONENT headers)


#
# Install libpyne
#
INSTALL(TARGETS pyne LIBRARY DESTINATION lib COMPONENT pyne)
