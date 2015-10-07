set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
 
# setup source files
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

#check operating system
IF(NOT ${PYNE_PLATFORM} STREQUAL "__LINUX__")
    set(BUILD_SPATIAL_SOLVER false)
ENDIF()

# setup source files
IF(BUILD_SPATIAL_SOLVER)
    set(PYNE_SRCS
        "data.cpp"
        "decay.cpp"
        "enrichment.cpp"
        "enrichment_cascade.cpp"
        "enrichment_symbolic.cpp"
        "jsoncpp.cpp"
        "jsoncustomwriter.cpp"
        "material.cpp"
        "nucname.cpp"
        "particle.cpp"
        "rxname.cpp"
        "tally.cpp"
        "utils.cpp"
        "endftod.f90"
        "${TRANSPORT_SPATIAL_METHODS_SRCS}"
        )
ELSE()
    set(PYNE_SRCS
        "data.cpp"
        "decay.cpp"
        "enrichment.cpp"
        "enrichment_cascade.cpp"
        "enrichment_symbolic.cpp"
        "jsoncpp.cpp"
        "jsoncustomwriter.cpp"
        "material.cpp"
        "nucname.cpp"
        "particle.cpp"
        "rxname.cpp"
        "tally.cpp"
        "utils.cpp"
        "endftod.f90"
        )
ENDIF()

if (MOAB_FOUND)
    set(PYNE_SRCS "${PYNE_SRCS}"
        "dagmc_bridge.cpp"
        "measure.cpp"
        "source_sampling.cpp"
        )
endif (MOAB_FOUND)

# compile and link library
add_library(pyne ${PYNE_SRCS})
target_link_libraries(pyne hdf5)
IF(BUILD_SPATIAL_SOLVER)
    target_link_libraries(pyne hdf5 blas lapack)
ENDIF(BUILD_SPATIAL_SOLVER)
if (MOAB_FOUND)
    target_link_libraries(pyne dagmc MOAB)
endif (MOAB_FOUND)

set_source_files_properties(decay.cpp PROPERTIES COMPILE_FLAGS 
    "-O0")

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
