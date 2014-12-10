set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(AHOT_SRCS
    "ahot/spatial_solvers_3d/source/ahotn_kernel_module.f90"
    "ahot/spatial_solvers_3d/source/readsrc.f90"
    "ahot/spatial_solvers_3d/source/angle.f90"
    "ahot/spatial_solvers_3d/source/readxs.f90"
    "ahot/spatial_solvers_3d/source/check.f90"
    "ahot/spatial_solvers_3d/source/read_inflow_ahotn.f90"
    "ahot/spatial_solvers_3d/source/dgfem_kernel.f90"
    "ahot/spatial_solvers_3d/source/read_inflow_dgfem.f90"
    "ahot/spatial_solvers_3d/source/echo.f90"
    "ahot/spatial_solvers_3d/source/solvar.f90"
    "ahot/spatial_solvers_3d/source/inner.f90"
    "ahot/spatial_solvers_3d/source/solve.f90"
    "ahot/spatial_solvers_3d/source/invar.f90"
    "ahot/spatial_solvers_3d/source/sweep_ahotn_l.f90"
    #"ahot/spatial_solvers_3d/source/main.f90"
    "ahot/spatial_solvers_3d/source/sweep_ahotn_nefd.f90"
    "ahot/spatial_solvers_3d/source/output.f90"
    "ahot/spatial_solvers_3d/source/sweep_dgfem.f90"
    "ahot/spatial_solvers_3d/source/output_phi.f90"
    "ahot/spatial_solvers_3d/source/timevar.f90"
    "ahot/spatial_solvers_3d/source/p.f90"
    "ahot/spatial_solvers_3d/source/version.f90"
    "ahot/spatial_solvers_3d/source/readmt.f90"
    )

# setup source files
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
    "${AHOT_SRCS}"
    )

if (MOAB_FOUND)
    set(PYNE_SRCS "${PYNE_SRCS}"
        "dagmc_bridge.cpp"
        "measure.cpp"
        "source_sampling.cpp"
        )
endif (MOAB_FOUND)

# compile and link library
add_library(pyne ${PYNE_SRCS})
target_link_libraries(pyne hdf5 blas lapack)
if (MOAB_FOUND)
    target_link_libraries(pyne dagmc MOAB)
endif (MOAB_FOUND)

set_source_files_properties(decay.cpp PROPERTIES COMPILE_FLAGS 
    "-O0 -time")

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
