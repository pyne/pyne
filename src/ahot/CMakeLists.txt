set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(AHOT_SRCS
    "spatial_solvers_3d/source/ahotn_kernel_module.f90"
    "spatial_solvers_3d/source/readsrc.f90"
    "spatial_solvers_3d/source/angle.f90"
    "spatial_solvers_3d/source/readxs.f90"
    "spatial_solvers_3d/source/check.f90"
    "spatial_solvers_3d/source/read_inflow_ahotn.f90"
    "spatial_solvers_3d/source/dgfem_kernel.f90"
    "spatial_solvers_3d/source/read_inflow_dgfem.f90"
    "spatial_solvers_3d/source/echo.f90"      
    "spatial_solvers_3d/source/solvar.f90"
    "spatial_solvers_3d/source/inner.f90"
    "spatial_solvers_3d/source/solve.f90"
    "spatial_solvers_3d/source/invar.f90"
    "spatial_solvers_3d/source/sweep_ahotn_l.f90"
    #"spatial_solvers_3d/source/main.f90"
    "spatial_solvers_3d/source/sweep_ahotn_nefd.f90"
    "spatial_solvers_3d/source/output.f90"
    "spatial_solvers_3d/source/sweep_dgfem.f90"
    "spatial_solvers_3d/source/output_phi.f90"
    "spatial_solvers_3d/source/timevar.f90"
    "spatial_solvers_3d/source/p.f90"
    "spatial_solvers_3d/source/version.f90"
    "spatial_solvers_3d/source/readmt.f90"
    )

# compile and link library
add_library(ahot ${AHOT_SRCS})
target_link_libraries(ahot blas lapack)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
