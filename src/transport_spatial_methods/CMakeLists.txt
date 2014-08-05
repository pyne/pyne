set(PYNE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" ${PYNE_INCLUDE_DIRS} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# setup source files
set(AHOT_SRCS
    "pyne/src/transport_spatial_methods /3d/ahotn_kernel_module.f90"
    "/3d/readsrc.f90"
    "3d/angle.f90"
    "3d/readxs.f90"
    "3d/check.f90"
    "3d/read_inflow_ahotn.f90"
    "3d/dgfem_kernel.f90"
    "3d/read_inflow_dgfem.f90"
    "3d/echo.f90"      
    "3d/solvar.f90"
    "3d/inner.f90"
    "3d/solve.f90"
    "3d/invar.f90"
    "3d/sweep_ahotn_l.f90"
    #"3d/main.f90"
    "3d/sweep_ahotn_nefd.f90"
    "3d/output.f90"
    "3d/sweep_dgfem.f90"
    "3d/output_phi.f90"
    "3d/timevar.f90"
    "3d/p.f90"
    "3d/version.f90"
    "3d/readmt.f90"
    )

# compile and link library
add_library(ahot ${AHOT_SRCS})
target_link_libraries(ahot blas lapack)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- Include paths for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
