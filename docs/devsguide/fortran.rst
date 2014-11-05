.. _devsguide_fortran:

=============================
Integrating Fortran with PyNE
=============================

Using the f2py tool distributed under numpy, PyNE can now support development
with fortran.  The Spatial Solver component of PyNE uses this in the manner
outlined below:

----------------------------------
Overview
----------------------------------

F2py is a tool that can generate both python and c interfaces for fortran codes.  Because PyNE uses the CMAKE build system, all f2py interface building is done via CMAKE.  When the initial build process occurs, CMAKE:
  1. Creates fortran objects for each of the required source files (as specified in the CMAKELISTS.TXT files)
  2. Wraps them with f2py, creating a shared object
  3. Links that shared object with the PyNE python interface

************************
Where Fortran Files Live
************************
* All fortran files currently live in the pyne/src folder.  For Spatial Solver, all files live in pyne/src/transport_spatial_methods.

************************
Spatial Solver Specifics
************************
* The Spatial Solver module has a fortran API, but it is not accessible through PyNE. The only
  accessible API is the wrapped and linked python interface that uses the fortran codes through
  f2py.  Pyne currently has no frameworks for a direct fortran API to be set up with.

*****************************
Steps for wrapping with CMAKE
*****************************
* All fortran files, except for the one being wrapped, must be declared in a CMakeLists.txt file located in the fortran source directory???
* Following the file declaration, the file grouping must be linked.  The following is the linking for the spatial solver module:
 
   # compile and link library

   add_library(ahot ${AHOT_SRCS})

   target_link_libraries(ahot blas lapack)
* F2py wrapping must then be done in pyne/pyne/CMakeLists.txt.  The spatial solver wrapping (and linking) is done as following:

add_custom_target(transport_spatial_methods ALL DEPENDS transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX})

add_custom_command(OUTPUT transport_spatial_methods${CMAKE_SHARED_LIBRARY_SUFFIX}

     COMMAND f2py -m transport_spatial_methods -I${PROJECT_BINARY_DIR}/src -L${PROJECT_BINARY_DIR}/src

                 --f90flags="-fdefault-real-8"

                 -c ${PROJECT_SOURCE_DIR}/src/transport_spatial_methods/3d/main.f90

                 -lpyne

                 #-lpyne --debug-capi

     DEPENDS pyne

     )

