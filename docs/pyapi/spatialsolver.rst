.. _pyne_spatialsolver:

==================================================
Spatialsolver Support -- :mod:`pyne.spatialsolver`
==================================================
Spatialsolver is a pyne module that contains seven neutron transport equation solvers.  Each
solver is its own unique nodal method.  The solvers included in this module are listed below.  The theory and methodology behind each
can be found in the :ref:`theory documentation <theorymanual_spatialsolver>`.

 #. **AHOTN-LN**: Arbitrarily higher order transport method of the nodal type 
    linear-nodal method
 #. **AHOTN-LL**:  Arbitrarily higher order transport method of the nodal type 
    linear-linear method
 #. **AHOTN-NEFD**: Arbitrarily higher order transport method of the nodal type
    that makes use of the unknown nodal flux moments (NEFD algorithm).
 #. **DGFEM-LD**: The Discontinuous Galerkin Finite Element Method (DGFEM) with
    a linear discontinuous (LD) approximation for angular flux.
 #. **DGFEM-DENSE**: The Discontinuous Galerkin Finite Element Method (DGFEM) 
    that use Lagrange polynomials to create a function space in each  dimension. 
 #. **DGFEM-LAGRANGE**:   The Discontinuous Galerkin Finite Element Method (DGFEM) 
    uses Lagrange polynomials to create a function space in each dimension.
 #. **SCTSTEP**: SCT Step algorithm uses a step approximation in all cells that
    are intersected by lines and planes of non-smoothness.

-----------------
Spatialsolver API
-----------------

.. automodule:: pyne.spatialsolver
    :members:

.. _source_file_formatting:

-----------------------------------
Source File Formatting
-----------------------------------

The spatial solver dictionary requires multiple input binary source files.  The required files are the following:
    1. XS file
    2. Source input file
    3. BC input file
    4. Quad file (optional)

Here is a brief description of how each should be (or is) formatted.


  (1.) XS file:
      The xs file contains information about the cross sections for materials used.  Each material should be assigned an ID,
      and the cross section data should be in the following format.  It should be saved as either an extensionless or .txt file.

        ! Cross section file

          ! Material # 1

          ! Group #1

          1.1          ! Total XS

          0.2         ! Scattering matrix

          ! Material 2

          . . .

          ! End Cross section file


  (2.) Source file:
      The source file is a file containing source information for each cell.  The formatting is dependant on the solver
      you select.

      For the AHOTN and DGFEM solvers, the source file should be formatted as following.  
      There should be ng * nx * ny * nz * lambda * lambda * lambda source entries present.
      We will refer to the index of each source value as (ng, nx, ny, nz, lambda_x, lambda_y, lambda_z).
      The source entries should be in the following order:

          (1,1,1,1,1,1,1)

          (1,1,1,1,1,1,2)

          (1,1,1,1,1,1,.)



          (1,1,1,1,1,2,1)

          (1,1,1,1,1,2,2)

          (1,1,1,1,1,2,.)

          (1,1,1,1,1,.,.)



          (1,1,1,1,2,1,1)

          (1,1,1,1,2,1,2)

          (1,1,1,1,2,1,.)

          (1,1,1,1,2,2,1)

          (1,1,1,1,2,2,2)

          (1,1,1,1,2,2,.)

          (1,1,1,1,2,.,.)

          (1,1,1,1,.,.,.)



          (1,1,1,2,1,1,1)

          (1,1,1,2,1,1,2)

          (1,1,1,2,1,1,.)

          (1,1,1,2,1,2,1)

          (1,1,1,2,1,2,2)

          (1,1,1,2,1,2,.)

          (1,1,1,2,1,.,.)

          (1,1,1,2,2,1,1)

          (1,1,1,2,2,1,2)

          (1,1,1,2,2,1,.)

          (1,1,1,2,2,2,1)

          (1,1,1,2,2,2,2)

          (1,1,1,2,2,2,.)

          (1,1,1,2,2,.,.)

          (1,1,1,.,.,.,.)


          ...

          (1,1,.,.,.,.,.)

          ...

          (1,.,.,.,.,.,.)

          ...

          (.,.,.,.,.,.,.)

      When being read in, they will be iterated over by the following loop:
	      DO g=1,ng
		       DO ix=1,nx
		          DO iy=1,ny
		             DO iz=1,nz
		                DO jx=0,lambda
		                   DO jy=0,lambda
		                      DO jz=0,lambda
		                         READ(12) s(jx,jy,jz,ix,iy,iz,g)


      For the SCT STEP solver, the source file should be formatted as following.  
      There should be ng * nx * ny * nz source entries present.
      We will refer to the index of each source value as (ng, nx, ny, nz).
      The source entries should be in the following order:

          (1,1,1,1)

          (1,1,1,2)

          (1,1,1,.)


          (1,1,2,1)

          (1,1,2,2)

          (1,1,2,.)

          (1,1,.,.)


          (1,2,1,1)

          (1,2,1,2)

          (1,2,1,.)

          (1,2,2,1)

          (1,2,2,2)

          (1,2,2,.)

          (1,2,.,.)

          (1,.,.,.)


          (2,1,1,1)

          (2,1,1,2)

          (2,1,1,.)

          (2,1,2,1)

          (2,1,2,2)

          (2,1,2,.)

          (2,1,.,.)

          (2,2,1,1)

          (2,2,1,2)

          (2,2,1,.)

          (2,2,2,1)

          (2,2,2,2)

          (2,2,2,.)

          (2,2,.,.)

          (.,.,.,.)

      When being read in, they will be iterated over by the following loop:
        DO g=1,ng
           DO ix=1,nx
              DO iy=1,ny
                 DO iz=1,nz
                     READ(12) s(ix,iy,iz,g,1,1,1)


  (3.) Boundary Condition file (only needed if one of the boundary conditions specified above is 2):
      The boundary condition file contains information about the incoming scalar flux on each face of each cell.
    
  (4.) Quadrature file (optional):
      If the quadrature type you selected was 2, a quadrature file is required for running the solver.  If the quadrature type is not 2, no quadrature file is necessary.

  (5.) Flux output file (output & optional):
      If a output file name was specified, the final flux will be printed to that file in the following format.  Note that all flux values will be printed as a fortran REAL, and the termination key will be 0.0d0 (to indicate the end of the flux info)

      AHOTN Solvers:  Unformatted file, with all mesh scalar flux values in the following order:

      NOTE ORDERING: flux(ng,nx,ny,nz,jx,jy,jz)

                                (0,0,0,0,0,0,0)

                                (0,0,0,0,0,0,1)

                                (0,0,0,0,0,1,0)

                                (0,0,0,0,0,1,1)

                                ...

                                (0,0,0,0,0,.,.)

                                ...

                                (0,0,0,0,.,.,.)

                                ...

                                (0,0,0,nz,.,.,.)

                                ...

                                (0,0,0,.,.,.,.)

                                ...

                                (0,0,ny,.,.,.,.)

                                ...

                                (0,0,.,.,.,.,.)

                                ...

                                (0,nz,.,.,.,.,.)

                                ...

                                (0,.,.,.,.,.,.)

                                ...

                                (ng,.,.,.,.,.,.)

                                ...


      DGFEM Solvers: Unformatted file, with all mesh scalar flux values in the following order:
             NOTE ORDERING: flux(ng,nx,ny,nz,spatial_order,spatial_order,spatial_order)

                                (0,0,0,0,0,0,0)

                                (0,0,0,0,0,0,1)

                                ...

                                (0,0,0,0,0,0,spatial_order)

                                (0,0,0,0,0,1,0)

                                (0,0,0,0,0,1,1)

                                ...

                                (0,0,0,0,0,1,spatial_order)

                                ...

                                (0,0,0,0,0,spatial_order,.)

                                ...

                                (0,0,0,0,spatial_order,.,.)

                                ...

                                (0,0,0,0,.,.,.)

                                ...

                                (0,0,0,nz,.,.,.)

                                ...

                                (0,0,0,.,.,.,.)

                                ...

                                (0,0,ny,.,.,.,.)

                                ...

                                (0,0,.,.,.,.,.)

                                ...

                                (0,nz,.,.,.,.,.)

                                ...

                                (0,.,.,.,.,.,.)

                                ...

                                (ng,.,.,.,.,.,.)

                                ...

      SCTSTEP Solvers: SCTSTEP Currently not supported.  Nothing will be printed to file (although file may still be created).


.. currentmodule:: pyne.spatialsolver

All functionality may be found in the ``spatialsolver`` package::

 from pyne import spatialsolver

.. _Spatialsolver: http://something.com/



