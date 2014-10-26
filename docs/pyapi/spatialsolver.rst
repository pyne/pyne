.. _pyne_spatialsolver:

==================================================
Spatialsolver Support -- :mod:`pyne.spatialsolver`
==================================================
Spatialsolver is a pyne module that contains seven neutron transport equation solvers.  Each
solver is its own unique nodal method.  The solvers included in this module are listed below.  The theory and methodology behind each
can be found in the pyne theory documentation.

 #. **AHOTN-LN**: Arbitrarily higher order transport method of the nodal type linear-nodal method
 #. **AHOTN-LL**:  Arbitrarily higher order transport method of the nodal type linear-linear method
 #. **AHOTN-NEFD**: Arbitrarily higher order transport method of the nodal type that makes use of the
    unknown nodal flux moments (NEFD algorithm).
 #. **DGFEM-LD**: The Discontinuous Galerkin Finite Element Method (DGFEM) with a linear discontinuous (LD)
    approximation for angular flux. (SEE PAGE 27 of thesis)
 #. **DGFEM-DENSE**: The Discontinuous Galerkin Finite Element Method (DGFEM) that use ??dense??? lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
 #. **DGFEM-LAGRANGE**:   The Discontinuous Galerkin Finite Element Method (DGFEM) that use lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
 #. **SCTSTEP**: SCT Step algorithm similar to Duo's SCT algorithm implemented in three dimensional Cartesian
    geometry.

-----------------------------------
Input Dictionary Entries
-----------------------------------
As these are complicated solvers, they require a large amount of input data supplied by the user.  The
format we choose to take all this information in by is with a python dictionary.   Of the many key-pair values listed below, most are required, but some are optional.  The optional entries will be overridden by default values if not present/not specified.

**Entry: Solver type (AHOTN, DGFEM or SCTSTEP)**::

  key: "solver"
  type: String
  ex: "AHOTN"
  default: no default
  Note:  The three supported "solver"s are the following:
    1.  "AHOTN"
    2.  "DGFEM"
    3.  "SCTSTEP"
 
**Entry: Spatial expansion order**::

  key: "solver_type"
  type: String
  ex: "LN"
  default: No default
  Note: This entry is dependent on the "solver" entry.
    For AHOTN solver, there exist the LN, LL and NEFD solver types
    For the DGFEM solvers, there exist the LD, DENSE and LAGRANGE solver types
    For the SCTSTEP solver, no solver_type key is required.  The key can be set to something or be left empty.

**Entry: Spatial expansion order (lambda; ahot spatial order, 0, 1, or 2)**::

  key: "spatial_order"
  type: Integer
  ex: 0
  default: 1

  The Spatial expansion order is the expansion order of the spatial moment.

**Entry: Quadrature order**::

  key: "quadrature_order"
  type: Integer
  ex: 4
  default: 4

  The quadrature order is the number of angles to be used per octet.  For N sets of angles, there will
  be (N * (N + 2) / 8) ordinates per octet. The quadrature order may only be an even number!

**Entry: Quadrature type:**::

  key: "quadrature_type"
  type: Integer
  ex: 1
  default: 1

  The quadrature type is the type of quadrature scheme the code should use.  The possibilities are as follow:
    1 - TWOTRAN
    2 - EQN
    3 - Read-in

**Entry: Number of Nodes in x, y, and z directions (nx/ny/nz)**::

  key: "nodes_xyz"
  type: Integer array
  ex: [4, 4, 4]
  default: No default
    
**Entry: Number of groups (ng)**::

 key: "num_groups"
 type: Integer
 ex: 1
 default: No default

**Entry: Number of Materials (nm)**::

 key: "num_materials"
 type: Integer
 ex: 1
 default: No default

**Entry: x-size of cells (dx)**::

 key: "x_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: y-size of cells (dy)**::

 key: "y_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: z-size of cells (dz)**::

 key: "z_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: x start and end boundary conditions**::
 key: "x_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default
 'x_boundary_conditions', 'y_boundary_conditions', and 'z_boundary_conditions' are the boundary conditions for each face of the cubic mesh. The entries are as following:
    x[0] = xsbc
    x[1] = xebc
    y[0] = ysbc
    y[1] = yebc
    y[0] = zsbc
    y[1] = zebc
    The following are supported boundary conditions: 
      0 - vacuum
      1 - reflective
      2 - fixed inflow

**Entry: y start and end boundary conditions**::

 key: "y_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default

**Entry: z start and end boundary conditions**::

 key: "z_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default

**Entry: Material info**::

 key: "material_id"
 type: Integer 3 Dimensional Array
 ex: [ [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] ]
 default: No default
 Note:  Dimensions must match cell spacing and ordering

_Note: we need to give directions about the ordering. RS
    
**Entry: "quadrature_file" [optional; only needed for quadrature type 2]**::

 type: string
 ex: 'quad_file'
 default: No default  
 note: See input file formatting notes in the Source File Formatting section.

**Entry: cross section info file name**::

 key: "xs_file"
 type: string
 default: 'xs_file'
 note: See input file formatting notes in the Source File Formatting section.


**Entry: source file name**::

  key: "source_input_file"roman atwood
  type: string
  default: 'src.dat'
  note: See input file formatting notes in the Source File Formatting section.

**Entry: boundary condition file name [optional]**::

 key: "bc_input_file"
 type: string
 default: No default
 note: See input file formatting notes in the Source File Formatting section.

**Entry: output file name [optional]**::

 key: "flux_output_file"
 type: string
 default: 'flux.out'

_note: need to add file format / contents description

**Entry: Convergence Criterion**::

 key: "convergence_criterion"
 type: float
 ex: 1.e-12
 default: 1.e-12
 The convergence criterion is the maximum allowed relative difference (df) in the flux value from one sweep to the next.  The more times the solver
 sweeps the more definite the solution will be.  As soon as the convergence criterion is met, the solver will stop sweeping and calculate
 the final flux for that cell. ??Correct??
 df is calculated as the following:  
 f = current flux value
 ct = convergence tolerance (key pair value "converge_tolerance"0
 f1 = flux value from one iteration prior
 If f - f1 > ct:
   df = absolute(f - f1) / f1
 Else
   df = absolute(f - f1)

**Entry: Maximum Iterations**::

 key: "max_iterations"
 type: int
 ex: 10000
 default: 10000

**Entry: Moments Converged**::

 key: "moments_converged" ??
 type: int
 ex: 0
 default: 0

_We need to provide a clear explanation of what this means. RS

**Entry: Tolerance**::

 key: "converge_tolerance"
 type: float
 ex: 1.e-10
 default: 1.e-10
 Converge tolerance is the tolerance that determines how the relative difference between flux iterations (df)
 will be calculated.  This effects when the solver will stop its mesh sweeps.  See the convergence criterion key
 value pair for more information.

**Entry: Flag for computing, presenting quadrant fluxes**::

 key: "quad_flux_print_flag"
 type: Integer
 ex: 0
 default: 0

**Entry: Flag for printing matrix file**::

 key: "matrix_print_flag"
 type: Integer
 ex: 0
 default: 0

**Entry: ITM solution flag**::

 key: "itm_direct_solution_flag" [only relevant if itm is selected]
 type: Integer
 ex: 0
 default: default of 0 if itm solution method

_We need to double check the meaning of this one. RS
_I have not added this yet.  As soon as we verify the meaning I can add it. JH
_great RS

-----------------------------------
Output Dictionary Entries
-----------------------------------
When ran, the solvers return a dictionary of useful solution data.  It contains the following key-pair entries:

**Entry: Flux output array**::
  key:  "flux"
  type: Double Array of 3 dimensions
  format: Flux output array is in following format.  Each cell in the array has a scalar flux, the integral of the angular flux over all angles in that cell.   The first index refers to the plane on the z axis, beginning at 0 with the lowest plane, and moving upwards to the highest plane on the mesh.  The second index is the row on the z plane, and the third index is the cell in the row.
  format examples: If you had a mesh with 4 by 4 by 4 cells extending in the x, y and z directions, then to get the following flux values, you would use the following index's:

  (1.) Scalar flux across top of cell 1,1,1:  flux_array[1][1][1]
       Geometric location of this cell:
          Plane: Bottom of cube
          Row: First y row (j) of cells
          Cell: First cell in x direction
  (2.) Scalar flux across top of cell 1,1,2:  flux_array[1][1][2]
       Geometric location of this cell:
          Plane: Bottom of cube
          Row: First y row (j) of cells
          Cell: Second cell in x direction
  (3.) Scalar flux across top of cell 1,2,1:  flux_array[1][2][1]
       Geometric location of this cell:
          Plane: Bottom of cube
          Row: Second y row (j) of cells
          Cell: First cell in x direction
  (4.) Scalar flux across top of cell 2,1,1:  flux_array[2][1][1]
       Geometric location of this cell:
          Plane: Top of one cell up from bottom of cube
          Row: First y row (j) of cells
          Cell: First cell in x direction

**Entry: Solver success code**::
  key:  "success"
  type: Integer
  format: 1 means yes, the solve succeeded.  0 means it failed.

**Entry: Raw system time of solver start **::
  time_start provides you with the system time when the solver began running.
  key:  "time_start"
  type: double
  format: system time

**Entry: Total run time **::
  total_time is the total time the solver took to solve.
  key:  "total_time"
  type: double
  format: system time

**Entry: Total print time **::
  print_time is the total time the solver took to print results.
  key:  "print_time"
  type: double
  format: system time

**Entry: Error Message **::
  If the solver fails, error_msg is a string describing why the solver failed.
  key:  "error_msg"
  type: String

-----------------------------------
Source File Formatting
-----------------------------------

The spatial solver dictionary requires multiple input binary source files.  The required files are the following:
    1. XS file
    2. Source input file
    3. BC input file
    4. Quad file (optional)

Here is a brief description of how each should be formatted.


  (1.) XS file:
      The xs file contains information about the cross sections for materials used.  Each material should be assigned an ID,
      and the cross section data should be in the following format.  It should be saved as either an extensionless or .txt file.

        ! Cross section file
          ! Material # 1
          ! Group #1
          1.1          ! Total XS
          0.2         ! Scattering matrix
          ! Material 2
          ...
        ! End Cross section file

  (2.) Source file:
      The source file is a file containing source information for each cell ????.  The formatting is dependant on the solver
      you select.

      For the AHOTN and DGFEM solvers, the source file should be formatted as following.  
      There should be ng * nx * ny * nz * lambda * lambda * lambda ?source? entries present.
      We will refer to the index of each ?source? value as (ng, nx, ny, nz, lambda_x, lambda_y, lambda_z).
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
      There should be ng * nx * ny * nz ?source? entries present.
      We will refer to the index of each ?source? value as (ng, nx, ny, nz).
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

! This line is supposed to be quadrature file name (qdfile) if you need one (type 2)            

.. currentmodule:: pyne.spatialsolver

All functionality may be found in the ``spatialsolver`` package::

 from pyne import spatialsolver

Spatialsolver API
-----------

.. automodule:: pyne.spatialsolver
    :members:

.. _Spatialsolver: http://something.com/



