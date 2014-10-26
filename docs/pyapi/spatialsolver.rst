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

-----------------------------------
Input Dictionary Entries
-----------------------------------
As these are complicated solvers, they require a large amount of input data 
supplied by the user.  This information needs to be entered as a Python 
dictionary. Of the many key-pair values listed below, most are required, but 
some are optional.  The optional entries will be overridden by default values 
if not present/specified.

**Entry: Solver type (AHOTN, DGFEM or SCTSTEP)**
::
  key: "solver"
  type: String
  ex: "AHOTN"
  default: no default
  Note:  The three supported "solver"s are the following:
    1.  "AHOTN"
    2.  "DGFEM"
    3.  "SCTSTEP"
 
**Entry: Spatial expansion order**
::
  key: "solver_type"
  type: String
  ex: "LN"
  default: No default
  Note: This entry is dependent on the "solver" entry.
    For AHOTN solver, there exist the "LN", "LL", and "NEFD" solver types
    For the DGFEM solvers, there exist the "LD", "DENSE", and "LAGRANGE" solver types
    For the SCTSTEP solver, the "solver_type" key is not used.

**Entry: Spatial expansion order (lambda; ahot spatial order, 0, 1, or 2)**
::
  key: "spatial_order"
  type: Integer
  ex: 0
  default: 1

  The Spatial expansion order is the expansion order of the spatial moment.

**Entry: Angular quadrature order**
::
  key: "quadrature_order"
  type: Integer
  ex: 4
  default: 4

  The angular quadrature order is the number of angles to be used per octant.  
  For N sets of angles, there will be (N * (N + 2) / 8) ordinates per octant. 
  The quadrature order may only be an even number!

**Entry: Quadrature type:**
::
  key: "quadrature_type"
  type: Integer
  ex: 1
  default: 1

  The quadrature type is the type of quadrature scheme the code uses.  
  The possibilities are:
    1 - TWOTRAN
    2 - EQN
    3 - Read-in

**Entry: Number of spatial nodes in x, y, and z directions (nx/ny/nz)**
::
  key: "nodes_xyz"
  type: Integer array
  ex: [4, 4, 4]
  default: No default
    
**Entry: Number of energy groups (ng)**
::
 key: "num_groups"
 type: Integer
 ex: 1
 default: No default

**Entry: Number of materials (nm)**
::
 key: "num_materials"
 type: Integer
 ex: 1
 default: No default

**Entry: x-size of cells (dx)**
::
 key: "x_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: y-size of cells (dy)**
::
 key: "y_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: z-size of cells (dz)**
::
 key: "z_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: x start and end boundary conditions**
::
 key: "x_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default
 'x_boundary_conditions' are the x boundary conditions for each face of the cubic mesh. 
 The entries are:
    x[0] = x starting bc; the left
    x[1] = x ending bc; the right
    The following are supported boundary conditions: 
      0 - vacuum
      1 - reflective
      2 - fixed inflow

**Entry: y start and end boundary conditions**
::
 key: "y_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default
 'y_boundary_conditions' are the y boundary conditions for each face of the cubic mesh. 
 The entries are:
    y[0] = y starting bc; the front
    y[1] = y ending bc; the back
    The following are supported boundary conditions: 
      0 - vacuum
      1 - reflective
      2 - fixed inflow

**Entry: z start and end boundary conditions**
::
 key: "z_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default
 'z_boundary_conditions' are the z boundary conditions for each face of the cubic mesh. 
 The entries are:
    z[0] = z starting bc; the bottom
    z[1] = z ending bc; the top
    The following are supported boundary conditions: 
      0 - vacuum
      1 - reflective
      2 - fixed inflow

**Entry: Material info**
::
 key: "material_id"
 type: Integer 3 dimensional array
 ex: [ [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] 
       [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] 
       [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] ]
 default: No default
 note: Dimensions must match cells such that there is one material number
       in each spatial cell. The cells are ordered as x, y, z.

    
**Entry: "quadrature_file" [optional; only needed for quadrature type 2]**
::
 type: String
 ex: 'quad_file'
 default: No default  
 note: See input file formatting notes in the Source File Formatting section.

**Entry: cross section info file name**
::
 key: "xs_file"
 type: String
 default: 'xs_file'
 note: See input file formatting notes in the Source File Formatting section.

**Entry: source file name**
::
  key: "source_input_file"
  type: String
  default: 'src.dat'
  note: See input file formatting notes in the Source File Formatting section.

**Entry: boundary condition file name [optional]**
::
 key: "bc_input_file"
 type: String
 default: No default
 note: See input file formatting notes in the Source File Formatting section.

**Entry: output file name [optional]**
::
 key: "flux_output_file"
 type: String
 default: 'flux.out'
 note: See input file formatting notes in the Source File Formatting section.

.. Josh: please to add file format / contents description

**Entry: Convergence Criterion**
::
 key: "convergence_criterion"
 type: float
 ex: 1.e-12
 default: 1.e-12
 To terminate the calculation, the convergence criterion 
 The convergence criterion is the cell-wise relative difference (df) in the flux value between solution iterates. 
 Once the flux in each cell satisfies this criterion, the solver will calculate the final flux.
 df is calculated as the following:  
 f = current flux value
 ct = convergence tolerance (value for this key, "converge_tolerance")
 f1 = flux value from the previous iteration
 If (f - f1) > ct:
   df = absolute(f - f1) / f1
 Else
   df = absolute(f - f1)

**Entry: Tolerance**
::
 key: "converge_tolerance"
 type: float
 ex: 1.e-10
 default: 1.e-10
 Converge tolerance is the tolerance that determines how the relative difference between flux iterations (df)
 will be calculated.  This effects when the solver will stop its mesh sweeps.  See the convergence criterion key
 value pair for more information.

**Entry: Maximum Iterations**
::
 key: "max_iterations"
 type: int
 ex: 10000
 default: 10000
 note: If this number of iterations is reached before the convergence criterion
       is satisfied, the calculation will terminate and report the current flux
       estimate.

**Entry: Moments Converged**
::
 key: "moments_converged" ??
 type: int
 ex: 0
 default: 0

.. Josh: Isn't this taken out now?

**Entry: Flag for computing, presenting quadrant fluxes**
::
 key: "quad_flux_print_flag"
 type: Integer
 ex: 0
 default: 0

.. Josh: Isn't this taken out now?

**Entry: Flag for printing matrix file**
::
 key: "matrix_print_flag"
 type: Integer
 ex: 0
 default: 0

.. Josh: Isn't this taken out now?

**Entry: ITM solution flag**
::
 key: "itm_direct_solution_flag" [only relevant if itm is selected]
 type: Integer
 ex: 0
 default: default of 0 if itm solution method

.. Josh: Isn't this taken out now?

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

**Entry: Solver success code**
::
  key:  "success"
  type: Integer
  format: 1 means yes, the solve succeeded.  0 means it failed.

**Entry: Raw system time of solver start**
::
  time_start provides you with the system time when the solver began running.
  key:  "time_start"
  type: double
  format: system time

**Entry: Total run time**
::
  total_time is the total time the solver took to solve.
  key:  "total_time"
  type: double
  format: system time

**Entry: Total print time**
::
  print_time is the total time the solver took to print results.
  key:  "print_time"
  type: double
  format: system time

**Entry: Error Message**
::
  If the solver fails, error_msg is a string describing why the solver failed.
  key:  "error_msg"
  type: String

.. _source_file_formatting:

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



