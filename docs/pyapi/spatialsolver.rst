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

**Entry: Method (meth): 0/1 => AHOT-N-SI/AHOT-N-ITM**::

  key: "spatial_method"   
  type: Integer
  ex: 0
  default: 0

**Entry: Quadrature order (can only be an even number)**::

  key: "quadrature_order"
  type: Integer
  ex: 4
  default: 4

**Entry: Qudrature type:**::

  key: "quadrature_type"
  type: Integer
  ex: 1
  default: 1

**Entry: Number of Nodes in x, y, and z directions (nx/ny/nz)**::

  key: "nodes_xyz"
  type: Integer array
  ex: [4, 4, 4]
  default: No default
    
**Entry: Number of groupds (ng)**::

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
 note: need to add file format / contents description

**Entry: cross section info file name**::

 key: "xs_file"
 type: string
 default: 'xs_file'

_note: need to add more specific file format / contents description, see conversation in working notes

**Entry: source file name**::

        key: "source_input_file"
        type: string
        default: 'src.dat'

_note: need to add file format / contents description

**Entry: boundary condition file name [optional]**::

 key: "bc_input_file"
 type: string
 default: No default

_note: need to add file format / contents description

_Rachel, how should we deal with these? They're both required for the example to run, but future users of pyne won't necessarily have them relevant to their problems. MM

_We need to figure out how to generate them so we can give directions for that. Sebastian can help us with this part RS

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

_We need to define what this criterion means. That is, what is it checking exactly so the user can pick a good value. Our default might change as well. RS

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

 key: "tolerance"
 type: float
 ex: 1.e-10
 default: 1.e-10

_We also need to define what this criterion means. That is, what is it checking exactly so the user can pick a good value. Our default might change as well. RS

_I would call this relative_convergence_criteria. I would format double, since long is an integer format, right? But double precision has ~16 significant digits of precision, so maybe we should have a warning if the user goes beyond that?  MM

**Entry: Frequency of Checking Solution**::

 key: "need_a_better_name"
 type: Integer
 ex:  0
 default: 0

**Entry: Tolerance of Check**::

 key: "need_a_better_name"
 type: float
 ex: 1.e-14
 default: 1.e-14

_Call it ichk_tolerance because it implies a relation/dependence on ichk, which it has. I think the double precision question (slash warning) is relevant here too. MM

_We need better names on both of these based on what they actually do. RS

_I left it as ichk and ichk_tolerance.  I'll try and think of something better, if one of you comes up with a better name in the meantime let me know JH

_Sounds good. RS

_The three entries below may not be necessary.
    
**Entry: Max moment to print**::

 key: "max_mom_printed"
 type: Integer
 ex: 0
 default: 0

**Entry: flag for moment sum**::

 key: "moment_sum_flag"
 type: Integer
 ex: 1
 default: 0

**Entry: flag for moment at a point**::

 key: "mom_at_a_pt_flag"
 type: Integer
 ex: 1
 default: 0

_I'm not sure if we actually need any of these, so we may want to take all or some of them out. We don't have to take them out of the code for now, but we can keep the defaults as off and not tell the user about them. RS
_I removed all three of these for now, but they are still passed in on the python side.  That way, later on we can go back and add them or chose to remove them once again. JH

_great RS

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

_I'm not sure if we actually need either of these ones either. We don't have to take them out of the code for now, but we can keep the defaults as off and not tell the user about them. RS

_i'm guessing one of these is a zero MM


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



-----------------------------------
Source File Formatting
-----------------------------------

The spatial solver dictionary requires multiple input binary source files.  The required files are the following:
    **BC input file
    **Source input file
    **...

Here is a breif description of how each should be formatted.

  (1.) Source file:
      The source file is a file contaning source information for each cell ????.  The formatting is dependant on the solver
      you select.

      For the AHOTN solvers, the source file should be formatted as following.  
      There should be ng * nx * ny * nz * lambda * lambda * lambda ?source? values present.
      We will refer to the index of each ?source? value as (ng, nx, ny, nz, lambda_x, lambda_y, lambda_z).
      Thus, the source values should be in the following order:
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

          (1,.,.,.,.,.,.)

      When being read in, they will be iterated over by the following loop:
	      DO g=1,ng
		       DO ix=1,nx
		          DO iy=1,ny
		             DO iz=1,nz
		                DO jx=0,lambda
		                   DO jy=0,lambda
		                      DO jz=0,lambda
		                         READ(12) s(jx,jy,jz,ix,iy,iz,g)

    


.. currentmodule:: pyne.spatialsolver

All functionality may be found in the ``spatialsolver`` package::

 from pyne import spatialsolver

Spatialsolver API
-----------

.. automodule:: pyne.spatialsolver
    :members:

.. _Spatialsolver: http://something.com/

