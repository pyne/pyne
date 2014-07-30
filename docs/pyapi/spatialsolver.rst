.. _pyne_spatialsolver:

==================================================
Spatialsolver Support -- :mod:`pyne.spatialsolver`
==================================================
Spatialsolver is a pyne module that contains seven neutron transport equation solvers.  Each
solver is its own unique nodal method.

The solvers included in this module are listed below.  The theory and methodology behind each
can be found in the pyne theory doccumentation.
 #. **AHOTN-LN**: Arbitrarily higher order transport method of the nodal type linear-nodal method
 #. **AHOTN-LL**:  Arbitrarily higher order transport method of the nodal type linear-linear method
 #. **AHOTN-NEFD**: Arbitrarily higher order transport method of the nodal type that makes use of the
    unknown nodal flux moments (NEFD algorithm).
 #. **DGFEM-LD**: The Discontinuous Galerkin Finite Element Method (DGFEM) with a linear discontinuous (LD)
    approximation for angular flux. (SEE PAGE 27 of thesis)
 #. **DGFEM-DENSE**: The Discontinuous Galerkin Finite Element Method (DGFEM) that use ??dense??? lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
Use Lagrange polynomials to create a function
space per dimension, i.e. in x, y and z direction.
 #. **DGFEM-LAGRANGE**:   The Discontinuous Galerkin Finite Element Method (DGFEM) that use lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
 #. **SCTSTEP**: SCT Step algorithm similar to Duo's SCT algorithm implemented in three dimensional Cartesian
    geometry.

-----------------------------------
Dictionary Entries
-----------------------------------
As these are complicated solvers, they require a large amount of input data supplied by the user.  The
format we choose to take all this information in by is with a python dictionary.   Of the many key-pair values listed below, most are required, but some are optional.  The optional entries will be overriden by default values if not present/not specified. 

**Entry: Spatial expansion order (lambda; ahot spatial order, 0, 1, or 2)**::
    key: "spatial_order"
    type: Integer
    ex: 0
    default: 1

**Entry: Method (meth): 0/1 => AHOT-N-SI/AHOT-N-ITM**
    key: "spatial_method"   
    type: Integer
    ex: 0
    default: 0

**Entry: Quadrature order (can only be an even number)**
    key: "quadrature_order"
    type: Integer
    ex: 4
    default: 4

**Entry: Qudrature type:**
    key: "quadrature_type"
    type: Integer
    ex: 1
    default: 1

**Entry: Number of Nodes in x, y, and z directions (nx/ny/nz)**
    key: "nodes_xyz"
    type: Integer array
    ex: [4, 4, 4]
    default: No default
	
**Entry: Number of groupds (ng)**
 key: "num_groups"
 type: Integer
 ex: 1
 default: No default

**Entry: Number of Materials (nm)**
 key: "num_materials"
 type: Integer
 ex: 1
 default: No default

**Entry: x-size of cells (dx)**
 key: "x_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: y-size of cells (dy)**
 key: "y_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: z-size of cells (dz)**
 key: "z_cells_widths"
 type: double array
 ex: [0.25, 0.25, 0.25, 0.25]
 default: No default

**Entry: x start and end boundary conditions**
 key: "x_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default

**Entry: y start and end boundary conditions**
 key: "y_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default

**Entry: z start and end boundary conditions**
 key: "z_boundry_conditions"
 type: Integer array
 ex: [2, 2]
 default: No default

**Entry: Material info**
 key: "material_id"
 type: Integer 3 Dimensional Array
 ex: [ [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] ]
 default: No default
 Note:  Dimensions must match cell spacing and ordering

_Note: we need to give directions about the ordering. -RS_
	
**Entry: "quadrature_file" [optional; only needed for quadrature type 2]**
 type: string
 ex: 'quad_file'
 default: No default  
 note: need to add file format / contents description

**Entry: cross section info file name**
 key: "xs_file"
 type: string
 default: 'xs_file'
_note: need to add more specific file format / contents description
_note: see conversation in working notes

**Entry: source file name**
|		key: "source_input_file"
|		type: string
|		default: 'src.dat'
_note: need to add file format / contents description

**Entry: boundary condition file name [optional]**
 key: "bc_input_file"
 type: string
 default: No default
_note: need to add file format / contents description

_Rachel, how should we deal with these? They're both required for the example to run, but future users of pyne won't necessarily have them relevant to their problems. -MM_

_We need to figure out how to generate them so we can give directions for that. Sebastian can help us with this part -RS_

**Entry: output file name [optional]**
 key: "flux_output_file"
 type: string
 default: 'flux.out'
_note: need to add file format / contents description

**Entry: Convergence Criterion**
 key: "convergence_criterion"
 type: float
 ex: 1.e-12
 default: 1.e-12

_We need to define what this criterion means. That is, what is it checking exactly so the user can pick a good value. Our default might change as well. -RS_

**Entry: Maximum Iterations**
 key: "max_iterations"
 type: int
 ex: 10000
 default: 10000

**Entry: Moments Converged**
 key: "moments_converged" ??
 type: int
 ex: 0
 default: 0

_We need to provide a clear explanation of what this means. -RS_

**Entry: Tolerence**
 key: "tolerence"
 type: float
 ex: 1.e-10
 default: 1.e-10

_We also need to define what this criterion means. That is, what is it checking exactly so the user can pick a good value. Our default might change as well. -RS_

_I would call this relative_convergence_criteria. I would format double, since long is an integer format, right? But double precision has ~16 significant digits of precision, so maybe we should have a warning if the user goes beyond that?-MM_

**Entry: Frequency of Checking Solution**
 key: "need_a_better_name"
 type: Integer
 ex:  0
 default: 0

**Entry: Tolerance of Check**
 key: "need_a_better_name"
 type: float
 ex: 1.e-14
 default: 1.e-14

_Call it ichk_tolerance because it implies a relation/dependence on ichk, which it has. I think the double precision question (slash warning) is relevant here too. -MM_

_We need better names on both of these based on what they actually do. -RS_

_I left it as ichk and ichk_tolerance.  I'll try and think of something better, if one of you comes up with a better name in the meantime let me know -JH_

_Sounds good. -RS_

_The three entries below may not be necessary.
	
**Entry: Max moment to print**
 key: "max_mom_printed"
 type: Integer
 ex: 0
 default: 0

**Entry: flag for moment sum**
 key: "moment_sum_flag"
 type: Integer
 ex: 1
 default: 0 

**Entry: flag for moment at a point**
 key: "mom_at_a_pt_flag"
 type: Integer
 ex: 1 
 default: 0

_I'm not sure if we actually need any of these, so we may want to take all or some of them out. We don't have to take them out of the code for now, but we can keep the defaults as off and not tell the user about them. -RS_
_I removed all three of these for now, but they are still passed in on the python side.  That way, later on we can go back and add them or chose to remove them once again. -JH_

_great -RS_

**Entry: Flag for computing, presenting quadrant fluxes**
 key: "quad_flux_print_flag"
 type: Integer
 ex: 0
 default: 0

**Entry: Flag for printing matrix file**
 key: "matrix_print_flag"
 type: Integer
 ex: 0
 default: 0

_I'm not sure if we actually need either of these ones either. We don't have to take them out of the code for now, but we can keep the defaults as off and not tell the user about them. -RS_

_i'm guessing one of these is a zero -MM_


**Entry: ITM solution flag**
 key: "itm_direct_solution_flag" [only relevant if itm is selected]
 type: Integer
 ex: 0
 default: default of 0 if itm solution method

_We need to double check the meaning of this one. -RS_
_I have not added this yet.  As soon as we verify the meaning I can add it. -JH_
_great -RS_

.. currentmodule:: pyne.spatialsolver

All functionality may be found in the ``spatialsolver`` package::

 from pyne import spatialsolver

Spatialsolver API
-----------

.. automodule:: pyne.spatialsolver
    :members:

.. _Spatialsolver: http://something.com/
