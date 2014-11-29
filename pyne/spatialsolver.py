"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
from __future__ import division

import sys
import os
import time

#Solver imports
import pyne.transport_spatial_methods as transport_spatial_methods

def solve(inputdict_unchecked):
    """ 
    As these are complicated solvers, they require a large amount of input data 
    supplied by the user.  This information needs to be entered as a Python 
    dictionary. Of the many key-pair values listed below, most are required, but 
    some are optional.  The optional entries will be overridden by default values 
    if not present/specified.

    **Input Dictionary Entry: Solver type (AHOTN, DGFEM or SCTSTEP)**
    ::  

      key: "solver"
      type: String
      ex: "AHOTN"
      default: no default
      Note:  The three supported "solver"s are the following:
        1.  "AHOTN"
        2.  "DGFEM"
        3.  "SCTSTEP"
     
    **Input Dictionary Entry: Spatial expansion order**
    ::
      
      key: "solver_type"
      type: String
      ex: "LN"
      default: No default
      Note: This entry is dependent on the "solver" entry.
        For AHOTN solver, there exist the "LN", "LL", and "NEFD" solver types
        For the DGFEM solvers, there exist the "LD", "DENSE", and "LAGRANGE" solver types
        For the SCTSTEP solver, the "solver_type" key is not used.

    **Input Dictionary Entry: Spatial expansion order (lambda; ahot spatial order, 0, 1, or 2)**
    ::

      key: "spatial_order"
      type: Integer
      ex: 0
      default: 1

      The Spatial expansion order is the expansion order of the spatial moment.

    **Input Dictionary Entry: Angular quadrature order**
    ::

      key: "quadrature_order"
      type: Integer
      ex: 4
      default: 4

      The angular quadrature order is the number of angles to be used per octant.  
      For N sets of angles, there will be (N * (N + 2) / 8) ordinates per octant. 
      The quadrature order may only be an even number!

    **Input Dictionary Entry: Quadrature type:**
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

    **Input Dictionary Entry: Number of spatial nodes in x, y, and z directions (nx/ny/nz)**
    ::

      key: "nodes_xyz"
      type: Integer array
      ex: [4, 4, 4]
      default: No default
        
    **Input Dictionary Entry: Number of energy groups (ng)**
    ::

     key: "num_groups"
     type: Integer
     ex: 1
     default: No default

    **Input Dictionary Entry: Number of materials (nm)**
    ::

     key: "num_materials"
     type: Integer
     ex: 1
     default: No default

    **Input Dictionary Entry: x-size of cells (dx)**
    ::

     key: "x_cells_widths"
     type: double array
     ex: [0.25, 0.25, 0.25, 0.25]
     default: No default

    **Input Dictionary Entry: y-size of cells (dy)**
    ::

     key: "y_cells_widths"
     type: double array
     ex: [0.25, 0.25, 0.25, 0.25]
     default: No default

    **Input Dictionary Entry: z-size of cells (dz)**
    ::

     key: "z_cells_widths"
     type: double array
     ex: [0.25, 0.25, 0.25, 0.25]
     default: No default

    **Input Dictionary Entry: x start and end boundary conditions**
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

    **Input Dictionary Entry: y start and end boundary conditions**
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

    **Input Dictionary Entry: z start and end boundary conditions**
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

    **Input Dictionary Entry: Material info**
    ::

     key: "material_id"
     type: Integer 3 dimensional array
     ex: [ [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] 
           [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] 
           [ [1 1 1 1] [1 1 1 1] [1 1 1 1] [1 1 1 1] ] ]
     default: No default
     note: Dimensions must match cells such that there is one material number
           in each spatial cell. The cells are ordered as x, y, z.

        
    **Input Dictionary Entry: "quadrature_file" [optional; only needed for quadrature type 2]**
    ::

     key: "quad_file"
     type: String
     ex: 'quad_file'
     default: No default  
     note: See input file formatting notes in the Quadrature File Formatting section.

    **Input Dictionary Entry: cross section info file name**
    ::

     key: "xs_file"
     type: String
     default: 'xs_file'
     note: See input file formatting notes in the Cross Section File Formatting section.

    **Input Dictionary Entry: source file name**
    ::

      key: "source_input_file"
      type: String
      default: 'src.dat'
      note: See input file formatting notes in the Source File Formatting section.

    **Input Dictionary Entry: boundary condition file name [optional]**
    ::

     key: "bc_input_file"
     type: String
     default: No default
     note: See input file formatting notes in the Boundry Condition File Formatting section.

    **Input Dictionary Entry: output file name [optional]**
    ::

     key: "flux_output_file"
     type: String
     default: 'flux.out'
     note: See input file formatting notes in the Flux Output (5.) File Formatting section.

    **Input Dictionary Entry: Convergence Criterion**
    ::

     key: "convergence_criterion"
     type: float
     ex: 1.e-5
     default: 1.e-5
     The solution is considered converged and the calculation completes when the flux
     in each cell at the current iteration is within "convergence_criterion" of the
     previous iterate. This is generally the relative difference, but in cases of 
     very small flux values the absolute difference is used instead (see the 
     Convergence Tolerance entry below).  

    **Input Dictionary Entry: Tolerance**
    ::

     key: "converge_tolerance"
     type: float
     ex: 1.e-10
     default: 1.e-10
     Converge tolerance is the tolerance that determines how the difference between
     flux iterates (df) that is used to determine convergence will be calculated. 
     df is calculated as follows:
       f = current flux value
       ct = convergence tolerance (value for this key, "converge_tolerance")
       f1 = flux value from the previous iteration
       If f1 > ct:
         df = absolute(f - f1) / f1
       Else
         df = absolute(f - f1)
     The idea is to use the absolute difference instead of the relative difference
     between iterates when the flux is very small to help avoid rounding error.

    **Input Dictionary Entry: Maximum Iterations**
    ::

     key: "max_iterations"
     type: int
     ex: 10000
     default: 10000
     note: If this number of iterations is reached before the convergence criterion
           is satisfied, the calculation will terminate and report the current flux
           estimate.

    **Input Dictionary Entry: Moments Converged**
    ::

     key: "moments_converged"
     type: int
     ex: 0
     default: 0
     Moments converged is the number of moments that should be converged upon for each quadrature in the
     solution space.  Value for moments converged must be in range [0, spatial_order_in].



-----------------------------------
Output Dictionary Entries
-----------------------------------
When run, the solvers return a dictionary of useful solution data.  It contains the following key-pair entries:

**Output Dictionary Entry: Flux output array**
::

  key:  "flux"
  type: Double Array of 3 dimensions
  format: Flux output array is in following format:
  Each cell in the array has a scalar flux, the integral of the angular
  flux over all angles in that cell.   The first index refers to the 
  plane on the z axis, beginning at 0 with the lowest plane, and moving
  upwards to the highest plane on the mesh.  The second index is the 
  row on the z plane, and the third index is the cell in the row.

  format examples: If you had a mesh with 4 by 4 by 4 cells extending
  in the x, y and z directions, then to get the following flux values,
  you would use the following indices:

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

**Output Dictionary Entry: Solver success code**
::

  key:  "success"
  type: Integer
  format: 1 means yes, the solve succeeded.  0 means it failed.

**Output Dictionary Entry: Raw system time of solver start**
::

  time_start provides you with the system time when the solver began running.
  key:  "time_start"
  type: double
  format: system time

**Output Dictionary Entry: Total run time**
::

  total_time is the total time the solver took to solve.
  key:  "total_time"
  type: double
  format: system time

**Output Dictionary Entry: Total print time**
::

  print_time is the total time the solver took to print results.
  key:  "print_time"
  type: double
  format: system time

**Output Dictionary Entry: Error Message**
::

  If the solver fails, error_msg is a string describing why the solver failed.
  key:  "error_msg"
  type: String

    """

    inputdict = _dict_complete(inputdict_unchecked)
    flux_output = []
    solver_output = {};
    fortran_returns = transport_spatial_methods.main("test title in",
        inputdict['solver'],
        inputdict['solver_type'],
        inputdict['spatial_order'], 
        #inputdict['spatial_method'],
        inputdict['angular_quadrature_order'],
        inputdict['angular_quadrature_type'],
        inputdict['nodes_xyz'][0],
        inputdict['nodes_xyz'][1],
        inputdict['nodes_xyz'][2],
        inputdict['num_groups'],
        inputdict['num_materials'],
        inputdict['x_cells_widths'],
        inputdict['y_cells_widths'],
        inputdict['z_cells_widths'],
        inputdict['x_boundry_conditions'][0],
        inputdict['x_boundry_conditions'][1],
        inputdict['y_boundry_conditions'][0],
        inputdict['y_boundry_conditions'][1],
        inputdict['z_boundry_conditions'][0],
        inputdict['z_boundry_conditions'][1],
        inputdict['material_id'],
        inputdict['quadrature_file'],
        inputdict['xs_file'],
        inputdict['source_input_file'],
        inputdict['bc_input_file'],
        inputdict['flux_output_file'],
        inputdict['convergence_criterion'],
        inputdict['max_iterations'],
        inputdict['moments_converged'],
        inputdict['converge_tolerence'],
        inputdict['max_mom_printed'],
        inputdict['moment_sum_flag'],
        inputdict['mom_at_a_pt_flag'],
        inputdict['quad_flux_print_flag'])

    solver_output['flux'] = fortran_returns[6].tolist()
    error_code = fortran_returns[7]
    tsolve = fortran_returns[8]
    ttosolve = fortran_returns[9]
    tend = fortran_returns[10]
  
    if(error_code == 0):
        solver_output['success'] = 1
        solver_output['time_start'] = tsolve
        solver_output['total_time'] = tsolve-ttosolve
        solver_output['print_time'] = tend-ttosolve
        solver_output['error_msg'] = 0
    else:
        solver_output['success'] = 0
        solver_output['error_msg'] = _error_to_string(error_code)
        print(solver_output['error_msg'])
    return solver_output

def _error_to_string(error_code):
    err_dictionary = {
        1001: "ERROR: Lambda must be equal to one.",
        1002: "ERROR: Illegal value for qdord. Must be greater than zero.",
        1003: "ERROR: Illegal value for the quadrature order. Even #s only.",
        1004: "ERROR: illegal entry for the qdfile name.",
        1005: "ERROR: illegal entry for lambda. Must be zero or greater.",
        1007: "ERROR: Illegal number of x cells. Must be positive.",
        1008: "ERROR: Illegal number of y cells. Must be positive.",
        1009: "ERROR: Illegal number of z cells. Must be positive.",
        1010: "ERROR: Illegal number of energy groups. Must be positive.",
        1011: "ERROR: Illegal number of materials. Must be positive.",
        1012: "ERROR: Illegal x cell dimension, dx. Must be positive.",
        1013: "ERROR: Illegal y cell dimension, dy. Must be positive.",
        1014: "ERROR: Illegal z cell dimension, dz. Must be positive.",
        1015: "ERROR: Illegal value in material map. Must be in [1, #materials].",
        1016: "ERROR: Illegal lower x BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1017: "ERROR: Illegal lower y BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1018: "ERROR: Illegal lower z BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1019: "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed.",
        1020: "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed.",
        1021: "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed.",
        1022: "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed.",
        1023: "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed.",
        1024: "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed.",
        1025: "ERROR: Illegal convergence criterion. Must be positive.",
        1026: "ERROR: Illegal max inner iterations, itmx. Must be positive.",
        1027: "ERROR: Illegal tolerance setting, tolr. Must be positive.",
        1028: "ERROR: Illegal value for moments to converge, iall. Must be in [0, lambda].",
        1029: "ERROR: Illegal value for solution check flag, iall. iall is 0 (skip check) or positive integer",
        1030: "ERROR: Illegal value for solution check tolerance, tchk. Must be positive.",
        1031: "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1.",
        1032: "ERROR: Illegal value for weight. Must be entered positive, less than 0.125.",
        1033: "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values.",
        1034: "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda.",
        1035: "ERROR: Illegal value for flag for moment summing. Must be 0 for off or 1 for on.",
        1036: "ERROR: Illegal value for flag for moment sum at non-center point. Must be 0/1=off/on.",
        1037: "ERROR: Illegal value for flag for printing average flux of the four quadrants. Must be 0/1 = off/on."
  };
    return err_dictionary[error_code] 


def _dict_complete(inputdict):

    formatted_dict = {}
    try:
        if((inputdict['solver'] == "AHOTN") or (inputdict['solver']=="DGFEM") or (inputdict['solver']=="SCTSTEP")):
            formatted_dict['solver'] = inputdict['solver']
        else:
            assert (0==1), "solver does not exist"
    except:
        assert  (0==1), "solver key does not exist"
    try:
        formatted_dict['solver_type'] = inputdict['solver_type']
    except:
        if(inputdict['solver'] == "AHOTN"):
            formatted_dict['solver_type'] = "LN"
        elif(inputdict['solver'] == "DGFEM"):
            formatted_dict['solver_type'] = "LD"

    assert 'nodes_xyz' in inputdict, 'nodes_xyz key not in dict'
    formatted_dict['nodes_xyz'] = inputdict['nodes_xyz']
    assert 'num_groups' in inputdict, 'num_groups key not in dict'
    formatted_dict['num_groups'] = inputdict['num_groups']
    assert 'num_materials' in inputdict, 'num_materials key not in dict'
    formatted_dict['num_materials'] = inputdict['num_materials']
    assert 'x_cells_widths' in inputdict, 'x_cells_widths not in dict'
    formatted_dict['x_cells_widths'] = inputdict['x_cells_widths']
    assert 'y_cells_widths' in inputdict, 'y_cells_widths not in dict'
    formatted_dict['y_cells_widths'] = inputdict['y_cells_widths']
    assert 'z_cells_widths' in inputdict, 'z_cells_widths not in dict'
    formatted_dict['z_cells_widths'] = inputdict['z_cells_widths']
    assert 'x_boundry_conditions' in inputdict, 'x_boundry_conditions not in dict'
    formatted_dict['x_boundry_conditions'] = inputdict['x_boundry_conditions']
    assert 'y_boundry_conditions' in inputdict, 'y_boundry_conditions not in dict'
    formatted_dict['y_boundry_conditions'] = inputdict['y_boundry_conditions']
    assert 'z_boundry_conditions' in inputdict, 'z_boundry_conditions not in dict'
    formatted_dict['z_boundry_conditions'] = inputdict['z_boundry_conditions']
    assert 'material_id' in inputdict, 'material_id not in dict'
    formatted_dict['material_id'] = inputdict['material_id']
    assert 'quadrature_file' in inputdict, 'quadrature_file not in dict'
    formatted_dict['quadrature_file'] = inputdict['quadrature_file']
    assert 'xs_file' in inputdict, 'xs_file not in dict'
    formatted_dict['xs_file'] = inputdict['xs_file']
    assert 'source_input_file' in inputdict, 'source_input_file not in dict'
    formatted_dict['source_input_file'] = inputdict['source_input_file']
    assert 'bc_input_file' in inputdict, 'bc_input_file not in dict'
    formatted_dict['bc_input_file'] = inputdict['bc_input_file']
    assert 'flux_output_file' in inputdict, 'flux_output_file not in dict'
 
    formatted_dict['spatial_order'] = inputdict.get('spatial_order', 1)
    formatted_dict['angular_quadrature_order'] = inputdict.get('angular_quadrature_order',4)
    formatted_dict['angular_quadrature_type'] = inputdict.get('angular_quadrature_type',1)
    formatted_dict['flux_output_file'] = inputdict['flux_output_file']  
    formatted_dict['convergence_criterion'] = inputdict.get('convergence_criterion',1e-5)
    formatted_dict['max_iterations'] = inputdict.get('max_iterations', 6000)
    formatted_dict['moments_converged'] = inputdict.get('moments_converged',0)
    formatted_dict['converge_tolerence'] = inputdict.get('converge_tolerence',1e-10)
    formatted_dict['max_mom_printed'] = 0
    formatted_dict['moment_sum_flag'] = 0
    formatted_dict['mom_at_a_pt_flag'] = 0
    formatted_dict['quad_flux_print_flag'] = 0

    return formatted_dict


