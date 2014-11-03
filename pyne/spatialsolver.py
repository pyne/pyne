"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
from __future__ import division

import sys
import os
import time

#Solver imports
#sys.path.append("../pyne_transport/pyne/fortran/spatial_solvers_3d/source")
import pyne.transport_spatial_methods as transport_spatial_methods

#imports being used for testing
#from dict_util import dict_complete

def solve(inputdict_unchecked):
  inputdict = dict_complete(inputdict_unchecked)
  flux_output = []
  solver_output = {
     };
  #if(inputdict['solver'] == "AHOTN" or inputdict['solver'] == "DGFEM"):
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
  inputdict['quad_flux_print_flag']#,
#    inputdict['out_dims']
    )
    #time.sleep(.5)
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
    solver_output['error_msg'] = error_toString(error_code)
    print(solver_output['error_msg'])
  return solver_output
  #elif(inputdict['solver'] == "SCTSTEP"):
  #  print("SCT-STEP NOT IMPLEMENTED YET...")
  #  return null
  #else:
    #Throw error
  #  print("Not a supported solver")

def error_toString(error_code):
  err_str = ""  
  
  if(error_code == 1001):
    err_str = "ERROR: Lambda must be equal to one."
  elif(error_code == 1002):
    err_str = "ERROR: Illegal value for qdord. Must be greater than zero."
  elif(error_code == 1003):
    err_str = "ERROR: Illegal value for the quadrature order. Even #s only."
  elif(error_code == 1004):
    err_str = "ERROR: illegal entry for the qdfile name."
  elif(error_code == 1005):
    err_str = "ERROR: illegal entry for lambda. Must be zero or greater."
  #Error code 1006 removed.  Leaving spot 1006 as a placeholder for new error
  elif(error_code == 1007):
    err_str = "ERROR: Illegal number of x cells. Must be positive."
  elif(error_code == 1008):
    err_str = "ERROR: Illegal number of y cells. Must be positive."
  elif(error_code == 1009):
    err_str = "ERROR: Illegal number of z cells. Must be positive."
  elif(error_code == 1010):
    err_str = "ERROR: Illegal number of energy groups. Must be positive."
  elif(error_code == 1011):
    err_str = "ERROR: Illegal number of materials. Must be positive."
  elif(error_code == 1012):
    err_str = "ERROR: Illegal x cell dimension, dx. Must be positive."
  elif(error_code == 1013):
    err_str = "ERROR: Illegal y cell dimension, dy. Must be positive."
  elif(error_code == 1014):
    err_str = "ERROR: Illegal z cell dimension, dz. Must be positive."
  elif(error_code == 1015):
    err_str = "ERROR: Illegal value in material map. Must be in [1, #materials]."
  elif(error_code == 1016):
    err_str = "ERROR: Illegal lower x BC. Must be 0-Vac, 1-Refl or 2-Fixed."
  elif(error_code == 1017):
    err_str = "ERROR: Illegal lower y BC. Must be 0-Vac, 1-Refl or 2-Fixed."
  elif(error_code == 1018):
    err_str = "ERROR: Illegal lower z BC. Must be 0-Vac, 1-Refl or 2-Fixed."
  elif(error_code == 1019):
    err_str = "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1020):
    err_str = "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1021):
    err_str = "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1022):
    err_str = "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1023):
    err_str = "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1024):
    err_str = "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed."
  elif(error_code == 1025):
    err_str = "ERROR: Illegal convergence criterion. Must be positive."
  elif(error_code == 1026):
    err_str = "ERROR: Illegal max inner iterations, itmx. Must be positive."
  elif(error_code == 1027):
    err_str = "ERROR: Illegal tolerance setting, tolr. Must be positive."
  elif(error_code == 1028):
    err_str = "ERROR: Illegal value for moments to converge, iall. Must be in [0, lambda]."
  elif(error_code == 1029):
    err_str = "ERROR: Illegal value for solution check flag, iall. iall is 0 (skip check) or positive integer"
  elif(error_code == 1030):
    err_str = "ERROR: Illegal value for solution check tolerance, tchk. Must be positive."
  elif(error_code == 1031):
    err_str = "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
  elif(error_code == 1032):
    err_str = "ERROR: Illegal value for weight. Must be entered positive, less than 0.125."
  elif(error_code == 1033):
    err_str = "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values."
  elif(error_code == 1034):
    err_str = "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda."
  elif(error_code == 1035):
    err_str = "ERROR: Illegal value for flag for moment summing. Must be 0 for off or 1 for on."
  elif(error_code == 1036):
    err_str = "ERROR: Illegal value for flag for moment sum at non-center point. Must be 0/1=off/on."
  elif(error_code == 1037):
    err_str = "ERROR: Illegal value for flag for printing average flux of the four quadrants. Must be 0/1 = off/on."

  return "Error code: %d" % (error_code,) + " " + err_str

def dict_complete(inputdict):

    warning_msg = "Input dictionary taking on default "

    formatted_dict = {}
    try:
      if((inputdict['solver'] == "AHOTN") or (inputdict['solver']=="DGFEM") or (inputdict['solver']=="SCTSTEP")):
        formatted_dict['solver'] = inputdict['solver']
      else:
        raise InputDictError("solver does not exist")
    except:
      raise InputDictError("solver")
    try:
      formatted_dict['solver_type'] = inputdict['solver_type']
    except:
      #raise InputDictError("solver_type")
       if(inputdict['solver'] == "AHOTN"):
        formatted_dict['solver_type'] = "LN"
       elif(inputdict['solver'] == "DGFEM"):
        formatted_dict['solver_type'] = "LD"
		
    try:
      formatted_dict['spatial_order'] = inputdict['spatial_order']
    except:
      formatted_dict['spatial_order'] = 1
      warn(warning_msg + " spatial_order value of 1")
#No longer needed.  Spatial method of 1 was never implemented!
    #try:
    #  formatted_dict['spatial_method'] = inputdict['spatial_method']
    #except:
    #  formatted_dict['spatial_method'] = 0
    #  warn(warning_msg + " spatial_method value of 0")
    try:
      formatted_dict['angular_quadrature_order'] = inputdict['angular_quadrature_order']
    except:
      formatted_dict['angular_quadrature_order'] = 4
      warn(warning_msg + " angular_quadrature_order value of 4")
    try:
      formatted_dict['angular_quadrature_type'] = inputdict['angular_quadrature_type']
    except:
      formatted_dict['qangular_uadrature_type'] = 1
      warn(warning_msg + " angular_quadrature_type value of 1")
    try:
      formatted_dict['nodes_xyz'] = inputdict['nodes_xyz']
    except:
      raise InputDictError("nodes_xyz")
    try:
      formatted_dict['num_groups'] = inputdict['num_groups']
    except:
      raise InputDictError("num_groups")
    try:
      formatted_dict['num_materials'] = inputdict['num_materials']
    except:
      raise InputDictError('num_materials')
    try:
      formatted_dict['x_cells_widths'] = inputdict['x_cells_widths']
    except:
      raise InputDictError("x_cells_widths")		
    try:
      formatted_dict['y_cells_widths'] = inputdict['y_cells_widths']
    except:
      raise InputDictError("y_cells_widths")	
    try:
      formatted_dict['z_cells_widths'] = inputdict['z_cells_widths']
    except:
      raise InputDictError("z_cells_widths")
    try:
      formatted_dict['x_boundry_conditions'] = inputdict['x_boundry_conditions']
    except:
      raise InputDictError("x_boundry_conditions")
    try:
      formatted_dict['y_boundry_conditions'] = inputdict['y_boundry_conditions']
    except:
      raise InputDictError("y_boundry_conditions")
    try:
      formatted_dict['z_boundry_conditions'] = inputdict['z_boundry_conditions']
    except:
      raise InputDictError("z_boundry_conditions")
    try:
      formatted_dict['material_id'] = inputdict['material_id']
    except:
      raise InputDictError("material_id file")
    try:
      formatted_dict['quadrature_file'] = inputdict['quadrature_file']
    except:
      raise InputDictError("quadrature_file")
    try:
		  formatted_dict['xs_file'] = inputdict['xs_file']
    except:
      raise InputDictError("xs_file")
    try:
      formatted_dict['source_input_file'] = inputdict['source_input_file']
    except:
      raise InputDictError("source_input_file")
    try:
      formatted_dict['bc_input_file'] = inputdict['bc_input_file']
    except:
      raise InputDictError("bc_input_file")
    try:
      formatted_dict['flux_output_file'] = inputdict['flux_output_file']
    except:
      raise InputDictError("flux_output_file")
    try:
      formatted_dict['convergence_criterion'] = inputdict['convergence_criterion']
    except:
      formatted_dict['convergence_criterion'] = 1e-5
      warn(warning_msg + " convergence_criterion value of 1e-5")
    try:
      formatted_dict['max_iterations'] = inputdict['max_iterations']
    except:
      formatted_dict['max_iterations'] = 6000
      warn(warning_msg + " max_iterations value of 6000")
    try:
      formatted_dict['moments_converged'] = inputdict['moments_converged']
    except:
      formatted_dict['moments_converged'] = 0
      warn(warning_msg + " moments_converged value of 0")
    try:
      formatted_dict['converge_tolerence'] = inputdict['converge_tolerence']
    except:
      formatted_dict['converge_tolerence'] = 1e-10
      warn(warning_msg + " converge_tolerence value of 1e-10")
			
    formatted_dict['max_mom_printed'] = 0
    formatted_dict['moment_sum_flag'] = 0
    formatted_dict['mom_at_a_pt_flag'] = 0
    formatted_dict['quad_flux_print_flag'] = 0
#!if solver == ahotn if solvertype == ln
#    formatted_dict['out_dims'] = [4,formatted_dict['nodes_xyz'][0],formatted_dict['nodes_xyz'][1],formatted_dict##['nodes_xyz'][2],formatted_dict['num_groups'],1,1];

    return formatted_dict


class InputDictError(Exception):
    global missing
    def __init__(self, missing=None):
        self.missing = missing

    def __str__(self):
        msg = "Input dictionary missing required key value pair: "
        if self.missing is not None:
            msg += ": " + missing 
        return msg

