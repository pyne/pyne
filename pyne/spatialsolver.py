"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
from __future__ import division

import sys
import os

#Solver imports
#sys.path.append("../fortran/spatial_solvers_3d/source")
sys.path.append("../pyne_transport/pyne/fortran/spatial_solvers_3d/source")
#sys.path.append("../../fortran/spatial_solvers_3d/source")

#goal_dir = os.path.join(os.getcwd(), "../../fortran/spatial_solvers_3d/source")
from main import main as main
#from ... fortran.spatial_solvers_3d.source.main import main as main

#imports being used for testing
#from dict_util import dict_complete


def solve(inputdict_unchecked):
	inputdict = dict_complete(inputdict_unchecked)
	if(inputdict['solver'] == "AHOTN" or inputdict['solver'] == "DGFEM"):
		main("test title in",
		inputdict['solver'],
		inputdict['solver_type'],
		inputdict['spatial_order'],
		inputdict['spatial_method'],
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
		inputdict['convergence_criterion'],
		inputdict['max_iterations'],
		inputdict['moments_converged'],
		inputdict['converge_tolerence'],
		inputdict['ichk'],
		inputdict['ichk_tolerence'],
		inputdict['max_mom_printed'],
		inputdict['moment_sum_flag'],
		inputdict['mom_at_a_pt_flag'],
		inputdict['quad_flux_print_flag'])
	elif(inputdict['solver'] == "SCT-STEP"):
		print("SCT-STEP NOT IMPLEMENTED YET...")
	else:
		#Throw error
		print("Not a supported solver")





def dict_complete(inputdict):

		warning_msg = "Input dictionary taking on default "

		formatted_dict = {}
		try:
			if((inputdict['solver'] == "AHOTN") or (inputdict['solver']=="DGFEM")):
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
		try:
			formatted_dict['spatial_method'] = inputdict['spatial_method']
		except:
			formatted_dict['spatial_method'] = 0
			warn(warning_msg + " spatial_method value of 0")
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
			formatted_dict['convergence_criterion'] = 1e-12
			warn(warning_msg + " convergence_criterion value of 1e-12")
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
		try:
			formatted_dict['ichk'] = inputdict['ichk']
		except:
			formatted_dict['ichk'] = 0
			warn(warning_msg + " ichk value of 0")
		try:
			formatted_dict['ichk_tolerence'] = inputdict['ichk_tolerence']
		except:
			formatted_dict['ichk_tolerence'] = 1e-14
			warn(warning_msg + " ichk_tolerence value of 1e-14")

		formatted_dict['max_mom_printed'] = 0
		formatted_dict['moment_sum_flag'] = 0
		formatted_dict['mom_at_a_pt_flag'] = 0
		formatted_dict['quad_flux_print_flag'] = 0

		return formatted_dict


class InputDictError(Exception):
    def __init__(self, missing=None):
        self.missing = missing

    def __str__(self):
        msg = "Input dictionary missing required key value pair: "
        if self.missing is not None:
            msg += ": " + missing
        return msg
