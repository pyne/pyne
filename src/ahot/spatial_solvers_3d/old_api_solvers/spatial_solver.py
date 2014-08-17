"""This module provides a way to access neutron spatial solver codes.  It primarily consists of wrapped fortran code's being called with python via f2py"""

#Usual imports
from source.main import main as main

#imports being used for testing
from dict_util import dict_complete

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
	
