"""This file provides utility functions for the dictionary that is passed into various solvers on in spatial_solver.  It checks that all required key-pair values are present, then fills in any optional key-pair values that are not present.  Warns the user for each overriden key-pair value."""

# Python imports 
from warnings import warn

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
				formatted_dict['solver_type'] = "NEED_DEFAULT..."
		
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
			formatted_dict['num_groups'] = 1
			warn(warning_msg + " num_groups value of 1")
		try:
			formatted_dict['num_materials'] = inputdict['num_materials']
		except:
			formatted_dict['num_materials'] = 1
			warn(warning_msg + " num_materials value of 1")
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
			formatted_dict['src_file'] = inputdict['src_file']
		except:
			raise InputDictError("src_file")
		try:
			formatted_dict['converge_critical'] = inputdict['converge_critical']
		except:
			formatted_dict['converge_critical'] = 1e-12
			warn(warning_msg + " converge_critical value of 1e-12")
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
		try:
			formatted_dict['momsum'] = inputdict['momsum']
		except:
			formatted_dict['momsum'] = 0
			warn(warning_msg + " momsum value of 0")
		try:
			formatted_dict['momp'] = inputdict['momp']
		except:
			formatted_dict['momp'] = 0
			warn(warning_msg + " momp value of 0")
		try:
			formatted_dict['mompt'] = inputdict['mompt']
		except:
			formatted_dict['mompt'] = 0
			warn(warning_msg + " mompt value of 0")
		try:
			formatted_dict['qdflx'] = inputdict['qdflx']
		except:
			formatted_dict['qdflx'] = 0
			warn(warning_msg + " qdflx value of 0")
			
		return formatted_dict


class InputDictError(Exception):
    def __init__(self, missing=None):
        self.missing = missing

    def __str__(self):
        msg = "Input dictionary missing required key value pair: "
        if self.missing is not None:
            msg += ": " + missing 
        return msg
