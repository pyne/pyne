from warnings import warn

def dict_complete:

		formatted_dict = {}
		try:
			formatted_dict['spatial_order'] = inputdict['spatial_order']
		except:
			formatted_dict['spatial_order'] = 1
		try:
			formatted_dict['spatial_method'] = inputdict['spatial_method']
		except:
			formatted_dict['spatial_method'] = 1
		try:
			formatted_dict['quadrature_order'] = inputdict['quadrature_order']
		except:
			formatted_dict['quadrature_order'] = 4
		try:
			formatted_dict['quadrature_type'] = inpudict['quadrature_type']
		except:
			formatted_dict['quadrature_type'] = 1
		try:
			formatted_dict['nodes_xyz'] = inputdict['nodes_xyz']
		except:
			return 'input dictionary missing nodes_xyz'
		try:
			formatted_dict['num_groups'] = inputdict['num_groups']
		except:
			formatted_dict['num_groups'] = 1
		try:
			formatted_dict['num_materials'] = inputdict['num_materials']
		except:
			formatted_dict['num_materials'] = 1
		try:
			formatted_dict['x_cells_widths'] = inputdict['x_cells_widths']
		except:
			return 'input dictionary missing x_cells_widths'		
		try:
			formatted_dict['y_cells_widths'] = inputdict['y_cells_widths']
		except:
			return 'input dictionary missing y_cells_widths'		
		try:
			formatted_dict['z_cells_widths'] = inputdict['z_cells_widths']
		except:
			return 'input dictionary missing z_cells_widths'
		try:
			formatted_dict['x_boundry_conditions'] = inputdict['x_boundry_conditions']
		except:
			return 'input dictionary missing x_boundry_conditions'
		try:
			formatted_dict['y_boundry_conditions'] = inputdict['y_boundry_conditions']
		except:
			return 'input dictionary missing y_boundry_conditions'
		try:
			formatted_dict['z_boundry_conditions'] = inputdict['z_boundry_conditions']
		except:
			return 'input dictionary missing z_boundry_conditions'
		try:
			formatted_dict['material_id'] = inputdict['material_id']
		except:
			return 'input dictionary missing material_id file'
		try:
			formatted_dict['quadrature_file'] = inputdict['quadrature_file']
		except:
			return 'input dictionary missing quadrature_file'
		try:
			formatted_dict['xs_file'] = inputdict['xs_file']
		except:
			return 'input dictionary missing xs_file'
		try:
			formatted_dict['converge_critical'] = inputdict['converge_critical']
		except:
			formatted_dict['converge_critical'] = 1e-12
		try:
			formatted_dict['max_iterations'] = inputdict['max_iterations']
		except:
			formatted_dict['max_iterations'] = 6000
		try:
			formatted_dict['moments_converged'] = inputdict['moments_converged']
		except:
			formatted_dict['moments_converged'] = 0
		try:
			formatted_dict['converge_tolerence'] = inputdict['converge_tolerence']
		except:
			formatted_dict['converge_tolerence'] = 1e-10
		try:
			formatted_dict['ichk'] = inputdict['ichk']
		except:
			formatted_dict['ichk'] = 0
		try:
			formatted_dict['ichk_tolerence'] = inputdict['ichk_tolerence']
		except:
			formatted_dict['ichck_tolerence'] = 1e-14
		try:
			formatted_dict['momsum'] = inputdict['momsum']
		except:
			formatted_dict['momsum'] = 0
		try:
			formatted_dict['mompt'] = inputdict['mompt']
		except:
			formatted_dict['mompt'] = 0
		try:
			formatted_dict['qdflx'] = inputdict['qdflx']
		except:
			formatted_dict['qdflx'] = 0
			
		return formatted_dict
