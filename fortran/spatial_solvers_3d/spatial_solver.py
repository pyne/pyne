from AHOTN.src_LN.input import input
from AHOTN.src_LN.solve import solve

def ahotn(inputdict):
	print "in ahotan"
	print inputdict['material_id']
	input("test title in",
	inputdict['spatial_order'],
	inputdict['spatial_method'],
	inputdict['quadrature_order'],
	inputdict['quadrature_type'],
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
	inputdict['srcfile'],
	inputdict['converge_critical'],
	inputdict['max_iterations'],
	inputdict['moments_converged'],
	inputdict['converge_tolerence'],
	inputdict['ichk'],
	inputdict['ichk_tolerence'],
	inputdict['momp'],
	inputdict['momsum'],
	inputdict['mompt'],
	inputdict['qdflx'])
	print "INPUT WORKED!"
	print inputdict['qdflx']
	solve()
