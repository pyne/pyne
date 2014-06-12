def populate():
		solver_dict = {
		'spatial_order':1,
		'spatial_method':1,
		'quadrature_order':4,
		'quadrature_type':1,
		'nodes_xyz':[4,4,4],
		'num_groups':1,
		'num_materials':1,
		'x_cells_widths':[0.25, 0.25, 0.25, 0.25],
	 	'y_cells_widths':[0.25, 0.25, 0.25, 0.25],
		'z_cells_widths':[0.25, 0.25, 0.25, 0.25],
		'x_boundry_conditions':[2,2],
		'y_boundry_conditions':[2,2],
		'z_boundry_conditions':[2,2],
		'material_id': [[[1]*4]*4]*4,
		'quadrature_file':'quad_file',
		'xs_file':'xs',
		'srcfile':'src_4.dat',
		'converge_critical':1.e-12,
		'max_iterations':6000,
		'moments_converged':0,
		'converge_tolerence':1.e-10,
		'ichk':0,
		'ichk_tolerence':1.e-14,
		'momp':0,
		'momsum':0,
		'mompt':0,
		'qdflx':0
		};
		return solver_dict
'''
To create dictionary:
from dictionary_populate import populate
c = populate()
print(c)
'''


'''

Wrapping correlation:

titlein : input string(len=80)
"test title in"
lambdain : input int
input['spatial_order']
methin : input int
input['spatial_method']
qdordin : input int
input['quadrature_order']
qdtypin : input int
input['quadrature_type']
nxin : input int
input['nodes_xyz'][0]
nyin : input int
input['nodes_xyz'][1]
nzin : input int
input['nodes_xyz'][2]
ngin : input int
input['num_groups']
nmin : input int
input['num_materials']
dxin : input rank-1 array('d') with bounds (f2py_dxin_d0)
input['x_cells_widths']
dyin : input rank-1 array('d') with bounds (f2py_dyin_d0)
input['y_cells_widths']
dzin : input rank-1 array('d') with bounds (f2py_dzin_d0)
input['z_cells_widths']
xsbcin : input int
input['x_boundry_conditions'][0]
xebcin : input int
input['x_boundry_conditions'][1]
ysbcin : input int
input['y_boundry_conditions'][0]
yebcin : input int
input['y_boundry_conditions'][1]
zsbcin : input int
input['z_boundry_conditions'][0]
zebcin : input int
input['z_boundry_conditions'][1]
matin : input rank-3 array('i') with bounds (f2py_matin_d0,f2py_matin_d1,f2py_matin_d2)
input['material_id']
qdfilein : input string(len=30)
input['quadrature_file']
xsfilein : input string(len=30)
input['xs_file']
srcfilein : input string(len=30)
input['srcfile'
errin : input float
input['converge_critical']
itmxin : input int
input['max_iterations']
iallin : input int
input['moments_converged']
tolrin : input float
input['converge_tolerence']
tchkin : input float
input['ichk']
ichkin : input int
input['ichk_tolerence']
mompin : input int
input['momp']
momsumin : input int
input['momsum']
momptin : input int
input['mompt']
qdflxin : input int
input['qdflx']

'''
