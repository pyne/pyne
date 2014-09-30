.. currentmodule:: pyne.spatialsolver

.. _usersguide_spatialsolver:

==========================
Spatial Solver Users Guide
==========================


--------------
Example of Use
--------------

import pyne.spatialsolver
solver_dict = {
    'solver':solverin,
    'solver_type':solvertypein,
    'spatial_order':1,
    'spatial_method':0,
    'angular_quadrature_order':4,
    'angular_quadrature_type':1,
    'nodes_xyz':[4,4,4],
    'num_groups':1,
    'num_materials':2,
    'x_cells_widths':[0.25, 0.25, 0.25, 0.25],
    'y_cells_widths':[0.25, 0.25, 0.25, 0.25],
    'z_cells_widths':[0.25, 0.25, 0.25, 0.25],
    'x_boundry_conditions':[2,2],
    'y_boundry_conditions':[2,2],
    'z_boundry_conditions':[2,2],
    'material_id':[
    [
    [1,2,1,2],
    [2,1,2,1],
    [1,2,1,2],
    [2,1,2,1],
    ],
    [
    [2,1,2,1],
    [1,2,1,2],
    [2,1,2,1],
    [1,2,1,2],
    ],
    [
    [1,2,1,2],
    [2,1,2,1],
    [1,2,1,2],
    [2,1,2,1],
    ],
    [
    [2,1,2,1],
    [1,2,1,2],
    [2,1,2,1],
    [1,2,1,2],
    ]
    ],
    'quadrature_file':'quad_file',
    'xs_file':'xs_alternating',
    'source_input_file':'src_4.dat',
    'bc_input_file':'bc_4.dat',
    'flux_output_file':'phi_4.ahot',
    'convergence_criterion':1.e-12,
    'max_iterations':6000,
    'moments_converged':0,
    'converge_tolerence':1.e-10,
    'ichk':0,
    'ichk_tolerence':1.e-14
    };

For meanings of the above keys, see the pyne.spatialsolver python api!

dict_results = pyne.spatialsolver.solve(a)
print(dict_results['flux']

For a more in depth example see the spatial solver section of the tutorials.
:ref:`pyne_spatialsolver`.
