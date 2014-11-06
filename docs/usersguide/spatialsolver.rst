.. currentmodule:: pyne.spatialsolver

.. _usersguide_spatialsolver:

==========================
Spatial Solver Users Guide
==========================



*****************************
Example Use of Spatial Solver
*****************************

.. code-block:: ipython

   In [1]: import pyne
   In [2]: import pyne.spatialsolver
   In [3]: input_dict = {'Name': 'Jane', 'Age': 27};
   In [4]: input_dict['solver'] = "AHOTN"
   In [5]: input_dict['solver_type'] = "LN" 
   In [6]: input_dict['spatial_order'] = 1
   In [7]: input_dict['spatial_method'] = 0
   In [8]: input_dict['angular_quadrature_order'] = 4
   In [9]: input_dict['angular_quadrature_type'] = 1
   In [10]: input_dict['nodes_xyz'] = [4,4,4]
   In [11]: input_dict['num_groups'] = 1
   In [12]: input_dict['num_materials'] = 1
   In [13]: input_dict['x_cells_widths'] = [0.25, 0.25, 0.25, 0.25]
   In [14]: input_dict['y_cells_widths'] = [0.25, 0.25, 0.25, 0.25]
   In [15]: input_dict['z_cells_widths'] = [0.25, 0.25, 0.25, 0.25]
   In [16]: input_dict['x_boundry_conditions'] = [2, 2]
   In [17]: input_dict['y_boundry_conditions'] = [2, 2]
   In [18]: input_dict['z_boundry_conditions'] = [2, 2]
   In [19]: input_dict['material_id'] = [ [ [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1] ], 
                                [ [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1] ],  
                                [ [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1] ],  
                                [ [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1] ] ]
   In [20]: input_dict['quadrature_file'] = 'quad_file'
   In [21]: input_dict['xs_file'] = 'xs_file'
   In [22]: input_dict['source_input_file'] = 'src_4.dat'
   In [23]: input_dict['bc_input_file'] = 'bc_4.dat'
   In [24]: input_dict['flux_output_file'] = 'phi_4.ahot'
   In [25]: input_dict['convergence_criterion'] = 1.e-12
   In [26]: input_dict['max_iterations'] = 6000
   In [27]: input_dict['moments_converged'] = 0
   In [28]: input_dict['converge_tolerence'] = 1.e-10
   In [29]: input_dict['ichk'] = 0
   In [30]: input_dict['ichk_tolerence'] = 1.e-14
   In [31]: dict_results = {}
   In [32]: dict_results = pyne.spatialsolver.solve(input_dict)
   In [33]: if(dict_results['success']):
   In [34]:   print('Yay, job ran succesfully!')
   In [35]: print(dict_results['success'])
   In [36]: print(dict_results['error_msg'])
   In [37]: print(dict_results['flux'])
   In [38]: print('Total solving time was: ')
   In [39]: print(dict_results['total_time'])
   In [40]: print('Solver call started at: ')
   In [41]: print(dict_results['time_start'])
   
 
For meanings of the above keys, see the pyne.spatialsolver python api!

For a more in depth example see the spatial solver section of the tutorials.
:ref:`pyne_spatialsolver`.
