EXAMPLE README:

NOTE:  This example is set up in a very very very poor fashon.
Currently, the pyne hookup is not complete, so while the built
pyne module is called, it relies on the pyne/fortran source folder
from which you built pyne from.  It also requires the compile script
in pyne/fortran/spatial_solvers_3d/source to be run. See below for
the old example.

To run the new pyne interface example, do the following:

(1) Copy this example folder to some location outside of the pyne build folder
(2) run the compile script in source
(3) CHANGE THE FOLLOWING LINE IN pyne/pyne/spatialsolver.py
		sys.path.append("../pyne_transport/pyne/fortran/spatial_solvers_3d/source")
		Change the "../pyne_transport..." to the relative path from the example folder
		you copied to the pyne folder that you built pyne with.
			For example, if your pyne folder (the folder that you ran setup.py from) is
			located at : desktop/pyne, and the copy of this example folder that you copied
			is located at: desktop/example, then change the path line in spatialsolver.py to:
				"../pyne/pyne/fortran/spatial_solvers_3d/source"
(4) Rebuild pyne (python setup.py install)
(5) Run the example by either:
	(a) Running the python test_script.py file located in the example folder
			**the example folder you copied to the location outside of the
	(b) Run the following in an interactive python session:
			$import pyne.spatialsolver
			$ from dictionary_populate_test import populate, populate_with_warnings

			$ a = populate("DGFEM","DENSE")

			#Call solve to excecute.  Output written to file called fort.8
 			$ pyne.spatialsolver.solve(a)

			Here are more supported solver configurations:
				#ALL supported configurations without warnings
				#a = populate("AHOTN","LN")
				#a = populate("AHOTN","LL")
				#a = populate("AHOTN","NEFD")
				#a = populate("DGFEM","LD")
				#a = populate("DGFEM","DENSE")
				#a = populate("DGFEM","LAGRANGE")
				#Supported configurations with ALL POSSIBLE warnings
				#a = populate_with_warnings("AHOTN")
				#a = populate_with_warnings("DGFEM")

Alternately, you could use the old example that does not use the new pyne interface.  Until
the pyne interface is finished, this is probably the best way to go.

Requirements:

	(1) In addition to the pyne dependencies, you must have the following packages installed:
			(a) gfortran compiler
			(b) f2py

	In order to run this example, you must have run the compile script located at:
		pyne/fortran/spatial_solvers_3d/source/compile
	by typing the command:
		./compile

There are two options for running this example:

(1) Cd into the directory pyne/pyne/spatialsolver.  Run the test_script.py python file.  It generates a dictionary
	with all needed key-pair entries, then runs the NEFD solver.  To change dictionary entries, alter
	the dictionary_populate_test.py file.
		To change the solver, change the "solver" entry.  Note that so far, only AHOTN and DGFEM are supported
		To change the solver type, change the "solver_type" entry.  See below for full list of solver "flavors"
	Once the test_script.py file is run, it will generate a file called fort.8.  This file contains the
		output.  There is no way to get this data back into python (will me implemented soon).

(2) In a python session, run the following commands:
			**Note: you must be in the pyne/fortran/spatial_solver_3d directory.  Integration to pyne is not
					supported yet, so this (and pyne/fortran/spatial_solver_3d/soure) is the only directory that
					has access to the solvers python interface.

		# Imports
		$ import spatial_solver
		$ from dictionary_populate_test import populate, populate_with_warnings

		# Creating a dictionary, a, then populating it with the test script
		$ solver_dict = populate("AHOTN","LN")

		# Call the solver
		$ spatial_solver.solve(solver_dict)

				#ALL supported configurations without warnings
				#a = populate("AHOTN","LN")
				#a = populate("AHOTN","LL")
				#a = populate("AHOTN","NEFD")
				#a = populate("DGFEM","LD")
				#a = populate("DGFEM","DENSE")
				#a = populate("DGFEM","LAGRANGE")

					# Alternately, you could use the test dictionary generate script to generate a dictionary with only the
					# required key-pair values.  The rest will be overriden by the solver; according warnings will be printed.
						#Supported configurations with ALL POSSIBLE warnings
						#a = populate_with_warnings("AHOTN")
						#a = populate_with_warnings("DGFEM")

				# You can also create your own dictionary (more realistic, although the test dictionary script used above
				# does just that).  Below is a commented out example:

				# solver_dict = {
				#'solver':'AHOTN',
				#'solver_type':'NEFD', # OR LN OR LL
				#'spatial_order':1,
				#'spatial_method':0,
				#'angular_quadrature_order':4,  #WORKS
				#'angular_quadrature_type':1,
				#'nodes_xyz':[4,4,4],
				#'num_groups':1,
				#'num_materials':1,
				#'x_cells_widths':[0.25, 0.25, 0.25, 0.25],
		 		#'y_cells_widths':[0.25, 0.25, 0.25, 0.25],
				#'z_cells_widths':[0.25, 0.25, 0.25, 0.25],
				#'x_boundry_conditions':[2,2],
				#'y_boundry_conditions':[2,2],
				#'z_boundry_conditions':[2,2],
				#'material_id': [[[1]*4]*4]*4,
				#'quadrature_file':'quad_file',
				#'xs_file':'xs',
				#'source_input_file':'src_4.dat',
				#'bc_input_file':'bc_4.dat',
				#'flux_output_file':'phi_4.ahot',
				#'convergence_criterion':1.e-12,
				#'max_iterations':6000,
				#'moments_converged':0,
				#'converge_tolerence':1.e-10,
				#'ichk':0,
				#'ichk_tolerence':1.e-14
				#};
