'''This is a sample test script to quickly verify the integrity of the solve as work is being done.  To run, 
simply run this file in python. '''

import spatial_solver
from dictionary_populate_test import populate, populate_with_warnings

a = populate("DGFEM","DENSE")
#a = populate_with_warnings()
spatial_solver.solve(a)
