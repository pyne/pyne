'''This is a sample test script to quickly verify the integrity of the solve as work is being done.  To run, 
simply run this file in python. '''

import spatial_solver
from dictionary_populate_test import populate, populate_with_warnings

a = populate("DGFEM","DENSE")

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

#Call solve to excecute.  Output written to file called fort.8
spatial_solver.solve(a)
