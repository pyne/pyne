import spatial_solver
from dictionary_populate_test import populate, populate_with_warnings

a = populate()
#a = populate_with_warnings()
spatial_solver.solve(a)
