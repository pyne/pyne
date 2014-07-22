'''AHOT Spatial Solver tests'''

import pyne.spatialsolver
from dictionary_populate_test import populate, populate_with_warnings

def test_ahotn_ln():
  a = populate("AHOTN","LN")
  pyne.spatialsolver.solve(a)

#NOTE: Before multiple tests can be ran, a python callback needs to be added in the wrapped
#solver codes so that the following tests dont begin until the prior test has terminated in
#the fortran runtime environment.

#def test_ahotn_ll():
#  a = populate("AHOTN","LL")
#  pyne.spatialsolver.solve(a)

#def test_ahotn_nefd():
#  a = populate("AHOTN","NEFD")
#  pyne.spatialsolver.solve(a)

#def test_dgfem_ld():
#  a = populate("DGFEM","LD")
#  pyne.spatialsolver.solve(a)

#def test_dgfem_dense():
#  a = populate("DGFEM","DENSE")
#  pyne.spatialsolver.solve(a)

#def test_dgfem_lagrange():
#  a = populate("DGFEM","LAGRANGE")
#  pyne.spatialsolver.solve(a)

#TODO:
#Add tests with supported warning configurations?
  #a = populate_with_warnings("AHOTN")
  #a = populate_with_warnings("DGFEM")
#Add result verifications to each solver?
