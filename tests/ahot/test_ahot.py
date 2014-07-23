'''AHOT Spatial Solver tests'''

import pyne.spatialsolver
from dictionary_populate_test import populate, populate_with_warnings

def test_ahotn_ln():
  a = populate("AHOTN","LN")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

#NOTE: Before multiple tests can be ran, a python callback needs to be added in the wrapped
#solver codes so that the following tests dont begin until the prior test has terminated in
#the fortran runtime environment.

def test_ahotn_ll():
  a = populate("AHOTN","LL")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

#def test_ahotn_nefd():
#  a = populate("AHOTN","NEFD")
#  dict_results = pyne.spatialsolver.solve(a)

#def test_dgfem_ld():
#  a = populate("DGFEM","LD")
#  dict_results = pyne.spatialsolver.solve(a)

#def test_dgfem_dense():
#  a = populate("DGFEM","DENSE")
#  dict_results = pyne.spatialsolver.solve(a)

#def test_dgfem_lagrange():
#  a = populate("DGFEM","LAGRANGE")
#  dict_results = pyne.spatialsolver.solve(a)

#TODO:
#Add tests with supported warning configurations?
  #a = populate_with_warnings("AHOTN")
  #a = populate_with_warnings("DGFEM")
#Add result verifications to each solver?
