'''AHOT Spatial Solver tests'''

#TODO:
#Add tests with supported warning configurations?
  #a = populate_with_warnings("AHOTN")
  #a = populate_with_warnings("DGFEM")
#Add result verifications to each solver?

import pyne.spatialsolver
from dictionary_populate_test import populate_simple, populate_simple_with_warnings, populate_intermediate_1

'''
def test_ahotn_ln():
  a = populate_simple("AHOTN","LN") 
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_ahotn_ll():
  a = populate_simple("AHOTN","LL")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)


def test_ahotn_nefd():
  a = populate_simple("AHOTN","NEFD")
  dict_results = pyne.spatialsolver.solve(a)
BAD

def test_dgfem_ld():
  a = populate_simple("DGFEM","LD")
  dict_results = pyne.spatialsolver.solve(a)

def test_dgfem_dense():
  a = populate_simple("DGFEM","DENSE")
  dict_results = pyne.spatialsolver.solve(a)

def test_dgfem_lagrange():
  a = populate_simple("DGFEM","LAGRANGE")
  dict_results = pyne.spatialsolver.solve(a)

def test_sct_step():
  a = populate_simple("SCTSTEP","anything")
  dict_results = pyne.spatialsolver.solve(a)
'''

def test_ahotn_ln_alternating():
  a = populate_intermediate_1("AHOTN", "LN")
  dict_results = pyne.spatialsolver.solve(a)

