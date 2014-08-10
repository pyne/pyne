'''AHOT Spatial Solver tests'''

#TODO:
#Add tests with supported warning configurations?
  #a = populate_with_warnings("AHOTN")
  #a = populate_with_warnings("DGFEM")

import numpy as np

import pyne.spatialsolver
from dictionary_populate_test import populate_simple, populate_simple_with_warnings, populate_intermediate_1

def test_ahotn_ln():
  a = populate_simple("AHOTN","LN") 
  dict_results = pyne.spatialsolver.solve(a)
  rounded_flux = np.around(dict_results['flux'], decimals=5)

  correct_flux =  [[[ 3.52650199,  3.09260257,  3.09260257,  3.52650199],
                    [ 3.09260257,  2.73209732,  2.73209732,  3.09260257],
                    [ 3.09260257,  2.73209732,  2.73209732,  3.09260257],
                    [ 3.52650199,  3.09260257,  3.09260257,  3.52650199],],

                   [[ 2.89021832,  2.61284811,  2.61284811,  2.89021832],
                    [ 2.61284811,  2.38571678,  2.38571678,  2.61284811],
                    [ 2.61284811,  2.38571678,  2.38571678,  2.61284811],
                    [ 2.89021832,  2.61284811,  2.61284811,  2.89021832],],

                   [[ 2.89021832,  2.61284811,  2.61284811,  2.89021832],
                    [ 2.61284811,  2.38571678,  2.38571678,  2.61284811],
                    [ 2.61284811,  2.38571678,  2.38571678,  2.61284811],
                    [ 2.89021832,  2.61284811,  2.61284811,  2.89021832],],

                   [[ 3.52650199,  3.09260257,  3.09260257,  3.52650199],
                    [ 3.09260257,  2.73209732,  2.73209732,  3.09260257],
                    [ 3.09260257,  2.73209732,  2.73209732,  3.09260257],
                    [ 3.52650199,  3.09260257,  3.09260257,  3.52650199]]]
  correct_flux_rounded = np.around(correct_flux, decimals=5)
  if (rounded_flux==correct_flux_rounded).all():
    print("flux's are equal!")
  else:
    raise AssertionError("Flux outputs are not equal for ahotn-ln example.  Check system setup.")
  print(dict_results)

def test_ahotn_ll():
  a = populate_simple("AHOTN","LL")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)


def test_ahotn_nefd():
  a = populate_simple("AHOTN","NEFD")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_dgfem_ld():
  a = populate_simple("DGFEM","LD")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_dgfem_dense():
  a = populate_simple("DGFEM","DENSE")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_dgfem_lagrange():
  a = populate_simple("DGFEM","LAGRANGE")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_sct_step():
  a = populate_simple("SCTSTEP","anything")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)

def test_ahotn_ln_alternating():
  a = populate_intermediate_1("AHOTN", "LN")
  dict_results = pyne.spatialsolver.solve(a)
  print(dict_results)


