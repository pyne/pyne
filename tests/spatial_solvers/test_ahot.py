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
  rounded_flux = np.around(dict_results['flux'], decimals=4)

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
  correct_flux_rounded = np.around(correct_flux, decimals=4)
  if (rounded_flux==correct_flux_rounded).all():
    print("flux's are equal!")
  else:
    raise AssertionError("Flux outputs are not equal for ahotn-ln example.  Check system setup.")
  #print(dict_results)

def test_ahotn_ll():
  a = populate_simple("AHOTN","LL")
  dict_results = pyne.spatialsolver.solve(a)
  rounded_flux = np.around(dict_results['flux'], decimals=4)

  correct_flux =  [[[3.478172,  3.178575,  3.178575,  3.478172],
                    [3.178575,  3.069365,  3.069365,  3.178575],
                    [3.178575,  3.069365,  3.069365,  3.178575],
                    [3.478172,  3.178575,  3.178575,  3.478172],],

                   [[2.855123,  2.666765,  2.666765,  2.855123],
                    [2.666765,  2.599730,  2.599730,  2.666765],
                    [2.666765,  2.599730,  2.599730,  2.666765],
                    [2.855123,  2.666765,  2.666765,  2.855123],],

                   [[2.855123,  2.666765,  2.666765,  2.855123],
                    [2.666765,  2.599730,  2.599730,  2.666765],
                    [2.666765,  2.599730,  2.599730,  2.666765],
                    [2.855123,  2.666765,  2.666765,  2.855123],],

                   [[3.478172,  3.178575,  3.178575,  3.478172],
                    [3.178575,  3.069365,  3.069365,  3.178575],
                    [3.178575,  3.069365,  3.069365,  3.178575],
                    [3.478172,  3.178575,  3.178575,  3.478172]]]

  correct_flux_rounded = np.around(correct_flux, decimals=4)
  if (rounded_flux==correct_flux_rounded).all():
    print("flux's are equal!")
  else:
    raise AssertionError("Flux outputs are not equal for ahotn-ll example.  Check system setup.")
  #print(dict_results)


def test_ahotn_nefd():
  a = populate_simple("AHOTN","NEFD")
  dict_results = pyne.spatialsolver.solve(a)
  rounded_flux = np.around(dict_results['flux'], decimals=4)

  correct_flux =  [[[3.524073,  3.091501,  3.091501,  3.524073],
                    [3.091501,  2.734906,  2.734906,  3.091501],
                    [3.091501,  2.734906,  2.734906,  3.091501],
                    [3.524073,  3.091501,  3.091501,  3.524073],],


                   [[2.888798,  2.612178,  2.612178,  2.888798],
                    [2.612178,  2.387341,  2.387341,  2.612178],
                    [2.612178,  2.387341,  2.387341,  2.612178],
                    [2.888798,  2.612178,  2.612178,  2.888798],],


                   [[2.888798,  2.612178,  2.612178,  2.888798],
                    [2.612178,  2.387341,  2.387341,  2.612178],
                    [2.612178,  2.387341,  2.387341,  2.612178],
                    [2.888798,  2.612178,  2.612178,  2.888798],],


                   [[3.524073,  3.091501,  3.091501,  3.524073],
                    [3.091501,  2.734906,  2.734906,  3.091501],
                    [3.091501,  2.734906,  2.734906,  3.091501],
                    [3.524073,  3.091501,  3.091501,  3.524073]]]

  correct_flux_rounded = np.around(correct_flux, decimals=4)
  if (rounded_flux==correct_flux_rounded).all():
    print("flux's are equal!")
  else:
    raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")
  #print(dict_results)

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


