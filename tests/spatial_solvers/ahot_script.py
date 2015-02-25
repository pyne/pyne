'''AHOT Spatial Solver tests'''

#TODO:
#Add tests with supported warning configurations?
  #a = populate_with_warnings("AHOTN")
  #a = populate_with_warnings("DGFEM")

import numpy as np

import pyne.spatialsolver
from .dictionary_populate_test import populate_simple, populate_simple_with_warnings, populate_intermediate_1

def test_ahotn_ln():
    a = populate_simple("AHOTN","LN")
    dict_results = pyne.spatialsolver.solve(a)

    if(dict_results['success'] == 0):
        raise AssertionError("Error: " + dict_results['error_msg'])


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

def test_dgfem_ld():
    a = populate_simple("DGFEM","LD")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[3.540511,  3.104096,  3.104096,  3.540511],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.540511,  3.104096,  3.104096,  3.540511],],

                    [[2.899079,  2.620152,  2.620152,  2.899079],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.899079,  2.620152,  2.620152,  2.899079],],

                    [[2.899079,  2.620152,  2.620152,  2.899079],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.899079,  2.620152,  2.620152,  2.899079],],

                    [[3.540511,  3.104096,  3.104096,  3.540511],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.540511,  3.104096,  3.104096,  3.540511]]]

    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_dgfem_dense():
    a = populate_simple("DGFEM","DENSE")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[3.540511,  3.104096,  3.104096,  3.540511],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.540511,  3.104096,  3.104096,  3.540511],],

                    [[2.899079,  2.620152,  2.620152,  2.899079],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.899079,  2.620152,  2.620152,  2.899079],],

                    [[2.899079,  2.620152,  2.620152,  2.899079],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.620152,  2.383940,  2.383940,  2.620152],
                    [2.899079,  2.620152,  2.620152,  2.899079],],

                    [[3.540511,  3.104096,  3.104096,  3.540511],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.104096,  2.730554,  2.730554,  3.104096],
                    [3.540511,  3.104096,  3.104096,  3.540511]]]

    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_dgfem_lagrange():
    a = populate_simple("DGFEM","LAGRANGE")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[3.536038,  3.096808,  3.096808,  3.536038],
                    [3.096808,  2.732475,  2.732475,  3.096808],
                    [3.096808,  2.732475,  2.732475,  3.096808],
                    [3.536038,  3.096808,  3.096808,  3.536038],],

                    [[2.896267,  2.615275,  2.615275,  2.896267],
                    [2.615275,  2.385484,  2.385484,  2.615275],
                    [2.615275,  2.385484,  2.385484,  2.615275],
                    [2.896267,  2.615275,  2.615275,  2.896267],],

                    [[2.896267,  2.615275,  2.615275,  2.896267],
                    [2.615275,  2.385484,  2.385484,  2.615275],
                    [2.615275,  2.385484,  2.385484,  2.615275],
                    [2.896267,  2.615275,  2.615275,  2.896267],],

                    [[3.536038,  3.096808,  3.096808,  3.536038],
                    [3.096808,  2.732475,  2.732475,  3.096808],
                    [3.096808,  2.732475,  2.732475,  3.096808],
                    [3.536038,  3.096808,  3.096808,  3.536038]]]

    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_sct_step():
    a = populate_simple("SCTSTEP","anything")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[3.273572,  2.948301,  2.948502,  3.291909],
                    [2.811363,  2.464789,  2.468086,  2.813676],
                    [2.921249,  2.576771,  2.593078,  2.919847],
                    [3.138840,  2.784381,  2.785791,  3.139999],],

                    [[2.466767,  2.188464,  2.191274,  2.465690],
                    [2.168904,  1.883310,  1.884325,  2.169292],
                    [2.181507,  1.891052,  1.895120,  2.178766],
                    [2.438198,  2.161378,  2.161873,  2.438270],],

                    [[2.429940,  2.143983,  2.143274,  2.427243],
                    [2.144259,  1.849312,  1.848996,  2.143790],
                    [2.142347,  1.843699,  1.841852,  2.140937],
                    [2.425510,  2.142483,  2.142357,  2.425371],],

                    [[3.091479,  2.729188,  2.728940,  3.091578],
                    [2.727627,  2.366091,  2.365882,  2.727488],
                    [2.726782,  2.365203,  2.364727,  2.726503],
                    [3.087793,  2.725209,  2.725085,  3.087700]]]

    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_ahotn_ln_alternating():
    a = populate_intermediate_1("AHOTN", "LN")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.302715,  2.230236,  1.817902,  2.952883],
                    [2.230236,  1.292285,  1.620001,  1.817902],
                    [1.817902,  1.620001,  1.292285,  2.230236],
                    [2.952883,  1.817902,  2.230236,  2.302715],],

                    [[2.289555,  1.443020,  1.762396,  1.811167],
                    [1.443020,  1.283541,  1.038793,  1.762396],
                    [1.762396,  1.038793,  1.283541,  1.443020],
                    [1.811167,  1.762396,  1.443020,  2.289555],],

                    [[1.811167,  1.762396,  1.443020,  2.289555],
                    [1.762396,  1.038793,  1.283541,  1.443020],
                    [1.443020,  1.283541,  1.038793,  1.762396],
                    [2.289555,  1.443020,  1.762396,  1.811167],],

                    [[2.952883,  1.817902,  2.230236,  2.302715],
                    [1.817902,  1.620001,  1.292285,  2.230236],
                    [2.230236,  1.292285,  1.620001,  1.817902],
                    [2.302715,  2.230236,  1.817902,  2.952883]]]


    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_ahotn_ll_alternating():
    a = populate_intermediate_1("AHOTN", "LL")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[1.685745,  2.553723,  1.536011,  2.927681],
                    [2.553723,  1.444449,  2.291995,  1.536011],
                    [1.536011,  2.291995,  1.444449,  2.553723],
                    [2.927681,  1.536011,  2.553723,  1.685745],],

                    [[2.351537,  1.291496,  2.101697,  1.394040],
                    [1.291496,  1.922387,  1.225058,  2.101697],
                    [2.101697,  1.225058,  1.922387,  1.291496],
                    [1.394040,  2.101697,  1.291496,  2.351537],],

                    [[1.394040,  2.101697,  1.291496,  2.351537],
                    [2.101697,  1.225058,  1.922387,  1.291496],
                    [1.291496,  1.922387,  1.225058,  2.101697],
                    [2.351537,  1.291496,  2.101697,  1.394040],],

                    [[2.927681,  1.536011,  2.553723,  1.685745],
                    [1.536011,  2.291995,  1.444449,  2.553723],
                    [2.553723,  1.444449,  2.291995,  1.536011],
                    [1.685745,  2.553723,  1.536011,  2.927681]]]


    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")


def test_ahotn_nefd_alternating():
    a = populate_intermediate_1("AHOTN", "NEFD")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.320847,  2.193170,  1.836823,  2.923995],
                    [2.193170,  1.310507,  1.568554,  1.836823],
                    [1.836823,  1.568554,  1.310507,  2.193170],
                    [2.923995,  1.836823,  2.193170,  2.320847],],

                    [[2.266863,  1.456056,  1.732060,  1.824538],
                    [1.456056,  1.241531,  1.049696,  1.732060],
                    [1.732060,  1.049696,  1.241531,  1.456056],
                    [1.824538,  1.732060,  1.456056,  2.266863],],

                    [[1.824538,  1.732060,  1.456056,  2.266863],
                    [1.732060,  1.049696,  1.241531,  1.456056],
                    [1.456056,  1.241531,  1.049696,  1.732060],
                    [2.266863,  1.456056,  1.732060,  1.824538],],

                    [[2.923995,  1.836823,  2.193170,  2.320847],
                    [1.836823,  1.568554,  1.310507,  2.193170],
                    [2.193170,  1.310507,  1.568554,  1.836823],
                    [2.320847,  2.193170,  1.836823,  2.923995]]]



    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_dgfem_ld_alternating():
    a = populate_intermediate_1("DGFEM", "LD")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.420725,  2.104426,  1.900442,  2.889886],
                    [2.104426,  1.299636,  1.433389,  1.900442],
                    [1.900442,  1.433389,  1.299636,  2.104426],
                    [2.889886,  1.900442,  2.104426,  2.420725],],

                    [[2.224013,  1.498666,  1.647904,  1.894524],
                    [1.498666,  1.119896,  1.039153,  1.647904],
                    [1.647904,  1.039153,  1.119896,  1.498666],
                    [1.894524,  1.647904,  1.498666,  2.224013],],

                    [[1.894524,  1.647904,  1.498666,  2.224013],
                    [1.647904,  1.039153,  1.119896,  1.498666],
                    [1.498666,  1.119896,  1.039153,  1.647904],
                    [2.224013,  1.498666,  1.647904,  1.894524],],

                    [[2.889886,  1.900442,  2.104426,  2.420725],
                    [1.900442,  1.433389,  1.299636,  2.104426],
                    [2.104426,  1.299636,  1.433389,  1.900442],
                    [2.420725,  2.104426,  1.900442,  2.889886]]]



    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_dgfem_dense_alternating():
    a = populate_intermediate_1("DGFEM", "DENSE")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.420725,  2.104426,  1.900442,  2.889886],
                    [2.104426,  1.299636,  1.433389,  1.900442],
                    [1.900442,  1.433389,  1.299636,  2.104426],
                    [2.889886,  1.900442,  2.104426,  2.420725],],

                    [[2.224013,  1.498666,  1.647904,  1.894524],
                    [1.498666,  1.119896,  1.039153,  1.647904],
                    [1.647904,  1.039153,  1.119896,  1.498666],
                    [1.894524,  1.647904,  1.498666,  2.224013],],

                    [[1.894524,  1.647904,  1.498666,  2.224013],
                    [1.647904,  1.039153,  1.119896,  1.498666],
                    [1.498666,  1.119896,  1.039153,  1.647904],
                    [2.224013,  1.498666,  1.647904,  1.894524],],

                    [[2.889886,  1.900442,  2.104426,  2.420725],
                    [1.900442,  1.433389,  1.299636,  2.104426],
                    [2.104426,  1.299636,  1.433389,  1.900442],
                    [2.420725,  2.104426,  1.900442,  2.889886]]]


    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_dgfem_lagrange_alternating():
    a = populate_intermediate_1("DGFEM", "LAGRANGE")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.403548,  2.135009,  1.885348,  2.906123],
                    [2.135009,  1.300693,  1.469197,  1.885348],
                    [1.885348,  1.469197,  1.300693,  2.135009],
                    [2.906123,  1.885348,  2.135009,  2.403548],],

                    [[2.241881,  1.486578,  1.673153,  1.882209],
                    [1.486578,  1.145347,  1.036189,  1.673153],
                    [1.673153,  1.036189,  1.145347,  1.486578],
                    [1.882209,  1.673153,  1.486578,  2.241881],],

                    [[1.882209,  1.673153,  1.486578,  2.241881],
                    [1.673153,  1.036189,  1.145347,  1.486578],
                    [1.486578,  1.145347,  1.036189,  1.673153],
                    [2.241881,  1.486578,  1.673153,  1.882209],],

                    [[2.906123,  1.885348,  2.135009,  2.403548],
                    [1.885348,  1.469197,  1.300693,  2.135009],
                    [2.135009,  1.300693,  1.469197,  1.885348],
                    [2.403548,  2.135009,  1.885348,  2.906123]]]

    correct_flux_rounded = np.around(correct_flux, decimals=4)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")

def test_sct_step_alternating():
    a = populate_intermediate_1("SCTSTEP", "anything")
    dict_results = pyne.spatialsolver.solve(a)
    rounded_flux = np.around(dict_results['flux'], decimals=4)

    correct_flux =  [[[2.103727,  2.129333,  1.775806,  2.709218],
                    [1.984849,  1.172710,  1.337597,  1.664623],
                    [1.757312,  1.459605,  1.282230,  2.107971],
                    [2.551582,  1.644416,  1.966496,  1.996478],],

                    [[1.909362,  1.216011,  1.443766,  1.521228],
                    [1.198507,  .8426090,  .7858172,  1.423269],
                    [1.435932,  .7960783,  .8584189,  1.209827],
                    [1.500600,  1.417286,  1.194468,  1.887075],],

                    [[1.497664,  1.410221,  1.186999,  1.881503],
                    [1.408052,  .7672912,  .8230592,  1.185632],
                    [1.186346,  .8224311,  .7656347,  1.407697],
                    [1.878868,  1.184635,  1.406690,  1.494015],],

                    [[2.519203,  1.608783,  1.927761,  1.963608],
                    [1.608023,  1.265341,  1.108607,  1.927101],
                    [1.9271,  1.108730,  1.265047,  1.608085],
                    [1.962463,  1.926423,  1.607454,  2.518035]]]



    correct_flux_rounded = np.around(correct_flux, decimals=4)
    print(correct_flux_rounded)
    print(rounded_flux)
    if (rounded_flux==correct_flux_rounded).all():
        print("flux's are equal!")
    else:
        raise AssertionError("Flux outputs are not equal for ahotn-nefd example.  Check system setup.")
