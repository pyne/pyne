#!/usr/bin/python

import nose
from nose.tools import assert_equal, assert_almost_equal
from pyne import material
from pyne.material import Material, MaterialLibrary

# global variable filename for material libarary
filename = 'sample_output.h5m'
"""
test Lead
"""
def test_material1():
    # load the material library
    output_lib = MaterialLibrary()
    output_lib.from_hdf5(filename,'/materials')
    # create the material for the test
    Lead = Material({822040000: 0.013999999999999999, 822060000:
                     0.24100000000000002, 822070000: 0.22100000000000003, 822080000: 0.524})
    for material in output_lib.iteritems():
        if material[1].metadata['original_name'] == 'Lead':
            assert_almost_equal(material[1].comp, Lead.comp, places=4)

"""
test Nitrogen
"""
def test_material2():
    # load the material library
    output_lib = MaterialLibrary()
    output_lib.from_hdf5(filename,'/materials')
    # create the material for the test
    Nitrogen = Material({70140000: 0.99636, 70150000: 0.00364})
    for material in output_lib.iteritems():
        if material[1].metadata['original_name'] == 'Nitrogen':
            assert_almost_equal(material[1].comp, Nitrogen.comp, places=4)

"""
test Mercury
"""
def test_material3():
    # load the material library
    output_lib = MaterialLibrary()
    output_lib.from_hdf5(filename,'/materials')
    # create the material for the test
    Mercury = Material({801960000: 0.0015, 801980000: 0.09970000000000001, 801990000: 0.16870000000000002,
                   802000000: 0.231, 802010000: 0.1318, 802020000: 0.2986, 802040000: 0.0687})

    for material in output_lib.iteritems():
        if material[1].metadata['original_name'] == 'Mercury':
            assert_almost_equal(material[1].comp, Mercury.comp, places=4)

"""
test Steel, Stainless 321
"""
def test_material4():
    # load the material library
    output_lib = MaterialLibrary()
    output_lib.from_hdf5(filename,'/materials')
    # create the material for the test
    Steel = Material({60120000: 0.10992222222222224, 60130000: 0.0011888888888888893, 140280000: 0.10247000000000002, 140290000: 0.005205555555555556, 140300000: 0.0034355555555555563, 150310000: 0.11111111111111112, 160320000: 0.10554444444444445, 160330000: 0.0008333333333333334, 160340000: 0.004722222222222223, 160360000: 1.1111111111111113e-05, 220460000: 0.009166666666666668, 220470000: 0.008266666666666669, 220480000: 0.08191111111111112, 220490000: 0.006011111111111112, 220500000:
                 0.005755555555555556, 240500000: 0.004827777777777778, 240520000: 0.09309888888888891, 240530000: 0.010556666666666667, 240540000: 0.002627777777777779, 250550000: 0.11111111111111112, 260540000: 0.006494444444444446, 260560000: 0.10194888888888891, 260570000: 0.0023544444444444455, 260580000: 0.0003133333333333333, 280580000: 0.07564111111111112, 280600000: 0.029136666666666672, 280610000: 0.0012665555555555557, 280620000: 0.004038444444444445, 280640000: 0.0010283333333333336})

    for material in output_lib.iteritems():
        if material[1].metadata['original_name'] == 'Steel, Stainless 321':
            assert_almost_equal(material[1].comp, Steel.comp, places=4)

