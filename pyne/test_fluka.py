#!/usr/bin/python

import os
import fluka

import nose.tools
from nose.tools import assert_almost_equal, assert_equal, assert_true, \
    assert_not_equal, assert_false, assert_raises
from nose.plugins.skip import SkipTest

# Mesh specific imports
try:
    from itaps import iMesh

    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest

    HAVE_PYTAPS = False
    pass

from pyne.mesh import Mesh, StatMesh, MeshError
from numpy.testing import assert_array_equal

# Test USR and USRBIN classes

def test_single_usrbin():
    """Test a usrbin file containing a single tally.
    """

    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_single.lis")
    expected_h5m = os.path.join(thisdir, "single_n.h5m")
    expected_sm = Mesh(mesh=expected_h5m, structured=True)

    usrbin_object = fluka.UsrbinFile(usrbin_file)

    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, -1.0, 1.0, 3.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -1.0, 1.0, 3.0]
    assert_equal(usrbin_object.tally['single_n'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['single_n'].name, 'single_n')
    assert_equal(usrbin_object.tally['single_n'].particle, '8')
    assert_equal(usrbin_object.tally['single_n'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['single_n'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['single_n'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
                            1.0077E-02, 5.1336E-03, 2.3750E-03, 2.0697E-02, 
                            8.7204E-03, 3.1160E-03, 1.1974E-02, 5.8080E-03, 
                            2.8636E-03, 2.1603E-02, 8.5671E-03, 4.0857E-03, 
                            1.4707E-01, 1.5037E-02, 4.1191E-03, 2.2340E-02, 
                            8.3678E-03, 2.9380E-03, 1.1004E-02, 5.2263E-03,
                            2.4421E-03, 2.2533E-02, 8.9809E-03, 2.8157E-03, 
                            1.1616E-02, 5.8871E-03, 2.6331E-03
                            ]
    expected_error_data = [
                            6.5873E+00, 3.9876E+00, 6.6591E+00, 5.1410E+00, 
                            6.3925E+00, 7.3896E+00, 4.4964E+00, 1.3437E+01, 
                            1.2783E+01, 3.9902E+00, 1.8128E+00, 7.2382E+00, 
                            4.9879E-01, 6.0769E+00, 5.5579E+00, 3.0970E+00, 
                            3.0329E+00, 9.5369E+00, 6.1279E+00, 1.0225E+01,
                            1.5106E+01, 4.1769E+00, 5.8046E+00, 1.3207E+01, 
                            5.9847E+00, 1.2240E+01, 1.5934E+01
                            ]

    for i,v_e in enumerate(usrbin_object.tally['single_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['single_n'].mesh.getTagHandle('part_data_8')[v_e]
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i,v_e in enumerate(usrbin_object.tally['single_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['single_n'].mesh.getTagHandle('error_data_8')[v_e]
        expected = expected_error_data[i]
        assert_equal(read, expected)

#def test_multiple_usrbin():

test_single_usrbin()
