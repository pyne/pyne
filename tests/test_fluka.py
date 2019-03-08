#!/usr/bin/python

import os
import nose.tools
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal

from pyne import fluka
# Mesh specific imports
from pyne.mesh import Mesh, StatMesh, MeshError, HAVE_PYMOAB

# Test Usrbin and UsrbinTally classes


def test_single_usrbin():
    """Test a usrbin file containing a single tally.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_single.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)

    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally['single_n'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['single_n'].name, 'single_n')
    assert_equal(usrbin_object.tally['single_n'].particle, '8')
    assert_equal(usrbin_object.tally['single_n'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['single_n'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['single_n'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.0984E-02, 4.1051E-03, 1.0636E-03, 2.1837E-02, 5.5610E-03, 1.9119E-03,
        1.0971E-02, 3.3943E-03, 1.2456E-03, 1.6615E-02, 2.9501E-03, 7.4597E-04,
        1.0395E-01, 6.1186E-03, 1.4997E-03, 1.7421E-02, 3.0824E-03, 7.3878E-04,
        1.8097E-02, 5.2532E-03, 2.1572E-03, 1.0465E-01, 6.2611E-03, 1.8829E-03,
        1.7323E-02, 5.5092E-03, 2.1418E-03
    ]
    expected_error_data = [
        5.0179E+00, 1.6521E+01, 1.3973E+01, 4.2025E+00, 8.1766E+00, 1.1465E+01,
        7.2005E+00, 1.0479E+01, 1.5640E+01, 5.5994E+00, 1.3275E+01, 2.7617E+01,
        7.3788E-01, 6.7200E+00, 1.9092E+01, 7.3670E+00, 1.3018E+01, 2.8866E+01,
        5.7221E+00, 1.5916E+01, 2.6001E+01, 8.3490E-01, 1.6715E+01, 1.2759E+01,
        5.0763E+00, 1.1420E+01, 1.0040E+01
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['single_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['single_n'].part_data_tag[v_e]
        assert_equal(
            usrbin_object.tally['single_n'].part_data_tag.name, 'part_data_8')
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['single_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['single_n'].error_data_tag[v_e]
        assert_equal(
            usrbin_object.tally['single_n'].error_data_tag.name, 'error_data_8')
        expected = expected_error_data[i]
        assert_equal(read, expected)


def test_multiple_usrbin():
    """Test a usrbin file containing multiple (two) tallies.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_multiple.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)

    # Tally #1:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally['multi_p'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['multi_p'].name, 'multi_p')
    assert_equal(usrbin_object.tally['multi_p'].particle, '7')
    assert_equal(usrbin_object.tally['multi_p'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['multi_p'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['multi_p'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        7.5083E-04, 1.7570E-04, 3.3361E-05, 1.1232E-03, 3.4735E-04, 1.5816E-04,
        6.2264E-04, 2.3071E-04, 8.3469E-05, 1.6700E-03, 4.1785E-04, 7.6990E-05,
        3.3842E-03, 9.2931E-04, 2.4958E-04, 1.0121E-03, 2.7993E-04, 6.1043E-05,
        7.7401E-04, 3.2480E-04, 9.3145E-06, 1.4245E-03, 4.3352E-04, 1.7392E-04,
        7.3166E-04, 2.4210E-04, 1.4804E-04
    ]
    expected_error_data = [
        2.2149E+01, 7.4509E+01, 1.0000E+02, 2.4621E+01, 4.6383E+01, 3.3621E+01,
        2.1616E+01, 7.5885E+01, 1.0000E+02, 2.0067E+01, 3.3654E+01, 6.1265E+01,
        1.8407E+01, 1.6239E+01, 5.2119E+01, 1.5791E+01, 3.8452E+01, 1.0000E+02,
        7.6577E+00, 3.5290E+01, 1.0000E+02, 8.3702E+00, 5.3283E+01, 6.2602E+01,
        1.1655E+01, 6.2289E+01, 6.7541E+01
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['multi_p'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['multi_p'].part_data_tag[v_e]
        assert_equal(
            usrbin_object.tally['multi_p'].part_data_tag.name, 'part_data_7')
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['multi_p'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['multi_p'].error_data_tag[v_e]
        assert_equal(
            usrbin_object.tally['multi_p'].error_data_tag.name, 'error_data_7')
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #2:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally['multi_n'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['multi_n'].name, 'multi_n')
    assert_equal(usrbin_object.tally['multi_n'].particle, '8')
    assert_equal(usrbin_object.tally['multi_n'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['multi_n'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['multi_n'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.0984E-02, 4.1051E-03, 1.0636E-03, 2.1837E-02, 5.5610E-03, 1.9119E-03,
        1.0971E-02, 3.3943E-03, 1.2456E-03, 1.6615E-02, 2.9501E-03, 7.4597E-04,
        1.0395E-01, 6.1186E-03, 1.4997E-03, 1.7421E-02, 3.0824E-03, 7.3878E-04,
        1.8097E-02, 5.2532E-03, 2.1572E-03, 1.0465E-01, 6.2611E-03, 1.8829E-03,
        1.7323E-02, 5.5092E-03, 2.1418E-03
    ]
    expected_error_data = [
        5.0179E+00, 1.6521E+01, 1.3973E+01, 4.2025E+00, 8.1766E+00, 1.1465E+01,
        7.2005E+00, 1.0479E+01, 1.5640E+01, 5.5994E+00, 1.3275E+01, 2.7617E+01,
        7.3788E-01, 6.7200E+00, 1.9092E+01, 7.3670E+00, 1.3018E+01, 2.8866E+01,
        5.7221E+00, 1.5916E+01, 2.6001E+01, 8.3490E-01, 1.6715E+01, 1.2759E+01,
        5.0763E+00, 1.1420E+01, 1.0040E+01
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['multi_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['multi_n'].part_data_tag[v_e]
        assert usrbin_object.tally['multi_n'].part_data_tag.name == 'part_data_8'
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['multi_n'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['multi_n'].error_data_tag[v_e]
        assert usrbin_object.tally['multi_n'].error_data_tag.name == 'error_data_8'
        expected = expected_error_data[i]
        assert_equal(read, expected)


def test_degenerate_usrbin():
    """Test usrbin file containing tallies with different number of bins in 
    each direction.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_degenerate.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)

    # Tally #1:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, 0.0, 3.0]
    expected_zbounds = [-3.0, 0.0]
    assert_equal(usrbin_object.tally['degen1'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['degen1'].name, 'degen1')
    assert_equal(usrbin_object.tally['degen1'].particle, '8')
    assert_equal(usrbin_object.tally['degen1'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['degen1'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['degen1'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        3.5279E-02, 4.7334E-03, 1.4458E-03,
        3.6242E-02, 4.6521E-03, 1.5292E-03
    ]
    expected_error_data = [
        1.2016E+00, 6.4313E+00, 7.7312E+00,
        2.0235E+00, 9.4199E+00, 8.0514E+00
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['degen1'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen1'].part_data_tag[v_e]
        assert usrbin_object.tally['degen1'].part_data_tag.name == 'part_data_8'
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['degen1'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen1'].error_data_tag[v_e]
        assert usrbin_object.tally['degen1'].error_data_tag.name == 'error_data_8'
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #2:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 1.5, 6.0]
    expected_ybounds = [-3.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally['degen2'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['degen2'].name, 'degen2')
    assert_equal(usrbin_object.tally['degen2'].particle, '8')
    assert_equal(usrbin_object.tally['degen2'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['degen2'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['degen2'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.1543E-02, 2.0295E-03, 3.2603E-02,
        1.4229E-03, 3.3492E-02, 2.7923E-03
    ]
    expected_error_data = [
        2.7321E+00, 5.2342E+00, 7.4679E-01,
        4.2862E+00, 1.3090E+00, 1.4151E+01
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['degen2'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen2'].part_data_tag[v_e]
        assert usrbin_object.tally['degen2'].part_data_tag.name == 'part_data_8'
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['degen2'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen2'].error_data_tag[v_e]
        assert usrbin_object.tally['degen2'].error_data_tag.name == 'error_data_8'
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #3:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -1.5, 0.0]
    assert_equal(usrbin_object.tally['degen3'].coord_sys, 'Cartesian')
    assert_equal(usrbin_object.tally['degen3'].name, 'degen3')
    assert_equal(usrbin_object.tally['degen3'].particle, '8')
    assert_equal(usrbin_object.tally['degen3'].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally['degen3'].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally['degen3'].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        5.8037E-03, 1.3260E-02, 5.6046E-03,
        7.9677E-03, 4.3111E-02, 8.1349E-03
    ]
    expected_error_data = [
        6.1913E+00, 2.3684E+00, 4.6124E+00,
        3.2523E+00, 1.3714E+00, 4.3161E+00
    ]

    for i, v_e in enumerate(
            usrbin_object.tally['degen3'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen3'].part_data_tag[v_e]
        assert usrbin_object.tally['degen3'].part_data_tag.name == 'part_data_8'
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
            usrbin_object.tally['degen3'].structured_iterate_hex("zyx")):
        read = usrbin_object.tally['degen3'].error_data_tag[v_e]
        assert usrbin_object.tally['degen3'].error_data_tag.name == 'error_data_8'
        expected = expected_error_data[i]
        assert_equal(read, expected)

# test file writing to catch upstream changes in mesh


def test_mesh_write():
    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_single.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)
    data = usrbin_object.tally['single_n']
    data.write_hdf5("test_fluka_data.h5m")
