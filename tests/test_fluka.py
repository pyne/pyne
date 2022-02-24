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
    """Test a usrbin file containing a single tally."""

    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_single.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)

    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally["single_n"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["single_n"].name, "single_n")
    assert_equal(usrbin_object.tally["single_n"].particle, "8")
    assert_equal(usrbin_object.tally["single_n"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["single_n"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["single_n"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.0984e-02,
        4.1051e-03,
        1.0636e-03,
        2.1837e-02,
        5.5610e-03,
        1.9119e-03,
        1.0971e-02,
        3.3943e-03,
        1.2456e-03,
        1.6615e-02,
        2.9501e-03,
        7.4597e-04,
        1.0395e-01,
        6.1186e-03,
        1.4997e-03,
        1.7421e-02,
        3.0824e-03,
        7.3878e-04,
        1.8097e-02,
        5.2532e-03,
        2.1572e-03,
        1.0465e-01,
        6.2611e-03,
        1.8829e-03,
        1.7323e-02,
        5.5092e-03,
        2.1418e-03,
    ]
    expected_error_data = [
        5.0179e00,
        1.6521e01,
        1.3973e01,
        4.2025e00,
        8.1766e00,
        1.1465e01,
        7.2005e00,
        1.0479e01,
        1.5640e01,
        5.5994e00,
        1.3275e01,
        2.7617e01,
        7.3788e-01,
        6.7200e00,
        1.9092e01,
        7.3670e00,
        1.3018e01,
        2.8866e01,
        5.7221e00,
        1.5916e01,
        2.6001e01,
        8.3490e-01,
        1.6715e01,
        1.2759e01,
        5.0763e00,
        1.1420e01,
        1.0040e01,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["single_n"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["single_n"].part_data_tag[v_e]
        assert_equal(usrbin_object.tally["single_n"].part_data_tag.name, "part_data_8")
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["single_n"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["single_n"].error_data_tag[v_e]
        assert_equal(
            usrbin_object.tally["single_n"].error_data_tag.name, "error_data_8"
        )
        expected = expected_error_data[i]
        assert_equal(read, expected)


def test_multiple_usrbin():
    """Test a usrbin file containing multiple (two) tallies."""

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
    assert_equal(usrbin_object.tally["multi_p"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["multi_p"].name, "multi_p")
    assert_equal(usrbin_object.tally["multi_p"].particle, "7")
    assert_equal(usrbin_object.tally["multi_p"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["multi_p"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["multi_p"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        7.5083e-04,
        1.7570e-04,
        3.3361e-05,
        1.1232e-03,
        3.4735e-04,
        1.5816e-04,
        6.2264e-04,
        2.3071e-04,
        8.3469e-05,
        1.6700e-03,
        4.1785e-04,
        7.6990e-05,
        3.3842e-03,
        9.2931e-04,
        2.4958e-04,
        1.0121e-03,
        2.7993e-04,
        6.1043e-05,
        7.7401e-04,
        3.2480e-04,
        9.3145e-06,
        1.4245e-03,
        4.3352e-04,
        1.7392e-04,
        7.3166e-04,
        2.4210e-04,
        1.4804e-04,
    ]
    expected_error_data = [
        2.2149e01,
        7.4509e01,
        1.0000e02,
        2.4621e01,
        4.6383e01,
        3.3621e01,
        2.1616e01,
        7.5885e01,
        1.0000e02,
        2.0067e01,
        3.3654e01,
        6.1265e01,
        1.8407e01,
        1.6239e01,
        5.2119e01,
        1.5791e01,
        3.8452e01,
        1.0000e02,
        7.6577e00,
        3.5290e01,
        1.0000e02,
        8.3702e00,
        5.3283e01,
        6.2602e01,
        1.1655e01,
        6.2289e01,
        6.7541e01,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["multi_p"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["multi_p"].part_data_tag[v_e]
        assert_equal(usrbin_object.tally["multi_p"].part_data_tag.name, "part_data_7")
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["multi_p"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["multi_p"].error_data_tag[v_e]
        assert_equal(usrbin_object.tally["multi_p"].error_data_tag.name, "error_data_7")
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #2:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 0.0, 3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally["multi_n"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["multi_n"].name, "multi_n")
    assert_equal(usrbin_object.tally["multi_n"].particle, "8")
    assert_equal(usrbin_object.tally["multi_n"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["multi_n"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["multi_n"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.0984e-02,
        4.1051e-03,
        1.0636e-03,
        2.1837e-02,
        5.5610e-03,
        1.9119e-03,
        1.0971e-02,
        3.3943e-03,
        1.2456e-03,
        1.6615e-02,
        2.9501e-03,
        7.4597e-04,
        1.0395e-01,
        6.1186e-03,
        1.4997e-03,
        1.7421e-02,
        3.0824e-03,
        7.3878e-04,
        1.8097e-02,
        5.2532e-03,
        2.1572e-03,
        1.0465e-01,
        6.2611e-03,
        1.8829e-03,
        1.7323e-02,
        5.5092e-03,
        2.1418e-03,
    ]
    expected_error_data = [
        5.0179e00,
        1.6521e01,
        1.3973e01,
        4.2025e00,
        8.1766e00,
        1.1465e01,
        7.2005e00,
        1.0479e01,
        1.5640e01,
        5.5994e00,
        1.3275e01,
        2.7617e01,
        7.3788e-01,
        6.7200e00,
        1.9092e01,
        7.3670e00,
        1.3018e01,
        2.8866e01,
        5.7221e00,
        1.5916e01,
        2.6001e01,
        8.3490e-01,
        1.6715e01,
        1.2759e01,
        5.0763e00,
        1.1420e01,
        1.0040e01,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["multi_n"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["multi_n"].part_data_tag[v_e]
        assert usrbin_object.tally["multi_n"].part_data_tag.name == "part_data_8"
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["multi_n"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["multi_n"].error_data_tag[v_e]
        assert usrbin_object.tally["multi_n"].error_data_tag.name == "error_data_8"
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
    assert_equal(usrbin_object.tally["degen1"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["degen1"].name, "degen1")
    assert_equal(usrbin_object.tally["degen1"].particle, "8")
    assert_equal(usrbin_object.tally["degen1"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["degen1"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["degen1"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        3.5279e-02,
        4.7334e-03,
        1.4458e-03,
        3.6242e-02,
        4.6521e-03,
        1.5292e-03,
    ]
    expected_error_data = [
        1.2016e00,
        6.4313e00,
        7.7312e00,
        2.0235e00,
        9.4199e00,
        8.0514e00,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["degen1"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen1"].part_data_tag[v_e]
        assert usrbin_object.tally["degen1"].part_data_tag.name == "part_data_8"
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["degen1"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen1"].error_data_tag[v_e]
        assert usrbin_object.tally["degen1"].error_data_tag.name == "error_data_8"
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #2:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 1.5, 6.0]
    expected_ybounds = [-3.0, 3.0]
    expected_zbounds = [-3.0, -2.0, -1.0, 0.0]
    assert_equal(usrbin_object.tally["degen2"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["degen2"].name, "degen2")
    assert_equal(usrbin_object.tally["degen2"].particle, "8")
    assert_equal(usrbin_object.tally["degen2"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["degen2"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["degen2"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        1.1543e-02,
        2.0295e-03,
        3.2603e-02,
        1.4229e-03,
        3.3492e-02,
        2.7923e-03,
    ]
    expected_error_data = [
        2.7321e00,
        5.2342e00,
        7.4679e-01,
        4.2862e00,
        1.3090e00,
        1.4151e01,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["degen2"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen2"].part_data_tag[v_e]
        assert usrbin_object.tally["degen2"].part_data_tag.name == "part_data_8"
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["degen2"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen2"].error_data_tag[v_e]
        assert usrbin_object.tally["degen2"].error_data_tag.name == "error_data_8"
        expected = expected_error_data[i]
        assert_equal(read, expected)

    # Tally #3:
    # Test UsrbinTally attributes
    expected_xbounds = [-3.0, 6.0]
    expected_ybounds = [-3.0, -1.0, 1.0, 3.0]
    expected_zbounds = [-3.0, -1.5, 0.0]
    assert_equal(usrbin_object.tally["degen3"].coord_sys, "Cartesian")
    assert_equal(usrbin_object.tally["degen3"].name, "degen3")
    assert_equal(usrbin_object.tally["degen3"].particle, "8")
    assert_equal(usrbin_object.tally["degen3"].x_bounds, expected_xbounds)
    assert_equal(usrbin_object.tally["degen3"].y_bounds, expected_ybounds)
    assert_equal(usrbin_object.tally["degen3"].z_bounds, expected_zbounds)

    # Test error and part data values match
    expected_part_data = [
        5.8037e-03,
        1.3260e-02,
        5.6046e-03,
        7.9677e-03,
        4.3111e-02,
        8.1349e-03,
    ]
    expected_error_data = [
        6.1913e00,
        2.3684e00,
        4.6124e00,
        3.2523e00,
        1.3714e00,
        4.3161e00,
    ]

    for i, v_e in enumerate(
        usrbin_object.tally["degen3"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen3"].part_data_tag[v_e]
        assert usrbin_object.tally["degen3"].part_data_tag.name == "part_data_8"
        expected = expected_part_data[i]
        assert_equal(read, expected)

    for i, v_e in enumerate(
        usrbin_object.tally["degen3"].structured_iterate_hex("zyx")
    ):
        read = usrbin_object.tally["degen3"].error_data_tag[v_e]
        assert usrbin_object.tally["degen3"].error_data_tag.name == "error_data_8"
        expected = expected_error_data[i]
        assert_equal(read, expected)


# test file writing to catch upstream changes in mesh


def test_mesh_write():
    if not HAVE_PYMOAB:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    usrbin_file = os.path.join(thisdir, "fluka_usrbin_single.lis")

    usrbin_object = fluka.Usrbin(usrbin_file)
    data = usrbin_object.tally["single_n"]
    data.write_hdf5("test_fluka_data.h5m")
