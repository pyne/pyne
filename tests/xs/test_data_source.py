from __future__ import unicode_literals, print_function
import os
import sys
import warnings
from io import StringIO

try:
    import urllib.request as urllib
except ImportError:
    import urllib

import numpy as np
import tables as tb

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_almost_equal,
    assert_true,
    assert_is,
)
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne.xs import data_source
from pyne.pyne_config import pyne_conf

sys.path.insert(0, os.path.dirname(__file__))
from utils import download_file

del sys.path[0]

download_file(
    "https://www-nds.iaea.org/wolfram/w180/beta3/W180.ace",
    "W180.ace",
    "5349513f196ad594d172bc6ea61dc382",
)

nuc_data = pyne_conf.NUC_DATA_PATH

# These tests require nuc_data
if not os.path.isfile(nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")

simpleds = data_source.SimpleDataSource()
cinderds = data_source.CinderDataSource()
eafds = data_source.EAFDataSource()


def test_cinder_E_g():
    if not cinderds.exists:
        return
    with tb.open_file(nuc_data, "r") as f:
        E_g = np.array(f.root.neutron.cinder_xs.E_g)
    assert_array_equal(E_g, cinderds.src_group_struct)


def test_cinder_sigma_f():
    if not cinderds.exists:
        return
    with tb.open_file(nuc_data, "r") as f:
        sigma_f_n_U235 = np.array(f.root.neutron.cinder_xs.fission[28]["xs"])
    obs = cinderds.reaction("U235", "fission")
    assert_array_equal(sigma_f_n_U235, obs)
    assert_equal(id(obs), id(cinderds.reaction("U235", "fission")))
    assert_equal(id(obs), id(cinderds.reaction(922350, "fission")))


def test_cinder_sigma_a():
    if not cinderds.exists:
        return
    with tb.open_file(nuc_data, "r") as f:
        sigma_a_n_H1 = np.array(f.root.neutron.cinder_xs.absorption[0]["xs"])
    obs = cinderds.reaction(10010, "absorption")
    assert_array_equal(sigma_a_n_H1, obs)
    assert_equal(id(obs), id(cinderds.reaction(10010, "absorption")))
    assert_equal(id(obs), id(cinderds.reaction("H1", "absorption")))


def test_cinder_sigma_f_n1():
    if not cinderds.exists:
        return
    observed = cinderds.reaction(922350, "fission")
    expected = np.array(
        [
            1.74780000e03,
            1.09570000e03,
            8.54720000e02,
            8.21910000e02,
            5.96110000e02,
            6.55820000e02,
            4.85430000e02,
            5.24960000e02,
            4.01070000e02,
            3.84060000e02,
            8.32680000e02,
            3.68510000e02,
            2.66930000e02,
            2.27710000e02,
            1.83750000e02,
            1.67020000e02,
            7.96280000e01,
            6.53830000e01,
            2.87850000e01,
            1.43510000e01,
            1.87710000e01,
            1.92710000e01,
            7.72680000e01,
            4.90740000e01,
            5.32240000e01,
            4.62680000e01,
            2.47770000e01,
            2.08130000e01,
            2.07720000e01,
            1.36800000e01,
            1.30990000e01,
            8.23490000e00,
            6.81700000e00,
            8.26300000e00,
            5.23320000e00,
            4.96880000e00,
            4.32240000e00,
            3.26220000e00,
            2.71850000e00,
            2.31530000e00,
            2.16830000e00,
            1.98670000e00,
            1.80300000e00,
            1.61870000e00,
            1.46980000e00,
            1.32110000e00,
            1.23810000e00,
            1.18940000e00,
            1.15190000e00,
            1.13810000e00,
            1.18470000e00,
            1.22020000e00,
            1.25640000e00,
            1.29290000e00,
            1.26850000e00,
            1.19820000e00,
            1.12000000e00,
            1.06560000e00,
            1.53220000e00,
            2.06170000e00,
            2.10070000e00,
            1.96770000e00,
            1.96770000e00,
        ][::-1]
    )
    assert_array_equal(observed, expected)


def test_cinder_sigma_f_n2():
    if not cinderds.exists:
        return
    observed = cinderds.reaction(10010, "fission")
    expected = None
    assert_array_equal(observed, expected)


def test_get_sigma_a_n1():
    # Test example with one entry
    if not cinderds.exists:
        return
    observed = cinderds.reaction(10010, "absorption")
    expected = np.array(
        [
            9.96360000e-01,
            6.07160000e-01,
            4.72250000e-01,
            3.99360000e-01,
            3.52130000e-01,
            3.18550000e-01,
            2.93030000e-01,
            2.69290000e-01,
            2.46390000e-01,
            2.27390000e-01,
            2.11380000e-01,
            1.95050000e-01,
            1.76480000e-01,
            1.50820000e-01,
            1.20850000e-01,
            9.42780000e-02,
            7.26880000e-02,
            5.65930000e-02,
            4.38660000e-02,
            3.42250000e-02,
            2.68640000e-02,
            2.08190000e-02,
            1.61410000e-02,
            1.27460000e-02,
            9.89720000e-03,
            7.70440000e-03,
            5.88870000e-03,
            4.59450000e-03,
            3.57770000e-03,
            2.78980000e-03,
            2.18600000e-03,
            1.74630000e-03,
            1.32340000e-03,
            1.13220000e-03,
            1.04010000e-03,
            9.54700000e-04,
            7.95210000e-04,
            6.03360000e-04,
            4.53900000e-04,
            3.60530000e-04,
            2.99660000e-04,
            2.36340000e-04,
            1.71240000e-04,
            1.19680000e-04,
            8.19370000e-05,
            5.64050000e-05,
            4.48700000e-05,
            3.98850000e-05,
            3.71800000e-05,
            3.54070000e-05,
            3.45980000e-05,
            3.44370000e-05,
            3.43310000e-05,
            3.43120000e-05,
            3.49470000e-05,
            3.58100000e-05,
            3.62750000e-05,
            3.60950000e-05,
            3.48160000e-05,
            2.97130000e-05,
            2.89150000e-05,
            2.67690000e-05,
            2.60400000e-05,
        ][::-1]
    )
    assert_array_equal(observed, expected)


def test_get_sigma_a_n2():
    # Test example with multiple entries but not that have reaction_type = 'c'
    if not cinderds.exists:
        return
    observed = cinderds.reaction(10020, "absorption")
    expected = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0011439,
            0.019188,
            0.047752,
            0.087012,
            0.18032,
            0.18847,
            0.20148,
            0.2053,
        ][::-1]
    )
    expected += np.array(
        [
            1.52240000e-03,
            9.30720000e-04,
            7.23200000e-04,
            6.13370000e-04,
            5.40020000e-04,
            4.89210000e-04,
            4.49150000e-04,
            4.13780000e-04,
            3.78790000e-04,
            3.48940000e-04,
            3.24070000e-04,
            2.99220000e-04,
            2.70960000e-04,
            2.31220000e-04,
            1.85480000e-04,
            1.44650000e-04,
            1.11350000e-04,
            8.69150000e-05,
            6.73760000e-05,
            5.26030000e-05,
            4.15620000e-05,
            3.20750000e-05,
            2.48770000e-05,
            1.95970000e-05,
            1.52790000e-05,
            1.17130000e-05,
            9.09390000e-06,
            7.05040000e-06,
            5.53650000e-06,
            4.39520000e-06,
            3.43240000e-06,
            2.68440000e-06,
            2.00390000e-06,
            1.70200000e-06,
            1.58790000e-06,
            1.46300000e-06,
            1.28160000e-06,
            1.09490000e-06,
            1.01380000e-06,
            1.04050000e-06,
            1.07440000e-06,
            1.17800000e-06,
            1.41250000e-06,
            1.77080000e-06,
            2.30510000e-06,
            2.96130000e-06,
            3.51820000e-06,
            3.92870000e-06,
            4.43470000e-06,
            4.95710000e-06,
            5.57910000e-06,
            6.32150000e-06,
            7.01920000e-06,
            7.69640000e-06,
            8.43000000e-06,
            9.14600000e-06,
            9.83610000e-06,
            1.03010000e-05,
            1.06530000e-05,
            9.44650000e-06,
            9.18780000e-06,
            8.38510000e-06,
            8.10500000e-06,
        ][::-1]
    )
    assert_array_equal(observed, expected)


def test_get_sigma_a_n3():
    # Test example with multiple entries including one that has reaction_type = 'c'
    if not cinderds.exists:
        return
    observed = cinderds.reaction(20030, "absorption")
    expected = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0043036,
            0.060056,
            0.11588,
            0.15215,
            0.15165,
            0.12801,
            0.112,
        ][::-1]
    )
    expected += np.array(
        [
            1.59870000e04,
            9.74200000e03,
            7.57740000e03,
            6.40770000e03,
            5.65000000e03,
            5.11100000e03,
            4.70140000e03,
            4.32040000e03,
            3.95290000e03,
            3.64790000e03,
            3.39100000e03,
            3.12890000e03,
            2.83090000e03,
            2.41910000e03,
            1.93810000e03,
            1.51180000e03,
            1.16550000e03,
            9.07320000e02,
            7.03200000e02,
            5.48580000e02,
            4.30550000e02,
            3.33630000e02,
            2.58630000e02,
            2.04220000e02,
            1.58550000e02,
            1.23410000e02,
            9.43150000e01,
            7.36280000e01,
            5.73090000e01,
            4.46360000e01,
            3.50230000e01,
            2.79830000e01,
            2.10090000e01,
            1.78020000e01,
            1.61630000e01,
            1.46700000e01,
            1.20900000e01,
            9.13540000e00,
            6.93580000e00,
            5.58270000e00,
            4.72330000e00,
            3.83980000e00,
            2.91150000e00,
            2.19690000e00,
            1.63240000e00,
            1.22460000e00,
            1.05240000e00,
            9.59460000e-01,
            9.20880000e-01,
            9.06540000e-01,
            8.90440000e-01,
            8.79710000e-01,
            8.66060000e-01,
            8.25060000e-01,
            7.46760000e-01,
            6.08280000e-01,
            4.47760000e-01,
            3.38630000e-01,
            2.50370000e-01,
            1.21070000e-01,
            1.11410000e-01,
            8.97170000e-02,
            8.20000000e-02,
        ][::-1]
    )
    expected += np.array(
        [
            1.40840000e-03,
            1.10660000e-03,
            8.04470000e-04,
            5.02330000e-04,
            2.00200000e-04,
            3.15200000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.10000000e-05,
            3.09990000e-05,
            3.09990000e-05,
            3.09980000e-05,
            3.09970000e-05,
            3.09950000e-05,
            3.09920000e-05,
            3.09860000e-05,
            3.09770000e-05,
            3.09630000e-05,
            3.09390000e-05,
            3.08990000e-05,
            3.08620000e-05,
            3.08380000e-05,
            3.08080000e-05,
            3.07250000e-05,
            3.05460000e-05,
            3.02520000e-05,
            2.99180000e-05,
            2.95930000e-05,
            2.89430000e-05,
            2.76470000e-05,
            2.54710000e-05,
            2.18830000e-05,
            1.59690000e-05,
            9.60350000e-06,
            3.53490000e-06,
            5.38530000e-07,
            1.74800000e-06,
            3.50630000e-06,
            5.50360000e-06,
            7.86570000e-06,
            1.11590000e-05,
            1.53870000e-05,
            2.08170000e-05,
            2.86970000e-05,
            3.76470000e-05,
            5.65300000e-05,
            8.97390000e-05,
            1.15640000e-04,
            1.34700000e-04,
            1.46310000e-04,
        ][::-1]
    )
    assert_array_equal(observed, expected)


def test_get_sigma_a_n4():
    if not cinderds.exists:
        return
    # Test that a zeros array is returned for an entry that is not in the table
    observed = cinderds.reaction(10420, "absorption")
    expected = None
    assert_array_equal(observed, expected)


def test_simple_E_g():
    if not simpleds.exists:
        return
    assert_array_equal(np.array([14.0, 1.0, 2.53e-8, 0.0]), simpleds.src_group_struct)


def test_simple_sigma_f():
    if not simpleds.exists:
        return
    sigma_f_n_U235 = np.array([2.056, 1.235, 584.4])
    obs = simpleds.reaction("U235", "fission")
    assert_array_equal(sigma_f_n_U235, obs)
    assert_equal(id(obs), id(simpleds.reaction("U235", "fission")))
    assert_equal(id(obs), id(simpleds.reaction(922350, "fission")))


def test_simple_sigma_a():
    if not simpleds.exists:
        return
    sigma_a_n_H1 = np.array([2.983e-5, 3.927e-5, 0.332])
    obs = simpleds.reaction(10010, "absorption")
    assert_array_almost_equal(sigma_a_n_H1, obs)
    assert_equal(id(obs), id(simpleds.reaction(10010, "absorption")))
    assert_equal(id(obs), id(simpleds.reaction("H1", "absorption")))


def test_simple_not_in_table():
    if not simpleds.exists:
        return
    exp = None
    obs = simpleds.reaction(10030, "z_4n")
    assert_equal(obs, exp)


def test_simple_not_a_rx():
    if not simpleds.exists:
        return
    exp = None
    obs = simpleds.reaction(10010, "42")
    assert_equal(obs, exp)


def test_simple_discretize_no_weights1():
    if not simpleds.exists or not cinderds.exists:
        return
    dst_g = cinderds.src_group_struct
    simpleds.dst_group_struct = dst_g
    obs = simpleds.discretize("U235", "fission")
    mask = obs[:-1] <= obs[1:]
    assert_true((obs[mask][:-1] <= obs[mask][1:]).all())
    assert_true((obs[~mask][:-1] >= obs[~mask][1:]).all())
    simpleds.dst_group_struct = None


def test_simple_discretize_no_weights2():
    if not simpleds.exists or not cinderds.exists:
        return
    dst_g = cinderds.src_group_struct
    simpleds.dst_group_struct = dst_g
    obs = simpleds.discretize(10010, "a")
    assert_true((obs[:-1] <= obs[1:]).all())
    simpleds.dst_group_struct = None


def test_simple_discretize_weights1():
    if not simpleds.exists or not cinderds.exists:
        return
    dst_g = cinderds.src_group_struct
    phi_g = np.ones(cinderds.src_ngroups, dtype=float)
    phi_g[:25] = 0.0
    simpleds.dst_group_struct = dst_g
    obs = simpleds.discretize("U235", "fission", dst_phi_g=phi_g)
    assert_true((obs[:-1] <= obs[1:]).all())
    simpleds.dst_group_struct = None


def test_shield_weights1():
    mat = {922350000: 0.5, 922380000: 0.5}
    simpleds.shield_weights(mat, 300)
    assert_array_equal(
        simpleds.slf_shld_wgts[922350000], simpleds.slf_shld_wgts[922380000]
    )


def test_eaf_E_g():
    if not eafds.exists:
        return
    with tb.open_file(nuc_data, "r") as f:
        E_g = np.array(f.root.neutron.eaf_xs.E_g)
    assert_array_equal(E_g, eafds.src_group_struct)


def test_eaf_valid_mtnum_RX():
    if not eafds.exists:
        return
    observed = eafds.reaction(250550, "22")
    expected = [
        4.48860e-02,
        2.89454e-02,
        2.36589e-02,
        1.61601e-02,
        6.04350e-03,
        3.65975e-03,
        2.11363e-03,
        1.14410e-03,
        7.26727e-04,
        4.03486e-04,
        6.35067e-05,
        3.15836e-05,
        4.62556e-06,
        1.31764e-06,
        1.09397e-06,
        8.81033e-07,
        6.73990e-07,
        4.81856e-07,
        2.99065e-07,
        1.25156e-07,
        5.38030e-09,
    ]
    while len(expected) < 175:
        expected.append(0.0)
    expected = np.array(expected)
    assert_array_equal(observed, expected)


def test_eaf_invalid_mtnum_RX():
    if not eafds.exists:
        return
    observed = eafds.reaction(250550, 24)
    expected = None
    assert_array_equal(observed, expected)


def test_eaf_valid_str_RX():
    if not eafds.exists:
        return
    observed = eafds.reaction(250550, "na")
    expected = [
        4.48860e-02,
        2.89454e-02,
        2.36589e-02,
        1.61601e-02,
        6.04350e-03,
        3.65975e-03,
        2.11363e-03,
        1.14410e-03,
        7.26727e-04,
        4.03486e-04,
        6.35067e-05,
        3.15836e-05,
        4.62556e-06,
        1.31764e-06,
        1.09397e-06,
        8.81033e-07,
        6.73990e-07,
        4.81856e-07,
        2.99065e-07,
        1.25156e-07,
        5.38030e-09,
    ]
    while len(expected) < 175:
        expected.append(0.0)
    expected = np.array(expected)
    assert_array_equal(observed, expected)


def test_eaf_invalid_str_RX():
    if not eafds.exists:
        return
    observed = eafds.reaction(10010, "fission")
    expected = None
    assert_array_equal(observed, expected)


def test_eaf_multiple_xs():
    # Currently no case where multiple rows from nuc_data should be combined...
    pass


sample_xs_openmc = StringIO(
    """<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="W-180.21c" awr="178.401" location="1" name="74180.21c" path="W180.ace" temperature="2.585e-08" zaid="74180"/>
  <ace_table alias="C-12.00c" awr="11.896900"location="1" name="6000.00c"  path="C012-n.ace" temperature="2.5263e-08" zaid="6000"/>
</cross_sections>
"""
)


def test_openmc():
    sample_xs_openmc.seek(0)
    ods = data_source.OpenMCDataSource(
        cross_sections=sample_xs_openmc, src_group_struct=np.logspace(1, -9, 11)
    )
    obs = ods.reaction("W180", 2)
    assert_equal(10, len(obs))

    obs = ods.reaction("W180", "total")
    assert_equal(10, len(obs))

    obs = ods.reaction("U-235", 42)
    assert_is(None, obs)

    # threshold reactoin
    obs = ods.reaction("W180", "z_3n")
    assert_equal(10, len(obs))
    assert_true(np.all(obs >= 0.0))


def test_openmc_bkg_none():
    C_12 = 60120000
    W_180 = 741800000
    sample_xs_openmc.seek(0)
    ods = data_source.OpenMCDataSource(
        cross_sections=sample_xs_openmc, src_group_struct=np.logspace(1, -9, 11)
    )
    expected = ods.reaction(C_12, "total", 300)
    atom_dens = {C_12: 1.0e22, W_180: 1.0e22}
    ods.atom_dens = atom_dens

    observed = ods.bkg_xs(W_180, 300)
    assert_array_almost_equal(0.0, observed)


def test_openmc_bkg():
    C = 60000000
    W_180 = 741800000
    sample_xs_openmc.seek(0)
    ods = data_source.OpenMCDataSource(
        cross_sections=sample_xs_openmc, src_group_struct=np.logspace(1, -9, 11)
    )
    expected = ods.reaction(C, "total", 300)
    atom_dens = {C: 1.0e22, W_180: 1.0e22}
    ods.atom_dens = atom_dens

    observed = ods.bkg_xs(W_180, 300)
    if expected is not None and observed is not None:
        assert_array_almost_equal(expected, observed)


def test_openmc_self_shielding1():
    C = 60000000
    W_180 = 741800000
    sample_xs_openmc.seek(0)
    ods = data_source.OpenMCDataSource(
        cross_sections=sample_xs_openmc, src_group_struct=np.logspace(1, -9, 11)
    )
    non_ss = ods.reaction(W_180, "gamma", 300)
    atom_dens = {C: 1.0e22, W_180: 1.0e22}
    ods.atom_dens = atom_dens

    observed = ods.reaction(W_180, "gamma", 300)
    assert_true(np.all(observed <= non_ss))


def test_open_self_shielding2():
    C = 60000000
    W_180 = 741800000
    sample_xs_openmc.seek(0)
    ods = data_source.OpenMCDataSource(
        cross_sections=sample_xs_openmc, src_group_struct=np.logspace(1, -9, 11)
    )
    non_ss = ods.reaction(W_180, "gamma", 300)
    atom_dens = {C: 1.0e22, W_180: 1.0e2}
    ods.atom_dens = atom_dens

    observed = ods.reaction(W_180, "gamma", 300)
    assert_array_almost_equal(observed, non_ss)
