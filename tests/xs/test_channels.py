import os

import numpy as np
import tables as tb

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_almost_equal,
    assert_true,
    assert_raises,
)
from numpy.testing import assert_array_equal, assert_array_almost_equal

import pyne
import pyne.data
import pyne.xs.models
from pyne.xs.cache import xs_cache
from pyne.xs.channels import (
    sigma_f,
    sigma_s_gh,
    sigma_s,
    sigma_a_reaction,
    metastable_ratio,
    sigma_a,
    chi,
    sigma_t,
    _atom_mass_channel,
)
from pyne.pyne_config import pyne_conf
from pyne.material import Material

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


def setup():
    np.seterr(all="ignore")
    xs_cache.clear()


#
# Test helper functions
#
def test_atom_weight_channel1():
    xs_cache["E_g"] = np.array([3, 2, 1.0])
    chanfunc = lambda nuc: np.array([1.0, nuc], float)

    # Test dict
    nucspec = {1: 1, 10: 2}
    obs = _atom_mass_channel(chanfunc, nucspec)
    exp = np.array([1.0, 7.0])
    assert_array_equal(obs, exp)

    # Test list of tuples
    nucspec = [(1, 1), (10, 2)]
    obs = _atom_mass_channel(chanfunc, nucspec)
    exp = np.array([1.0, 7.0])
    assert_array_equal(obs, exp)

    # test material
    h2o = Material({10010000: 0.11191487328888054, 80160000: 0.8880851267111195})
    obs = _atom_mass_channel(chanfunc, h2o)
    exp = np.array([1.0, 33393333.333333336])
    assert_array_almost_equal(obs, exp)


def test_atom_weight_channel2():
    xs_cache["E_g"] = np.logspace(-6, 1, 10)[::-1]
    exp = (sigma_t("H1") * 2.0 + sigma_t("O16")) / 3.0

    # Test dict
    nucspec = {"H1": 2.0, "O16": 1.0}
    obs = _atom_mass_channel(sigma_t, nucspec)
    assert_array_almost_equal(obs, exp)

    # Test list of tuples
    nucspec = [("H1", 2.0), ("O16", 1.0)]
    obs = _atom_mass_channel(sigma_t, nucspec)
    assert_array_equal(obs, exp)

    # test material
    h2o = Material({10010000: 0.11191487328888054, 80160000: 0.8880851267111195})
    obs = _atom_mass_channel(sigma_t, h2o)
    assert_array_almost_equal(obs, exp)


#
# Test channels, make sure they run rather than test their values.
# This is OK since the underlying functions are very well tested.
#


def test_sigma_f():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    sig_f = sigma_f("U238", group_struct=E_g)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f("U238", group_struct=E_g)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f("U235")
    observed = (0.0 <= sig_f).all()
    assert_true(observed)


def test_sigma_s_gh():
    # Tests stub
    b = pyne.data.b("H1")
    aw = pyne.data.atomic_mass("H1")
    E = np.logspace(-6, 1, 10)[::-1]
    E_centers = (E[1:] + E[:-1]) / 2.0
    expected = np.diag(pyne.xs.models.sigma_s(E_centers, b, aw, 600.0))
    observed = sigma_s_gh("H1", 600.0, group_struct=E)
    assert_array_equal(expected, observed)


def test_sigma_s():
    E_g = np.logspace(-6, 1, 10)[::-1]
    expected = sigma_s_gh("H1", 600.0, E_g).sum(axis=1)
    observed = sigma_s("H1", 600.0, E_g)
    assert_array_equal(expected, observed)


def test_sigma_a_reaction():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    sig_rx = sigma_a_reaction("U238", "z_2n", group_struct=E_g)
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)

    sig_rx = sigma_a_reaction("U238", "gamma", group_struct=E_g)
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)

    sig_rx = sigma_a_reaction("H1", "gamma")
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)


def test_metastable_ratio():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    ms_rx = metastable_ratio("U238", "z_2n", group_struct=E_g)
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)

    ms_rx = metastable_ratio("U238", "gamma", group_struct=E_g)
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)

    ms_rx = metastable_ratio("H1", "gamma")
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)


def test_sigma_a():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    sig_a = sigma_a("U238", group_struct=E_g)
    observed = (0.0 <= sig_a).all()
    assert_true(observed)

    sig_a = sigma_a("U238", group_struct=E_g)
    observed = (0.0 <= sig_a).all()
    assert_true(observed)

    sig_a = sigma_a("U235")
    observed = (0.0 <= sig_a).all()
    assert_true(observed)


def test_chi():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    c = chi("U238", group_struct=E_g)
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi("U238", group_struct=E_g)
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi("U235")
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi("H1")
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 0.0)


def test_sigma_t():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    sig_t = sigma_t("U238", 600.0, E_g)
    observed = (0.0 <= sig_t).all()
    assert_true(observed)
    expected = sigma_a("U238", 600.0, E_g) + sigma_s("U238", 600.0, E_g)
    assert_array_almost_equal(sig_t, expected)

    sig_t = sigma_t("U238", 600.0, E_g)
    observed = (0.0 <= sig_t).all()
    assert_true(observed)
    expected = sigma_a("U238", 600.0, E_g) + sigma_s("U238", 600.0, E_g)
    assert_array_almost_equal(sig_t, expected)

    sig_t = sigma_t("U235")
    observed = (0.0 <= sig_t).all()
    assert_true(observed)
    expected = sigma_a("U235", 600.0) + sigma_s("U235", 600.0, E_g)
    assert_array_almost_equal(sig_t, expected)
