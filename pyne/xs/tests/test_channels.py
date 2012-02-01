import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal, assert_almost_equal, assert_true, \
                       assert_raises
from numpy.testing import assert_array_equal, assert_array_almost_equal

import pyne.data
import pyne.xs.models
from pyne.xs.cache import xs_cache
from pyne.xs.channels import sigma_f, sigma_s_gh, sigma_s, sigma_a_reaction, \
    metastable_ratio, sigma_a, chi
from pyne.pyne_config import pyne_conf



#
# Test channels, make sure they run rather than test their values.
# This is OK since the underlying functions are very well tested.
#

def test_sigma_f():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    sig_f = sigma_f('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f('U235')
    observed = (0.0 <= sig_f).all()
    assert_true(observed)


def test_sigma_s_gh():
    # Tests stub
    b = pyne.data.b('H1')
    aw = pyne.data.nuc_weight('H1')
    E = np.logspace(-6, 1, 10)[::-1]
    E_centers = (E[1:] + E[:-1]) / 2.0
    expected = np.diag(pyne.xs.models.sigma_s(E_centers, b, aw, 600.0))
    observed = sigma_s_gh('H1', 600.0, E_g=E)
    assert_array_equal(expected, observed)


def test_sigma_s():
    E_g = np.logspace(-6, 1, 10)[::-1]
    expected = sigma_s_gh('H1', 600.0, E_g=E_g).sum(axis=1)
    observed = sigma_s('H1', 600.0, E_g=E_g)
    assert_array_equal(expected, observed)


def test_sigma_a_reaction():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    sig_rx = sigma_a_reaction('U238', '2n', E_g, E_n, phi_n)
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)

    sig_rx = sigma_a_reaction('U238', 'gamma', E_g, E_n, phi_n)
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)

    sig_rx = sigma_a_reaction('H1', 'g')
    observed = (0.0 <= sig_rx).all()
    assert_true(observed)


def test_metastable_ratio():
    # Hide warnings from numpy
    np.seterr(divide='ignore')

    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    ms_rx = metastable_ratio('U238', '2n', E_g, E_n, phi_n)
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)

    ms_rx = metastable_ratio('U238', 'gamma', E_g, E_n, phi_n)
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)

    ms_rx = metastable_ratio('H1', 'g')
    observed = (0.0 <= ms_rx).all()
    assert_true(observed)


def test_sigma_a():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    sig_a = sigma_a('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_a).all()
    assert_true(observed)

    sig_a = sigma_a('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_a).all()
    assert_true(observed)

    sig_a = sigma_a('U235')
    observed = (0.0 <= sig_a).all()
    assert_true(observed)


def test_chi():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    c = chi('U238', E_g, E_n, phi_n)
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi('U238', E_g, E_n, phi_n)
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi('U235')
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 1.0)

    c = chi('H1')
    observed = (0.0 <= c).all()
    assert_true(observed)
    assert_almost_equal(c.sum(), 0.0)
