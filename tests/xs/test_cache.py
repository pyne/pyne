import os

import numpy as np
import tables as tb

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_almost_equal,
    assert_true,
    assert_false,
)
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.xs import data_source
from pyne.xs.cache import xs_cache
from pyne.pyne_config import pyne_conf

nuc_data = pyne_conf.NUC_DATA_PATH
cinderds = data_source.CinderDataSource()

# These tests require nuc_data
if not os.path.isfile(nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


def test_xs_cache_sigma_f_n():
    xs_cache.clear()
    if not cinderds.exists:
        return

    with tb.open_file(nuc_data, "r") as f:
        sigma_f_n_U235 = np.array(f.root.neutron.cinder_xs.fission[28]["xs"])
    from_cache = xs_cache[922350, "fiss"]

    assert_not_equal(id(sigma_f_n_U235), id(from_cache))
    assert_equal(id(from_cache), id(xs_cache[922350, "fiss"]))
    assert_array_equal(sigma_f_n_U235, xs_cache[922350, "fiss"])


def test_xs_cache_sigma_a_n():
    xs_cache.clear()
    if not cinderds.exists:
        return

    with tb.open_file(nuc_data, "r") as f:
        sigma_a_n_H1 = np.array(f.root.neutron.cinder_xs.absorption[0]["xs"])
    from_cache = xs_cache[10010, "abs"]

    assert_not_equal(id(sigma_a_n_H1), id(from_cache))
    assert_equal(id(from_cache), id(xs_cache[10010, "abs"]))
    assert_array_equal(sigma_a_n_H1, xs_cache[10010, "abs"])


def test_xs_cache_set_E_g():
    xs_cache.clear()

    # Add an energy stucture
    xs_cache["E_g"] = [10.0, 1.0]
    E_g = xs_cache["E_g"]

    # Assert that the cache is working
    assert_equal(E_g.shape, (2,))
    assert_equal(id(E_g), id(xs_cache["E_g"]))

    # Assert that the cache has been reloaded
    xs_cache["E_g"] = [10.0, 2.0, 1.0]
    assert_not_equal(id(E_g), id(xs_cache["E_g"]))
    assert_equal(len(E_g), 2)
    assert_equal(len(xs_cache["E_g"]), 3)

    # Test auto-clearing
    xs_cache[10010, "fiss"]
    assert_true((10010, "fiss") in xs_cache)
    xs_cache["E_g"] = [10.0, 8.0, 2.0, 1.0]
    assert_true(xs_cache["phi_g"] is None)
    assert_false((10010, "fiss") in xs_cache)


def text_xs_get_reaction():
    xs_cache.clear()
    assert_raises(KeyError, xs_cache[10010, 1089, 300])


def test_xs_cache_get_phi_g():
    xs_cache.clear()
    xs_cache["E_g"] = np.array([1e-8, 5.0, 10.0])
    xs_cache["phi_g"] = [1.0, 1.0]
    phi_g = xs_cache["phi_g"]
    expected = np.array([1.0, 1.0])
    assert_array_equal(phi_g, expected)
