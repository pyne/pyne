"""chainsolve transmutation tests."""
import os
import nose
import warnings

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_almost_equal,
    assert_true,
    assert_false,
    assert_is,
    with_setup,
    assert_less,
)

from numpy.testing import dec, assert_array_equal

import numpy as np
import tables as tb
from scipy import linalg

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne import nuc_data
from pyne import nucname as nn
from pyne import data
from pyne.material import Material
from pyne.transmute.chainsolve import Transmuter

tm = None


def setup():
    global tm
    tm = Transmuter()


def teardown():
    global tm
    del tm


def test_check_phi():
    """Tests the _check_phi function"""
    numeaf = 175

    def set_phi(f):
        tm.phi = f

    # First check that None is properly converted
    tm._phi = None
    assert_is(tm.phi, None)
    tm.phi = np.ones(numeaf)
    assert_array_equal(tm.phi, np.ones(numeaf))
    # Check that incorrect number of entries raises an exception
    assert_raises(ValueError, set_phi, np.ones((50, 1)))
    # Check that a negative entry raises an exception
    x = np.ones(numeaf)
    x[123] = -1
    assert_raises(ValueError, set_phi, x)


def test_grow_matrix1():
    "Tests correct implementation of the _grow_matrix function"
    prod = 0.1848
    dest = 1.337
    orig = np.array([[-0.5, 0.0, 0.0], [0.25, -0.3, 0.0], [0.0, 0.123, -1.2]])
    exp = np.array(
        [
            [-0.5, 0.0, 0.0, 0.0],
            [0.25, -0.3, 0.0, 0.0],
            [0.0, 0.123, -1.2, 0.0],
            [0.0, 0.0, 0.1848, -1.337],
        ]
    )
    obs = tm._grow_matrix(orig, prod, dest)
    assert_array_equal(exp, obs)


def test_grow_matrix2():
    prod = 0.1848
    dest = 1.337
    orig = np.array([[-1.0]])
    exp = np.array([[-1.0, 0.0], [0.1848, -1.337]])
    obs = tm._grow_matrix(orig, prod, dest)
    assert_array_equal(exp, obs)


@with_setup(None, lambda: os.remove("log.txt") if os.path.exists("log.txt") else None)
def test_tree_log():
    "Tests corret implementation of the _log_tree() function"
    filename = "log.txt"
    tm.log = open(filename, "w")
    d0 = 0
    d1 = 1
    d2 = 2
    d11 = 1
    d20 = 0
    nuc0 = nn.id("O16")
    nuc1 = nn.id("O17")
    nuc2 = nn.id("O18")
    nuc11 = nn.id("He4")
    nuc20 = nn.id("C12")
    N0 = 123.456
    N1 = 12.3456
    N2 = 1.23456
    N11 = 1111.0
    N20 = 12.0
    exp = (
        "--> O16 123.456\n"
        "   |--> O17 12.3456\n"
        "   |   |--> O18 1.23456\n"
        "   |--> He4 1111.0\n"
        "--> C12 12.0\n"
    )
    with open(filename, "w") as tree:
        tm._log_tree(d0, nuc0, N0)
        tm._log_tree(d1, nuc1, N1)
        tm._log_tree(d2, nuc2, N2)
        tm._log_tree(d11, nuc11, N11)
        tm._log_tree(d20, nuc20, N20)
    tm.log.close()
    tm.log = None
    with open(filename, "r") as f:
        obs = f.read()
    # print repr(exp)
    # print repr(obs)
    # print obs == exp
    assert_equal(exp, obs)


def test_zero_flux():
    """Tests correct implementation of a transmutation with zero flux on
    an isotope with a zero decay-constant."""
    inp = Material({"FE56": 1.0}, mass=1.0)
    obs = tm.transmute(inp, t=100.0, tol=1e-7)
    assert_almost_equal(obs["FE56"], 1.0)


def test_root_decrease():
    "Tests that the root isotope is not being skipped"
    phi = 1e12 * np.ones(175)
    inp = Material({"FE56": 1.0}, mass=1.0)
    obs = tm.transmute(inp, t=100.0, phi=phi, tol=1e-7)
    assert_less(obs["FE56"], 1.0)


def test_tm171_decay():
    "Tests if decay is properly implemented"
    t_sim = 1.2119e8  # Run for 3.843 years (approx 2 half lives)
    lamb = data.decay_const("TM171")
    exp = np.exp(-1 * lamb * t_sim)
    inp = Material({"TM171": 1.0}, mass=1.0)
    obs = tm.transmute(inp, t=t_sim, phi=0.0, tol=1e-7)
    assert_almost_equal(exp, obs["TM171"], 12)


#
# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
