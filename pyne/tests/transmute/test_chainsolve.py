"""chainsolve transmutation tests."""

import os
import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_is

from numpy.testing import dec, assert_array_equal

import numpy  as np
import tables as tb
from scipy import linalg

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
    orig = np.array([[-0.5,0.,0.],
                     [0.25,-0.3,0.],
                     [0.,0.123,-1.2]])
    exp = np.array([[-0.5,0.,0.,0.],
                    [0.25,-0.3,0.,0.],
                    [0.,0.123,-1.2,0.],
                    [0.,0.,0.1848,-1.337]])
    obs = tm._grow_matrix(orig, prod, dest)
    assert_array_equal(exp, obs)

def test_grow_matrix2():
    prod = 0.1848
    dest = 1.337
    orig = np.array([[-1.]])
    exp = np.array([[-1.,0.],
                    [0.1848,-1.337]])
    obs = tm._grow_matrix(orig, prod, dest)
    assert_array_equal(exp, obs)

"""\

def test_tree_log():
    "Tests corret implementation of the _tree_log() function"
    filename = 'testTreeFile'
    d0 = 0
    d1 = 1
    d2 = 2
    d11 = 1
    d20 = 0
    nuc0 = nn.zzaaam('O16')
    nuc1 = nn.zzaaam('O17')
    nuc2 = nn.zzaaam('O18')
    nuc11 = nn.zzaaam('He4')
    nuc20 = nn.zzaaam('C12')
    N0 = 123.456
    N1 = 12.3456
    N2 = 1.23456
    N11 = 1111.
    N20 = 12.
    temp = ['--> O16 123.456\n']
    temp.append('   |--> O17 12.3456\n')
    temp.append('   |   |--> O18 1.23456\n')
    temp.append('   |--> HE4 1111.0\n')
    temp.append('--> C12 12.0\n')
    with open(filename, 'w') as tree:
        tm._tree_log(d0, nuc0, N0, tree)
        tm._tree_log(d1, nuc1, N1, tree)
        tm._tree_log(d2, nuc2, N2, tree)
        tm._tree_log(d11, nuc11, N11, tree)
        tm._tree_log(d20, nuc20, N20, tree)
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        assert_equal(lines[i], temp[i])
    os.remove(filename)

def test_zero_flux():
    "Tests correct implementation of a transmutation with zero flux on an isotope with a zero decay-constant."
    nuc = nn.zzaaam('FE56')
    t_sim = 100.
    phi = None
    tree = None
    tol = 1e-7
    out = tm.transmute_core(nuc, t_sim, phi, tree, tol)
    assert_equal(out[nuc], 1)

def test_root_decrease():
    "Tests that the root isotope is not being skipped"
    nuc = nn.zzaaam('FE56')
    t_sim = 100.
    phi = np.zeros((175,1))
    for i in np.arange(phi.shape[0]):
        phi[i] = 1e12
    tree = None
    tol = 1e-7
    out = tm.transmute_core(nuc, t_sim, phi, tree, tol)
    assert_true(out[nuc] < 1)

def test_trans_v_transCore():
    "Tests that transmute_core and transmute agree"
    nuc = nn.zzaaam('FE56')
    t_sim = 100.
    phi = np.zeros((175,1))
    for i in np.arange(phi.shape[0]):
        phi[i] = 1.0E+12
    out_core = tm.transmute_core(nuc, t_sim, phi)
    inp = {nuc: 1.0}
    out = tm.transmute(inp, t_sim, phi)
    for key in out.keys():
        assert_true(key in out_core.keys())
        assert_equal(out[key], out_core[key])

def test_trans_v_transSpat():
    "Tests that transmute and transmute_spatial agree"
    nuc = nn.zzaaam('FE56')
    t_sim = 100.
    phi = np.zeros((175,1))
    for i in np.arange(phi.shape[0]):
        phi[i] = 1.0E+12
    inp = {nuc: 1.0}
    inp2 = {nuc: 2.0}
    space = {1: (phi, inp), 2: (phi, inp2)}
    out = tm.transmute(inp, t_sim, phi)
    space_out = tm.transmute_spatial(space, t_sim)
    for key in out.keys():
        assert_true(key in space_out[1].keys())
        assert_true(key in space_out[2].keys())
        assert_equal(out[key], space_out[1][key])
        assert_equal(2*out[key], space_out[2][key])

def test_tm171_decay():
    "Tests if decay is properly implemented"
    nuc = nn.zzaaam('TM171')
    t_sim = 1.2119E+8 # Run for 3.843 years (approx 2 half lives)
    out = tm.transmute_core(nuc, t_sim, None)
    # With no flux, the result should be pure decay
    assert_true(nuc in out.keys())
    tm_res = out[nuc]
    lamb = data.decay_const(nuc)
    analytical = np.exp(-1*lamb*t_sim)
    assert_equal(tm_res, analytical)
"""

#
# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
