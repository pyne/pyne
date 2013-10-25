"""Chainsolve transmutation tests."""

import nose
import os

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

from numpy.testing import dec

import numpy  as np
import scipy  as sp
import tables as tb
from scipy import linalg

from pyne.material import Material
from pyne import nuc_data
from pyne import nucname as nn
from pyne import data
from pyne.transmute.chainsolve import Transmuter


def test_check_phi():
    """Tests the _check_phi function"""
    eaf_numEntries = 175
    # First check that None is properly converted
    nothing = None
    phi = tm._check_phi(None)
    assert_equal(phi.ndim, 2)
    assert_equal(phi.shape[0], eaf_numEntries)
    assert_equal(phi.shape[1], 1)
    for entry in phi:
        assert_equal(entry, 0)
    # Check that incorrect shape raises an exception
    phi = np.ones((175))
    assert_raises(ValueError, tm._check_phi, phi)
    # Check that incorrect number of entries raises an exception
    phi = np.ones((50,1))
    assert_raises(ValueError, tm._check_phi, phi)
    # Check that a negative entry raises an exception
    phi = np.ones((175,1))
    phi[123] = -1
    assert_raises(ValueError, tm._check_phi, phi)

def test_get_daughters():
    """Tests correct application of the _get_daughters function"""
    nuc = nn.zzaaam("O16")
    with tb.openFile(nuc_data) as data_table:
        daughtersTest = [row['daughter'] for row in \
            data_table.root.neutron.eaf_xs.eaf_xs.where('nuc_zz == nuc')]
        for i in range(len(daughtersTest)):
            daughtersTest[i] = tm._convert_eaf(daughtersTest[i])
        daughter_dict = tm._get_daughters(nuc, data_table)
    for daugh in daughter_dict.keys():
        assert(daugh in daughtersTest)

def test_convert_eaf():
    """Tests conversion of EAF formatted nuc strings to zzaaam format"""
    nuc1 = nn.zzaaam('O16')
    test1 = tm._convert_eaf('O 16')
    assert_equal(nuc1, test1)
    nuc2 = nn.zzaaam('AU196')
    test2 = tm._convert_eaf('AU196G')
    assert_equal(nuc2, test2)
    nuc3 = nn.zzaaam('AU196M')
    test3 = tm._convert_eaf('AU196M1')
    assert_equal(nuc3, test3)
    nuc4 = nn.zzaaam('AU196') + 2
    test4 = tm._convert_eaf('AU196M2')
    assert_equal(nuc4, test4)
    nuc5 = nn.zzaaam('MG28')
    test5 = tm._convert_eaf('MG28')
    assert_equal(nuc5, test5)

def test_get_decay():
    """Tests correct implementation of the _get_decay function"""
    manual_dict = {}
    nuc = nn.zzaaam('O16')
    children = data.decay_children(nuc)
    for child in children:
        branch = data.branch_ratio(nuc,child)
        manual_dict[child] = branch
    transmute_dict = tm._get_decay(nuc)
    assert_equal(manual_dict, transmute_dict)

def test_grow_matrix():
    """Tests correct implementation of the _grow_matrix function"""
    orig = np.array([[-0.5,0.,0.],
                     [0.25,-0.3,0.],
                     [0.,0.123,-1.2]])
    prod = 0.1848
    dest = 1.337
    manual = np.array([[-0.5,0.,0.,0.],
                       [0.25,-0.3,0.,0.],
                       [0.,0.123,-1.2,0.],
                       [0.,0.,0.1848,-1.337]])
    method = tm._grow_matrix(orig, prod, dest)
    assert_true(np.array_equal(manual, method))
    orig = np.array([[-1.]])
    manual = np.array([[-1.,0.],
                       [0.1848,-1.337]])
    method = tm._grow_matrix(orig, prod, dest)
    assert_true(np.array_equal(manual, method))


def test_tree_log():
    """Tests corret implementation of the _tree_log() function"""
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
    """Tests correct implementation of a transmutation with zero flux on an
    isotope with a zero decay-constant."""
    nuc = nn.zzaaam('FE56')
    t_sim = 100.
    phi = None
    tree = None
    tol = 1e-7
    out = tm.transmute_core(nuc, t_sim, phi, tree, tol)
    assert_equal(out[nuc], 1)

def test_root_decrease():
    """Tests that the root isotope is not being skipped"""
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
    """Tests that transmute_core and transmute agree"""
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
    """Tests that transmute and transmute_spatial agree"""
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
    """Tests if decay is properly implemented"""
    nuc = nn.zzaaam('TM171')
    t_sim = 1.2119E+8 # Run for 3.843 years (approx 2 half lives)
    out = tm.transmute_core(nuc, t_sim, None)
    # With no flux, the result should be pure decay
    assert_true(nuc in out.keys())
    tm_res = out[nuc]
    lamb = data.decay_const(nuc)
    analytical = np.exp(-1*lamb*t_sim)
    assert_equal(tm_res, analytical)

#
# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
