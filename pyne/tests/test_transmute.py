"""Transmute tests"""

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
from pyne import nucname
from pyne import data
from pyne import transmute
from pyne.xs.data_source import EAF_RX


"""Tests _solve_decay_matrix calculating the matrix exponential"""
def test_expm():
    A = np.zeros((3,3))
    t = 1.
    np.fill_diagonal(A,[-1,-2,-3])
    transA = transmute._matrix_exp(A, t)
    eA = linalg.expm2(A)
    tayA = linalg.expm3(A,q=20)
    # Check the equivalence of _solve_decay_matrix
    assert_true(np.array_equal(eA, transA))


"""Ensures the flux vector format is completing"""
"""
def test_phi_size_check():
    phiPass = np.arange(175)
    phiPass = phiPass.reshape((175,1))
    assert_true(np.array_equal(phiPass, transmute._format_phi(phiPass)))
    phi = np.arange(100)
    phi2 = phi.reshape((100,1))
    out = np.append(phi2, np.zeros((75,1)), 0)
    assert_true(np.array_equal(out, transmute._format_phi(phi)))
"""


"""Tests correct application of the _get_daughters function"""
def test_get_daughters():
    nuc = nucname.zzaaam("O16")
    data_table = tb.openFile(nuc_data)
    daughtersTest = [row['daughter'] for row in \
        data_table.root.neutron.eaf_xs.eaf_xs.where('nuc_zz == nuc')]
    data_table.close()
    for i in range(len(daughtersTest)):
        daughtersTest[i] = transmute._convert_eaf(daughtersTest[i])
    daughter_dict = transmute._get_daughters(nuc)
    for daugh in daughter_dict.keys():
        assert(daugh in daughtersTest)


"""Tests conversion of EAF formatted nuc strings to zzaaam format"""
def test_convert_eaf():
    nuc1 = nucname.zzaaam('O16')
    test1 = transmute._convert_eaf('O 16')
    assert_equal(nuc1, test1)
    nuc2 = nucname.zzaaam('AU196')
    test2 = transmute._convert_eaf('AU196G')
    assert_equal(nuc2, test2)
    nuc3 = nucname.zzaaam('AU196M')
    test3 = transmute._convert_eaf('AU196M1')
    assert_equal(nuc3, test3)
    nuc4 = nucname.zzaaam('AU196') + 2
    test4 = transmute._convert_eaf('AU196M2')
    assert_equal(nuc4, test4)

"""Tests correct implementation of the _get_decay function"""
def test_get_decay():
    manual_dict = {}
    nuc = nucname.zzaaam('O16')
    children = data.decay_children(nuc)
    for child in children:
        branch = data.branch_ratio(nuc,child)
        manual_dict[child] = branch
    transmute_dict = transmute._get_decay(nuc)
    assert_equal(manual_dict, transmute_dict)

"""Tests correct implementation of the _grow_matrix function"""
def test_grow_matrix():
    orig = np.array([[-0.5,0.,0.],
                     [0.25,-0.3,0.],
                     [0.,0.123,-1.2]])
    prod = 0.1848
    dest = 1.337
    manual = np.array([[-0.5,0.,0.,0.],
                       [0.25,-0.3,0.,0.],
                       [0.,0.123,-1.2,0.],
                       [0.,0.,0.1848,-1.337]])
    method = transmute._grow_matrix(orig, prod, dest)
    assert_true(np.array_equal(manual, method))
    orig = np.array([[-1.]])
    manual = np.array([[-1.,0.],
                       [0.1848,-1.337]])
    method = transmute._grow_matrix(orig, prod, dest)
    assert_true(np.array_equal(manual, method))

"""Tests correct implementation of the _check_tol function"""
def test_check_tol():
    current = 1848.3
    tol = 1337.5
    N_ini = 1
    assert_true(transmute._check_tol(current, tol))

"""Tests corret implementation of the _tree_log() function"""
def test_tree_log():
    filename = 'testTreeFile'
    d0 = 0
    d1 = 1
    d2 = 2
    d11 = 1
    d20 = 0
    nuc0 = nucname.zzaaam('O16')
    nuc1 = nucname.zzaaam('O17')
    nuc2 = nucname.zzaaam('O18')
    nuc11 = nucname.zzaaam('He4')
    nuc20 = nucname.zzaaam('C12')
    N0 = 123.456
    N1 = 12.3456
    N2 = 1.23456
    N11 = 1111.
    N20 = 12.
    temp = ['--> O16 (123.456)\n']
    temp.append('   |--> O17 (12.3456)\n')
    temp.append('   |   |--> O18 (1.23456)\n')
    temp.append('   |--> HE4 (1111.0)\n')
    temp.append('--> C12 (12.0)\n')
    with open(filename, 'w') as tree:
        transmute._tree_log(d0, nuc0, N0, tree)
        transmute._tree_log(d1, nuc1, N1, tree)
        transmute._tree_log(d2, nuc2, N2, tree)
        transmute._tree_log(d11, nuc11, N11, tree)
        transmute._tree_log(d20, nuc20, N20, tree)
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        assert_equal(lines[i], temp[i])
    os.remove(filename)

"""Tests correct implementation of a transmutation with zero flux on an
isotope with a zero decay-constant."""
def test_zero_flux():
    nuc = nucname.zzaaam('FE56')
    t_sim = 100.
    phi = None
    tree = None
    tol = 1e-7
    out = transmute.transmute_core(nuc, t_sim, phi, tree, tol)
    assert_equal(out[nuc], 1)

#
# Run as script
#
if __name__ == "__main__":
    nose.main()
