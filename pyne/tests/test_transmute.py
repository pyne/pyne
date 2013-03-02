"""Transmute tests"""

import nose

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
    eA = linalg.expm(A)
    tayA = linalg.expm3(A,q=20)
    # Check the equivalence of _solve_decay_matrix
    assert_true(np.array_equal(eA, transA))


"""Ensures the flux vector size check is completing"""
def test_phi_size_check():
    nuc = nucname.zzaaam("U235")
    phi = np.arange(1.,175.)
    t_sim = 1.0
    tol = 1.0e-15
    assert_raises(SystemExit,transmute.decay,nuc,phi,t_sim,tol)


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

#
# Run as script
#
if __name__ == "__main__":
    nose.main()
