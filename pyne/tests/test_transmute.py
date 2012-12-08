"""Transmute tests"""

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

from numpy.testing import dec

from pyne.material import Material
import numpy  as np
import scipy  as sp
import tables as tb
from pyne import transmute, nucname, data, nuc_data
from scipy import linalg


"""Example test for future reference"""
"""def test_decay1():
    mat = Material({'C14': 1.0, 'N14': 0.0, 'C13': 0.0})
    obs = transmute.decay(mat, 1.0)

    assert False
"""

"""Tests _solve_decay_matrix calculating the matrix exponential"""
def test_expm():
    A = np.zeros((3,3))
    np.fill_diagonal(A,[-1,-2,-3])
    transA = transmute._solve_decay_matrix(A)
    eA = linalg.expm(A)
    tayA = linalg.expm3(A,q=20)
    # Check the equivalence of _solve_decay_matrix
    assert_equal(eA.all(),transA.all())
    # Check convergence to high level Taylor Series solution
    assert_almost_equal(tayA.all(),transA.all())


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
    daughters = transmute._get_daughters(nuc)
    assert_equal(daughters, daughtersTest)
    data_table.close()


#
# Run as script
#
if __name__ == "__main__":
    nose.main()
