"""Transmute tests"""

from unittest import TestCase
import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
from pyne.material import Material
import numpy  as np
import scipy  as sp
import tables as tb
from pyne import transmute
from scipy import linalg


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

#
# Run as script
#
if __name__ == "__main__":
    nose.main()
