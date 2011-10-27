import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal, assert_almost_equal, assert_true
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.xs.models import partial_energy_matrix, partial_energy_matrix_inc
from pyne.pyne_config import pyne_conf

nuc_data = pyne_conf.NUC_DATA_PATH


#
# Test Partial Energy Matrix
#

def test_partial_energy_matrix_inc1():
    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc_inc2():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[1.0, 0.0], 
                         [0.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[0.5, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 0.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_inc(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)    

