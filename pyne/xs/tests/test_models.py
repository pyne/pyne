import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal, assert_almost_equal, assert_true, \
                       assert_raises
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.xs.models import partial_energy_matrix, partial_energy_matrix_mono, chi, alpha, k, \
    m_n, beta, alpha_at_theta_0, alpha_at_theta_pi, one_over_gamma_squared, E_prime_min, sigma_s_const, \
    phi_g
from pyne.pyne_config import pyne_conf

nuc_data = pyne_conf.NUC_DATA_PATH


#
# Test Partial Energy Matrix
#

def test_partial_energy_matrix_inc1():
    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc2():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 0.0], 
                         [0.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[0.5, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 0.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_inc6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_dec1():
    E_g = np.array([10.0, 0.0])
    E_n = np.array([10.0, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_dec2():
    E_g = np.array([10.0, 5.0, 0.0])
    E_n = np.array([10.0, 5.0, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0, 0.0], 
                         [0.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_dec3():
    E_g = np.array([7.5, 5.0, 1.25])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[0.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 0.5]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_dec4():
    E_g = np.array([10.0, 5.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix_dec5():
    E_g = np.array([10.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0, 1.0, 0.4, 0.0], 
                         [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)  


def test_partial_energy_matrix_dec6():
    E_g = np.array([8.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[0.2, 1.0, 0.4, 0.0], 
                         [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix1():
    # tests dispach to inc
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix2():
    E_g = np.array([8.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix(E_g, E_n)

    expected = np.array([[0.2, 1.0, 0.4, 0.0], 
                         [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix3():
    # tests monotonicity
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0][::-1])

    assert_raises(ValueError,  partial_energy_matrix, E_g, E_n)


#
# Test Group Collapse
#

def test_phi_g1():
    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    phi_n = np.ones(1)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0])

    assert_array_equal(observed, expected)


def test_phi_g2():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    phi_n = np.ones(2)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(observed, expected)


def test_phi_g3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.5, 1.0])

    assert_array_equal(observed, expected)

def test_phi_g4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([2.0, 2.0])

    assert_array_equal(observed, expected)


def test_phi_g5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 2.4])

    assert_array_equal(observed, expected)


def test_phi_g6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 1.6])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(phi_g, expected)

def test_phi_g6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 1.6])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(observed, expected)


def test_phi_g7():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([0.0, 2.0, 1.0, 0.5])

    observed = phi_g(E_g, E_n, phi_n)

    expected = np.array([1.2, 1.9])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(observed, expected)




#
# Test physical models
#

def test_chi():
    assert_equal(chi(0.0), 0.0)
    assert_equal(chi(1.0), 0.453 * np.exp(-1.036) * np.sinh(np.sqrt(2.29)))
    assert_equal(chi(10.0), 0.453 * np.exp(-10.36) * np.sinh(np.sqrt(22.9)))

    e = np.arange(50, dtype=float)
    assert_array_equal(chi(e), 0.453 * np.exp(-1.036*e) * np.sinh(np.sqrt(2.29*e)))


def test_alpha():
    assert_equal(1.0 / k, alpha(0.0, 1.0, 0.0, m_n, 1.0))

    assert_almost_equal(1.0, (1.0 / (12.0*k)) / alpha(0.0, 1.0, 0.0, 12*m_n, 1.0))

    assert_almost_equal(1.0, (1.5 / (12.0*k)) / alpha(0.5, 1.0, np.pi/2, 12*m_n, 1.0))

    assert_almost_equal(1.0, (1.5 / (12.0*k*2.0)) / alpha(0.5, 1.0, np.pi/2, 12*m_n, 2.0))

    assert_almost_equal(1.0, ((1.5 - 2*np.sqrt(0.5)*np.cos(np.pi/4)) / (12.0*k*2.0)) / alpha(0.5, 1.0, np.pi/4, 12*m_n, 2.0))


def test_beta():
    assert_equal((1.0 / k), beta(2.0, 1.0, 1.0))

    assert_equal((2.0 / k), beta(3.0, 1.0, 1.0))

    assert_equal((2.0 / (2.0 * k)), beta(3.0, 1.0, 2.0))

    assert_equal((-1.0 / (2.0 * k)), beta(0.0, 1.0, 2.0))


def test_alpha_at_theta_0():
    E_prime = np.linspace(0.5, 1.5, 101)
    E = np.linspace(0.75, 1.25, 101)

    M_A = np.linspace(1, 300, 101)
    T = np.linspace(1, 1800, 101)

    assert_array_almost_equal(alpha_at_theta_0(E_prime, 1.0, 12*m_n, 2.0),  alpha(E_prime, 1.0, 0.0, 12*m_n, 2.0))

    assert_array_almost_equal(alpha_at_theta_0(E_prime, E, 12*m_n, 2.0),  alpha(E_prime, E, 0.0, 12*m_n, 2.0))

    assert_array_almost_equal(alpha_at_theta_0(E_prime, E, M_A, 2.0),  alpha(E_prime, E, 0.0, M_A, 2.0))

    assert_array_almost_equal(alpha_at_theta_0(E_prime, E, M_A, T),  alpha(E_prime, E, 0.0, M_A, T))


def test_alpha_at_theta_pi():
    E_prime = np.linspace(0.5, 1.5, 101)
    E = np.linspace(0.75, 1.25, 101)

    M_A = np.linspace(1, 300, 101)
    T = np.linspace(1, 1800, 101)

    assert_array_almost_equal(alpha_at_theta_pi(E_prime, 1.0, 12*m_n, 2.0),  alpha(E_prime, 1.0, np.pi, 12*m_n, 2.0))

    assert_array_almost_equal(alpha_at_theta_pi(E_prime, E, 12*m_n, 2.0),  alpha(E_prime, E, np.pi, 12*m_n, 2.0))

    assert_array_almost_equal(alpha_at_theta_pi(E_prime, E, M_A, 2.0),  alpha(E_prime, E, np.pi, M_A, 2.0))

    assert_array_almost_equal(alpha_at_theta_pi(E_prime, E, M_A, T),  alpha(E_prime, E, np.pi, M_A, T))


def test_one_over_gamma_squared():
    E = np.logspace(-9, 2, 101)
    rcf = one_over_gamma_squared(E)
    expected = 1.0 - 2.0 * E / (931.46 * 1.0086649159700001)
    assert_array_almost_equal(rcf, expected)


def test_E_prime_min():
    assert_equal(0.0, E_prime_min(1.0, m_n))
    assert_equal(1.0/9.0, E_prime_min(1.0, 2.0*m_n))
    assert_equal(4.0/9.0, E_prime_min(2.0, 2.0*m_n))


def test_sigma_s_const():
    assert_equal(sigma_s_const(0.0), 0.0)
    assert_equal(sigma_s_const(0.5), 1E24 * np.pi)
    assert_equal(sigma_s_const(1.0), 4E24 * np.pi)
