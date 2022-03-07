import os

import numpy as np
import tables as tb

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_almost_equal,
    assert_true,
    assert_raises,
)
from numpy.testing import assert_array_equal, assert_array_almost_equal

import pyne
from pyne.xs.models import (
    partial_energy_matrix,
    partial_energy_matrix_mono,
    chi,
    alpha,
    k,
    m_n,
    beta,
    alpha_at_theta_0,
    alpha_at_theta_pi,
    one_over_gamma_squared,
    E_prime_min,
    sigma_s_const,
    sigma_s,
    phi_g,
    group_collapse,
    thermspect,
    fastspect,
)
from pyne.pyne_config import pyne_conf

nuc_data = pyne_conf.NUC_DATA_PATH

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")

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

    expected = np.array([[1.0, 0.0], [0.0, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_inc3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[0.5, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_inc4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_inc5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], [0.0, 0.4, 1.0, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_inc6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix_mono(E_g, E_n, 1)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], [0.0, 0.4, 1.0, 0.2]])

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

    expected = np.array([[1.0, 0.0], [0.0, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_dec3():
    E_g = np.array([7.5, 5.0, 1.25])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.5]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_dec4():
    E_g = np.array([10.0, 5.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_dec5():
    E_g = np.array([10.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[1.0, 1.0, 0.4, 0.0], [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix_dec6():
    E_g = np.array([8.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix_mono(E_g, E_n, -1)

    expected = np.array([[0.2, 1.0, 0.4, 0.0], [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix1():
    # tests dispach to inc
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix2():
    E_g = np.array([8.0, 4.0, 0.0])
    E_n = np.array([10.0, 7.5, 5.0, 2.5, 0.0])

    pem = partial_energy_matrix(E_g, E_n)

    expected = np.array([[0.2, 1.0, 0.4, 0.0], [0.0, 0.0, 0.6, 1.0]])

    assert_array_equal(pem, expected)


def test_partial_energy_matrix3():
    # tests monotonicity
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0][::-1])

    assert_raises(ValueError, partial_energy_matrix, E_g, E_n)


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


def test_group_collapse1():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([0.0, 2.0, 1.0, 0.5])
    sigma_n = np.array([1.0, 2.0, 3.0, 4.0])

    expected = np.array([2.0, 5.0 / 1.9])

    # First way of calling
    observed = group_collapse(sigma_n, phi_n, E_g=E_g, E_n=E_n)
    assert_array_almost_equal(observed, expected)

    # Second method of calling
    p_g = phi_g(E_g, E_n, phi_n)
    pem = partial_energy_matrix(E_g, E_n)
    observed = group_collapse(sigma_n, phi_n, phi_g=p_g, partial_energies=pem)
    assert_array_almost_equal(observed, expected)

    # bad call
    assert_raises(ValueError, group_collapse, sigma_n, phi_n)


def test_wgt_group_collapse1():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([0.0, 2.0, 1.0, 0.5])
    sigma_n = np.array([1.0, 2.0, 3.0, 4.0])
    wgts = np.array([0.00001, 0.00001, 0.00001, 0.00001])

    observed = group_collapse(sigma_n, phi_n, E_g=E_g, E_n=E_n, weights=wgts)
    expected = group_collapse(sigma_n, phi_n, E_g=E_g, E_n=E_n)
    assert_array_almost_equal(observed, expected)


#
# Test physical models
#


def test_chi():
    assert_equal(chi(0.0), 0.0)
    assert_equal(chi(1.0), 0.453 * np.exp(-1.036) * np.sinh(np.sqrt(2.29)))
    assert_equal(chi(10.0), 0.453 * np.exp(-10.36) * np.sinh(np.sqrt(22.9)))

    e = np.arange(50, dtype=float)
    assert_array_equal(chi(e), 0.453 * np.exp(-1.036 * e) * np.sinh(np.sqrt(2.29 * e)))


def test_alpha():
    assert_equal(1.0 / k, alpha(0.0, 1.0, 0.0, m_n, 1.0))

    assert_almost_equal(1.0, (1.0 / (12.0 * k)) / alpha(0.0, 1.0, 0.0, 12 * m_n, 1.0))

    assert_almost_equal(
        1.0, (1.5 / (12.0 * k)) / alpha(0.5, 1.0, np.pi / 2, 12 * m_n, 1.0)
    )

    assert_almost_equal(
        1.0, (1.5 / (12.0 * k * 2.0)) / alpha(0.5, 1.0, np.pi / 2, 12 * m_n, 2.0)
    )

    assert_almost_equal(
        1.0,
        ((1.5 - 2 * np.sqrt(0.5) * np.cos(np.pi / 4)) / (12.0 * k * 2.0))
        / alpha(0.5, 1.0, np.pi / 4, 12 * m_n, 2.0),
    )


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

    assert_array_almost_equal(
        alpha_at_theta_0(E_prime, 1.0, 12 * m_n, 2.0),
        alpha(E_prime, 1.0, 0.0, 12 * m_n, 2.0),
    )

    assert_array_almost_equal(
        alpha_at_theta_0(E_prime, E, 12 * m_n, 2.0),
        alpha(E_prime, E, 0.0, 12 * m_n, 2.0),
    )

    assert_array_almost_equal(
        alpha_at_theta_0(E_prime, E, M_A, 2.0), alpha(E_prime, E, 0.0, M_A, 2.0)
    )

    assert_array_almost_equal(
        alpha_at_theta_0(E_prime, E, M_A, T), alpha(E_prime, E, 0.0, M_A, T)
    )


def test_alpha_at_theta_pi():
    E_prime = np.linspace(0.5, 1.5, 101)
    E = np.linspace(0.75, 1.25, 101)

    M_A = np.linspace(1, 300, 101)
    T = np.linspace(1, 1800, 101)

    assert_array_almost_equal(
        alpha_at_theta_pi(E_prime, 1.0, 12 * m_n, 2.0),
        alpha(E_prime, 1.0, np.pi, 12 * m_n, 2.0),
    )

    assert_array_almost_equal(
        alpha_at_theta_pi(E_prime, E, 12 * m_n, 2.0),
        alpha(E_prime, E, np.pi, 12 * m_n, 2.0),
    )

    assert_array_almost_equal(
        alpha_at_theta_pi(E_prime, E, M_A, 2.0), alpha(E_prime, E, np.pi, M_A, 2.0)
    )

    assert_array_almost_equal(
        alpha_at_theta_pi(E_prime, E, M_A, T), alpha(E_prime, E, np.pi, M_A, T)
    )


def test_one_over_gamma_squared():
    E = np.logspace(-9, 2, 101)
    rcf = one_over_gamma_squared(E)
    expected = 1.0 - 2.0 * E / (931.46 * 1.0086649159700001)
    assert_array_almost_equal(rcf, expected)


def test_E_prime_min():
    assert_equal(0.0, E_prime_min(1.0, m_n))
    assert_equal(1.0 / 9.0, E_prime_min(1.0, 2.0 * m_n))
    assert_equal(4.0 / 9.0, E_prime_min(2.0, 2.0 * m_n))


def test_sigma_s_const():
    assert_equal(sigma_s_const(0.0), 0.0)
    assert_equal(sigma_s_const(0.5), 1e24 * np.pi)
    assert_equal(sigma_s_const(1.0), 4e24 * np.pi)


def test_sigma_s():
    # Probably could use some more testing
    E = np.logspace(-9, 2, 101)

    sig_s = sigma_s(E)
    assert_true((0.0 <= sig_s).all())
    assert_true((sig_s[1:] <= sig_s[:-1]).all())

    sig_s = sigma_s(E, 12.0, 13.0, 1900.0)
    assert_true((0.0 <= sig_s).all())
    assert_true((sig_s[1:] <= sig_s[:-1]).all())


def test_thermspect():
    e1 = [1.0e-6]
    e2 = [1.0]
    e3 = [0.9375e3, 1.5e4]
    assert_array_equal(np.asarray([1.0]), thermspect(np.asarray(e1)))
    assert_array_equal(np.asarray([1.0]), thermspect(np.asarray(e2)))
    assert_array_equal(np.asarray([0.8, 0.2]), thermspect(np.asarray(e3)))


def test_thermspect2():
    e1 = [0.9375e3, 1.5e4]
    phi = thermspect(np.asarray(e1), T=400.0, lower=1.0e4)
    assert_array_equal(np.asarray([0.0, 1.0]), phi)


def test_fastspect():
    e1 = np.asarray([1.0e-6])
    e2 = np.asarray([1.0])
    e3 = np.asarray([1.291e-5, 10])
    uno = np.asarray([1.0])
    assert_array_equal(uno, fastspect(e1))
    assert_array_equal(uno, fastspect(e2))
    v1 = np.asarray([0.04872966, 0.95127034])
    assert_array_almost_equal(v1, fastspect(e3))


def test_fastspect2():
    e1 = np.asarray([1.291e-5, 10])
    phi1 = fastspect(e1, T=1000.0, lower=1.0e-5)
    v1 = np.asarray([0.0032472, 0.9967528])
    assert_array_almost_equal(v1, phi1)
