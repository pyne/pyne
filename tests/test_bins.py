"""PyNE bins tests"""
import os
import warnings

import nose

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_in,
    assert_true,
)

from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne import bins
import numpy as np


def test_ninespace():
    obs = bins.ninespace(0.9, 0.9999, 4)
    exp = np.array([0.9, 0.99, 0.999, 0.9999])
    assert_array_almost_equal(obs, exp)


def test_stair_step():
    x = [0.1, 1.0, 10.0, 100.0]
    y = [2.0, 3.0, 4.0]

    xobs, yobs = bins.stair_step(x, y)
    xexp = [0.1, 1.0, 1.0, 10.0, 10.0, 100.0]
    yexp = [2.0, 2.0, 3.0, 3.0, 4.0, 4.0]

    assert_equal(len(xobs), len(yobs))
    assert_equal(len(yobs), 2 * len(y))
    assert_array_almost_equal(xobs, xexp)
    assert_array_almost_equal(yobs, yexp)


def check_pointwise_linear_collapse(x_g, x, y, exp):
    obs = bins.pointwise_linear_collapse(x_g, x, y)
    assert_equal(len(exp), len(obs))
    assert_equal(len(x_g) - 1, len(obs))
    assert_array_almost_equal(exp, obs)


def test_pointwise_linear_collapse():
    cases = [
        [
            np.array([0.0, 1.0, 2.0]),
            np.linspace(0.0, 2.0, 101),
            np.ones(101),
            [1.0, 1.0],
        ],
        [
            np.array([0.0, 1.0, 2.0]),
            np.linspace(0.0, 2.0, 101),
            np.linspace(0.0, 2.0, 101),
            [0.5, 1.5],
        ],
        [np.array([-1.0, 0.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [3.5]],
        [np.array([0.0, 1.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [4.5]],
        [np.array([1.0, 2.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [5.5]],
        [np.array([0.0, 1.0]), np.array([0.0, 2.0]), np.array([0.0, 2.0]), [0.5]],
        [np.array([0.5, 1.5]), np.array([0.0, 2.0]), np.array([0.0, 2.0]), [1.0]],
        [
            np.array([0.0, 1.0]),
            np.array([0.0, 0.75]),
            np.array([0.0, 1.0]),
            [2.0 / 3.0],
        ],
        [np.array([1.0, 2.0]), np.array([0.0, 2.0]), np.array([4.0, 6.0]), [5.5]],
    ]
    for x_g, x, y, exp in cases:
        yield check_pointwise_linear_collapse, x_g, x, y, exp
        yield check_pointwise_linear_collapse, x_g[::-1], x[::-1], y[::-1], exp[::-1]


def test_point_collapse():
    """Verify all types of interpolation work for pointwise_collapse."""

    x_g = np.array([1.0, 10.0])
    x = np.array([1.0, 2.0, 20.0, np.power(10.0, 2.6)])
    y = np.array([100.0, 100.0, 1000.0, 1.0])

    lin_lin = bins.pointwise_collapse(x_g, x, y)
    lin_log = bins.pointwise_collapse(x_g, x, y, logy=True)
    log_lin = bins.pointwise_collapse(x_g, x, y, logx=True)
    log_log = bins.pointwise_collapse(x_g, x, y, log=True)

    # Verified via hand calcs
    exp_lin_lin = np.array([277.77777778])
    exp_lin_log = np.array([157.59080145])
    exp_log_lin = np.array([319.85158013])
    exp_log_log = np.array([175.50097498])

    assert np.allclose(lin_lin, exp_lin_lin)
    assert np.allclose(lin_log, exp_lin_log)
    assert np.allclose(log_lin, exp_log_lin)
    assert np.allclose(log_log, exp_log_log)


def test_point_collase_raises():

    x_g = np.array([0.1, 1.0, 2.0])
    x = np.array([0.1, 1.0, 2.0])
    y = np.array([0.1, 1.0, 2.0])

    with assert_raises(ValueError):
        bins.pointwise_collapse(x_g, np.array([2.0, 1.0, 2.0]), y)

    with assert_raises(ValueError):
        bins.pointwise_collapse(np.array([-1, 1.0, 2.0]), x, y, logx=True)

    with assert_raises(ValueError):
        bins.pointwise_collapse(x_g, x, np.array([0, 1.0, 2.0]), log=True)


if __name__ == "__main__":
    nose.runmodule()
