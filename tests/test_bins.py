"""PyNE bins tests"""
import os
import warnings

import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in, assert_true

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
    yexp = [2.0, 2.0, 3.0, 3.0,  4.0,  4.0]

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
        [np.array([0.0, 1.0, 2.0]), np.linspace(0.0, 2.0, 101), 
         np.ones(101), [1.0, 1.0]],
        [np.array([0.0, 1.0, 2.0]), np.linspace(0.0, 2.0, 101), 
         np.linspace(0.0, 2.0, 101), [0.5, 1.5]],
        [np.array([-1.0, 0.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [3.5]],
        [np.array([0.0, 1.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [4.5]],
        [np.array([1.0, 2.0]), np.array([-1.0, 2.0]), np.array([3.0, 6.0]), [5.5]],
        [np.array([0.0, 1.0]), np.array([0.0, 2.0]), np.array([0.0, 2.0]), [0.5]],
        [np.array([0.5, 1.5]), np.array([0.0, 2.0]), np.array([0.0, 2.0]), [1.0]],
        [np.array([0.0, 1.0]), np.array([0.0, 0.75]), np.array([0.0, 1.0]), [2.0/3.0]],
        [np.array([1.0, 2.0]), np.array([0.0, 2.0]), np.array([4.0, 6.0]), [5.5]],
        ]
    for x_g, x, y, exp in cases:
        yield check_pointwise_linear_collapse, x_g, x, y, exp
        yield check_pointwise_linear_collapse, x_g[::-1], x[::-1], y[::-1], exp[::-1]

def test_linear_collapse():

    x_g = np.array([0.00, 0.01, 0.03, 0.07, 0.10, 0.11, 0.12, 0.21, 0.41, 0.42, 0.5, 0.8, 0.91,
           0.92, 1.1, 1.2])
    x = np.array([0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
    y = np.array([0.85, 0.75, 0.63, 0.23, 0.26, 0.24, 0.38, 0.56, 0.80, 0.77, 0.97])

    lin_lin = bins.pointwise_collapse(x_g, x, y, (False, False))
    #lin_log = bins.pointwise_collapse(x_g, x, y, (False, True))
    #log_lin = bins.pointwise_collapse(x_g, x, y, (True, False))
    #log_log = bins.pointwise_collapse(x_g, x, y, (True, True))

    exp_lin_lin = np.array([0.8499998972, 0.8499980839, 0.8285738601,
                            0.7736293872, 0.7438675679, 0.7319192405, 0.6704822882, 0.3486381283,
                            0.2570217582, 0.2482952453, 0.4622711939, 0.7847809729, 0.8002478617,
                            0.8983716176, 0.9699962467])
    exp_lin_log = np.array([0.8499998972, 0.8499980839, 0.8276281356,
                            0.7726665077, 0.743357055, 0.7305636342, 0.6675312678, 0.3275902384,
                            0.256920334, 0.2481426929, 0.4563733894, 0.7846051604, 0.7973832307,
                            0.8945183403, 0.9699962467])
    exp_log_lin = np.array([0.8499998972, 0.8499980839, 0.8229431289,
                            0.7687535679, 0.7414321701, 0.7257565196, 0.6632651436, 0.3403854023,
                            0.2567265603, 0.2478736595, 0.4645827573, 0.784555933, 0.8016253225,
                            0.9001178593, 0.9699962467])
    exp_log_log = np.array([0.8499998972, 0.8499980839, 0.8218771941,
                            0.7679207585, 0.7407403275, 0.7240572047, 0.6601639557, 0.3194561454,
                            0.2566165951, 0.2477226101, 0.4586655975, 0.784373088, 0.7986550456,
                            0.8962838244, 0.9699962467])
    print(lin_lin)
    print(exp_lin_lin)
    assert(np.allclose(lin_lin, exp_lin_lin))
    #assert(np.allclose(lin_log, exp_lin_log))
    #assert(np.allclose(log_lin, exp_log_lin))
    #assert(np.allclose(log_log, exp_log_log))

if __name__ == "__main__":
    nose.runmodule()

