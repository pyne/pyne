"""PyNE bins tests"""
import os

import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in, \
                       assert_true

from numpy.testing import assert_array_equal, assert_array_almost_equal

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

if __name__ == "__main__":
    nose.main()

