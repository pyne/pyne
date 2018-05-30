"""PyNE vta module tests"""
import os

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
                       assert_in, assert_true

from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne import vta
import numpy as np


def test_ray_voxel_traverse():
    x_bound = [0, 1, 2, 3, 4, 5, 6]
    y_bound = [0, 1, 2, 3, 4, 5, 6]
    z_bound = [0, 1]
    bounds = [x_bound, y_bound, z_bound]
    A = vta.Point(3, 5.5, 0.5)
    B = vta.Point(0.5, 0.5, 0.5)
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = [[2, 5, 0], [2, 4, 0], [2, 3, 0], [1, 3, 0],
               [1, 2, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0]]
    for i in range(len(exp_ans)):
        assert_array_equal(results[i], exp_ans[i])

def test_facet_voxel_rraverse():
    x_bound = [0, 1, 2, 3, 4, 5, 6]
    y_bound = [0, 1, 2, 3, 4, 5, 6]
    z_bound = [0, 1]
    bounds = [x_bound, y_bound, z_bound]
    A = vta.Point(0.5, 0.5, 0.5)
    B = vta.Point(5.5, 0.5, 0.5)
    C = vta.Point(3, 5.5, 0.5)
    results = vta._facet_voxel_traverse(A, B, C, bounds)
    exp_ans = [[2, 5, 0], [2, 4, 0], [2, 3, 0], [1, 3, 0],
               [1, 2, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],
               [2, 2, 0], [1, 0, 0], [2, 1, 0], [2, 0, 0],
               [3, 5, 0], [3, 4, 0], [3, 3, 0], [3, 2, 0],
               [3, 1, 0], [3, 0, 0], [4, 2, 0], [4, 1, 0],
               [4, 0, 0], [4, 3, 0], [5, 1, 0], [5, 0, 0]]
    for i in range(len(exp_ans)):
        assert_array_equal(results[i], exp_ans[i])


if __name__ == "__main__":
    nose.runmodule()

