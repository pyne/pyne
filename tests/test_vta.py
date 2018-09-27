"""PyNE vta module tests"""
import os

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
                       assert_in, assert_true

from numpy.testing import assert_array_equal, assert_array_almost_equal

import numpy as np
from pyne import vta


def test_ray_voxel_traverse_2D():
    x_bound = [0, 1, 2, 3, 4, 5, 6]
    y_bound = [0, 1, 2, 3, 4, 5, 6]
    z_bound = [0, 1]
    bounds = [x_bound, y_bound, z_bound]
    A = np.array([3, 5.5, 0.5])
    B = np.array([0.5, 0.5, 0.5])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(2, 5, 0), (2, 4, 0), (2, 3, 0), (1, 3, 0),
               (1, 2, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0)])
#    for i in range(len(exp_ans)):
#        assert_array_equal(results[i], exp_ans[i])
    assert_equal(exp_ans, results)

def test_ray_voxel_traverse_3D():
    """
    This unit test uses a cubit with the size of 2 * 2 * 2.
    Start point: (0.3, 0.4, 0.5), end point: (2, 2, 2).
    Direction: [0.62, 0.58, 0.54]
    Ray traverse history:
    Point: [0.3, 0.4, 0.5] -> [0.88, 0.94, 1] -> [0.94, 1, 1.06] -> [1, 1.06, 1.11]
    Voxel: (0, 0, 0)       -> (0, 0, 1)       -> (0, 1, 1)       -> (1, 1, 1)
    """
    x_bound = [0, 1, 2]
    y_bound = [0, 1, 2]
    z_bound = [0, 1, 2]
    bounds = [x_bound, y_bound, z_bound]
    A = np.array([0.3, 0.4, 0.5])
    B = np.array([2, 2, 2])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
#    for i in range(len(exp_ans)):
#        assert_array_equal(results[i], exp_ans[i])
    assert_equal(exp_ans, results)

def test_facet_voxel_rraverse():
    x_bound = [0, 1, 2, 3, 4, 5, 6]
    y_bound = [0, 1, 2, 3, 4, 5, 6]
    z_bound = [0, 1]
    bounds = [x_bound, y_bound, z_bound]
    A = np.array([0.5, 0.5, 0.5])
    B = np.array([5.5, 0.5, 0.5])
    C = np.array([3, 5.5, 0.5])
#    import pdb; pdb.set_trace()
    results = vta._facet_voxel_traverse(A, B, C, bounds)
    exp_ans = set([(2, 5, 0), (2, 4, 0), (2, 3, 0), (1, 3, 0),
               (1, 2, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0),
               (2, 2, 0), (1, 0, 0), (2, 1, 0), (2, 0, 0),
               (3, 5, 0), (3, 4, 0), (3, 3, 0), (3, 2, 0),
               (3, 1, 0), (3, 0, 0), (4, 2, 0), (4, 1, 0),
               (4, 0, 0), (4, 3, 0), (5, 1, 0), (5, 0, 0)])
#    for i in range(len(exp_ans)):
#        assert_array_equal(results[i], exp_ans[i])
    assert_equal(exp_ans, results)


if __name__ == "__main__":
    nose.runmodule()

