"""PyNE vta module tests"""
import os

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
                       assert_in, assert_true

from numpy.testing import assert_array_equal, assert_array_almost_equal

import numpy as np
from pyne import vta

def test_find_voxel_idx_1d():
    bounds_1d = np.array([1.0, 2.0, 3.0])

    # tests for right dir
    vec_1d = 0.5
    # left boundary right dir
    cor = 1.0
    exp_ans = 0
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)
    # middle, right dir
    cor = 1.5
    exp_ans = 0
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)
    # right boundary, right dir
    cor = 3.0
    exp_ans = -1
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)

    # tests for left dir
    vec_1d = -0.5
    # left boundary, left dir
    cor = 1.0
    exp_ans = -1
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)
    # middle, left dir
    cor = 1.5
    exp_ans = 0
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)
    # right boundary, left dir
    cor = 3.0
    exp_ans = 1
    result = vta._find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    assert_equal(exp_ans, result)

def test_ray_voxel_traverse_2D():
    bounds = np.array([[0, 1, 2, 3, 4, 5, 6],
                       [0, 1, 2, 3, 4, 5, 6],
                       [0, 1]])
    A = np.array([3, 5.5, 0.5])
    B = np.array([0.5, 0.5, 0.5])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(2, 5, 0), (2, 4, 0), (2, 3, 0), (1, 3, 0),
               (1, 2, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0)])
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
    bounds = np.array([[0, 1, 2],
                       [0, 1, 2],
                       [0, 1, 2]])
    A = np.array([0.3, 0.4, 0.5])
    B = np.array([2, 2, 2])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
    assert_equal(exp_ans, results)

    # reverse the ray direction
    A = np.array([2, 2, 2])
    B = np.array([0.3, 0.4, 0.5])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
    assert_equal(exp_ans, results)

    # set start and end outside the mesh
    A = np.array([-0.32, -0.18, -0.04])
    B = np.array([2.62, 2.58, 2.54])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
    assert_equal(exp_ans, results)

    # ray has no intersection with the mesh
    A = np.array([-1.0, -1.0, -1.0])
    B = np.array([3.0, -1.0, -1.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

    # ray on the outside boundary surface of the mesh
    A = np.array([0.0, 0.0, 0.0])
    B = np.array([0.0, 0.0, 2.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

    # ray inside the mesh, on the boundary of x and y direction
    A = np.array([1.0, 1.0, 0.0])
    B = np.array([1.0, 1.0, 2.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

    # ray on boundary in x direction
    A = np.array([1.0, 0.5, 0.0])
    B = np.array([1.0, 0.5, 2.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

    # ray on boundary in y direction
    A = np.array([0.5, 1.0, 0.0])
    B = np.array([0.5, 1.0, 2.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

    # ray on boundary in z direction
    A = np.array([0.5, 0.0, 1.0])
    B = np.array([0.5, 2.0, 1.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
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
    assert_equal(exp_ans, results)

    # reverse the ray direction
    A = np.array([2, 2, 2])
    B = np.array([0.3, 0.4, 0.5])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
    assert_equal(exp_ans, results)

    # set start and end outside the mesh
    A = np.array([-0.32, -0.18, -0.04])
    B = np.array([2.62, 2.58, 2.54])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set([(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])
    assert_equal(exp_ans, results)

    # ray has no intersection with the mesh
    A = np.array([-1.0, -1.0, -1.0])
    B = np.array([3.0, -1.0, -1.0])
    results = vta._ray_voxel_traverse(bounds, A, B)
    exp_ans = set()
    assert_equal(exp_ans, results)

def test_facet_voxel_rraverse():
    bounds = np.array([[0, 1, 2, 3, 4, 5, 6],
                       [0, 1, 2, 3, 4, 5, 6],
                       [0, 1]])
    A = np.array([0.5, 0.5, 0.5])
    B = np.array([5.5, 0.5, 0.5])
    C = np.array([3, 5.5, 0.5])
    results = vta._facet_voxel_traverse(A, B, C, bounds)
    exp_ans = set([(2, 5, 0), (2, 4, 0), (2, 3, 0), (1, 3, 0),
               (1, 2, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0),
               (2, 2, 0), (1, 0, 0), (2, 1, 0), (2, 0, 0),
               (3, 5, 0), (3, 4, 0), (3, 3, 0), (3, 2, 0),
               (3, 1, 0), (3, 0, 0), (4, 2, 0), (4, 1, 0),
               (4, 0, 0), (4, 3, 0), (5, 1, 0), (5, 0, 0)])
    assert_equal(exp_ans, results)

def test_is_tri_aabb_intersects():
    box_center = np.array([0.0, 0.0, 0.0])
    box_extends = np.array([1.0, 1.0, 1.0])

    # triangle with all points inside the cube
    triangle = np.array([[0.0, 0.0, 0.0],
                         [0.5, 0.2, 0.3],
                         [0.3, 0.2, 0.4]])
    exp_ans = True
    result = vta._is_tri_intersects_box(triangle, box_center, box_extends)
    assert_equal(exp_ans, result)

    # triangle with one point inside the cube
    triangle = np.array([[0.0, 0.0, 0.0],
                         [2.0, 0.0, 0.0],
                         [2.0, 3.0, 0.0]])
    exp_ans = True
    result = vta._is_tri_intersects_box(triangle, box_center, box_extends)
    assert_equal(exp_ans, result)

    # triangle with all points outside the cube, but intersects with cube
    triangle = np.array([[-2.0, -2.0, -2.0],
                         [2.0, 2.0, 2.0],
                         [2.0, 2.0, 0.0]])
    exp_ans = True
    result = vta._is_tri_intersects_box(triangle, box_center, box_extends)
    assert_equal(exp_ans, result)

    # triangle with points outside the cube, but does not intersect with cube
    triangle = np.array([[2.0, 0.0, 0.0],
                         [2.0, 2.0, 0.0],
                         [3.0, 0.0, 0.0]])
    exp_ans = False
    result = vta._is_tri_intersects_box(triangle, box_center, box_extends)
    assert_equal(exp_ans, result)

if __name__ == "__main__":
    nose.runmodule()

