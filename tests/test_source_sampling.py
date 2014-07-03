from __future__ import print_function

import warnings
import itertools

from operator import itemgetter
from nose.tools import assert_true, assert_equal, assert_raises, with_setup, \
    assert_is, assert_is_instance, assert_in, assert_not_in, assert_almost_equal
from random import uniform

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
try:
    from itaps import iBase, iMesh, iMeshExtensions
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)

from pyne.mesh import Mesh, IMeshTag
from pyne.source_sampling import Sampler

def test_analog_single_hex():

    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]], mats = None)
    m.src = IMeshTag(1, float)
    m.src[0] = 1.0
    m.bias = IMeshTag(1, float)
    m.bias[0] = 1.0
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", np.array([0, 1]), False)

    num_samples = 100000
    num_divs = 2
    samples = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s[4], 1.0)
        samples[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs), 
                int(s[3]*num_divs)] += 1.0/num_samples
    
    for s in samples.flat:
        assert(abs(s*(num_divs)**4 - 1) < 0.3)

def test_analog_multiple_hex():
    m = Mesh(structured=True, structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 1]], mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = np.array([[1,1],[1,1],[1,1],[1,1]])
    m.mesh.save("sampling_mesh2.h5m")
    sampler = Sampler("sampling_mesh2.h5m", "src", np.array([0,0.5,1]), False)

    num_samples = 10000
    num_divs = 2
    samples = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s[4], 1.0)
        samples[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs), 
                int(s[3]*num_divs)] += 1.0/num_samples
    
    for s in samples.flat:
        assert(abs(s*(num_divs)**4 - 1) < 0.1)


def test_analog_single_tet():
    mesh = iMesh.Mesh()
    v1 = [0, 0, 0]
    v2 = [1, 0, 0]
    v3 = [0, 1, 0]
    v4 = [0, 0, 1]
    verts = mesh.createVtx([v1, v2, v3, v4])
    mesh.createEnt(iMesh.Topology.tetrahedron, verts)
    m = Mesh(structured=False, mesh=mesh)
    m.src = IMeshTag(1, float)
    m.src[:] = np.array([1])
    m.mesh.save("tet.h5m")
    center = m.ve_center(list(m.iter_ve())[0])

    subtets = [[center, v1, v2, v3], 
               [center, v1, v2, v4], 
               [center, v1, v3, v4], 
               [center, v2, v3, v4]]

    sampler = Sampler("tet.h5m", "src", np.array([0, 1]), False)
    num_samples = 10000
    samples = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s[4], 1.0)
        for i, tet in enumerate(subtets):
            if point_in_tet(tet, [s[0], s[1], s[2]]):
                samples[i] += 1.0/num_samples
    
    for s in samples:
        for a in samples:
            print(a)
        assert(abs(s*4 - 1) < 0.1)

def point_in_tet(t, p):

    matricies = [
    np.array( [[t[0][0], t[0][1], t[0][2], 1],
              [t[1][0], t[1][1], t[1][2], 1],
              [t[2][0], t[2][1], t[2][2], 1],
              [t[3][0], t[3][1], t[3][2], 1]]),
    np.array( [[p[0], p[1], p[2], 1],
              [t[1][0], t[1][1], t[1][2], 1],
              [t[2][0], t[2][1], t[2][2], 1],
              [t[3][0], t[3][1], t[3][2], 1]]),
    np.array( [[t[0][0], t[0][1], t[0][2], 1],
              [p[0], p[1], p[2], 1],
              [t[2][0], t[2][1], t[2][2], 1],
              [t[3][0], t[3][1], t[3][2], 1]]),
    np.array( [[t[0][0], t[0][1], t[0][2], 1],
              [t[1][0], t[1][1], t[1][2], 1],
              [p[0], p[1], p[2], 1],
              [t[3][0], t[3][1], t[3][2], 1]]),
    np.array( [[t[0][0], t[0][1], t[0][2], 1],
              [t[1][0], t[1][1], t[1][2], 1],
              [t[2][0], t[2][1], t[2][2], 1],
              [p[0], p[1], p[2], 1]])]

    determinates =[np.linalg.det(x) for x in matricies]
    print(determinates[0] - np.sum(determinates[1:]))
    return all(x >= 0 for x in determinates) or all(x < 0 for x in determinates)
