from __future__ import print_function

import os
import warnings
import itertools

from operator import itemgetter
from nose.tools import assert_equal, with_setup, assert_almost_equal
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
from pyne.source_sampling import Sampler, AliasTable

def try_rm_file(filename):
    return lambda: os.remove(filename) if os.path.exists(filename) else None

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex():
    """This test tests that particles of sampled evenly within the phase-space 
    of a single mesh volume element with one energy group in an analog sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then 
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into. 
    """
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]], 
             mats = None)
    m.src = IMeshTag(1, float)
    m.src[0] = 1.0
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", np.array([0, 1]), False)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s[4], 1.0) # analog: all weights must be one
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs), 
              int(s[3]*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_multiple_hex():
    """This test tests that particle are sampled uniformly from a uniform source
    defined on eight mesh volume elements in two energy groups. This is done
    using the exact same method ass test_analog_multiple_hex.
    """
    m = Mesh(structured=True, 
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]], 
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = np.ones(shape=(8,2))
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", np.array([0, 0.5, 1]), False)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s[4], 1.0)
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs), 
              int(s[3]*num_divs)] += score
    
    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

@with_setup(None, try_rm_file('tet.h5m'))
def test_analog_single_tet():
    """This test tests uniform sampling within a single tetrahedron. This is
    done by dividing the tetrahedron in 4 smaller tetrahedrons and ensuring
    that each sub-tet is sampled equally.
    """
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
    num_samples = 5000
    score = 1.0/num_samples
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s[4], 1.0)
        for i, tet in enumerate(subtets):
            if point_in_tet(tet, [s[0], s[1], s[2]]):
                tally[i] += score
                break
    
    for t in tally:
        assert(abs(t - 0.25)/0.25 < 0.2)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_uniform():
    """This test tests that the uniform biasing scheme:
    1. Samples space uniformly. This is checked using the same method
       described in test_analog_single_hex().
    2. Adjusts weights accordingly. Sample calculations are provided in Case 1
       in the Theory Manual.
    """
    m = Mesh(structured=True, 
             structured_coords=[[0, 2.5, 10], [0, 10], [0, 10]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 3.0], [1.0, 4.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, True)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s[0] < 2.5:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.8) # hand calcs
            else:
              assert_almost_equal(s[4], 1.2) # hand calcs
        else:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.4) # hand calcs
            else:
              assert_almost_equal(s[4], 1.6) # hand calcs

        tally[int(s[0]*num_divs)/10, 
              int(s[1]*num_divs)/10, 
              int(s[2]*num_divs)/10]  += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias():
    """This test tests that a user-specified biasing scheme:
    1. Samples space uniformly according to the scheme.
    2. Adjusts weights accordingly. Sample calculations are provided in Case 2
       in the Theory Manual.
    """
    m = Mesh(structured=True, 
             structured_coords=[[0, 2.5, 10], [0, 10], [0, 10]], 
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 3.0], [1.0, 4.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.bias = IMeshTag(2, float)
    m.bias[:] = [[8.0, 3.0], [5.0, 8.0]]
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, "bias")

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s[0] < 2.5:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.625) # hand calcs
              tally[0] += score
            else:
              assert_almost_equal(s[4], 2.5) # hand calcs
              tally[1] += score
        else:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.5) # hand calcs
              tally[2] += score
            else:
              assert_almost_equal(s[4], 1.25) # hand calcs
              tally[3] += score

    expected_tally = [0.16, 0.06, 0.3, 0.48]
    for a, b in zip(tally, expected_tally):
       assert(abs(a-b)/b < 0.25)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias_spacial():
    """This test tests a user-specified biasing scheme for which the only 1
    bias group is supplied for a source distribution containing two energy 
    groups. This bias group is applied to both energy groups. In this test,
    the user-supplied bias distribution that was choosen, correspondes to 
    uniform sampling, so that results can be checked against Case 1 in the
    theory manual.
    """
    m = Mesh(structured=True, 
             structured_coords=[[0, 2.5, 10], [0, 10], [0, 10]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 3.0], [1.0, 4.0]]
    m.bias = IMeshTag(2, float)
    m.bias[:] = [1, 1]
    e_bounds = np.array([0, 0.5, 1.0])
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, "bias")

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s[0] < 2.5:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.8) # hand calcs
            else:
              assert_almost_equal(s[4], 1.2) # hand calcs
        else:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 0.4) # hand calcs
            else:
              assert_almost_equal(s[4], 1.6) # hand calcs

        tally[int(s[0]*num_divs)/10, 
              int(s[1]*num_divs)/10, 
              int(s[2]*num_divs)/10]  += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

def test_alias_table():
    """This tests that the AliasTable class produces samples in the ratios
    consistant with the supplied PDF.
    """
    pdf = np.array([0.1, 0.2, 0.7])
    at = AliasTable(pdf)
    num_samples = 50000
    score = 1.0/num_samples
    tally = np.zeros(shape=(3))

    for i in range(num_samples):
        s = at.sample_pdf(uniform(0, 1), uniform(0,1))
        tally[s] += score

    for i in range(0, 3):
       print(tally[i])
       assert(abs(tally[i] - pdf[i])/pdf[i] < 0.05)

def point_in_tet(t, p):
    """ This function determines if some point <p> lies within some tetrahedron
    <t> using the method described here:
    http://steve.hollasch.net/cgindex/geometry/ptintet.html
    """
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
    return all(x >= 0 for x in determinates) or all(x < 0 for x in determinates)
