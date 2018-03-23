import os
import warnings
import itertools

from operator import itemgetter
from nose.tools import assert_equal, with_setup, assert_almost_equal, assert_raises
from random import uniform, seed

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
try:
    from itaps import iBase, iMesh, iMeshExtensions
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne.mesh import Mesh, IMeshTag
from pyne.source_sampling import Sampler, AliasTable

# Define modes
DEFAULT_ANALOG = 0
DEFAULT_UNIFORM = 1
DEFAULT_USER = 2
SUBVOXEL_ANALOG = 3
SUBVOXEL_UNIFORM = 4
SUBVOXEL_USER = 5

def try_rm_file(filename):
    return lambda: os.remove(filename) if os.path.exists(filename) else None

@with_setup(None, try_rm_file('tet.h5m'))
def test_single_tet_tag_names_map():
    """This test tests uniform sampling within a single tetrahedron. This is
    done by dividing the tetrahedron in 4 smaller tetrahedrons and ensuring
    that each sub-tet is sampled equally.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.bias = IMeshTag(2, float)
    m.bias[:] = [[1.0, 2.0], [3.0, 3.0]]
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)

    # right condition
    tag_names = {"src_tag_name": "src"}
    e_bounds = np.array([0, 1])
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_ANALOG)

    # src_tag_name not given
    tag_names = {}
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_ANALOG)
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_UNIFORM)

    # bias_tag_name not given
    tag_names = {"src_tag_name": "src"}
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_USER)

    # subvoxel r2s source.h5m used for r2s calculation
    cell_fracs = np.zeros(2, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0), (1, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_ANALOG)
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_UNIFORM)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name":"cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    assert_raises(ValueError, Sampler, filename, tag_names, e_bounds, DEFAULT_USER)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex():
    """This test tests that particles of sampled evenly within the phase-space 
    of a single mesh volume element with one energy group in an analog sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then 
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into. 
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]], 
             mats = None)
    m.src = IMeshTag(1, float)
    m.src[0] = 1.0
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), DEFAULT_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0) # analog: all weights must be one
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs), 
              int(s.e*num_divs)] += score

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
    seed(1953)
    m = Mesh(structured=True, 
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]], 
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = np.ones(shape=(8,2))
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src"}
    sampler = Sampler(filename, tag_names, np.array([0, 0.5, 1]), DEFAULT_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s.w, 1.0)
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs), 
              int(s.e*num_divs)] += score
    
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
    seed(1953)
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
    filename = "tet.h5m"
    m.mesh.save(filename)
    center = m.ve_center(list(m.iter_ve())[0])

    subtets = [[center, v1, v2, v3], 
               [center, v1, v2, v4], 
               [center, v1, v3, v4], 
               [center, v2, v3, v4]]
    tag_names = {"src_tag_name": "src"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), DEFAULT_ANALOG)
    num_samples = 5000
    score = 1.0/num_samples
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s.w, 1.0)
        for i, tet in enumerate(subtets):
            if point_in_tet(tet, [s.x, s.y, s.z]):
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
    seed(1953)
    m = Mesh(structured=True, 
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src"}
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_UNIFORM)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4)) # number of phase space groups
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3.0:
            assert_almost_equal(s.w, 0.7) # hand calcs
        else:
            assert_almost_equal(s.w, 2.8) # hand calcs

        spatial_tally[int(s.x*num_divs/3.5), 
                      int(s.y*num_divs/1.0), 
                      int(s.z*num_divs/1.0)]  += score

        if s.x < 3 and s.e < 0.5:
            e_tally[0] += score
        elif s.x < 3 and s.e > 0.5:
            e_tally[1] += score
        if s.x > 3 and s.e < 0.5:
            e_tally[2] += score
        if s.x > 3 and s.e > 0.5:
            e_tally[3] += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(spatial_tally, i)[j,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

    expected_e_tally = [4./7, 2./7, 3./28, 1./28] # hand calcs
    for i in range(4):
        assert(abs(e_tally[i] - expected_e_tally[i]) \
               /expected_e_tally[i] < 0.1)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_single_subvoxel_analog():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element (also a sub-voxel) with one energy group
    in an analog sampling scheme. This done by dividing each dimension
    (x, y, z, E) in half, then sampling particles and tallying on the basis of
    which of the 2^4 = 16 regions of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(1, float)
    m.src[0] = 1.0
    cell_fracs = np.zeros(1, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), SUBVOXEL_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0) # analog: all weights must be one
        assert_equal(s.c, 11) # analog: the cell number
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_multiple_subvoxel_analog():
    """This test tests that particles of sampled analog within the phase-space
    of a single mesh volume element but multiple sub-voxels with one energy
    group in an analog sampling scheme. Then sampling particles and tallying
    the particles and check the probability of particles born in each
    sub-voxel and the cell_number.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(3, float)
    m.src[:] = np.empty(shape=(1, 3))
    m.src[0] = [0, 0.2, 0.8]
    cell_fracs = np.zeros(3, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.3, 0.0), (0, 12, 0.3, 0.0), (0, 13, 0.4, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), SUBVOXEL_ANALOG)
    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0) # analog: all weights must be one
        if s.c == 11:
            tally[0] += score
        elif s.c == 12:
            tally[1] += score
        elif s.c == 13:
            tally[2] += score

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.158)/0.158 < 0.05)
    assert(abs(tally[2] - 0.842)/0.842 < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_multiple_hex_multiple_subvoxel_analog():
    """This test tests that particle are sampled analog from a uniform source
    defined on eight mesh volume elements in two energy groups.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = np.ones(shape=(8,2))
    cell_fracs = np.zeros(8, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 1, 1.0, 0.0), (1, 2, 1.0, 0.0), (2, 3, 1.0, 0.0),
                     (3, 4, 1.0, 0.0), (4, 5, 1.0, 0.0), (5, 6, 1.0, 0.0),
                     (6, 7, 1.0, 0.0), (7, 8, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 0.5, 1]), SUBVOXEL_ANALOG)
    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s.w, 1.0)
        assert_equal(s.c, 4*int(s.x*num_divs) + 2*int(s.y*num_divs)
                     + int(s.z*num_divs) + 1)
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_subvoxel_uniform():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an uniform sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(1, float)
    m.src[0] = 1.0
    cell_fracs = np.zeros(1, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), SUBVOXEL_UNIFORM)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0) # analog: all weights must be one
        assert_equal(s.c, 11) # analog: the cell number
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
               int(s.e*num_divs)] += score

     # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
     # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_multiple_subvoxel_uniform():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an uniform sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(3, float)
    m.src[:] = np.empty(shape=(1, 3))
    m.src[0] = [0, 0.2, 0.8]
    cell_fracs = np.zeros(3, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.3, 0.0), (0, 12, 0.3, 0.0), (0, 13, 0.4, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), SUBVOXEL_UNIFORM)
    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.c == 11:
            tally[0] += score
        if s.c == 12:
            tally[1] += score
            assert(abs(s.w - 0.369)/0.369 < 0.05) # analog: all weights must be one
        if s.c == 13:
            tally[2] += score
            assert(abs(s.w - 1.475)/1.475 < 0.05)

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.428) < 0.05)
    assert(abs(tally[2] - 0.572) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_multiple_hex_multiple_subvoxel_uniform():
    """This test tests that particle are sampled uniformly from a uniform source
    defined on eight mesh volume elements in two energy groups.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = np.empty(shape=(8,2), dtype=float)
    m.src[:] = [[0,0], [1,0], [0,0], [2,0],
                [0,0], [3,0], [0,0], [4,0]]
    cell_fracs = np.zeros(8, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 0, 1.0, 0.0), (1, 1, 1.0, 0.0), (2, 2, 1.0, 0.0),
                     (3, 3, 1.0, 0.0), (4, 4, 1.0, 0.0), (5, 5, 1.0, 0.0),
                     (6, 6, 1.0, 0.0), (7, 7, 1.0, 0.0)]
    empty_cells = [0, 2, 4, 6]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 0.5, 1]), SUBVOXEL_UNIFORM)
    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 8
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        # check the cell_number
        assert_equal(s.c, 4*int(s.x*num_divs) + 2*int(s.y*num_divs)
                     + int(s.z*num_divs))
        # check the weight of each subvoxel
        if s.c not in empty_cells:
            # weight for cell 1, 3, 5, 7 should be: 0.4, 0.8, 1.2, 1.6
            exp_w = (s.c + 1) / 2 * 0.4
            out_w = s.w
            assert(abs(out_w - exp_w)/exp_w < 0.05) # hand calculate
        # count the tally
        tally[s.c] += score

    # check the real sample rate
    for i, item in enumerate(tally):
        if i not in empty_cells:
            assert(abs(item - 0.25)/0.25 < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias():
    """This test tests that a user-specified biasing scheme:
    1. Samples space uniformly according to the scheme.
    2. Adjusts weights accordingly. Sample calculations are provided in Case 2
       in the Theory Manual.
    """
    seed(1953)
    m = Mesh(structured=True, 
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]], 
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.bias = IMeshTag(2, float)
    m.bias[:] = [[1.0, 2.0], [3.0, 3.0]]
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_USER)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3:
            if s.e < 0.5:
              assert_almost_equal(s.w, 1.6) # hand calcs
              tally[0] += score
            else:
              assert_almost_equal(s.w, 0.4) # hand calcs
              tally[1] += score
        else:
            if s.e < 0.5:
              assert_almost_equal(s.w, 2.4) # hand calcs
              tally[2] += score
            else:
              assert_almost_equal(s.w, 0.8) # hand calcs
              tally[3] += score

    expected_tally = [0.25, 0.5, 0.125, 0.125] # hand calcs
    for a, b in zip(tally, expected_tally):
       assert(abs(a-b)/b < 0.25)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias_spatial():
    """This test tests a user-specified biasing scheme for which the only 1
    bias group is supplied for a source distribution containing two energy 
    groups. This bias group is applied to both energy groups. In this test,
    the user-supplied bias distribution that was choosen, correspondes to 
    uniform sampling, so that results can be checked against Case 1 in the
    theory manual.
    """
    seed(1953)
    m = Mesh(structured=True, 
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    m.bias = IMeshTag(1, float)
    m.bias[:] = [1, 1]
    e_bounds = np.array([0, 0.5, 1.0])
    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_USER)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4)) # number of phase space groups
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3.0:
            assert_almost_equal(s.w, 0.7) # hand calcs
        else:
            assert_almost_equal(s.w, 2.8) # hand calcs

        spatial_tally[int(s.x*num_divs/3.5), 
                      int(s.y*num_divs/1.0), 
                      int(s.z*num_divs/1.0)]  += score

        if s.x < 3 and s.e < 0.5:
            e_tally[0] += score
        elif s.x < 3 and s.e > 0.5:
            e_tally[1] += score
        if s.x > 3 and s.e < 0.5:
            e_tally[2] += score
        if s.x > 3 and s.e > 0.5:
            e_tally[3] += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(spatial_tally, i)[j,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

    expected_e_tally = [4./7, 2./7, 3./28, 1./28] # hand calcs
    for i in range(4):
        assert(abs(e_tally[i] - expected_e_tally[i])
               /expected_e_tally[i] < 0.1)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_1():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length of 1.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats = None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = IMeshTag(4, float)
    m.src[:] = np.empty(shape=(2, 4), dtype=float)
    m.src[:] = [[0.05, 0.05, 0.10, 0.10],
                [0.15, 0.15, 0.20, 0.20]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = IMeshTag(1, float)
    m.bias[:] = [[0.4], [0.6]]

    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.c/10 - 1
        cid = s.c%10 - 1
        eid = 0 if s.e < 0.5 else 1;
        # check the cell_number
        if s.x < 0.5:
            assert(s.c in [11, 12])
        if s.x > 0.5:
            assert(s.c in [21, 22])
        # check the weight of each subvoxel
        if vid == 0:
            assert(abs(s.w - 0.746) / 0.746 < 0.05)
        if vid == 1:
            assert(abs(s.w - 1.163) / 1.163 < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    # exp_tally calculated by hand
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.067, 0.067],
                     [0.133, 0.133]],
                    [[0.129, 0.129],
                     [0.171, 0.171]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]) / exp_tally[v, c, e] < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_max_num_cells_num_e_groups():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length
    of max_num_cells*num_e_group.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats = None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = IMeshTag(4, float)
    m.src[:] = np.empty(shape=(2,4), dtype=float)
    m.src[:] = [[0.125, 0.125, 0.125, 0.125],
                [0.125, 0.125, 0.125, 0.125]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = IMeshTag(4, float)
    m.bias[:] = [[0.125, 0.125, 0.1, 0.15], [0.1, 0.1, 0.15, 0.15]]

    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt[:] = [[[1.0, 1.0], [1.25, 0.83]], [[1.25, 1.25], [0.83, 0.83]]]
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.c/10 - 1
        cid = s.c%10 - 1
        eid = 0 if s.e < 0.5 else 1;
        # check the cell_number
        if s.x < 0.5:
            assert(s.c in [11, 12])
        if s.x > 0.5:
            assert(s.c in [21, 22])
        # check the weight of each subvoxel
        assert(abs(s.w - exp_wgt[vid, cid, eid]) / exp_wgt[vid, cid, eid] < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.125, 0.125], [0.100, 0.150]],
                    [[0.100, 0.100], [0.150, 0.150]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]) / exp_tally[v, c, e] < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_e_groups():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length
    of energy groups.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats = None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = IMeshTag(4, float)
    m.src[:] = np.empty(shape=(2,4), dtype=float)
    m.src[:] = [[0.05, 0.05, 0.10, 0.10],
                [0.15, 0.15, 0.20, 0.20]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = IMeshTag(2, float)
    m.bias[:] = [[0.1, 0.3], [0.2, 0.4]]

    filename = "sampling_mesh.h5m"
    m.mesh.save(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.c/10 - 1
        cid = s.c%10 - 1
        eid = 0 if s.e < 0.5 else 1;
        # check the cell_number
        if s.x < 0.5:
            assert(s.c in [11, 12])
        if s.x > 0.5:
            assert(s.c in [21, 22])
        # check the weight of each subvoxel
        if vid == 0 and eid == 0:
            assert(abs(s.w - 1.5)/1.5 < 0.05)
        if vid == 0 and eid == 1:
            assert(abs(s.w - 0.5)/0.5 < 0.05)
        if vid == 1 and eid == 0:
            assert(abs(s.w - 1.75)/1.75 < 0.05)
        if vid == 1 and eid == 1:
            assert(abs(s.w - 0.875)/0.875 < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.0333, 0.1000],
                     [0.0667, 0.2000]],
                    [[0.0857, 0.1714],
                     [0.1143, 0.2286]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]) / exp_tally[v, c, e] < 0.05)

def test_alias_table():
    """This tests that the AliasTable class produces samples in the ratios
    consistant with the supplied PDF.
    """
    seed(1953)
    pdf = np.array([0.1, 0.2, 0.7])
    at = AliasTable(pdf)
    num_samples = 50000
    score = 1.0/num_samples
    tally = np.zeros(shape=(3))

    for i in range(num_samples):
        s = at.sample_pdf(uniform(0, 1), uniform(0,1))
        tally[s] += score

    for i in range(0, 3):
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

def _creat_mesh_via_num_ve(num_ve):
    """
    This function used to creat mesh from number of voxels
    ----------
    num_ve : int. Number of voxels

    return : mesh. MOAB mesh.
    """
    x_coords = [v*1.0/(num_ve) for v in range(num_ve)] + [1.0]
    mesh = Mesh(structured=True,
             structured_coords=[x_coords, [0, 1], [0, 1]],
             mats = None)
    return mesh

def source_sampling_test_template(constrcut_paras=None, exp_answers=None):
    """
    This function serve as a template for all source_sampling test cases.
    It constrcut Sampler from construct_paras.
    And then perform a standardized sampling and tally,
    Finally, it compares tallied results with exp_answers.

    Assumptions:
        * Construct_paras and exp_answers are two maps.
        * Use unit cube for all the meshes.
        * Use structured meshes for all the tests.
        * filename will always be: "sampling_mesh.h5m"
        * Energy have only two options:
            - [0.0, 1.0]
            - [0.0, 0.5, 1.0]
        * Voxel number of meshes could be:
            - 1 voxel 1 subvoxel -> Single voxel single subvoxel
            - 1 voxel 2 subvoxel -> Single voxel multiple subvoxel
            - 2 voxel 2 subvoxel -> Multiple voxel multiple subvoxel (subvoxel == voxel)
            - 2 voxel 4 subvoxel -> Multiple voxel multiple subvoxel (subvoxel != voxel)
    
    Under these assumptions:
        * Mesh could be derived from cell_fracs
        * e_bounds could be derived from src_tag
        * required construct_paras contain:
            - mode
            - cell_fracs
            - src_tag
            - bias_tag (optional for USER mode)

    Check items:
        * position for each particle, and distribution
        * energy for each particle, and distribution
        * weight for each particle, and distribution
        * cell_number for each particle, and distribution
    """
    sub_mode_r2s = (0, 1, 2)
    sub_mode_subvoxel = (3, 4, 5)
    avail_mode = (0, 1, 2, 3, 4, 5)
    # input check
    if construct_paras["mode"] == None:
        raise ValueError("mode must be given")
    elif construct_paras["mode"] not in avail_mode:
        raise ValueError("mode must be in (0, 1, 2, 3, 4, 5)")
    if construct_paras["cell_fracs"] == None:
        raise ValueError("cell_fracs must be given")
    if construct_paras["src_tag"] == None:
        raise ValueError("src_tag must be given")
    # set initial value for input parameters
    mode = construct_paras["mode"]
    cell_fracs_list = construct_paras["cell_fracs"]
    src_tag = construct_paras["src_tag"]
    bias_tag = construct_paras["bias_tag"]
    # set up cell_fracs
    cell_fracs = np.zeros(len(cell_fracs_list),
                          dtype=[('idx', np.int64),
                                 ('cell', np.int64),
                                 ('vol_frac', np.float64),
                                 ('rel_error', np.float64)])
    cell_fracs[:] = cell_fracs_list
    # set up mesh
    if mode in sub_mode_r2s:
        # DEFAULT r2s
        num_ve = len(cell_fracs)
        m = _creat_mesh_via_num_ve(num_ve)
    else:
        # SUBVOXEL r2s
        num_sve = len(cell_fracs)
        num_ve = len(set(cell_fracs['idx']))
        m = _creat_mesh_via_num_ve(num_ve)
    # set up e_bounds
    num_e_groups = len(src_tag[0])




