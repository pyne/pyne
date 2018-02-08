import os
import warnings
import itertools

from operator import itemgetter
from nose.tools import assert_equal, with_setup, assert_almost_equal
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
from pyne.source_sampling import Sampler, AliasTable, PySourceParticle

def try_rm_file(filename):
    return lambda: os.remove(filename) if os.path.exists(filename) else None

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an analog sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 16 regions
    of phase space the particle is born into.
    """
    seed(1953)
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
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        c = src.get_src_c()
        assert_equal(s[4], 1.0) # analog: all weights must be one
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs),
              int(s[3]*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex_subvoxel_analog():
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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 1]), False)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        assert_equal(s[4], 1.0) # analog: all weights must be one
        assert_equal(cell_number, 11) # analog: the cell number
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs),
              int(s[3]*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 1]), True)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        assert_equal(s[4], 1.0) # analog: all weights must be one
        assert_equal(cell_number, 11) # analog: the cell number
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs),
              int(s[3]*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j,:,:,:]) - 0.5) < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex_multiple_subvoxel():
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

    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 1]), False)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3

    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        assert_equal(s[4], 1.0) # analog: all weights must be one
        if int(cell_number) == 11:
            tally[0] += score
        if int(cell_number) == 12:
            tally[1] += score
        if int(cell_number) == 13:
            tally[2] += score

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.158)/0.158 < 0.05)
    assert(abs(tally[2] - 0.842)/0.842 < 0.05)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex_multiple_subvoxel_uniform():
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

    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 1]), True)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3

    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        if int(cell_number) == 11:
            tally[0] += score
        if int(cell_number) == 12:
            tally[1] += score
            assert(abs(s[4] - 0.369)/0.369 < 0.05) # analog: all weights must be one
        if int(cell_number) == 13:
            tally[2] += score
            assert(abs(s[4] - 1.475)/1.475 < 0.05)

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.428) < 0.05)
    assert(abs(tally[2] - 0.572) < 0.05)

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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", np.array([0, 0.5, 1]), False)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        assert_equal(s[4], 1.0)
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs),
              int(s[3]*num_divs)] += score

    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_multiple_hex_subvoxel():
    """This test tests that particle are sampled analog from a uniform source
    defined on eight mesh volume elements in two energy groups. This is done
    using the exact same method as test_analog_multiple_hex_subvoxel.
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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 0.5, 1]), False)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        assert_equal(s[4], 1.0)
        assert_equal(cell_number, 4*int(s[0]*num_divs) + 2*int(s[1]*num_divs)
                     + int(s[2]*num_divs) + 1)
        tally[int(s[0]*num_divs), int(s[1]*num_divs), int(s[2]*num_divs),
              int(s[3]*num_divs)] += score

    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j,:,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_multiple_hex_subvoxel_uniform():
    """This test tests that particle are sampled uniformly from a uniform source
    defined on eight mesh volume elements in two energy groups. This is done
    using the exact same method ass test_analog_multiple_hex_subvoxel.
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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      np.array([0, 0.5, 1]), True)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 8
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        # check the cell_number
        assert_equal(cell_number, 4*int(s[0]*num_divs) + 2*int(s[1]*num_divs)
                     + int(s[2]*num_divs))
        # check the weight of each subvoxel
        if cell_number not in empty_cells:
            # weight for cell 1, 3, 5, 7 should be: 0.4, 0.8, 1.2, 1.6
            exp_w = (cell_number + 1) / 2 * 0.4
            out_w = s[4]
            assert(abs(out_w - exp_w)/exp_w < 0.05) # hand calculate
        # count the tally
        tally[cell_number] += score

    # check the real sample rate
    for i, item in enumerate(tally):
        if i not in empty_cells:
            assert(abs(item - 0.25)/0.25 < 0.05)

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
    m.src[:] = np.empty(shape=(2,4), dtype=float)
    m.src[:] = [[0.05, 0.05, 0.10, 0.10],
                [0.15, 0.15, 0.20, 0.20]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = IMeshTag(1, float)
    m.bias[:] = [[0.4], [0.6]]

    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      e_bounds, "bias")

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        vid = cell_number/10 - 1
        cid = cell_number%10 - 1
        eid = 0 if s[3] < 0.5 else 1;
        wgt = s[4]
        # check the cell_number
        if s[0] < 0.5:
            assert(cell_number in [11, 12])
        if s[0] > 0.5:
            assert(cell_number in [21, 22])
        # check the weight of each subvoxel
        if vid == 0:
            assert(abs(wgt - 0.746) / 0.746 < 0.05)
        if vid == 1:
            assert(abs(wgt - 1.163) / 1.163 < 0.05)
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

    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      e_bounds, "bias")

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt[:] = [[[1.0, 1.0], [1.25, 0.83]], [[1.25, 1.25], [0.83, 0.83]]]
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        vid = cell_number/10 - 1
        cid = cell_number%10 - 1
        eid = 0 if s[3] < 0.5 else 1;
        wgt = s[4]
        # check the cell_number
        if s[0] < 0.5:
            assert(cell_number in [11, 12])
        if s[0] > 0.5:
            assert(cell_number in [21, 22])
        # check the weight of each subvoxel
        assert(abs(wgt - exp_wgt[vid, cid, eid]) / exp_wgt[vid, cid, eid] < 0.05)
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

    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", "cell_number", "cell_fracs",
                      e_bounds, "bias")

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
        cell_number = src.get_src_c()
        vid = cell_number/10 - 1
        cid = cell_number%10 - 1
        eid = 0 if s[3] < 0.5 else 1;
        wgt = s[4]
        # check the cell_number
        if s[0] < 0.5:
            assert(cell_number in [11, 12])
        if s[0] > 0.5:
            assert(cell_number in [21, 22])
        # check the weight of each subvoxel
        if vid == 0 and eid == 0:
            assert(abs(wgt - 1.5)/1.5 < 0.05)
        if vid == 0 and eid == 1:
            assert(abs(wgt - 0.5)/0.5 < 0.05)
        if vid == 1 and eid == 0:
            assert(abs(wgt - 1.75)/1.75 < 0.05)
        if vid == 1 and eid == 1:
            assert(abs(wgt - 0.875)/0.875 < 0.05)
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
        src = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        s = src.get_src_xyzew()
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
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats = None)
    m.src = IMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, True)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4)) # number of phase space groups
    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        if s[0] < 3.0:
            assert_almost_equal(s[4], 0.7) # hand calcs
        else:
            assert_almost_equal(s[4], 2.8) # hand calcs

        spatial_tally[int(s[0]*num_divs/3.5),
                      int(s[1]*num_divs/1.0),
                      int(s[2]*num_divs/1.0)]  += score

        if s[0] < 3 and s[3] < 0.5:
            e_tally[0] += score
        elif s[0] < 3 and s[3] > 0.5:
            e_tally[1] += score
        if s[0] > 3 and s[3] < 0.5:
            e_tally[2] += score
        if s[0] > 3 and s[3] > 0.5:
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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, "bias")

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        if s[0] < 3:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 1.6) # hand calcs
              tally[0] += score
            else:
              assert_almost_equal(s[4], 0.4) # hand calcs
              tally[1] += score
        else:
            if s[3] < 0.5:
              assert_almost_equal(s[4], 2.4) # hand calcs
              tally[2] += score
            else:
              assert_almost_equal(s[4], 0.8) # hand calcs
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
    m.mesh.save("sampling_mesh.h5m")
    sampler = Sampler("sampling_mesh.h5m", "src", e_bounds, "bias")

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4)) # number of phase space groups
    for i in range(num_samples):
        src = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        s = src.get_src_xyzew()
        if s[0] < 3.0:
            assert_almost_equal(s[4], 0.7) # hand calcs
        else:
            assert_almost_equal(s[4], 2.8) # hand calcs

        spatial_tally[int(s[0]*num_divs/3.5),
                      int(s[1]*num_divs/1.0),
                      int(s[2]*num_divs/1.0)]  += score

        if s[0] < 3 and s[3] < 0.5:
            e_tally[0] += score
        elif s[0] < 3 and s[3] > 0.5:
            e_tally[1] += score
        if s[0] > 3 and s[3] < 0.5:
            e_tally[2] += score
        if s[0] > 3 and s[3] > 0.5:
            e_tally[3] += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(spatial_tally, i)[j,:,:])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

    expected_e_tally = [4./7, 2./7, 3./28, 1./28] # hand calcs
    for i in range(4):
        assert(abs(e_tally[i] - expected_e_tally[i])
               /expected_e_tally[i] < 0.1)

def test_alias_table():
    """This tests that the AliasTable class produces samples in the ratios
    consistent with the supplied PDF.
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

