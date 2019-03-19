"""alara module tests"""
import os
import nose
import subprocess

from nose.tools import assert_almost_equal
from nose.tools import assert_equal, assert_true, with_setup
from nose.plugins.skip import SkipTest
import numpy as np
from numpy.testing import assert_array_equal
import tables as tb
import warnings
import filecmp

# mesh specific imports
from pyne.mesh import HAVE_PYMOAB
from pyne.mesh import Mesh, StatMesh, MeshError
from pyne.alara import mesh_to_fluxin, photon_source_to_hdf5, \
    photon_source_hdf5_to_mesh, mesh_to_geom, num_density_to_mesh, \
    irradiation_blocks, record_to_geom, phtn_src_energy_bounds
from pyne.material import Material
from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)


thisdir = os.path.dirname(__file__)


def test_write_fluxin_single():
    """This function tests the flux_mesh_to_fluxin function for a single energy
    group case.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    output_name = "fluxin.out"
    forward_fluxin = os.path.join(thisdir, "files_test_alara",
                                  "fluxin_single_forward.txt")
    output = os.path.join(os.getcwd(), output_name)

    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    tag_flux = flux_mesh.tag(name="flux", size=1, dtype=float)
    flux_data = [1, 2, 3, 4]
    ves = flux_mesh.structured_iterate_hex("xyz")
    for i, ve in enumerate(ves):
        flux_mesh.flux[i] = flux_data[i]

    # test forward writting
    mesh_to_fluxin(flux_mesh, "flux", output_name, False)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)


def test_write_fluxin_multiple():
    """This function tests the flux_mesh_to_fluxin function for a multiple
    energy group case.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    output_name = "fluxin.out"
    forward_fluxin = os.path.join(thisdir, "files_test_alara",
                                  "fluxin_multiple_forward.txt")
    reverse_fluxin = os.path.join(thisdir, "files_test_alara",
                                  "fluxin_multiple_reverse.txt")
    output = os.path.join(os.getcwd(), output_name)

    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1], [0, 1]])
    flux_data = [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14]]
    flux_mesh.tag("flux", flux_data, 'nat_mesh', size=7, dtype=float)

    # test forward writting
    mesh_to_fluxin(flux_mesh, "flux", output_name, False)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)

    # test reverse writting
    mesh_to_fluxin(flux_mesh, "flux", output_name, True)

    with open(output) as f:
        written = f.readlines()

    with open(reverse_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)


def test_write_fluxin_multiple_subvoxel():
    """This function tests the flux_mesh_to_fluxin function for a multiple
    energy group case under sub-voxel r2s.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    output_name = "fluxin_subvoxel.out"
    forward_fluxin = os.path.join(thisdir, "files_test_alara",
                                  "fluxin_multiple_forward_subvoxel.txt")
    reverse_fluxin = os.path.join(thisdir, "files_test_alara",
                                  "fluxin_multiple_reverse_subvoxel.txt")
    output = os.path.join(os.getcwd(), output_name)

    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    flux_data = [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14],
                 [15, 16, 17, 18, 19, 20, 21],
                 [22, 23, 24, 25, 26, 27, 28]]
    flux_mesh.tag("flux", flux_data, 'nat_mesh', size=7, dtype=float)

    cell_fracs = np.zeros(6, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])

    cell_fracs[:] = [(0, 11, 1, 0.0), (1, 11, 0.5, 0.0), (1, 12, 0.5, 0.0),
                     (2, 11, 0.5, 0.0), (2, 13, 0.5, 0.0), (3, 13, 1, 0.0)]

    cell_mats = {11: Material({'H': 1.0}, density=1.0),
                 12: Material({'He': 1.0}, density=1.0),
                 13: Material({}, density=0.0, metadata={'name': 'void'})}

    # test forward writting
    mesh_to_fluxin(flux_mesh, "flux", output_name, False, True, cell_fracs,
                   cell_mats)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)

    # test reverse writting
    mesh_to_fluxin(flux_mesh, "flux", output_name, True, True, cell_fracs,
                   cell_mats)

    with open(output) as f:
        written = f.readlines()

    with open(reverse_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)


def test_photon_source_to_hdf5():
    """Tests the function photon_source_to_hdf5.
    """
    filename = os.path.join(thisdir, "files_test_alara", "phtn_src")
    photon_source_to_hdf5(filename, chunkshape=(10,))
    assert_true(os.path.exists(filename + '.h5'))

    with tb.open_file(filename + '.h5') as h5f:
        obs = h5f.root.data[:]

    with open(filename, 'r') as f:
        lines = f.readlines()
        count = 0
        old = ""
        for i, row in enumerate(obs):
            ls = lines[i].strip().split('\t')
            if ls[0] != 'TOTAL' and old == 'TOTAL':
                count += 1

            assert_equal(count, row['idx'])
            assert_equal(ls[0].strip(), row['nuc'].decode())
            assert_equal(ls[1].strip(), row['time'].decode())
            assert_array_equal(np.array(ls[2:], dtype=np.float64),
                               row['phtn_src'])
            old = ls[0]

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')


def test_photon_source_hdf5_to_mesh():
    """Tests the function photon source_h5_to_mesh."""

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara", "phtn_src")
    photon_source_to_hdf5(filename, chunkshape=(10,))
    assert_true(os.path.exists(filename + '.h5'))

    mesh = Mesh(structured=True,
                structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])

    tags = {('1001', 'shutdown'): 'tag1', ('TOTAL', '1.0 h'): 'tag2'}
    photon_source_hdf5_to_mesh(mesh, filename + '.h5', tags)

    # create lists of lists of expected results
    tag1_answers = [[1] + [0] * 41, [2] + [0] * 41,
                    [3] + [0] * 41, [4] + [0] * 41]
    tag2_answers = [[5] + [0] * 41, [6] + [0] * 41,
                    [7] + [0] * 41, [8] + [0] * 41]

    ves = list(mesh.structured_iterate_hex("xyz"))
    for i, ve in enumerate(ves):
        assert_array_equal(mesh.tag1[ve], tag1_answers[i])
        assert_array_equal(mesh.tag2[ve], tag2_answers[i])

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')


def test_photon_source_hdf5_to_mesh_subvoxel():
    """Tests the function photon source_h5_to_mesh
    under sub-voxel r2s condition."""

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara", "phtn_src")
    photon_source_to_hdf5(filename, chunkshape=(10,))
    assert_true(os.path.exists(filename + '.h5'))
    sub_voxel = True
    mesh = Mesh(structured=True,
                structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    cell_fracs = np.zeros(6, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])

    cell_fracs[:] = [(0, 11, 1.0, 0.0), (1, 11, 0.5, 0.0), (1, 12, 0.5, 0.0),
                     (2, 11, 0.5, 0.0), (2, 13, 0.5, 0.0), (3, 13, 1.0, 0.0)]

    cell_mats = {11: Material({'H': 1.0}, density=1.0),
                 12: Material({'He': 1.0}, density=1.0),
                 13: Material({}, density=0.0, metadata={'name': 'void'})}
    mesh.tag_cell_fracs(cell_fracs)
    tags = {('1001', 'shutdown'): 'tag1', ('TOTAL', '1 h'): 'tag2'}

    photon_source_hdf5_to_mesh(mesh, filename + '.h5', tags,
                               sub_voxel=sub_voxel, cell_mats=cell_mats)

    # create lists of lists of expected results
    tag1_answers = [[1.0] + [0.0] * 41 + [0.0] * 42,
                    [2.0] + [0.0] * 41 + [3.0] + [0.0] * 41,
                    [4.0] + [0.0] * 41 + [0.0] * 42,
                    [0.0] * 42 * 2]
    tag2_answers = [[5.0] + [0.0] * 41 + [0.0] * 42,
                    [6.0] + [0.0] * 41 + [7.0] + [0.0] * 41,
                    [8.0] + [0.0] * 41 + [0.0] * 42,
                    [0.0] * 42 * 2]

    for i, _, ve in mesh:
        assert_array_equal(mesh.tag1[ve], tag1_answers[i])
        assert_array_equal(mesh.tag2[ve], tag2_answers[i])

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')


def test_photon_source_hdf5_to_mesh_subvoxel_size1():
    """Tests the function photon source_h5_to_mesh
    under sub-voxel r2s condition."""

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara", "phtn_src")
    photon_source_to_hdf5(filename, chunkshape=(10,))
    assert_true(os.path.exists(filename + '.h5'))
    sub_voxel = True
    mesh = Mesh(structured=True,
                structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])

    cell_fracs[:] = [(0, 11, 1.0, 0.0), (1, 12, 1.0, 0.0),
                     (2, 13, 1.0, 0.0), (3, 14, 1.0, 0.0)]

    cell_mats = {11: Material({'H': 1.0}, density=1.0),
                 12: Material({'He': 1.0}, density=1.0),
                 13: Material({'He': 1.0}, density=1.0),
                 14: Material({}, density=0.0, metadata={'name': 'void'})}
    mesh.tag_cell_fracs(cell_fracs)
    tags = {('1001', 'shutdown'): 'tag1', ('TOTAL', '1 h'): 'tag2'}

    photon_source_hdf5_to_mesh(mesh, filename + '.h5', tags,
                               sub_voxel=sub_voxel, cell_mats=cell_mats)

    # create lists of lists of expected results
    tag1_answers = [[1.0] + [0.0] * 41,
                    [2.0] + [0.0] * 41,
                    [3.0] + [0.0] * 41,
                    [0.0] * 42]
    tag2_answers = [[5.0] + [0.0] * 41,
                    [6.0] + [0.0] * 41,
                    [7.0] + [0.0] * 41,
                    [0.0] * 42]

    for i, _, ve in mesh:
        assert_array_equal(mesh.tag1[ve], tag1_answers[i])
        assert_array_equal(mesh.tag2[ve], tag2_answers[i])

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')


def test_record_to_geom():

    if not HAVE_PYMOAB:
        raise SkipTest

    expected_geom = os.path.join(thisdir, "files_test_alara",
                                 "alara_record_geom.txt")
    expected_matlib = os.path.join(thisdir, "files_test_alara",
                                   "alara_record_matlib.txt")
    geom = os.path.join(os.getcwd(), "alara_record_geom")
    matlib = os.path.join(os.getcwd(), "alara_record_matlib")
    cell_fracs = np.zeros(11, dtype=[('idx', np.int64),
                                     ('cell', np.int64),
                                     ('vol_frac', np.float64),
                                     ('rel_error', np.float64)])

    cell_mats = {11: Material({'H1': 1.0, 'K39': 1.0}, density=1.1,
                              metadata={'name': 'fake_mat'}),
                 12: Material({'H1': 0.1, 'O16': 1.0}, density=1.2,
                              metadata={'name': 'water'}),
                 13: Material({'He4': 42.0}, density=1.3,
                              metadata={'name': 'helium'}),
                 14: Material({}, density=0.0, metadata={'name': 'void'}),
                 15: Material({}, density=0.0, metadata={'name': 'void'}),
                 16: Material({'H1': 1.0, 'K39': 1.0}, density=1.1,
                              metadata={'name': 'fake_mat'})}

    cell_fracs[:] = [(0, 11, 0.55, 0.0), (0, 12, 0.45, 0.0), (1, 11, 0.2, 0.0),
                     (1, 12, 0.3, 0.0), (1, 13, 0.5, 0.0), (2, 11, 0.15, 0.0),
                     (2, 14, 0.01, 0.0), (2, 15, 0.04, 0.0), (2, 16, 0.8, 0.0),
                     (3, 11, 0.55, 0.0), (3, 12, 0.45, 0.0)]

    m = Mesh(structured_coords=[[-1, 0, 1], [-1, 0, 1], [0, 1]],
             structured=True, mats=None)

    record_to_geom(m, cell_fracs, cell_mats, geom, matlib)

    assert(filecmp.cmp(geom, expected_geom))
    if os.path.isfile(geom):
        os.remove(geom)

    assert(filecmp.cmp(matlib, expected_matlib))
    if os.path.isfile(matlib):
        os.remove(matlib)


def test_record_to_geom_subvoxel():

    if not HAVE_PYMOAB:
        raise SkipTest

    expected_geom = os.path.join(thisdir, "files_test_alara",
                                 "alara_record_geom_subvoxel.txt")
    expected_matlib = os.path.join(thisdir, "files_test_alara",
                                   "alara_record_matlib_subvoxel.txt")
    geom = os.path.join(os.getcwd(), "alara_record_geom")
    matlib = os.path.join(os.getcwd(), "alara_record_matlib")
    cell_fracs = np.zeros(11, dtype=[('idx', np.int64),
                                     ('cell', np.int64),
                                     ('vol_frac', np.float64),
                                     ('rel_error', np.float64)])

    cell_mats = {11: Material({'H1': 1.0, 'K39': 1.0}, density=1.1,
                              metadata={'name': 'fake_mat'}),
                 12: Material({'H1': 0.1, 'O16': 1.0}, density=1.2,
                              metadata={'name': 'water'}),
                 13: Material({'He4': 42.0}, density=1.3,
                              metadata={'name': 'helium'}),
                 14: Material({}, density=0.0, metadata={'name': 'void'}),
                 15: Material({}, density=0.0, metadata={'name': 'void'}),
                 16: Material({'H1': 1.0, 'K39': 1.0}, density=1.1,
                              metadata={'name': 'fake_mat'})}

    cell_fracs[:] = [(0, 11, 0.55, 0.0), (0, 12, 0.45, 0.0), (1, 11, 0.2, 0.0),
                     (1, 12, 0.3, 0.0), (1, 13, 0.5, 0.0), (2, 11, 0.15, 0.0),
                     (2, 14, 0.01, 0.0), (2, 15, 0.04, 0.0), (2, 16, 0.8, 0.0),
                     (3, 11, 0.55, 0.0), (3, 12, 0.45, 0.0)]

    m = Mesh(structured_coords=[[-1, 0, 1], [-1, 0, 1], [0, 1]],
             structured=True, mats=None)

    record_to_geom(m, cell_fracs, cell_mats, geom, matlib, sub_voxel=True)

    assert(filecmp.cmp(geom, expected_geom))
    if os.path.isfile(geom):
        os.remove(geom)

    assert(filecmp.cmp(matlib, expected_matlib))
    if os.path.isfile(matlib):
        os.remove(matlib)


def test_mesh_to_geom():

    if not HAVE_PYMOAB:
        raise SkipTest

    expected_geom = os.path.join(thisdir, "files_test_alara", "alara_geom.txt")
    expected_matlib = os.path.join(thisdir, "files_test_alara",
                                   "alara_matlib.txt")
    geom = os.path.join(os.getcwd(), "alara_geom")
    matlib = os.path.join(os.getcwd(), "alara_matlib")

    mats = {
        0: Material({'H1': 1.0, 'K39': 1.0}, density=1.1),
        1: Material({'H1': 0.1, 'O16': 1.0}, density=1.2),
        2: Material({'He4': 42.0}, density=1.3),
        3: Material({'Tm171': 171.0}, density=1.4),
    }
    m = Mesh(structured_coords=[[-1, 0, 1], [-1, 0, 1], [0, 1]], structured=True,
             mats=mats)
    mesh_to_geom(m, geom, matlib)

    with open(expected_geom) as f:
        written = f.readlines()

    with open(geom) as f:
        expected = f.readlines()

    assert_equal(written, expected)

    if os.path.isfile(geom):
        os.remove(geom)

    with open(expected_matlib) as f:
        written = f.readlines()

    with open(matlib) as f:
        expected = f.readlines()

    assert_equal(written, expected)

    if os.path.isfile(matlib):
        os.remove(matlib)


def test_num_den_to_mesh_shutdown():

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara",
                            "num_density_output.txt")
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1, 2]])
    with open(filename) as f:
        lines = f.readlines()
    num_density_to_mesh(lines, 'shutdown', m)

    # expected composition results:
    exp_comp_0 = {10010000: 5.3390e+19,
                  10020000: 3.0571e+17,
                  10030000: 1.2082e+12,
                  20030000: 7.4323e+09,
                  20040000: 7.1632e+02}
    exp_comp_1 = {10010000: 4.1240e+13,
                  10020000: 4.7443e+11,
                  10030000: 2.6627e+13,
                  20030000: 8.3547e+10,
                  20040000: 2.6877e+19}

    # actual composition results
    act_comp_0 = m.mats[0].to_atom_frac()
    act_comp_1 = m.mats[1].to_atom_frac()

    assert_equal(len(exp_comp_0), len(act_comp_0))
    for key, value in exp_comp_0.iteritems():
        assert_almost_equal(value/act_comp_0[key], 1.0, 15)

    assert_equal(len(exp_comp_1), len(act_comp_1))
    for key, value in exp_comp_1.iteritems():
        assert_almost_equal(value/act_comp_1[key], 1.0, 15)

    # compare densities
    exp_density_0 = 8.96715E-05
    exp_density_1 = 1.785214E-04

    assert_almost_equal(exp_density_0, m.mats[0].density)
    assert_almost_equal(exp_density_1, m.mats[1].density)


def test_num_den_to_mesh_stdout():

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara",
                            "num_density_output.txt")
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1, 2]])

    p = subprocess.Popen(["cat", filename], stdout=subprocess.PIPE)
    lines, err = p.communicate()

    num_density_to_mesh(lines.split('\n'), 'shutdown', m)

    # expected composition results:
    exp_comp_0 = {10010000: 5.3390e+19,
                  10020000: 3.0571e+17,
                  10030000: 1.2082e+12,
                  20030000: 7.4323e+09,
                  20040000: 7.1632e+02}
    exp_comp_1 = {10010000: 4.1240e+13,
                  10020000: 4.7443e+11,
                  10030000: 2.6627e+13,
                  20030000: 8.3547e+10,
                  20040000: 2.6877e+19}

    # actual composition results
    act_comp_0 = m.mats[0].to_atom_frac()
    act_comp_1 = m.mats[1].to_atom_frac()

    assert_equal(len(exp_comp_0), len(act_comp_0))
    for key, value in exp_comp_0.iteritems():
        assert_almost_equal(value/act_comp_0[key], 1.0, 15)

    assert_equal(len(exp_comp_1), len(act_comp_1))
    for key, value in exp_comp_1.iteritems():
        assert_almost_equal(value/act_comp_1[key], 1.0, 15)

    # compare densities
    exp_density_0 = 8.96715E-05
    exp_density_1 = 1.785214E-04

    assert_almost_equal(exp_density_0, m.mats[0].density)
    assert_almost_equal(exp_density_1, m.mats[1].density)


def test_num_den_to_mesh_1_y():

    if not HAVE_PYMOAB:
        raise SkipTest

    filename = os.path.join(thisdir, "files_test_alara",
                            "num_density_output.txt")
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1, 2]])
    num_density_to_mesh(filename, '1 y', m)

    # expected results:
    exp_comp_0 = {10010000: 5.3390e+19,
                  10020000: 3.0571e+17,
                  10030000: 1.1424e+12,
                  20030000: 7.3260e+10,
                  20040000: 7.1632e+02}
    exp_comp_1 = {10010000: 4.1240e+13,
                  10020000: 4.7443e+11,
                  10030000: 2.5176e+13,
                  20030000: 1.5343e+12,
                  20040000: 2.6877e+19}

    # actual results
    act_comp_0 = m.mats[0].to_atom_frac()
    act_comp_1 = m.mats[1].to_atom_frac()

    assert_equal(len(exp_comp_0), len(act_comp_0))
    for key, value in exp_comp_0.iteritems():
        assert_almost_equal(value/act_comp_0[key], 1.0, 15)

    assert_equal(len(exp_comp_1), len(act_comp_1))
    for key, value in exp_comp_1.iteritems():
        assert_almost_equal(value/act_comp_1[key], 1.0, 15)

    # compare densities
    exp_density_0 = 8.96715E-05
    exp_density_1 = 1.78521E-04
    assert_almost_equal(exp_density_0, m.mats[0].density)
    assert_almost_equal(exp_density_1, m.mats[1].density)


def test_irradiation_blocks():

    # actual results
    act = irradiation_blocks("matlib", "isolib",
                             "FEINDlib CINDER CINDER90 THERMAL",
                             ["1 h", "0.5 y"], "fluxin.out", "1 y",
                             output="number_density")

    exp = ("material_lib matlib\n"
           "element_lib isolib\n"
           "data_library FEINDlib CINDER CINDER90 THERMAL\n"
           "\n"
           "cooling\n"
           "    1 h\n"
           "    0.5 y\n"
           "end\n"
           "\n"
           "flux flux_1 fluxin.out 1.0 0 default\n"
           "schedule simple_schedule\n"
           "    1 y flux_1 pulse_once 0 s\n"
           "end\n"
           "\n"
           "pulsehistory pulse_once\n"
           "    1 0.0 s\n"
           "end\n"
           "\n"
           "output zone\n"
           "    units Ci cm3\n"
           "    number_density\n"
           "end\n"
           "\n"
           "truncation 1e-12\n"
           "impurity 5e-06 0.001\n"
           "dump_file dump_file\n")

    assert_equal(act, exp)


def test_phtn_src_energy_bounds():

    input_file = os.path.join(thisdir, "files_test_alara",
                              "alara_other_blocks.txt")
    e_bounds = phtn_src_energy_bounds(input_file)
    expected_e_bounds = [0, 1.00E4, 2.00E4, 5.00E4, 1.00E5, 2.00E5, 3.00E5,
                         4.00E5, 6.00E5, 8.00E5, 1.00E6, 1.22E6, 1.44E6, 1.66E6,
                         2.00E6, 2.50E6, 3.00E6, 4.00E6, 5.00E6, 6.50E6, 8.00E6,
                         1.00E7, 1.20E7, 1.40E7, 2.00E7]

    assert_array_equal(e_bounds, expected_e_bounds)
