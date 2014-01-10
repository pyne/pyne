"""alara module tests"""
import os
import nose

from nose.tools import assert_almost_equal
from nose.tools import assert_equal, assert_true, with_setup
from numpy.testing import assert_array_equal
import numpy as np
import tables as tb

# mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest
    HAVE_PYTAPS = False
    pass

from pyne.mesh import Mesh, StatMesh, MeshError
from pyne.alara import flux_mesh_to_fluxin, photon_source_to_hdf5, \
    photon_source_hdf5_to_mesh

thisdir = os.path.dirname(__file__)

def test_write_fluxin_single():
    """This function tests the flux_mesh_to_fluxin function for a single energy
    group case.
    """
    output_name = "fluxin.out"
    forward_fluxin = os.path.join(thisdir, 
                     "files_test_alara/fluxin_single_forward.txt")
    output = os.path.join(os.getcwd(), output_name)

    flux_mesh = Mesh(structured=True, structured_coords=[[0,1,2],[0,1,2],[0,1]])
    tag_flux = flux_mesh.mesh.createTag("flux", 1, float)
    flux_data = [1, 2, 3, 4]
    ves = flux_mesh.structured_iterate_hex("xyz")
    for i, ve in enumerate(ves):
        tag_flux[ve] = flux_data[i]

    # test forward writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, False)

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
    output_name = "fluxin.out"
    forward_fluxin = os.path.join(thisdir, 
                     "files_test_alara/fluxin_multiple_forward.txt")
    reverse_fluxin = os.path.join(thisdir, 
                     "files_test_alara/fluxin_multiple_reverse.txt")
    output = os.path.join(os.getcwd(), output_name)

    flux_mesh = Mesh(structured=True, structured_coords=[[0,1,2],[0,1],[0,1]])
    tag_flux = flux_mesh.mesh.createTag("flux", 7, float)
    flux_data = [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14]]
    ves = flux_mesh.structured_iterate_hex("xyz")
    for i, ve in enumerate(ves):
        tag_flux[ve] = flux_data[i]

    # test forward writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, False)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    if os.path.isfile(output):
        os.remove(output)

    # test reverse writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, True)

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
    
    with tb.openFile(filename + '.h5') as h5f:
        obs = h5f.root.data[:]
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        count = 0
        for i, row in enumerate(obs):
            ls = lines[i].strip().split('\t')
            if ls[0]  != 'TOTAL' and lines[i-1].strip().split('\t')[0] == 'TOTAL':
                if i != 0:
                    count += 1

            assert_equal(count, row['ve_idx'])
            assert_equal(ls[0].strip(), row['nuc'])
            assert_equal(ls[1].strip(), row['time'])
            assert_array_equal(np.array(ls[2:], dtype=np.float64), row['phtn_src'])

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')

def test_photon_source_hdf5_to_mesh():
    """Tests the function photon source_h5_to_mesh"""

    filename = os.path.join(thisdir, "files_test_alara", "phtn_src") 
    photon_source_to_hdf5(filename, chunkshape=(10,))
    assert_true(os.path.exists(filename + '.h5'))

    mesh = Mesh(structured=True, 
                structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])

    tags = {('h-1', 'shutdown') : 'tag1', ('TOTAL', '1 h') : 'tag2'}
    photon_source_hdf5_to_mesh(mesh, filename + '.h5', tags)

    # create lists of lists of expected results
    tag1_answers = [[1] + [0] * 41, [2] + [0] * 41, 
                    [3] + [0] * 41, [4] + [0] * 41] 
    tag2_answers = [[5] + [0] * 41, [6] + [0] * 41, 
                    [7] + [0] * 41, [8] + [0] * 41] 

    ves = list(mesh.structured_iterate_hex("xyz"))
    for i, ve in enumerate(ves):
	assert_array_equal(mesh.mesh.getTagHandle("tag1")[ve], tag1_answers[i])
        assert_array_equal(mesh.mesh.getTagHandle("tag2")[ve], tag2_answers[i])

    if os.path.isfile(filename + '.h5'):
        os.remove(filename + '.h5')
