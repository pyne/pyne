"""alara module tests"""
import os
import nose

from nose.tools import assert_almost_equal
from nose.tools import assert_equal
from numpy.testing import assert_array_equal

# mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest
    HAVE_PYTAPS = False
    pass

from pyne.mesh import Mesh, StatMesh, MeshError
from pyne.alara import flux_mesh_to_fluxin

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
    ves = flux_mesh.structured_iterate_hex("zyx")
    for i, ve in enumerate(ves):
        tag_flux[ve] = flux_data[i]

    # test forward writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, False)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
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
    ves = flux_mesh.structured_iterate_hex("zyx")
    for i, ve in enumerate(ves):
        tag_flux[ve] = flux_data[i]

    # test forward writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, False)

    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
    os.remove(output)

    # test reverse writting
    flux_mesh_to_fluxin(flux_mesh, "flux", output_name, True)

    with open(output) as f:
        written = f.readlines()

    with open(reverse_fluxin) as f:
        expected = f.readlines()

    assert_equal(written, expected)
