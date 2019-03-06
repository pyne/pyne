"""
Tests for PyNE variance_reduction module.
"""
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_almost_equal
from nose.plugins.skip import SkipTest
import warnings

from pyne.mesh import Mesh, NativeMeshTag, MeshError, HAVE_PYMOAB
if not HAVE_PYMOAB:
    raise SkipTest
try:
    from pyne import mcnp
except ImportError:
    raise SkipTest
from pyne.variancereduction import cadis, magic
from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)


def test_cadis_single_e():
    """Test single energy group cadis case"""
    adj_flux_tag = "adj_flux"
    q_tag = "q"
    ww_tag = "ww"
    q_bias_tag = "q_bias"

    # create meshes
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    adj_flux_mesh = Mesh(structured=True, structured_coords=coords)
    q_mesh = Mesh(structured=True, structured_coords=coords)
    ww_mesh = Mesh(structured=True, structured_coords=coords)
    q_bias_mesh = Mesh(structured=True, structured_coords=coords)

    # add tags to meshes
    adj_flux_mesh.adj_flux = NativeMeshTag(1, float)
    q_mesh.q = NativeMeshTag(1, float)

    # create data for input meshes
    adj_flux_data = [1.1, 1.3, 1.5, 1.7]
    q_data = [2.9, 2.6, 2.4, 2.2]

    # data data to mesh
    adj_flux_mesh.adj_flux[:] = adj_flux_data
    q_mesh.q[:] = q_data

    # run CADIS
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)

    # checkout output meshes
    expected_ww = [0.3995338, 0.33806706, 0.29299145, 0.258521908]
    expected_q_bias = [0.04652859, 0.04929988, 0.052508751, 0.0545507858]
    ww_mesh.ww = NativeMeshTag(1, float)
    q_bias_mesh.q_bias = NativeMeshTag(1, float)

    assert_array_almost_equal(ww_mesh.ww[:], expected_ww[:])

    assert_array_almost_equal(q_bias_mesh.q_bias[:], expected_q_bias[:])


def test_cadis_multiple_e():
    """Test multiple energy group CADIS case"""

    adj_flux_tag = "adj_flux"
    q_tag = "q"
    ww_tag = "ww"
    q_bias_tag = "q_bias"

    # create meshes
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    adj_flux_mesh = Mesh(structured=True, structured_coords=coords)
    q_mesh = Mesh(structured=True, structured_coords=coords)
    ww_mesh = Mesh(structured=True, structured_coords=coords)
    q_bias_mesh = Mesh(structured=True, structured_coords=coords)

    # add tags to input meshes
    adj_flux_mesh.adj_flux = NativeMeshTag(2, float)
    q_mesh.q = NativeMeshTag(2, float)

    # create data for input meshes
    adj_flux_data = [[1.1, 1.2], [1.3, 1.4], [0.0, 1.6], [1.7, 1.9]]
    q_data = [[2.9, 2.8], [2.6, 2.5], [2.4, 2.2], [2.9, 0.0]]

    # data data to mesh
    adj_flux_mesh.adj_flux[:] = adj_flux_data
    q_mesh.q[:] = q_data

    # run cadis
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)

    # expected results
    expected_q_bias = [[0.0306200806, 0.0322518718], [0.0324438472, 0.0335956998],
                       [0.0, 0.0337876752], [0.0473219428, 0.0]]
    expected_ww = [[0.3208302538, 0.2940943993], [0.2714717532, 0.2520809137],
                   [0.0, 0.2205707995], [0.2075960465, 0.1857438311]]

    ww_mesh.ww = NativeMeshTag(2, float)
    q_bias_mesh.q_bias = NativeMeshTag(2, float)

    assert_array_almost_equal(ww_mesh.ww[:], expected_ww[:])
    assert_array_almost_equal(q_bias_mesh.q_bias[:], expected_q_bias[:])


def test_magic_below_tolerance():
    """Test MAGIC case when all flux errors are below the default tolerance"""

    # create mesh
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    flux_data = [1.2, 3.3, 1.6, 1.7]
    flux_error = [0.11, 0.013, 0.14, 0.19]
    tally = Mesh(structured=True, structured_coords=coords)

    tally.particle = "neutron"
    tally.e_bounds = [0.0, 0.5, 1.0]
    tally.n_total_flux = NativeMeshTag(1, float)
    tally.n_total_flux[:] = flux_data

    tally.n_rel_error = NativeMeshTag(1, float)
    tally.n_rel_error[:] = flux_error

    magic(tally, "n_total_flux", "n_rel_error")

    expected_ww = [0.181818182, 0.5, 0.2424242, 0.2575757576]

    assert_array_almost_equal(tally.ww_x[:], expected_ww[:])


def test_magic_multi_bins():
    """Test multiple energy bins MAGIC case"""

    # create a basic meshtally
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    flux_data = [[1.2, 3.3], [1.6, 1.7], [1.5, 1.4], [2.6, 1.0]]
    flux_error = [[0.11, 0.013], [0.14, 0.19], [0.02, 0.16], [0.04, 0.09]]
    tally = Mesh(structured=True, structured_coords=coords)

    tally.particle = "neutron"
    tally.e_bounds = [0.0, 0.5, 1.0]
    tally.n_flux = NativeMeshTag(2, float)
    tally.n_flux[:] = flux_data

    tally.n_rel_error = NativeMeshTag(2, float)
    tally.n_rel_error[:] = flux_error

    tolerance = 0.15
    null_value = 0.001

    magic(tally, "n_flux", "n_rel_error",
          tolerance=tolerance, null_value=null_value)

    expected_ww = [[0.2307692308, 0.5],
                   [0.3076923077, 0.001],
                   [0.2884615385, 0.001],
                   [0.5, 0.15151515]]

    assert_array_almost_equal(tally.ww_x[:], expected_ww[:])


def test_magic_e_total():
    """Test total energy group MAGIC case"""

    # create mesh
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    flux_data = [1.2, 3.3, 1.6, 1.7]
    flux_error = [0.11, 0.013, 0.14, 0.19]
    tally = Mesh(structured=True, structured_coords=coords)

    tally.particle = "neutron"
    tally.e_bounds = [0.0, 0.5, 1.0]
    tally.n_total_flux = NativeMeshTag(1, float)
    tally.n_total_flux[:] = flux_data

    tally.n_rel_error = NativeMeshTag(1, float)
    tally.n_rel_error[:] = flux_error

    tolerance = 0.15
    null_value = 0.001

    magic(tally, "n_total_flux", "n_rel_error",
          tolerance=tolerance, null_value=null_value)

    expected_ww = [0.181818182, 0.5, 0.2424242, 0.001]

    assert_array_almost_equal(tally.ww_x[:], expected_ww[:])


def test_magic_single_e():
    """Test a single energy group MAGIC case"""

    # create mesh
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    flux_data = [1.2, 3.3, 1.6, 1.7]
    flux_error = [0.11, 0.013, 0.14, 0.19]
    tally = Mesh(structured=True, structured_coords=coords)

    tally.particle = "neutron"
    tally.e_bounds = [0.0, 1.0]
    tally.n_flux = NativeMeshTag(1, float)
    tally.n_flux[:] = flux_data

    tally.n_rel_error = NativeMeshTag(1, float)
    tally.n_rel_error[:] = flux_error

    tolerance = 0.15
    null_value = 0.001

    magic(tally, "n_flux", "n_rel_error",
          tolerance=tolerance, null_value=null_value)

    expected_ww = [0.181818182, 0.5, 0.2424242, 0.001]

    assert_array_almost_equal(tally.ww_x[:], expected_ww[:])
