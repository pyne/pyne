"""
Tests for PyNE variance_reduction module.
"""
import warnings

try:
    from itaps import iBase, iMesh, iMeshExtensions
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest
from nose.tools import assert_almost_equal
from numpy.testing import assert_array_almost_equal

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne.variancereduction import cadis
from pyne.mesh import Mesh
from pyne.mesh import MeshError

def test_cadis_single_e():
    """Test single energy group cadis case"""
    adj_flux_tag = "adj_flux"
    q_tag = "q"
    ww_tag = "ww"
    q_bias_tag= "q_bias"

    #create meshes
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    adj_flux_mesh = Mesh(structured=True, structured_coords=coords)
    q_mesh = Mesh(structured=True, structured_coords=coords)
    ww_mesh = Mesh(structured=True, structured_coords=coords)
    q_bias_mesh = Mesh(structured=True, structured_coords=coords)
    
    #add tags to input meshes
    tag_adj_flux = adj_flux_mesh.mesh.createTag(adj_flux_tag, 1, float)
    tag_q = q_mesh.mesh.createTag(q_tag, 1, float)

    #create data for input meshes
    adj_flux_data = [1.1, 1.3, 1.5, 1.7]
    q_data = [2.9, 2.6, 2.4, 2.2]

    #data data to mesh
    adj_flux_ves = list(adj_flux_mesh.mesh.iterate(iBase.Type.region, 
                                                   iMesh.Topology.all))
    for i, adj_flux_ve in enumerate(adj_flux_ves):
        tag_adj_flux[adj_flux_ve] = adj_flux_data[i]

    q_ves = list(q_mesh.mesh.iterate(iBase.Type.region, 
                                         iMesh.Topology.all))
    for i, q_ve in enumerate(q_ves):
        tag_q[q_ve] = q_data[i]

    #run CADIS
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)
    
    #checkout output meshes
    expected_ww = [1.4535005225, 1.3717948718, 1.287962963, 1.2397504456]
    expected_q_bias = [0.2293314162, 0.2429906542, 0.2588066139, 0.2688713156]
    
    ww_ves = list(ww_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all))
    for i, ww_ve in enumerate(ww_ves):
        assert_almost_equal(ww_mesh.mesh.getTagHandle(ww_tag)[ww_ve], 
                            expected_ww[i])

    q_bias_ves = list(q_bias_mesh.mesh.iterate(iBase.Type.region, 
                                         iMesh.Topology.all))
    for i, q_bias_ve in enumerate(q_bias_ves):
        assert_almost_equal(
            q_bias_mesh.mesh.getTagHandle(q_bias_tag)[q_bias_ve], 
            expected_q_bias[i])
            

def test_cadis_multiple_e():
    """Test multiple energy group CADIS case"""

    adj_flux_tag = "adj_flux"
    q_tag = "q"
    ww_tag = "ww"
    q_bias_tag= "q_bias"

    #create meshes
    coords = [[0, 1, 2], [-1, 3, 4], [10, 12]]
    adj_flux_mesh = Mesh(structured=True, structured_coords=coords)
    q_mesh = Mesh(structured=True, structured_coords=coords)
    ww_mesh = Mesh(structured=True, structured_coords=coords)
    q_bias_mesh = Mesh(structured=True, structured_coords=coords)
    
    #add tags to input meshes
    tag_adj_flux = adj_flux_mesh.mesh.createTag(adj_flux_tag, 2, float)
    tag_q = q_mesh.mesh.createTag(q_tag, 2, float)

    #create data for input meshes
    adj_flux_data = [[1.1, 1.2], [1.3, 1.4], [1.5, 1.6], [1.7, 1.8]]
    q_data = [[2.9, 2.8], [2.6, 2.5], [2.4, 2.2], [2.2, 2.0]]

    #data data to mesh
    adj_flux_ves = list(adj_flux_mesh.mesh.iterate(iBase.Type.region, 
                                                       iMesh.Topology.all))
    for i, adj_flux_ve in enumerate(adj_flux_ves):
        tag_adj_flux[adj_flux_ve] = adj_flux_data[i]

    q_ves = list(q_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all))
    for i, q_ve in enumerate(q_ves):
        tag_q[q_ve] = q_data[i]

    #run cadis
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)
    
    #expected results
    expected_ww = [[1.4535005225, 1.3869047619], 
                   [1.3717948718, 1.3314285714], 
                   [1.287962963, 1.3238636364], 
                   [1.2397504456, 1.2944444444]]
    expected_q_bias = [[0.2293314162, 0.2403433476], 
                       [0.2429906542, 0.2503576538],
                       [0.2588066139, 0.251788269],
                       [0.2688713156, 0.2575107296]]
    
    ww_ves = list(ww_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all))
    for i, ww_ve in enumerate(ww_ves):
        assert_array_almost_equal(ww_mesh.mesh.getTagHandle(ww_tag)[ww_ve], 
                                  expected_ww[i])

    q_bias_ves = list(q_bias_mesh.mesh.iterate(iBase.Type.region, 
                                               iMesh.Topology.all))

    for i, q_bias_ve in enumerate(q_bias_ves):
        assert_array_almost_equal(
            q_bias_mesh.mesh.getTagHandle(q_bias_tag)[q_bias_ve],
            expected_q_bias[i])
