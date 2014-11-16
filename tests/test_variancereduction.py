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

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne.variancereduction import cadis, magic

from pyne.mesh import Mesh, IMeshTag
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
    
    #add tags to meshes
    adj_flux_mesh.adj_flux = IMeshTag(1, float)
    q_mesh.q = IMeshTag(1, float)

    #create data for input meshes
    adj_flux_data = [1.1, 1.3, 1.5, 1.7]
    q_data = [2.9, 2.6, 2.4, 2.2]

    #data data to mesh
    adj_flux_mesh.adj_flux[:] = adj_flux_data
    q_mesh.q[:] = q_data

    #run CADIS
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)
    
    #checkout output meshes
    expected_ww = [0.3995338, 0.33806706, 0.29299145, 0.258521908]
    expected_q_bias = [0.04652859, 0.04929988, 0.052508751, 0.0545507858]
    ww_mesh.ww = IMeshTag(1, float)
    q_bias_mesh.q_bias = IMeshTag(1, float)
    
    assert_array_almost_equal(ww_mesh.ww[:], expected_ww[:])

    assert_array_almost_equal(q_bias_mesh.q_bias[:], expected_q_bias[:])
            

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
    adj_flux_mesh.adj_flux = IMeshTag(2, float)
    q_mesh.q = IMeshTag(2, float)

    #create data for input meshes
    adj_flux_data = [[1.1, 1.2], [1.3, 1.4], [0.0, 1.6], [1.7, 1.9]]
    q_data = [[2.9, 2.8], [2.6, 2.5], [2.4, 2.2], [2.9, 0.0]]

    #data data to mesh
    adj_flux_mesh.adj_flux[:] = adj_flux_data
    q_mesh.q[:] = q_data

    #run cadis
    cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag,
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5)
    
    #expected results
    expected_q_bias = [[0.0306200806, 0.0322518718], [0.0324438472, 0.0335956998],
                       [0.0, 0.0337876752], [0.0473219428, 0.0]]
    expected_ww = [[0.3208302538, 0.2940943993], [0.2714717532, 0.2520809137],
                   [0.0, 0.2205707995], [0.2075960465, 0.1857438311]]

    ww_mesh.ww = IMeshTag(2, float)
    q_bias_mesh.q_bias = IMeshTag(2, float)
    
    assert_array_almost_equal(ww_mesh.ww[:], expected_ww[:])
    assert_array_almost_equal(q_bias_mesh.q_bias[:], expected_q_bias[:])
    

def test_magic_multi_bins():
    # create a basic meshtally
    meshtal_file = "mcnp_meshtal_single_meshtal.txt"
    tags = {4: ["n_result", "n_rel_error",
                "n_total_result", "n_total_rel_error"]}
    meshtal_object = mcnp.Meshtal(meshtal_file, tags, meshes_have_mats=True)
    
    results = magic(meshtal_object.tally[4], "n_result", 
                    "n_rel_error", .2, .001)
                
    expected_ww_x = [   [ 0.05572753,  0.08699954,  0.07214425]
                        [ 0.08686583,  0.11118486,  0.09164014]
                        [ 0.001     ,  0.06850455,  0.05955106]
                        [ 0.12385143,  0.12588822,  0.11486921]
                        [ 0.14062005,  0.15200565,  0.15422776]
                        [ 0.08140096,  0.09324095,  0.09507091]
                        [ 0.11855448,  0.13473888,  0.13326227]
                        [ 0.17357025,  0.22207127,  0.20034488]
                        [ 0.11799773,  0.14092871,  0.10970505]
                        [ 0.0959014 ,  0.10666259,  0.11924539]
                        [ 0.12186914,  0.17309702,  0.15669681]
                        [ 0.10607864,  0.09856309,  0.09243499]
                        [ 0.06814074,  0.08970813,  0.0709805 ]
                        [ 0.10852473,  0.10442963,  0.09596184]
                        [ 0.001     ,  0.06864924,  0.06446685]
                        [ 0.10736275,  0.10860951,  0.10390776]
                        [ 0.1054862 ,  0.13453259,  0.14436399]
                        [ 0.001     ,  0.09935915,  0.09152904]
                        [ 0.16389565,  0.20567791,  0.17598681]
                        [ 0.27533702,  0.2921259 ,  0.29065098]
                        [ 0.1191101 ,  0.13152554,  0.14642511]
                        [ 0.18992019,  0.21840044,  0.22218371]
                        [ 0.5       ,  0.5       ,  0.5       ]
                        [ 0.12792488,  0.17321354,  0.17124345]
                        [ 0.15380799,  0.20500506,  0.18252061]
                        [ 0.28461651,  0.29838497,  0.29824387]
                        [ 0.11898439,  0.12983076,  0.13957817]
                        [ 0.08328065,  0.09690795,  0.10125554]
                        [ 0.16736073,  0.16926813,  0.14419627]
                        [ 0.001     ,  0.09319319,  0.08569872]
                        [ 0.001     ,  0.07557973,  0.07155787]
                        [ 0.06218994,  0.09975598,  0.09463895]
                        [ 0.001     ,  0.06085632,  0.06879198]
                        [ 0.09145113,  0.11695828,  0.11041626]
                        [ 0.13902165,  0.15554516,  0.15806517]
                        [ 0.09289564,  0.09457757,  0.09430651]
                        [ 0.10877774,  0.13962074,  0.13614011]
                        [ 0.15946076,  0.18927357,  0.19282319]
                        [ 0.001     ,  0.10141351,  0.10644415]
                        [ 0.10346137,  0.11518805,  0.11624401]
                        [ 0.16234664,  0.16910338,  0.15386744]
                        [ 0.10064509,  0.09295108,  0.09645042]
                        [ 0.05356004,  0.07904283,  0.07449233]
                        [ 0.0862368 ,  0.09445628,  0.09213172]
                        [ 0.001     ,  0.06796732,  0.06162933]]
    
    assert_equal(results, expected_ww_x[:])
    
def test_magic_e_total():
        # create a basic meshtally
    meshtal_file = "mcnp_meshtal_single_meshtal.txt"
    tags = {4: ["n_result", "n_rel_error",
                "n_total_result", "n_total_rel_error"]}
    meshtal_object = mcnp.Meshtal(meshtal_file, tags, meshes_have_mats=True)
    
    results = magic(meshtal_object.tally[4], "n_total_result", 
                    "n_total_rel_error", .2, .001)
                
    expected_ww_x = [   0.0727709 ,  0.09303473,  0.05973245,  0.11605089,  
                        0.15358968,  0.09446185,  0.13288133,  0.20116779,  
                        0.11247224,  0.11745254,  0.15682268,  0.09338495,
                        0.07237543,  0.09706171,  0.06452718,  0.10439888,  
                        0.14226393,  0.09148097,  0.17794095,  0.2902493 ,  
                        0.14431317,  0.22078935,  0.5       ,  0.16993262,
                        0.1833384 ,  0.29779333,  0.13810424,  0.10030041,  
                        0.14697729,  0.08678777,  0.07130232,  0.09394698,  
                        0.06765712,  0.11029486,  0.15721929,  0.09428047,
                        0.13549031,  0.19141024,  0.10496243,  0.11572691,
                        0.15536779,  0.09631373,  0.07414555,  0.09211715,  
                        0.06157482]

    assert_equal(results, expected_ww_x[:])

