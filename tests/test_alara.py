"""alara module tests"""
import os
import unittest
import nose

from nose.tools import assert_almost_equal
from nose.tools import assert_equal
from numpy.testing import assert_array_equal

thisdir = os.path.dirname(__file__)
# mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest
    HAVE_PYTAPS = False
    pass

from pyne.mesh import Mesh, StatMesh, MeshError

def test_write_fluxin_single():
    mesh = Mesh(structured=True, structured_coords=[[0,1,2],[0,1,2],[0,1]])
    flux_tag = mesh.createTag("flux", 1, float)
    flux_data = [1, 2, 3, 4]
    
