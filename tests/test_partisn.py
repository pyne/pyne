"""
Tests for PyNE partisn module.
"""
import warnings

try:
    from itaps import iBase, iMesh, iMeshExtensions
    HAVE_PYTAPS = True
except ImportError:
    HAVE_PYTAPS = False
    from nose.plugins.skip import SkipTest
    raise SkipTest

if HAVE_PYTAPS:
    from pyne.mesh import Mesh

from nose.tools import assert_almost_equal
from numpy.testing import assert_array_almost_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import partisn


def test_get_material_lib_with_names():
    pass
    

def test_get_material_lib_no_names():
    pass


def test_nucid_to_xs_with_names():
    """Test the _nucid_to_xs function given a nuc_names dictionary.
    """
    
    # Create mat_lib and nuc_names
    mat_lib = {'mat:M1': {290630000: 0.058, 290650000: 0.026}, 
               'mat:M2': {10010000: 0.067, 80160000: 0.033}}
    
    nuc_names = {}
    nuc_names[290630000] = 'cu63'
    nuc_names[290650000] = 'cu65'
    nuc_names[10010000] = 'h1'
    nuc_names[80160000] = 'o16'
    
    # Expected mat_xs_lib:
    mat_xs_names_expected = {'mat:M1': {'cu63': 0.058, 'cu65': 0.026}, 
                             'mat:M2': {'h1': 0.067, 'o16': 0.033}}
    
    # Call function
    mat_xs_names = partisn._nucid_to_xs(mat_lib, nuc_names = nuc_names)
    
    # Check that material dictionaries are equal:
    assert(mat_xs_names == mat_xs_names_expected)

    
def test_nucid_to_xs_no_names():
    """Test the _nucid_to_xs function without a nuc_names dictionary.
    """
    
    # Create mat_lib
    mat_lib = {'mat:M1': {290630000: 0.058, 290650000: 0.026}, 
               'mat:M2': {10010000: 0.067, 80160000: 0.033}}
    
    # Expected mat_xs_lib:
    mat_xs_names_expected = {'mat:M1': {u'Cu63': 0.058, u'Cu65': 0.026}, 
                             'mat:M2': {u'H1': 0.067, u'O16': 0.033}}
    
    # Call function
    mat_xs_names = partisn._nucid_to_xs(mat_lib)
    
    # Check that material dictionaries are equal:
    assert(mat_xs_names == mat_xs_names_expected)


def test_get_xs_names():
    """Test the _get_xs_names function.
    """
    
    # Create dictionary to pass
    mat_xs_names = {'mat:M1': {'cu63': 0.058, 'cu65': 0.026}, 
                    'mat:M2': {'h1': 0.067, 'o16': 0.033}}
    
    # Expected xs_names
    xs_names_expected = ['cu63', 'cu65', 'h1', 'o16']
    
    # Call function
    xs_names = partisn._get_xs_names(mat_xs_names)
    
    # Check that lists are equal
    assert(sorted(xs_names_expected) == sorted(xs_names))
    

def test_get_coord_sys_1D():
    """Test the _get_coord_sys function for a 1D mesh.
    """
    
    # Create mesh
    xvals = [0.0, 2.0]
    yvals = [0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'SLAB'
    bounds_expected = {'z':[-1.0, 0.0, 1.0]}
    
    # Call function
    igeom, bounds = partisn._get_coord_sys(mesh)
    
    # Test that output is equal
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def test_get_coord_sys_2D():
    """Test the _get_coord_sys function for a 2D mesh.
    """
    
    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'X-Y'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'z':[-1.0, 0.0, 1.0]}
    
    # Call function
    igeom, bounds = partisn._get_coord_sys(mesh)
    
    # Test that output is equal
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def test_get_coord_sys_3D():
    """Test the _get_coord_sys function for a 3D mesh.
    """
    
    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'X-Y-Z'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'y':[-3.0, 0.0, 3.0], 
                        'z':[-1.0, 0.0, 1.0]}
    
    # Call function
    igeom, bounds = partisn._get_coord_sys(mesh)
    
    # Test that output is equal
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def test_get_zones():
    pass


def test_write_partisn_input_1D():
    pass
    

def test_write_partisn_input_2D():
    pass


def test_write_partisn_input_3D():
    pass


def test_format_repeated_vector():
    """Test the format_repeated_vector function.
    """
    # Create vector
    vector = [1.0, 1, 1, 0, 0, 1, 2, 2.0, 2.1]
    
    # Expected string
    string_expected = "3R 1.0 2R 0 1 2R 2 2.1 "
    
    # Call function
    string = partisn.format_repeated_vector(vector)
    
    # Test if equal
    assert(string == string_expected)
