"""
Tests for PyNE partisn module.
"""
import warnings
import os
import numpy as np
import filecmp
from nose.tools import assert_almost_equal
from numpy.testing import assert_array_almost_equal
from pyne import partisn
from pyne.utils import QAWarning
import multiprocessing
import unittest

try:
    from itaps import iBase, iMesh, iMeshExtensions
    HAVE_PYTAPS = True
except ImportError:
    HAVE_PYTAPS = False
    from nose.plugins.skip import SkipTest
    raise SkipTest

if HAVE_PYTAPS:
    from pyne.mesh import Mesh

warnings.simplefilter("ignore", QAWarning)


def test_get_material_lib_with_names():
    """Test get_material_lib with a provided nuc_names list.
    """
    
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # nuc_names list
    names = {}
    names[800000000] = 'hg'
    names[20030000] = 'he3'
    names[20040000] = 'he4'
    
    mat_lib = partisn._get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, nuc_names=names)
    mat_lib_expected = {u'mat:Mercury':{800000000:4.066613534078662e-2}, 
                        u'mat:Helium, Natural':{20040000:2.4975599277878773e-05, 
                                                20030000:4.4414858514189387e-11}}
    assert(mat_lib == mat_lib_expected)


def test_get_material_lib_no_names():
    """Test get_material_lib without a provided nuc_names list.
    """
    
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    mat_lib = partisn._get_material_lib(hdf5, data_hdf5path, nuc_hdf5path)
    mat_lib_expected = {'mat:Mercury': {802020000:1.2060451913893048e-02, 
                                        802000000:9.423512145483618e-03,
                                        802010000:5.3498985962366465e-03, 
                                        801960000:6.24414427454006e-05, 
                                        802040000:2.7475463582858147e-03, 
                                        801980000:4.108325935058038e-03, 
                                        801990000:6.916609590819954e-03},
                        'mat:Helium, Natural':{20040000:2.4975599277878773e-05, 
                                            20030000:4.4414858514189387e-11}}
    assert(mat_lib == mat_lib_expected)


def test_nucid_to_xs_with_names():
    """Test the _nucid_to_xs function given a nuc_names dictionary.
    """
    mat_lib = {'mat:M1': {290630000: 0.058, 290650000: 0.026}, 
            'mat:M2': {10010000: 0.067, 80160000: 0.033}}
    
    nuc_names = {}
    nuc_names[290630000] = 'cu63'
    nuc_names[290650000] = 'cu65'
    nuc_names[10010000] = 'h1'
    nuc_names[80160000] = 'o16'
    
    mat_xs_names_expected = {'M1': {'cu63': 0.058, 'cu65': 0.026}, 
                            'M2': {'h1': 0.067, 'o16': 0.033}}
    mat_xs_names = partisn._nucid_to_xs(mat_lib, nuc_names = nuc_names)
    assert(mat_xs_names == mat_xs_names_expected)

    
def test_nucid_to_xs_no_names():
    """Test the _nucid_to_xs function without a nuc_names dictionary.
    """
    mat_lib = {'mat:M1': {290630000: 0.058, 290650000: 0.026}, 
            'mat:M2': {10010000: 0.067, 80160000: 0.033}}
    mat_xs_names_expected = {'M1': {u'Cu63': 0.058, u'Cu65': 0.026}, 
                            'M2': {u'H1': 0.067, u'O16': 0.033}}
    mat_xs_names = partisn._nucid_to_xs(mat_lib)
    assert(mat_xs_names == mat_xs_names_expected)


def test_get_xs_names():
    """Test the _get_xs_names function.
    """
    # Create dictionary to pass
    mat_xs_names = {'M1': {'cu63': 0.058, 'cu65': 0.026}, 
                    'M2': {'h1': 0.067, 'o16': 0.033}}
    xs_names_expected = ['cu63', 'cu65', 'h1', 'o16']
    xs_names = partisn._get_xs_names(mat_xs_names)
    assert(sorted(xs_names_expected) == sorted(xs_names))
    

def test_get_coord_sys_1D():
    """Test the _get_coord_sys function for a 1D mesh.
    """
    # Create mesh
    xvals = [0.0, 2.0]
    yvals = [0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    pn = 2
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'slab'
    bounds_expected = {'z':[-1.0, 0.0, 1.0]}
    
    igeom, bounds, nmq = partisn._get_coord_sys(mesh, pn)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)
    assert(nmq == 3)


def test_get_coord_sys_2D():
    """Test the _get_coord_sys function for a 2D mesh.
    """
    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    pn = 2
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'x-y'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'z':[-1.0, 0.0, 1.0]}
    
    igeom, bounds, nmq = partisn._get_coord_sys(mesh, pn)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)
    assert(nmq == 6)


def test_get_coord_sys_3D():
    """Test the _get_coord_sys function for a 3D mesh.
    """
    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    pn = 2
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Expected values
    igeom_expected = 'x-y-z'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'y':[-3.0, 0.0, 3.0], 
                        'z':[-1.0, 0.0, 1.0]}
    
    igeom, bounds, nmq = partisn._get_coord_sys(mesh, pn)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)
    assert(nmq == 9)

   
def get_zones_no_void():
    """Test the _get_zones function if no void is in the meshed area.
    """
    # hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # Load dagmc geometry
    from pyne import dagmc
    dagmc.load(hdf5)
    
    # mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    # more inputs
    bounds = {'x':xvals, 'y':yvals, 'z':zvals}
    num_rays = 144
    grid = True
    
    voxel_zones, zones = partisn._get_zones(mesh, hdf5, bounds, num_rays, grid)    

    voxel_zones_expected = np.array([[1, 1, 1, 1],
                                    [2, 2, 2, 2],
                                    [3, 3, 3, 3]])
    zones_expected = {1:{'vol_frac':[1.0], 'mat':[u'HeliumNatural']}, 
                    2:{'vol_frac':[0.5, 0.5], 'mat':[u'HeliumNatural', u'Mercury']}, 
                    3:{'vol_frac':[1.0], 'mat':[u'Mercury']}}
    
    vz_tf = False
    z_tf = False
    if voxel_zones.all() == voxel_zones_expected.all():
        vz_tf = True
    if zones == zones_expected:
        z_tf = True
    #assert(voxel_zones.all() == voxel_zones_expected.all())
    #assert(zones == zones_expected)
    return [vz_tf, z_tf]


def test_get_zones_no_void():
    """Test the _get_zones function if no void is in the meshed area.
    """
    p = multiprocessing.Pool()
    r = p.apply_async(get_zones_no_void)
    assert(r.get() == [True, True])

    
def get_zones_with_void():
    """Test the _get_zones function if a void is present.
    """
    # hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    from pyne import dagmc
    dagmc.load(hdf5)
    
    # mesh
    xvals = [-5., 0., 10., 15., 15.1]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    # more inputs
    bounds = {'x':xvals, 'y':yvals, 'z':zvals}
    num_rays = 400
    grid = True
    
    voxel_zones, zones = partisn._get_zones(mesh, hdf5, bounds, num_rays, grid)
    
    # expected results
    voxel_zones_expected = np.array([[1, 1, 1, 1],
                                    [2, 2, 2, 2],
                                    [3, 3, 3, 3],
                                    [0, 0, 0, 0]])
    zones_expected = {1:{'vol_frac':[1.0], 'mat':[u'HeliumNatural']}, 
                    2:{'vol_frac':[0.5, 0.5], 'mat':[u'HeliumNatural', u'Mercury']}, 
                    3:{'vol_frac':[1.0], 'mat':[u'Mercury']}}
    
    vz_tf = False
    z_tf = False
    if voxel_zones.all() == voxel_zones_expected.all():
        vz_tf = True
    if zones == zones_expected:
        z_tf = True
    #assert(voxel_zones.all() == voxel_zones_expected.all())
    #assert(zones == zones_expected)
    return [vz_tf, z_tf]
    
    
def test_get_zones_with_void():
    """Test the _get_zones function if a void is present.
    """
    p = multiprocessing.Pool()
    r = p.apply_async(get_zones_with_void)
    assert(r.get() == [True, True])


def test_check_fine_mesh_total_true():
    """Check that if fine mesh is less than 7, warning is issued.
    """
    block01 = {'it':2, 'jt':4}
    warn_out = partisn._check_fine_mesh_total(block01)
    assert(warn_out == True)
    

def test_check_fine_mesh_total_false():
    """Check that if fine mesh is greater than 7, warning is not issued.
    """
    block01 = {'it':2, 'jt':4, 'kt':4}
    warn_out = partisn._check_fine_mesh_total(block01)
    assert(warn_out == False)
    

def write_partisn_input_1D():
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 5.]
    zvals = [-5., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_1D.inp'
    
    # other inputs
    ngroup = 5
    pn = 2
    
    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_1D_expected.inp'
    
    partisn.write_partisn_input(mesh, hdf5, ngroup, pn, 
        data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path, 
        input_file=input_file, num_rays=100, grid=True)
    
    out = filecmp.cmp(input_file, file_expected)
    assert(out == True)
    return out
        
        
def test_write_partisn_input_1D():
    """Test full input file creation for 1D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_1D)
    assert(r.get() == True)


def write_partisn_input_2D():
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_2D.inp'
    
    # other inputs
    ngroup = 5
    pn = 2
    
    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_2D_expected.inp'
    
    partisn.write_partisn_input(mesh, hdf5, ngroup, pn, 
        data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path, 
        input_file=input_file, num_rays=100, grid=True)
    
    out = filecmp.cmp(input_file, file_expected)
    return out

        
def test_write_partisn_input_2D():
    """Test full input file creation for 2D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_2D)
    assert(r.get() == True)


def write_partisn_input_3D():
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_3D.inp'
    
    # other inputs
    ngroup = 5
    pn = 2
    
    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_3D_expected.inp'
    
    partisn.write_partisn_input(mesh, hdf5, ngroup, pn, 
        data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path, 
        input_file=input_file, num_rays=100, grid=True)
    
    out = filecmp.cmp(input_file, file_expected)
    return out


def test_write_partisn_input_3D():
    """Test full input file creation for 3D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_3D)
    assert(r.get() == True)


def write_partisn_input_with_names_dict():
    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'
    
    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 5.]
    zvals = [-5., 5.]
    mesh=Mesh(structured_coords=[xvals, yvals, zvals], structured=True, 
                structured_ordering='xyz')
    
    # nuc_names list
    names = {}
    names[800000000] = 'hg'
    names[20030000] = 'he3'
    names[20040000] = 'he4'
    
    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_nucnames.inp'
    
    # other inputs
    ngroup = 5
    pn = 2
    
    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_nucnames_expected.inp'
    
    partisn.write_partisn_input(mesh, hdf5, ngroup, pn, 
        data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path, 
        input_file=input_file, num_rays=100, grid=True, names_dict=names)
    
    out = filecmp.cmp(input_file, file_expected)
    return out
    
    
def test_write_partisn_input_with_names_dict():
    """Test full input file creation for 1D case with a names_dict provided
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_with_names_dict)
    assert(r.get() == True)
    
    
def test_format_repeated_vector():
    """Test the format_repeated_vector function.
    """
    vector = [1.0, 1, 1, 0, 0, 1, 2, 2.0, 2.1]
    string_expected = "3R 1.0 2R 0 1 2R 2 2.1 "
    string = partisn.format_repeated_vector(vector)

    assert(string == string_expected)


def test_strip_mat_name():
    """Test the strip_mat_name function.
    """
    name = 'mat:Helium, Natural'
    mat_name_expected = 'HeliumNatural'
    mat_name = partisn.strip_mat_name(name)
    
    assert(mat_name_expected == mat_name)

