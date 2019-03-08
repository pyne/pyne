"""
Tests for PyNE partisn module.
"""
import warnings
import os
import numpy as np
from numpy.testing import assert_array_almost_equal
import filecmp
from nose.tools import assert_almost_equal
from nose.plugins.skip import SkipTest
import multiprocessing
import unittest

from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB
from pyne.utils import QAWarning
from pyne import partisn

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

    mat_lib, unique_names = partisn._get_material_lib(
        hdf5, data_hdf5path, nuc_hdf5path, nuc_names=names)
    mat_lib_expected = {u'MERCURY1': {800000000: 4.066613534078662e-2},
                        u'HELIUMNA': {20040000: 2.4975599277878773e-05,
                                      20030000: 4.4414858514189387e-11}}
    expected_unique_names = {
        'mat:Helium, Natural': 'HELIUMNA', 'mat:Mercury': 'MERCURY1'}
    assert(unique_names == expected_unique_names)
    assert(mat_lib == mat_lib_expected)


def test_get_material_lib_no_names():
    """Test get_material_lib without a provided nuc_names list.
    """

    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    mat_lib, unique_names = partisn._get_material_lib(
        hdf5, data_hdf5path, nuc_hdf5path)
    mat_lib_expected = {'MERCURY1': {802020000: 1.2060451913893048e-02,
                                     802000000: 9.423512145483618e-03,
                                     802010000: 5.3498985962366465e-03,
                                     801960000: 6.24414427454006e-05,
                                     802040000: 2.7475463582858147e-03,
                                     801980000: 4.108325935058038e-03,
                                     801990000: 6.916609590819954e-03},
                        'HELIUMNA': {20040000: 2.4975599277878773e-05,
                                     20030000: 4.4414858514189387e-11}}
    expected_unique_names = {
        'mat:Helium, Natural': 'HELIUMNA', 'mat:Mercury': 'MERCURY1'}
    assert(unique_names == expected_unique_names)
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

    mat_xs_names_expected = {'mat:M1': {'cu63': 0.058, 'cu65': 0.026},
                             'mat:M2': {'h1': 0.067, 'o16': 0.033}}
    mat_xs_names = partisn._nucid_to_xs(mat_lib, nuc_names=nuc_names)
    print(mat_xs_names)
    assert(mat_xs_names == mat_xs_names_expected)


def test_nucid_to_xs_no_names():
    """Test the _nucid_to_xs function without a nuc_names dictionary.
    """
    mat_lib = {'mat:M1': {290630000: 0.058, 290650000: 0.026},
               'mat:M2': {10010000: 0.067, 80160000: 0.033}}
    mat_xs_names_expected = {'mat:M1': {u'Cu63': 0.058, u'Cu65': 0.026},
                             'mat:M2': {u'H1': 0.067, u'O16': 0.033}}
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

    if not HAVE_PYMOAB:
        raise SkipTest

    # Create mesh
    xvals = [0.0, 2.0]
    yvals = [0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Expected values
    igeom_expected = 'slab'
    bounds_expected = {'z': [-1.0, 0.0, 1.0]}

    igeom, bounds = partisn._get_coord_sys(mesh)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def test_get_coord_sys_2D():
    """Test the _get_coord_sys function for a 2D mesh.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Expected values
    igeom_expected = 'x-y'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'z': [-1.0, 0.0, 1.0]}

    igeom, bounds = partisn._get_coord_sys(mesh)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def test_get_coord_sys_3D():
    """Test the _get_coord_sys function for a 3D mesh.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    # Create mesh
    xvals = [-1.0, 0.0, 2.0]
    yvals = [-3.0, 0.0, 3.0]
    zvals = [-1.0, 0.0, 1.0]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Expected values
    igeom_expected = 'x-y-z'
    bounds_expected = {'x': [-1.0, 0.0, 2.0], 'y': [-3.0, 0.0, 3.0],
                       'z': [-1.0, 0.0, 1.0]}

    igeom, bounds = partisn._get_coord_sys(mesh)
    assert(igeom == igeom_expected)
    assert(bounds == bounds_expected)


def get_zones_no_void():
    """Test the _get_zones function if no void is in the meshed area.
    """

    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Load dagmc geometry
    dagmc.load(hdf5)

    # mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')
    # more inputs
    bounds = {'x': xvals, 'y': yvals, 'z': zvals}
    num_rays = 144
    grid = True
    dg = None
    mat_assigns = None

    unique_names = {'mat:Helium, Natural': 'HELIUMNA',
                    'mat:Mercury': 'MERCURY1'}
    voxel_zones, zones = partisn._get_zones(
        mesh, hdf5, bounds, num_rays, grid, dg, mat_assigns, unique_names)

    voxel_zones_expected = np.array([[1, 1, 1, 1],
                                     [2, 2, 2, 2],
                                     [3, 3, 3, 3]])
    zones_expected = {1: {'vol_frac': [1.0], 'mat': [u'HELIUMNA']},
                      2: {'vol_frac': [0.5, 0.5], 'mat': [u'HELIUMNA', u'MERCURY1']},
                      3: {'vol_frac': [1.0], 'mat': [u'MERCURY1']}}

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
    p.close()
    p.join()
    assert(r.get() == [True, True])


def get_zones_iteration_order():
    """Test that _get_zones function gives results in zyx order.
    """
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/fractal_box.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Load dagmc geometry
    dagmc.load(hdf5)

    bounds = [-5., 0., 5.]
    sc = [bounds]*3
    mesh = Mesh(structured_coords=sc, structured=True)
    bounds = {'x': bounds, 'y': bounds, 'z': bounds}
    num_rays = 9
    grid = True
    dg = None
    mat_assigns = None
    unique_names = {'mat:m1': 'M1', 'mat:m2': 'M2',
                    'mat:m3': 'M3', 'mat:m4': 'M4'}

    voxel_zones, zones = partisn._get_zones(
        mesh, hdf5, bounds, num_rays, grid, dg, mat_assigns, unique_names)

    voxel_zones_expected = np.array([[1, 2],
                                     [1, 4],
                                     [1, 3],
                                     [1, 4]])
    zones_expected = {1: {'vol_frac': [1.0], 'mat': ['M2']},
                      2: {'vol_frac': [1.0], 'mat': ['M4']},
                      3: {'vol_frac': [1.0], 'mat': ['M1']},
                      4: {'vol_frac': [1.0], 'mat': ['M3']}}

    vz_tf = voxel_zones.all() == voxel_zones_expected.all()
    z_tf = zones == zones_expected
    return [vz_tf, z_tf]


def test_get_zones_iteration_order():
    """Test that _get_zones function gives results in zyx order.
    """
    p = multiprocessing.Pool()
    r = p.apply_async(get_zones_iteration_order)
    p.close()
    p.join()
    assert(r.get() == [True, True])


def get_zones_with_void():
    """Test the _get_zones function if a void is present.
    """
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    dagmc.load(hdf5)

    # mesh
    xvals = [-5., 0., 10., 15., 15.1]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')
    # more inputs
    bounds = {'x': xvals, 'y': yvals, 'z': zvals}
    num_rays = 400
    grid = True
    dg = None
    mat_assigns = None

    unique_names = {'mat:Helium, Natural': 'HELIUMNA',
                    'mat:Mercury': 'MERCURY1'}
    voxel_zones, zones = partisn._get_zones(
        mesh, hdf5, bounds, num_rays, grid, dg, mat_assigns, unique_names)

    # expected results
    voxel_zones_expected = np.array([[1, 1, 1, 1],
                                     [2, 2, 2, 2],
                                     [3, 3, 3, 3],
                                     [0, 0, 0, 0]])
    zones_expected = {1: {'vol_frac': [1.0], 'mat': [u'HELIUMNA']},
                      2: {'vol_frac': [0.5, 0.5], 'mat': [u'HELIUMNA', u'MERCURY1']},
                      3: {'vol_frac': [1.0], 'mat': [u'MERCURY1']}}

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
    p.close()
    p.join()
    assert(r.get() == [True, True])


def test_check_fine_mesh_total_true():
    """Check that if fine mesh is less than 7, warning is issued.
    """
    block01 = {'it': 2, 'jt': 4}
    with warnings.catch_warnings(record=True) as w:
        partisn._check_fine_mesh_total(block01)
        assert(len(w) == 1)


def test_check_fine_mesh_total_false():
    """Check that if fine mesh is greater than 7, warning is not issued.
    """
    block01 = {'it': 2, 'jt': 4, 'kt': 4}
    with warnings.catch_warnings(record=True) as w:
        partisn._check_fine_mesh_total(block01)
    assert(len(w) == 0)


def write_partisn_input_1D():
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 5.]
    zvals = [-5., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_1D.inp'

    # other inputs
    ngroup = 5

    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_1D_expected.inp'

    partisn.write_partisn_input(mesh, hdf5, ngroup,
                                data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path,
                                input_file=input_file, num_rays=100, grid=True)

    out = filecmp.cmp(input_file, file_expected)
    os.remove(input_file)
    assert(out == True)
    return out


def test_write_partisn_input_1D():
    """Test full input file creation for 1D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_1D)
    p.close()
    p.join()
    assert(r.get() == True)


def write_partisn_input_2D():
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_2D.inp'

    # other inputs
    ngroup = 5

    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_2D_expected.inp'

    partisn.write_partisn_input(mesh, hdf5, ngroup,
                                data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path,
                                input_file=input_file, num_rays=100, grid=True)

    out = filecmp.cmp(input_file, file_expected)
    os.remove(input_file)
    return out


def test_write_partisn_input_2D():
    """Test full input file creation for 2D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_2D)
    p.close()
    p.join()
    assert(r.get() == True)


def write_partisn_input_3D():
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 0., 5.]
    zvals = [-5., 0., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
                structured_ordering='xyz')

    # Path for output file
    input_file = THIS_DIR + '/files_test_partisn/partisn_3D.inp'

    # other inputs
    ngroup = 5

    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_3D_expected.inp'

    partisn.write_partisn_input(mesh, hdf5, ngroup,
                                data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path,
                                input_file=input_file, num_rays=100, grid=True)

    out = filecmp.cmp(input_file, file_expected)
    os.remove(input_file)
    return out


def test_write_partisn_input_3D():
    """Test full input file creation for 3D case
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_3D)
    p.close()
    p.join()
    assert(r.get() == True)


def write_partisn_input_with_names_dict():
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    # Path to hdf5 test file
    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = THIS_DIR + '/files_test_partisn/partisn_test_geom.h5m'
    data_hdf5path = '/materials'
    nuc_hdf5path = '/nucid'

    # Create mesh
    xvals = [-5., 0., 10., 15.]
    yvals = [-5., 5.]
    zvals = [-5., 5.]
    mesh = Mesh(structured_coords=[xvals, yvals, zvals], structured=True,
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

    # expected output file
    file_expected = THIS_DIR + '/files_test_partisn/partisn_nucnames_expected.inp'

    partisn.write_partisn_input(mesh, hdf5, ngroup,
                                data_hdf5path=data_hdf5path, nuc_hdf5path=nuc_hdf5path,
                                input_file=input_file, num_rays=100, grid=True, names_dict=names)

    out = filecmp.cmp(input_file, file_expected)
    os.remove(input_file)
    return out


def test_write_partisn_input_with_names_dict():
    """Test full input file creation for 1D case with a names_dict provided
    """
    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_with_names_dict)
    p.close()
    p.join()
    assert(r.get() == True)


def write_partisn_input_options():
    if not HAVE_PYMOAB:
        raise SkipTest

    """Test PARTISN input file creation with a slew of keyword arguments
    """

    THIS_DIR = os.path.dirname(os.path.realpath(__file__))
    hdf5 = os.path.join(THIS_DIR, 'files_test_partisn',
                        'partisn_test_geom.h5m')
    input_file = os.path.join(
        THIS_DIR, 'files_test_partisn', 'partisn_options.inp')
    file_expected = os.path.join(
        THIS_DIR, 'files_test_partisn', 'partisn_options_expected.inp')

    sc = [-5., 0., 10., 15.], [-5., 5.], [-5., 5.]
    mesh = Mesh(structured_coords=sc, structured=True)
    ngroup = 66

    dg = np.zeros(4, dtype=[('idx', np.int64),
                            ('cell', np.int64),
                            ('vol_frac', np.float64),
                            ('rel_error', np.float64)])
    dg[:] = [(0, 1, 1.0, 0.0), (1, 1, 0.5, 0.04714045207910317),
             (1, 2, 0.5, 0.04714045207910317), (2, 2, 1.0, 0.0)]
    mat_assigns = {1: 'mat:Helium, Natural', 2: 'mat:Mercury',
                   5: 'mat:Graveyard', 6: u'mat:Vacuum'}

    cards = {"block1": {"isn": 6,
                        "maxscm": 3000000,
                        "maxlcm": 6000000,
                        },
             "block2": {"hello": "from block2"},
             "block3": {"lib": "xsf21-71",
                        "lng": 175,
                        "maxord": 5,
                        "ihm": 227,
                        "iht": 10,
                        "ihs": 11,
                        "ifido": 1,
                        "ititl": 1,
                        "i2lp1": 0,
                        "savbxs": 1,
                        "kwikrd": 1
                        },
             "block4": {"hello": "from block4"},
             "block5": {"source": "<this is a dummy source>"}
             }

    with warnings.catch_warnings(record=True) as w:
        partisn.write_partisn_input(mesh, hdf5, ngroup, input_file=input_file,
                                    dg=dg, mat_assigns=mat_assigns, fine_per_coarse=3,
                                    cards=cards, num_rays=9)  # include num_rays to get warning

    # verify we get a warning from including num_rays and dg
    out1 = len(w) == 1
    out2 = filecmp.cmp(input_file, file_expected)
    os.remove(input_file)
    return out1 and out2


def test_write_partisn_input_options():
    """Test full input file creation for 1D case with a lot of key work args
    """

    p = multiprocessing.Pool()
    r = p.apply_async(write_partisn_input_options)
    p.close()
    p.join()
    assert(r.get() == True)


def test_format_repeated_vector():
    """Test the format_repeated_vector function.
    """
    vector = [1.0, 1, 1, 0, 0, 1, 2, 2.0, 2.1]
    string_expected = "3R 1.0 2R 0 1 2R 2 2.1 "
    string = partisn.format_repeated_vector(vector)

    assert(string == string_expected)


def test_mesh_to_isotropic_source():
    """Test isotropic SOURCF generation.
    """
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    m = Mesh(structured=True, structured_coords=[range(5), range(5), range(5)])
    m.src = NativeMeshTag(4, float)
    # These source values were carefully choosen so that:
    # 1. The iteration order could be visually checked based on RTFLUX output using
    #    the resulting SOURCF card.
    # 2. The repeation capability (i.e. 3R 0 = 0 0 0) could be tested.
    m.src[:] = [[0, 0, 0, 100], [0, 0, 0, 0], [0, 0, 0, 6], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [
                    0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 100, 0], [0, 0, 0, 5], [0, 0, 0, 6], [
                    0, 0, 0, 0], [0, 0, 0, 8], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [
                    0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 100, 0, 0], [0, 0, 0, 5], [0, 0, 0, 0], [
                    0, 0, 0, 7], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [
                    0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [100, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [
                    0, 0, 0, 7], [0, 0, 0, 8], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [
                    0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    out = partisn.mesh_to_isotropic_source(m, "src")
    exp = ("sourcf= 1.00000E+02 3R 0; 0 8.00000E+00 0 8.00000E+00; 4R 0; 4R 0; 0 5.00000E+00\n"
           "5.00000E+00 0; 4R 0; 4R 0; 4R 0; 6.00000E+00 6.00000E+00 2R 0; 4R 0; 4R 0; 4R 0;\n"
           "2R 0 7.00000E+00 7.00000E+00; 4R 0; 4R 0; 4R 0; 0 1.00000E+02 2R 0; 4R 0; 4R 0;\n"
           "4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 2R\n"
           "0 1.00000E+02 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R\n"
           "0; 4R 0; 4R 0; 4R 0; 4R 0; 3R 0 1.00000E+02; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0;\n"
           "4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0; 4R 0;")

    assert(out == exp)


def test_isotropic_vol_source():
    """Test isotropic volumetric source generation from DAGMC geometry.
    """
    try:
        from pyne import dagmc
    except:
        raise SkipTest

    if not HAVE_PYMOAB:
        raise SkipTest

    sc = np.linspace(-25, 25, 6)
    m = Mesh(structured=True, structured_coords=[sc, sc, sc])

    cells = [14, 15]
    spectra = [[0.1, 0.1, 0.1, 0.7], [0.3, 0.3, 0.3, 0.1]]
    intensities = [1, 2]

    dg, s = partisn.isotropic_vol_source("files_test_partisn/source_boxes.h5m",
                                         m, cells, spectra, intensities,
                                         num_rays=4, tag_name="src", grid=True)
    m.src = NativeMeshTag(4, float)
    data = m.src[:]

    # setup expected result, confirmed by hand calcs and inspection
    exp = np.zeros(shape=((len(sc) - 1)**3, len(spectra[0])))
    exp[22, :] = [0.025, 0.025,  0.025,  0.175]
    exp[62, :] = [0.075, 0.075, 0.075, 0.025]
    exp[63, :] = [0.075, 0.075, 0.075, 0.025]
    exp[87, :] = [0.075, 0.075, 0.075, 0.025]
    exp[88, :] = [0.075, 0.075, 0.075, 0.025]

    assert(np.allclose(data, exp))
