from __future__ import print_function
import sys
import unittest
import os.path
import warnings
from nose.tools import assert_equal, assert_almost_equal, assert_raises, assert_true
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_equal
import imp
import multiprocessing

try:
    from itaps import iMesh
    HAVE_IMESH = True
except ImportError:
    HAVE_IMESH = False

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

try:
    from pyne.mesh import Mesh
    # See if dagmc module exists but do not import it
    pyne_info = imp.find_module('pyne')
    pyne_mod = imp.load_module('pyne', *pyne_info)
    imp.find_module('dagmc', pyne_mod.__path__)
except ImportError:
    raise SkipTest

if sys.version_info[0] < 3:
    STRING_TYPES = (basestring, str, unicode)
else:
    STRING_TYPES = (str,)

path = os.path.join(os.path.dirname(__file__), 'unitbox.h5m')

# use extra underscore to ensure this function is first in alpabetical
# sorted order, because it must run before the others.
def load():
    from pyne import dagmc
    # FIXME laoding causes infor to be printied to stderr (or stdout).
    dagmc.load(path)

def metadata():
    from pyne import dagmc
    dagmc.load(path)
    
    rets = [dagmc.volume_is_graveyard(x) for x in range(1, 5)]
    assert_equal(rets, [True, False, False, False])

    rets = [dagmc.volume_is_implicit_complement(x) for x in range(1, 5)]
    assert_equal(rets, [False, False, False, True])

    md = dagmc.volume_metadata(2)
    assert_equal(md['material'], 5)
    assert_equal(md['rho'], 0.5)
    assert_true(all(x in md for x in ['material', 'rho', 'imp']))

def versions():
    from pyne import dagmc
    dagmc.load(path)
    
    returned = dagmc.versions()
    assert_equal(len(returned), 2)
    assert_true(isinstance(returned[0], STRING_TYPES))
    assert_true(isinstance(returned[1], int))

def list_functions():
    from pyne import dagmc
    dagmc.load(path)
    
    surfs = dagmc.get_surface_list()
    assert_equal(set(surfs), set(range(1, 19)))

    vols = dagmc.get_volume_list()
    assert_equal(set(vols), set(range(1, 5)))

def boundary():
    from pyne import dagmc
    dagmc.load(path)
    
    low, high = dagmc.volume_boundary(2)
    for i in range(0, 3):
        assert_true(low[i] <= -1.0)
        assert_true(high[i] >= 1.0)

def pt_in_vol():
    from pyne import dagmc
    dagmc.load(path)
    
    # there are 4 volumes; (0,0,.2) is in volume 2
    rets = [dagmc.point_in_volume(x, [0, 0, .2]) for x in range(1, 5)]
    assert_equal(rets, [False, True, False, False])

    # (1.1,0,0) is in volume 3
    rets = [dagmc.point_in_volume(x, [1.1, 0, 0]) for x in range(1, 5)]
    assert_equal(rets, [False, False, True, False])

def find_volume():
    from pyne import dagmc
    dagmc.load(path)
    
    vol = dagmc.find_volume([0, 0, 0])
    assert_equal(vol, 2)

    vol = dagmc.find_volume([.9, .9, .9])
    assert_equal(vol, 2)

    # boundary case -- point [1,.1,.1] is on surface between vols 2 and 3
    # the behavior on boundaries will vary with uvw, but ensure that
    # only the two volumes the touch the location are ever returned.
    for uvw in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]:
        vol = dagmc.find_volume([1, .1, .1], uvw)
        assert_true(vol in (2, 3))

    # boundary case-- exiting volume 3 => in volume 3
    vol = dagmc.find_volume([1, .1, .1], [-1, 0, 0])
    assert_equal(vol, 3)

    vol = dagmc.find_volume([1.1, 0, 0])
    assert_equal(vol, 3)

def one_ray():
    from pyne import dagmc
    dagmc.load(path)
    
    fromcenter = dagmc.fire_one_ray(2, [0, 0, 0], [1, 0, 0])
    assert_almost_equal(fromcenter[1], 1.0)

    fromhalf = dagmc.fire_one_ray(2, [.5, .5, .5], [0, -1, 0])
    assert_almost_equal(fromhalf[1], 1.5)

    # the following is actually a misuse of ray_fire because it should only
    # be called from inside a volume, but at least verify that a ray is
    # lost if it 1) starts outside a volume and 2) does not intersect the
    # volume
    fromoutside = dagmc.fire_one_ray(2, [0, 1.1, 0], [-1, 0, 0])
    assert_equal(fromoutside, None)

    fromvol3 = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, -1, 0])
    assert_almost_equal(fromvol3[1], 0.1)

    fromvol3 = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, 1, 0])
    assert_almost_equal(fromvol3[1], 3.056921938)

def failures():
    from pyne import dagmc
    dagmc.load(path)
    
    assert_raises(Exception, dagmc.point_in_volume, [100, (0, 0, 0)])
    assert_raises(Exception, dagmc.point_in_volume, [1, (0, 0, 0, 0)])
    assert_raises(Exception, dagmc.fire_one_ray, [2, (0, 0, 0), 1])

def ray_iterator():
    from pyne import dagmc
    dagmc.load(path)
    
    start = [-2, 0, 0]
    startvol = dagmc.find_volume(start)
    assert_equal(startvol, 3)
    direction = [1, 0, 0]

    expected_vols = [2, 3, 1, 4]
    expected_dists = [1, 2, 3.156921938, 0]

    for i, (vol, dist, surf) in enumerate(
            dagmc.ray_iterator(startvol, start, direction)):
        assert_equal(expected_vols[i], vol)
        if expected_dists[i] != 0:
            assert_almost_equal(expected_dists[i], dist)
    assert_equal(i, 3)

    for i, (vol, dist, surf) in enumerate(
            dagmc.ray_iterator(startvol, start, direction, dist_limit=4)):
        assert_equal(expected_vols[i], vol)
        if expected_dists[i] != 0:
            assert_almost_equal(expected_dists[i], dist)
    assert_equal(i, 1)

def ray_story():
    from pyne import dagmc
    dagmc.load(path)
    
    # Run with `nosetests -s` to see output printed on stdout
    dagmc.tell_ray_story((0, 0, 0), (1, 1, 0))
    dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0))
    # tests edge behavior
    dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0), dist_limit=4)
    dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0), dist_limit=4.1)

def util_graveyard_bound():
    from pyne import dagmc
    dagmc.load(path)
    
    lo, hi = dagmc.find_graveyard_inner_box()

    grave_diam = 4.15692194
    for i in range(0, 3):
        assert_almost_equal(lo[i], -grave_diam)
        assert_almost_equal(hi[i], grave_diam)

def util_matlist():
    from pyne import dagmc
    dagmc.load(path)
    
    mats = dagmc.get_material_set()
    assert_equal(set((0, 5)), mats)

    mats = dagmc.get_material_set(with_rho=True)
    assert_equal(set([(0, 0.0), (5, 0.5)]), mats)

def discretize_geom_rand():
    from pyne import dagmc
    dagmc.load(path)
    
    coords = [-4, -1, 1, 4]
    mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
    results = dagmc.discretize_geom(mesh, num_rays=50)
    
    assert_equal(len(results), (len(coords) - 1)**3)

    for res in results:
        if res['idx'] != 13:
            assert_equal(res['cell'], 3)
        else:
            assert_equal(res['cell'], 2)

        assert_almost_equal(res['vol_frac'], 1.0)

def discretize_geom_grid():
    from pyne import dagmc
    dagmc.load(path)
    
    coords = [-4, -1, 1, 4]
    mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
    results = dagmc.discretize_geom(mesh, num_rays=49, grid=True)
    
    assert_equal(len(results), (len(coords) - 1)**3)

    for res in results:
        if res['idx'] != 13:
            assert_equal(res['cell'], 3)
        else:
            assert_equal(res['cell'], 2)

        assert_almost_equal(res['vol_frac'], 1.0)

def discretize_geom_mix():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    coords2 = [0, 2]
    mesh = Mesh(structured=True, 
                structured_coords=[coords2, coords, coords])
    results1 = dagmc.discretize_geom(mesh, num_rays=100, grid=True)

    assert_equal(results1[0]['cell'], 2)
    assert_almost_equal(results1[0]['vol_frac'], 0.5)

    assert_equal(results1[1]['cell'], 3)
    assert_almost_equal(results1[1]['vol_frac'], 0.5)

    
    # To to make sure standard error decreases with increasing rays
    results2 = dagmc.discretize_geom(mesh, num_rays=625, grid=True)
    assert(results2[0]['rel_error'] < results1[0]['rel_error'])
    assert(results2[1]['rel_error'] < results1[1]['rel_error'])

def discretize_non_square():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    mesh = Mesh(structured=True, 
                structured_coords=[coords, coords, coords])
    assert_raises(ValueError, dagmc.discretize_geom, mesh, num_rays=3, grid=True)

def discretize_geom_centers():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    coords2 = [0, 2, 4]
    mesh = Mesh(structured=True, 
                structured_coords=[coords2, coords, coords])
    #  explicitly set to unstructured to trigger the ve centers sampling
    #  method
    mesh.structured = False
    res = dagmc.discretize_geom(mesh)
    assert_array_equal(res["idx"], [0, 1])
    assert_array_equal(res["cell"], [2, 3])
    assert_array_equal(res["vol_frac"], [1.0, 1.0])
    assert_array_equal(res["rel_error"], [1.0, 1.0])

    #  ensure kwargs are not accepted for unstructured mesh
    assert_raises(ValueError, dagmc.discretize_geom, mesh, num_rays=3, grid=True)

def cells_at_ve_centers():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    coords2 = [0, 2, 4]
    mesh = Mesh(structured=True, 
                structured_coords=[coords2, coords, coords])
    cells = dagmc.cells_at_ve_centers(mesh)
    assert_array_equal(cells, [2, 3])

def test__load():
    p = multiprocessing.Pool()
    results = p.apply_async(load)
    r = results.get()
    
def test_metadata():
    p = multiprocessing.Pool()
    results = p.apply_async(metadata)
    r = results.get()

def test_versions():
    p = multiprocessing.Pool()
    results = p.apply_async(versions)
    r = results.get()

def test_versions():
    p = multiprocessing.Pool()
    results = p.apply_async(versions)
    r = results.get()

def test_list_functions():
    p = multiprocessing.Pool()
    results = p.apply_async(list_functions)
    r = results.get()

def test_boundary():
    p = multiprocessing.Pool()
    results = p.apply_async(boundary)
    r = results.get()

def test_pt_in_vol():
    p = multiprocessing.Pool()
    results = p.apply_async(pt_in_vol)
    r = results.get()

def test_find_volume():
    p = multiprocessing.Pool()
    results = p.apply_async(find_volume)
    r = results.get()

def test_one_ray():
    p = multiprocessing.Pool()
    results = p.apply_async(one_ray)
    r = results.get()

def test_failures():
    p = multiprocessing.Pool()
    results = p.apply_async(failures)
    r = results.get()

def test_ray_iterator():
    p = multiprocessing.Pool()
    results = p.apply_async(ray_iterator)
    r = results.get()

def test_ray_story():
    p = multiprocessing.Pool()
    results = p.apply_async(ray_story)
    r = results.get()

def test_util_graveyard_bounds():
    p = multiprocessing.Pool()
    results = p.apply_async(util_graveyard_bound)
    r = results.get()

def test_util_matlist():
    p = multiprocessing.Pool()
    results = p.apply_async(util_matlist)
    r = results.get()

def test_discretize_geom_rand():
    """The 14th (index 13) mesh volume element fully contains volume 2. Use 
    random sampling.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_rand)
    r = results.get()

def test_discretize_geom_grid():
    """The 14th (index 13) mesh volume element fully contains volume 2. Use 
    grid sampling.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_grid)
    r = results.get()

def test_discretize_geom_mix():
    """Single mesh volume element that is a 50:50 split of geometry volumes
    2 and 3.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_mix)
    r = results.get()

def test_descritize_non_square():
    """Test to make sure requesting a grid with a num_rays that is not a
    perfect square raises ValueError.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(discretize_non_square)
    r = results.get()

def test_discretize_geom_centers():
    """Test that unstructured mesh is sampled by mesh ve centers.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_centers)
    r = results.get()

def test_cells_at_ve_centers():
    """Test that a mesh with one ve in cell 2 and one ve in cell 3 produces
    correct results.
    """
    if not HAVE_IMESH:
        raise SkipTest
        
    p = multiprocessing.Pool()
    results = p.apply_async(cells_at_ve_centers)
    r = results.get()

