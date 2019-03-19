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
import numpy as np

from pyne.mesh import HAVE_PYMOAB

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

    rets1 = [dagmc.volume_is_graveyard(x) for x in range(1, 5)]
    rets2 = [dagmc.volume_is_implicit_complement(x) for x in range(1, 5)]
    md = dagmc.volume_metadata(2)

    return [rets1, rets2, md]

def versions():
    from pyne import dagmc
    dagmc.load(path)

    returned = dagmc.versions()
    return returned

def list_functions():
    from pyne import dagmc
    dagmc.load(path)

    surfs = dagmc.get_surface_list()
    vols = dagmc.get_volume_list()

    return [surfs, vols]

def boundary():
    from pyne import dagmc
    dagmc.load(path)

    low, high = dagmc.volume_boundary(2)
    return [low, high]

def pt_in_vol():
    from pyne import dagmc
    dagmc.load(path)

    # there are 4 volumes; (0,0,.2) is in volume 2
    rets1 = [dagmc.point_in_volume(x, [0, 0, .2]) for x in range(1, 5)]

    # (1.1,0,0) is in volume 3
    rets2 = [dagmc.point_in_volume(x, [1.1, 0, 0]) for x in range(1, 5)]

    return [rets1, rets2]

def find_volume():
    from pyne import dagmc
    dagmc.load(path)

    vol1 = dagmc.find_volume([0, 0, 0])
    vol2 = dagmc.find_volume([.9, .9, .9])

    # boundary case -- point [1,.1,.1] is on surface between vols 2 and 3
    # the behavior on boundaries will vary with uvw, but ensure that
    # only the two volumes the touch the location are ever returned.
    vols = []
    for uvw in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]:
        vols.append(dagmc.find_volume([1, .1, .1], uvw))

    # boundary case-- exiting volume 3 => in volume 3
    vol3 = dagmc.find_volume([1, .1, .1], [-1, 0, 0])
    vol4 = dagmc.find_volume([1.1, 0, 0])

    return [vol1, vol2, vol3, vol4, vols]

def one_ray():
    from pyne import dagmc
    dagmc.load(path)

    fromcenter = dagmc.fire_one_ray(2, [0, 0, 0], [1, 0, 0])
    fromhalf = dagmc.fire_one_ray(2, [.5, .5, .5], [0, -1, 0])

    # the following is actually a misuse of ray_fire because it should only
    # be called from inside a volume, but at least verify that a ray is
    # lost if it 1) starts outside a volume and 2) does not intersect the
    # volume
    fromoutside = dagmc.fire_one_ray(2, [0, 1.1, 0], [-1, 0, 0])

    fromvol3a = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, -1, 0])
    fromvol3b = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, 1, 0])

    return [fromcenter, fromhalf, fromoutside, fromvol3a, fromvol3b]

def failures():
    from pyne import dagmc
    dagmc.load(path)

    # if Exception is raised, then None is returned, else nothing is returned
    # and test will fail
    i = assert_raises(Exception, dagmc.point_in_volume, [100, (0, 0, 0)])
    j = assert_raises(Exception, dagmc.point_in_volume, [1, (0, 0, 0, 0)])
    k = assert_raises(Exception, dagmc.fire_one_ray, [2, (0, 0, 0), 1])

    return [i, j, k]

def ray_iterator():
    from pyne import dagmc
    dagmc.load(path)

    start = [-2, 0, 0]
    startvol = dagmc.find_volume(start)
    direction = [1, 0, 0]

    expected_vols = [2, 3, 1, 4]
    expected_dists = [1, 2, 3.156921938, 0]

    vols1 = []
    dists1 = []
    surfs1 = []
    for i, (vol, dist, surf) in enumerate(
            dagmc.ray_iterator(startvol, start, direction)):
        vols1.append(vol)
        dists1.append(dist)
        surfs1.append(surf)

    vols2 = []
    dists2 = []
    surfs2 = []
    for j, (vol, dist, surf) in enumerate(
            dagmc.ray_iterator(startvol, start, direction, dist_limit=4)):
        vols2.append(vol)
        dists2.append(dist)
        surfs2.append(surf)

    return [startvol, vols1, dists1, surfs1, i, vols2, dists2, surfs2, j]

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
    return [lo, hi]

def util_matlist():
    from pyne import dagmc
    dagmc.load(path)

    mats1 = dagmc.get_material_set()
    mats2 = dagmc.get_material_set(with_rho=True)

    return [mats1, mats2]

def discretize_geom_rand():
    from pyne import dagmc
    dagmc.load(path)

    coords = [-4, -1, 1, 4]
    mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
    results = dagmc.discretize_geom(mesh, num_rays=50)

    return [results, coords]

def discretize_geom_grid():
    from pyne import dagmc
    dagmc.load(path)

    coords = [-4, -1, 1, 4]
    mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
    results = dagmc.discretize_geom(mesh, num_rays=49, grid=True)

    return [results, coords]

def discretize_geom_mix():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    coords2 = [0, 2]
    mesh = Mesh(structured=True,
                structured_coords=[coords2, coords, coords])
    results1 = dagmc.discretize_geom(mesh, num_rays=100, grid=True)

    # To to make sure standard error decreases with increasing rays
    results2 = dagmc.discretize_geom(mesh, num_rays=625, grid=True)

    return [results1, results2]

def discretize_non_square():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    mesh = Mesh(structured=True,
                structured_coords=[coords, coords, coords])
    # if assert_raises is True, then k will be None, else a fail will occur
    k = assert_raises(ValueError, dagmc.discretize_geom, mesh, num_rays=3, grid=True)

    return k

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

    #  ensure kwargs are not accepted for unstructured mesh
    # if assert_raises is True, then k will be None, else a fail will occur
    k = assert_raises(ValueError, dagmc.discretize_geom, mesh, num_rays=3, grid=True)

    return [res, k]

def cells_at_ve_centers():
    from pyne import dagmc
    dagmc.load(path)

    coords = [0, 1]
    coords2 = [0, 2, 4]
    mesh = Mesh(structured=True,
                structured_coords=[coords2, coords, coords])
    cells = dagmc.cells_at_ve_centers(mesh)

    return cells

def test__load():
    p = multiprocessing.Pool()
    results = p.apply_async(load)
    p.close()
    p.join()
    r = results.get()

def test_metadata():
    p = multiprocessing.Pool()
    results = p.apply_async(metadata)
    p.close()
    p.join()
    r = results.get()

    rets1 = r[0]
    rets2 = r[1]
    md = r[2]

    assert_equal(rets1, [True, False, False, False])
    assert_equal(rets2, [False, False, False, True])
    assert_equal(md['material'], 5)
    assert_equal(md['rho'], 0.5)
    assert_true(all(x in md for x in ['material', 'rho', 'imp']))

def test_versions():
    p = multiprocessing.Pool()
    results = p.apply_async(versions)
    p.close()
    p.join()
    returned = results.get()

    assert_equal(len(returned), 2)
    assert_true(isinstance(returned[0], STRING_TYPES))
    assert_true(isinstance(returned[1], int))

def test_list_functions():
    p = multiprocessing.Pool()
    results = p.apply_async(list_functions)
    p.close()
    p.join()
    r = results.get()
    surfs = r[0]
    vols = r[1]
    assert_equal(set(surfs), set(range(1, 19)))
    assert_equal(set(vols), set(range(1, 5)))

def test_boundary():
    p = multiprocessing.Pool()
    results = p.apply_async(boundary)
    p.close()
    p.join()
    r = results.get()
    low = r[0]
    high = r[1]
    for i in range(0, 3):
        assert_true(low[i] <= -1.0)
        assert_true(high[i] >= 1.0)

def test_pt_in_vol():
    p = multiprocessing.Pool()
    results = p.apply_async(pt_in_vol)
    p.close()
    p.join()
    r = results.get()
    rets1 = r[0]
    rets2 = r[1]
    assert_equal(rets1, [False, True, False, False])
    assert_equal(rets2, [False, False, True, False])

def test_find_volume():
    p = multiprocessing.Pool()
    results = p.apply_async(find_volume)
    p.close()
    p.join()
    r = results.get()

    vol1 = r[0]
    vol2 = r[1]
    vol3 = r[2]
    vol4 = r[3]
    vols = r[4]

    assert_equal(vol1, 2)
    assert_equal(vol2, 2)
    assert_equal(vol3, 3)
    assert_equal(vol4, 3)

    for vol in vols:
        assert_true(vol in (2, 3))

def test_one_ray():
    p = multiprocessing.Pool()
    results = p.apply_async(one_ray)
    p.close()
    p.join()
    r = results.get()

    fromcenter = r[0]
    fromhalf = r[1]
    fromoutside = r[2]
    fromvol3a = r[3]
    fromvol3b = r[4]

    assert_almost_equal(fromcenter[1], 1.0)
    assert_almost_equal(fromhalf[1], 1.5)
    assert_equal(fromoutside, None)
    assert_almost_equal(fromvol3a[1], 0.1)
    assert_almost_equal(fromvol3b[1], 3.056921938)

def test_failures():
    p = multiprocessing.Pool()
    results = p.apply_async(failures)
    p.close()
    p.join()
    r = results.get()

    i = r[0]
    j = r[1]
    k = r[2]

    assert_equal(i, None)
    assert_equal(j, None)
    assert_equal(k, None)

def test_ray_iterator():
    p = multiprocessing.Pool()
    results = p.apply_async(ray_iterator)
    p.close()
    p.join()
    r = results.get()

    startvol = r[0]
    vols1 = r[1]
    dists1 = r[2]
    surfs1 = r[3]
    i = r[4]
    vols2 = r[5]
    dists2 = r[6]
    surfs2 = r[7]
    j = r[8]

    expected_vols = [2, 3, 1, 4]
    expected_dists = [1, 2, 3.156921938, 0]

    assert_equal(startvol, 3)

    for k, vol in enumerate(vols1):
        assert_equal(expected_vols[k], vol)
        if expected_dists[k] != 0:
            assert_almost_equal(expected_dists[k], dists1[k])
    assert_equal(i, 3)

    for k, vol in enumerate(vols2):
        assert_equal(expected_vols[k], vol)
        if expected_dists[k] != 0:
            assert_almost_equal(expected_dists[k], dists2[k])
    assert_equal(j, 1)

def test_ray_story():
    p = multiprocessing.Pool()
    results = p.apply_async(ray_story)
    p.close()
    p.join()
    r = results.get()

def test_util_graveyard_bounds():
    p = multiprocessing.Pool()
    results = p.apply_async(util_graveyard_bound)
    p.close()
    p.join()
    r = results.get()

    lo = r[0]
    hi = r[1]

    grave_diam = 4.15692194
    for i in range(0, 3):
        assert_almost_equal(lo[i], -grave_diam)
        assert_almost_equal(hi[i], grave_diam)

def test_util_matlist():
    p = multiprocessing.Pool()
    results = p.apply_async(util_matlist)
    p.close()
    p.join()
    r = results.get()

    mats1 = r[0]
    mats2 = r[1]
    assert_equal(set((0, 5)), mats1)
    assert_equal(set([(0, 0.0), (5, 0.5)]), mats2)

def test_discretize_geom_rand():
    """The 14th (index 13) mesh volume element fully contains volume 2. Use
    random sampling.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    r = p.apply_async(discretize_geom_rand)
    p.close()
    p.join()
    rr = r.get()

    results = rr[0]
    coords = rr[1]

    assert_equal(len(results), (len(coords) - 1)**3)

    for res in results:
        if res['idx'] != 13:
            assert_equal(res['cell'], 3)
        else:
            assert_equal(res['cell'], 2)

        assert_almost_equal(res['vol_frac'], 1.0)

def test_discretize_geom_grid():
    """The 14th (index 13) mesh volume element fully contains volume 2. Use
    grid sampling.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    rr = p.apply_async(discretize_geom_grid)
    p.close()
    p.join()
    r = rr.get()

    results = r[0]
    coords = r[1]

    assert_equal(len(results), (len(coords) - 1)**3)

    for res in results:
        if res['idx'] != 13:
            assert_equal(res['cell'], 3)
        else:
            assert_equal(res['cell'], 2)

        assert_almost_equal(res['vol_frac'], 1.0)

def test_discretize_geom_mix():
    """Single mesh volume element that is a 50:50 split of geometry volumes
    2 and 3.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_mix)
    p.close()
    p.join()
    r = results.get()

    results1 = r[0]
    results2 = r[1]

    assert_equal(results1[0]['cell'], 2)
    assert_almost_equal(results1[0]['vol_frac'], 0.5)
    assert_equal(results1[1]['cell'], 3)
    assert_almost_equal(results1[1]['vol_frac'], 0.5)

    # To to make sure standard error decreases with increasing rays
    assert(results2[0]['rel_error'] < results1[0]['rel_error'])
    assert(results2[1]['rel_error'] < results1[1]['rel_error'])

def test_descritize_non_square():
    """Test to make sure requesting a grid with a num_rays that is not a
    perfect square raises ValueError.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    results = p.apply_async(discretize_non_square)
    r = results.get()
    assert_equal(r, None)

def test_discretize_geom_centers():
    """Test that unstructured mesh is sampled by mesh ve centers.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    results = p.apply_async(discretize_geom_centers)
    p.close()
    p.join()
    r = results.get()

    res = r[0]
    k = r[1]

    assert_array_equal(res["idx"], [0, 1])
    assert_array_equal(res["cell"], [2, 3])
    assert_array_equal(res["vol_frac"], [1.0, 1.0])
    assert_array_equal(res["rel_error"], [1.0, 1.0])
    assert_equal(k, None)

def test_cells_at_ve_centers():
    """Test that a mesh with one ve in cell 2 and one ve in cell 3 produces
    correct results.
    """
    if not HAVE_PYMOAB:
        raise SkipTest

    p = multiprocessing.Pool()
    results = p.apply_async(cells_at_ve_centers)
    p.close()
    p.join()
    cells = results.get()

    assert_array_equal(cells, [2, 3])

def cell_material_assignments():
    from pyne import dagmc
    path = os.path.join(os.path.dirname(__file__), 'files_test_dagmc',
                        'three_blocks.h5m')
    c = dagmc.cell_material_assignments(path)
    r = []
    r.append(c[1] == 'mat:m1/rho:1.0')
    r.append(c[2] == 'mat:m2')
    r.append(c[3] == 'mat:m2/rho:3.0')
    r.append(c[6] == 'mat:graveyard')
    return np.all(r)

def test_cell_material_assignments():
    """Test cell_material_assigments().
    """
    if not HAVE_PYMOAB:
        raise SkipTest
    p = multiprocessing.Pool()
    r = p.apply_async(cell_material_assignments)
    p.close()
    p.join()
    assert_true(r.get())

def cell_materials():
    from pyne import dagmc
    path = os.path.join(os.path.dirname(__file__), 'files_test_dagmc',
                        'three_blocks.h5m')
    c = dagmc.cell_materials(path)
    r = []
    r.append(c[1].comp == {10010000: 1.0})
    r.append(c[1].density == 1.0)
    r.append(c[1].metadata['name'] == 'mat:m1/rho:1.0')
    r.append(c[2].comp == {20040000: 1.0})
    r.append(c[2].density == 2.0)
    r.append(c[2].metadata['name'] == 'mat:m2')
    r.append(c[3].comp == {20040000: 1.0})
    r.append(c[3].density == 3.0)
    r.append(c[3].metadata['name'] == 'mat:m2/rho:3.0')
    r.append(c[6].comp == {})
    r.append(c[6].density == 0)
    r.append(c[6].metadata['name'] == 'void')
    return np.all(r)

def test_cell_materials():
    """Test cell_materials().
    """
    if not HAVE_PYMOAB:
        raise SkipTest
    p = multiprocessing.Pool()
    r = p.apply_async(cell_materials)
    p.close()
    p.join()
    assert_true(r.get())
