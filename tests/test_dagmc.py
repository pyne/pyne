from __future__ import print_function
import unittest
import os.path
from nose.tools import assert_equal, assert_almost_equal, assert_raises

try:
    from pyne import dagmc
    from pyne.mesh import Mesh
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest


class TestDagmcWithUnitbox(unittest.TestCase):

    # use extra underscore to ensure this function is first in alpabetical
    # sorted order, because it must run before the others.
    def test__load(self):
        # FIXME laoding causes infor to be printied to stderr (or stdout).
        path = os.path.join(os.path.dirname(__file__), 'unitbox.h5m')
        dagmc.load(path)

    def test_metadata(self):

        rets = [dagmc.volume_is_graveyard(x) for x in range(1, 5)]
        self.assertEqual(rets, [True, False, False, False])

        rets = [dagmc.volume_is_implicit_complement(x) for x in range(1, 5)]
        self.assertEqual(rets, [False, False, False, True])

        md = dagmc.volume_metadata(2)
        self.assertEqual(md['material'], 5)
        self.assertEqual(md['rho'], 0.5)
        self.assertTrue(all(x in md for x in ['material', 'rho', 'imp']))

    def test_versions(self):
        returned = dagmc.versions()
        self.assertEqual(len(returned), 2)
        self.assertTrue(isinstance(returned[0], str))
        self.assertTrue(isinstance(returned[1], int))

    def test_list_functions(self):

        surfs = dagmc.get_surface_list()
        self.assertEqual(set(surfs), set(range(1, 19)))

        vols = dagmc.get_volume_list()
        self.assertEqual(set(vols), set(range(1, 5)))

    def test_boundary(self):
        low, high = dagmc.volume_boundary(2)
        for i in range(0, 3):
            self.assertTrue(low[i] <= -1.0)
            self.assertTrue(high[i] >= 1.0)

    def test_pt_in_vol(self):

        # there are 4 volumes; (0,0,.2) is in volume 2
        rets = [dagmc.point_in_volume(x, [0, 0, .2]) for x in range(1, 5)]
        self.assertEqual(rets, [False, True, False, False])

        # (1.1,0,0) is in volume 3
        rets = [dagmc.point_in_volume(x, [1.1, 0, 0]) for x in range(1, 5)]
        self.assertEqual(rets, [False, False, True, False])

    def test_find_volume(self):

        vol = dagmc.find_volume([0, 0, 0])
        self.assertEqual(vol, 2)

        vol = dagmc.find_volume([.9, .9, .9])
        self.assertEqual(vol, 2)

        # boundary case -- point [1,.1,.1] is on surface between vols 2 and 3
        # the behavior on boundaries will vary with uvw, but ensure that
        # only the two volumes the touch the location are ever returned.
        for uvw in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]:
            vol = dagmc.find_volume([1, .1, .1], uvw)
            self.assertTrue(vol in (2, 3))

        # boundary case-- exiting volume 3 => in volume 3
        vol = dagmc.find_volume([1, .1, .1], [-1, 0, 0])
        self.assertEqual(vol, 3)

        vol = dagmc.find_volume([1.1, 0, 0])
        self.assertEqual(vol, 3)

    def test_one_ray(self):

        fromcenter = dagmc.fire_one_ray(2, [0, 0, 0], [1, 0, 0])
        self.assertAlmostEqual(fromcenter[1], 1.0)

        fromhalf = dagmc.fire_one_ray(2, [.5, .5, .5], [0, -1, 0])
        self.assertAlmostEqual(fromhalf[1], 1.5)

        # the following is actually a misuse of ray_fire because it should only
        # be called from inside a volume, but at least verify that a ray is
        # lost if it 1) starts outside a volume and 2) does not intersect the
        # volume
        fromoutside = dagmc.fire_one_ray(2, [0, 1.1, 0], [-1, 0, 0])
        self.assertEqual(fromoutside, None)

        fromvol3 = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, -1, 0])
        self.assertAlmostEqual(fromvol3[1], 0.1)

        fromvol3 = dagmc.fire_one_ray(3, [0, 1.1, 0], [0, 1, 0])
        self.assertAlmostEqual(fromvol3[1], 3.056921938)

    def test_failures(self):

        self.assertRaises(Exception, dagmc.point_in_volume, [100, (0, 0, 0)])
        self.assertRaises(Exception, dagmc.point_in_volume, [1, (0, 0, 0, 0)])
        self.assertRaises(Exception, dagmc.fire_one_ray, [2, (0, 0, 0), 1])

    def test_ray_iterator(self):

        start = [-2, 0, 0]
        startvol = dagmc.find_volume(start)
        self.assertEqual(startvol, 3)
        direction = [1, 0, 0]

        expected_vols = [2, 3, 1, 4]
        expected_dists = [1, 2, 3.156921938, 0]

        for i, (vol, dist, surf) in enumerate(
                dagmc.ray_iterator(startvol, start, direction)):
            self.assertEqual(expected_vols[i], vol)
            if expected_dists[i] != 0:
                self.assertAlmostEqual(expected_dists[i], dist)
        self.assertEqual(i, 3)

        for i, (vol, dist, surf) in enumerate(
                dagmc.ray_iterator(startvol, start, direction, dist_limit=4)):
            self.assertEqual(expected_vols[i], vol)
            if expected_dists[i] != 0:
                self.assertAlmostEqual(expected_dists[i], dist)
        self.assertEqual(i, 1)

    def test_ray_story(self):
        # Run with `nosetests -s` to see output printed on stdout
        dagmc.tell_ray_story((0, 0, 0), (1, 1, 0))
        dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0))
        # tests edge behavior
        dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0), dist_limit=4)
        dagmc.tell_ray_story((-3, 0, 0), (1, 0, 0), dist_limit=4.1)

    def test_util_graveyard_bound(self):

        lo, hi = dagmc.find_graveyard_inner_box()

        grave_diam = 4.15692194
        for i in range(0, 3):
            self.assertAlmostEqual(lo[i], -grave_diam)
            self.assertAlmostEqual(hi[i], grave_diam)

    def test_util_matlist(self):

        mats = dagmc.get_material_set()
        self.assertEqual(set((0, 5)), mats)

        mats = dagmc.get_material_set(with_rho=True)
        self.assertEqual(set([(0, 0.0), (5, 0.5)]), mats)

    def test_discretize_geom_rand(self):
        """The 14th (index 13) mesh volume element fully contains volume 2. Use 
        random sampling.
        """
        coords = [-4, -1, 1, 4]
        mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
        results = dagmc.discretize_geom(mesh, 50)
        
        assert_equal(len(results), (len(coords) - 1)**3)

        for res in results:
            if res['idx'] != 13:
                assert_equal(res['cell'], 3)
            else:
                assert_equal(res['cell'], 2)

            assert_almost_equal(res['vol_frac'], 1.0)

    def test_discretize_geom_grid(self):
        """The 14th (index 13) mesh volume element fully contains volume 2. Use 
        grid sampling.
        """
        coords = [-4, -1, 1, 4]
        mesh = Mesh(structured=True, structured_coords=[coords, coords, coords])
        results = dagmc.discretize_geom(mesh, 49, grid=True)
        
        assert_equal(len(results), (len(coords) - 1)**3)

        for res in results:
            if res['idx'] != 13:
                assert_equal(res['cell'], 3)
            else:
                assert_equal(res['cell'], 2)

            assert_almost_equal(res['vol_frac'], 1.0)
 
    def test_discretize_geom_mix(self):
        """Single mesh volume element that is a 50:50 split of geometry volumes
        2 and 3.
        """
        coords = [0, 1]
        coords2 = [0, 2]
        mesh = Mesh(structured=True, 
                    structured_coords=[coords2, coords, coords])
        results1 = dagmc.discretize_geom(mesh, 100, grid=True)

        assert_equal(results1[0]['cell'], 2)
        assert_almost_equal(results1[0]['vol_frac'], 0.5)

        assert_equal(results1[1]['cell'], 3)
        assert_almost_equal(results1[1]['vol_frac'], 0.5)

        
        # To to make sure standard error decreases with increasing rays
        results2 = dagmc.discretize_geom(mesh, 625, grid=True)
        assert(results2[0]['rel_error'] < results1[0]['rel_error'])
        assert(results2[1]['rel_error'] < results1[1]['rel_error'])

        
    def test_discretize_non_square(self):
        """Test to make sure requesting a grid with a num_rays that is not a
        perfect square raises ValueError.
        """
        coords = [0, 1]
        mesh = Mesh(structured=True, 
                    structured_coords=[coords, coords, coords])
        assert_raises(ValueError, dagmc.discretize_geom, mesh, 3, grid=True)

if __name__ == "__main__":
    import nose
    nose.runmodule()
