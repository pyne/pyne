from itaps import iBase, iMesh, iMeshExtensions

import numpy as np
import unittest
import os.path
import itertools
from operator import itemgetter

from pyne.mesh import Mesh, MeshError

class MeshTest(unittest.TestCase):

    def setUp(self):
        self.mesh = iMesh.Mesh()

    def test_unstructured_mesh_from_file(self):
        filename = os.path.join(os.path.dirname(__file__), 
                                "files_mesh_test/unstr.h5m")
        sm = Mesh(mesh_file=filename)

    def test_unstructured_mesh_from_instance(self):
        filename = os.path.join(os.path.dirname(__file__), 
                                "files_mesh_test/unstr.h5m")
        hello = iMesh.Mesh()
        self.mesh.load(filename)
        sm = Mesh(mesh = self.mesh)

    def test_structured_mesh_from_coords(self):
        sm = Mesh(structured_coords = [range(1,5), range(1,4), range(1,3)], \
                  structured=True)
        self.assertTrue(all(sm.dims == [0, 0, 0, 3, 2, 1]))

    def test_create_by_set(self):
        a =  self.mesh.createStructuredMesh( [0,0,0,1,1,1], 
                                             i=[1,2], j=[1,2], k=[1,2], 
                                             create_set=True )
        sm = Mesh(mesh = self.mesh, structured_set=a, structured=True)
        self.assertTrue(all(sm.dims == [0, 0, 0, 1, 1, 1]))

    def test_create_by_file(self):
        filename = os.path.join(os.path.dirname(__file__), 
                                "files_mesh_test/grid543.h5m")
        sm = Mesh(mesh_file = filename, structured=True)
        self.assertTrue(all(sm.dims == [1, 11, -5, 5, 14, -3]))

        # This mesh is interesting because the i/j/k space is not numbered from
        # zero. Check that divisions are correct

        self.assertEqual(sm.structured_get_divisions("x"), range(1,6) )
        self.assertEqual(sm.structured_get_divisions("y"), [1.0, 5.0, 10.0, 15.0] )
        self.assertEqual(sm.structured_get_divisions("z"), [-10.0, 2.0, 12.0] )

        # loading a test file without structured mesh metadata should raise an 
        # error
        filename2 = os.path.join(os.path.dirname(__file__), 
                                 "files_mesh_test/no_str_mesh.h5m")
        self.assertRaises(iBase.TagNotFoundError, Mesh, mesh_file = filename2, 
                           structured=True)

    def test_structured_get_hex(self):
        # mesh with valid i values 0-4, j values 0-3, k values 0-2
        sm = Mesh(structured_coords = [range(11,16), range(21,25), range(31,34)], 
                  structured=True)
        def check(e):
            self.assertTrue(isinstance(e, iBase.Entity))
        check(sm.structured_get_hex(0, 0, 0))
        check(sm.structured_get_hex(1, 1, 1))
        check(sm.structured_get_hex(3, 0, 0))
        check(sm.structured_get_hex(3, 2, 1))

        self.assertRaises(MeshError, sm.structured_get_hex,-1,-1,-1)
        self.assertRaises(MeshError, sm.structured_get_hex, 4, 0, 0)
        self.assertRaises(MeshError, sm.structured_get_hex, 0, 3, 0)
        self.assertRaises(MeshError, sm.structured_get_hex, 0, 0, 2)

    def test_hex_volume(self):

        sm = Mesh(structured_coords = [[0,1,3], [-3,-2,0], [12,13,15]], 
                  structured=True)
        self.assertEqual(sm.structured_get_hex_volume(0,0,0), 1)
        self.assertEqual(sm.structured_get_hex_volume(1,0,0), 2)
        self.assertEqual(sm.structured_get_hex_volume(0,1,0), 2)
        self.assertEqual(sm.structured_get_hex_volume(1,1,0), 4)
        self.assertEqual(sm.structured_get_hex_volume(1,1,1), 8)

        ijk_all = itertools.product(*([[0,1]]*3))

        for V, ijk in itertools.izip_longest(sm.structured_iterate_hex_volumes(), 
                                             ijk_all):
            self.assertEqual(V, sm.structured_get_hex_volume(*ijk))


    def test_structured_get_vertex(self):
        # mesh with valid i values 0-4, j values 0-3, k values 0-2
        x_range = np.array(range(10,15),dtype=np.float64)
        y_range = np.array(range(21,24),dtype=np.float64)
        z_range = np.array(range(31,33),dtype=np.float64)

        sm = Mesh(structured_coords=[x_range, y_range, z_range], structured=True)

        for i,x in enumerate(x_range):
            for j,y in enumerate(y_range):
                for k,z in enumerate(z_range):
                    print i, j, k
                    vtx = sm.structured_get_vertex(i,j,k)
                    vcoord = sm.mesh.getVtxCoords(vtx)
                    self.assertTrue(all(vcoord == [x,y,z]))

    def test_get_divs(self):
        x = [1, 2.5, 4, 6.9]
        y = [-12, -10, -.5]
        z = [100, 200]

        sm = Mesh(structured_coords=[x, y, z], structured=True)

        self.assertEqual(sm.structured_get_divisions("x"), x)
        self.assertEqual(sm.structured_get_divisions("y"), y)
        self.assertEqual(sm.structured_get_divisions("z"), z)

class StrMeshIterateTest(unittest.TestCase):

    def setUp(self):
        self.mesh = iMesh.Mesh()
        #i = 0,1,2,3, j = 0,1,2, k = 0,1
        self.sm = Mesh(structured=True,\
                       structured_coords =[range(10,15), range(21,25), range(31,34)])   
        self.I = range(0,4)
        self.J = range(0,3)
        self.K = range(0,2)

    def test_bad_iterates(self):

        self.assertRaises(MeshError, self.sm.structured_iterate_hex, "abc")
        self.assertRaises(TypeError, self.sm.structured_iterate_hex, 12)
        self.assertRaises(MeshError, self.sm.structured_iterate_hex, "xxyz")
        self.assertRaises(MeshError, self.sm.structured_iterate_hex, "yyx")
        self.assertRaises(MeshError, self.sm.structured_iterate_hex, "xyz", z=[0,1,2])

    def test_iterate_3d(self):        
        # use izip_longest in the lockstep iterations below; this will catch any
        # situations where one iterator turns out to be longer than expected.
        izip = itertools.izip_longest

        it = self.sm.structured_set.iterate( iBase.Type.region, 
                                     iMesh.Topology.hexahedron )

        print "testing zyx"

        # Test the zyx order, which is default; it should be equivalent
        # to the standard imesh iterator
        for it_x, sm_x in izip( it, self.sm.structured_iterate_hex() ):
            self.assertEqual( it_x, sm_x )

        print "testing xyz"

        all_indices_zyx = itertools.product( self.I, self.J, self.K )
        # Test the xyz order, the default from original mmGridGen
        for ijk_index, sm_x in izip( all_indices_zyx, 
                                     self.sm.structured_iterate_hex("xyz") ):
            self.assertEqual( self.sm.structured_get_hex(*ijk_index), sm_x )

        def tuple_sort( collection, indices ):
            # sorting function for order test
            def t( tup ):
                # sort this 3-tuple according to the order of x, y, and z in 
                #indices
                return ( tup["xyz".find(indices[0])]*100 +
                         tup["xyz".find(indices[1])]*10 +
                         tup["xyz".find(indices[2])] )
            return sorted( collection, key = t )

        def test_order( order, *args,  **kw ):
            print "testing",order
            all_indices = itertools.product(*args)
            for ijk_index, sm_x in izip( tuple_sort(all_indices, order),
                                         self.sm.structured_iterate_hex(order,**kw) ):
                self.assertEqual(self.sm.structured_get_hex(*ijk_index), sm_x)

        test_order( "yxz", self.I, self.J, self.K )
        test_order( "yzx", self.I, self.J, self.K )
        test_order( "xzy", self.I, self.J, self.K )
        test_order( "zxy", self.I, self.J, self.K )

        # Specify z=[1] to iterator
        test_order("xyz", self.I, self.J, [1], z=[1])
        # Specify y=2 to iterator
        test_order("zyx", self.I, [2], self.K, y=2)
        # specify x and y both to iterator
        test_order("yzx", [1,2,3],self.J[:-1], self.K, y=self.J[:-1], x=[1,2,3])

    def test_iterate_2d(self):
        def test_order(iter1, iter2):
            for i1, i2 in itertools.izip_longest(iter1, iter2):
                self.assertEqual( i1, i2 )

        test_order(self.sm.structured_iterate_hex("yx"), 
                   self.sm.structured_iterate_hex("zyx", z=[0]))
        test_order(self.sm.structured_iterate_hex("yx",z=1), 
                   self.sm.structured_iterate_hex("zyx",z=[1]))
        test_order(self.sm.structured_iterate_hex("yx",z=1), 
                   self.sm.structured_iterate_hex("yzx",z=[1]))
        test_order(self.sm.structured_iterate_hex("zy",x=[3]), 
                   self.sm.structured_iterate_hex("zxy",x=3))

        # Cannot iterate over multiple z's without specifing z order
        self.assertRaises(MeshError, self.sm.structured_iterate_hex, "yx", z=[0,1])

    def test_iterate_1d(self):
        
        def test_equal( ijk_list, miter ):
            for ijk, i in itertools.izip_longest( ijk_list, miter ):
                self.assertEqual( self.sm.structured_get_hex(*ijk), i )

        test_equal( [[0,0,0],[0,0,1]], 
                    self.sm.structured_iterate_hex("z") )

        test_equal( [[0,1,1],[0,2,1]],
                    self.sm.structured_iterate_hex("y", y=[1,2], z=1) )

        test_equal( [[2,0,0],[2,1,0],[2,2,0]],
                    self.sm.structured_iterate_hex("y", x=2) ) 
        test_equal( [[0,0,0],[1,0,0],[2,0,0]], 
            self.sm.structured_iterate_hex("x", x=[0,1,2]) )

    def test_vtx_iterator(self):
        
        #use vanilla izip as we"ll test using non-equal-length iterators
        izip = itertools.izip

        sm = self.sm
        it = sm.structured_set.iterate(iBase.Type.vertex, iMesh.Topology.point)

        # test the default order
        for (it_x, sm_x) in itertools.izip_longest(it, 
                                                  sm.structured_iterate_vertex("zyx")):
            self.assertEqual(it_x,sm_x)

        # Do the same again, but use an arbitrary kwarg to structured_iterate_vertex to
        # prevent optimization from kicking in
        it.reset()
        for (it_x, sm_x) in itertools.izip_longest(it,
                                    sm.structured_iterate_vertex("zyx", no_opt=True) ):
            self.assertEqual(it_x,sm_x)

        it.reset()
        for (it_x, sm_x) in izip(it, sm.structured_iterate_vertex("yx",z=sm.dims[2])):
            self.assertEqual(it_x,sm_x)

        it.reset()
        for (it_x, sm_x) in izip( it, sm.structured_iterate_vertex("x")):
            self.assertEqual(it_x,sm_x)

    
class StrPerfTest(unittest.TestCase):

    # Give this class the perf attribute; this slow test may be skipped
    # using nosetests -a '!perf'
    perf = True

    def test_large_iterator(self):

        print "building large mesh"
        big = Mesh(structured_coords=[range(1,100), range(101,200), range(201,300)], 
                   structured = True)
        print "iterating (1)"
        for i in big.structured_iterate_hex():
            pass
        print "iterating (2)"
        for i in big.structured_iterate_hex( "yzx" ):
            pass

