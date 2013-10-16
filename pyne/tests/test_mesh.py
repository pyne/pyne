import os.path
import itertools
from operator import itemgetter
from nose.tools import assert_true, assert_equal, assert_raises

import numpy as np
from itaps import iBase, iMesh, iMeshExtensions
from pyne.mesh import Mesh, StatMesh, MeshError

#############################################
#Test unstructured mesh functionality
#############################################

def test_unstructured_mesh_from_file():
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/unstr.h5m")
    sm = Mesh(mesh_file=filename)

def test_unstructured_mesh_from_instance():
    filename = os.path.join(os.path.dirname(__file__), 
                            "files_mesh_test/unstr.h5m")
    mesh = iMesh.Mesh()
    mesh.load(filename)
    sm = Mesh(mesh = mesh)

#############################################
#Test structured mesh functionality
#############################################

def test_structured_mesh_from_coords():
    sm = Mesh(structured_coords = [range(1,5), range(1,4), range(1,3)], \
              structured=True)
    assert_true(all(sm.dims == [0, 0, 0, 3, 2, 1]))

def test_create_by_set():
    mesh = iMesh.Mesh()
    a = mesh.createStructuredMesh([0,0,0,1,1,1], i=[1,2], j=[1,2], k=[1,2], 
                                 create_set=True)
    sm = Mesh(mesh = mesh, structured_set=a, structured=True)
    assert_true(all(sm.dims == [0, 0, 0, 1, 1, 1]))

def test_create_by_file():
    filename = os.path.join(os.path.dirname(__file__), 
                            "files_mesh_test/grid543.h5m")
    sm = Mesh(mesh_file = filename, structured=True)
    assert_true(all(sm.dims == [1, 11, -5, 5, 14, -3]))

    # This mesh is interesting because the i/j/k space is not numbered from
    # zero. Check that divisions are correct

    assert_equal(sm.structured_get_divisions("x"), range(1,6))
    assert_equal(sm.structured_get_divisions("y"), [1.0, 5.0, 10.0, 15.0] )
    assert_equal(sm.structured_get_divisions("z"), [-10.0, 2.0, 12.0] )

    # loading a test file without structured mesh metadata should raise an 
    # error
    filename2 = os.path.join(os.path.dirname(__file__), 
                             "files_mesh_test/no_str_mesh.h5m")
    assert_raises(iBase.TagNotFoundError, Mesh, mesh_file = filename2, 
                       structured=True)

def test_structured_get_hex():
    # mesh with valid i values 0-4, j values 0-3, k values 0-2
    sm = Mesh(structured_coords = [range(11,16), range(21,25), range(31,34)], 
              structured=True)
    def check(e):
        assert_true(isinstance(e, iBase.Entity))
    check(sm.structured_get_hex(0, 0, 0))
    check(sm.structured_get_hex(1, 1, 1))
    check(sm.structured_get_hex(3, 0, 0))
    check(sm.structured_get_hex(3, 2, 1))

    assert_raises(MeshError, sm.structured_get_hex,-1,-1,-1)
    assert_raises(MeshError, sm.structured_get_hex, 4, 0, 0)
    assert_raises(MeshError, sm.structured_get_hex, 0, 3, 0)
    assert_raises(MeshError, sm.structured_get_hex, 0, 0, 2)

def test_hex_volume():

    sm = Mesh(structured_coords = [[0,1,3], [-3,-2,0], [12,13,15]], 
              structured=True)
    assert_equal(sm.structured_get_hex_volume(0,0,0), 1)
    assert_equal(sm.structured_get_hex_volume(1,0,0), 2)
    assert_equal(sm.structured_get_hex_volume(0,1,0), 2)
    assert_equal(sm.structured_get_hex_volume(1,1,0), 4)
    assert_equal(sm.structured_get_hex_volume(1,1,1), 8)

    ijk_all = itertools.product(*([[0,1]]*3))

    for V, ijk in itertools.izip_longest(sm.structured_iterate_hex_volumes(), 
                                         ijk_all):
        assert_equal(V, sm.structured_get_hex_volume(*ijk))


def test_structured_get_vertex():
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
                assert_true(all(vcoord == [x,y,z]))

def test_get_divs():
    x = [1, 2.5, 4, 6.9]
    y = [-12, -10, -.5]
    z = [100, 200]

    sm = Mesh(structured_coords=[x, y, z], structured=True)

    assert_equal(sm.structured_get_divisions("x"), x)
    assert_equal(sm.structured_get_divisions("y"), y)
    assert_equal(sm.structured_get_divisions("z"), z)

#############################################
#Test mesh arithematic for Mesh and StatMesh
#############################################
def test_add_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh1.add(mesh2)

    exp_res = [2.1, 4.2, 6.3, 8.4]

    for i, vol in enumerate(list(mesh1.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8


def test_op_add_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh3 = mesh1 + mesh2

    exp_res = [2.1, 4.2, 6.3, 8.4]

    for i, vol in enumerate(list(mesh3.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_subtract_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh1.sub(mesh2)

    exp_res = [-0.1, -0.2, -0.3, -0.4]

    for i, vol in enumerate(list(mesh1.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_op_subtract_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh3 = mesh1 - mesh2

    exp_res = [-0.1, -0.2, -0.3, -0.4]

    for i, vol in enumerate(list(mesh3.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_multiply_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh1.mul(mesh2)

    exp_res = [1.1, 4.4, 9.9, 17.6]

    for i, vol in enumerate(list(mesh1.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_op_multiply_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh3 = mesh1 * mesh2

    exp_res = [1.1, 4.4, 9.9, 17.6]

    for i, vol in enumerate(list(mesh3.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_divide_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh1.div(mesh2)

    exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]

    for i, vol in enumerate(list(mesh1.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_op_divide_mesh():
    mesh1 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    mesh2 = Mesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    mesh3 = mesh1/mesh2

    exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]

    for i, vol in enumerate(list(mesh3.structured_iterate_hex("xyz"))):
        assert abs(mesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8

def test_add_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh1.add(statmesh2)

    exp_res = [2.1, 4.2, 6.3, 8.4]
    exp_err = [0.070790803558659549, 0.1415816071173191, 0.21237241067597862, 0.28316321423463819]

    for i, vol in enumerate(list(statmesh1.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol] - exp_err[i]) < 1E-8

def test_op_add_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh3 = statmesh1 + statmesh2

    exp_res = [2.1, 4.2, 6.3, 8.4]
    exp_err = [0.070790803558659549, 0.1415816071173191, 0.21237241067597862, 0.28316321423463819]

    for i, vol in enumerate(list(statmesh3.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol] - exp_err[i]) <1E-8

def test_subtract_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh1.sub(statmesh2)

    exp_res = [-0.1, -0.2, -0.3, -0.4]
    exp_err = [-1.4866068747, -2.9732137495, -4.4598206242, -5.9464274989]

    for i, vol in enumerate(list(statmesh1.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol]- exp_err[i]) < 1E-8

def test_op_subtract_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh3 = statmesh1 - statmesh2

    exp_res = [-0.1, -0.2, -0.3, -0.4]
    exp_err = [-1.4866068747, -2.9732137495, -4.4598206242, -5.9464274989]

    for i, vol in enumerate(list(statmesh3.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol] - exp_err[i]) < 1E-8

def test_multiply_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh1.mul(statmesh2)

    exp_res = [1.1, 4.4, 9.9, 17.6]
    exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249,]

    for i, vol in enumerate(list(statmesh1.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol]- exp_err[i]) < 1E-8

def test_op_multiply_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh3 = statmesh1 * statmesh2

    exp_res = [1.1, 4.4, 9.9, 17.6]
    exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249,]

    for i, vol in enumerate(list(statmesh3.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol] - exp_err[i]) < 1E-8

def test_divide_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh1.div(statmesh2)

    exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]
    exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249]

    for i, vol in enumerate(list(statmesh1.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol]- exp_err[i]) < 1E-8

def test_op_divide_statmesh():
    statmesh1 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_1.h5m", 
                         structured=True)
    statmesh2 = StatMesh(mesh_file="files_mesh_test/test_mesh_2x2x1_2.h5m", 
                         structured=True)
    statmesh3 = statmesh1/statmesh2

    exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]
    exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249]

    for i, vol in enumerate(list(statmesh3.structured_iterate_hex("xyz"))):
        assert abs(statmesh1.mesh.getTagHandle("flux")[vol] - exp_res[i]) < 1E-8
        assert abs(statmesh1.mesh.getTagHandle("flux_error")[vol] - exp_err[i]) < 1E-8


#############################################
#Test Structured mesh iteration functionality
#############################################

def test_bad_iterates():
    sm = Mesh(structured=True,
             structured_coords =[range(10,15), range(21,25), range(31,34)])

    assert_raises(MeshError, sm.structured_iterate_hex, "abc")
    assert_raises(TypeError, sm.structured_iterate_hex, 12)
    assert_raises(MeshError, sm.structured_iterate_hex, "xxyz")
    assert_raises(MeshError, sm.structured_iterate_hex, "yyx")
    assert_raises(MeshError, sm.structured_iterate_hex, "xyz", z=[0,1,2])

def test_iterate_3d():        
    # use izip_longest in the lockstep iterations below; this will catch any
    # situations where one iterator turns out to be longer than expected.
    sm = Mesh(structured=True,
             structured_coords =[range(10,15), range(21,25), range(31,34)])
    I = range(0,4)
    J = range(0,3)
    K = range(0,2)
    izip = itertools.izip_longest

    it = sm.structured_set.iterate(iBase.Type.region, 
                                 iMesh.Topology.hexahedron)

    # Test the zyx order, which is default; it should be equivalent
    # to the standard imesh iterator
    for it_x, sm_x in izip(it, sm.structured_iterate_hex()):
        assert_equal(it_x, sm_x)

    #testing xyz

    all_indices_zyx = itertools.product(I, J, K)
    # Test the xyz order, the default from original mmGridGen
    for ijk_index, sm_x in izip(all_indices_zyx, 
                                 sm.structured_iterate_hex("xyz")):
        assert_equal(sm.structured_get_hex(*ijk_index), sm_x )

    def _tuple_sort(collection, indices ):
        # sorting function for order test
        def t(tup):
            # sort this 3-tuple according to the order of x, y, and z in 
            #indices
            return (tup["xyz".find(indices[0])]*100 +
                    tup["xyz".find(indices[1])]*10 +
                    tup["xyz".find(indices[2])])
        return sorted(collection, key = t)

    def test_order(order, *args,  **kw):
        all_indices = itertools.product(*args)
        for ijk_index, sm_x in izip(_tuple_sort(all_indices, order),
                                     sm.structured_iterate_hex(order,**kw)):
            assert_equal(sm.structured_get_hex(*ijk_index), sm_x)

    test_order("yxz", I, J, K)
    test_order("yzx", I, J, K)
    test_order("xzy", I, J, K)
    test_order("zxy", I, J, K)

    # Specify z=[1] to iterator
    test_order("xyz", I, J, [1], z=[1])
    # Specify y=2 to iterator
    test_order("zyx", I, [2], K, y=2)
    # specify x and y both to iterator
    test_order("yzx", [1,2,3],J[:-1], K, y=J[:-1], x=[1,2,3])

def test_iterate_2d():
    sm = Mesh(structured=True,
             structured_coords =[range(10,15), range(21,25), range(31,34)])
    def test_order(iter1, iter2):
        for i1, i2 in itertools.izip_longest(iter1, iter2):
            assert_equal(i1, i2)

    test_order(sm.structured_iterate_hex("yx"), 
               sm.structured_iterate_hex("zyx", z=[0]))
    test_order(sm.structured_iterate_hex("yx",z=1), 
               sm.structured_iterate_hex("zyx",z=[1]))
    test_order(sm.structured_iterate_hex("yx",z=1), 
               sm.structured_iterate_hex("yzx",z=[1]))
    test_order(sm.structured_iterate_hex("zy",x=[3]), 
               sm.structured_iterate_hex("zxy",x=3))

    # Cannot iterate over multiple z's without specifing z order
    assert_raises(MeshError, sm.structured_iterate_hex, "yx", z=[0,1])

def test_iterate_1d():
    sm = Mesh(structured=True,
             structured_coords =[range(10,15), range(21,25), range(31,34)])
    
    def test_equal(ijk_list, miter):
        for ijk, i in itertools.izip_longest(ijk_list, miter):
            assert_equal(sm.structured_get_hex(*ijk), i)

    test_equal([[0,0,0],[0,0,1]], 
                sm.structured_iterate_hex("z"))

    test_equal([[0,1,1],[0,2,1]],
                sm.structured_iterate_hex("y", y=[1,2], z=1))

    test_equal([[2,0,0],[2,1,0],[2,2,0]],
                sm.structured_iterate_hex("y", x=2)) 
    test_equal([[0,0,0],[1,0,0],[2,0,0]], 
        sm.structured_iterate_hex("x", x=[0,1,2]))

def test_vtx_iterator():
    #use vanilla izip as we"ll test using non-equal-length iterators
    izip = itertools.izip

    sm = Mesh(structured=True,
             structured_coords =[range(10,15), range(21,25), range(31,34)])
    it = sm.structured_set.iterate(iBase.Type.vertex, iMesh.Topology.point)

    # test the default order
    for (it_x, sm_x) in itertools.izip_longest(it, 
                                           sm.structured_iterate_vertex("zyx")):
        assert_equal(it_x,sm_x)

    # Do the same again, but use an arbitrary kwarg to structured_iterate_vertex 
    # to prevent optimization from kicking in
    it.reset()
    for (it_x, sm_x) in itertools.izip_longest(it,
                            sm.structured_iterate_vertex("zyx", no_opt=True)):
        assert_equal(it_x,sm_x)

    it.reset()
    for (it_x, sm_x) in izip(it, 
                               sm.structured_iterate_vertex("yx",z=sm.dims[2])):
        assert_equal(it_x,sm_x)

    it.reset()
    for (it_x, sm_x) in izip(it, sm.structured_iterate_vertex("x")):
        assert_equal(it_x,sm_x)


def test_large_iterator():
    #Test performance with large mesh
    print "building large mesh"
    big = Mesh(structured_coords=[range(1,100), range(101,200), range(201,300)], 
               structured = True)
    print "iterating (1)"
    for i in big.structured_iterate_hex():
        pass
    print "iterating (2)"
    for i in big.structured_iterate_hex("yzx"):
        pass
