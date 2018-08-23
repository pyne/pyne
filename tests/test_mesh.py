
from __future__ import print_function
import os
import sys
import time
import shutil
import warnings
import itertools
from operator import itemgetter
from nose.tools import assert_true, assert_equal, assert_raises, with_setup, \
    assert_is, assert_is_instance, assert_in, assert_not_in, assert_almost_equal

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
# try:
#     from itaps import iBase, iMesh, iMeshExtensions
# except ImportError:
#     from nose.plugins.skip import SkipTest
#     raise SkipTest

from pyne.mesh import Mesh, StatMesh, MeshError, meshset_iterate

from pymoab import core, hcoord, scd, types

def try_rm_file(filename):
    return lambda: os.remove(filename) if os.path.exists(filename) else None

def gen_mesh(mats=()):
    mesh_1 = Mesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True,
                  structured_ordering='zyx', mats=mats)
    volumes1 = list(mesh_1.structured_iterate_hex("xyz"))
    flux_tag = mesh_1.mesh.tag_get_handle("flux", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
    flux_data = [1.0, 2.0, 3.0, 4.0]
    mesh_1.mesh.tag_set_data(flux_tag, volumes1, flux_data)
    return mesh_1

#############################################
#Test unstructured mesh functionality
#############################################

def test_unstructured_mesh_from_file():
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/unstr.h5m")
    sm = Mesh(mesh=filename)

def test_unstructured_mesh_from_instance():
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/unstr.h5m")
    mesh = core.Core()
    mesh.load_file(filename)
    sm = Mesh(mesh=mesh)

def test_elem_volume():
    """Test the get_elem_volume method"""
    # Test tet elements recognition and calculation
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/unstr.h5m")
    tetmesh = Mesh(mesh=filename)
    vols = list()
    for __, __, ve in tetmesh:
        vols.append(tetmesh.elem_volume(ve))

    assert_almost_equal(np.min(vols), 0.13973, places=5)
    assert_almost_equal(np.max(vols), 2.1783, places=4)
    assert_almost_equal(np.mean(vols), 0.52702, places=5)

    # Test hex elements recognition and calculation
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/grid543.h5m")
    mesh = Mesh(mesh=filename)
    vols = list()
    for __, __, ve in mesh:
        vols.append(mesh.elem_volume(ve))
    assert_almost_equal(np.mean(vols), 51.3333, places=4)

def test_ve_center():
    m = Mesh(structured=True, structured_coords=[[-1, 3, 5], [-1, 1], [-1, 1]])
    exp_centers = [(1, 0, 0), (4, 0, 0)]
    for i, mat, ve in m:
        assert_equal(m.ve_center(ve), exp_centers[i])

#############################################
#Test structured mesh functionality
#############################################

def test_structured_mesh_from_coords():
    sm = Mesh(structured_coords = [range(1,5), range(1,4), range(1,3)], \
              structured=True)
    assert_true(sm.dims == [0, 0, 0, 3, 2, 1])
    assert_array_equal(sm.structured_coords, [range(1,5), range(1,4), range(1,3)])
    assert_equal(sm.structured_ordering, 'xyz')

def test_create_by_set():
    mesh = core.Core()
    scdi = scd.ScdInterface(mesh)
    low = hcoord.HomCoord([0,0,0])
    high = hcoord.HomCoord([1,1,1])
    structured_coords = [[1,2],[1,2],[1,2]]
    coords = []
    for x in range(low[0], high[0]+1):
        for y in range(low[1], high[1]+1):
            for z in range(low[2], high[2]+1):
                coords.append(structured_coords[0][x])
                coords.append(structured_coords[1][y])
                coords.append(structured_coords[2][z])
    a = scdi.construct_box(low, high, coords).box_set()   
    sm = Mesh(mesh=mesh, structured_set=a, structured=True)
    assert_true(sm.dims == [0, 0, 0, 1, 1, 1])

        
def test_create_by_file():
    filename = os.path.join(os.path.dirname(__file__),
                            "files_mesh_test/")
    filename += "grid543.h5m"
    sm = Mesh(mesh=filename, structured=True)
    assert_true(sm.dims == [1, 11, -5, 5, 14, -3])
    # # This mesh is interesting because the i/j/k space is not numbered from
    # # zero. Check that divisions are correct

    assert_equal(sm.structured_get_divisions("x"), [float(i) for i in range(1,6)])
    assert_equal(sm.structured_get_divisions("y"), [1.0, 5.0, 10.0, 15.0] )
    assert_equal(sm.structured_get_divisions("z"), [-10.0, 2.0, 12.0] )

    assert_equal(sm.structured_coords[0], range(1,6))
    assert_equal(sm.structured_coords[1], [1.0, 5.0, 10.0, 15.0] )
    assert_equal(sm.structured_coords[2], [-10.0, 2.0, 12.0] )

    # loading a test file without structured mesh metadata should raise an
    # error
    filename2 = os.path.join(os.path.dirname(__file__),
                              "files_mesh_test/no_str_mesh.h5m")
    assert_raises(RuntimeError, Mesh, mesh=filename2,
                        structured=True)

def test_structured_get_hex():
    # mesh with valid i values 0-4, j values 0-3, k values 0-2
    sm = Mesh(structured_coords = [range(11,16), range(21,25), range(31,34)],
              structured=True)
    def check(e):
        assert_true(isinstance(e, long))
    check(sm.structured_get_hex(0, 0, 0))
    check(sm.structured_get_hex(1, 1, 1))
    check(sm.structured_get_hex(3, 0, 0))
    check(sm.structured_get_hex(3, 2, 1))

    assert_raises(MeshError, sm.structured_get_hex,-1,-1,-1)
    assert_raises(MeshError, sm.structured_get_hex, 4, 0, 0)
    assert_raises(MeshError, sm.structured_get_hex, 0, 3, 0)
    assert_raises(MeshError, sm.structured_get_hex, 0, 0, 2)

def test_structured_hex_volume():

    sm = Mesh(structured_coords = [[0,1,3], [-3,-2,0], [12,13,15]],
              structured=True)
    assert_equal(sm.structured_hex_volume(0,0,0), 1)
    assert_equal(sm.structured_hex_volume(1,0,0), 2)
    assert_equal(sm.structured_hex_volume(0,1,0), 2)
    assert_equal(sm.structured_hex_volume(1,1,0), 4)
    assert_equal(sm.structured_hex_volume(1,1,1), 8)

    ijk_all = itertools.product(*([[0,1]]*3))

    for V, ijk in itertools.izip_longest(sm.structured_iterate_hex_volumes(),
                                         ijk_all):
        assert_equal(V, sm.structured_hex_volume(*ijk))

def test_structured_get_vertex():
    # mesh with valid i values 0-4, j values 0-3, k values 0-2
    x_range = np.array(range(10,15),dtype=np.float64)
    y_range = np.array(range(21,24),dtype=np.float64)
    z_range = np.array(range(31,33),dtype=np.float64)

    sm = Mesh(structured_coords=[x_range, y_range, z_range], structured=True)

    for i,x in enumerate(x_range):
        for j,y in enumerate(y_range):
            for k,z in enumerate(z_range):
                print("{0} {1} {2}".format(i, j, k))
                vtx = sm.structured_get_vertex(i,j,k)
                vcoord = sm.mesh.get_coords(vtx)
                assert_true(all(vcoord == [x,y,z]))

def test_get_divs():
    x = [1, 2.5, 4, 6.9]
    y = [-12, -10, -.5]
    z = [100, 200]

    sm = Mesh(structured_coords=[x, y, z], structured=True)

    assert_equal(sm.structured_get_divisions("x"), x)
    assert_equal(sm.structured_get_divisions("y"), y)
    assert_equal(sm.structured_get_divisions("z"), z)

def test_iter_structured_idx():
    m = gen_mesh() # gen_mesh uses zyx order
    xyz_idx = [0, 2, 1, 3] # expected results in xyz order
    for n, i in enumerate(m.iter_structured_idx('xyz')):
        assert_equal(i, xyz_idx[n])

#############################################
#Test mesh arithmetic for Mesh and StatMesh
#############################################

class TestArithmetic():

    def arithmetic_mesh_setup(self):
        self.mesh_1 = Mesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True)
        volumes1 = list(self.mesh_1.structured_iterate_hex("xyz"))
        volumes2 = list(self.mesh_1.structured_iterate_hex("xyz"))
        flux_tag = self.mesh_1.mesh.tag_get_handle("flux", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        flux_data = [1.0, 2.0, 3.0, 4.0]
        self.mesh_1.mesh.tag_set_data(flux_tag, volumes1, flux_data)

        self.mesh_2 = Mesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True)
        volumes1 = list(self.mesh_2.structured_iterate_hex("xyz"))
        volumes2 = list(self.mesh_2.structured_iterate_hex("xyz"))
        flux_tag = self.mesh_2.mesh.tag_get_handle("flux", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        flux_data = [1.1, 2.2, 3.3, 4.4]
        self.mesh_2.mesh.tag_set_data(flux_tag, volumes1, flux_data)

    def arithmetic_statmesh_setup(self):
        self.statmesh_1 = StatMesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True)
        volumes1 = list(self.statmesh_1.structured_iterate_hex("xyz"))
        volumes2 = list(self.statmesh_1.structured_iterate_hex("xyz"))
        flux_tag = self.statmesh_1.mesh.tag_get_handle("flux", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        error_tag = self.statmesh_1.mesh.tag_get_handle("flux_rel_error", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        flux_data = [1.0, 2.0, 3.0, 4.0]
        error_data = [0.1, 0.2, 0.3, 0.4]
        self.statmesh_1.mesh.tag_set_data(flux_tag, volumes1, flux_data)
        self.statmesh_1.mesh.tag_set_data(error_tag, volumes2, error_data)

        self.statmesh_2 = StatMesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True)
        volumes1 = list(self.statmesh_2.structured_iterate_hex("xyz"))
        volumes2 = list(self.statmesh_2.structured_iterate_hex("xyz"))
        flux_tag = self.statmesh_2.mesh.tag_get_handle("flux", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        error_tag = self.statmesh_2.mesh.tag_get_handle("flux_rel_error", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing = True)
        flux_data = [1.1, 2.2, 3.3, 4.4]
        error_data = [0.1, 0.2, 0.3, 0.4]
        self.statmesh_2.mesh.tag_set_data(flux_tag, volumes1, flux_data)
        self.statmesh_2.mesh.tag_set_data(error_tag, volumes2, error_data)

    def test_add_mesh(self):
        self.arithmetic_mesh_setup()
        self.mesh_1 += self.mesh_2
        exp_res = [2.1, 4.2, 6.3, 8.4]
        obs_res = [self.mesh_1.mesh.tag_get_data(self.mesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.mesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)
 
    def test_subtract_mesh(self):
        self.arithmetic_mesh_setup()
        self.mesh_1 -= self.mesh_2
        exp_res = [-0.1, -0.2, -0.3, -0.4]
        obs_res = [self.mesh_1.mesh.tag_get_data(self.mesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.mesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)

    def test_multiply_mesh(self):
        self.arithmetic_mesh_setup()
        self.mesh_1 *= self.mesh_2
        exp_res = [1.1, 4.4, 9.9, 17.6]
        obs_res = [self.mesh_1.mesh.tag_get_data(self.mesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.mesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)

    def test_divide_mesh(self):
        self.arithmetic_mesh_setup()
        self.mesh_1 /= self.mesh_2
        exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]
        obs_res = [self.mesh_1.mesh.tag_get_data(self.mesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.mesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)

    def test_add_statmesh(self):
        self.arithmetic_statmesh_setup()
        self.statmesh_1 += self.statmesh_2
        exp_res = [2.1, 4.2, 6.3, 8.4]
        exp_err = [0.070790803558659549, 0.1415816071173191,
                   0.21237241067597862, 0.28316321423463819]
        obs_res = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        obs_err = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux_rel_error"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)
        assert_array_almost_equal(exp_err, obs_err)

    def test_subtract_statmesh(self):
        self.arithmetic_statmesh_setup()
        self.statmesh_1 -= self.statmesh_2
        exp_res = [-0.1, -0.2, -0.3, -0.4]
        exp_err = [-1.4866068747, -2.9732137495, -4.4598206242, -5.9464274989]
        obs_res = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        obs_err = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux_rel_error"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)
        assert_array_almost_equal(exp_err, obs_err)

    def test_multiply_statmesh(self):
        self.arithmetic_statmesh_setup()
        self.statmesh_1 *= self.statmesh_2
        exp_res = [1.1, 4.4, 9.9, 17.6]
        exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249,]
        obs_res = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        obs_err = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux_rel_error"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)
        assert_array_almost_equal(exp_err, obs_err)

    def test_divide_statmesh(self):
        self.arithmetic_statmesh_setup()
        self.statmesh_1 /= self.statmesh_2
        exp_res = [0.9090909091, 0.9090909091, 0.9090909091, 0.9090909091]
        exp_err = [0.1414213562, 0.2828427125, 0.4242640687, 0.5656854249]
        obs_res = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        obs_err = [self.statmesh_1.mesh.tag_get_data(self.statmesh_1.mesh.tag_get_handle("flux_rel_error"), vol, flat = True)[0]
                   for vol in self.statmesh_1.structured_iterate_hex("xyz")]
        assert_array_almost_equal(exp_res, obs_res)
        assert_array_almost_equal(exp_err, obs_err)

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

    it = meshset_iterate(sm.mesh, sm.structured_set, types.MBHEX)

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
    it = meshset_iterate(sm.mesh, sm.structured_set, types.MBVERTEX)

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

# @with_setup(None, try_rm_file('test_matlib.h5m'))
# @with_setup(None, try_rm_file('test_matlib2.h5m'))
# def test_matlib():
#     mats = {
#         0: Material({'H1': 1.0, 'K39': 1.0}, density=1.1),
#         1: Material({'H1': 0.1, 'O16': 1.0}, density=2.2),
#         2: Material({'He4': 42.0}, density=3.3),
#         3: Material({'Tm171': 171.0}, density=4.4),
#         }
#     m = gen_mesh(mats=mats)
#     for i, ve in enumerate(m.mesh.iterate(iBase.Type.region, iMesh.Topology.all)):
#         assert_is(m.mats[i], mats[i])
#         assert_equal(m.mesh.getTagHandle('idx')[ve], i)

#     m.write_hdf5('test_matlib.h5m')
#     shutil.copy('test_matlib.h5m', 'test_matlib2.h5m')
#     m2 = Mesh(mesh='test_matlib2.h5m')  # MOAB fails to flush
#     for i, mat, ve in m2:
#         assert_equal(len(mat.comp), len(mats[i].comp))
#         for key in mats[i].iterkeys():
#             assert_equal(mat.comp[key], mats[i].comp[key])
#         assert_equal(mat.density, mats[i].density)
#         assert_equal(m2.idx[i], i)

@with_setup(None, try_rm_file('test_no_matlib.h5m'))
def test_no_matlib():
    m = gen_mesh(mats=None)
    m.write_hdf5('test_no_matlib.h5m')

if __name__ == "__main__":

    
    def trace(frame, event, arg):
        print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
        return trace
    
    sys.settrace(trace)
    test_unstructured_mesh_from_file()
