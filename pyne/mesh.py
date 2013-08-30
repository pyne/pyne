import math
import itertools
from collections import namedtuple
from collections import Iterable

from itaps import iMesh
from itaps import iBase
from itaps import iMeshExtensions

class MeshError(Exception):
    pass

class StrMeshError(Exception):
    pass

class Mesh(object):

    def __init__(self, mesh=None, mesh_file=None, structured=False, str_coords=None, str_set=None):
        """
        Unstructured mesh instantiation:
             - From iMesh instance by specifying: <mesh>
             - From mesh file by specifying: <mesh_file>

        Structured mesh instantiation:
            - From iMesh instance with exactly 1 entity set (with BOX_DIMS tag)
              by specifying <mesh> and structured = True.
            - From mesh file with exactly 1 entity set (with BOX_DIMS tag) by
              specifying <mesh_file> and structured = True.
            - From an imesh instance with multiple entity sets by specifying 
              <mesh>, <str_set>, structured=True.
            - From coordinates by specifying <str_coords>, structured=True, and 
              optional preexisting iMesh instance <mesh>
        """
        if mesh:
            self.mesh = mesh
        else: 
            self.mesh = iMesh.Mesh()

        self.structured = structured

        #Unstructured mesh cases
        if not self.structured:
            #Error if structured arguments are passed
            if str_coords or str_set:
                MeshError("Structured mesh arguments should not be present for\
                            unstructured Mesh instantiation.")

            #From imesh instance
            if mesh and not mesh_file:
                pass
            #From file
            elif mesh_file and not mesh:
                self.mesh = iMesh.Mesh().load(mesh_file)
            else:
                raise MeshError("To instantiate structured mesh object, must " \
                                 + "supply exactly 1 of the following: "\
                                 + "<mesh>, <mesh_file>.")

        #structured mesh cases
        elif self.structured:
            #From mesh or mesh_file
            if (mesh or mesh_file) and not str_coords and not str_set:
                if mesh_file:
                    self.mesh.load(mesh_file)
                try:
                    self.mesh.getTagHandle("BOX_DIMS")
                except iBase.TagNotFoundError as e:
                    print "BOX_DIMS not found on iMesh instance"
                    raise e

                count = 0
                for ent_set in self.mesh.rootSet.getEntSets():
                    try:
                        self.mesh.getTagHandle("BOX_DIMS")[ent_set]
                    except iBase.TagNotFoundError:
                        pass
                    else:
                        self.str_set = ent_set
                        count += 1

                if count == 0:
                    raise MeshError('Found no structured meshes in \
                                        file {0}'.format(mesh_file))
                elif count > 1:
                    raise MeshError("Found {0} structured meshes.".format(count)
                                 + " Instantiate individually using from_ent_set()")
                                       
            elif not mesh and not mesh_file and str_coords and not str_set:
                extents = [0, 0, 0] + [len(x) - 1 for x in str_coords]
                self.str_set = self.mesh.createStructuredMesh(
                     extents, i=str_coords[0], j=str_coords[1], k=str_coords[2],
                     create_set=True)

            #From mesh and str_set:
            elif mesh and not mesh_file and not str_coords and str_set:
                 self.str_set = str_set
                

            else:
                raise MeshError("For structured mesh instantiation, need to \
                                   supply exactly one of the following:\n \
                                   A. iMesh instance\n\
                                   B. Mesh file\n\
                                   C. Mesh coordinates\n\
                                   D. Structured entity set AND iMesh instance")

            self.dims = self.mesh.getTagHandle('BOX_DIMS')[self.str_set]
            self.vertex_dims = list(self.dims[0:3]) + [x + 1 for x in self.dims[3:6]]                           


    # A six-element tuple corresponding to the BOX_DIMS tag on the
    # structured mesh.  See the MOAB library's metadata-info.doc file.
    #extents_tuple = namedtuple('extents',
    #                           ('imin', 'jmin', 'kmin',
    #                            'imax', 'jmax', 'kmax'))


    def str_get_vertex(self, i, j, k):
        """Return the (i,j,k)'th vertex in the mesh"""
        n = _str_find_idx(self.vertex_dims, (i, j, k))
        return _str_stepIter(self.str_set.iterate(iBase.Type.vertex, iMesh.Topology.point), n)


    def str_get_hex(self, i, j, k):
        """Return the (i,j,k)'th hexahedron in the mesh"""
        n = _str_find_idx(self.dims, (i, j, k))
        return _str_stepIter(self.str_set.iterate(iBase.Type.region, iMesh.Topology.hexahedron), n)


    def str_get_hex_volume(self, i, j, k):
        """Return the volume of the (i,j,k)'th hexahedron in the mesh"""
        v = list(self.str_iterate_vertex(x=[i, i + 1],
                                 y=[j, j + 1],
                                 z=[k, k + 1]))
        coord = self.mesh.getVtxCoords(v)
        dx = coord[1][0] - coord[0][0]
        dy = coord[2][1] - coord[0][1]
        dz = coord[4][2] - coord[0][2]
        return dx * dy * dz


    def str_iterate_hex(self, order='zyx', **kw):
        """Get an iterator over the hexahedra of the mesh

        The order argument specifies the iteration order.  It must be a string
        of 1-3 letters from the set (x,y,z).  The rightmost letter is the axis
        along which the iteration will advance the most quickly.  Thus 'zyx' --
        x coordinates changing fastest, z coordinates changing least fast-- is
        the default, and is identical to the order that would be given by the
        str_set.iterate() function.

        When a dimension is absent from the order, iteration will proceed over
        only the column in the mesh that has the lowest corresonding (i/j/k)
        coordinate.  Thus, with order 'xy,' iteration proceeds over the i/j
        plane of the structured mesh with the smallest k coordinate.

        Specific slices can be specified with keyword arguments:

        Keyword args::

          x: specify one or more i-coordinates to iterate over.
          y: specify one or more j-coordinates to iterate over.
          z: specify one or more k-coordinates to iterate over.

        Examples::

          str_iterate_hex(): equivalent to iMesh iterator over hexes in mesh
          str_iterate_hex( 'xyz' ): iterate over entire mesh, with k-coordinates
                               changing fastest, i-coordinates least fast.
          str_iterate_hex( 'yz', x=3 ): Iterate over the j-k plane of the mesh
                                   whose i-coordinate is 3, with k values
                                   changing fastest.
          str_iterate_hex( 'z' ): Iterate over k-coordinates, with i=dims.imin
                             and j=dims.jmin
          str_iterate_hex( 'yxz', y=(3,4) ): Iterate over all hexes with
                                        j-coordinate = 3 or 4.  k-coordinate
                                        values change fastest, j-values least
                                        fast.
        """

        # special case: zyx order is the standard pytaps iteration order,
        # so we can save time by simply returning a pytaps iterator
        # if no kwargs were specified
        if order == 'zyx' and not kw:
            return self.str_set.iterate(iBase.Type.region,
                                       iMesh.Topology.hexahedron)

        indices, ordmap = _str_IterSetup(self.dims, order, **kw)
        return _str_Iter(indices, ordmap, self.dims, self.str_set.iterate(iBase.Type.region, iMesh.Topology.hexahedron))


    def str_iterate_vertex(self, order='zyx', **kw):
        """Get an iterator over the vertices of the mesh

        See str_iterate_hex() for an explanation of the order argument and the
        available keyword arguments.
        """

        #special case: zyx order without kw is equivalent to pytaps iterator
        if order == 'zyx' and not kw:
            return self.str_set.iterate(iBase.Type.vertex, iMesh.Topology.point)

        indices, ordmap = _str_IterSetup(self.vertex_dims, order, **kw)
        return _str_Iter(indices, ordmap, self.vertex_dims, self.str_set.iterate(iBase.Type.vertex, iMesh.Topology.point))


    def str_iterate_hex_volumes(self, order='zyx', **kw):
        """Get an iterator over the volumes of the mesh hexahedra

        See str_iterate_hex() for an explanation of the order argument and the
        available keyword arguments.
        """

        indices, _ = _str_IterSetup(self.dims, order, **kw)
        # Use an inefficient but simple approach: call str_get_hex_volume()
        # on each required i,j,k pair.  
        # A better implementation would only make one call to getVtxCoords.
        for A in itertools.product(*indices):
            # the ordmap returned from _str_IterSetup maps to kji/zyx ordering,
            # but we want ijk/xyz ordering, so create the ordmap differently.
            ordmap = [order.find(L) for L in 'xyz']
            ijk = [A[ordmap[x]] for x in range(3)]
            yield self.str_get_hex_volume(*ijk)


    def str_get_divisions(self, dim):
        """Get the mesh divisions on a given dimension

        Given a dimension 'x', 'y', or 'z', return a list of the mesh vertices
        along that dimension
        """
        if len(dim) == 1 and dim in 'xyz':
            idx = 'xyz'.find(dim)
            return [self.mesh.getVtxCoords(i)[idx]
                    for i in self.str_iterate_vertex(dim)]
        else:
            raise StrMeshError('Invalid dimension: '+str(dim))


##########################
# private helper functions
##########################

def _str_find_idx(dims, ijk):
    """Helper method fo str_get_vertex and str_get_hex

    For tuple (i,j,k), return the number N in the appropriate iterator.
    """
    dim0 = [0] * 3
    for i in xrange(0, 3):
        if (dims[i] > ijk[i] or dims[i + 3] <= ijk[i]):
            raise StrMeshError(str(ijk) + ' is out of bounds')
        dim0[i] = ijk[i] - dims[i]
    i0, j0, k0 = dim0
    n = (((dims[4] - dims[1]) * (dims[3] - dims[0]) * k0) +
         ((dims[3] - dims[0]) * j0) +
         i0)
    return n


def _str_stepIter(it, n):
    """Helper method for str_get_vertex and str_get_hex

    Return the nth item in the iterator"""
    it.step(n)
    r = it.next()
    it.reset()
    return r


def _str_IterSetup(dims, order, **kw):
    """Setup helper function for StrMesh iterator functions

    Given dims and the arguments to the iterator function, return
    a list of three lists, each being a set of desired coordinates,
    with fastest-changing coordinate in the last column),
    and the ordmap used by _str_Iter to reorder each coodinate to (i,j,k)
    """
    # a valid order has the letters 'x', 'y', and 'z'
    # in any order without duplicates
    if not (len(order) <= 3 and
            len(set(order)) == len(order) and
            all([a in 'xyz' for a in order])):
        raise StrMeshError('Invalid iteration order: ' + str(order))

    # process kw for validity
    spec = {}
    for idx, d in enumerate('xyz'):
        if d in kw:
            spec[d] = kw[d]
            if not isinstance(spec[d], Iterable):
                spec[d] = [spec[d]]
            if not all(x in range(dims[idx], dims[idx + 3])
                    for x in spec[d]):
                raise StrMeshError( \
                        'Invalid iterator kwarg: {0}={1}'.format(d, spec[d]))
            if d not in order and len(spec[d]) > 1:
                raise StrMeshError('Cannot iterate over' + str(spec[d]) +
                                   'without a proper iteration order')
        if d not in order:
            order = d + order
            spec[d] = spec.get(d, [dims[idx]])

    # get indices and ordmap
    indices = []
    for L in order:
        idx = 'xyz'.find(L)
        indices.append(spec.get(L, xrange(dims[idx], dims[idx + 3])))

    ordmap = ['zyx'.find(L) for L in order]

    return indices, ordmap


def _str_Iter(indices, ordmap, dims, it):
    """Iterate over the indices lists, yielding _str_stepIter(it) for each.
    """
    d = [0, 0, 1]
    d[1] = (dims[3] - dims[0])
    d[0] = (dims[4] - dims[1]) * d[1]
    mins = [dims[2], dims[1], dims[0]]
    offsets = ([(a - mins[ordmap[x]]) * d[ordmap[x]]
                for a in indices[x]]
                for x in range(3))
    for ioff, joff, koff in itertools.product(*offsets):
        yield _str_stepIter(it, (ioff + joff + koff))
