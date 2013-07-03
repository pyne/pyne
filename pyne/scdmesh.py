import itertools
from collections import namedtuple, Iterable

from itaps import iBase, iMesh, iMeshExtensions


class ScdMeshError(Exception):
    pass


class ScdMesh:
    """A structured mesh in the spirit of MOAB's ScdMesh interface.
    
    Public member variables::
        
        self.imesh -- the iMesh instance in which this 
                      structured mesh resides
        self.dims -- A namedtuple indicating the minimum and maximum
                     (i,j,k) coordinates of this structured mesh.
    """


    # A six-element tuple corresponding to the BOX_DIMS tag on the
    # structured mesh.  See the MOAB library's metadata-info.doc file.
    extents_tuple = namedtuple('extents',
                               ('imin', 'jmin', 'kmin',
                                'imax', 'jmax', 'kmax'))

    def __init__(self, x_points, y_points, z_points, imesh=None, **kw):
        """Construct a ScdMesh from given x, y, and z coordinates.

        Parameters
        ----------
        x_points, y_points, z_points : list of floats
            List of points.
            The i,j,k extents of the structured mesh will be numbered from 0.
        imesh : iMesh.Mesh object
            Parameter can be used to designate an iMesh instance
            in which to place the mesh; if None, a new instance will be 
            created.
        keyword arguments:
            bdtag : iMesh.Tag handle
                if the BOX_DIMS tag has already been looked up, it may
                be passed thus
            _scdset : iMesh.EntitySet handle (?)
                an existing scdset to use; clients should use 
                fromEntSet() instead of using this parameter.
        """
        if imesh:
            self.imesh = imesh
        else:
            self.imesh = iMesh.Mesh()
        if x_points is None and y_points is None and z_points is None:
            self.scdset = kw['_scdset']
        else:
            extents = [0, 0, 0]
            extents.extend([len(x)-1 for x in [x_points, y_points, z_points]])
            self.scdset = self.imesh.createStructuredMesh(
                    extents, i=x_points, j=y_points, k=z_points,
                    create_set=True)
        bdtag = kw.get('bdtag', self.imesh.getTagHandle('BOX_DIMS'))
        self.dims = ScdMesh.extents_tuple(*bdtag[self.scdset])
        vdims_incr = list(self.dims[0:3]) + [x + 1 for x in self.dims[3:6]]
        self.vdims = ScdMesh.extents_tuple(*vdims_incr)
        # access iterators for hexes and vertices; used by getHex()/getVtx()
        self.hexit = self.scdset.iterate(iBase.Type.region,
                                         iMesh.Topology.hexahedron)
        self.vtxit = self.scdset.iterate(iBase.Type.vertex,
                                         iMesh.Topology.point)

    @classmethod
    def fromFile(cls, filename, imesh=None):
        """Load structured meshes from a file

        Returns one strutured mesh if the file contains one mesh, or a list if
        the file contains multiple meshes.
        """
        retlist = []
        if not imesh:
            imesh = iMesh.Mesh()
        imesh.rootSet.load(filename)
        for eset in imesh.rootSet.getEntSets():
            try:
                bdtag = imesh.getTagHandle('BOX_DIMS')
                bdtag[eset]
                retlist.append(cls(None, None, None, imesh,
                                   _scdset=eset, bdtag=bdtag))
            except iBase.TagNotFoundError:
                pass
        if not retlist:
            raise ScdMeshError('Found no structured meshes in file '+filename)
        if len(retlist) == 1:
            return retlist[0]
        else:
            return retlist

    @classmethod
    def fromEntSet(cls, imesh, eset):
        """Constructor function: create a ScdMesh from an existing entity set

        The eset parameter must be a structured mesh set with the BOX_DIMS
        tag set on it.
        """
        m = cls(None, None, None, imesh, _scdset=eset)
        return m

    def getVtx(self, i, j, k):
        """Return the (i,j,k)'th vertex in the mesh"""
        n = _dimConvert(self.vdims, (i, j, k))
        return _stepIter(self.vtxit, n)

    def getHex(self, i, j, k):
        """Return the (i,j,k)'th hexahedron in the mesh"""
        n = _dimConvert(self.dims, (i, j, k))
        return _stepIter(self.hexit, n)

    def getHexVolume(self, i, j, k):
        """Return the volume of the (i,j,k)'th hexahedron in the mesh"""
        v = list(self.iterateVtx(x=[i, i + 1],
                                 y=[j, j + 1],
                                 z=[k, k + 1]))
        coord = self.imesh.getVtxCoords(v)
        dx = coord[1][0] - coord[0][0]
        dy = coord[2][1] - coord[0][1]
        dz = coord[4][2] - coord[0][2]
        return dx * dy * dz

    def iterateHex(self, order='zyx', **kw):
        """Get an iterator over the hexahedra of the mesh

        The order argument specifies the iteration order.  It must be a string
        of 1-3 letters from the set (x,y,z).  The rightmost letter is the axis
        along which the iteration will advance the most quickly.  Thus 'zyx' --
        x coordinates changing fastest, z coordinates changing least fast-- is
        the default, and is identical to the order that would be given by the
        scdset.iterate() function.

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

          iterateHex(): equivalent to iMesh iterator over hexes in mesh
          iterateHex( 'xyz' ): iterate over entire mesh, with k-coordinates
                               changing fastest, i-coordinates least fast.
          iterateHex( 'yz', x=3 ): Iterate over the j-k plane of the mesh
                                   whose i-coordinate is 3, with k values
                                   changing fastest.
          iterateHex( 'z' ): Iterate over k-coordinates, with i=dims.imin
                             and j=dims.jmin
          iterateHex( 'yxz', y=(3,4) ): Iterate over all hexes with
                                        j-coordinate = 3 or 4.  k-coordinate
                                        values change fastest, j-values least
                                        fast.
        """

        # special case: zyx order is the standard pytaps iteration order,
        # so we can save time by simply returning a pytaps iterator
        # if no kwargs were specified
        if order == 'zyx' and not kw:
            return self.scdset.iterate(iBase.Type.region,
                                       iMesh.Topology.hexahedron)

        indices, ordmap = _scdIterSetup(self.dims, order, **kw)
        return _scdIter(indices, ordmap, self.dims, self.hexit)

    def iterateVtx(self, order='zyx', **kw):
        """Get an iterator over the vertices of the mesh

        See iterateHex() for an explanation of the order argument and the
        available keyword arguments.
        """

        #special case: zyx order without kw is equivalent to pytaps iterator
        if order == 'zyx' and not kw:
            return self.scdset.iterate(iBase.Type.vertex, iMesh.Topology.point)

        indices, ordmap = _scdIterSetup(self.vdims, order, **kw)
        return _scdIter(indices, ordmap, self.vdims, self.vtxit)

    def iterateHexVolumes(self, order='zyx', **kw):
        """Get an iterator over the volumes of the mesh hexahedra

        See iterateHex() for an explanation of the order argument and the
        available keyword arguments.
        """

        indices, _ = _scdIterSetup(self.dims, order, **kw)
        # Use an inefficient but simple approach: call getHexVolume()
        # on each required i,j,k pair.  
        # A better implementation would only make one call to getVtxCoords.
        for A in itertools.product(*indices):
            # the ordmap returned from _scdIterSetup maps to kji/zyx ordering,
            # but we want ijk/xyz ordering, so create the ordmap differently.
            ordmap = [order.find(L) for L in 'xyz']
            ijk = [A[ordmap[x]] for x in range(3)]
            yield self.getHexVolume(*ijk)

    def getDivisions(self, dim):
        """Get the mesh divisions on a given dimension

        Given a dimension 'x', 'y', or 'z', return a list of the mesh vertices
        along that dimension
        """
        if len(dim) == 1 and dim in 'xyz':
            idx = 'xyz'.find(dim)
            return [self.imesh.getVtxCoords(i)[idx]
                    for i in self.iterateVtx(dim)]
        else:
            raise ScdMeshError('Invalid dimension: '+str(dim))


def _dimConvert(dims, ijk):
    """Helper method fo getVtx and getHex

    For tuple (i,j,k), return the number N in the appropriate iterator.
    """
    dim0 = [0] * 3
    for i in xrange(0, 3):
        if (dims[i] > ijk[i] or dims[i + 3] <= ijk[i]):
            raise ScdMeshError(str(ijk) + ' is out of bounds')
        dim0[i] = ijk[i] - dims[i]
    i0, j0, k0 = dim0
    n = (((dims[4] - dims[1]) * (dims[3] - dims[0]) * k0) +
         ((dims[3] - dims[0]) * j0) +
         i0)
    return n


def _stepIter(it, n):
    """Helper method for getVtx and getHex

    Return the nth item in the iterator"""
    it.step(n)
    r = it.next()
    it.reset()
    return r


def _scdIterSetup(dims, order, **kw):
    """Setup helper function for ScdMesh iterator functions

    Given dims and the arguments to the iterator function, return
    a list of three lists, each being a set of desired coordinates,
    with fastest-changing coordinate in the last column),
    and the ordmap used by _scdIter to reorder each coodinate to (i,j,k)
    """
    # a valid order has the letters 'x', 'y', and 'z'
    # in any order without duplicates
    if not (len(order) <= 3 and
            len(set(order)) == len(order) and
            all([a in 'xyz' for a in order])):
        raise ScdMeshError('Invalid iteration order: ' + str(order))

    # process kw for validity
    spec = {}
    for idx, d in enumerate('xyz'):
        if d in kw:
            spec[d] = kw[d]
            if not isinstance(spec[d], Iterable):
                spec[d] = [spec[d]]
            if not all(x in range(dims[idx], dims[idx + 3])
                    for x in spec[d]):
                raise ScdMeshError( \
                        'Invalid iterator kwarg: {0}={1}'.format(d, spec[d]))
            if d not in order and len(spec[d]) > 1:
                raise ScdMeshError('Cannot iterate over' + str(spec[d]) +
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


def _scdIter(indices, ordmap, dims, it):
    """Iterate over the indices lists, yielding _stepIter(it) for each"""
    d = [0, 0, 1]
    d[1] = (dims[3] - dims[0])
    d[0] = (dims[4] - dims[1]) * d[1]
    mins = [dims[2], dims[1], dims[0]]
    offsets = ([(a - mins[ordmap[x]]) * d[ordmap[x]]
                for a in indices[x]]
                for x in range(3))
    for ioff, joff, koff in itertools.product(*offsets):
        yield _stepIter(it, (ioff + joff + koff))
