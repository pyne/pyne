import copy
import itertools
from collections import namedtuple, Iterable, Sequence

import numpy as np

try:
    from itaps import iMesh, iBase, iMeshExtensions
except ImportError:
    raise ImportError("The mesh module requires imports of iMesh, iBase, and"
          " iMeshExtensions from PyTAPS")

from pyne.material import Material, MaterialLibrary

# dictionary of lamba functions for mesh arithmetic
_ops = {"+": lambda val_1, val_2: (val_1 + val_2), 
        "-": lambda val_1, val_2: (val_1 - val_2),
        "*": lambda val_1, val_2: (val_1 * val_2),
        "/": lambda val_1, val_2: (val_1 / val_2)}

err__ops = {"+": lambda val_1, val_2, val_1_err, val_2_err: \
                 (1/(val_1 + val_2)*np.sqrt((val_1*val_1_err)**2 \
                  + (val_2*val_2_err)**2)), 
            "-": lambda val_1, val_2, val_1_err, val_2_err: \
                 (1/(val_1 - val_2)*np.sqrt((val_1*val_1_err)**2 \
                 + (val_2*val_2_err)**2)),
            "*": lambda val_1, val_2, val_1_err, val_2_err: \
                 (np.sqrt(val_1_err**2 + val_2_err**2)),
                  "/": lambda val_1, val_2, val_1_err, val_2_err: \
                 (np.sqrt(val_1_err**2 + val_2_err**2))}

_INTEGRAL_TYPES = (int, np.integer, np.bool_)

class Tag(object):
    """A mesh tag, which acts as a descriptor on the mesh.  This dispatches
    access to intrinsic material properties, the iMesh.Mesh tags, and material
    metadata attributes.
    """

    def __init__(self, mesh=None, name=None, doc=None):
        """Parameters
        ----------
        mesh : Mesh, optional
            The PyNE mesh to tag.
        name : str, optional
            The name of the tag.
        doc : str, optional
            Documentation string for the tag.

        """
        if mesh is None or name is None:
            self._lazy_args = {'mesh': mesh, 'name': name, 'doc': doc}
            return 
        self.mesh = mesh
        self.name = name
        mesh.tags[name] = self
        if doc is None:
            doc = "the {0!r} tag".format(name)
        self.__doc__ = doc
        if hasattr(self, '_lazy_args'):
            del self._lazy_args

    def __str__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    def __repr__(self):
        return "{0}(name={1!r}, doc={2!r})".format(self.__class__.__name__, self.name,
                                                   self.__doc__)

    def __get__(self, mesh, objtype=None):
        return self

    def __set__(self, mesh, value):
        if not isinstance(value, Tag):
            raise AttributeError("can't set tag from non-tag objects, "
                                 "got {0}".format(type(value)))
        if self.name != value.name:
            raise AttributeError("tags names must match, found "
                                 "{0} and {1}".format(self.name, value.name))
        self[:] = value[:]

    def __delete__(self, mesh):
        del self[:]

class MaterialPropertyTag(Tag):
    """A mesh tag which looks itself up as a material property (attribute).
    This makes the following expressions equivalent for a given material property
    name::

        mesh.name[i] == mesh.mats[i].name

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.
    """

    def __getitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return getattr(mats[key], name)
        elif isinstance(key, slice):
            return np.array([getattr(mats[i],name) for i in range(*key.indices(size))])
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return np.array([getattr(mats[i], name) for i, b in enumerate(key) if b])
        elif isinstance(key, Iterable):
            return np.array([getattr(mats[i], name) for i in key])
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __setitem__(self, key, value):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            setattr(mats[key], name, value)
        elif isinstance(key, slice):
            idx = range(*key.indices(size))
            if isinstance(value, Sequence) and len(value) == len(idx):
                for i, v in zip(idx, value):
                    setattr(mats[i], name, v)
            else:
                for i in idx:
                    setattr(mats[i], name, value) 
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            idx = np.where(key)[0]
            if isinstance(value, Sequence) and len(value) == key.sum():
                for i, v in zip(idx, value):
                    setattr(mats[i], name, v)
            else:
                for i in idx:
                    setattr(mats[i], name, value) 
        elif isinstance(key, Iterable):
            if isinstance(value, Sequence) and len(value) == len(key):
                for i, v in zip(key, value):
                    setattr(mats[i], name, v)
            else:
                for i in key:
                    setattr(mats[i], name, value) 
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __delitem__(self, key):
        msg = "the material property tag {0!r} may not be deleted".format(self.name)
        raise AttributeError(msg)


class MaterialMethodTag(Tag):
    """A mesh tag which looks itself up by calling a material method which takes 
    no arguments.  This makes the following expressions equivalent for a given 
    material method name::

        mesh.name[i] == mesh.mats[i].name()

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.
    """

    def __getitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return getattr(mats[key], name)()
        elif isinstance(key, slice):
            return np.array([getattr(mats[i], name)() for i in \
                                                          range(*key.indices(size))])
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return np.array([getattr(mats[i], name)() for i, b in enumerate(key) if b])
        elif isinstance(key, Iterable):
            return np.array([getattr(mats[i], name)() for i in key])
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __setitem__(self, key, value):
        msg = "the material method tag {0!r} may not be set".format(self.name)
        raise AttributeError(msg)

    def __delitem__(self, key):
        msg = "the material method tag {0!r} may not be deleted".format(self.name)
        raise AttributeError(msg)


class MetadataTag(Tag):
    """A mesh tag which looks itself up as a material metadata attribute.
    Tags of this are untyped and may have any size.  Use this for catch-all tags.
    This makes the following expressions equivalent for a given material property
    name::

        mesh.name[i] == mesh.mats[i].attrs['name']

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.  
    """

    def __getitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return mats[key].attrs[name]
        elif isinstance(key, slice):
            return [mats[i].attrs[name] for i in range(*key.indices(size))]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return [mats[i].attrs[name] for i, b in enumerate(key) if b]
        elif isinstance(key, Iterable):
            return [mats[i].attrs[name] for i in key]
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __setitem__(self, key, value):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            mats[key].attrs[name] = value
        elif isinstance(key, slice):
            idx = range(*key.indices(size))
            if isinstance(value, Sequence) and len(value) == len(idx):
                for i, v in zip(idx, value):
                    mats[i].attrs[name] = v
            else:
                for i in idx:
                    mats[i].attrs[name] = value
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            idx = np.where(key)[0]
            if isinstance(value, Sequence) and len(value) == key.sum():
                for i, v in zip(idx, value):
                    mats[i].attrs[name] = v
            else:
                for i in idx:
                    mats[i].attrs[name] = value
        elif isinstance(key, Iterable):
            if isinstance(value, Sequence) and len(value) == len(key):
                for i, v in zip(key, value):
                    mats[i].attrs[name] = v
            else:
                for i in key:
                    mats[i].attrs[name] = value
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __delitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            del mats[key].attrs[name]
        elif isinstance(key, slice):
            for i in range(*key.indices(size)):
                del mats[i].attrs[name]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            for i, b in enumerate(key): 
                if b:
                    del mats[i].attrs[name] 
        elif isinstance(key, Iterable):
            for i in key:
                del mats[i].attrs[name]
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

class IMeshTag(Tag):
    """A mesh tag which looks itself up as a tag on the iMesh.Mesh instance.
    This makes the following expressions equivalent for a given iMesh.Mesh tag
    name::

        mesh.name[i] == mesh.mesh.getTagHandle(name)[list(mesh.mesh.iterate(
                                iBase.Type.region, iMesh.Topology.all))[i]]

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.
    """

    def __init__(self, size=1, dtype='f8', mesh=None, name=None, doc=None):
        """Parameters
        ----------
        size : int, optional
            The number of elements of type dtype that this tag stores.
        dtype : np.dtype or similar, optional
            The data type of this tag from int, float, and byte. See PyTAPS
            tags for more details.
        mesh : Mesh, optional
            The PyNE mesh to tag.
        name : str, optional
            The name of the tag.
        doc : str, optional
            Documentation string for the tag.

        """
        super(IMeshTag, self).__init__(mesh=mesh, name=name, doc=doc)
        if mesh is None or name is None:
            self._lazy_args['size'] = size
            self._lazy_args['dtype'] = dtype
            return 
        try:
            self.tag = self.mesh.mesh.getTagHandle(self.name)
        except iBase.TagNotFoundError: 
            self.tag = self.mesh.mesh.createTag(self.name, size, dtype)

    def __delete__(self, mesh):
        super(IMeshTag, self).__delete__(mesh)
        self.mesh.mesh.destroyTag(self.name, force=True)

    def __getitem__(self, key):
        m = self.mesh.mesh
        size = len(self.mesh)
        mtag = self.tag
        miter = m.iterate(iBase.Type.region, iMesh.Topology.all)
        if isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError("key index {0} greater than the size of the "
                                 "mesh {1}".format(key, size))
            for i_ve in zip(range(key+1), miter):
                pass
            return mtag[i_ve[1]]
        elif isinstance(key, slice):
            return mtag[list(miter)[key]]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return mtag[[ve for b, ve in zip(key, miter) if b]]
        elif isinstance(key, Iterable):
            ves = list(miter)
            return mtag[[ves[i] for i in key]]
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __setitem__(self, key, value):
        m = self.mesh.mesh
        size = len(self.mesh)
        mtag = self.tag
        miter = m.iterate(iBase.Type.region, iMesh.Topology.all)
        if isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError("key index {0} greater than the size of the "
                                 "mesh {1}".format(key, size))
            for i_ve in zip(range(key+1), miter):
                pass
            mtag[i_ve[1]] = value
        elif isinstance(key, slice):
            idx = range(*key.indices(size))
            if not (isinstance(value, Sequence) and len(value) == len(idx)):
                value = [value] * len(idx)
            mtag[list(miter)[key]] = value
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            ntrues = key.sum()
            if not (isinstance(value, Sequence) and len(value) == ntrues):
                value = [value] * ntrues                
            mtag[[ve for b, ve in zip(key, miter) if b]] = value
        elif isinstance(key, Iterable):
            if not (isinstance(value, Sequence) and len(value) == len(key)):
                value = [value] * len(key)
            ves = list(miter)
            mtag[[ves[i] for i in key]] = value
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __delitem__(self, key):
        m = self.mesh.mesh
        size = len(self.mesh)
        mtag = self.tag
        miter = m.iterate(iBase.Type.region, iMesh.Topology.all)
        if isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError("key index {0} greater than the size of the "
                                 "mesh {1}".format(key, size))
            for i_ve in zip(range(key+1), miter):
                pass
            del mtag[i_ve[1]]
        elif isinstance(key, slice):
            del mtag[list(miter)[key]]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            del mtag[[ve for b, ve in zip(key, miter) if b]]
        elif isinstance(key, Iterable):
            ves = list(miter)
            del mtag[[ves[i] for i in key]]
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

class ComputedTag(Tag):
    '''A mesh tag which looks itself up by calling a function (or other callable)
    with the following signature::

        def f(mesh, i):
            """mesh is a pyne.mesh.Mesh() object and i is the volume element index
            to compute.
            """
            # ... do some work ...
            return anything_you_want

    This makes the following expressions equivalent for a given computed tag
    name::

        mesh.name[i] == f(mesh, i)

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.

    Notes
    -----
    The results of computed tags are not stored and the function object itself 
    is also not persisted.  Therefore, you must manually re-tag the mesh with
    the desired functions each session.

    '''

    def __init__(self, f, mesh=None, name=None, doc=None):
        """Parameters
        ----------
        f : callable object
            The function that performs the computation.
        mesh : Mesh, optional
            The PyNE mesh to tag.
        name : str, optional
            The name of the tag.
        doc : str, optional
            Documentation string for the tag.

        """
        doc = doc or f.__doc__
        super(ComputedTag, self).__init__(mesh=mesh, name=name, doc=doc)
        if mesh is None or name is None:
            self._lazy_args['f'] = f
            return 
        self.f = f

    def __getitem__(self, key):
        m = self.mesh
        f = self.f
        size = len(m)
        if isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError("key index {0} greater than the size of the "
                                 "mesh {1}".format(key, size))
            return f(m, key)
        elif isinstance(key, slice):
            return [f(m, i) for i in range(*key.indices(size))]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return [f(m, i) for i, b in enumerate(key) if b]
        elif isinstance(key, Iterable):
            return [f(m, i) for i in key]
        else:
            raise TypeError("{0} is not an int, slice, mask, "
                            "or fancy index.".format(key))        

    def __setitem__(self, key, value):
        msg = "the computed tag {0!r} may not be set".format(self.name)
        raise AttributeError(msg)

    def __delitem__(self, key):
        msg = "the computed tag {0!r} may not be deleted".format(self.name)
        raise AttributeError(msg)



class MeshError(Exception):
    """Errors related to instantiating mesh objects and utilizing their methods.
    """
    pass

class Mesh(object):        
    """This class houses any iMesh instance and contains methods for various mesh
    operations. Special methods exploit the properties of structured mesh.

    Attributes
    ----------
    mesh : iMesh instance
    mesh_file : string
        File name of file containing iMesh instance.
    structured : bool
        True for structured mesh.
    structured_coords : list of lists
        A list containing lists of x_points, y_points and z_points that make up
        a structured mesh. 
    structured_set : iMesh entity set handle
        A preexisting structured entity set on an iMesh instance with a
        "BOX_DIMS" tag.

        Unstructured mesh instantiation:
             - From iMesh instance by specifying: <mesh>
             - From mesh file by specifying: <mesh_file>

        Structured mesh instantiation:
            - From iMesh instance with exactly 1 entity set (with BOX_DIMS tag)
              by specifying <mesh> and structured = True.
            - From mesh file with exactly 1 entity set (with BOX_DIMS tag) by
              specifying <mesh_file> and structured = True.
            - From an imesh instance with multiple entity sets by specifying 
              <mesh>, <structured_set>, structured=True.
            - From coordinates by specifying <structured_coords>,
              structured=True, and optional preexisting iMesh instance <mesh>

        The "BOX_DIMS" tag on iMesh instances containing structured mesh is
        a vector of floats it the following form:
        [i_min, j_min, k_min, i_max, j_max, k_max]
        where each value is a volume element index number. Typically volume 
        elements should be indexed from 0. The "BOX_DIMS" information is stored
        in self.dims.
    mats : MaterialLibrary or dict or Materials or None
        This is a mapping of volume element handles to Material objects.

    """

    

    def __init__(self, mesh=None, mesh_file=None, structured=False, \
                 structured_coords=None, structured_set=None, mats=None):
        if mesh:
            self.mesh = mesh
        else: 
            self.mesh = iMesh.Mesh()

        self.structured = structured

        #Unstructured mesh cases
        if not self.structured:
            #Error if structured arguments are passed
            if structured_coords or structured_set:
                MeshError("Structured mesh arguments should not be present for\
                            unstructured Mesh instantiation.")

            #From imesh instance
            if mesh and not mesh_file:
                pass
            #From file
            elif mesh_file and not mesh:
                self.mesh.load(mesh_file)
                self.mats = MaterialLibrary(mesh_file)
            else:
                raise MeshError("To instantiate unstructured mesh object, "
                                 "must supply exactly 1 of the following: "
                                 "<mesh>, <mesh_file>.")

        #structured mesh cases
        elif self.structured:
            #From mesh or mesh_file
            if (mesh or mesh_file) and not structured_coords \
                                   and not structured_set:
                if mesh_file:
                    self.mesh.load(mesh_file)
                    self.mats = MaterialLibrary(mesh_file)
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
                        self.structured_set = ent_set
                        count += 1

                if count == 0:
                    raise MeshError("Found no structured meshes in "
                                    "file {0}".format(mesh_file))
                elif count > 1:
                    raise MeshError("Found {0} structured meshes."
                                    " Instantiate individually using"
                                    " from_ent_set()".format(count))
            # from coordinates                       
            elif not mesh and not mesh_file and structured_coords \
                                            and not structured_set:
                extents = [0, 0, 0] + [len(x) - 1 for x in structured_coords]
                self.structured_set = self.mesh.createStructuredMesh(
                     extents, i=structured_coords[0], j=structured_coords[1], 
                     k=structured_coords[2], create_set=True)

            #From mesh and structured_set:
            elif mesh and not mesh_file and not structured_coords \
                                        and structured_set:
                try:
                    self.mesh.getTagHandle("BOX_DIMS")[structured_set]
                except iBase.TagNotFoundError as e:
                    print("Supplied entity set does not contain BOX_DIMS tag")
                    raise e

                self.structured_set = structured_set
            else:
                raise MeshError("For structured mesh instantiation, need to"
                                "supply exactly one of the following:\n"
                                "A. iMesh instance\n"
                                "B. Mesh file\n"
                                "C. Mesh coordinates\n"
                                "D. Structured entity set AND iMesh instance")

            self.dims = self.mesh.getTagHandle("BOX_DIMS")[self.structured_set]
            self.vertex_dims = list(self.dims[0:3]) \
                               + [x + 1 for x in self.dims[3:6]]
        # sets mats
        if mats is None:
            mats = MaterialLibrary()
        elif not isinstance(mats, MaterialLibrary):
            mats = MaterialLibrary(mats)
        self.mats = mats

        # tag with volume id and ensure mats exist.
        ves = list(self.mesh.iterate(iBase.Type.region, iMesh.Topology.all))
        tags = self.mesh.getAllTags(ves[0])
        tags = set(tag.name for tag in tags)
        if 've_idx' in tags:
            tag_ve_idx = self.mesh.getTagHandle('ve_idx')
        else:
            tag_ve_idx = self.mesh.createTag('ve_idx', 1, int)
        for i, ve in enumerate(ves):
            tag_ve_idx[ve] = i
            if i not in mats:
                mats[i] = Material()
        self._len = i + 1

        # Default tags
        self.tags = {}
        # metadata tags, these should come first so they don't accidentally 
        # overwite hard coded tag names.
        metatagnames = set()
        for mat in mats.values():
            metatagnames.update(mat.attrs.keys())
        for name in metatagnames:
            setattr(self, name, MetadataTag(mesh=self, name=name))
        # iMesh.Mesh() tags
        tagnames = set()
        for ve in ves:
            tagnames.update(t.name for t in self.mesh.getAllTags(ve))
        for name in tagnames:
            setattr(self, name, IMeshTag(mesh=self, name=name)) 
        # Material property tags
        self.atoms_per_mol = MaterialPropertyTag(mesh=self, name='atoms_per_mol', 
                                                 doc='Number of atoms per molecule')
        self.attrs = MaterialPropertyTag(mesh=self, name='attrs', 
                        doc='metadata attributes, stored on the material')
        self.comp = MaterialPropertyTag(mesh=self, name='comp', 
                doc="normalized composition mapping from nuclides to mass fractions")
        self.mass = MaterialPropertyTag(mesh=self, name='mass', 
                                        doc='the mass of the material')
        self.density = MaterialPropertyTag(mesh=self, name='density', 
                                           doc='the density [g/cc]')
        # Material method tags
        methtagnames = ('expand_elements', 'mass_density', 'molecular_weight', 
                        'mult_by_mass', 'number_density', 'sub_act', 'sub_fp', 
                        'sub_lan', 'sub_ma', 'sub_tru', 'to_atom_frac')
        for name in methtagnames:
            doc = "see Material.{0}() for more information".format(name)
            setattr(self, name, MaterialMethodTag(mesh=self, name=name, doc=doc))
        
        
    def __len__(self):
        return self._len

    def __iter__(self):
        """Iterates through the mesh and at each step yield the volume element 
        index i, the material mat, and the volume element itself ve.
        """
        mats = self.mats
        for i, ve in enumerate(self.mesh.iterate(iBase.Type.region, 
                                                 iMesh.Topology.all)):
            yield i, mats[i], ve

    def __contains__(self, i):
        return i < len(self)

    def __setattr__(self, name, value):
        if isinstance(value, Tag) and hasattr(value, '_lazy_args'):
            # some 1337 1Azy 3\/a1
            kwargs = value._lazy_args
            kwargs['mesh'] = self if kwargs['mesh'] is None else kwargs['mesh']
            kwargs['name'] = name if kwargs['name'] is None else kwargs['name']
            value = type(value)(**kwargs)
        super(Mesh, self).__setattr__(name, value)

    def tag(self, name, value=None, tagtype=None, doc=None, size=None, dtype=None):
        """Adds a new tag to the mesh, guessing the approriate place to store the
        data.

        Parameters
        ----------
        name : str
            The tag name
        value : optional
            The value to initialize the tag with, skipped if None.
        tagtype : Tag or str, optional
            The type of tag this should be any of the following classes or 
            strings are accepted: IMeshTag, MetadataTag, ComputedTag, 'imesh', 
            'metadata', or 'computed'.
        doc : str, optional
            The tag documentation string.
        size : int, optional
            The size of the tag. This only applies to IMeshTags.
        dtype : numpy dtype, optional
            The data type of the tag. This only applies to IMeshTags. See PyTAPS
            for more details.

        """
        if name in self.tags:
            raise KeyError('{0} tag already exists on the mesh'.format(name))
        if tagtype is None:
            if callable(value):
                tagtype = ComputedTag
            elif size is None and dtype is not None:
                size = 1
                tagtype = IMeshTag
            elif size is not None and dtype is None:
                dtype = 'f8'
                tagtype = IMeshTag
            elif value is None:
                size = 1
                value = 0.0
                dtype = 'f8'
                tagtype = IMeshTag
            elif isinstance(value, float):
                size = 1
                dtype = 'f8'
                tagtype = IMeshTag
            elif isinstance(value, int):
                size = 1
                dtype = 'i'
                tagtype = IMeshTag
            elif isinstance(value, str):
                tagtype = MetadataTag                
            elif isinstance(value, Sequence):
                raise ValueError('ambiguous tag {0!r} creation when value is a '
                        'sequence, please set tagtype, size, or dtype'.format(name))
            else:
                tagtype = MetadataTag
        if tagtype is IMeshTag or tagtype.lower() == 'imesh':
            t = IMeshTag(size=size, dtype=dtype, mesh=self, name=name, doc=doc)
        elif tagtype is MetadataTag or tagtype.lower() == 'metadata':
            t = MetadataTag(mesh=self, name=name, doc=doc)
        elif tagtype is ComputedTag or tagtype.lower() == 'computed':
            t = ComputedTag(f=value, mesh=self, name=name, doc=doc)
        else:
            raise ValueError('tagtype {0} not valid'.format(tagtype))
        if value is not None and tagtype is not ComputedTag:
            t[:] = value
        setattr(self, name, t)


#    def __add__(self, other):
#        """Adds the common tags of other and returns a new mesh object.
#        """
#        tags = self.common_ve_tags(other)
#        return self._do_op(other, tags, "+", in_place=False)
#
#    def __sub__(self, other):
#        """Subtracts the common tags of other and returns a new mesh object.
#        """
#        tags = self.common_ve_tags(other)
#        return self._do_op(other, tags, "-", in_place=False)
#
#    def __mul__(self, other):
#        """Multiplies the common tags of other and returns a new mesh object.
#        """
#        tags = self.common_ve_tags(other)
#        return  self._do_op(other, tags, "*", in_place=False)
#
#    def __div__(self, other):
#        """Adds the common tags of other and returns a new mesh object.
#        """
#        tags = self.common_ve_tags(other)
#        return self._do_op(other, tags, "/", in_place=False)


    def add(self, other):
        """Adds the common tags of other to the mesh object.
        """
        tags = self.common_ve_tags(other)
        self._do_op(other, tags, "+")

    def sub(self, other):
        """Substracts the common tags of other to the mesh object.
        """
        tags = self.common_ve_tags(other)
        self._do_op(other, tags, "-")

    def mul(self, other):
        """Multiplies the common tags of other to the mesh object.
        """
        tags = self.common_ve_tags(other)
        self._do_op(other, tags, "*")

    def div(self, other):
        """Divides the common tags of other to the mesh object.
        """
        tags = self.common_ve_tags(other)
        self._do_op(other, tags, "/")

    def _do_op(self, other, tags, op, in_place=True):
        """Private function to do mesh +, -, *, /.
        """
        # Exclude error tags in a case a StatMesh is mistakenly initialized as a
        # Mesh object.
        tags = set(tag for tag in tags if not tag.endswith('_error'))

        if in_place:
            mesh_1 = self
        else:
            mesh_1 = copy.copy(self)

        for tag in tags:
            for ve_1, ve_2 in \
                zip(zip(iter(mesh_1.mesh.iterate(iBase.Type.region, iMesh.Topology.all))),
                    zip(iter(other.mesh.iterate(iBase.Type.region, iMesh.Topology.all)))):
                self.mesh.getTagHandle(tag)[ve_1] = \
                    _ops[op](mesh_1.mesh.getTagHandle(tag)[ve_1], 
                            other.mesh.getTagHandle(tag)[ve_2])

        return mesh_1


    def common_ve_tags(self, other):
        """Returns the volume element tags in common between self and other.
        """
        self_tags = self.mesh.getAllTags(list(self.mesh.iterate(
                                     iBase.Type.region, iMesh.Topology.all))[0])
        other_tags = other.mesh.getAllTags(list(other.mesh.iterate(iBase.Type.region, 
                                           iMesh.Topology.all))[0])
        self_tags = set(x.name for x in self_tags)
        other_tags = set(x.name for x in other_tags)
        intersect = self_tags & other_tags
        intersect.discard('ve_idx')
        return intersect
                           
    def __copy__(self):
        #first copy full imesh instance
        imesh_copy = iMesh.Mesh()

        #now create Mesh objected from copied iMesh instance
        mesh_copy = Mesh(mesh=imesh_copy, structured=copy.copy(self.structured))
        return mesh_copy

    #Structured methods:
    def structured_get_vertex(self, i, j, k):
        """Return the handle for (i,j,k)'th vertex in the mesh"""
        self._structured_check()
        n = _structured_find_idx(self.vertex_dims, (i, j, k))
        return _structured_step_iter(
            self.structured_set.iterate(iBase.Type.vertex, 
                                        iMesh.Topology.point), n)


    def structured_get_hex(self, i, j, k):
        """Return the handle for the (i,j,k)'th hexahedron in the mesh"""
        self._structured_check()
        n = _structured_find_idx(self.dims, (i, j, k))
        return _structured_step_iter(
            self.structured_set.iterate(iBase.Type.region, 
                                 iMesh.Topology.hexahedron), n)


    def structured_get_hex_volume(self, i, j, k):
        self._structured_check()
        """Return the volume of the (i,j,k)'th hexahedron in the mesh"""
        v = list(self.structured_iterate_vertex(x=[i, i + 1],
                                 y=[j, j + 1],
                                 z=[k, k + 1]))
        coord = self.mesh.getVtxCoords(v)
        dx = coord[1][0] - coord[0][0]
        dy = coord[2][1] - coord[0][1]
        dz = coord[4][2] - coord[0][2]
        return dx * dy * dz


    def structured_iterate_hex(self, order="zyx", **kw):
        """Get an iterator over the hexahedra of the mesh

        The order argument specifies the iteration order.  It must be a string
        of 1-3 letters from the set (x,y,z).  The rightmost letter is the axis
        along which the iteration will advance the most quickly.  Thus "zyx" --
        x coordinates changing fastest, z coordinates changing least fast-- is
        the default, and is identical to the order that would be given by the
        structured_set.iterate() function.

        When a dimension is absent from the order, iteration will proceed over
        only the column in the mesh that has the lowest corresonding (i/j/k)
        coordinate.  Thus, with order "xy," iteration proceeds over the i/j
        plane of the structured mesh with the smallest k coordinate.

        Specific slices can be specified with keyword arguments:

        Keyword args::

          x: specify one or more i-coordinates to iterate over.
          y: specify one or more j-coordinates to iterate over.
          z: specify one or more k-coordinates to iterate over.

        Examples::

          structured_iterate_hex(): equivalent to iMesh iterator over hexes in mesh
          structured_iterate_hex("xyz"): iterate over entire mesh, with k-coordinates
                                         changing fastest, i-coordinates least fast.
          structured_iterate_hex("yz", x=3): Iterate over the j-k plane of the mesh
                                             whose i-coordinate is 3, with k values
                                             changing fastest.
          structured_iterate_hex("z"): Iterate over k-coordinates, with i=dims.imin
                             and j=dims.jmin
          structured_iterate_hex("yxz", y=(3,4)): Iterate over all hexes with
                                        j-coordinate = 3 or 4.  k-coordinate
                                        values change fastest, j-values least
                                        fast.
        """
        self._structured_check()

        # special case: zyx order is the standard pytaps iteration order,
        # so we can save time by simply returning a pytaps iterator
        # if no kwargs were specified
        if order == "zyx" and not kw:
            return self.structured_set.iterate(iBase.Type.region,
                                       iMesh.Topology.hexahedron)

        indices, ordmap = _structured_iter_setup(self.dims, order, **kw)
        return _structured_iter(indices, ordmap, self.dims, 
            self.structured_set.iterate(iBase.Type.region, 
                                        iMesh.Topology.hexahedron))


    def structured_iterate_vertex(self, order="zyx", **kw):
        """Get an iterator over the vertices of the mesh

        See structured_iterate_hex() for an explanation of the order argument and the
        available keyword arguments.
        """
        self._structured_check()
        #special case: zyx order without kw is equivalent to pytaps iterator
        if order == "zyx" and not kw:
            return self.structured_set.iterate(iBase.Type.vertex,
                                               iMesh.Topology.point)

        indices, ordmap = _structured_iter_setup(self.vertex_dims, order, **kw)
        return _structured_iter(indices, ordmap, self.vertex_dims, 
                self.structured_set.iterate(iBase.Type.vertex, 
                                            iMesh.Topology.point))


    def structured_iterate_hex_volumes(self, order="zyx", **kw):
        """Get an iterator over the volumes of the mesh hexahedra

        See structured_iterate_hex() for an explanation of the order argument and the
        available keyword arguments.
        """
        self._structured_check()
        indices, _ = _structured_iter_setup(self.dims, order, **kw)
        # Use an inefficient but simple approach: call structured_get_hex_volume()
        # on each required i,j,k pair.  
        # A better implementation would only make one call to getVtxCoords.
        for A in itertools.product(*indices):
            # the ordmap returned from _structured_iter_setup maps to kji/zyx ordering,
            # but we want ijk/xyz ordering, so create the ordmap differently.
            ordmap = [order.find(L) for L in "xyz"]
            ijk = [A[ordmap[x]] for x in range(3)]
            yield self.structured_get_hex_volume(*ijk)


    def structured_get_divisions(self, dim):
        """Get the mesh divisions on a given dimension

        Given a dimension "x", "y", or "z", return a list of the mesh vertices
        along that dimension.
        """
        self._structured_check()
        if len(dim) == 1 and dim in "xyz":
            idx = "xyz".find(dim)
            return [self.mesh.getVtxCoords(i)[idx]
                    for i in self.structured_iterate_vertex(dim)]
        else:
            raise MeshError("Invalid dimension: {0}".format(str(dim)))

    def _structured_check(self):
        if not self.structured:
            raise MeshError("Structured mesh methods cannot be called from "\
                            "unstructured mesh instances.")

    def write_hdf5(self, filename):
        """Writes the mesh to an hdf5 file."""
        self.mesh.save(filename)
        self.mats.write_hdf5(filename)

######################################################
# private helper functions for structured mesh methods
######################################################

def _structured_find_idx(dims, ijk):
    """Helper method fo structured_get_vertex and structured_get_hex.

    For tuple (i,j,k), return the number N in the appropriate iterator.
    """
    dim0 = [0] * 3
    for i in xrange(0, 3):
        if (dims[i] > ijk[i] or dims[i + 3] <= ijk[i]):
            raise MeshError(str(ijk) + " is out of bounds")
        dim0[i] = ijk[i] - dims[i]
    i0, j0, k0 = dim0
    n = (((dims[4] - dims[1]) * (dims[3] - dims[0]) * k0) +
         ((dims[3] - dims[0]) * j0) +
         i0)
    return n


def _structured_step_iter(it, n):
    """Helper method for structured_get_vertex and structured_get_hex

    Return the nth item in the iterator."""
    it.step(n)
    r = it.next()
    it.reset()
    return r


def _structured_iter_setup(dims, order, **kw):
    """Setup helper function for StrMesh iterator functions

    Given dims and the arguments to the iterator function, return
    a list of three lists, each being a set of desired coordinates,
    with fastest-changing coordinate in the last column),
    and the ordmap used by _structured_iter to reorder each coodinate to (i,j,k).
    """
    # a valid order has the letters "x", "y", and "z"
    # in any order without duplicates
    if not (len(order) <= 3 and
            len(set(order)) == len(order) and
            all([a in "xyz" for a in order])):
        raise MeshError("Invalid iteration order: " + str(order))

    # process kw for validity
    spec = {}
    for idx, d in enumerate("xyz"):
        if d in kw:
            spec[d] = kw[d]
            if not isinstance(spec[d], Iterable):
                spec[d] = [spec[d]]
            if not all(x in range(dims[idx], dims[idx + 3])
                    for x in spec[d]):
                raise MeshError( \
                        "Invalid iterator kwarg: {0}={1}".format(d, spec[d]))
            if d not in order and len(spec[d]) > 1:
                raise MeshError("Cannot iterate over" + str(spec[d]) +
                                   "without a proper iteration order")
        if d not in order:
            order = d + order
            spec[d] = spec.get(d, [dims[idx]])

    # get indices and ordmap
    indices = []
    for L in order:
        idx = "xyz".find(L)
        indices.append(spec.get(L, xrange(dims[idx], dims[idx + 3])))

    ordmap = ["zyx".find(L) for L in order]
    return indices, ordmap


def _structured_iter(indices, ordmap, dims, it):
    """Iterate over the indices lists, yielding _structured_step_iter(it) for each.
    """
    d = [0, 0, 1]
    d[1] = (dims[3] - dims[0])
    d[0] = (dims[4] - dims[1]) * d[1]
    mins = [dims[2], dims[1], dims[0]]
    offsets = ([(a - mins[ordmap[x]]) * d[ordmap[x]]
                for a in indices[x]]
                for x in range(3))
    for ioff, joff, koff in itertools.product(*offsets):
        yield _structured_step_iter(it, (ioff + joff + koff))


class StatMesh(Mesh):
    def __init__(self, mesh=None, mesh_file=None, structured=False,
                 structured_coords=None, structured_set=None):

        super(StatMesh, self).__init__(mesh=mesh, mesh_file=mesh_file, 
              structured=structured, structured_coords=structured_coords, 
              structured_set=structured_set)

    def _do_op(self, other, tags, op, in_place=True):
        """Private function to do mesh +, -, *, /. Called by operater overloading
        functions.
        """
        # Exclude error tags because result and error tags are treated simotaneously
        # so there is not need to include both in the tag list to iterate through.
        tags = set(tag for tag in tags if not tag.endswith('_error'))

        if in_place:
            mesh_1 = self
        else:
            mesh_1 = copy.copy(self)

        for tag in tags:
            for ve_1, ve_2 in \
                zip(zip(iter(mesh_1.mesh.iterate(iBase.Type.region, iMesh.Topology.all))),
                    zip(iter(other.mesh.iterate(iBase.Type.region, iMesh.Topology.all)))):

                mesh_1.mesh.getTagHandle(tag + "_error")[ve_1] = err__ops[op](
                    mesh_1.mesh.getTagHandle(tag)[ve_1], 
                    other.mesh.getTagHandle(tag)[ve_2], 
                    mesh_1.mesh.getTagHandle(tag + "_error")[ve_1], 
                    other.mesh.getTagHandle(tag + "_error")[ve_2])

                mesh_1.mesh.getTagHandle(tag)[ve_1] = \
                    _ops[op](mesh_1.mesh.getTagHandle(tag)[ve_1], 
                            other.mesh.getTagHandle(tag)[ve_2])

        return mesh_1

