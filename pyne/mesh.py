from __future__ import print_function, division
from future.utils import implements_iterator
from pyne.material import Material, MultiMaterial
from pyne.material_library import MaterialLibrary
import sys
import copy
import itertools

try:
    from collections.abc import Iterable, Sequence
except ImportError:
    from collections import Iterable, Sequence
from warnings import warn
from pyne.utils import QA_warn, check_iterable

import numpy as np
import tables as tb

QA_warn(__name__)

try:
    from pymoab import core as mb_core, hcoord, scd, types
    from pymoab.rng import subtract
    from pymoab.tag import Tag
    from pymoab.types import _eh_py_type, _TAG_TYPE_STRS

    HAVE_PYMOAB = True

except ImportError:
    HAVE_PYMOAB = False
    warn(
        "The PyMOAB optional dependency could not be imported. "
        "Some aspects of the mesh module may be incomplete.",
        ImportWarning,
    )


_BOX_DIMS_TAG_NAME = "BOX_DIMS"

if sys.version_info[0] > 2:
    basestring = str

# dictionary of lamba functions for mesh arithmetic
_ops = {
    "+": lambda val_1, val_2: (val_1 + val_2),
    "-": lambda val_1, val_2: (val_1 - val_2),
    "*": lambda val_1, val_2: (val_1 * val_2),
    "/": lambda val_1, val_2: (val_1 / val_2),
}

err__ops = {
    "+": lambda val_1, val_2, val_1_err, val_2_err: (
        1
        / (val_1 + val_2)
        * np.sqrt((val_1 * val_1_err) ** 2 + (val_2 * val_2_err) ** 2)
    ),
    "-": lambda val_1, val_2, val_1_err, val_2_err: (
        1
        / (val_1 - val_2)
        * np.sqrt((val_1 * val_1_err) ** 2 + (val_2 * val_2_err) ** 2)
    ),
    "*": lambda val_1, val_2, val_1_err, val_2_err: (
        np.sqrt(val_1_err**2 + val_2_err**2)
    ),
    "/": lambda val_1, val_2, val_1_err, val_2_err: (
        np.sqrt(val_1_err**2 + val_2_err**2)
    ),
}

_INTEGRAL_TYPES = (int, np.integer, np.bool_)
_SEQUENCE_TYPES = (Sequence, np.ndarray)


class Tag(object):
    """A mesh tag, which acts as a descriptor on the mesh.  This dispatches
    access to intrinsic material properties, PyMOAB tags, and material
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
            self._lazy_args = {"mesh": mesh, "name": name, "doc": doc}
            return
        self.mesh = mesh
        self.name = name
        mesh.tags[name] = self
        if doc is None:
            doc = "the {0!r} tag".format(name)
        self.__doc__ = doc
        if hasattr(self, "_lazy_args"):
            del self._lazy_args

    def __str__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    def __repr__(self):
        return "{0}(name={1!r}, doc={2!r})".format(
            self.__class__.__name__, self.name, self.__doc__
        )

    def __get__(self, mesh, objtype=None):
        return self

    def __set__(self, mesh, value):
        if not isinstance(value, Tag):
            raise AttributeError(
                "can't set tag from non-tag objects, " "got {0}".format(type(value))
            )
        if self.name != value.name:
            raise AttributeError(
                "tags names must match, found "
                "{0} and {1}".format(self.name, value.name)
            )
        self[:] = value[:]

    def __delete__(self, mesh):
        del mesh.tags[self.name]
        del self[:]

    def delete(self, mesh=None):
        if mesh == None:
            mesh = self.mesh
        del self[:]
        del mesh.tags[self.name]


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
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return getattr(mats[int(key)], name)
        elif isinstance(key, slice):
            return np.array([getattr(mats[i], name) for i in range(*key.indices(size))])
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length of the mesh.")
            return np.array([getattr(mats[i], name) for i, b in enumerate(key) if b])
        elif isinstance(key, Iterable):
            return np.array([getattr(mats[i], name) for i in key])
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __setitem__(self, key, value):
        name = self.name
        mats = self.mesh.mats
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            setattr(mats[key], name, value)
        elif isinstance(key, slice):
            idx = range(*key.indices(size))
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == len(idx):
                for i, v in zip(idx, value):
                    setattr(mats[i], name, v)
            else:
                for i in idx:
                    setattr(mats[i], name, value)
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match " "the length of the mesh.")
            idx = np.where(key)[0]
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == key.sum():
                for i, v in zip(idx, value):
                    setattr(mats[i], name, v)
            else:
                for i in idx:
                    setattr(mats[i], name, value)
        elif isinstance(key, Iterable):
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == len(key):
                for i, v in zip(key, value):
                    setattr(mats[i], name, v)
            else:
                for i in key:
                    setattr(mats[i], name, value)
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __delitem__(self, key):
        msg = ("the material property tag {0!r} may " "not be deleted").format(
            self.name
        )
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
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return getattr(mats[key], name)()
        elif isinstance(key, slice):
            return np.array(
                [getattr(mats[i], name)() for i in range(*key.indices(size))]
            )
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the " "length of the mesh.")
            return np.array([getattr(mats[i], name)() for i, b in enumerate(key) if b])
        elif isinstance(key, Iterable):
            return np.array([getattr(mats[i], name)() for i in key])
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __setitem__(self, key, value):
        msg = "the material method tag {0!r} may not be set".format(self.name)
        raise AttributeError(msg)

    def __delitem__(self, key):
        msg = ("the material method tag {0!r} may not be " "deleted").format(self.name)
        raise AttributeError(msg)


class MetadataTag(Tag):
    """A mesh tag which looks itself up as a material metadata attribute.
    Tags of this are untyped and may have any size.  Use this for catch-all
    tags. This makes the following expressions equivalent for a given material
    property.

    name::

        mesh.name[i] == mesh.mats[i].metadata['name']

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.
    """

    def __getitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            return mats[key].metadata[name]
        elif isinstance(key, slice):
            return [mats[i].metadata[name] for i in range(*key.indices(size))]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            return [mats[i].metadata[name] for i, b in enumerate(key) if b]
        elif isinstance(key, Iterable):
            return [mats[i].metadata[name] for i in key]
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __setitem__(self, key, value):
        name = self.name
        mats = self.mesh.mats
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            mats[key].metadata[name] = value
        elif isinstance(key, slice):
            idx = range(*key.indices(size))
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == len(idx):
                for i, v in zip(idx, value):
                    mats[i].metadata[name] = v
            else:
                for i in idx:
                    mats[i].metadata[name] = value
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            idx = np.where(key)[0]
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == key.sum():
                for i, v in zip(idx, value):
                    mats[i].metadata[name] = v
            else:
                for i in idx:
                    mats[i].metadata[name] = value
        elif isinstance(key, Iterable):
            if isinstance(value, _SEQUENCE_TYPES) and len(value) == len(key):
                for i, v in zip(key, value):
                    mats[i].metadata[name] = v
            else:
                for i in key:
                    mats[i].metadata[name] = value
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __delitem__(self, key):
        name = self.name
        mats = self.mesh.mats
        if mats is None:
            RuntimeError("Mesh.mats is None, please add a MaterialLibrary.")
        size = len(self.mesh)
        if isinstance(key, _INTEGRAL_TYPES):
            del mats[key].metadata[name]
        elif isinstance(key, slice):
            for i in range(*key.indices(size)):
                del mats[i].metadata[name]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            for i, b in enumerate(key):
                if b:
                    del mats[i].metadata[name]
        elif isinstance(key, Iterable):
            for i in key:
                del mats[i].metadata[name]
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )


class NativeMeshTag(Tag):
    """A mesh tag which looks itself up as a tag on a PyMOAB core instance.
    This makes the following expressions equivalent for a given PyNE/PyNE Mesh
    tag name::

        mesh.tag_name[i] == mesh.mesh.tag_get_data(mesh.mesh.tag_get_handle(name),
                                               mesh.mesh.get_entities_by_type(
                                               mesh.mesh.get_root_set(),
                                               types.MBHEX))[i]

    It also adds slicing, fancy indexing, boolean masking, and broadcasting
    features to this process.
    """

    def __init__(
        self,
        size=1,
        dtype="f8",
        default=0.0,
        mesh=None,
        name=None,
        doc=None,
        storage_type=None,
    ):
        """Parameters
        ----------
        size : int, optional
            The number of elements of type dtype that this tag stores.
        dtype : np.dtype or similar, optional
            The data type of this tag from int, float, and byte. See PyMOAB
            tags for more details.
        default : dtype or None, optional
            The default value to fill this tag with upon creation. If None,
            then the tag is created empty.
        mesh : Mesh, optional
            The PyNE mesh to tag.
        name : str, optional
            The name of the tag.
        doc : str, optional
            Documentation string for the tag.
        storage_type: str, optional
            MOAB tag storage type (MB_TAG_DENSE, MB_TAG_SPARSE, etc.)
            in advanced use of the database, this flag controls how this tag's
            data is stored in memory. The supported storage types are:
            MB_TYPE_SPARSE -  sparse tags are stored as a list of (entity
                handle, tag value) tuples, one list per sparse tag, sorted by
                entity handle.
            dense - MB_TAG_DENSE tag, values are stored in arrays which match
                arrays of contiguous entity handles. Dense tags are more
                efficient in both storage and memory if large numbers of
                entities are assigned the same tag. Storage for a given dense
                tag is not allocated until a tag value is set on an entity, at
                which point memory allocation for the dense tag occurs for all
                entities.
        """

        super(NativeMeshTag, self).__init__(mesh=mesh, name=name, doc=doc)

        if mesh is None or name is None:
            self._lazy_args["size"] = size
            self._lazy_args["dtype"] = dtype
            self._lazy_args["default"] = default
            return
        self.size = size
        self.dtype = dtype
        self.pymbtype = types.pymoab_data_type(self.dtype)
        self.default = default
        if storage_type is None:
            self.storage_type = types.MB_TAG_DENSE
        else:
            self.storage_type = storage_type
        # if the tag already exists, pick up its properties
        try:
            self.tag = self.mesh.mesh.tag_get_handle(self.name)
            self.size = self.tag.get_length()
            self.dtype = self.tag.get_dtype()
            self.pymbtype = types.pymoab_data_type(self.dtype)
            self.default = self.tag.get_default_value()
        except RuntimeError:
            self.tag = self.mesh.mesh.tag_get_handle(
                self.name,
                self.size,
                self.pymbtype,
                self.storage_type,
                create_if_missing=True,
                default_value=default,
            )
            if default is not None:
                self[:] = default

    def __delete__(self, mesh):
        super(NativeMeshTag, self).__delete__(mesh)
        mesh.mesh.tag_delete(self.tag)

    def delete(self, mesh=None):
        if mesh == None:
            mesh = self.mesh
        super(NativeMeshTag, self).delete()
        mesh.mesh.tag_delete(self.tag)

    def _collect_iterables(self, key, miter):
        ves = list(miter)
        list_of_ves = []
        # support either indexes or entityhandles
        for k in key:
            if isinstance(k, _eh_py_type):
                list_of_ves.append(k)
            elif isinstance(k, _INTEGRAL_TYPES):
                list_of_ves.append(ves[k])
            else:
                raise TypeError(
                    "{0} contains invalid element references "
                    "(non-ints, non-handles)".format(key)
                )
        return list_of_ves

    def __getitem__(self, key):
        m = self.mesh.mesh
        size = len(self.mesh)
        mtag = self.tag
        miter = self.mesh.iter_ve()
        # special case, get data on mesh's root set
        if isinstance(key, Mesh) and key == self.mesh:
            return self.mesh.mesh.tag_get_data(
                self.tag, self.mesh.mesh.get_root_set(), flat=True
            )
        elif isinstance(key, _eh_py_type):
            return self.mesh.mesh.tag_get_data(self.tag, key, flat=True)
        elif isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError(
                    "key index {0} greater than the size of the "
                    "mesh {1}".format(key, size)
                )
            for i_ve in zip(range(key + 1), miter):
                pass
            return self.mesh.mesh.tag_get_data(self.tag, i_ve[1], flat=True)
        elif isinstance(key, slice):
            flat = True if self.size == 1 else False
            ents = list(miter)[key]
            data = self.mesh.mesh.tag_get_data(self.tag, ents, flat=flat)
            return data
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            return self.mesh.mesh.tag_get_data(
                self.tag, [ve for b, ve in zip(key, miter) if b], flat=True
            )
        elif isinstance(key, Iterable):
            ves_to_get = self._collect_iterables(key, miter)
            return self.mesh.mesh.tag_get_data(self.tag, ves_to_get, flat=True)
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __setitem__(self, key, value):
        # get value into canonical form
        tsize = self.size
        value = np.asarray(value, self.tag.get_dtype())
        value = np.atleast_1d(value) if tsize == 1 else np.atleast_2d(value)
        # set up mesh to be iterated over
        m = self.mesh.mesh
        msize = len(self.mesh)
        mtag = self.tag
        miter = self.mesh.iter_ve()
        # special case: tag the mesh's root set
        if isinstance(key, Mesh) and key == self.mesh:
            self.mesh.mesh.tag_set_data(self.tag, self.mesh.mesh.get_root_set(), value)
        elif isinstance(key, _eh_py_type):
            self.mesh.mesh.tag_set_data(self.tag, key, value)
        elif isinstance(key, _INTEGRAL_TYPES):
            if key >= msize:
                raise IndexError(
                    "key index {0} greater than the size of the "
                    "mesh {1}".format(key, msize)
                )
            for i_ve in zip(range(key + 1), miter):
                pass
            self.mesh.mesh.tag_set_data(
                self.tag, i_ve[1], value if tsize == 1 else value[0]
            )
        elif isinstance(key, slice):
            key = list(miter)[key]
            v = np.empty((len(key), tsize), self.tag.get_dtype())
            if tsize == 1 and len(value.shape) == 1:
                v.shape = (len(key),)
            v[...] = value
            self.mesh.mesh.tag_set_data(mtag, key, v)
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != msize:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            key = [ve for b, ve in zip(key, miter) if b]
            v = np.empty((len(key), tsize), self.tag.get_dtype())
            if tsize == 1 and len(value.shape) == 1:
                v.shape = (len(key),)
            v[...] = value
            self.mesh.mesh.tag_set_data(mtag, key, v)
        elif isinstance(key, Iterable):
            if tsize != 1 and len(value) != len(key):
                v = np.empty((len(key), tsize), self.tag.get_dtype())
                v[...] = value
                value = v
            ves_to_tag = self._collect_iterables(key, miter)
            self.mesh.mesh.tag_set_data(mtag, ves_to_tag, value)
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __delitem__(self, key):
        m = self.mesh.mesh
        size = len(self.mesh)
        mtag = self.tag
        miter = self.mesh.iter_ve()
        # special case, look up mesh's root set
        if isinstance(key, Mesh) and key == self.mesh:
            self.mesh.mesh.tag_delete_data(mtag, self.mesh.mesh.get_root_set())
        elif isinstance(key, _eh_py_type):
            self.mesh.mesh.tag_delete_data(mtag, key)
        elif isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError(
                    "key index {0} greater than the size of the "
                    "mesh {1}".format(key, size)
                )
            for i_ve in zip(range(key + 1), miter):
                pass
            self.mesh.mesh.tag_delete_data(mtag, i_ve[1])
        elif isinstance(key, slice):
            self.mesh.mesh.tag_delete_data(mtag, list(miter)[key])
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the " "length of the mesh.")
            self.mesh.mesh.tag_delete_data(mtag, [ve for b, ve in zip(key, miter) if b])
        elif isinstance(key, Iterable):
            ves_to_del = self._collect_iterables(key, miter)
            self.mesh.mesh.tag_delete_data(mtag, ves_to_del)
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def expand(self):
        """This function creates a group of scalar tags from a vector tag. For
        a vector tag named <tag_name> of length N, scalar tags in the form:

        <tag_name>_000, <tag_name>_001, <tag_name>_002... <tag_name>_N

        are created and the data is tagged accordingly.
        """

        if self.size < 2:
            raise TypeError("Cannot expand a tag that is already a scalar.")
        if self.storage_type == types.MB_TAG_SPARSE:
            raise TypeError("Expansion of sparse tags is not implemented.")
        for j in range(self.size):
            data = [x[j] for x in self[:]]
            tag = self.mesh.mesh.tag_get_handle(
                "{0}_{1:03d}".format(self.name, j),
                1,
                self.pymbtype,
                storage_type=types.MB_TAG_DENSE,
                create_if_missing=True,
            )
            self.mesh.mesh.tag_set_data(tag, list(self.mesh.iter_ve()), data)

    def __add__(self, addend):
        """Addition operator for NativeMeshTag. Adds self to the addend using
        NumPy slicing if applicable. Addend can be any of the following types:
        ndarray, list, NativeMeshTag. Shapes must be the same for addition to work."""

        return self._do_op(addend, "+")

    def __sub__(self, subtrahend):
        """Subtraction operator for NativeMeshTag. Subtracts subtrahend from self using
        NumPy slicing, if applicable. Subtrahend can be any of the following types:
        ndarray, list, NativeMeshTag. Shapes must be the same for subtraction to work."""

        return self._do_op(subtrahend, "-")

    def __mul__(self, multiplier):
        """Multiplication operator for NativeMeshTag. Multiplies self by the multiplier using
        NumPy slicing if applicable. Multiplier can be any of the following types: int, float,
        ndarray, list, NativeMeshTag. Throws error if shapes are incorrect"""

        return self._do_op(multiplier, "*")

    def __truediv__(self, divisor):
        """Division operator for NativeMeshTag. Divides self by the divisor using
        NumPy slicing, if applicable. Divisor can be any of the following types: int, float,
        ndarray, list, NativeMeshTag. Throws error if shapes are incorrect."""

        return self._do_op(divisor, "/")

    def _do_op(self, other, op):
        """Helper function for NativeMeshTag operators"""
        if isinstance(other, NativeMeshTag):
            return self._do_native_op(other, op)
        elif isinstance(other, list) or isinstance(other, np.ndarray):
            return self._do_list_op(other, op)
        elif isinstance(other, int) or isinstance(other, float):
            return self._do_scalar_op(other, op)
        else:
            raise TypeError("Incorrect operand type provided")

    def _do_native_op(self, other, op):
        """Helper function for the NativeMeshTag operators when the operation is
        being done on two NativeMeshTag objects"""

        if self.size == other.size:
            try:
                return _ops[op](self[:], other[:])
            except:
                raise ValueError("Data must have the same shape")
        elif other.size > 1:
            try:
                return _ops[op](self[:][:, None], other[:])
            except:
                raise ValueError("Incompatible shape for scalar or vector data")
        elif self.size > 1:
            try:
                return _ops[op](self[:], other[:][None, :])
            except:
                raise ValueError("Incompatible shape for vector or scalar data")

    def _do_list_op(self, other, op):
        """Helper function for the NativeMeshTag operators when the operation is
        being done on a list or ndarray object"""

        try:
            return _ops[op](self[:], other)
        except:
            raise ValueError("Incompatible array or list shape")

    def _do_scalar_op(self, other, op):
        """Helper function for the NativeMeshTag operators when the operation is
        being done on a float or int"""

        return _ops[op](self[:], other)


class ComputedTag(Tag):
    '''A mesh tag which looks itself up by calling a function (or other callable)
    with the following signature::

        def f(mesh, i):
            """mesh is a pyne.mesh.Mesh() object and i is the volume element
            index to compute.
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
            self._lazy_args["f"] = f
            return
        self.f = f

    def __getitem__(self, key):
        m = self.mesh
        f = self.f
        size = len(m)
        if isinstance(key, _INTEGRAL_TYPES):
            if key >= size:
                raise IndexError(
                    "key index {0} greater than the size of the "
                    "mesh {1}".format(key, size)
                )
            return f(m, key)
        elif isinstance(key, slice):
            return [f(m, i) for i in range(*key.indices(size))]
        elif isinstance(key, np.ndarray) and key.dtype == np.bool:
            if len(key) != size:
                raise KeyError("boolean mask must match the length " "of the mesh.")
            return [f(m, i) for i, b in enumerate(key) if b]
        elif isinstance(key, Iterable):
            return [f(m, i) for i in key]
        else:
            raise TypeError(
                "{0} is not an int, slice, mask, " "or fancy index.".format(key)
            )

    def __setitem__(self, key, value):
        msg = "the computed tag {0!r} may not be set".format(self.name)
        raise AttributeError(msg)

    def __delitem__(self, key):
        msg = "the computed tag {0!r} may not be deleted".format(self.name)
        raise AttributeError(msg)


class MeshError(Exception):
    """Errors related to instantiating mesh objects and utilizing their methods."""

    pass


class Mesh(object):
    """This class houses a PyMOAB core instance and contains methods for various
    mesh operations. Special methods exploit the properties of structured mesh.

    Attributes
    ----------
    mesh : PyMOAB core instance
    structured : bool
        True for structured mesh.
    structured_coords : list of lists
        A list containing lists of x_points, y_points and z_points that make up
        a structured mesh.
    structured_ordering : str
        A three character string denoting the iteration order of the mesh (e.g.
        'xyz', meaning z changest fastest, then y, then x.)
    """

    def __init__(
        self,
        mesh=None,
        structured=False,
        structured_coords=None,
        structured_set=None,
        structured_ordering="xyz",
        mats=(),
    ):
        """Parameters
        ----------
        mesh : PyMOAB core instance or str, optional
            Either a PyMOAB core instance or a file name of a PyMOAB mesh file.
        structured : bool, optional
            True for structured mesh.
        structured_coords : list of lists, optional
            A list containing lists of x_points, y_points and z_points
            that make up a structured mesh.
        structured_set : PyMOAB entity set handle, optional
            A preexisting structured entity set on an PyMOAB core instance with
            a "BOX_DIMS" tag.
        structured_ordering : str, optional
            A three character string denoting the iteration order of the mesh
            (e.g. 'xyz', meaning z changest fastest, then y, then x.)
        mats : MaterialLibrary or dict or Materials or None, optional
            This is a mapping of volume element handles to Material objects.
            If mats is None, then no empty materials are created for the mesh.

            Unstructured mesh instantiation:
                 - From PyMOAB core instance by specifying: <mesh>
                 - From mesh file by specifying: <mesh_file>

            Structured mesh instantiation:
                - From PyMOAB core instance with exactly 1 entity set
                  (with BOX_DIMS tag) by specifying <mesh> and structured = True.
                - From mesh file with exactly 1 entity set (with BOX_DIMS tag)
                  by specifying <mesh_file> and structured = True.
                - From a PyMOAB instance with multiple entity sets by
                  specifying <mesh>, <structured_set>, structured=True.
                - From coordinates by specifying <structured_coords>,
                  structured=True, and optional pre-existing PyMOAB core
                  instance <mesh>

            The "BOX_DIMS" tag on PyMOAB core instances containing structured
            mesh is a vector of floats in the following form:
            [i_min, j_min, k_min, i_max, j_max, k_max]
            where each value is a volume element index number. Typically volume
            elements should be indexed from 0. The "BOX_DIMS" information is
            stored in self.dims.
        """

        # if Mesh is made and no parameters passed, raise MeshError
        if (
            (mesh is None)
            and (not structured)
            and (structured_coords is None)
            and (structured_set is None)
            and (structured_ordering == "xyz")
            and (mats == ())
        ):
            raise MeshError(
                "Trivial mesh instantiation detected. "
                "For structured mesh instantiation, "
                "supply exactly one of the following:\n"
                "A. PyMOAB instance\n"
                "B. Mesh file\n"
                "C. Mesh coordinates\n"
                "D. Structured entity set AND PyMOAB instance"
            )
        if mesh is None:
            self.mesh = mb_core.Core()
        elif isinstance(mesh, basestring):
            self.mesh = mb_core.Core()
            self.mesh.load_file(mesh)
        else:
            self.mesh = mesh

        self.structured = structured

        if self.structured:
            self.scd = scd.ScdInterface(self.mesh)
            self.structured_coords = structured_coords
            self.structured_ordering = structured_ordering
            # if a MOAB mesh instance exists and no structured coords
            # or structured set is provided, search for a single
            # structured set
            if (mesh is not None) and not structured_coords and not structured_set:
                # check for the structured box tag on the instance
                try:
                    box_tag = self.mesh.tag_get_handle(_BOX_DIMS_TAG_NAME)
                except RuntimeError as e:
                    print("BOX_DIMS not found on MOAB mesh instance")
                    raise e

                # find all entity sets with the structured box tag
                count = 0
                root_set = self.mesh.get_root_set()
                for ent_set in self.mesh.get_entities_by_type(
                    root_set, types.MBENTITYSET
                ):
                    try:
                        self.mesh.tag_get_data(box_tag, ent_set)
                    except RuntimeError:
                        pass
                    else:
                        self.structured_set = ent_set
                        count += 1

                if count == 0:
                    raise MeshError(
                        "Found no structured meshes in " "file {0}".format(mesh)
                    )
                elif count > 1:
                    raise MeshError(
                        "Found {0} structured meshes."
                        " Instantiate individually using"
                        " from_ent_set()".format(count)
                    )

            # from coordinates
            elif (mesh is None) and structured_coords and not structured_set:
                # check for single vertex coordinates here? it seems we only support volumetric mesh -PCS
                extents = [0, 0, 0] + [len(x) - 1 for x in structured_coords]
                low = hcoord.HomCoord([0, 0, 0])
                high = hcoord.HomCoord([len(x) - 1 for x in structured_coords])
                # get coordinates as array
                xs = np.asarray(structured_coords[0])
                ys = np.asarray(structured_coords[1])
                zs = np.asarray(structured_coords[2])
                # generate array
                coords = np.empty((xs.size * ys.size * zs.size * 3,), dtype=np.float64)
                # set mesh values
                coords[0::3] = np.tile(xs, ys.size * zs.size)
                coords[1::3] = np.tile(
                    np.repeat(
                        ys,
                        xs.size,
                    ),
                    zs.size,
                )
                coords[2::3] = np.repeat(zs, xs.size * ys.size)
                # construct the structured mesh
                scd_box = self.scd.construct_box(low, high, coords)
                self.structured_set = scd_box.box_set()

            # from mesh and structured_set:
            elif not structured_coords and structured_set:
                # check for the structured box tag on the instance
                try:
                    box_tag = self.mesh.tag_get_handle(_BOX_DIMS_TAG_NAME)
                except types.MB_TAG_NOT_FOUND as e:
                    print("BOX_DIMS not found on MOAB mesh instance")
                    raise e
                # check that the structured_set found is tagged as a structured set
                try:
                    self.mesh.tag_get_data(box_tag, structured_set)
                except:
                    print("Supplied entity set does not contain BOX_DIMS tag")
                    raise e

                self.structured_set = structured_set
            else:
                raise MeshError(
                    "For structured mesh instantiation, need to"
                    "supply exactly one of the following:\n"
                    "A. PyMOAB instance\n"
                    "B. Mesh file\n"
                    "C. Mesh coordinates\n"
                    "D. Structured entity set AND PyMOAB instance"
                )

            self.dims = list(
                self.mesh.tag_get_data(
                    self.mesh.tag_get_handle(_BOX_DIMS_TAG_NAME),
                    self.structured_set,
                    flat=True,
                )
            )

            self.vertex_dims = list(self.dims[0:3]) + [x + 1 for x in self.dims[3:6]]

            if self.structured_coords is None:
                self.structured_coords = [
                    self.structured_get_divisions("x"),
                    self.structured_get_divisions("y"),
                    self.structured_get_divisions("z"),
                ]
        else:
            # Unstructured mesh cases
            # Error if structured arguments are passed
            if structured_coords or structured_set:
                MeshError(
                    "Structured mesh arguments should not be present for\
                            unstructured Mesh instantiation."
                )

        # sets mats
        mats_in_mesh_file = False
        if isinstance(mesh, basestring) and len(mats) == 0:
            with tb.open_file(mesh) as h5f:
                if "/mat_name" in h5f:
                    mats_in_mesh_file = True
                    mat_path = "/mat_name"
                elif ("/materials" in h5f) or ("/material_library/materials" in h5f):
                    mats_in_mesh_file = True
                    mat_path = "/materials"

            if mats_in_mesh_file:
                mats = MaterialLibrary(mesh, datapath=mat_path)

        if mats is None:
            pass
        elif len(mats) == 0 and not mats_in_mesh_file:
            mats = MaterialLibrary()
        elif not isinstance(mats, MaterialLibrary):
            mats = MaterialLibrary(mats)

        self.mats = mats

        # tag with volume id and ensure mats exist.
        ves = list(self.iter_ve())
        tags = self.mesh.tag_get_tags_on_entity(ves[0])
        tags = set(tag.get_name() for tag in tags)
        if "idx" in tags:
            tag_idx = self.mesh.tag_get_handle("idx")
        else:
            tag_idx = self.mesh.tag_get_handle(
                "idx",
                1,
                types.MB_TYPE_INTEGER,
                types.MB_TAG_DENSE,
                create_if_missing=True,
            )
        # tag elements with index
        idxs = np.arange(0, len(ves))
        self.mesh.tag_set_data(tag_idx, ves, idxs)
        # check for and populate materials
        if mats is not None:
            for i in range(len(ves)):
                if i not in mats:
                    mats[i] = Material()
        self._len = len(ves)

        # Default tags
        self.tags = {}
        if mats is not None:
            # metadata tags, these should come first so they don't accidentally
            # overwite hard coded tag names.
            metatagnames = set()
            for mat in mats.values():
                metatagnames.update(mat.metadata.keys())
            for name in metatagnames:
                setattr(self, name, MetadataTag(mesh=self, name=name))
        # PyMOAB tags
        tagnames = set()
        for ve in ves:
            tagnames.update(t.get_name() for t in self.mesh.tag_get_tags_on_entity(ve))
        for name in tagnames:
            setattr(self, name, NativeMeshTag(mesh=self, name=name))

        if mats is not None:
            # Material property tags
            self.atoms_per_molecule = MaterialPropertyTag(
                mesh=self, name="atoms_per_molecule", doc="Number of atoms per molecule"
            )
            self.metadata = MaterialPropertyTag(
                mesh=self,
                name="metadata",
                doc="metadata attributes, stored on the material",
            )
            self.comp = MaterialPropertyTag(
                mesh=self,
                name="comp",
                doc="normalized composition mapping from nuclides to " "mass fractions",
            )
            self.mass = MaterialPropertyTag(
                mesh=self, name="mass", doc="the mass of the material"
            )
            self.density = MaterialPropertyTag(
                mesh=self, name="density", doc="the density [g/cc]"
            )
            # Material method tags
            methtagnames = (
                "expand_elements",
                "mass_density",
                "molecular_mass",
                "mult_by_mass",
                "number_density",
                "sub_act",
                "sub_fp",
                "sub_lan",
                "sub_ma",
                "sub_tru",
                "to_atom_frac",
            )
            for name in methtagnames:
                doc = "see Material.{0}() for more information".format(name)
                setattr(self, name, MaterialMethodTag(mesh=self, name=name, doc=doc))

    def get_all_tags(self):
        return [
            getattr(self, t) for t in dir(self) if isinstance(getattr(self, t), Tag)
        ]

    def __len__(self):
        return self._len

    def __iter__(self):
        """Iterates through the mesh and at each step yield the volume element
        index i, the material mat, and the volume element itself ve.
        """

        mats = self.mats
        if mats is None:
            for i, ve in enumerate(self.iter_ve()):
                yield i, None, ve
        else:
            for i, ve in enumerate(self.iter_ve()):
                yield i, mats[i], ve

    def iter_ve(self):
        """Returns an iterator that yields on the volume elements."""
        if self.structured:
            return self.structured_iterate_hex(self.structured_ordering)
        else:
            return iter(
                list(
                    self.mesh.get_entities_by_dimension(
                        self.mesh.get_root_set(), 3, True
                    )
                )
            )

    def __contains__(self, i):
        return i < len(self)

    def __setattr__(self, name, value):
        if isinstance(value, Tag) and hasattr(value, "_lazy_args"):
            # some 1337 1Azy 3\/a1
            kwargs = value._lazy_args
            kwargs["mesh"] = self if kwargs["mesh"] is None else kwargs["mesh"]
            kwargs["name"] = name if kwargs["name"] is None else kwargs["name"]
            value = type(value)(**kwargs)
        super(Mesh, self).__setattr__(name, value)

    def tag(
        self,
        name,
        value=None,
        tagtype=None,
        doc=None,
        size=None,
        dtype=None,
        storage_type=None,
    ):
        """Adds a new tag to the mesh, guessing the approriate place to store
        the data.

        Parameters
        ----------
        name : str
            The tag name
        value : optional
            The value to initialize the tag with, skipped if None.
        tagtype : Tag or str, optional
            The type of tag this should be any of the following classes or
            strings are accepted: NativeMeshTag, MetadataTag, ComputedTag,
            'nat_mesh', 'metadata', or 'computed'.
        doc : str, optional
            The tag documentation string.
        size : int, optional
            The size of the tag. This only applies to NativeMeshTags.
        dtype : numpy dtype, optional
            The data type of the tag. This only applies to NativeMeshTags. See PyMOAB
            for more details.
        storage_type: str, optional
            MOAB tag storage type (MB_TAG_DENSE, MB_TAG_SPARSE, etc.)
            in advanced use of the database, this flag controls how this tag's
            data is stored in memory, the supported storage types are:
            sparse - MB_TAG_SPARSE tag, values are stored as a list of (entity
                handle, tag value) tuples, one list per sparse tag, sorted by
                entity handle.
            dense - MB_TAG_DENSE tag, values are stored in arrays which match
                arrays of contiguous entity handles. Dense tags are more
                efficient in both storage and memory if large numbers of
                entities are assigned the same tag. Storage for a given dense
                tag is not allocated until a tag value is set on an entity, at
                which point memory allocation for the dense tag occurs for all
                entities.
        """

        if name in self.tags:
            raise KeyError("{0} tag already exists on the mesh".format(name))
        if tagtype is None:
            if callable(value):
                tagtype = ComputedTag
            elif size is None and dtype is not None:
                size = 1
                tagtype = NativeMeshTag
            elif size is not None and dtype is None:
                dtype = "f8"
                tagtype = NativeMeshTag
            elif value is None:
                size = 1
                value = 0.0
                dtype = "f8"
                tagtype = NativeMeshTag
            elif isinstance(value, float):
                size = 1
                dtype = "f8"
                tagtype = NativeMeshTag
            elif isinstance(value, int):
                size = 1
                dtype = "i"
                tagtype = NativeMeshTag
            elif isinstance(value, str):
                tagtype = MetadataTag
            elif isinstance(value, _SEQUENCE_TYPES):
                raise ValueError(
                    "ambiguous tag {0!r} creation when value is a"
                    " sequence, please set tagtype, size, "
                    "or dtype".format(name)
                )
            else:
                tagtype = MetadataTag

        if tagtype is NativeMeshTag or tagtype.lower() == "nat_mesh":
            if storage_type is None or storage_type.lower() == "dense":
                storage_type = types.MB_TAG_DENSE
                default_value = 0.0
            elif storage_type.lower() == "sparse":
                storage_type = types.MB_TAG_SPARSE
                default_value = None
            else:
                raise ValueError("storage_type {0} not valid".format(storage_type))
            t = NativeMeshTag(
                size=size,
                dtype=dtype,
                mesh=self,
                name=name,
                doc=doc,
                storage_type=storage_type,
                default=default_value,
            )
        elif storage_type is not None:
            raise ValueError("storage_type works for only NativeMeshTag")
        elif tagtype is MetadataTag or tagtype.lower() == "metadata":
            t = MetadataTag(mesh=self, name=name, doc=doc)
        elif tagtype is ComputedTag or tagtype.lower() == "computed":
            t = ComputedTag(f=value, mesh=self, name=name, doc=doc)
        else:
            raise ValueError("tagtype {0} not valid".format(tagtype))
        if value is not None and tagtype is not ComputedTag:
            t[:] = value

        setattr(self, name, t)

    def get_tag(self, tag_name):
        return getattr(self, tag_name)

    def delete_tag(self, tag):
        if isinstance(tag, Tag):
            tag_name = tag.name
            tag_handle = tag
        elif isinstance(tag, str):
            tag_name = tag
            tag_handle = self.tags[tag_name]
        else:
            raise ValueError("{0} is neither a Tag object nor a string".format(tag))

        tag_handle.delete()

        if hasattr(self, tag_name):
            delattr(self, tag_name)

    def __iadd__(self, other):
        """Adds the common tags of other to the mesh object."""
        tags = self.common_ve_tags(other)
        return self._do_op(other, tags, "+")

    def __isub__(self, other):
        """Substracts the common tags of other to the mesh object."""
        tags = self.common_ve_tags(other)
        return self._do_op(other, tags, "-")

    def __imul__(self, other):
        """Multiplies the common tags of other to the mesh object."""
        tags = self.common_ve_tags(other)
        return self._do_op(other, tags, "*")

    def __idiv__(self, other):
        """Divides the common tags of other to the mesh object."""
        tags = self.common_ve_tags(other)
        return self._do_op(other, tags, "/")

    def __itruediv__(self, other):
        """Divides the common tags of other to the mesh object."""
        tags = self.common_ve_tags(other)
        return self._do_op(other, tags, "/")

    def _do_op(self, other, tags, op, in_place=True):
        """Private function to do mesh +, -, *, /."""
        # Exclude error tags in a case a StatMesh is mistakenly initialized as
        # a Mesh object.
        tag_list = []

        if isinstance(tags, str):
            tag_list = [tags]
        else:
            tag_list = tags

        tag_list = set(tag for tag in tag_list if not tag.endswith("_error"))

        if in_place:
            mesh_1 = self
        else:
            mesh_1 = copy.copy(self)
        for tag in tag_list:
            for ve_1, ve_2 in zip(
                zip(
                    iter(
                        meshset_iterate(
                            mesh_1.mesh, mesh_1.structured_set, types.MBMAXTYPE, dim=3
                        )
                    )
                ),
                zip(
                    iter(
                        meshset_iterate(
                            other.mesh, other.structured_set, types.MBMAXTYPE, dim=3
                        )
                    )
                ),
            ):
                mesh_1_tag = mesh_1.mesh.tag_get_handle(tag)
                other_tag = other.mesh.tag_get_handle(tag)
                val = _ops[op](
                    mesh_1.mesh.tag_get_data(mesh_1_tag, ve_1, flat=True),
                    other.mesh.tag_get_data(other_tag, ve_2, flat=True),
                )
                # removed the [0] from the ops call
                mesh_1.mesh.tag_set_data(mesh_1_tag, ve_1, val)
        return mesh_1

    def common_ve_tags(self, other):
        """Returns the volume element tags in common between self and other."""
        self_it = MeshSetIterator(
            self.mesh, self.structured_set, types.MBMAXTYPE, dim=3
        )
        self_tags = self.mesh.tag_get_tags_on_entity(next(self_it))
        other_it = MeshSetIterator(
            other.mesh, other.structured_set, types.MBMAXTYPE, dim=3
        )
        other_tags = other.mesh.tag_get_tags_on_entity(next(other_it))
        self_tags = set(x.get_name() for x in self_tags)
        other_tags = set(x.get_name() for x in other_tags)
        intersect = self_tags & other_tags
        intersect.discard("idx")
        return intersect

    def __copy__(self):
        # first copy full pymoab instance
        pymb_copy = mb_core.Core()

        # now create Mesh objected from copied PyMOAB instance
        mesh_copy = Mesh(mesh=pymb_copy, structured=copy.copy(self.structured))
        return mesh_copy

    # Non-structured volume methods
    def elem_volume(self, ve):
        """Get the volume of a hexahedral or tetrahedral volume element

        Approaches are adapted from MOAB's measure.cpp.

        Parameters
        ----------
        ve : PyMOAB EntitySet handle
            A volume element

        Returns
        -------
        .. : float
            Element's volume. Returns None if volume is not a hex or tet.
        """

        coord = self.mesh.get_coords(self.mesh.get_connectivity(ve)).reshape(-1, 3)
        num_coords = coord.shape[0]

        if num_coords == 4:
            return abs(np.linalg.det(coord[:-1] - coord[1:])) / 6.0
        elif num_coords == 8:
            b = coord[
                np.array(
                    [
                        [0, 1, 3, 4],
                        [7, 3, 6, 4],
                        [4, 5, 1, 6],
                        [1, 6, 3, 4],
                        [2, 6, 3, 1],
                    ]
                )
            ]
            return np.sum(np.abs(np.linalg.det(b[:, :-1] - b[:, 1:]))) / 6.0
        else:
            return None

    def ve_center(self, ve):
        """Finds the point at the center of any tetrahedral or hexahedral mesh
        volume element.

        Parameters
        ----------
        ve : PyMOAB EntitySet handle
           Any mesh volume element.

        Returns
        -------
        center : tuple
           The (x, y, z) coordinates of the center of the mesh volume element.
        """

        ve_handle = _eh_py_type(ve)
        coords = self.mesh.get_coords(self.mesh.get_connectivity(ve_handle)).reshape(
            -1, 3
        )
        center = tuple([np.mean(coords[:, x]) for x in range(3)])
        return center

    # Structured methods:
    def structured_get_vertex(self, i, j, k):
        """Return the handle for (i,j,k)'th vertex in the mesh"""
        self._structured_check()
        n = _structured_find_idx(self.vertex_dims, (i, j, k))
        return _structured_step_iter(
            meshset_iterate(self.mesh, self.structured_set, entity_type=types.MBVERTEX),
            n,
        )

    def structured_get_hex(self, i, j, k):
        """Return the handle for the (i,j,k)'th hexahedron in the mesh"""
        self._structured_check()
        n = _structured_find_idx(self.dims, (i, j, k))
        return _structured_step_iter(
            meshset_iterate(self.mesh, self.structured_set, types.MBHEX, 3), n
        )

    def structured_hex_volume(self, i, j, k):
        """Return the volume of the (i,j,k)'th hexahedron in the mesh"""
        self._structured_check()
        handle = self.structured_get_hex(i, j, k)
        h = self.mesh.get_connectivity(handle)
        coord = self.mesh.get_coords(list(h))
        coord = coord.reshape(8, 3)
        # assumes a "well-behaved" hex element
        dx = max(coord[:, 0]) - min(coord[:, 0])
        dy = max(coord[:, 1]) - min(coord[:, 1])
        dz = max(coord[:, 2]) - min(coord[:, 2])
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

          structured_iterate_hex(): equivalent to mehset_iterator over hexes
                                    in mesh
          structured_iterate_hex("xyz"): iterate over entire mesh, with
                                         k-coordinates changing fastest,
                                         i-coordinates least fast.
          structured_iterate_hex("yz", x=3): Iterate over the j-k plane of the
                                             mesh whose i-coordinate is 3, with
                                             k values changing fastest.
          structured_iterate_hex("z"): Iterate over k-coordinates, with
                                       i=dims.imin and j=dims.jmin
          structured_iterate_hex("yxz", y=(3,4)): Iterate over all hexes with
                                        j-coordinate = 3 or 4.  k-coordinate
                                        values change fastest, j-values least
                                        fast.
        """

        self._structured_check()

        # special case: zyx order is the standard pytaps iteration order,
        # so we can save time by simply returning an iterator
        # if no kwargs were specified
        if order == "zyx" and not kw:
            return meshset_iterate(
                self.mesh, self.structured_set, entity_type=types.MBHEX, dim=3
            )

        indices, ordmap = _structured_iter_setup(self.dims, order, **kw)
        return _structured_iter(
            indices,
            ordmap,
            self.dims,
            meshset_iterate(
                self.mesh, self.structured_set, entity_type=types.MBHEX, dim=3
            ),
        )

    def structured_iterate_vertex(self, order="zyx", **kw):
        """Get an iterator over the vertices of the mesh

        See structured_iterate_hex() for an explanation of the order argument
        and the available keyword arguments.
        """

        self._structured_check()
        # special case: zyx order without kw is equivalent to an iterator
        if order == "zyx" and not kw:
            return meshset_iterate(
                self.mesh, self.structured_set, entity_type=types.MBVERTEX
            )

        indices, ordmap = _structured_iter_setup(self.vertex_dims, order, **kw)
        return _structured_iter(
            indices,
            ordmap,
            self.vertex_dims,
            meshset_iterate(self.mesh, self.structured_set, entity_type=types.MBVERTEX),
        )

    def structured_iterate_hex_volumes(self, order="zyx", **kw):
        """Get an iterator over the volumes of the mesh hexahedra

        See structured_iterate_hex() for an explanation of the order argument
        and the available keyword arguments.
        """

        self._structured_check()
        indices, _ = _structured_iter_setup(self.dims, order, **kw)
        # Use an inefficient but simple approach: call structured_hex_volume()
        # on each required i,j,k pair.
        # A better implementation would only make one call to getVtxCoords.
        for A in itertools.product(*indices):
            # the ordmap returned from _structured_iter_setup maps to kji/zyx
            # ordering, but we want ijk/xyz ordering, so create the ordmap
            # differently.
            ordmap = [order.find(L) for L in "xyz"]
            ijk = [A[ordmap[x]] for x in range(3)]
            yield self.structured_hex_volume(*ijk)

    def iter_structured_idx(self, order=None):
        """Return an iterater object of volume element indexes (idx) for any
        iteration order. Note that idx is assigned upon instantiation in the
        order of the structured_ordering attribute. This method is meant to be
        used when the order argument is different from structured_ordering.
        When they are the same, the iterator (0, 1, 2, ... N-1) is returned.

        Parameters
        ----------
        order : str, optional
            The requested iteration order (e.g. 'zyx').
        """

        self._structured_check()
        if not order:
            order = self.structured_ordering

        ves = self.structured_iterate_hex(order)
        tag = self.mesh.tag_get_handle("idx")
        for ve in ves:
            yield self.mesh.tag_get_data(tag, ve, flat=True)[0]

    def structured_get_divisions(self, dim):
        """Get the mesh divisions on a given dimension

        Given a dimension "x", "y", or "z", return a list of the mesh vertices
        along that dimension.
        """

        self._structured_check()

        ## sometimes the dim is the ascii of the 'x', 'y', 'z'
        if len(dim) == 1 and dim in "xyz":
            idx = "xyz".find(dim)
            return [
                self.mesh.get_coords(v)[idx]
                for v in self.structured_iterate_vertex(dim)
            ]

        else:
            raise MeshError("Invalid dimension: {0}".format(str(dim)))

    def _structured_check(self):
        if not self.structured:
            raise MeshError(
                "Structured mesh methods cannot be called from "
                "unstructured mesh instances."
            )

    def write_hdf5(self, filename, write_mats=True):
        """Writes the mesh to an hdf5 file."""
        self.mesh.write_file(filename)
        if write_mats and self.mats is not None:
            self.mats.write_hdf5(filename)

    def cell_fracs_to_mats(self, cell_fracs, cell_mats):
        """This function uses the output from dagmc.discretize_geom() and
        a mapping of geometry cells to Materials to assign materials
        to each mesh volume element.

        Parameters
        ----------
        cell_fracs : structured array
            The output from dagmc.discretize_geom(). A sorted, one dimensional
            array, each entry containing the following fields:

                :idx: int
                    The volume element index.
                :cell: int
                    The geometry cell number.
                :vol_frac: float
                    The volume fraction of the cell withing the mesh ve.
                :rel_error: float
                    The relative error associated with the volume fraction.

            The array must be sorted with respect to both idx and cell, with
            cell changing fastest.
        cell_mats : dict
            Maps geometry cell numbers to Material objects that represent what
            material each cell is made of.
        """

        for i in range(len(self)):
            mat_col = {}  # Collection of materials in the ith ve.
            for row in cell_fracs[cell_fracs["idx"] == i]:
                mat_col[cell_mats[row["cell"]]] = row["vol_frac"]

            mixed = MultiMaterial(mat_col)
            self.mats[i] = mixed.mix_by_volume()

    def tag_cell_fracs(self, cell_fracs):
        """This function uses the output from dagmc.discretize_geom() and
        a mapping of geometry cells to set the cell_fracs_tag.

        Parameters
        ----------
        cell_fracs : structured array
            The output from dagmc.discretize_geom(). A sorted, one dimensional
            array, each entry containing the following fields:

                :idx: int
                    The volume element index.
                :cell: int
                    The geometry cell number.
                :vol_frac: float
                    The volume fraction of the cell withing the mesh ve.
                :rel_error: float
                    The relative error associated with the volume fraction.

            The array must be sorted with respect to both idx and cell, with
            cell changing fastest.
        """

        num_vol_elements = len(self)
        # sort cell_fracs
        cell_fracs = _cell_fracs_sort_vol_frac_reverse(cell_fracs)
        # Find the maximum cell number in a voxel
        max_num_cells = -1
        for i in range(num_vol_elements):
            max_num_cells = max(max_num_cells, len(cell_fracs[cell_fracs["idx"] == i]))

        # create tag frame with default value
        cell_largest_frac_number = [-1] * num_vol_elements
        cell_largest_frac = [0.0] * num_vol_elements
        voxel_cell_number = np.empty(shape=(num_vol_elements, max_num_cells), dtype=int)
        voxel_cell_fracs = np.empty(
            shape=(num_vol_elements, max_num_cells), dtype=float
        )
        voxel_cell_number.fill(-1)
        voxel_cell_fracs.fill(0.0)

        # set the data
        for i in range(num_vol_elements):
            for (cell, row) in enumerate(cell_fracs[cell_fracs["idx"] == i]):
                voxel_cell_number[i, cell] = row["cell"]
                voxel_cell_fracs[i, cell] = row["vol_frac"]
            # cell_largest_frac_tag
            cell_largest_frac[i] = max(voxel_cell_fracs[i, :])
            largest_index = list(voxel_cell_fracs[i, :]).index(cell_largest_frac[i])
            cell_largest_frac_number[i] = int(voxel_cell_number[i, largest_index])

        # create the tags
        self.tag(
            name="cell_number",
            value=voxel_cell_number,
            doc="cell numbers of the voxel, -1 used to fill vacancy",
            tagtype=NativeMeshTag,
            size=max_num_cells,
            dtype=int,
        )
        self.tag(
            name="cell_fracs",
            value=voxel_cell_fracs,
            tagtype=NativeMeshTag,
            doc="volume fractions of each cell in the "
            "voxel, 0.0 used to fill vacancy",
            size=max_num_cells,
            dtype=float,
        )
        self.tag(
            name="cell_largest_frac_number",
            value=cell_largest_frac_number,
            tagtype=NativeMeshTag,
            doc="cell number of the cell with largest volume fraction in " "the voxel",
            size=1,
            dtype=int,
        )
        self.tag(
            name="cell_largest_frac",
            value=cell_largest_frac,
            tagtype=NativeMeshTag,
            doc="cell fraction of the cell with " "largest cell volume fraction",
            size=1,
            dtype=float,
        )


class StatMesh(Mesh):
    """This class extends the basic Mesh class by modifying the standard
    mathematical operations that are performed on multiple meshes.

    A StatMesh assumes that each value being operated upon also has a
    statistical error associaed with it, and forces operations on the
    statistical error as well.  For any tag with name `tag_name` the StatMesh
    assumes that there is also a tag with name `tag_name_rel_error`.

    For example, when to quantities are added together, c = a + b, the
    statistical error of c is found by combining the statistical errors of a
    and b.
    """

    def __init__(
        self,
        mesh=None,
        structured=False,
        structured_coords=None,
        structured_set=None,
        mats=(),
    ):

        super(StatMesh, self).__init__(
            mesh=mesh,
            structured=structured,
            structured_coords=structured_coords,
            structured_set=structured_set,
            mats=mats,
        )

    def _do_op(self, other, tags, op, in_place=True):
        """Private function to do mesh +, -, *, /. Called by operater
        overloading functions.
        """

        # Exclude error tags because result and error tags are treated
        # simultaneously so there is not need to include both in the tag
        # list to iterate through.
        error_suffix = "_rel_error"

        tags = set(tag for tag in tags if not tag.endswith(error_suffix))

        if in_place:
            mesh_1 = self
        else:
            mesh_1 = copy.copy(self)

        for tag in tags:
            for ve_1, ve_2 in zip(
                zip(
                    iter(
                        meshset_iterate(
                            mesh_1.mesh, mesh_1.structured_set, types.MBMAXTYPE, dim=3
                        )
                    )
                ),
                zip(
                    iter(
                        meshset_iterate(
                            other.mesh, other.structured_set, types.MBMAXTYPE, dim=3
                        )
                    )
                ),
            ):

                mesh_1_err_tag = mesh_1.mesh.tag_get_handle(tag + error_suffix)
                other_err_tag = other.mesh.tag_get_handle(tag + error_suffix)
                mesh_1_tag = mesh_1.mesh.tag_get_handle(tag)
                other_tag = other.mesh.tag_get_handle(tag)

                mesh_1_val = mesh_1.mesh.tag_get_data(mesh_1_tag, ve_1, flat=True)[0]
                other_val = other.mesh.tag_get_data(other_tag, ve_2, flat=True)[0]
                mesh_1_err = mesh_1.mesh.tag_get_data(mesh_1_err_tag, ve_1, flat=True)[
                    0
                ]
                other_err = other.mesh.tag_get_data(other_err_tag, ve_2, flat=True)[0]

                new_err_val = err__ops[op](mesh_1_val, other_val, mesh_1_err, other_err)
                mesh_1.mesh.tag_set_data(mesh_1_err_tag, ve_1, new_err_val)

                new_val = _ops[op](mesh_1_val, other_val)
                mesh_1.mesh.tag_set_data(mesh_1_tag, ve_1, new_val)

        return mesh_1


class MeshTally(StatMesh):
    """This class stores all information from a single mesh tally that
    exists within some meshtal or state point file. Header information is
    stored as attributes and the "mesh" attribute is a PyNE mesh object tagged
    with all result and relative error data. This class inherits from StatMesh,
    exposing all statistical mesh manipulation methods.

    Attributes
    ----------
    tally_number : int
        The tally number.
        For mesh tally from MCNP, it must end with 4 (e.g. 4, 14, 214).
        For mesh tally from OpenMC, it could be any int.
    particle : string
        Either "neutron" for a neutron mesh tally or "photon" for a photon mesh
        tally.
    dose_response : bool
        True if the tally is modified by a dose response function.
    x_bounds : tuple of floats
        The locations of mesh vertices in the x direction.
    y_bounds : tuple of floats
        The locations of mesh vertices in the y direction.
    z_bounds : tuple of floats
        The locations of mesh vertices in the z direction.
    dims : list
        Dimensions of the mesh.
    num_ves : int
        Number of volume elements.
    e_bounds : tuple of floats
        The minimum and maximum bounds for energy bins.
    num_e_groups: int
        Number of energy groups.
    mesh :
        An PyMOAB core instance tagged with all results and
        relative errors.
    tag_names : iterable
        Four strs that specify the tag names for the results, relative errors,
        total results, and total relative error.

    Notes
    -----
    All Mesh/StatMesh attributes are also present via a super() call to
    StatMesh.__init__().
    """

    def __init__(self):
        """Create an empty MeshTally object and set default values."""
        if not HAVE_PYMOAB:
            raise NotImplementedError(
                "PyMOAB is not available, " "unable to create Meshtally Mesh."
            )

        self.tally_number = None
        self.particle = "neutron"
        self.set_default_tag_names()

    @property
    def tag_names(self):
        return self._tag_names

    @tag_names.setter
    def tag_names(self, tag_names):
        _check_meshtally_tag_names(tag_names)
        self._tag_names = tag_names

    def set_default_tag_names(self):
        """Set default tag_names according to particle type."""
        self.tag_names = (
            "{0}_result".format(self.particle),
            "{0}_result_rel_error".format(self.particle),
            "{0}_result_total".format(self.particle),
            "{0}_result_total_rel_error".format(self.particle),
        )

    def tag_flux_error_from_tally_results(self, result, rel_err, res_tot, rel_err_tot):
        """
        This function uses the output tally result, rel_err, res_tot and the
        rel_err_tot to set the flux and error tags.

        Parameters
        ----------
        result : numpy array
            This numpy array contains the flux data read from MCNP meshtal or
            OpenMC state point file.
            The shape of this numpy array is (ves, num_e_groups).
        rel_err : numpy array
            This numpy array contains the relative error data read from MCNP
            meshtal or OpenMC state point file. The shape of this numpy array
            is (num_ves, num_e_groups).
        res_tot : list
            The total results.
        rel_err_tot : list
            Relative error of total results.
        """

        num_ves = len(self)
        self.tag(
            name=self.tag_names[0],
            value=result,
            doc="{0} flux".format(self.particle),
            tagtype=NativeMeshTag,
            size=self.num_e_groups,
            dtype=float,
        )
        # set result_rel_error tag
        self.tag(
            name=self.tag_names[1],
            value=rel_err,
            doc="{0} flux relative error".format(self.particle),
            tagtype=NativeMeshTag,
            size=self.num_e_groups,
            dtype=float,
        )
        # set result_total tag
        self.tag(
            name=self.tag_names[2],
            value=res_tot,
            doc="total {0} flux".format(self.particle),
            tagtype=NativeMeshTag,
            size=1,
            dtype=float,
        )
        # set result_total_rel_error tag
        self.tag(
            name=self.tag_names[3],
            value=rel_err_tot,
            doc="total {0} flux relative error".format(self.particle),
            tagtype=NativeMeshTag,
            size=1,
            dtype=float,
        )


######################################################
# private helper functions for structured mesh methods
######################################################


def _structured_find_idx(dims, ijk):
    """Helper method fo structured_get_vertex and structured_get_hex.

    For tuple (i,j,k), return the number N in the appropriate iterator.
    """

    dim0 = [0] * 3
    for i in range(0, 3):
        if dims[i] > ijk[i] or dims[i + 3] <= ijk[i]:
            raise MeshError(str(ijk) + " is out of bounds")
        dim0[i] = ijk[i] - dims[i]
    i0, j0, k0 = dim0
    n = (
        ((dims[4] - dims[1]) * (dims[3] - dims[0]) * k0)
        + ((dims[3] - dims[0]) * j0)
        + i0
    )
    return n


def _structured_step_iter(it, n):
    """Helper method for structured_get_vertex and structured_get_hex

    Return the nth item in the iterator.
    """

    it.step(n)
    r = next(it)
    it.reset()
    return r


def _structured_iter_setup(dims, order, **kw):
    """Setup helper function for StrMesh iterator functions

    Given dims and the arguments to the iterator function, return
    a list of three lists, each being a set of desired coordinates,
    with fastest-changing coordinate in the last column), and the
    ordmap used by _structured_iter to reorder each coodinate to (i,j,k).
    """

    # a valid order has the letters "x", "y", and "z"
    # in any order without duplicates
    if not (
        len(order) <= 3
        and len(set(order)) == len(order)
        and all([a in "xyz" for a in order])
    ):
        raise MeshError("Invalid iteration order: " + str(order))

    # process kw for validity
    spec = {}
    for idx, d in enumerate("xyz"):
        if d in kw:
            spec[d] = kw[d]
            if not isinstance(spec[d], Iterable):
                spec[d] = [spec[d]]
            if not all(x in range(dims[idx], dims[idx + 3]) for x in spec[d]):
                raise MeshError("Invalid iterator kwarg: " "{0}={1}".format(d, spec[d]))
            if d not in order and len(spec[d]) > 1:
                raise MeshError(
                    "Cannot iterate over"
                    + str(spec[d])
                    + "without a proper iteration order"
                )
        if d not in order:
            order = d + order
            spec[d] = spec.get(d, [dims[idx]])

    # get indices and ordmap
    indices = []
    for L in order:
        idx = "xyz".find(L)
        indices.append(spec.get(L, range(dims[idx], dims[idx + 3])))

    ordmap = ["zyx".find(L) for L in order]
    return indices, ordmap


def _structured_iter(indices, ordmap, dims, it):
    """Iterate over the indices lists, yielding _structured_step_iter(it) for
    each.
    """

    d = [0, 0, 1]
    d[1] = dims[3] - dims[0]
    d[0] = (dims[4] - dims[1]) * d[1]
    mins = [dims[2], dims[1], dims[0]]
    offsets = (
        [(a - mins[ordmap[x]]) * d[ordmap[x]] for a in indices[x]] for x in range(3)
    )
    for ioff, joff, koff in itertools.product(*offsets):
        yield _structured_step_iter(it, (ioff + joff + koff))


if HAVE_PYMOAB:

    def mesh_iterate(mesh, mesh_type=3, topo_type=types.MBMAXTYPE):
        return meshset_iterate(mesh, 0, topo_type, mesh_type, recursive=True)

    def meshset_iterate(
        pymb,
        meshset=0,
        entity_type=types.MBMAXTYPE,
        dim=-1,
        arr_size=1,
        recursive=False,
    ):
        return MeshSetIterator(pymb, meshset, entity_type, dim, arr_size, recursive)


class MeshSetIterator(object):
    def __init__(self, inst, meshset, entity_type, dim=-1, arr_size=1, recursive=False):
        self.pymb = inst
        self.meshset = meshset
        self.ent_type = entity_type
        self.dimension = dim
        self.arr_size = arr_size
        self.recur = recursive
        self.reset()

    def reset(self):
        # if a specific dimension is requested, return only that dimension
        if self.ent_type != types.MBMAXTYPE:
            ents = self.pymb.get_entities_by_type(
                self.meshset, self.ent_type, self.recur
            )
        # if a specific type is requested, return only that type
        elif self.dimension != -1:
            ents = self.pymb.get_entities_by_dimension(
                self.meshset, self.dimension, self.recur
            )
        # otherwise return everything
        else:
            ents = self.pymb.get_entities_by_handle(self.meshset, self.recur)

        self.pos = 0
        self.size = len(ents)
        self.entities = ents

    def __iter__(self):
        return self

    def __next__(self):
        if self.pos >= self.size:
            raise StopIteration
        else:
            self.pos, value = self.pos + 1, self.entities[self.pos]
            return value

    # for Python2 compatability
    def next(self):
        return self.__next__()

    def step(self, num_steps):
        self.pos += int(num_steps)  # protecting from Python 3 auto-promotion
        at_end = False
        if self.pos >= self.size:
            self.pos = self.size - 1
            at_end = True
        return at_end


def _cell_fracs_sort_vol_frac_reverse(cell_fracs):
    """
    Sort cell_fracs according to the order of increasing idx and decreasing
    with vol_frac.

    Parameters
    ----------
    cell_fracs : structured array
            The output from dagmc.discretize_geom(). A sorted, one dimensional
            array, each entry containing the following fields:

                :idx: int
                    The volume element index.
                :cell: int
                    The geometry cell number.
                :vol_frac: float
                    The volume fraction of the cell withing the mesh ve.
                :rel_error: float
                    The relative error associated with the volume fraction.

            The array must be sorted with respect to both idx and cell, with
            cell changing fastest.

    Returns
    -------
    cell_fracs : structured array
        Sorted cell_fracs.
    """

    # sort ascending along idx and vol_frac
    # ndarray.sort can't sort using desending sequence.
    # Multiply the vol_frac to -1.0 to sort the vol_frac in reverse order.
    cell_fracs["vol_frac"] *= -1.0
    cell_fracs.sort(order=["idx", "vol_frac"])
    cell_fracs["vol_frac"] *= -1.0
    return cell_fracs


def _check_meshtally_tag_names(tag_names):
    """Make sure tag_names is an iterable of 4 strings."""
    # check iterable
    if not check_iterable(tag_names):
        raise ValueError("The given tag_names is not an Iterable.")

    # check length of 4
    if len(tag_names) != 4:
        raise ValueError("The length of tag_names is not 4.")

    # check content strings
    for item in tag_names:
        if not isinstance(item, str):
            raise ValueError("The content of tag_names " "should be strings")

    # tag_names should be a string with length of 4
    if isinstance(tag_names, str):
        raise ValueError("The tag_names should not be a single string")

    return True
