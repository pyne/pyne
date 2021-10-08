"""Python wrapper for material library."""
from __future__ import division, unicode_literals
import pyne.material as material
import pyne.data as data
import pyne.nucname as nucname
import pyne.jsoncpp as jsoncpp
import pyne.stlcontainers as conv
from pyne cimport cpp_material_library
from pyne cimport cpp_material
import tables as tb
import sys
import os
from pyne.utils import QA_warn
import numpy as np

# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from libc.stdlib cimport malloc, free
from libcpp.unordered_map cimport unordered_map as cpp_umap
from libcpp.memory cimport shared_ptr

from libcpp.vector cimport vector as cpp_vector
from libcpp cimport bool as cpp_bool
# from cython.operator cimport reference as ref
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Python imports
import collections
try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
cimport numpy as np
if sys.version_info[0] >= 3:
    # Python2 basestring is now Python3 string
    basestring = str


# local imports
cimport pyne.stlcontainers as conv

cimport cpp_jsoncpp
cimport pyne.jsoncpp as jsoncpp

cimport pyne.nucname as nucname

cimport pyne.data as data

cimport pyne.material as material


QA_warn(__name__)

_INTEGRAL_TYPES = (int, np.integer, np.bool_)


#
#  Material Library
#

cdef class _MaterialLibrary:
    """This class allows the definitions of a set of materials, stored by name
    """

    def __cinit__(self, lib=None, datapath="/materials"):
        """
        Parameters
        ----------
        lib : dict-like, str, or None, optional
            The data to initialize the material library with.  If this is a
            string, it is interpreted as a path to a file (hdf5 or json format supported).
        datapath : str, optional
            The path in the hierarchy to the data table in an HDF5 file.
            The path in the hierarchy to the nuclide array in an HDF5 file.
        """
        if lib != None:
            if sys.version_info[0] >= 3 and isinstance(lib, bytes):
                lib = lib.decode()
            if (isinstance(lib, basestring) or isinstance(lib, unicode)) \
                    and not isinstance(lib, collectionsAbc.Mapping):
                # Python2: basestring = (std + unicode)
                c_filename = lib.encode('UTF-8')
                c_datapath = datapath.encode('UTF-8')
                self._inst = new cpp_material_library.MaterialLibrary(c_filename, c_datapath)
            elif isinstance(lib, collectionsAbc.Mapping) or isinstance(lib, collectionsAbc.Sequence):
                self._inst = new cpp_material_library.MaterialLibrary()
                for key in sorted(lib.keys()):
                    mat = lib[key]
                    self.__setitem__(key, material.ensure_material(mat))
            else:
                msg = "Could not initialize library with lib type {0!r}"
                raise TypeError(msg.format(type(lib)))
        else:
            self._inst = new cpp_material_library.MaterialLibrary()

#   def __dealloc__(self):
#        """MaterialLibrary C++ destructor."""
#        del self._inst

    def from_hdf5(self, filename, datapath):
        """Loads data from an HDF5 file into this material library.
        Parameters
        ----------
        file : str
            A path to an HDF5 file.
        datapath : str, optional
            The path in the hierarchy to the data table in an HDF5 file.
            The path in the hierarchy to the nuclide array in an HDF5 file.
        protocol : int, optional
            Specifies the protocol to use to read in the data.  Different
            protocols are used to represent different internal structures in
            the HDF5 file.

        """
        c_filename = filename.encode('UTF-8')
        c_datapath = datapath.encode('UTF-8')
        self._inst.from_hdf5(c_filename, c_datapath)

    def write_hdf5(self, filename, datapath="/materials", h5_overwrite=False):
        """Writes this material library to an HDF5 file.
        Parameters
        ----------
        filename : str
            A path to an HDF5 file.
        datapath : str, optional
            The path in the hierarchy to the data table in an HDF5 file.
        """
        c_filename = filename.encode('UTF-8')
        c_datapath = datapath.encode('UTF-8')
        self._inst.write_hdf5(c_filename, c_datapath, h5_overwrite)

    def remove_material(self, mat_name):
        """Remove a Material from this material library.
        Parameters
        ----------
        mat_name : str Name of the material be removed from this material library
        """
        if isinstance(mat_name, basestring):
            c_matname = mat_name.encode('UTF-8')
            self._inst.del_material(< std_string > c_matname)
        else:
            raise TypeError("the argument must be a string (material name) but is a "
                            "{0}".format(type(mat_name)))

    def merge(self, mat_library):
        """Merge a material library into this material library.
        Parameters
        ----------
        mat_library : MaterialLibrary
            Material Library to merge
        """

        if isinstance(mat_library, _MaterialLibrary):
            self._inst.merge(< cpp_material_library.MaterialLibrary*>( < _MaterialLibrary > mat_library)._inst)
        else:
            raise TypeError("the material library must be a MaterialLibrary but is a "
                            "{0}".format(type(mat_library)))

    cdef cpp_set[std_string] get_keylist(self):
        return self._inst.get_keylist()

    cdef cpp_set[int] get_nuclist(self):
        return self._inst.get_nuclist()

    def load_json(self, json):
        """load_json(json)
        Loads a JSON instance into this Material.

        Parameters
        ----------
        json : jsoncpp.Value
            An object-type JSON value.

        """
        self._inst.load_json(deref(( < jsoncpp.Value > json)._inst))

    def dump_json(self):
        """dump_json()
        Dumps the material to a JSON object.

        Returns
        -------
        val : jsoncpp.Value
            An object-type JSON value.

        """
        cdef jsoncpp.Value val = jsoncpp.Value(view=False)
        val._inst[0] = self._inst.dump_json()
        return val

    def from_json(self, filename):
        """from_json(char * filename)
        Initialize a Material object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to text file that contains the data to read in.

        """
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        self._inst.from_json(c_filename)

    def write_json(self, filename):
        """write_json(filename)
        Writes the material to a JSON file.

        Parameters
        ----------
        filename : str
            Path to text file to write the data to.  If the file already
            exists, it will be overwritten.

        """
        filename = filename.encode()
        self._inst.write_json(filename)

    def write_openmc(self, filename):
        """write_openmc(filename)
        Writes the material library in an OpenMC XML format.

        Parameters
        ----------
        filename : str
            Path to text file to write the data to.  If the file already
            exists, it will be overwritten.

        """
        filename = filename.encode()
        self._inst.write_openmc(filename)

    def __setitem__(self, key, value):
        """Add a Material to this material library, if the material key already exist it will be overwritten.
         Parameters
         ----------
         key : str or int (converted to str)
             key to register the material, if material has no name attribute
             will be added as the name metadata of the material
         mat : Material
             PyNE material object be added to this material library
         """
        cdef cpp_pair[std_string, cpp_material.Material] item
        value_proxy = material.Material(
            value, free_mat=not isinstance(value, material._Material))
        self._inst.add_material(ensure_material_key(key), deref(( < material._Material > value_proxy).mat_pointer))

    def __getitem__(self, key):
        cdef shared_ptr[cpp_material.Material] c_mat
        c_mat = self._inst.get_material_ptr(< std_string > ensure_material_key(key))

        # build a PyNE Material object from the cpp_material
        py_mat = material.Material(free_mat=False)
        ( < material._Material > py_mat).mat_pointer = c_mat.get()
        return py_mat

    def __len__(self):
        return self._inst.size()

    def __delitem__(self, key):
        c_key = ensure_material_key(key)
        self._inst.del_material(< std_string > c_key)

    def __iter__(self):
        mat_lib_dict = matlib_to_dict_str_matp(deref(self._inst))
        self._iter_mat_lib = mat_lib_dict
        mat_lib_iter = iter(mat_lib_dict)
        return mat_lib_iter


class MaterialLibrary(_MaterialLibrary, collectionsAbc.MutableMapping):
    """The material library is a collection of unique keys mapped to
    Material objects.
    """

    def __repr__(self):
        libs = ["{0!r}={1!r}".format(k, m) for k, m in self.items()]
        libs = "{" + ", ".join(libs) + "}"
        return "pyne.material.MaterialLibrary({0})".format(libs)


def ensure_material_key(key):
    if isinstance(key, basestring):
        key = key.encode('UTF-8')
    elif isinstance(key, _INTEGRAL_TYPES):
        key = str(key).encode('UTF-8')
    return key


# Python dict to u_map<string, Material *>
cdef cpp_umap[std_string, shared_ptr[cpp_material.Material]] dict_to_map_str_matp(dict pydict):
    cdef shared_ptr[cpp_material.Material] cpp_matp
    cdef cpp_umap[std_string, shared_ptr[cpp_material.Material]] cppmap = cpp_umap[std_string, shared_ptr[cpp_material.Material]]()
    cdef cpp_pair[std_string, shared_ptr[cpp_material.Material]] item

    for key, value in pydict.items():
        py_mat = material.Material(free_mat=False)
        py_mat = value
        cpp_matp = shared_ptr[cpp_material.Material](( < material._Material > py_mat).mat_pointer)
        #cppmap[std_string(key)] = cpp_matp
        item = cpp_pair[std_string, shared_ptr[cpp_material.Material]](std_string( < char * > key),         cpp_matp)
        cppmap.insert(item)

    return cppmap

# MaterialLibrary to python dict
cdef dict matlib_to_dict_str_matp(cpp_material_library.MaterialLibrary cpp_mat_lib):
    pydict = {}
    cdef cpp_umap[std_string, shared_ptr[cpp_material.Material]].iterator mapiter = cpp_mat_lib.begin()

    while mapiter != cpp_mat_lib.end():
        py_mat = material.Material(free_mat=False)
        ( < material._Material > py_mat).mat_pointer = (deref(mapiter).second).get()
        pydict[ < char * > deref(mapiter).first.c_str()] = py_mat
        inc(mapiter)

    return pydict
