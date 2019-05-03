"""Python wrapper for material library."""
from __future__ import division, unicode_literals

# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from libc.stdlib cimport malloc, free
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from libcpp cimport bool as cpp_bool
# from cython.operator cimport reference as ref
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Python imports
import collections
cimport numpy as np
import numpy as np
from warnings import warn
from pyne.utils import QAWarning
import os
import sys
if sys.version_info[0] >= 3:
    # Python2 basestring is now Python3 string
    basestring = str

import tables as tb

# local imports
from pyne cimport cpp_material
from pyne cimport cpp_material_library
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

cimport cpp_jsoncpp
cimport pyne.jsoncpp as jsoncpp
import pyne.jsoncpp as jsoncpp

cimport pyne.nucname as nucname
import pyne.nucname as nucname

cimport pyne.data as data
import pyne.data as data

cimport pyne.material as material
import pyne.material as material


warn(__name__ + " is not yet QA compliant.", QAWarning)

# Maximum 32-bit signed int
DEF INT_MAX = 2147483647


#
#  Material Library
#

cdef class _MaterialLibrary:
    """This class allows the definitions of a set of material, stored by names
    """


    def __cinit__(self, lib=None, datapath="/materials", nucpath="/nucid"):
        """Parameters
        ----------
        lib : dict-like, str, or None, optional
            The data to intialize the material library with.  If this is a
            string, it is interpreted as a path to a file.
        datapath : str, optional
            The path in the heirarchy to the data table in an HDF5 file.
        nucpath : str, optional
            The path in the heirarchy to the nuclide array in an HDF5 file.

        """
        if lib != None:
            if sys.version_info[0] >= 3 and isinstance(lib, bytes):
                lib = lib.decode()
            if isinstance(lib, collections.Mapping):
                self._inst = new cpp_material_library.MaterialLibrary()
                list_ = []
                for key in sorted(lib.keys()):
                    mat = lib[key]
                    self.__setitem__(key, material.ensure_material(mat))
            elif isinstance(lib, unicode):
                c_filename = lib.encode('utf-8')
                c_datapath = datapath.encode('utf-8')
                self._inst = new cpp_material_library.MaterialLibrary(c_filename, c_datapath)
            elif isinstance(lib, basestring):
                c_filename = lib
                c_datapath = datapath
                self._inst = new cpp_material_library.MaterialLibrary(c_filename, c_datapath)
            elif isinstance(lib, collections.Sequence):
                self._inst = new cpp_material_library.MaterialLibrary()
                for key in sorted(lib.keys()):
                    mat = lib[key]
                    self.__setitem__(key, material.ensure_material(mat))
            else:
                msg = "Could not initialize library with lib type {0!r}"
                raise TypeError(msg.format(type(lib)))
        else:
            self._inst = new cpp_material_library.MaterialLibrary()

    def __dealloc__(self):
        """MaterialLibrary C++ destructor."""
        del self._inst

    def from_hdf5(self, filename, datapath="/materials", nucpath="", protocol=1):
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        cdef char * c_datapath
        datapath_bytes = datapath.encode('UTF-8')
        c_datapath = datapath_bytes
        cdef char * c_nucpath
        nucpath_bytes = nucpath.encode('UTF-8')
        c_nucpath = nucpath_bytes
        self._inst.from_hdf5(c_filename, c_datapath, c_nucpath, protocol)


    def write_hdf5(self, filename, datapath="/materials", nucpath="/nucid"):
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        cdef char * c_datapath
        datapath_bytes = datapath.encode('UTF-8')
        c_datapath = datapath_bytes
        cdef char * c_nucpath
        nucpath_bytes = nucpath.encode('UTF-8')
        c_nucpath = nucpath_bytes
        self._inst.write_hdf5(c_filename, c_datapath, c_nucpath)

    def add_material(self, mat):
        cdef std_string c_matname
        if isinstance(mat, material.Material):
            value_proxy = material.Material(
                mat, free_mat=not isinstance(mat, material._Material))
            self._inst.add_material( (<material._Material> value_proxy).mat_pointer[0])
        else:
            raise TypeError("the material must be a material or a stri but is a "
                            "{0}".format(type(mat)))

    def del_material(self, mat):
        cdef std_string c_matname
        if isinstance(mat, material._Material):
            c_matname = std_string(< char * > mat).mat_pointer.metadata["name"] 
        elif isinstance(mat, basestring):
            c_matname = std_string(< char * > mat)
        else:
            raise TypeError("the material must be a material or a stri but is a "
                            "{0}".format(type(mat)))
        self._inst.del_material(c_matname)

    def get_material(self, matname):
        # Get the correct cpp_material
        cdef cpp_material.Material c_mat
        cdef std_string c_matname
        c_matname = matname
        c_mat = self._inst.get_material(c_matname);
        # build a PyNE Material object form the cpp_material
        py_mat = material.Material(free_mat = False)
        (< material._Material > py_mat).mat_pointer = new cpp_material.Material(c_mat.comp, c_mat.mass, c_mat.density, c_mat.atoms_per_molecule, c_mat.metadata)        
        return py_mat

    def merge(self, mat_library):
        if isinstance(mat_library, _MaterialLibrary):
            self._inst.merge(mat_library._inst)
        else:
            raise TypeError("the material library must be a MaterialLibrary but is a "
                            "{0}".format(type(mat_library)))

    cdef cpp_set[std_string] get_matlist(self):
        return self._inst.get_matlist()

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
        self._inst.load_json(deref((<jsoncpp.Value> json)._inst))

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

    def __setitem__(self, key, value):
        if isinstance(key, int):
            key = str(key)

        value.metadata["name"] = key.encode('utf-8')
        value_proxy = material.Material(value, free_mat=not isinstance(value, material._Material))
        self._inst.add_material( (<material._Material> value_proxy).mat_pointer[0])

    def __getitem__(self, key):
        if isinstance(key, basestring):
            key = key.encode('UTF-8')
        elif isinstance(key, int):
            key = str(key).encode('UTF-8')
        py_mat = self.get_material(key)
        return py_mat
    
    def __len__(self):
        return self.mat_library.size()

    def __delitem__(self, key):
        if isinstance(key, basestring):
            key = key.encode('UTF-8')
        elif isinstance(key, int):
            key = str(key).encode('UTF-8')
        self.del_material(key)

    def __iter__(self):
        mat_lib_dict = map_to_dict_str_matp(self._inst.get_mat_library())
        self._iter_mat_lib = mat_lib_dict
        mat_lib_iter = iter(mat_lib_dict)
        return mat_lib_iter


class MaterialLibrary(_MaterialLibrary, collections.MutableMapping):
    """Material composed of nuclides.

    Parameters
    ----------

    """
    def __str__(self):
        return " "

    def __repr__(self):
        return ""

        
        
        
# <string, Material *>

cdef cpp_map[std_string, matp] dict_to_map_str_matp(dict pydict):
    cdef material._Material pymat
    cdef cpp_material.Material * cpp_matp
    cdef cpp_map[std_string, matp] cppmap = cpp_map[std_string, matp]()
    cdef cpp_pair[std_string, matp] item

    for key, value in pydict.items():
        pymat = value
        cpp_matp = pymat.mat_pointer
        #cppmap[std_string(key)] = cpp_matp
        item = cpp_pair[std_string, matp](std_string(< char * > key), cpp_matp)
        cppmap.insert(item)

    return cppmap


cdef dict map_to_dict_str_matp(cpp_map[std_string, matp] cppmap):
    pydict = {}
    cdef material._Material pymat
    cdef cpp_map[std_string, matp].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pymat = material.Material()
        pymat.mat_pointer[0] = deref(deref(mapiter).second)
        pydict[< char * > deref(mapiter).first.c_str()] = pymat
        inc(mapiter)

    return pydict

