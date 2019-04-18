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

cdef class MaterialLibrary:
    """This class allows the definitions of a set of material, stored by names
    """

    def __cinit(self, *args, **kwargs):
        """MaterialLibrary C++ default constructor."""
        self._inst = new cpp_material_library.MaterialLibrary()
    def __init__(self, lib=None, datapath="/materials", nucpath="/nucid"):
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
        if sys.version_info[0] >=3 and isinstance(lib, bytes):
            lib = lib.decode()
        if isinstance(lib, collections.Mapping):
            self._inst = new cpp_material_library.MaterialLibrary()
            for key, mat in lib.items():
                print("here?", flush=True)
                self.__setitem__(key, material.ensure_material(mat))
        elif isinstance(lib, basestring):
            c_filename = std_string( < char * > data)
            c_datapath = std_string( < char * > datapath)
            self._inst = new cpp_material_library.MaterialLibrary(c_filename, c_datapath)
        elif isinstance(lib, collections.Sequence):
            self._inst = new cpp_material_library.MaterialLibrary()
            for key, mat in lib:
                self.__setitem__(key, material.ensure_material(mat))
        else:
            msg = "Could not initialize library with lib type {0!r}"
            raise TypeError(msg.format(type(lib)))

    def __dealloc__(self):
        """MaterialLibrary C++ destructor."""
        del self._inst

    def from_hdf5(self, filename, datapath="/materials", protocol=1):
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        cdef char * c_datapath
        datapath_bytes = datapath.encode('UTF-8')
        c_datapath = datapath_bytes

        self._inst.from_hdf5(c_filename, c_datapath, protocol)

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
        cdef material._Material value_proxy
        if isinstance(mat, material._Material):
            value_proxy = material.Material(mat, free_mat=not isinstance(mat, material._Material))
            (<cpp_material_library.MaterialLibrary *> self._inst).add_material(
                    value_proxy.mat_pointer[0])
        else:
            raise TypeError("the material must be a material or a stri but is a "
                            "{0}".format(type(mat)))

    def del_material(self, mat):
        cdef std_string c_matname
        if isinstance(mat, material._Material):
            c_matname = std_string( < char * > mat).mat_pointer.metadata["name"] 
        elif isinstance(mat, basestring):
            c_matname = std_string( < char * > mat)
        else:
            raise TypeError("the material must be a material or a stri but is a "
                            "{0}".format(type(mat)))
        self._inst.del_material(c_matname)

    def get_material(self, matname):
        # Get the correct cpp_material
        cdef cpp_material.Material c_mat
        cdef std_string c_matname
        

        c_matname = std_string( < char * > matname)
        c_mat = self._inst.get_material(c_matname)

        # build a PyNE Material object form the cpp_material
        cdef jsoncpp.Value metadata = jsoncpp.Value()
        metadata.__set_instance__(c_mat.metadata)
        #del metadata._inst
        #metadata._inst = &c_mat.metadata
        
        py_mat = material._Material(    
            c_mat.comp,
            c_mat.mass,
            c_mat.density,
            c_mat.atoms_per_molecule,
            metadata)
        return py_mat

    def merge(self, mat_library):
        if isinstance(mat_library, MaterialLibrary):
            self._inst.merge(mat_library._inst)
        else:
            raise TypeError("the material library must be a MaterialLibrary but is a "
                            "{0}".format(type(mat_library)))

    cdef cpp_set[std_string] get_matlist(self):
        return self._inst.get_matlist()

    cdef cpp_set[int] get_nuclist(self):
        return self._inst.get_nuclist()

    def __setitem__(self, key, value):
        value.metadata["name"] = key.encode('utf-8')
        print("here1", flush=True) 
        cdef material._Material value_proxy
        print("here2", flush=True) 
        value_proxy = material.Material(value, free_mat=not isinstance(value, material._Material))
        print("here3", flush=True) 
        (<cpp_material_library.MaterialLibrary *> self._inst).add_material(
                value_proxy.mat_pointer[0])

    def __getitem__(self, key):
        return self.get_material(key.encode('utf-8'))

    def __delitem__(self, key):
        self.del_material(key.encode('utf-8'))


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
        item = cpp_pair[std_string, matp](std_string( < char * > key), cpp_matp)
        cppmap.insert(item)

    return cppmap


cdef dict map_to_dict_str_matp(cpp_map[std_string, matp] cppmap):
    pydict = {}
    cdef material._Material pymat
    cdef cpp_map[std_string, matp].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pymat = material.Material()
        pymat.mat_pointer[0] = deref(deref(mapiter).second)
        pydict[ < char * > deref(mapiter).first.c_str()] = pymat
        inc(mapiter)

    return pydict

