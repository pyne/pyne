"""Python wrapper for material library."""
from __future__ import division, unicode_literals

# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from libc.stdlib cimport malloc, free
from libcpp.unordered_map cimport unordered_map as cpp_umap
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

_INTEGRAL_TYPES = (int, np.integer, np.bool_)
 

#
#  Material Library
#

cdef class _MaterialLibrary:
    """This class allows the definitions of a set of materials, stored by name
    """


    def __cinit__(self, lib=None, datapath="/materials", nucpath="/nucid"):
        """
        Parameters
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
            elif isinstance(lib, basestring) or isinstance(lib, unicode): 
                # Python2: basestring = (std + unicode)
                c_filename = lib.encode('UTF-8')
                c_datapath = datapath.encode('UTF-8')
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
        """Loads data from an HDF5 file into this material library.
        Parameters
        ----------
        file : str
            A path to an HDF5 file.
        datapath : str, optional
            The path in the heirarchy to the data table in an HDF5 file.
        nucpath : str, optional
            The path in the heirarchy to the nuclide array in an HDF5 file.
        protocol : int, optional
            Specifies the protocol to use to read in the data.  Different
            protocols are used to represent different internal structures in
            the HDF5 file.

        """ 
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
        """Writes this material library to an HDF5 file.
        Parameters
        ----------
        filename : str
            A path to an HDF5 file.
        datapath : str, optional
            The path in the heirarchy to the data table in an HDF5 file.
        nucpath : str, optional
            The path in the heirarchy to the nuclide array in an HDF5 file.
        """
        
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

    def add_material(self, key, mat):
        """Add a Material to this material library.
        Parameters
        ----------
        key : str or int (converted to str)
            key to regisgter the material, if material has no name attribute
            will be added as the name metadata of the material
        mat : Material
            PyNE material object be added to this material library
        """
        
        cdef std_string c_matname
        if isinstance(key, basestring):
            key = key.encode('UTF-8')
        elif isinstance(key, _INTEGRAL_TYPES):
            key = str(key).encode('UTF-8')
        if isinstance(mat, material.Material):
            value_proxy = material.Material(
                mat, free_mat=not isinstance(mat, material._Material))
            self._inst.add_material( key, (<material._Material> value_proxy).mat_pointer[0])
        else:
            raise TypeError("Material must be a PyNE Material or a string but is a "
                            "{0}".format(type(mat)))

    def remove_material(self, mat):
        """Remove a Material from this material library.
        Parameters
        ----------
        mat : Material or str
            PyNE material object or material name be removed from this material library
            if a material is provided it needs a name metadata as material are
            removed according to their names
        """

        cdef std_string c_matname
        if isinstance(mat, basestring):
            c_matname = std_string(< char * > mat)
        else:
            raise TypeError("the argument must be a string (material name) but is a "
                            "{0}".format(type(mat)))
        self._inst.del_material(c_matname)

    def get_material(self, key):
        """Get a Material from this material library.
        Parameters
        ----------
        key : str or int (converted to str)
            key of the material to return 
        """
        # Get the correct cpp_material
        cdef cpp_material.Material c_mat
        cdef std_string c_matname
        if isinstance(key, _INTEGRAL_TYPES):
            c_matname = str(key).encode('UTF-8')
        else:
            c_matname = key
        c_mat = self._inst.get_material(c_matname)

        # build a PyNE Material object form the cpp_material
        py_mat = material.Material(free_mat = False)
        (< material._Material > py_mat).mat_pointer = new cpp_material.Material(c_mat.comp, c_mat.mass, c_mat.density, c_mat.atoms_per_molecule, c_mat.metadata)        
        return py_mat

    def merge(self, mat_library):
        """Merge a material library into this material library.
        Parameters
        ----------
        mat_library : MaterialLibrary
            Material Library to merge 
        """

        if isinstance(mat_library, _MaterialLibrary):
            self._inst.merge(mat_library._inst)
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
        cdef cpp_pair[std_string, cpp_material.Material] item
        if isinstance(key, _INTEGRAL_TYPES):
            key = str(key).encode('UTF-8')
        else:
            key = key.encode('UTF-8') 
        value_proxy = material.Material(value, free_mat=not isinstance(value, material._Material))
        self._inst.add_material(key, deref((<material._Material> value_proxy).mat_pointer))


    def __getitem__(self, key):
        cdef cpp_material.Material* c_mat
        cdef std_string c_matname
        if isinstance(key, basestring):
            key = key.encode('UTF-8')
        elif isinstance(key, _INTEGRAL_TYPES):
            key = str(key).encode('UTF-8')
        
        c_matname = key
        c_mat = self._inst.get_material_ptr(<std_string>c_matname)

        # build a PyNE Material object form the cpp_material
        py_mat = material.Material(free_mat = False)
        (< material._Material > py_mat).mat_pointer = c_mat
        return py_mat
    
    def __len__(self):
        if hasattr(self, 'mat_library'):
            return self.mat_library.size()
        else:
            return 0

    def __delitem__(self, key):
        if isinstance(key, basestring):
            key = key.encode('UTF-8')
        elif isinstance(key, _INTEGRAL_TYPES):
            key = str(key).encode('UTF-8')
        self.del_material(key)

    def __iter__(self):
        mat_lib_dict = map_to_dict_str_matp(self._inst.get_mat_library())
        self._iter_mat_lib = mat_lib_dict
        mat_lib_iter = iter(mat_lib_dict)
        return mat_lib_iter


class MaterialLibrary(_MaterialLibrary, collections.MutableMapping):
    """The material library is a collection of unique keys mapped to
    Material objects.

    """
    def __repr__(self):
        libs = ["{0!r}={1!r}".format(k, m) for k, m in self.items()]
        libs = "{" + ", ".join(libs) + "}"
        return "pyne.material.MaterialLibrary({0})".format(libs)
        
        
# Python dict to u_map<string, Material *> 
cdef cpp_umap[std_string, cpp_material.Material*] dict_to_map_str_matp(dict pydict):
    cdef cpp_material.Material * cpp_matp
    cdef cpp_umap[std_string, matp ] cppmap = cpp_umap[std_string, matp ]()
    cdef cpp_pair[std_string, cpp_material.Material* ] item

    for key, value in pydict.items():
        py_mat = material.Material(free_mat = False)
        py_mat = value
        cpp_matp = (< material._Material > py_mat).mat_pointer
        #cppmap[std_string(key)] = cpp_matp
        item = cpp_pair[std_string, matp](std_string(< char * > key), cpp_matp)
        cppmap.insert(item)

    return cppmap

# u_map<string, Material *> to python dict
cdef dict map_to_dict_str_matp(cpp_umap[std_string, cpp_material.Material*] cppmap):
    pydict = {}
    cdef cpp_umap[std_string, cpp_material.Material*].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        py_mat = material.Material(free_mat = False)
        (< material._Material > py_mat).mat_pointer = deref(mapiter).second
        pydict[< char * > deref(mapiter).first.c_str()] = py_mat
        inc(mapiter)

    return pydict
