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

    def __cinit(self, filename, datapath="/materials"):
        """MaterialLibrary C++ constructor."""
        cdef std_string c_filename
        cdef std_string c_datapath

        c_filename = std_string(< char * > filename)
        c_datapath = std_string(< char * > datapath)

        self._inst = new cpp_material_library.MaterialLibrary(c_filename, c_datapath)

    def __dealloc__(self):
        """MaterialLibrary C++ destructor."""
        del self._inst

    def from_hdf5(self, filename, datapath="/materials", protocol=1):
        cdef std_string c_filename
        cdef std_string c_datapath

        c_filename = std_string(< char * > filename)
        c_datapath = std_string(< char * > datapath)

        self._inst.from_hdf5(c_filename, c_datapath, protocol)

    def write_hdf5(self, filename, datapath="/materials", nucpath="/nucid",
            chunksize=100):
        cdef char * c_filename
        filename_bytes = filename.encode('UTF-8')
        c_filename = filename_bytes
        cdef char * c_datapath
        datapath_bytes = datapath.encode('UTF-8')
        c_datapath = datapath_bytes
        cdef char * c_nucpath
        nucpath_bytes = nucpath.encode('UTF-8')
        c_nucpath = nucpath_bytes
        self._inst.write_hdf5(c_filename, c_datapath, c_nucpath, chunksize)

    def add_material(self, mat):
        cdef std_string c_matname
        if isinstance(mat, material._Material):
            self._inst.add_material( < cpp_material.Material > ( < material._Material > mat).mat_pointer[0])
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

    def get_material(self, mat):
        # Get the correct cpp_material
        cdef cpp_material.Material c_mat
        cdef std_string c_matname
        cdef jsoncpp.Value metadata

        c_matname = std_string(< char * > mat)
        c_mat = self._inst.get_material(c_matname)

        # build a PyNE Material object form the cpp_material
        metadata = jsoncpp.Value(mat.metadata)
        py_mat = material._Material(
            c_mat.comp,
            c_mat.mass,
            c_mat.density,
            c_mat.atoms_per_molecule,
            metadata)
        return py_mat

    cdef cpp_set[std_string] get_matlist(self):
        return self._inst.get_matlist()

    cdef cpp_set[int] get_nuclist(self):
        return self._inst.get_nuclist()


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
        item = cpp_pair[std_string, matp](std_string(<char *> key), cpp_matp)
        cppmap.insert(item)

    return cppmap


cdef dict map_to_dict_str_matp(cpp_map[std_string, matp] cppmap):
    pydict = {}
    cdef material._Material pymat
    cdef cpp_map[std_string, matp].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pymat = material.Material()
        pymat.mat_pointer[0] = deref(deref(mapiter).second)
        pydict[<char *> deref(mapiter).first.c_str()] = pymat
        inc(mapiter)

    return pydict

# cdef class _OldMaterialLibrary(object):
#
#    def __init__(self, lib=None, datapath="/materials", nucpath="/nucid"):
#        """Parameters
#        ----------
#        lib : dict-like, str, or None, optional
#            The data to intialize the material library with.  If this is a
#            string, it is interpreted as a path to a file.
#        datapath : str, optional
#            The path in the heirarchy to the data table in an HDF5 file.
#        nucpath : str, optional
#            The path in the heirarchy to the nuclide array in an HDF5 file.
#
#        """
#        if sys.version_info[0] >=3 and isinstance(lib, bytes):
#            lib = lib.decode()
#        cdef dict _lib = {}
#        if lib is None:
#            self._lib = _lib
#        elif isinstance(lib, collections.Mapping):
#            for key, mat in lib.items():
#                _lib[key] = ensure_material(mat)
#            self._lib = _lib
#        elif isinstance(lib, basestring):
#            self._lib = _lib
#            if lib.endswith('.json') or lib.endswith('.js'):
#                self.from_json(lib)
#            if lib.endswith('.h5') or lib.endswith('.hdf5') \
#                                   or lib.endswith('.h5m'):
#                self.from_hdf5(lib, datapath=datapath, nucpath=nucpath)
#        elif isinstance(lib, collections.Sequence):
#            for key, mat in lib:
#                _lib[key] = ensure_material(mat)
#            self._lib = _lib
#        else:
#            msg = "Could not initialize library with lib type {0!r}"
#            raise TypeError(msg.format(type(lib)))
#
#    def __contains__(self, key):
#        return key in self._lib
#
#    def __len__(self):
#        return len(self._lib)
#
#    def __iter__(self):
#        return iter(self._lib)
#
#    def __getitem__(self, key):
#        return self._lib[key]
#
#    def __setitem__(self, key, value):
#        self._lib[key] = ensure_material(value)
#
#    def __delitem__(self, key):
#        del self._lib[key]
#
#    def from_json(self, file):
#        """Loads data from a JSON file into this material library.
#
#        Parameters
#        ----------
#        file : str
#            A path to a JSON file.
#
#        """
#        cdef std_string s
#        cdef bint opened_here = False
#        cdef cpp_jsoncpp.Value jsonlib
#        cdef cpp_jsoncpp.Reader reader = cpp_jsoncpp.Reader()
#        cdef int i
#        cdef std_string key
#        cdef cpp_vector[std_string] keys
#        cdef _Material mat
#        cdef dict _lib = (<_OldMaterialLibrary> self)._lib
#        if isinstance(file, basestring):
#            file = open(file, 'r')
#            opened_here = True
#        fstr = file.read()
#        if isinstance(fstr, str):
#            fstr = fstr.encode()
#        s = std_string(<char *> fstr)
#        if opened_here:
#            file.close()
#        reader.parse(s, jsonlib)
#        keys = jsonlib.getMemberNames()
#        for i in range(len(keys)):
#            mat = Material()
#            key = keys[i]
#            (<_Material> mat).mat_pointer.load_json(jsonlib[key])
#            _lib[bytes(key.c_str()).decode()] = mat
#
#    def write_json(self, file):
#        """Writes this material library to a JSON file.
#
#        Parameters
#        ----------
#        file : str
#            A path to a JSON file.
#
#        """
#        cdef std_string s
#        cdef std_string skey
#        cdef bint opened_here = False
#        cdef cpp_jsoncpp.Value jsonlib = cpp_jsoncpp.Value(cpp_jsoncpp.objectValue)
#        cdef cpp_jsoncpp.StyledWriter writer = cpp_jsoncpp.StyledWriter()
#        for key, mat in self._lib.items():
#            key = key.encode()
#            skey = std_string(<char *> key)
#            jsonlib[skey] = (<_Material> mat).mat_pointer.dump_json()
#        s = writer.write(jsonlib)
#        if isinstance(file, basestring):
#            file = open(file, 'w')
#            opened_here = True
#        file.write(bytes(s).decode())
#        if opened_here:
#            file.close()
#
#    def from_hdf5(self, file, datapath="/materials", nucpath="/nucid"):
#        """Loads data from an HDF5 file into this material library.
#
#        Parameters
#        ----------
#        file : str
#            A path to an HDF5 file.
#        datapath : str, optional
#            The path in the heirarchy to the data table in an HDF5 file.
#        nucpath : str, optional
#            The path in the heirarchy to the nuclide array in an HDF5 file.
#
#        """
#        cdef std_string s
#        cdef cpp_jsoncpp.Reader reader = cpp_jsoncpp.Reader()
#        cdef cpp_jsoncpp.Value attribs
#        cdef int i
#        cdef _Material mat
#        cdef dict _lib = (<_OldMaterialLibrary> self)._lib
#        cdef np.ndarray mattable
#        with tb.open_file(file, 'r') as f:
#            matstable = f.get_node(datapath)[:]
#            nucs = f.get_node(nucpath)[:]
#            matsmetadata = f.get_node(datapath + '_metadata').read()
#        for i in range(len(matstable)):
#            row = matstable[i]
#            comp = dict((<int> k, v) for k, v in zip(nucs, row[3]) if v != 0.0)
#            mat = Material(comp, mass=row[0], density=row[1],
#                                    atoms_per_molecule=row[2])
#            strmetadata = "".join(map(chr, matsmetadata[i]))
#            strmetadata = strmetadata.encode()
#            s = std_string(<char *> strmetadata)
#            attribs = cpp_jsoncpp.Value()
#            reader.parse(s, attribs)
#            (<_Material> mat).mat_pointer.metadata = attribs
#            if "name" in mat.metadata:
#                name = mat.metadata["name"]
#            else:
#                name = "_" + str(i)
#            _lib[name] = mat
#
#    def write_hdf5(self, filename, datapath="/materials", nucpath="/nucid"):
#        """Writes this material library to an HDF5 file.
#
#        Parameters
#        ----------
#        filename : str
#            A path to an HDF5 file.
#        datapath : str, optional
#            The path in the heirarchy to the data table in an HDF5 file.
#        nucpath : str, optional
#            The path in the heirarchy to the nuclide array in an HDF5 file.
#
#        """
#        cdef _Material mat
#        cdef dict _lib = (<_OldMaterialLibrary> self)._lib
#        cdef set nucids = set()
#        for mat in _lib.values():
#            nucids.update(mat.comp.keys())
#        with tb.open_file(filename, 'a') as f:
#            nucgrp, nucdsname = os.path.split(nucpath)
#            f.create_array(nucgrp, nucdsname, np.array(sorted(nucids)),
#                          createparents=True)
#        for key, mat in _lib.items():
#            if "name" not in mat.metadata:
#                mat.metadata["name"] = key
#            mat.write_hdf5(filename, datapath=datapath, nucpath=nucpath)
#
# class OldMaterialLibrary(_OldMaterialLibrary, collections.MutableMapping):
#    """The material library is a collection of unique keys mapped to
#    Material objects.  This is useful for organization and declaring
#    prefernces between several sources (multiple libraries).
#    """
#    def __repr__(self):
#        libs = ["{0!r}={1!r}".format(k, m) for k, m in self.items()]
#        libs = "{" + ", ".join(libs) + "}"
#        return "pyne.material.OldMaterialLibrary({0})".format(libs)
#
#ensure_material = lambda m: m if isinstance(m, Material) else Material(m)
