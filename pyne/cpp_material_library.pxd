"""C++ wrapper for material class."""
from libcpp.string cimport string as std_string
from libcpp.set cimport set as std_set
from libcpp.map cimport map

from pyne cimport cpp_material 

cdef extern from "material_pyne.h" namespace "pyne":
    # Cython does not allow for typdef'ing tamplated types :(
    #ctypedef map[int, double] comp_map
    #ctypedef map[int, double].iterator comp_iter

    cdef cppclass MaterialLibrary:
        # Constuctors
        MaterialLibrary()
        MaterialLibrary(std_string) except +
        MaterialLibrary(std_string, std_string) except +

        # Attributes
        map[std_string, cpp_material] materials
        std_set[std_string] material_name
        std_set[int] nuclides


        # Methods
        void from_hdf5(char *) except +
        void from_hdf5(char *, char *) except +
        void from_hdf5(char *, char *, int) except +

        void write_hdf5(char *) except +
        void write_hdf5(char *, char *) except +
        void write_hdf5(char *, char *, char *) except +
        void write_hdf5(char *, char *, char *, int) except +
        
        void add_material(cpp_material) except +
        
        void del_material(cpp_material) except +
        void del_material(std_string) except +
        
        cpp_material get_material(cpp_material) except +
        cpp_material get_material(std_string) except +
        
        std_set[std_string] get_matlist() except +
        std_set[int] get_nuclist() except +
