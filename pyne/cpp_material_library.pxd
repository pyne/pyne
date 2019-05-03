"""C++ wrapper for material class."""
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.set cimport set as std_set
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool

cimport cpp_jsoncpp

from pyne cimport cpp_material 

cdef extern from "material_library.h" namespace "pyne":
    # Cython does not allow for typdef'ing tamplated types :(
    #ctypedef map[int, double] comp_map
    #ctypedef map[int, double].iterator comp_iter

    cdef cppclass MaterialLibrary:
        # Constuctors
        MaterialLibrary()
        MaterialLibrary(std_string) except +
        MaterialLibrary(std_string, std_string) except +

        # Attributes
        map[std_string, cpp_material.Material*] material_library
        std_set[std_string] material_name
        std_set[int] nuclides


        # Methods
        void from_hdf5(char*) except +
        void from_hdf5(char*, char*) except +
        void from_hdf5(char*, char*, char*) except +
        void from_hdf5(char*, char*, char*, int) except +
        void from_hdf5(char*, char*, int) except +
        
        void load_json(cpp_jsoncpp.Value) except +
        cpp_jsoncpp.Value dump_json() except +
        void from_json(char *) except +
        void write_json(char *) except +

        void write_hdf5(char*, char*, char*) except +
        
        void add_material(cpp_material.Material) except +
        void add_material(char*, cpp_material.Material) except +
        
        void del_material(cpp_material.Material) except +
        void del_material(std_string) except +
        void merge(std_string) except +
        
        cpp_material.Material get_material(std_string) except +
        map[std_string, cpp_material.Material*] get_mat_library() except + 
        std_set[std_string] get_matlist() except +
        std_set[int] get_nuclist() except +
