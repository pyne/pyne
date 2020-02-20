"""C++ wrapper for material class."""
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.set cimport set as std_set
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool

cimport cpp_jsoncpp

from pyne cimport cpp_material 

cdef extern from "material_library.h" namespace "pyne":

    cdef cppclass MaterialLibrary:
        # Constuctors
        MaterialLibrary()
        MaterialLibrary(std_string) except +
        MaterialLibrary(std_string, std_string) except +

        # Attributes
        unordered_map[std_string, cpp_material.Material*] material_library
        std_set[std_string] material_name
        std_set[int] nuclides
        vector[std_string] name_order

        # Methods
        void from_hdf5(std_string, std_string) except +
        
        void load_json(cpp_jsoncpp.Value) except +
        cpp_jsoncpp.Value dump_json() except +
        void from_json(std_string) except +
        void write_json(std_string) except +

        void write_hdf5(std_string, std_string) except +
        
        void add_material(cpp_material.Material) except +
        void add_material(std_string, cpp_material.Material) except +
        
        void del_material(cpp_material.Material) except +
        void del_material(std_string) except +
        void merge(MaterialLibrary*) except +
        void replace(int, cpp_material.Material) except +

        cpp_material.Material get_material(std_string) except +
        cpp_material.Material* get_material_ptr( std_string ) except +
        
        unordered_map[std_string, cpp_material.Material*] get_mat_library() except + 
        std_set[std_string] get_keylist() except +
        std_set[int] get_nuclist() except +
