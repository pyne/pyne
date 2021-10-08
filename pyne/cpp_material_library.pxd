"""C++ wrapper for material class."""
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.set cimport set as std_set
from libcpp.unordered_map cimport unordered_map as cpp_umap
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool
from libcpp.memory cimport shared_ptr

cimport cpp_jsoncpp

from pyne cimport cpp_material

cdef extern from "material_library.h" namespace "pyne":

    cdef cppclass MaterialLibrary:
        # Constuctors
        MaterialLibrary()
        MaterialLibrary(std_string) except +
        MaterialLibrary(std_string, std_string) except +

        # Methods
        void from_hdf5(std_string, std_string) except +

        void load_json(cpp_jsoncpp.Value) except +
        cpp_jsoncpp.Value dump_json() except +
        void from_json(std_string) except +
        void write_json(std_string) except +
        void write_openmc(std_string) except +
        void write_hdf5(std_string, std_string, bool) except +

        void add_material(cpp_material.Material) except +
        void add_material(std_string, cpp_material.Material) except +

        void del_material(cpp_material.Material) except +
        void del_material(std_string) except +
        void merge(MaterialLibrary*) except +

        cpp_material.Material get_material(std_string) except +
        shared_ptr[cpp_material.Material] get_material_ptr( std_string ) except +

        cpp_umap[std_string, shared_ptr[cpp_material.Material]] get_mat_library() except +
        std_set[std_string] get_keylist() except +
        std_set[int] get_nuclist() except +
        int size() except +
        cpp_umap[std_string, shared_ptr[cpp_material.Material]].iterator begin() except +
        cpp_umap[std_string, shared_ptr[cpp_material.Material]].iterator end() except +

