# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.string cimport string as std_string
from libcpp.unordered_map cimport unordered_map as cpp_umap
from libcpp.set cimport set as cpp_set
from cython import pointer
from libcpp.memory cimport shared_ptr


import collections

# Local imports
cimport cpp_material
cimport pyne.stlcontainers as conv
cimport cpp_material_library

# Dictionary - Map Converters
ctypedef cpp_material.Material * matp

cdef class _MaterialLibrary:
    cdef cpp_material_library.MaterialLibrary * _inst
    cdef cpp_set[std_string] get_keylist(self)
    cdef cpp_set[int] get_nuclist(self)

cdef cpp_umap[std_string, shared_ptr[cpp_material.Material]] dict_to_map_str_matp(dict)
cdef dict matlib_to_dict_str_matp(cpp_material_library.MaterialLibrary)
