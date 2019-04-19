# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.string cimport string as std_string
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython import pointer

import collections

# Local imports
cimport cpp_material
cimport pyne.stlcontainers as conv
cimport cpp_material_library

# Dictionary - Map Converters
ctypedef cpp_material.Material * matp

cdef class _MaterialLibrary:
    cdef cpp_material_library.MaterialLibrary * _inst
    cdef cpp_set[std_string] get_matlist(self)
    cdef cpp_set[int] get_nuclist(self)

cdef cpp_map[std_string, matp] dict_to_map_str_matp(dict)
cdef dict map_to_dict_str_matp(cpp_map[std_string, matp])
cdef dict map_to_dict_str_mat(cpp_map[std_string, cpp_material.Material])
