# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from cython import pointer

# Local imports
cimport std
cimport cpp_material
cimport pyne.stlconverters as conv


cdef cpp_map[int, double] dict_to_comp(dict)

cdef class _Material:
    cdef cpp_material.Material * mat_pointer

# Dictionary - Map Converters
ctypedef cpp_material.Material * matp
cdef cpp_map[std.string, matp] dict_to_map_str_matp(dict)
cdef dict map_to_dict_str_matp(cpp_map[std.string, matp])
