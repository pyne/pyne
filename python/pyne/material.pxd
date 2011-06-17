# Cython imports
from libcpp.map cimport map as cpp_map
from cython import pointer


# Local imports
cimport std
cimport cpp_material

cdef class Material:
    cdef cpp_material.Material * mat_pointer


# Dictionary - Map Converters
ctypedef cpp_material.Material * matp
cdef cpp_map[std.string, matp] dict_to_map_str_msp(dict)
cdef dict map_to_dict_str_matp(cpp_map[std.string, matp])
