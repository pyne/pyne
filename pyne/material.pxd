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
    cdef public bint _free_mat

#
# Material containers
#
ctypedef cpp_material.Material * matp

# Dictionary - Map Converters
cdef cpp_map[std.string, matp] dict_to_map_str_matp(dict)
cdef dict map_to_dict_str_matp(cpp_map[std.string, matp])

# (Str, Material)
cdef class MapIterStrMaterial(object):
    cdef cpp_map[std.string, matp].iterator * iter_now
    cdef cpp_map[std.string, matp].iterator * iter_end
    cdef void init(MapIterStr, cpp_map[std.string, matp] *)

cdef class _MapStrMaterial:
    cdef cpp_map[std.string, matp] * map_ptr
    cdef public bint _free_map
    cdef dict _cache
