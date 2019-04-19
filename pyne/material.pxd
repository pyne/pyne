# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.string cimport string as std_string
from libcpp.map cimport map as cpp_map
from cython import pointer

import collections

# Local imports
cimport cpp_material
cimport pyne.stlcontainers as conv

cdef cpp_map[int, double] dict_to_comp(dict)

cdef class _Material:
    cdef cpp_material.Material * mat_pointer
    cdef public bint _free_mat
#
# Material containers
#
ctypedef cpp_material.Material * matp


# (Str, Material)
cdef class MapIterStrMaterial(object):
    cdef cpp_map[std_string, matp].iterator * iter_now
    cdef cpp_map[std_string, matp].iterator * iter_end
    cdef void init(MapIterStr, cpp_map[std_string, matp] *)

cdef class _MapStrMaterial:
    cdef cpp_map[std_string, matp] * map_ptr
    cdef public bint _free_map
    cdef dict _cache

#cdef class _MaterialLibrary(object):
#    cdef dict _lib
