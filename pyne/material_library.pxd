# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.string cimport string as std_string
from libcpp.map cimport map as cpp_map
from cython import pointer

import collections

# Local imports
cimport cpp_material_librar_library
cimport pyne.stlcontainers as conv

cdef cpp_map[int, double] dict_to_comp(dict)

cdef class _MaterialLibrary:
    cdef cpp_material_library.MaterialLibrary * matlibrary_pointer
    cdef public bint _free_mat

#
# MaterialLibrary containers
#
ctypedef cpp_material_library.MaterialLibrary * matp

# Dictionary - Map Converters
cdef cpp_map[std_string, matp] dict_to_map_str_matp(dict)
cdef dict map_to_dict_str_matp(cpp_map[std_string, matp])

# (Str, Material)
cdef class MapIterStrMaterialLibrary(object):
    cdef cpp_map[std_string, matp].iterator * iter_now
    cdef cpp_map[std_string, matp].iterator * iter_end
    cdef void init(MapIterStr, cpp_map[std_string, matp] *)

cdef class _MapStrMaterialLibrary:
    cdef cpp_map[std_string, matp] * map_ptr
    cdef public bint _free_map
    cdef dict _cache

