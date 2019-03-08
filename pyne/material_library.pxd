# Cython imports
from cython import pointer

from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string

cimport cpp_material
cimport cpp_material_library

cdef class MaterialLibrary:
    cdef cpp_material_library.MaterialLibrary * _inst
    cdef cpp_set[std_string] get_matlist(self)
    cdef cpp_set[int] get_nuclist(self)

