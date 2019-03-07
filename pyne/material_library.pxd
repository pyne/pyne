# Cython imports
from cython import pointer

cimport cpp_material
cimport cpp_material_library

cdef class MaterialLibrary:
    cdef cpp_material_library.MaterialLibrary * _inst
    cdef cpp_set[std_string] get_matlist(self)

