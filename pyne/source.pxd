from libcpp.string cimport string as std_string

cimport cpp_source

cdef class PointSource:
    cdef cpp_source.PointSource * _inst
    cdef public bint _free_inst
