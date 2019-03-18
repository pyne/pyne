# Local imports
cimport cpp_jsoncpp


#
# Json containers
#

cdef class Value:
    cdef cpp_jsoncpp.Value * _inst
    cdef public bint _view

    cdef __set_instance__(self, cpp_jsoncpp.Value new_inst)

cdef class Reader:
    cdef cpp_jsoncpp.Reader * _inst

cdef class FastWriter:
    cdef cpp_jsoncpp.FastWriter * _inst

cdef class StyledWriter:
    cdef cpp_jsoncpp.StyledWriter * _inst

cdef class CustomWriter:
    cdef cpp_jsoncpp.CustomWriter * _inst
