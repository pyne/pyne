# Local imports
cimport cpp_jsoncpp


#
# Json containers
#

cdef class Value:
    cdef cpp_jsoncpp.Value * _inst
    cdef public bint _view

cdef class Reader:
    cdef cpp_jsoncpp.Reader * _inst
