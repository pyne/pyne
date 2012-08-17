# Local imports
cimport cpp_jsoncpp


#
# Json containers
#

cdef class Value:
    cdef cpp_jsoncpp.Value * _inst
