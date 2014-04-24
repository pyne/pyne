from warnings import warn

# Local imports
cimport cpp_jsoncpp


warn(__name__ + " is not yet V&V compliant.", ImportWarning)

#
# Json containers
#

cdef class Value:
    cdef cpp_jsoncpp.Value * _inst
    cdef public bint _view

cdef class Reader:
    cdef cpp_jsoncpp.Reader * _inst

cdef class FastWriter:
    cdef cpp_jsoncpp.FastWriter * _inst

cdef class StyledWriter:
    cdef cpp_jsoncpp.StyledWriter * _inst
