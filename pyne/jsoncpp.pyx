"""Python wrapper for jsoncpp."""
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free

# Python imports
import collections

# local imports
cimport std
cimport cpp_jsoncpp


cdef class Value:
    def __cinit__(self):
        """Value C++ constuctor."""
        self._inst = new cpp_jsoncpp.Value()

    def __dealloc__(self):
        """Value C++ destructor."""
        del self._inst


    #
    # Class Attributes
    #

