"""C++ wrapper for extra types header."""
cimport std


cdef extern from "extra_types.h" namespace "extra_types":

    ctypedef struct complex_t:
        double re
        double im
