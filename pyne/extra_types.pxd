"""C++ wrapper for extra types header."""

ctypedef unsigned int uint

cdef extern from "extra_types.h" namespace "extra_types":

    ctypedef struct complex_t:
        double re
        double im

cdef complex_t py2c_complex(object pyv)
