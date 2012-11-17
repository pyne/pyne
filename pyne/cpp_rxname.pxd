"""C++ wrapper for rxname library."""
from libcpp.map cimport map
from libcpp.set cimport set

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string

cdef extern from "rxname.h" namespace "pyne::rxname":
    # Conversion dictionaries
    map[std_string, std_string] labels

    # Elemental string sets
    set[std_string] names

