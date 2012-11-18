"""C++ wrapper for rxname library."""
from libcpp.map cimport map
from libcpp.set cimport set

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport extra_types

cdef extern from "rxname.h" namespace "pyne::rxname":
    # names sets
    set[std_string] names

    # Conversion dictionaries
    map[extra_types.uint, std_string] labels
    map[extra_types.uint, std_string] id_name
    map[std_string, extra_types.uint] name_id


