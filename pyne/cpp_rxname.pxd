"""C++ wrapper for rxname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libc.string cimport const_char

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
    map[extra_types.uint, std_string] id_name
    map[std_string, extra_types.uint] name_id
    map[extra_types.uint, extra_types.uint] id_mt
    map[extra_types.uint, extra_types.uint] mt_id
    map[extra_types.uint, std_string] labels
    map[extra_types.uint, std_string] docs

    extra_types.uint hash(std_string) except +
    extra_types.uint hash(const_char *) except +

    std_string name(int) except +
    std_string name(extra_types.uint) except +
    std_string name(char *) except +
    std_string name(std_string) except +
    std_string name(int, int) except +
    std_string name(int, std_string) except +
    std_string name(std_string, int) except +
    std_string name(std_string, std_string) except + 
    std_string name(int, int, std_string) except +
    std_string name(int, std_string, std_string) except +
    std_string name(std_string, int, std_string) except +
    std_string name(std_string, std_string, std_string) except + 

    extra_types.uint id(int) except +
    extra_types.uint id(extra_types.uint) except +
    extra_types.uint id(char *) except +
    extra_types.uint id(std_string) except +
    extra_types.uint id(int, int) except +
    extra_types.uint id(int, std_string) except +
    extra_types.uint id(std_string, int) except +
    extra_types.uint id(std_string, std_string) except + 
    extra_types.uint id(int, int, std_string) except +
    extra_types.uint id(int, std_string, std_string) except +
    extra_types.uint id(std_string, int, std_string) except +
    extra_types.uint id(std_string, std_string, std_string) except + 

    extra_types.uint mt(int) except +
    extra_types.uint mt(extra_types.uint) except +
    extra_types.uint mt(char *) except +
    extra_types.uint mt(std_string) except +
    extra_types.uint mt(int, int) except +
    extra_types.uint mt(int, std_string) except +
    extra_types.uint mt(std_string, int) except +
    extra_types.uint mt(std_string, std_string) except + 
    extra_types.uint mt(int, int, std_string) except +
    extra_types.uint mt(int, std_string, std_string) except +
    extra_types.uint mt(std_string, int, std_string) except +
    extra_types.uint mt(std_string, std_string, std_string) except + 

    std_string label(int) except +
    std_string label(extra_types.uint) except +
    std_string label(char *) except +
    std_string label(std_string) except +
    std_string label(int, int) except +
    std_string label(int, std_string) except +
    std_string label(std_string, int) except +
    std_string label(std_string, std_string) except + 
    std_string label(int, int, std_string) except +
    std_string label(int, std_string, std_string) except +
    std_string label(std_string, int, std_string) except +
    std_string label(std_string, std_string, std_string) except + 
