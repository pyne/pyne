"""C++ wrapper for rxname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libc.string cimport const_char
from libcpp.string cimport string as std_string
cimport extra_types

cdef extern from "rxname.h" namespace "pyne::rxname":
    # names sets
    set[std_string] names

    # Conversion dictionaries
    map[std_string, extra_types.uint32] altnames
    map[extra_types.uint32, std_string] id_name
    map[std_string, extra_types.uint32] name_id
    map[extra_types.uint32, extra_types.uint32] id_mt
    map[extra_types.uint32, extra_types.uint32] mt_id
    map[extra_types.uint32, std_string] labels
    map[extra_types.uint32, std_string] docs

    extra_types.uint32 hash(std_string) except +
    extra_types.uint32 hash(const_char *) except +

    std_string name(int) except +
    std_string name(extra_types.uint32) except +
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

    extra_types.uint32 id(int) except +
    extra_types.uint32 id(extra_types.uint) except +
    extra_types.uint32 id(char *) except +
    extra_types.uint32 id(std_string) except +
    extra_types.uint32 id(int, int) except +
    extra_types.uint32 id(int, std_string) except +
    extra_types.uint32 id(std_string, int) except +
    extra_types.uint32 id(std_string, std_string) except + 
    extra_types.uint32 id(int, int, std_string) except +
    extra_types.uint32 id(int, std_string, std_string) except +
    extra_types.uint32 id(std_string, int, std_string) except +
    extra_types.uint32 id(std_string, std_string, std_string) except + 

    extra_types.uint32 mt(int) except +
    extra_types.uint32 mt(extra_types.uint) except +
    extra_types.uint32 mt(char *) except +
    extra_types.uint32 mt(std_string) except +
    extra_types.uint32 mt(int, int) except +
    extra_types.uint32 mt(int, std_string) except +
    extra_types.uint32 mt(std_string, int) except +
    extra_types.uint32 mt(std_string, std_string) except + 
    extra_types.uint32 mt(int, int, std_string) except +
    extra_types.uint32 mt(int, std_string, std_string) except +
    extra_types.uint32 mt(std_string, int, std_string) except +
    extra_types.uint32 mt(std_string, std_string, std_string) except + 

    std_string label(int) except +
    std_string label(extra_types.uint32) except +
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

    std_string doc(int) except +
    std_string doc(extra_types.uint32) except +
    std_string doc(char *) except +
    std_string doc(std_string) except +
    std_string doc(int, int) except +
    std_string doc(int, std_string) except +
    std_string doc(std_string, int) except +
    std_string doc(std_string, std_string) except + 
    std_string doc(int, int, std_string) except +
    std_string doc(int, std_string, std_string) except +
    std_string doc(std_string, int, std_string) except +
    std_string doc(std_string, std_string, std_string) except + 

    int child(int, int, std_string) except +
    int child(int, std_string, std_string) except +
    int child(std_string, int, std_string) except +
    int child(std_string, std_string, std_string) except + 

    int parent(int, int, std_string) except +
    int parent(int, std_string, std_string) except +
    int parent(std_string, int, std_string) except +
    int parent(std_string, std_string, std_string) except + 
