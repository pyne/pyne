"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std

cdef extern from "../cpp/data.h" namespace "pyne":
    # nuc_weight Functions
    map[int, double] nuc_weight_map
    double nuc_weight(int) except +
    double nuc_weight(char *) except +
    double nuc_weight(std.string) except +

