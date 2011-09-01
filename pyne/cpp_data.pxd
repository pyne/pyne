"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std
cimport extra_types

cdef extern from "../cpp/data.h" namespace "pyne":
    # nuc_weight functions
    map[int, double] nuc_weight_map
    double nuc_weight(int) except +
    double nuc_weight(char *) except +
    double nuc_weight(std.string) except +


    # Scattering length functions
    map[int, extra_types.complex_t] b_coherent_map
    extra_types.complex_t b_coherent(int) except +
    extra_types.complex_t b_coherent(char *) except +
    extra_types.complex_t b_coherent(std.string) except +

    map[int, extra_types.complex_t] b_incoherent_map
    extra_types.complex_t b_incoherent(int) except +
    extra_types.complex_t b_incoherent(char *) except +
    extra_types.complex_t b_incoherent(std.string) except +

    map[int, double] b_map
    double b(int) except +
    double b(char *) except +
    double b(std.string) except +


    # decay data functions
    map[int, double] half_life_map
    double half_life(int) except +
    double half_life(char *) except +
    double half_life(std.string) except +

    map[int, double] decay_const_map
    double decay_const(int) except +
    double decay_const(char *) except +
    double decay_const(std.string) except +


