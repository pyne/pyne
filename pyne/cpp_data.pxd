"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string

cimport extra_types

cdef extern from "data.h" namespace "pyne":
    # atomic_mass functions
    map[int, double] atomic_mass_map
    double atomic_mass(int) except +
    double atomic_mass(char *) except +
    double atomic_mass(std_string) except +

    # natural_abund functions
    map[int, double] natural_abund_map
    double natural_abund(int) except +
    double natural_abund(char *) except +
    double natural_abund(std_string) except +

    # Scattering length functions
    map[int, extra_types.complex_t] b_coherent_map
    extra_types.complex_t b_coherent(int) except +
    extra_types.complex_t b_coherent(char *) except +
    extra_types.complex_t b_coherent(std_string) except +

    map[int, extra_types.complex_t] b_incoherent_map
    extra_types.complex_t b_incoherent(int) except +
    extra_types.complex_t b_incoherent(char *) except +
    extra_types.complex_t b_incoherent(std_string) except +

    map[int, double] b_map
    double b(int) except +
    double b(char *) except +
    double b(std_string) except +


    # decay data functions
    map[int, double] half_life_map
    double half_life(int) except +
    double half_life(char *) except +
    double half_life(std_string) except +

    map[int, double] decay_const_map
    double decay_const(int) except +
    double decay_const(char *) except +
    double decay_const(std_string) except +


