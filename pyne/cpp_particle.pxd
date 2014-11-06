"""C++ wrapper for particle library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libc.string cimport const_char
from libcpp.string cimport string as std_string
cimport extra_types

cdef extern from "particle.h" namespace "pyne::particle":
    # names sets
    set[std_string] names
    set[int] pdc_nums

    # Conversion dictionaries
    map[std_string, int] altnames
    map[int, std_string] id_name
    map[std_string, int] name_id
    map[std_string, std_string] docs

    # functions
    # is_hydrogen
    bint is_hydrogen(int) except +
    bint is_hydrogen(char *) except +
    bint is_hydrogen(std_string) except +
    # is_heavy_ion
    bint is_heavy_ion(int) except +
    bint is_heavy_ion(char *) except +
    bint is_heavy_ion(std_string) except +
    # is_valid
    bint is_valid(int) except +
    bint is_valid(char *) except +
    bint is_valid(std_string) except +
    # pdc_number
    int id(int) except +
    int id(char *) except +
    int id(std_string) except +
    # mcnp
    std_string mcnp(int) except +
    std_string mcnp(char *) except +
    std_string mcnp(std_string) except +
    # mcnp6
    std_string mcnp6(int) except +
    std_string mcnp6(char *) except +
    std_string mcnp6(std_string) except +
    # fluka
    std_string fluka(int) except +
    std_string fluka(char *) except +
    std_string fluka(std_string) except +
    # geant4
    std_string geant4(int) except +
    std_string geant4(char *) except +
    std_string geant4(std_string) except +
    # name
    std_string name(int) except +
    std_string name(char *) except +
    std_string name(std_string) except +
    # describe
    std_string describe(int) except +
    std_string describe(char *) except +
    std_string describe(std_string) except +

