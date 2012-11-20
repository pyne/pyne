"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string

cdef extern from "nucname.h" namespace "pyne::nucname":
    # Conversion dictionaries
    map[std_string, int] name_zz
    map[int, std_string] zz_name

    # Elemental string sets
    set[std_string] LAN
    set[std_string] ACT
    set[std_string] TRU
    set[std_string] MA
    set[std_string] FP

    # Elemental integer sets
    set[int] lan
    set[int] act
    set[int] tru
    set[int] ma
    set[int] fp

    # Current Form
    std_string current_form(int) except +
    std_string current_form(char *) except +
    std_string current_form(std_string) except +

    # zzaaam Functions
    int zzaaam(int) except +
    int zzaaam(char *) except +
    int zzaaam(std_string) except +

    # name Functions
    std_string name(int) except +
    std_string name(char *) except +
    std_string name(std_string) except +

    # MCNP Functions 
    int mcnp(int) except + 
    int mcnp(char *) except + 
    int mcnp(std_string) except + 

    # Serpent Functions
    std_string serpent(int) except +
    std_string serpent(char *) except +
    std_string serpent(std_string) except +

    # NIST Functions
    std_string nist(int) except +
    std_string nist(char *) except +
    std_string nist(std_string) except +

    # Cinder Functions 
    int cinder(int) except + 
    int cinder(char *) except + 
    int cinder(std_string) except + 

