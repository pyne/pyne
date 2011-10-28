"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std

cdef extern from "nucname.h" namespace "pyne::nucname":
    # Conversion dictionaries
    map[std.string, int] name_zz
    map[int, std.string] zz_name

    # Elemental string sets
    set[std.string] LAN
    set[std.string] ACT
    set[std.string] TRU
    set[std.string] MA
    set[std.string] FP

    # Elemental integer sets
    set[int] lan
    set[int] act
    set[int] tru
    set[int] ma
    set[int] fp

    # Current Form
    std.string current_form(int) except +
    std.string current_form(char *) except +
    std.string current_form(std.string) except +

    # zzaaam Functions
    int zzaaam(int) except +
    int zzaaam(char *) except +
    int zzaaam(std.string) except +

    # name Functions
    std.string name(int) except +
    std.string name(char *) except +
    std.string name(std.string) except +

    # MCNP Functions 
    int mcnp(int) except + 
    int mcnp(char *) except + 
    int mcnp(std.string) except + 

    # Serpent Functions
    std.string serpent(int) except +
    std.string serpent(char *) except +
    std.string serpent(std.string) except +

    # NIST Functions
    std.string nist(int) except +
    std.string nist(char *) except +
    std.string nist(std.string) except +

    # Cinder Functions 
    int cinder(int) except + 
    int cinder(char *) except + 
    int cinder(std.string) except + 

