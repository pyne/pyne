"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std

cdef extern from "../isoname.h" namespace "isoname":
    # Conversion dictionaries
    map[std.string, int] LLzz
    map[int, std.string] zzLL

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
    std.string CurrentForm(int) except +
    std.string CurrentForm(std.string) except +

    # LLAAAM_2_* Functions
    int LLAAAM_2_zzaaam(std.string) except +
    int LLAAAM_2_MCNP(std.string) except +

    # zzaaam_2_* Functions
    std.string zzaaam_2_LLAAAM(int) except +
    int zzaaam_2_MCNP(int) except +

    # MCNP_2_* Functions 
    int MCNP_2_zzaaam(int) except + 
    std.string MCNP_2_LLAAAM(int) except + 

    # mixed_2_*_ Functions
    int mixed_2_zzaaam(std.string) except +
    int mixed_2_zzaaam(int) except +
    std.string mixed_2_LLAAAM(std.string) except +
    std.string mixed_2_LLAAAM(int) except +
    int mixed_2_MCNP(std.string) except +
    int mixed_2_MCNP(int) except +

    # Helper Functions
    double nuc_weight_zzaaam(int) except +
    double nuc_weight(int) except +
    double nuc_weight(std.string) except +

