"""C++ wrapper for nucname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string as std_string

cdef extern from "nucname.h" namespace "pyne::nucname":
    # Conversion dictionaries
    map[std_string, int] name_zz
    map[int, std_string] zz_name
    map[int, int] state_id_map

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

    # isnuclide
    bint isnuclide(int) except +
    bint isnuclide(char *) except +
    bint isnuclide(std_string) except +

    # iselement
    bint iselement(int) except +
    bint iselement(char *) except +
    bint iselement(std_string) except +

    # id functions
    int id(int) except +
    int id(char *) except +
    int id(std_string) except +

    # name Functions
    std_string name(int) except +
    std_string name(char *) except +
    std_string name(std_string) except +

    # znum functions
    int znum(int) except +
    int znum(char *) except +
    int znum(std_string) except +

    # anum functions
    int anum(int) except +
    int anum(char *) except +
    int anum(std_string) except +

    # znum functions
    int snum(int) except +
    int snum(char *) except +
    int snum(std_string) except +

    # zzaaam Functions
    int zzaaam(int) except +
    int zzaaam(char *) except +
    int zzaaam(std_string) except +
    int zzaaam_to_id(int) except +
    int zzaaam_to_id(char *) except +
    int zzaaam_to_id(std_string) except +

    # zzzaaa Functions
    int zzzaaa(int) except +
    int zzzaaa(char *) except +
    int zzzaaa(std_string) except +
    int zzzaaa_to_id(int) except +
    int zzzaaa_to_id(char *) except +
    int zzzaaa_to_id(std_string) except +

    # zzllaaam Functions
    std_string zzllaaam(int) except +
    std_string zzllaaam(char *) except +
    std_string zzllaaam(std_string) except +
    int zzllaaam_to_id(char *) except +
    int zzllaaam_to_id(std_string) except +

    # MCNP Functions 
    int mcnp(int) except + 
    int mcnp(char *) except + 
    int mcnp(std_string) except + 
    int mcnp_to_id(int) except +
    int mcnp_to_id(char *) except +
    int mcnp_to_id(std_string) except +

    # OpenMC Functions 
    std_string openmc(int) except + 
    std_string openmc(char *) except + 
    std_string openmc(std_string) except + 
    int openmc_to_id(char *) except +
    int openmc_to_id(std_string) except +

    # FLUKA Functions 
    std_string fluka(int) except + 
    int fluka_to_id(char *) except + 
    int fluka_to_id(std_string) except + 

    # Serpent Functions
    std_string serpent(int) except +
    std_string serpent(char *) except +
    std_string serpent(std_string) except +
    #int serpent_to_id(int) except +
    int serpent_to_id(char *) except +
    int serpent_to_id(std_string) except +

    # NIST Functions
    std_string nist(int) except +
    std_string nist(char *) except +
    std_string nist(std_string) except +
    #int nist_to_id(int) except +
    int nist_to_id(char *) except +
    int nist_to_id(std_string) except +

    # Cinder Functions 
    int cinder(int) except + 
    int cinder(char *) except + 
    int cinder(std_string) except + 
    int cinder_to_id(int) except +
    int cinder_to_id(char *) except +
    int cinder_to_id(std_string) except +

    # ALARA Functions
    std_string alara(int) except +
    std_string alara(char *) except +
    std_string alara(std_string) except +
    #int alara_to_id(int) except +
    int alara_to_id(char *) except +
    int alara_to_id(std_string) except +

    # SZA Functions
    int sza(int) except +
    int sza(char *) except +
    int sza(std_string) except +
    int sza_to_id(int) except +
    int sza_to_id(char *) except +
    int sza_to_id(std_string) except +

    # Groundstate Functions
    int groundstate(int) except +
    int groundstate(char *) except +
    int groundstate(std_string) except +
    
    # State id Functions
    int state_id_to_id(int state) except +
    int id_to_state_id(int nuc_id) except +

    # ENSDF id Functions
    int ensdf_to_id(char *) except +
