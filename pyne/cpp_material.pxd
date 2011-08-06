"""C++ wrapper for material class."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../cpp/material.h" namespace "pyne":
    # Cython does not allow for typdef'ing tamplated types :( 
    #ctypedef map[int, double] comp_map
    #ctypedef map[int, double].iterator comp_iter

    cdef cppclass Material:
        # Constuctors
        Material()
        Material(map[int, double], float, std.string) except +
        Material(char *, float, std.string) except +

        # Attributes
        map[int, double] comp
        double mass
        std.string name

        # Methods
        void norm_comp() except +
        void load_from_hdf5(char *, char *, int) except +
        void load_from_text(char *) except +

        void normalize() except +
        map[int, double] mult_by_mass() except +
        double molecular_weight() except +

        # Substream Methods
        Material sub_mat(set[int], std.string) except +
        Material set_mat(set[int], double, std.string) except +
        Material del_mat(set[int], std.string) except +

        Material sub_range(int, int, std.string) except +
        Material set_range(int, int, double, std.string) except +
        Material del_range(int, int, std.string) except +

        Material sub_u(std.string) except +
        Material sub_pu(std.string) except +
        Material sub_lan(std.string) except +
        Material sub_act(std.string) except +
        Material sub_tru(std.string) except +
        Material sub_ma(std.string) except +
        Material sub_fp(std.string) except +

        # Operator Overloads
        Material operator+(double) except +
        Material operator+(Material) except +
        Material operator*(double) except +
        Material operator/(double) except +
