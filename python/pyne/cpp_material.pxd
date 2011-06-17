"""C++ wrapper for material class."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../../cpp/material.h" namespace pyne:
    ctypedef map[int, double] comp_map
    ctypedef map[int, double].iterator comp_iter

    cdef cppclass Material:
        # Constuctors
        Material()
        Material(comp_map, float, std.string) except +
        Material(char *, float, std.string) except +

        # Attributes
        comp_map comp
        double mass
        std.string name

        # Methods
        void norm_comp() except +
        void load_from_hdf5(char *, char *, int) except +
        void load_from_text(char *) except +

        void normalize() except +
        map[int, double] mult_by_mass() except +
        double atomic_weight() except +

        # Substream Methods
        Material sub_mat(set[int], std.string) except +
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
