"""C++ wrapper for material class."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std

cimport cpp_jsoncpp 

cdef extern from "material.h" namespace "pyne":
    # Cython does not allow for typdef'ing tamplated types :( 
    #ctypedef map[int, double] comp_map
    #ctypedef map[int, double].iterator comp_iter

    cdef cppclass Material:
        # Constuctors
        Material()
        Material(map[int, double]) except +
        Material(map[int, double], double) except +
        Material(map[int, double], double, double) except +
        Material(map[int, double], double, double, double, cpp_jsoncpp.Value) except +
        Material(char *) except +
        Material(char *, double) except +
        Material(char *, double, double) except +
        Material(char *, double, double, double) except +
        Material(char *, double, double, double, cpp_jsoncpp.Value) except +

        # Attributes
        map[int, double] comp
        double mass
        double density
        double atoms_per_mol
        cpp_jsoncpp.Value attrs

        # Methods
        void norm_comp() except +
        void from_hdf5(char *, char *) except +
        void from_hdf5(char *, char *, int) except +
        void from_hdf5(char *, char *, int, int) except +

        void write_hdf5(char *, char *, char *) except +
        void write_hdf5(char *, char *, char *, float) except +
        void write_hdf5(char *, char *, char *, float, int) except +

        void from_text(char *) except +

        void write_text(char *) except +

        void normalize() except +
        map[int, double] mult_by_mass() except +
        double molecular_weight(double) except +

        # Substream Methods
        Material sub_mat(set[int]) except +
        Material set_mat(set[int], double) except +
        Material del_mat(set[int]) except +

        Material sub_range(int, int) except +
        Material set_range(int, int, double) except +
        Material del_range(int, int) except +

        Material sub_u() except +
        Material sub_pu() except +
        Material sub_lan() except +
        Material sub_act() except +
        Material sub_tru() except +
        Material sub_ma() except +
        Material sub_fp() except +

        # Atom frac member functions
        map[int, double] to_atom_frac() except +
        void from_atom_frac(map[int, double]) except +

        # Operator Overloads
        Material operator+(double) except +
        Material operator+(Material) except +
        Material operator*(double) except +
        Material operator/(double) except +
