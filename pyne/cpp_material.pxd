"""C++ wrapper for material class."""
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.set cimport set as std_set
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool

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
        Material(map[int, double], double, double, double) except +
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
        double atoms_per_molecule
        cpp_jsoncpp.Value metadata

        # Methods
        void norm_comp() except +
        std_string openmc(std_string, int) except +
        std_string mcnp(std_string) except +
        std_string mcnp(std_string, bool) except +
        std_string get_uwuw_name() except +
        std_string phits(std_string) except +
        std_string phits(std_string, bool) except +
        std_string fluka(int, std_string) except +
        bool not_fluka_builtin(std_string) except +
        std_string fluka_material_str(int) except +
        std_string fluka_material_component(int, int, std_string) except +
        std_string fluka_material_line(int, double, int, std_string) except +
        std_string fluka_format_field(float) except +
        std_string fluka_compound_str(int, std_string) except +
        void from_hdf5(char *, char *) except +
        void from_hdf5(char *, char *, int) except +
        void from_hdf5(char *, char *, int, int) except +

        void deprecated_write_hdf5(char *, char *, char *, float, int) except +
        void write_hdf5(char *, char *, float, int) except +

        void from_text(char *) except +

        void write_text(char *) except +

        void load_json(cpp_jsoncpp.Value) except +
        cpp_jsoncpp.Value dump_json() except +
        void from_json(char *) except +
        void write_json(char *) except +

        void normalize() except +
        map[int, double] mult_by_mass() except +
        map[int, double] activity() except +
        map[int, double] decay_heat() except +
        map[int, double] dose_per_g(std_string, int) except +
        double molecular_mass() except +
        double molecular_mass(double) except +
        Material expand_elements(std_set[int]) except +
        Material collapse_elements(std_set[int]) except +
        double mass_density() except +
        double mass_density(double) except +
        double mass_density(double, double) except +
        double number_density() except +
        double number_density(double) except +
        double number_density(double, double) except +

        # Substream Methods
        Material sub_mat(set[int]) except +
        Material set_mat(set[int], double) except +
        Material del_mat(set[int]) except +

        Material sub_range(int, int) except +
        Material set_range(int, int, double) except +
        Material del_range(int, int) except +

        Material sub_elem(int) except +
        Material sub_lan() except +
        Material sub_act() except +
        Material sub_tru() except +
        Material sub_ma() except +
        Material sub_fp() except +

        # Atom frac member functions
        map[int, double] to_atom_frac() except +
        void from_atom_frac(map[int, double]) except +

        map[int, double] to_atom_dens() except +


        vector[pair[double, double]] gammas() except +
        vector[pair[double, double]] xrays() except +
        vector[pair[double, double]] photons(bool) except +
        void from_activity(map[int, double]) except +

        Material decay(double) except +
        Material cram(vector[double]) except +
        Material cram(vector[double], int) except +

        # Operator Overloads
        Material operator+(double) except +
        Material operator+(Material) except +
        Material operator*(double) except +
        Material operator/(double) except +
