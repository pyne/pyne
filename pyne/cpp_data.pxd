"""C++ wrapper for data library."""
from libcpp.string cimport string as std_string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.utility cimport pair
from libcpp cimport bool
from libcpp.vector cimport vector

cimport extra_types

cdef extern from "data.h" namespace "pyne":
    # Mathematical constants
    double pi
    double N_A
    double barns_per_cm2
    double cm2_per_barn
    double sec_per_day

    # hash map and initialization function
    map[std_string, std_string] data_checksums
    
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

    # fission product data
    map[pair[int, int], double] wimsdfpy_data
    double fpyield(pair[int, int], int, bool) except +
    double fpyield(int, int, int, bool) except +
    double fpyield(char *, char *, int, bool) except +
    double fpyield(std_string, std_string, int, bool) except +

    # decay data functions
    map[int, double] half_life_map
    double half_life(int) except +
    double half_life(char *) except +
    double half_life(std_string) except +

    map[int, double] decay_const_map
    double decay_const(int) except +
    double decay_const(char *) except +
    double decay_const(std_string) except +

    map[pair[int, int], double] branch_ratio_map
    double branch_ratio(pair[int, int]) except +
    double branch_ratio(int, int) except +
    double branch_ratio(char *, char *) except +
    double branch_ratio(std_string, std_string) except +

    map[int, double] state_energy_map
    double state_energy(int) except +
    double state_energy(char *) except +
    double state_energy(std_string) except +

    map[int, set[int]] decay_children_map
    set[int] decay_children(int) except +
    set[int] decay_children(char *) except +
    set[int] decay_children(std_string) except +

    int metastable_id(int, int) except +
    int metastable_id(int) except +

    cdef struct decay_struct:
        int parent
        int daughter
        char decay[5]
        double half_life
        double half_life_error
        double branch_ratio
        double photon_branch_ratio
        double photon_branch_ratio_error
        double beta_branch_ratio
        double beta_branch_ratio_error

    pair[double,double] decay_half_life(pair[int, int]) except +
    vector[pair[double, double]] decay_half_lifes(int) except +
    double decay_branch_ratio(pair[int, int] from_to) except +
    vector[double] decay_branch_ratios(int) except +
    pair[double,double] decay_photon_branch_ratio(pair[int, int] from_to) except +
    vector[pair[double, double]] decay_photon_branch_ratios(int) except +
    pair[double,double] decay_beta_branch_ratio(pair[int, int] from_to) except +
    vector[pair[double, double]] decay_beta_branch_ratios(int) except +

    cdef struct gamma_struct:
        double energy
        double energy_err
        double photon_intensity
        double photon_intensity_err
        double conv_intensity
        double conv_intensity_err
        double total_intensity
        double total_intensity_err
        int from_nuc
        int to_nuc
        int parent_nuc
        double k_conv_e
        double l_conv_e
        double m_conv_e

    vector[pair[double, double]] gamma_energy(int parent) except +
    vector[pair[double, double]] gamma_photon_intensity(int parent) except +
    vector[pair[double, double]] gamma_conversion_intensity(int parent) except +
    vector[pair[double, double]] gamma_total_intensity(int parent) except +
    vector[pair[int, int]] gamma_from_to(int parent) except +
    vector[pair[int, int]] gamma_from_to(double energy, double error) except +
    vector[int] gamma_parent(double energy, double error) except +

    cdef struct alpha_struct:
        double energy
        double intensity
        int from_nuc
        int to_nuc

    vector[double] alpha_energy(int parent) except +
    vector[double] alpha_intensity(int parent) except +
    vector[int] alpha_parent(double energy, double error) except +
    vector[int] alpha_daughter(double energy, double error) except +
    vector[int] alpha_daughter(int parent) except +

    cdef struct beta_struct:
        double endpoint_energy
        double avg_energy
        double intensity
        int from_nuc
        int to_nuc

    vector[double] beta_endpoint_energy(int parent) except +
    vector[double] beta_average_energy(int parent) except +
    vector[double] beta_intensity(int parent) except +
    vector[int] beta_parent(double energy, double error) except +
    vector[int] beta_daughter(double energy, double error) except +
    vector[int] beta_daughter(int parent) except +


    cdef struct ecbp_struct:
        double endpoint_energy
        double avg_energy
        double beta_plus_intensity
        double ec_intensity
        int from_nuc
        int to_nuc
        double k_conv_e
        double l_conv_e
        double m_conv_e

    vector[double] ecbp_endpoint_energy(int parent) except +
    vector[double] ecbp_average_energy(int parent) except +
    vector[double] ec_intensity(int parent) except +
    vector[double] bp_intensity(int parent) except +
    vector[int] ecbp_parent(double energy, double error) except +
    vector[int] ecbp_daughter(double energy, double error) except +
    vector[int] ecbp_daughter(int parent) except +
