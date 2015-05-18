"""C++ wrapper for data library."""
from libcpp cimport bool
from libcpp.utility cimport pair
from libcpp.string cimport string as std_string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport extra_types

cdef extern from "data.h" namespace "pyne":
    # Mathematical constants
    double pi
    double N_A
    double barns_per_cm2
    double cm2_per_barn
    double sec_per_day
    double MeV_per_K
    double MeV_per_MJ
    double Bq_per_Ci 
    double Ci_per_Bq
   # hash map and initialization function
    map[std_string, std_string] data_checksums

    # atomic_mass functions
    map[int, double] atomic_mass_map
    double atomic_mass(int) except +
    double atomic_mass(char *) except +
    double atomic_mass(std_string) except +

    # simple_xs functions
    map[std_string, map[int, map[int, double]]] simple_xs_map
    double simple_xs(int nuc, int rx_id, std_string energy) except +
    double simple_xs(int nuc, std_string rx_id, std_string energy) except +
    double simple_xs(std_string nuc, int rx_id, std_string energy) except +
    double simple_xs(std_string nuc, std_string rx_id, std_string energy) except +

    # natural_abund functions
    map[int, double] natural_abund_map
    double natural_abund(int) except +
    double natural_abund(char *) except +
    double natural_abund(std_string) except +

    # q_val functions
    map[int, double] q_val_map
    double q_val(int) except +
    double q_val(char *) except +
    double q_val(std_string) except +

    map[int, double] gamma_frac_map
    double gamma_frac(int) except +
    double gamma_frac(char *) except +
    double gamma_frac(std_string) except +

    # ext_air_dose functions
    double ext_air_dose(int, int) except +
    double ext_air_dose(char *, int) except +
    double ext_air_dose(std_string, int) except +
    double dose_ratio(int, int) except +
    double dose_ratio(char *, int) except +
    double dose_ratio(std_string, int) except +

    # ext_soil_dose functions
    double ext_soil_dose(int, int) except +
    double ext_soil_dose(char *, int) except +
    double ext_soil_dose(std_string, int) except +

    # ingest_dose functions
    double ingest_dose(int, int) except +
    double ingest_dose(char *, int) except +
    double ingest_dose(std_string, int) except +
    double dose_fluid_frac(int, int) except +
    double dose_fluid_frac(char *, int) except +
    double dose_fluid_frac(std_string, int) except +

    # inhale_dose functions
    double inhale_dose(int, int) except +
    double inhale_dose(char *, int) except +
    double inhale_dose(std_string, int) except +
    std_string dose_lung_model(int, int) except +
    std_string dose_lung_model(char *, int) except +
    std_string dose_lung_model(std_string, int) except +

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

    # atomic data functions
    vector[pair[double, double]] calculate_xray_data(int, double,
                                                     double) except +

    # decay data functions
    double half_life(int) except +
    double half_life(char *) except +
    double half_life(std_string) except +

    double decay_const(int) except +
    double decay_const(char *) except +
    double decay_const(std_string) except +

    double branch_ratio(pair[int, int]) except +
    double branch_ratio(int, int) except +
    double branch_ratio(char *, char *) except +
    double branch_ratio(std_string, std_string) except +

    double state_energy(int) except +
    double state_energy(char *) except +
    double state_energy(std_string) except +

    set[int] decay_children(int) except +
    set[int] decay_children(char *) except +
    set[int] decay_children(std_string) except +

    int metastable_id(int, int) except +
    int metastable_id(int) except +

    int id_from_level(int, double) except +
    int id_from_level(int, double, std_string) except +

    #ENSDF data functions
    pair[double,double] decay_half_life(pair[int, int]) except +
    vector[pair[double, double]] decay_half_lifes(int) except +
    pair[double,double] decay_branch_ratio(pair[int, int] from_to) except +
    vector[double] decay_branch_ratios(int) except +
    pair[double,double] decay_photon_branch_ratio(pair[int, int] from_to) except +
    vector[pair[double, double]] decay_photon_branch_ratios(int) except +
    pair[double,double] decay_beta_branch_ratio(pair[int, int] from_to) except +
    vector[pair[double, double]] decay_beta_branch_ratios(int) except +
    vector[int] decay_data_children(int) except +

    vector[pair[double, double]] gamma_energy(int parent) except +
    vector[pair[double, double]] gamma_energy(double energy,
                                              double error) except +
    vector[pair[double, double]] gamma_photon_intensity(int parent) except +
    vector[pair[double, double]] gamma_photon_intensity(double energy,
                                                        double error) except +
    vector[pair[double, double]] gamma_conversion_intensity(int parent) except +
    vector[pair[double, double]] gamma_total_intensity(int parent) except +
    vector[pair[int, int]] gamma_from_to(int parent) except +
    vector[pair[int, int]] gamma_from_to(double energy, double error) except +
    vector[pair[int, int]] gamma_parent_child(double energy, double error) except +
    vector[int] gamma_parent(double energy, double error) except +
    vector[int] gamma_child(double energy, double error) except +
    vector[int] gamma_child(int parent) except +
    vector[pair[double, double]] gamma_xrays(int parent) except +

    vector[double] alpha_energy(int parent) except +
    vector[double] alpha_intensity(int parent) except +
    vector[int] alpha_parent(double energy, double error) except +
    vector[int] alpha_child(double energy, double error) except +
    vector[int] alpha_child(int parent) except +

    vector[double] beta_endpoint_energy(int parent) except +
    vector[double] beta_average_energy(int parent) except +
    vector[double] beta_intensity(int parent) except +
    vector[int] beta_parent(double energy, double error) except +
    vector[int] beta_child(double energy, double error) except +
    vector[int] beta_child(int parent) except +

    vector[double] ecbp_endpoint_energy(int parent) except +
    vector[double] ecbp_average_energy(int parent) except +
    vector[double] ec_intensity(int parent) except +
    vector[double] bp_intensity(int parent) except +
    vector[int] ecbp_parent(double energy, double error) except +
    vector[int] ecbp_child(double energy, double error) except +
    vector[int] ecbp_child(int parent) except +
    vector[pair[double, double]] ecbp_xrays(int parent) except +
