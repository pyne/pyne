/// \file data.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Implements basic nuclear data functions.

#ifndef PYNE_TEWK4A7VOFFLHDDXD5ZZ7KPXEQ
#define PYNE_TEWK4A7VOFFLHDDXD5ZZ7KPXEQ
#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "hdf5.h"
#include "hdf5_hl.h"

#ifndef PYNE_IS_AMALGAMATED
#include "h5wrap.h"
#include "extra_types.h"
#include "pyne.h"
#include "nucname.h"
#endif

namespace pyne
{
  /// \name Mathematical and Physical Constants
  /// \{
  extern const double pi;   ///< pi = 3.14159265359
  extern const double N_A;  ///< Avogadro's Number
  extern const double barns_per_cm2;  ///< barns per cm^2
  extern const double cm2_per_barn;   ///< cm^2 per barn
  extern const double sec_per_day;    ///< seconds per day
  /// \}

  extern std::string NUC_DATA_PATH; ///< Path to the nuc_data.h5 file.

  // Mapping from nodes in nuc_data.h5 to hashes of nodes
  extern std::map<std::string, std::string> data_checksums; 
  
  /// \name Atomic Mass Data
  /// \{

  /// Mapping from nuclides in id form to their atomic masses.
  extern std::map<int, double> atomic_mass_map;

  /// a struct matching the atomic_mass table in nuc_data.h5.
  typedef struct atomic_mass_struct {
    int nuc;      ///< nuclide in id form
    double mass;  ///< nuclide atomic mass [amu]
    double error; ///< error in atomic mass [amu]
    double abund; ///< natural abundance of nuclide [atom fraction]
  } atomic_mass_struct; 

  // Loads preset dataset hashes into memory.
  std::map<std::string, std::string> get_data_checksums();
  
  /// Loads the atomic mass and natural abundance data from the nuc_data.h5 file
  /// into memory.
  void _load_atomic_mass_map(); 

  /// \brief Returns the atomic mass of a nuclide \a nuc.  
  /// 
  /// This function will first try to find the atomic mass data in the atomic_mass_map.
  /// If this map is empty, it will load the data from disk.  If the nuclide is in an
  /// excited state and not found in the map, it will give the value for the ground
  /// state nuclide.  If the nuclide simply cannot be found, the A number is returned.
  double atomic_mass(int nuc);
  /// Returns the atomic mass of a nuclide \a nuc.
  double atomic_mass(char * nuc); 
  /// Returns the atomic mass of a nuclide \a nuc.
  double atomic_mass(std::string nuc); 
  /// \}


  /// \name Natural Abundance Data
  /// \{

  /// Mapping from nuclides in id form to their natural abundances.
  extern std::map<int, double> natural_abund_map;

  /// \brief Returns the natural abundance of a nuclide \a nuc.  
  /// 
  /// This follows the same the basic rules for finding or computing the natural
  /// abundances as the atomic_mass() functions do.  However, if the nuclide cannot
  /// be found, the default value returned is 0.0.
  double natural_abund(int nuc);
  /// Returns the natural abundance of a nuclide \a nuc. 
  double natural_abund(char * nuc); 
  /// Returns the natural abundance of a nuclide \a nuc. 
  double natural_abund(std::string nuc); 
  /// \}



  /// \name Q_value Data
  /// \{

  /// Mapping from nuclides in id form to their q_values.
  extern std::map<int, double> q_val_map;

  /// a struct matching the q_value table in nuc_data.h5.
  typedef struct q_val_struct {
    int nuc;          ///< nuclide in id form
    double q_val;      ///< nuclide q_value [MeV/fission]
    double gamma_frac; ///< fraction of q that comes from gammas
  } q_val_struct; 

  /// Loads the q_value data from the nuc_data.h5 file into memory.
  void _load_q_val_map(); 

  /// \brief Returns the q_value of a nuclide \a nuc.  
  /// 
  /// This function will first try to find the q_value data in the q_val_map.
  /// If this map is empty, it will load the data from disk. If the nuclide simply 
  /// cannot be found, the default value returned is 0.0.
  double q_val(int nuc);
  /// Returns the q_value of a nuclide \a nuc.
  double q_val(char * nuc); 
  /// Returns the q_value of a nuclide \a nuc.
  double q_val(std::string nuc); 
  /// \}



  /// \name Gamma Fraction of Q Values Data
  /// \{

  /// Mapping from nuclides in id form to the fraction of Q that comes from gammas.
  extern std::map<int, double> gamma_frac_map;

  /// \brief Returns the fraction of Q that comes from gammas of a nuclide \a nuc.  
  /// 
  /// This follows the same the basic rules for finding or computing the fraction
  /// of Q that comes from gammas as the q_val() functions do.
  /// If the nuclide cannot be found, the default value returned is 0.0.
  double gamma_frac(int nuc);
  /// Returns the gamma_frac of a nuclide \a nuc. 
  double gamma_frac(char * nuc); 
  /// Returns the gamma_frac of a nuclide \a nuc. 
  double gamma_frac(std::string nuc); 
  /// \}



  /// \name Dose Factor Data
  /// \{

  /// Mapping from nuclides in id form to their dose factors (4 types).
  extern std::map<int, double> ext_air_df_map;
  extern std::map<int, double> ext_soil_df_map;
  extern std::map<int, double> ingest_df_map;
  extern std::map<int, double> inhale_df_map;
  extern std::map<int, double> ratio_map;
  extern std::map<int, double> fluid_frac_map;
  extern std::map<int, std::string> lung_mod_map;
  
  /// A struct matching the dose factor table in nuc_data.h5.
  typedef struct df_struct {
    int nuc;            ///< nuclide in id form
    double ext_air_df;  ///< nuclide ext_air dose factor [mrem/h per Ci/m^3]
    double ratio;       ///< ratio of external air df to df due to inhalation
    double ext_soil_df; ///< nuclide ext_soil dose factor [mrem/h per Ci/m^2]
    double ingest_df;   ///< nuclide dose factor due to ingestion [mrem/pCi]
    double fluid_frac;  ///< fraction of activity abosorbed in body fluids
    double inhale_df;   ///< model of lung used (time of biological half life -- D, W, or Y)
    std::string lung_mod;      ///< nuclide dose factor [mrem/h per Ci/m^3]
  } df_struct;

  /// Loads the dose factor data from the nuc_data.h5 file into memory.
  ///
  void _load_df_map(const char * source_path);

  /// \brief Returns the dose factors of a nuclide \a nuc.  
  /// 
  /// These functions will first try to find the dose factor data in the df_maps.
  /// If the maps are empty, it will load the data from disk. If the nuclide simply 
  /// cannot be found, the default value returned is 0.0.
  double ext_air_df(int nuc, int source, bool get_error);
  double ext_soil_df(int nuc, int source, bool get_error);
  double ingest_df(int nuc, int source, bool get_error);
  double inhale_df(int nuc, int source, bool get_error);
  /// Returns the dose factor of a nuclide \a nuc.
  double ext_air_df(char * nuc, int source, bool get_error); 
  double ext_soil_df(char * nuc, int source, bool get_error); 
  double ingest_df(char * nuc, int source, bool get_error); 
  double inhale_df(char * nuc, int source, bool get_error);
  /// Returns the dose factor of a nuclide \a nuc.
  double ext_air_df(std::string nuc, int source, bool get_error);
  double ext_soil_df(std::string nuc, int source, bool get_error);
  double ingest_df(std::string nuc, int source, bool get_error);
  double inhale_df(std::string nuc, int source, bool get_error);
  
  /// \brief Returns the additional info of a nuclide \a nuc.  
  /// 
  /// This follows the same the basic rules for finding or computing the dose factor
  /// as the df() functions do.
  /// If the nuclide cannot be found, the default value returned is 0.0.
  double ratio(int nuc, int source, bool get_error);
  double fluid_frac(int nuc, int source, bool get_error);
  std::string lung_mod(int nuc, int source, bool get_error);
  /// Returns the requested info of a nuclide \a nuc. 
  double ratio(char * nuc, int source, bool get_error);
  double fluid_frac(char * nuc, int source, bool get_error);
  std::string lung_mod(char * nuc, int source, bool get_error); 
  /// Returns the requested info of a nuclide \a nuc. 
  double ratio(std::string nuc, int source, bool get_error);
  double fluid_frac(std::string nuc, int source, bool get_error);
  std::string lung_mod(std::string nuc, int source, bool get_error);
  /// \}


  
  /// \name Scattering Length Data
  /// \{

  /// Mapping from nuclides in id form to their coherent scattering length.
  extern std::map<int, xd_complex_t> b_coherent_map;
  /// Mapping from nuclides in id form to their incoherent scattering length.
  extern std::map<int, xd_complex_t> b_incoherent_map;
  /// Mapping from nuclides in id form to their scattering length.
  extern std::map<int, double> b_map;

  /// a struct matching the '/neutron/scattering_lengths' table in nuc_data.h5.
  typedef struct scattering_lengths_struct {
    int nuc;  ///< nuclide in id form
    xd_complex_t b_coherent;  ///< coherent scattering length [cm]
    xd_complex_t b_incoherent;  ///< incoherent scattering length [cm]
    double xs_coherent;   ///< coherent scattering cross section
    double xs_incoherent; ///< incoherent scattering cross section
    double xs;            ///< scattering cross section
  } scattering_lengths_struct; 

  /// Loads the scattering length data from the nuc_data.h5 file into memory.
  void _load_scattering_lengths();

  /// \brief Finds the coherent scattering length [cm] for a nuclide \a nuc.
  ///
  /// This function works by first checking the b_coherent_map.  If the map is empty
  /// then the data is read in from disk.  If no data is found than a value from a
  /// nuclide with the same A number is returned instead.  If none of these exist,
  /// then the value of a nuclide with the same Z number is used.  If none of these
  /// work then 0.0 is returned.
  xd_complex_t b_coherent(int nuc);
  /// Finds the coherent scattering length [cm] for a nuclide \a nuc.
  xd_complex_t b_coherent(char * nuc);
  /// Finds the coherent scattering length [cm] for a nuclide \a nuc.
  xd_complex_t b_coherent(std::string nuc);

  /// \brief Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  ///
  /// This function works in the same way that b_coherent() does.
  xd_complex_t b_incoherent(int nuc);
  /// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  xd_complex_t b_incoherent(char * nuc);
  /// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  xd_complex_t b_incoherent(std::string nuc);

  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(int nuc);
  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(char * nuc);
  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(std::string nuc);
  /// \}


  /// \name Fission Product Yield Data
  /// \{

  /// Mapping from nuclides in id form to their scattering length.
  extern std::map<std::pair<int, int>, double> wimsdfpy_data;

  /// a struct matching the '/neutron/wimsd_fission_product' table in nuc_data.h5.
  typedef struct wimsdfpy_struct {
    int from_nuc;  ///< from nuclide in id form
    int to_nuc;  ///< from nuclide in id form
    double yields; ///< fission product yield, fraction [unitless]
  } wimsdfpy_struct; 

  /// Loads the WIMSD fission product yield data from the nuc_data.h5 file into memory.
  void _load_wimsdfpy();

  /// a struct matching the '/neutron/nds_fission_product' table in nuc_data.h5
  typedef struct ndsfpy_struct {
    int from_nuc;
    int to_nuc;
    double yield_thermal;
    double yield_thermal_err;
    double yield_fast;
    double yield_fast_err;
    double yield_14MeV;
    double yield_14MeV_err;
  } ndsfpy_struct;

  /// a struct for the nds data for fpyield
  typedef struct ndsfpysub_struct {
    double yield_thermal;
    double yield_thermal_err;
    double yield_fast;
    double yield_fast_err;
    double yield_14MeV;
    double yield_14MeV_err;
  } ndsfpysub_struct;


  extern std::map<std::pair<int, int>, ndsfpysub_struct> ndsfpy_data;

  /// Loads the NDS fission product yield data from the nuc_data.h5 file into memory.
  void _load_ndsfpy();

  /// \brief Returns the fission product yield for a parent/child nuclide pair
  ///
  /// This function works by first checking the fission yield data.  If this is
  /// empty it loads the data from disk. If the parent/child nuclide pair
  /// is still not found, then the process is assumed to be impossible
  /// and 0.0 is returned. The data source is determined by the type value
  /// as follows: 0 WIMS, 1 thermal NDS, 2 fast NDS, 3 14 MeV NDS.
  /// negative type values return error for that data type.
  double fpyield(std::pair<int, int> from_to, int source, bool get_error);
  /// Returns the fission product yield for a parent/child nuclide pair
  double fpyield(int from_nuc, int to_nuc, int source, bool get_error);
  /// Returns the fission product yield for a parent/child nuclide pair
  double fpyield(char * from_nuc, char * to_nuc, int source, bool get_error);
  /// Returns the fission product yield for a parent/child nuclide pair
  double fpyield(std::string from_nuc, std::string to_nuc, int source, bool get_error);

  /// \}


  /// \name Decay Data
  /// \{

  /// Mapping from nuclides in id form to their half lives [s].
  extern std::map<int, double> half_life_map;
  /// Mapping from nuclides in id form to their decay constants [1/s].
  extern std::map<int, double> decay_const_map;
  /// Mapping from parent/child nuclide pairs in id form to their 
  /// branch ratio [fraction].
  extern std::map<std::pair<int, int>, double> branch_ratio_map;
  /// Mapping from nuclides in id form to their excitation energy [MeV].
  extern std::map<int, double> state_energy_map;
  /// Mapping from nuclides in id form to its decay children, if any.
  extern std::map<int, std::set<int> > decay_children_map;

  /// a struct matching the '/decay/half_life_decay' table in nuc_data.h5.
  typedef struct half_life_decay_struct {
    int from_nuc; ///< parent species in id form
    int to_nuc;   ///< child species in id form
    double level; ///< decay level [MeV]
    double half_life;     ///< species half life [s]
    double decay_const;   ///< decay constant [1/s]
    double branch_ratio;  ///< decay branch ratio [fraction]
  } half_life_decay_struct;

  /// Loads the half life data from the nuc_data.h5 file into memory.
  void _load_half_life_decay();

  /// \brief Returns the half life for a nuclide \a nuc.
  ///
  /// This function works by first checking the half_life_map.  If this is empty it
  /// loads the data from disk.  If the nuclide is still not found, then the species
  /// is assumed to be stable and infinity is returned.
  double half_life(int nuc);
  /// Returns the half life for a nuclide \a nuc.
  double half_life(char * nuc);
  /// Returns the half life for a nuclide \a nuc.
  double half_life(std::string nuc);

  /// \brief Returns the decay constant for a nuclide \a nuc.
  ///
  /// This function works by first checking the decay_const_map.  If this is empty it
  /// loads the data from disk.  If the nuclide is still not found, then the species
  /// is assumed to be stable and 0.0 is returned.
  double decay_const(int nuc);
  /// Returns the decay constant for a nuclide \a nuc.
  double decay_const(char * nuc);
  /// Returns the decay constant for a nuclide \a nuc.
  double decay_const(std::string nuc);

  /// \brief Returns the branch ratio for a parent/child nuclide pair.
  ///
  /// This function works by first checking the branch_ratio_map.  If this is empty it
  /// loads the data from disk.  If the parent/child nuclide pair is still not found, 
  /// then the decay is assumed to be impossible and 0.0 is returned.
  double branch_ratio(std::pair<int, int> from_to);
  /// Returns the branch ratio for a parent/child nuclide pair.
  double branch_ratio(int from_nuc, int to_nuc);
  /// Returns the branch ratio for a parent/child nuclide pair.
  double branch_ratio(char * from_nuc, char * to_nuc);
  /// Returns the branch ratio for a parent/child nuclide pair.
  double branch_ratio(std::string from_nuc, std::string to_nuc);

  /// \brief Returns the excitation energy [MeV] of a \a nuc in a given state.
  ///
  /// This function works by first checking the state_energy_map.  If this is empty it
  /// loads the data from disk.  If the nuclide is still not found, then the species
  /// is assumed to be in a ground state and 0.0 is returned.
  double state_energy(int nuc);
  /// Returns the excitation energy [MeV] of a \a nuc in a given state.
  double state_energy(char * nuc);
  /// Returns the excitation energy [MeV] of a \a nuc in a given state.
  double state_energy(std::string nuc);

  /// \brief Returns a set of decay children of a \a nuc.
  ///
  /// This function works by first checking decay_chidlren_map.  If this is empty it
  /// loads the data from disk.  If the nuclide is still not found, then the species
  /// is assumed to be stable and an empty set is returned.
  std::set<int> decay_children(int nuc);
  /// Returns the decay constant for a nuclide \a nuc.
  std::set<int> decay_children(char * nuc);
  /// Returns the decay constant for a nuclide \a nuc.
  std::set<int> decay_children(std::string nuc);

  /// Loads the level data from the nuc_data.h5 file into memory.
  void _load_level_data();

  /// \brief Returns the nuc_id of a metastable state
  ///
  /// This function looks through the level map for a given input nuc_id to find the
  /// nuc_id corresponding to the level
  int metastable_id(int nuc, int m);
  /// Assumes the first metastable state is the desired one
  int metastable_id(int nuc);

  /// a struct matching the '/decay/level_list' table in nuc_data.h5.
  typedef struct level_struct{
    int nuc_id;
    double level;
    double half_life;
    int metastable;
  } level_struct;

  /// Mapping from nuclides in id form to a struct containing data associated
  /// with that level.
  extern std::map<int, level_struct> level_data;

  template<typename T> void _load_data();

  /// a struct matching the '/decay/decays' table in nuc_data.h5.
  typedef struct decay_struct{
    int parent;
    int child;
    char decay[5];
    double half_life;
    double half_life_error;
    double branch_ratio;
    double photon_branch_ratio;
    double photon_branch_ratio_error;
    double beta_branch_ratio;
    double beta_branch_ratio_error;
  } decay_struct;

  /// Loads the decay data from the nuc_data.h5 file into memory.
  template<> void _load_data<decay_struct>();
  /// Mapping from a pair of nuclides in id form to a struct containing data
  /// associated with the decay from the first to the second
  extern std::map<std::pair<int, int>, decay_struct> decay_data;

  template<typename T, typename U> T data_access(std::pair<int, int> from_to, 
    size_t valoffset, std::map<std::pair<int, int>, U> &data);
  template<typename T, typename U> std::vector<T> data_access(int parent, 
    size_t valoffset, std::map<std::pair<int, int>, U> &data);
  //
  //void decay_data(std::pair<int, int> from_to, decay_struct *data);
  std::pair<double, double> decay_half_life(std::pair<int,int>);
  std::vector<std::pair<double, double> > decay_half_lifes(int);
  double decay_branch_ratio(std::pair<int,int>);
  std::vector<double> decay_branch_ratios(int parent);
  std::pair<double, double> decay_photon_branch_ratio(std::pair<int,int>);
  std::vector<std::pair<double, double> >decay_photon_branch_ratios(int parent);
  std::pair<double, double> decay_beta_branch_ratio(std::pair<int,int>);
  std::vector<std::pair<double, double> >decay_beta_branch_ratios(int parent);


  /// a struct matching the '/decay/gammas' table in nuc_data.h5.
  typedef struct gamma_struct{
    int from_nuc;
    int to_nuc;
    int parent_nuc;
    double energy;
    double energy_err;
    double photon_intensity;
    double photon_intensity_err;
    double conv_intensity;
    double conv_intensity_err;
    double total_intensity;
    double total_intensity_err;
    double k_conv_e;
    double l_conv_e;
    double m_conv_e;
  } gamma_struct;

  /// Loads the gamma ray data from the nuc_data.h5 file into memory.
  template<> void _load_data<gamma_struct>();

  class swapmapcompare{
    public:
        bool operator()(const std::pair<int, double>& lhs,
                        const std::pair<int, double>& rhs) const;
  };
  /// A vector of structs containing gamma ray data for access in memory
  //extern std::vector<gamma_struct> gamma_data;
  template<typename T, typename U> std::vector<T> data_access(double emin,
    double emax, size_t valoffset, std::map<std::pair<int, double>, U>  &data);
  template<typename T, typename U> std::vector<T> data_access(int parent,
    size_t valoffset, std::map<std::pair<int, double>, U>  &data);
  extern std::map<std::pair<int, double>, gamma_struct> gamma_data;

  //returns a list of gamma decay energies from input parent nuclide
  std::vector<std::pair<double, double> > gamma_energy(int parent);
  //returns a list of gamma photon intensities from input parent nuclide
  std::vector<std::pair<double, double> > gamma_photon_intensity(int parent);
  //returns a list of gamma conversion intensities from input parent nuclide
  std::vector<std::pair<double, double> > gamma_conversion_intensity(int parent);
  //returns a list of gamma total intensities from input parent nuclide
  std::vector<std::pair<double, double> > gamma_total_intensity(int parent);
  //returns a list of pairs of excited state transitions from an input parent nuclide
  std::vector<std::pair<int, int> > gamma_from_to(int parxent);
  //returns a list of pairs of excited state transitions from an decay energy
  std::vector<std::pair<int, int> > gamma_from_to(double energy, double error);
  //returns a list of parent nuclides associated with an input decay energy
  std::vector<int> gamma_parent(double energy, double error);

  /// a struct matching the '/decay/alphas' table in nuc_data.h5.
  typedef struct alpha_struct{
    int from_nuc;
    int to_nuc;
    double energy;
    double intensity;
  } alpha_struct;

  /// Loads the alpha decay data from the nuc_data.h5 file into memory.
  template<> void _load_data<alpha_struct>();

  /// A vector of structs containing alpha data for access in memory
  extern std::map<std::pair<int, double>, alpha_struct> alpha_data;
  
  //returns a list of alpha decay energies from input parent nuclide
  std::vector<double > alpha_energy(int parent);
  //returns a list of alpha decay intensities from input parent nuclide
  std::vector<double> alpha_intensity(int parent);
  //returns a list of alpha decay parents from input decay energy range
  std::vector<int> alpha_parent(double energy, double error);
  //returns a list of alpha decay children from input decay energy range
  std::vector<int> alpha_child(double energy, double error);
  //returns a list of alpha decay children from input parent nuclide
  std::vector<int> alpha_child(int parent);

  /// a struct matching the '/decay/betas' table in nuc_data.h5.
  typedef struct beta_struct{
    int from_nuc;
    int to_nuc;
    double endpoint_energy;
    double avg_energy;
    double intensity;
  } beta_struct;

  /// Loads the beta decay data from the nuc_data.h5 file into memory.
  template<> void _load_data<beta_struct>();

  /// A vector of structs containing beta data for access in memory
  extern std::map<std::pair<int, double>, beta_struct> beta_data;
  //returns a list of beta decay endpoint energies from input parent nuclide
  std::vector<double > beta_endpoint_energy(int parent);
  //returns a list of beta decay average energies from input parent nuclide
  std::vector<double > beta_average_energy(int parent);
  //returns a list of beta decay intensities from input parent nuclide
  std::vector<double> beta_intensity(int parent);
  //returns a list of beta decay parents from input decay energy range
  std::vector<int> beta_parent(double energy, double error);
  //returns a list of beta decay children from input decay energy range
  std::vector<int> beta_child(double energy, double error);
  //returns a list of beta decay children from input parent nuclide
  std::vector<int> beta_child(int parent);

  /// A struct matching the '/decay/ecbp' table in nuc_data.h5.
  typedef struct ecbp_struct{
    int from_nuc;
    int to_nuc;
    double endpoint_energy;
    double avg_energy;
    double beta_plus_intensity;
    double ec_intensity;
    double k_conv_e;
    double l_conv_e;
    double m_conv_e;
  } ecbp_struct;

  /// A vector of structs containing ecbp data for access in memory
  extern std::map<std::pair<int, double>, ecbp_struct> ecbp_data;

  /// Loads the electron capture and beta plus decay data from the
  /// nuc_data.h5 file into memory.
  template<> void _load_data<ecbp_struct>();
  ///returns a list of electron capture/ beta plus decay endpoint energies from
  ///input parent nuclide
  std::vector<double > ecbp_endpoint_energy(int parent);
  //returns a list of electron capture/ beta plus decay average energies from 
  //input parent nuclide
  std::vector<double > ecbp_average_energy(int parent);
  //returns a list of electron capture decay intensities from input parent 
  //nuclide
  std::vector<double> ec_intensity(int parent);
  //returns a list of beta plus decay intensities from input parent nuclide
  std::vector<double> bp_intensity(int parent);
  //returns a list of electron capture /beta plus decay parents from input 
  //decay energy range
  std::vector<int> ecbp_parent(double energy, double error);
  //returns a list of electron capture /beta plus decay children from input
  //decay energy range
  std::vector<int> ecbp_child(double energy, double error);
  //returns a list of electron capture /beta plus decay children from input 
  //parent nuclide
  std::vector<int> ecbp_child(int parent);
  /// \}
}

#endif
