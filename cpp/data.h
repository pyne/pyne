/// \file data.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Implements basic nuclear data functions.

#ifndef _PYNE_DATA_
#define _PYNE_DATA_
#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "h5wrap.h"

#include "extra_types.h"

#include "pyne.h"
#include "nucname.h"

namespace pyne
{
  extern std::string NUC_DATA_PATH; ///< Path to the nuc_data.h5 file.

  /// \name Atomic Mass Data
  /// \{

  /// Mapping from nuclides in id form to their atomic masses.
  extern std::map<int, double> atomic_mass_map;

  /// a struct matching the atomic_weight table in nuc_data.h5.
  typedef struct atomic_weight_struct {
    int nuc;      ///< nuclide in id form
    double mass;  ///< nuclide atomic mass [amu]
    double error; ///< error in atomic mass [amu]
    double abund; ///< natural abundance of nuclide [atom fraction]
  } atomic_weight_struct; 

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

  /// \name Scattering Length Data
  /// \{

  /// Mapping from nuclides in id form to their coherent scattering length.
  extern std::map<int, extra_types::complex_t> b_coherent_map;
  /// Mapping from nuclides in id form to their incoherent scattering length.
  extern std::map<int, extra_types::complex_t> b_incoherent_map;
  /// Mapping from nuclides in id form to their scattering length.
  extern std::map<int, double> b_map;

  /// a struct matching the '/neutron/scattering_lengths' table in nuc_data.h5.
  typedef struct scattering_lengths_struct {
    int nuc;  ///< nuclide in id form
    extra_types::complex_t b_coherent;  ///< coherent scattering length [cm]
    extra_types::complex_t b_incoherent;  ///< incoherent scattering length [cm]
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
  extra_types::complex_t b_coherent(int nuc);
  /// Finds the coherent scattering length [cm] for a nuclide \a nuc.
  extra_types::complex_t b_coherent(char * nuc);
  /// Finds the coherent scattering length [cm] for a nuclide \a nuc.
  extra_types::complex_t b_coherent(std::string nuc);

  /// \brief Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  ///
  /// This function works in the same way that b_coherent() does.
  extra_types::complex_t b_incoherent(int nuc);
  /// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  extra_types::complex_t b_incoherent(char * nuc);
  /// Finds the incoherent scattering length [cm] for a nuclide \a nuc.
  extra_types::complex_t b_incoherent(std::string nuc);

  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(int nuc);
  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(char * nuc);
  /// Computes the scattering length [cm] from the coherent and incoherent components.
  double b(std::string nuc);
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

  /// a struct matching the '/atomic_decay' table in nuc_data.h5.
  typedef struct atomic_decay_struct {
    int from_nuc; ///< parent species in id form
    double level; ///< decay level [MeV]
    int to_nuc;   ///< child species in id form
    double half_life;     ///< species half life [s]
    double decay_const;   ///< decay constant [1/s]
    double branch_ratio;  ///< decay branch ratio [fraction]
  } atomic_decay_struct;

  /// Loads the decay data from the nuc_data.h5 file into memory.
  void _load_atomic_decay();

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

  /// \}
}

#endif
