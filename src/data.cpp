// Implements basic nuclear data functions.
#ifndef PYNE_IS_AMALGAMATED
#include "data.h"
#include "atomic_data.h"
#endif

//
// Math Helpers
//

const double pyne::pi = 3.14159265359;
const double pyne::N_A = 6.0221415e+23;
const double pyne::barns_per_cm2 = 1e24;
const double pyne::cm2_per_barn = 1e-24;
const double pyne::sec_per_day = 24.0 * 3600.0;
const double pyne::MeV_per_K = 8.617343e-11;
const double pyne::MeV_per_MJ = 6.2415096471204E+18;
const double pyne::Bq_per_Ci = 3.7e10;
const double pyne::Ci_per_Bq = 2.7027027e-11;


/********************************/
/*** data_checksums Functions ***/
/********************************/

std::map<std::string, std::string> pyne::get_data_checksums() {
    std::map<std::string, std::string> temp_map;
    // Initialization of dataset hashes
    temp_map["/atomic_mass"]="10edfdc662e35bdfab91beb89285efff";
    temp_map["/material_library"]="8b10864378fbd88538434679acf908cc";
    temp_map["/neutron/eaf_xs"]="29622c636c4a3a46802207b934f9516c";
    temp_map["/neutron/scattering_lengths"]="a24d391cc9dc0fc146392740bb97ead4";
    temp_map["/neutron/simple_xs"]="3d6e086977783dcdf07e5c6b0c2416be";
    temp_map["/decay"]="4f41f3e46f4306cc44449f08a20922e0";
    temp_map["/dose_factors"]="dafa32c24b2303850a0bebdf3e6b122e";
    return temp_map;
}

std::map<std::string, std::string> pyne::data_checksums =
  pyne::get_data_checksums();

/*****************************/
/*** atomic_mass Functions ***/
/*****************************/
std::map<int, double> pyne::atomic_mass_map = std::map<int, double>();

void pyne::_load_atomic_mass_map() {
  // Loads the important parts of atomic_wight table into atomic_mass_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH)) {
    pyne::_load_atomic_mass_map_memory();
    return;
  }

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic_mass_data));
  H5Tinsert(desc, "nuc",   HOFFSET(atomic_mass_data, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "mass",  HOFFSET(atomic_mass_data, mass),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "error", HOFFSET(atomic_mass_data, error), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "abund", HOFFSET(atomic_mass_data, abund), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atomic_mass_set = H5Dopen2(nuc_data_h5, "/atomic_mass", H5P_DEFAULT);
  hid_t atomic_mass_space = H5Dget_space(atomic_mass_set);
  int atomic_mass_length = H5Sget_simple_extent_npoints(atomic_mass_space);

  // Read in the data
  atomic_mass_data * atomic_mass_array = new atomic_mass_data[atomic_mass_length];
  H5Dread(atomic_mass_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_mass_array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(atomic_mass_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of structs, put it in the map
  for(int n = 0; n < atomic_mass_length; n++) {
    atomic_mass_map.insert(std::pair<int, double>(atomic_mass_array[n].nuc, \
                                                  atomic_mass_array[n].mass));
    natural_abund_map.insert(std::pair<int, double>(atomic_mass_array[n].nuc, \
                                                    atomic_mass_array[n].abund));
  }

  delete[] atomic_mass_array;
}


double pyne::atomic_mass(int nuc) {
  // Find the nuclide's mass in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = atomic_mass_map.find(nuc);
  nuc_end = atomic_mass_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end) {
    return (*nuc_iter).second;
  }

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (atomic_mass_map.empty()) {
    // Don't fail if we can't load the library
    _load_atomic_mass_map();
    return atomic_mass(nuc);
  }

  double aw;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state mass...not strictly true, but good guess.
  if (0 < nucid%10000) {
    aw = atomic_mass((nucid/10000)*10000);
    if (atomic_mass_map.count(nuc) != 1) {
      atomic_mass_map.insert(std::pair<int, double>(nuc, aw));
    }
    return aw;
  };

  // Finally, if none of these work,
  // take a best guess based on the
  // aaa number.
  aw = (double) ((nucid/10000)%1000);
  if (atomic_mass_map.count(nuc) != 1) {
    atomic_mass_map.insert(std::pair<int, double>(nuc, aw));
  }
  return aw;
}


double pyne::atomic_mass(char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
}


double pyne::atomic_mass(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
}


/*******************************/
/*** natural_abund functions ***/
/*******************************/

std::map<int, double> pyne::natural_abund_map = std::map<int, double>();

double pyne::natural_abund(int nuc) {
  // Find the nuclide's natural abundance
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = natural_abund_map.find(nuc);
  nuc_end = natural_abund_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (natural_abund_map.empty()) {
    // Don't fail if we can't load the library
      _load_atomic_mass_map();
      return natural_abund(nuc);
  }

  double na;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state abundance...not strictly true, but good guess.
  if (0 < nucid%10000) {
    na = natural_abund((nucid/10000)*10000);
    natural_abund_map[nuc] = na;
    return na;
  }

  // Finally, if none of these work,
  // take a best guess based on the
  // aaa number.
  na = 0.0;
  natural_abund_map[nuc] = na;
  return na;
}


double pyne::natural_abund(char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
}


double pyne::natural_abund(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
}



/*************************/
/*** Q_value Functions ***/
/*************************/

void pyne::_load_q_val_map() {
  // Loads the important parts of q_value table into q_value_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(q_val_data));
  H5Tinsert(desc, "nuc", HOFFSET(q_val_data, nuc),  H5T_NATIVE_INT);
  H5Tinsert(desc, "q_val", HOFFSET(q_val_data, q_val), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "gamma_frac", HOFFSET(q_val_data, gamma_frac), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t q_val_set = H5Dopen2(nuc_data_h5, "/decay/q_values", H5P_DEFAULT);
  hid_t q_val_space = H5Dget_space(q_val_set);
  int q_val_length = H5Sget_simple_extent_npoints(q_val_space);

  // Read in the data
  q_val_data * q_val_array = new q_val_data[q_val_length];
  H5Dread(q_val_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, q_val_array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(q_val_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of structs, put it in the map
  for(int n = 0; n < q_val_length; n++) {
    q_val_map[q_val_array[n].nuc] = q_val_array[n].q_val;
    gamma_frac_map[q_val_array[n].nuc] = q_val_array[n].gamma_frac;
  }

  delete[] q_val_array;
}

std::map<int, double> pyne::q_val_map = std::map<int, double>();

double pyne::q_val(int nuc) {
  // Find the nuclide's q_val in MeV/fission
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = q_val_map.find(nuc);
  nuc_end = q_val_map.end();

  // First check if we already have the nuc q_val in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the nuc_data.h5 if the map is empty.
  if (q_val_map.empty()) {
      _load_q_val_map();
      return q_val(nuc);
  };

  double qv;
  int nucid = nucname::id(nuc);
  if (nucid != nuc)
    return q_val(nucid);

  // If nuclide is not found, return 0
  qv = 0.0;
  q_val_map[nuc] = qv;
  return qv;
}


double pyne::q_val(const char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return q_val(nuc_zz);
}


double pyne::q_val(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return q_val(nuc_zz);
}


/****************************/
/*** gamma_frac functions ***/
/****************************/

std::map<int, double> pyne::gamma_frac_map = std::map<int, double>();

double pyne::gamma_frac(int nuc) {
  // Find the nuclide's fraction of Q that comes from gammas
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = gamma_frac_map.find(nuc);
  nuc_end = gamma_frac_map.end();

  // First check if we already have the gamma_frac in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from nuc_data.h5 if the map is empty.
  if (gamma_frac_map.empty()) {
      _load_q_val_map();
      return gamma_frac(nuc);
  }

  double gf;
  int nucid = nucname::id(nuc);
  if (nucid != nuc)
    return gamma_frac(nucid);

  // If nuclide is not found, return 0
  gf = 0.0;
  gamma_frac_map[nucid] = gf;
  return gf;
}


double pyne::gamma_frac(const char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return gamma_frac(nuc_zz);
}


double pyne::gamma_frac(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return gamma_frac(nuc_zz);
}


/*****************************/
/*** Dose Factor Functions ***/
/*****************************/
/***************************************************************************
This data is from: [Exposure Scenarios and Unit Dose Factors for the Hanford
Immobilized Low-Activity Tank Waste Performance Assessment, ref.
HNF-SD-WM-TI-707 Rev. 1 December 1999] Appendix O of HNF-5636 [DATA PACKAGES
FOR THE HANFORD IMMOBILIZED LOW-ACTIVITY TANK WASTE PERFORMANCE ASSESSMENT:
2001 VERSION]

Liability Disclaimer: The PyNE Development Team shall not be liable for any
loss or injury resulting from decisions made with this data.
**************************************************************************/

void pyne::_load_dose_map(std::map<int, dose>& dm, std::string source_path) {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Defining string type for lung model data
  hid_t string_type_;
  string_type_ = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_, 1);
  H5Tset_strpad(string_type_, H5T_STR_NULLPAD);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(dose));
  status = H5Tinsert(desc, "nuc", HOFFSET(dose, nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "ext_air_dose", HOFFSET(dose, ext_air_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ratio", HOFFSET(dose, ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ext_soil_dose", HOFFSET(dose, ext_soil_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ingest_dose", HOFFSET(dose, ingest_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "fluid_frac", HOFFSET(dose, fluid_frac), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "inhale_dose", HOFFSET(dose, inhale_dose), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "lung_mod", HOFFSET(dose, lung_mod), string_type_);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Convert source_path to proper format for HD5open
  const char * c = source_path.c_str();

  // Open the data set
  hid_t dose_set = H5Dopen2(nuc_data_h5, c, H5P_DEFAULT);
  hid_t dose_space = H5Dget_space(dose_set);
  int dose_length = H5Sget_simple_extent_npoints(dose_space);

  // Read in the data
  dose * dose_array = new dose[dose_length];
  H5Dread(dose_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, dose_array);

  // Put array of structs in the map
  for (int n = 0; n < dose_length; n++) {
    dm[dose_array[n].nuc] = dose_array[n];
  }

  // Close the nuc_data library
  H5Dclose(dose_set);
  H5Tclose(string_type_);
  H5Fclose(nuc_data_h5);

  delete[] dose_array;
}

///
/// Functions for Source Location in nuc_data.h5
/// and related Map Pointers
///

std::string source_string(int source) {
  std::string source_location;
  if (source == 1) {
    source_location = "/dose_factors/DOE";
  } else if (source == 2) {
    source_location = "/dose_factors/GENII";
  } else {
    source_location = "/dose_factors/EPA";
  }
  return source_location;
}

std::map<int, pyne::dose>& dose_source_map(int source) {
  std::map<int, pyne::dose>* dm;
  if (source == 1) {
    dm = &pyne::doe_dose_map;
  } else if (source == 2) {
    dm = &pyne::genii_dose_map;
  } else {
    dm = &pyne::epa_dose_map;
  }
  if (dm->empty()) {
      std::string source_path = source_string(source);
      _load_dose_map(*dm, source_path);
  }
  return *dm;
}

std::map<int, pyne::dose> pyne::epa_dose_map;
std::map<int, pyne::dose> pyne::doe_dose_map;
std::map<int, pyne::dose> pyne::genii_dose_map;

///
/// Functions for External Air and
/// Ratio of External Air to Inhalation Dose Factors
///

/// External Air
double pyne::ext_air_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].ext_air_dose;
  } else {
    return -1;
  }
}

double pyne::ext_air_dose(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_air_dose(nuc_zz, source);
}


double pyne::ext_air_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_air_dose(nuc_zz, source);
}

/// Dose Ratio
double pyne::dose_ratio(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].ratio;
  } else {
    return -1;
  }
}

double pyne::dose_ratio(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_ratio(nuc_zz, source);
}


double pyne::dose_ratio(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_ratio(nuc_zz, source);
}

///
/// Function for External Soil Dose Factors
///

double pyne::ext_soil_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].ext_soil_dose;
  } else {
    return -1;
  }
}

double pyne::ext_soil_dose(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_soil_dose(nuc_zz, source);
}


double pyne::ext_soil_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ext_soil_dose(nuc_zz, source);
}

///
/// Functions for Ingestion Dose Factors and
/// Fraction of activity that is absorbed by body fluids
///

/// Ingestion
double pyne::ingest_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].ingest_dose;
  } else {
    return -1;
  }
}

double pyne::ingest_dose(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ingest_dose(nuc_zz, source);
}

double pyne::ingest_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return ingest_dose(nuc_zz, source);
}

/// Fluid Fraction
double pyne::dose_fluid_frac(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].fluid_frac;
  } else {
    return -1;
  }
}

double pyne::dose_fluid_frac(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_fluid_frac(nuc_zz, source);
}

double pyne::dose_fluid_frac(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_fluid_frac(nuc_zz, source);
}

///
/// Functions for Inhalation Dose Factors and
/// Lung Model used to obtain dose factors
///

/// Inhalation
double pyne::inhale_dose(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return dm[nucid].inhale_dose;
  } else {
    return -1;
  }
}

double pyne::inhale_dose(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return inhale_dose(nuc_zz, source);
}

double pyne::inhale_dose(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return inhale_dose(nuc_zz, source);
}

/// Lung Model
std::string pyne::dose_lung_model(int nuc, int source) {
  std::map<int, pyne::dose>& dm = dose_source_map(source);
  int nucid = nucname::id(nuc);

  if (dm.count(nucid)==1) {
    return std::string(1, dm[nucid].lung_mod);
  } else {
    return "Nada";
  }
}

std::string pyne::dose_lung_model(const char * nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_lung_model(nuc_zz, source);
}


std::string pyne::dose_lung_model(std::string nuc, int source) {
  int nuc_zz = nucname::id(nuc);
  return dose_lung_model(nuc_zz, source);
}


/***********************************/
/*** scattering length functions ***/
/***********************************/
std::map<int, xd_complex_t> pyne::b_coherent_map = std::map<int, xd_complex_t>();
std::map<int, xd_complex_t> pyne::b_incoherent_map = std::map<int, xd_complex_t>();
std::map<int, double> pyne::b_map = std::map<int, double>();


void pyne::_load_scattering_lengths() {
  // Loads the important parts of atomic_wight table into atomic_mass_map
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(scattering_lengths));
  status = H5Tinsert(desc, "nuc", HOFFSET(scattering_lengths, nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "b_coherent", HOFFSET(scattering_lengths, b_coherent),
                      h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "b_incoherent", HOFFSET(scattering_lengths, b_incoherent),
                      h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "xs_coherent", HOFFSET(scattering_lengths, xs_coherent),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs_incoherent", HOFFSET(scattering_lengths, xs_incoherent),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs", HOFFSET(scattering_lengths, xs), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t scat_len_set = H5Dopen2(nuc_data_h5, "/neutron/scattering_lengths", H5P_DEFAULT);
  hid_t scat_len_space = H5Dget_space(scat_len_set);
  int scat_len_length = H5Sget_simple_extent_npoints(scat_len_space);

  // Read in the data
  scattering_lengths * scat_len_array = new scattering_lengths[scat_len_length];
  status = H5Dread(scat_len_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, scat_len_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(scat_len_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for(int n = 0; n < scat_len_length; n++) {
    b_coherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_coherent;
    b_incoherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_incoherent;
  }

  delete[] scat_len_array;
}



//
// Coherent functions
//


xd_complex_t pyne::b_coherent(int nuc) {
  // Find the nuclide's bound scattering length in cm
  std::map<int, xd_complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_coherent_map.find(nuc);
  nuc_end = b_coherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (b_coherent_map.empty()) {
    _load_scattering_lengths();
    return b_coherent(nuc);
  }

  xd_complex_t bc;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (anum == nucname::anum((*nuc_iter).first)) {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    }
    nuc_iter++;
  }

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (znum == nucname::znum((*nuc_iter).first)) {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    }
    nuc_iter++;
  }

  // Finally, if none of these work,
  // just return zero...
  bc.re = 0.0;
  bc.im = 0.0;
  b_coherent_map[nuc] = bc;
  return bc;
}


xd_complex_t pyne::b_coherent(char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
}


xd_complex_t pyne::b_coherent(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
}



//
// Incoherent functions
//


xd_complex_t pyne::b_incoherent(int nuc) {
  // Find the nuclide's bound inchoherent scattering length in cm
  std::map<int, xd_complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_incoherent_map.find(nuc);
  nuc_end = b_incoherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (b_incoherent_map.empty()) {
    _load_scattering_lengths();
    return b_incoherent(nuc);
  }

  xd_complex_t bi;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (anum == nucname::anum((*nuc_iter).first)) {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    }
    nuc_iter++;
  }

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end) {
    if (znum == nucname::znum((*nuc_iter).first)) {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    }
    nuc_iter++;
  }

  // Finally, if none of these work,
  // just return zero...
  bi.re = 0.0;
  bi.im = 0.0;
  b_incoherent_map[nuc] = bi;
  return bi;
}


xd_complex_t pyne::b_incoherent(char * nuc) {
  return b_incoherent(nucname::id(nuc));
}


xd_complex_t pyne::b_incoherent(std::string nuc) {
  return b_incoherent(nucname::id(nuc));
}



//
// b functions
//

double pyne::b(int nuc) {
  // Find the nuclide's bound scattering length in cm
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = b_map.find(nuc);
  nuc_end = b_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, calculate the value from coherent and incoherent lengths
  xd_complex_t bc = b_coherent(nuc);
  xd_complex_t bi = b_incoherent(nuc);

  double b_val = sqrt(bc.re*bc.re + bc.im*bc.im + bi.re*bi.re + bi.im*bi.im);

  return b_val;
}


double pyne::b(char * nuc) {
  int nucid = nucname::id(nuc);
  return b(nucid);
}


double pyne::b(std::string nuc) {
  int nucid = nucname::id(nuc);
  return b(nucid);
}



//
// Fission Product Yield Data
//
std::map<std::pair<int, int>, double> pyne::wimsdfpy_data = \
  std::map<std::pair<int, int>, double>();

void pyne::_load_wimsdfpy() {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(wimsdfpy));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(wimsdfpy, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(wimsdfpy, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yields", HOFFSET(wimsdfpy, yields),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t wimsdfpy_set = H5Dopen2(nuc_data_h5, "/neutron/wimsd_fission_products",
                                H5P_DEFAULT);
  hid_t wimsdfpy_space = H5Dget_space(wimsdfpy_set);
  int wimsdfpy_length = H5Sget_simple_extent_npoints(wimsdfpy_space);

  // Read in the data
  wimsdfpy * wimsdfpy_array = new wimsdfpy[wimsdfpy_length];
  status = H5Dread(wimsdfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, wimsdfpy_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(wimsdfpy_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for(int n=0; n < wimsdfpy_length; n++) {
    wimsdfpy_data[std::make_pair(wimsdfpy_array[n].from_nuc,
      wimsdfpy_array[n].to_nuc)] = wimsdfpy_array[n].yields;
  }

  delete[] wimsdfpy_array;
}


std::map<std::pair<int, int>, pyne::ndsfpysub> pyne::ndsfpy_data = \
  std::map<std::pair<int, int>, pyne::ndsfpysub>();

void pyne::_load_ndsfpy() {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ndsfpy));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ndsfpy, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ndsfpy, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yield_thermal", HOFFSET(ndsfpy, yield_thermal),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_thermal_err", HOFFSET(ndsfpy, yield_thermal_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast", HOFFSET(ndsfpy, yield_fast),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast_err", HOFFSET(ndsfpy, yield_fast_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV", HOFFSET(ndsfpy, yield_14MeV),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV_err", HOFFSET(ndsfpy, yield_14MeV_err),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t ndsfpy_set = H5Dopen2(nuc_data_h5, "/neutron/nds_fission_products",
                                H5P_DEFAULT);
  hid_t ndsfpy_space = H5Dget_space(ndsfpy_set);
  int ndsfpy_length = H5Sget_simple_extent_npoints(ndsfpy_space);

  // Read in the data
  ndsfpy * ndsfpy_array = new ndsfpy[ndsfpy_length];
  status = H5Dread(ndsfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ndsfpy_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(ndsfpy_set);
  status = H5Fclose(nuc_data_h5);

  ndsfpysub ndsfpysub_temp;

  // Ok now that we have the array of structs, put it in the maps
  for(int n=0; n < ndsfpy_length; n++) {
    ndsfpysub_temp.yield_thermal = ndsfpy_array[n].yield_thermal;
    ndsfpysub_temp.yield_thermal_err = ndsfpy_array[n].yield_thermal_err;
    ndsfpysub_temp.yield_fast = ndsfpy_array[n].yield_fast;
    ndsfpysub_temp.yield_fast_err = ndsfpy_array[n].yield_fast_err;
    ndsfpysub_temp.yield_14MeV = ndsfpy_array[n].yield_14MeV;
    ndsfpysub_temp.yield_14MeV_err = ndsfpy_array[n].yield_14MeV_err;
    ndsfpy_data[std::make_pair(ndsfpy_array[n].from_nuc,
      ndsfpy_array[n].to_nuc)] = ndsfpysub_temp;
  }



  delete[] ndsfpy_array;
}

double pyne::fpyield(std::pair<int, int> from_to, int source, bool get_error) {
  // Note that this may be expanded eventually to include other
  // sources of fission product data.

  // Find the parent/child pair branch ratio as a fraction
  if (source == 0) {
    std::map<std::pair<int, int>, double>::iterator fpy_iter, fpy_end;
    fpy_iter = wimsdfpy_data.find(from_to);
    fpy_end = wimsdfpy_data.end();
    if (fpy_iter != fpy_end)
        //if (get_error == true) return 0;
        return (*fpy_iter).second;
  } else {
    std::map<std::pair<int, int>, ndsfpysub>::iterator fpy_iter, fpy_end;
    fpy_iter = ndsfpy_data.find(from_to);
    fpy_end = ndsfpy_data.end();
    if (fpy_iter != fpy_end) {
        switch (source) {
          case 1:
            if (get_error)
                return (*fpy_iter).second.yield_thermal_err;
            return (*fpy_iter).second.yield_thermal;
            break;
          case 2:
            if (get_error)
                return (*fpy_iter).second.yield_fast_err;
            return (*fpy_iter).second.yield_fast;
            break;
          case 3:
            if (get_error)
                return (*fpy_iter).second.yield_14MeV_err;
            return (*fpy_iter).second.yield_14MeV;
            break;
        }
    }
  }


  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if ((source == 0 ) && (wimsdfpy_data.empty())) {
    _load_wimsdfpy();
    return fpyield(from_to, 0, get_error);
  }else if (ndsfpy_data.empty()) {
    _load_ndsfpy();
    return fpyield(from_to, source, get_error);
  }

  // Finally, if none of these work,
  // assume the value is stable
  double fpy = 0.0;
  wimsdfpy_data[from_to] = fpy;
  return fpy;
}

double pyne::fpyield(int from_nuc, int to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}

double pyne::fpyield(char * from_nuc, char * to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}

double pyne::fpyield(std::string from_nuc, std::string to_nuc, int source,
                     bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc),
                                     nucname::id(to_nuc)), source, get_error);
}


/***********************/
/*** decay functions ***/
/***********************/

//
// Data access tools
//

bool pyne::swapmapcompare::operator()(const std::pair<int, double>& lhs,
const std::pair<int, double>& rhs) const {
    return lhs.second<rhs.second || (!(rhs.second<lhs.second) &&
      lhs.first<rhs.first);
}

template<typename T, typename U> std::vector<T> pyne::data_access(
double energy_min, double energy_max, size_t valoffset, std::map<std::pair<int,
double>, U>  &data) {
  typename std::map<std::pair<int, double>, U, swapmapcompare>::iterator
    nuc_iter, nuc_end, it;
  std::map<std::pair<int, double>, U, swapmapcompare> dc(data.begin(),
    data.end());
  std::vector<T> result;
  if (energy_max < energy_min){
    double temp = energy_max;
    energy_max = energy_min;
    energy_min = temp;
  }
  nuc_iter = dc.lower_bound(std::make_pair(0, energy_min));
  nuc_end = dc.upper_bound(std::make_pair(9999999999, energy_max));
  T *ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it!= nuc_end; ++it){
    ret = (T *)((char *)&(it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<T, U>(energy_min, energy_max, valoffset, data);
  }
  return result;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
double min, double max, size_t valoffset,
std::map<std::pair<int, double>, U>  &data) {
  typename std::map<std::pair<int, double>, U>::iterator nuc_iter, nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent,min));
  nuc_end = data.upper_bound(std::make_pair(parent,max));
  T *ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it!= nuc_end; ++it){
    ret = (T *)((char *)&(it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<T, U>(parent, min, max, valoffset, data);
  }
  return result;
}

template<typename T, typename U> T pyne::data_access(std::pair<int, int>
from_to, size_t valoffset, std::map<std::pair<int, int>, U> &data) {
  typename std::map<std::pair<int, int>, U>::iterator nuc_iter, nuc_end;

  nuc_iter = data.find(from_to);
  nuc_end = data.end();
  T *ret;
  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end){
    ret = (T *)((char *)&(nuc_iter->second) + valoffset);
    return *ret;
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<T, U>(from_to, valoffset, data);
  }
  // This is okay for now because we only return ints and doubles
  return 0;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
size_t valoffset, std::map<std::pair<int, int>, U> &data){
  typename std::map<std::pair<int, int>, U>::iterator nuc_iter, nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent,0));
  nuc_end = data.upper_bound(std::make_pair(parent,9999999999));
  T *ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it!= nuc_end; ++it){
    ret = (T *)((char *)&(it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<T, U>(parent, valoffset, data);
  }
  return result;
}

template<typename T, typename U> std::vector<T> pyne::data_access(int parent,
size_t valoffset, std::map<std::pair<int, unsigned int>, U> &data){
  typename std::map<std::pair<int, unsigned int>, U>::iterator nuc_iter,
   nuc_end, it;
  std::vector<T> result;
  nuc_iter = data.lower_bound(std::make_pair(parent,0));
  nuc_end = data.upper_bound(std::make_pair(parent,UINT_MAX));
  T *ret;
  // First check if we already have the nuc in the map
  for (it = nuc_iter; it!= nuc_end; ++it){
    ret = (T *)((char *)&(it->second) + valoffset);
    result.push_back(*ret);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<T, U>(parent, valoffset, data);
  }
  return result;
}

template<typename U> double pyne::data_access(int nuc,
size_t valoffset, std::map<int, U> &data){
  typename std::map<int, U>::iterator nuc_iter,
   nuc_end;
  nuc_iter = data.find(nuc);
  nuc_end = data.end();
  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end){
    return *(double *)((char *)&(nuc_iter->second) + valoffset);
  }
  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (data.empty())
  {
    _load_data<U>();
    return data_access<U>(nuc, valoffset, data);
  }
  throw pyne::nucname::NotANuclide(nuc, "");
}


//
// Load atomic data
//

std::map<int, pyne::atomic> pyne::atomic_data_map;

template<> void pyne::_load_data<pyne::atomic>() {
  // Loads the atomic table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic));
  status = H5Tinsert(desc, "z", HOFFSET(atomic, z),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "k_shell_fluor", HOFFSET(atomic, k_shell_fluor),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_fluor_error", HOFFSET(atomic, k_shell_fluor_error),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_shell_fluor", HOFFSET(atomic, l_shell_fluor),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_shell_fluor_error", HOFFSET(atomic, l_shell_fluor_error),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "prob", HOFFSET(atomic, prob),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_be", HOFFSET(atomic, k_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_shell_be_err", HOFFSET(atomic, k_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "li_shell_be", HOFFSET(atomic, li_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "li_shell_be_err", HOFFSET(atomic, li_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "mi_shell_be", HOFFSET(atomic, mi_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "mi_shell_be_err", HOFFSET(atomic, mi_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ni_shell_be", HOFFSET(atomic, ni_shell_be),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ni_shell_be_err", HOFFSET(atomic, ni_shell_be_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_to_ka", HOFFSET(atomic, kb_to_ka),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_to_ka_err", HOFFSET(atomic, kb_to_ka_err),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_to_ka1", HOFFSET(atomic, ka2_to_ka1),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_to_ka1_err", HOFFSET(atomic, ka2_to_ka1_err),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_auger", HOFFSET(atomic, l_auger),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_auger", HOFFSET(atomic, k_auger),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka1_x_ray_en", HOFFSET(atomic, ka1_x_ray_en),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka1_x_ray_en_err", HOFFSET(atomic, ka1_x_ray_en_err),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_x_ray_en", HOFFSET(atomic, ka2_x_ray_en),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ka2_x_ray_en_err", HOFFSET(atomic, ka2_x_ray_en_err),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "kb_x_ray_en", HOFFSET(atomic, kb_x_ray_en),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_x_ray_en", HOFFSET(atomic, l_x_ray_en),
                      H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);
  // Open the data set
  hid_t atomic_set = H5Dopen2(nuc_data_h5, "/decay/atomic", H5P_DEFAULT);
  hid_t atomic_space = H5Dget_space(atomic_set);
  int atomic_length = H5Sget_simple_extent_npoints(atomic_space);

  // Read in the data
  atomic * atomic_array = new atomic[atomic_length];
  status = H5Dread(atomic_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   atomic_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(atomic_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < atomic_length; ++i) {
      atomic_data_map[atomic_array[i].z] = atomic_array[i];
  }

  delete[] atomic_array;

}

std::vector<std::pair<double, double> >
  pyne::calculate_xray_data(int z, double k_conv, double l_conv) {
  double xk = 0;
  double xka = 0;
  double xka1 = 0;
  double xka2 = 0;
  double xkb = 0;
  double xl = 0;
  if (!isnan(k_conv)) {
    xk = data_access<atomic> (z, offsetof(atomic, k_shell_fluor),
     atomic_data_map)*k_conv;
    xka = xk / (1.0 + data_access<atomic> (z, offsetof(atomic,
     kb_to_ka), atomic_data_map));
    xka1 = xka / (1.0 + data_access<atomic> (z, offsetof(atomic,
     ka2_to_ka1), atomic_data_map));
    xka2 = xka - xka1;
    xkb = xk - xka;
    if (!isnan(l_conv)) {
        xl = (l_conv + k_conv*data_access<atomic> (z, offsetof(atomic,
     prob), atomic_data_map))*data_access<atomic> (z, offsetof(atomic,
     l_shell_fluor), atomic_data_map);
    }
  } else if (!isnan(l_conv)) {
    xl = l_conv*data_access<atomic> (z, offsetof(atomic,
     l_shell_fluor), atomic_data_map);
  }
  std::vector<std::pair<double, double> > result;
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
     ka1_x_ray_en), atomic_data_map),xka1));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
     ka2_x_ray_en), atomic_data_map),xka2));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
     kb_x_ray_en), atomic_data_map),xkb));
  result.push_back(std::make_pair(data_access<atomic> (z, offsetof(atomic,
     l_x_ray_en), atomic_data_map),xl));

  return result;
}


//
// Load level data
//

std::map<std::pair<int,double>, pyne::level_data> pyne::level_data_lvl_map;
std::map<std::pair<int,unsigned int>,
  pyne::level_data> pyne::level_data_rx_map;


template<> void pyne::_load_data<pyne::level_data>()
{

  // Loads the level table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(level_data));
  status = H5Tinsert(desc, "nuc_id", HOFFSET(level_data, nuc_id),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "rx_id", HOFFSET(level_data, rx_id),
                     H5T_NATIVE_UINT);
  status = H5Tinsert(desc, "half_life", HOFFSET(level_data, half_life),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "level", HOFFSET(level_data, level),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(level_data, branch_ratio),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "metastable", HOFFSET(level_data, metastable),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "special", HOFFSET(level_data, special),
                      H5T_C_S1);
  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);
  // Open the data set
  hid_t level_set = H5Dopen2(nuc_data_h5, "/decay/level_list", H5P_DEFAULT);
  hid_t level_space = H5Dget_space(level_set);
  int level_length = H5Sget_simple_extent_npoints(level_space);

  // Read in the data
  level_data * level_array = new level_data[level_length];
  status = H5Dread(level_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   level_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(level_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < level_length; ++i) {
    if (level_array[i].rx_id == 0)
      level_data_lvl_map[std::make_pair(level_array[i].nuc_id,
                                        level_array[i].level)] = level_array[i];
    else
      level_data_rx_map[std::make_pair(level_array[i].nuc_id,
                                       level_array[i].rx_id)] = level_array[i];
  }

  delete[] level_array;
}

//
// level id
//
int pyne::id_from_level(int nuc, double level, std::string special) {
  int nostate = (nuc / 10000) * 10000;
  if (level_data_lvl_map.empty()) {
    _load_data<level_data>();
  }

  std::map<std::pair<int, double>, level_data>::iterator nuc_lower, nuc_upper;

  nuc_lower = level_data_lvl_map.lower_bound(std::make_pair(nostate, 0.0));
  nuc_upper = level_data_lvl_map.upper_bound(std::make_pair(nostate+9999,
                                             DBL_MAX));
  double minv = DBL_MAX;
  //by default return input nuc_id with level stripped
  int ret_id = nuc;
  for (std::map<std::pair<int, double>, level_data>::iterator it=nuc_lower;
       it!=nuc_upper; ++it) {
    if ((std::abs(level - it->second.level) < minv) &&
    ((char)it->second.special == special.c_str()[0]) &&
    !isnan(it->second.level)) {
      minv = std::abs(level - it->second.level);
      ret_id = it->second.nuc_id;
    }
  }
  // This value was chosen so important transitions in U-235 are not missed
  if (minv > 3.0)
    return -nuc;
  else
    return ret_id;
}

int pyne::id_from_level(int nuc, double level){
    return id_from_level(nuc, level, " ");
}
//
// Metastable id data
//

int pyne::metastable_id(int nuc, int m) {
  int nostate = (nuc / 10000) * 10000;
  if (m==0) return nostate;
  if (level_data_lvl_map.empty()) {
    _load_data<level_data>();
  }

  std::map<std::pair<int, double>, level_data>::iterator nuc_lower, nuc_upper;

  nuc_lower = level_data_lvl_map.lower_bound(std::make_pair(nostate, 0.0));
  nuc_upper = level_data_lvl_map.upper_bound(std::make_pair(nostate+9999,
                                                            DBL_MAX));
  for (std::map<std::pair<int, double>, level_data>::iterator it=nuc_lower;
  it!=nuc_upper; ++it) {
    if (it->second.metastable == m)
        return it->second.nuc_id;
  }

  return -1;
}

int pyne::metastable_id(int nuc) {
  return metastable_id(nuc, 1);
}

//
// Decay children data
//


std::set<int> pyne::decay_children(int nuc) {
  // make sure spontaneous fission data is loaded
  if (wimsdfpy_data.empty())
    _load_wimsdfpy();

  std::vector<unsigned int> part = data_access<unsigned int, level_data>(nuc,
    offsetof(level_data, rx_id), level_data_rx_map);
  std::vector<unsigned int>::iterator it = part.begin();
  std::set<int> result;
  for (; it != part.end(); ++it) {
    switch (*it) {
      case 36125: {
        // internal transition, rx == 'it'
        result.insert((nuc /10000) * 10000);
        break;
      }
      case 36565:
      case 1794828612:
      {
        // spontaneous fission, rx == 'sf', 36565
        // beta- & spontaneous fission, rx == 'b-sf', 1794828612
        std::map<std::pair<int, int>, double>::iterator sf = wimsdfpy_data.begin();
        for (; sf != wimsdfpy_data.end(); ++sf)
          if (sf->first.first == nuc)
            result.insert(sf->first.second);
        break;
      }
      default: {
        result.insert((rxname::child(nuc, *it, "decay") /10000) * 10000);
      }
    }
  }
  return result;
}

std::set<int> pyne::decay_children(char * nuc)
{
  return decay_children(nucname::id(nuc));
}

std::set<int> pyne::decay_children(std::string nuc)
{
  return decay_children(nucname::id(nuc));
}

//
// Excitation state energy data
//

double pyne::state_energy(int nuc)
{
  std::vector<double> result = data_access<double, level_data>(nuc, 0.0,
  DBL_MAX, offsetof(level_data, level), level_data_lvl_map);
  if (result.size() == 1)
    return result[0]/1000.0;
  return 0.0;
}

double pyne::state_energy(char * nuc)
{
  return state_energy(nucname::id(nuc));
}


double pyne::state_energy(std::string nuc)
{
  return state_energy(nucname::id(nuc));
}


//
// Decay constant data
//

double pyne::decay_const(int nuc)
{
    std::vector<double> result = data_access<double, level_data>(nuc, 0.0,
      DBL_MAX, offsetof(level_data, half_life), level_data_lvl_map);
    if (result.size() == 1) {
        return log(2.0)/result[0];
    }
    return 0.0;
}


double pyne::decay_const(char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
}


double pyne::decay_const(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
}


//
// Half-life data
//
double pyne::half_life(int nuc) {
  std::vector<double> result = data_access<double, level_data>(
      nuc, 0.0, DBL_MAX, offsetof(level_data, half_life), level_data_lvl_map);
  if (result.size() == 1) {
    return result[0];
  }
  double num = 1.0;
  double den = 0.0;
  return num / den;
}

double pyne::half_life(char * nuc) {
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
}

double pyne::half_life(std::string nuc) {
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
}


//
// Branch ratio data
//
double pyne::branch_ratio(std::pair<int, int> from_to) {
  using std::vector;
  using pyne::nucname::groundstate;
  // make sure spontaneous fission data is loaded
  if (wimsdfpy_data.empty())
    _load_wimsdfpy();

  vector<unsigned int> part1 = \
    data_access<unsigned int, level_data>(from_to.first, offsetof(level_data, rx_id),
                                          level_data_rx_map);
  vector<double> part2 = \
    data_access<double, level_data>(from_to.first, offsetof(level_data, branch_ratio),
                                    level_data_rx_map);
  double result = 0.0;
  if ((from_to.first == from_to.second) && (decay_const(from_to.first) == 0.0))
    return 1.0;
  for (vector<unsigned int>::size_type i = 0; i < part1.size(); ++i) {
    if ((part1[i] == 36125) &&
        (groundstate(from_to.first) == groundstate(from_to.second)) &&
        (from_to.second % 10000 == 0)) {
      // internal transition, rx == 'it', 36125
      result = part2[i] * 0.01;
      break;
    } else if (part1[i] == 36565 || part1[i] == 1794828612) {
      // spontaneous fission, rx == 'sf', 36565
      // beta- & spontaneous fission, rx == 'b-sf', 1794828612
      result += part2[i] * 0.01 * wimsdfpy_data[from_to];
    } else if ((part1[i] != 0) && (groundstate(rxname::child(from_to.first,
                                   part1[i], "decay")) == from_to.second)) {
      result += part2[i] * 0.01;
    }
  }
  return result;
}

double pyne::branch_ratio(int from_nuc, int to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

double pyne::branch_ratio(char * from_nuc, char * to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

double pyne::branch_ratio(std::string from_nuc, std::string to_nuc) {
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc),
                                          nucname::id(to_nuc)));
}

std::map<std::pair<int, int>, pyne::decay> pyne::decay_data = \
  std::map<std::pair<int, int>, pyne::decay>();

template<> void pyne::_load_data<pyne::decay>() {

  // Loads the decay table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(decay));
  status = H5Tinsert(desc, "parent", HOFFSET(decay, parent),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "child", HOFFSET(decay, child),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "decay", HOFFSET(decay, decay),
                     H5T_NATIVE_UINT);
  status = H5Tinsert(desc, "half_life", HOFFSET(decay, half_life),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "half_life_error", HOFFSET(decay,
                     half_life_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(decay, branch_ratio),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio_error", HOFFSET(decay, branch_ratio_error),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio", HOFFSET(decay,
                     photon_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio_err", HOFFSET(decay,
                     photon_branch_ratio_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio", HOFFSET(decay,
                     beta_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio_err", HOFFSET(decay,
                     beta_branch_ratio_error), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t decay_set = H5Dopen2(nuc_data_h5, "/decay/decays", H5P_DEFAULT);
  hid_t decay_space = H5Dget_space(decay_set);
  int decay_length = H5Sget_simple_extent_npoints(decay_space);

  // Read in the data
  decay * decay_array = new decay[decay_length];
  status = H5Dread(decay_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   decay_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(decay_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < decay_length; ++i) {
    decay_data[std::make_pair(decay_array[i].parent, decay_array[i].child)] = \
      decay_array[i];
  }
  delete[] decay_array;
}


std::vector<int> pyne::decay_data_children(int parent) {
  std::vector<int> result = data_access<int, decay>(parent,
    offsetof(decay, child), decay_data);
  return result;
}

std::pair<double, double> pyne::decay_half_life(std::pair<int, int> from_to){
  return std::make_pair(data_access<double, decay>(from_to, offsetof(
   decay, half_life), decay_data), data_access<double, decay>(
   from_to, offsetof(decay, half_life_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_half_lifes(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
    offsetof(decay, half_life), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
    offsetof(decay, half_life_error), decay_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::pair<double, double> pyne::decay_branch_ratio(std::pair<int, int> from_to) {
  return std::make_pair(data_access<double, decay>(from_to, offsetof(decay,
    branch_ratio), decay_data),data_access<double, decay>(from_to, offsetof(decay,
    branch_ratio_error), decay_data));
}

std::vector<double> pyne::decay_branch_ratios(int parent) {
  return data_access<double, decay>(parent, offsetof(decay,
    branch_ratio), decay_data);
}

std::pair<double, double> pyne::decay_photon_branch_ratio(std::pair<int,int>
from_to) {
  return std::make_pair(data_access<double, decay>(from_to,
    offsetof(decay, photon_branch_ratio), decay_data),
    data_access<double, decay>(from_to, offsetof(decay,
    photon_branch_ratio_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_photon_branch_ratios(
int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
    offsetof(decay, photon_branch_ratio), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
    offsetof(decay, photon_branch_ratio_error), decay_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::pair<double, double> pyne::decay_beta_branch_ratio(std::pair<int,int>
from_to) {
  return std::make_pair(data_access<double, decay>(from_to,
    offsetof(decay, beta_branch_ratio), decay_data),
    data_access<double, decay>(from_to, offsetof(decay,
    beta_branch_ratio_error), decay_data));
}

std::vector<std::pair<double, double> >pyne::decay_beta_branch_ratios(
int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, decay>(parent,
    offsetof(decay, beta_branch_ratio), decay_data);
  std::vector<double> part2 = data_access<double, decay>(parent,
    offsetof(decay, beta_branch_ratio_error), decay_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::map<std::pair<int, double>, pyne::gamma> pyne::gamma_data;

template<> void pyne::_load_data<pyne::gamma>() {

  // Loads the gamma table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(gamma));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(gamma, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(gamma, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "parent_nuc", HOFFSET(gamma, parent_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "child_nuc", HOFFSET(gamma, child_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "energy", HOFFSET(gamma, energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "energy_err", HOFFSET(gamma, energy_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity", HOFFSET(gamma,
                     photon_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity_err", HOFFSET(gamma,
                     photon_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity", HOFFSET(gamma,
                     conv_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity_err", HOFFSET(gamma,
                     conv_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity", HOFFSET(gamma,
                     total_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity_err", HOFFSET(gamma,
                     total_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(gamma, k_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(gamma, l_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(gamma, m_conv_e),
                     H5T_NATIVE_DOUBLE);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t gamma_set = H5Dopen2(nuc_data_h5, "/decay/gammas", H5P_DEFAULT);
  hid_t gamma_space = H5Dget_space(gamma_set);
  int gamma_length = H5Sget_simple_extent_npoints(gamma_space);

  // Read in the data
  gamma * gamma_array = new gamma[gamma_length];
  status = H5Dread(gamma_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   gamma_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(gamma_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < gamma_length; ++i) {
    if ((gamma_array[i].parent_nuc != 0) && !isnan(gamma_array[i].energy))
      gamma_data[std::make_pair(gamma_array[i].parent_nuc,
        gamma_array[i].energy)] = gamma_array[i];
  }
  delete[] gamma_array;
}

std::vector<std::pair<double, double> > pyne::gamma_energy(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, energy), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, energy_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_energy(double energy,
double error) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(energy+error,
    energy-error, offsetof(gamma, energy), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(energy+error,
    energy-error, offsetof(gamma, energy_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_photon_intensity(
int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, photon_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, photon_intensity_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_photon_intensity(
double energy, double error) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(energy+error,
    energy-error, offsetof(gamma, photon_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(energy+error,
    energy-error, offsetof(gamma, photon_intensity_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_conversion_intensity(
int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, conv_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, conv_intensity_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::gamma_total_intensity(
int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<double> part1 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, total_intensity), gamma_data);
  std::vector<double> part2 = data_access<double, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, total_intensity_err), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<int, int> > pyne::gamma_from_to(int parent) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, to_nuc), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<std::pair<int, int> > pyne::gamma_from_to(double energy,
double error) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(energy+error,
    energy-error, offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(energy+error,
    energy-error, offsetof(gamma, to_nuc), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}


std::vector<std::pair<int, int> > pyne::gamma_parent_child(double energy,
double error) {
  std::vector<std::pair<int, int> > result;
  std::vector<int> part1 = data_access<int, gamma>(energy+error,
    energy-error, offsetof(gamma, parent_nuc), gamma_data);
  std::vector<int> part2 = data_access<int, gamma>(energy+error,
    energy-error, offsetof(gamma, child_nuc), gamma_data);
  for(int i = 0; i < part1.size(); ++i){
    result.push_back(std::make_pair(part1[i],part2[i]));
  }
  return result;
}

std::vector<int> pyne::gamma_parent(double energy, double error) {
  return data_access<int, gamma>(energy+error, energy-error,
    offsetof(gamma, parent_nuc), gamma_data);
}

std::vector<int> pyne::gamma_child(double energy, double error) {
  return data_access<int, gamma>(energy+error, energy-error,
  offsetof(gamma, child_nuc), gamma_data);
}

std::vector<int> pyne::gamma_child(int parent) {
  return data_access<int, gamma>(parent, 0.0, DBL_MAX,
  offsetof(gamma, child_nuc), gamma_data);
}

std::vector<std::pair<double, double> > pyne::gamma_xrays(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, k_conv_e), gamma_data);
  std::vector<double> l_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, l_conv_e), gamma_data);
  std::vector<int> children = data_access<int, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, from_nuc), gamma_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
    decay_photon_branch_ratios(parent);
  for(int i = 0; i < k_list.size(); ++i){
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
          k_list[i]*decay_br[j].first, l_list[i]*decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  return result;
}

std::map<std::pair<int, double>, pyne::alpha> pyne::alpha_data;

template<> void pyne::_load_data<pyne::alpha>() {

  // Loads the alpha table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(alpha));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(alpha, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(alpha, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "energy", HOFFSET(alpha, energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(alpha, intensity),
                     H5T_NATIVE_DOUBLE);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t alpha_set = H5Dopen2(nuc_data_h5, "/decay/alphas", H5P_DEFAULT);
  hid_t alpha_space = H5Dget_space(alpha_set);
  int alpha_length = H5Sget_simple_extent_npoints(alpha_space);

  // Read in the data
  alpha * alpha_array = new alpha[alpha_length];
  status = H5Dread(alpha_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   alpha_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(alpha_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < alpha_length; ++i) {
    if ((alpha_array[i].from_nuc != 0) && !isnan(alpha_array[i].energy))
      alpha_data[std::make_pair(alpha_array[i].from_nuc, alpha_array[i].energy)]
        = alpha_array[i];
  }
  delete[] alpha_array;
}

std::vector<double > pyne::alpha_energy(int parent){
  return data_access<double, alpha>(parent, 0.0, DBL_MAX,
                     offsetof(alpha,energy), alpha_data);
}
std::vector<double> pyne::alpha_intensity(int parent){
  return data_access<double, alpha>(parent, 0.0, DBL_MAX,
                     offsetof(alpha,intensity), alpha_data);
}

std::vector<int> pyne::alpha_parent(double energy, double error) {
  return data_access<int, alpha>(energy+error, energy-error,
                     offsetof(alpha, from_nuc), alpha_data);
}

std::vector<int> pyne::alpha_child(double energy, double error) {
  return data_access<int, alpha>(energy+error, energy-error,
                     offsetof(alpha, to_nuc), alpha_data);
}

std::vector<int> pyne::alpha_child(int parent){
  return data_access<int, alpha>(parent, 0.0, DBL_MAX,
                     offsetof(alpha, to_nuc), alpha_data);
}

std::map<std::pair<int, double>, pyne::beta> pyne::beta_data;

template<> void pyne::_load_data<pyne::beta>() {

  // Loads the beta table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(beta));
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(beta,
                     endpoint_energy), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(beta, avg_energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(beta, intensity),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(beta, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(beta, to_nuc),
                     H5T_NATIVE_INT);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t beta_set = H5Dopen2(nuc_data_h5, "/decay/betas", H5P_DEFAULT);
  hid_t beta_space = H5Dget_space(beta_set);
  int beta_length = H5Sget_simple_extent_npoints(beta_space);

  // Read in the data
  beta * beta_array = new beta[beta_length];
  status = H5Dread(beta_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, beta_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(beta_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < beta_length; ++i) {
    if ((beta_array[i].from_nuc != 0) && !isnan(beta_array[i].avg_energy))
      beta_data[std::make_pair(beta_array[i].from_nuc, beta_array[i].avg_energy)]
        = beta_array[i];
  }
  delete[] beta_array;
}

std::vector<double > pyne::beta_endpoint_energy(int parent){
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                     offsetof(beta, endpoint_energy), beta_data);
}

std::vector<double > pyne::beta_average_energy(int parent){
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                     offsetof(beta, avg_energy), beta_data);
}

std::vector<double> pyne::beta_intensity(int parent){
  return data_access<double, beta>(parent, 0.0, DBL_MAX,
                     offsetof(beta, intensity), beta_data);
}

std::vector<int> pyne::beta_parent(double energy, double error) {
  return data_access<int, beta>(energy+error, energy-error,
                     offsetof(beta, from_nuc), beta_data);
}

std::vector<int> pyne::beta_child(double energy, double error) {
  return data_access<int, beta>(energy+error, energy-error,
                     offsetof(beta, to_nuc), beta_data);
}

std::vector<int> pyne::beta_child(int parent){
  return data_access<int, beta>(parent, 0.0, DBL_MAX,
                     offsetof(beta, to_nuc),beta_data);
}


std::map<std::pair<int, double>, pyne::ecbp> pyne::ecbp_data;

template<> void pyne::_load_data<pyne::ecbp>() {

  // Loads the ecbp table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ecbp));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ecbp, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ecbp, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(ecbp,
                     endpoint_energy),H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(ecbp, avg_energy),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_plus_intensity", HOFFSET(ecbp,
                     beta_plus_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ec_intensity", HOFFSET(ecbp, ec_intensity),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(ecbp, k_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(ecbp, l_conv_e),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(ecbp, m_conv_e),
                     H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY,
                              H5P_DEFAULT);

  // Open the data set
  hid_t ecbp_set = H5Dopen2(nuc_data_h5, "/decay/ecbp", H5P_DEFAULT);
  hid_t ecbp_space = H5Dget_space(ecbp_set);
  int ecbp_length = H5Sget_simple_extent_npoints(ecbp_space);

  // Read in the data
  ecbp * ecbp_array = new ecbp[ecbp_length];
  status = H5Dread(ecbp_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecbp_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(ecbp_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < ecbp_length; ++i) {
    if ((ecbp_array[i].from_nuc != 0) && !isnan(ecbp_array[i].avg_energy))
      ecbp_data[std::make_pair(ecbp_array[i].from_nuc, ecbp_array[i].avg_energy)]
        = ecbp_array[i];
  }
  delete[] ecbp_array;
}

std::vector<double > pyne::ecbp_endpoint_energy(int parent){
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp,endpoint_energy), ecbp_data);
}

std::vector<double > pyne::ecbp_average_energy(int parent){
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp, avg_energy), ecbp_data);
}

std::vector<double> pyne::ec_intensity(int parent){
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp, ec_intensity), ecbp_data);
}

std::vector<double> pyne::bp_intensity(int parent){
  return data_access<double, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp, beta_plus_intensity), ecbp_data);
}

std::vector<int> pyne::ecbp_parent(double energy, double error) {
  return data_access<int, ecbp>(energy+error, energy-error,
                     offsetof(ecbp, from_nuc), ecbp_data);
}

std::vector<int> pyne::ecbp_child(double energy, double error) {
  return data_access<int, ecbp>(energy+error, energy-error,
                     offsetof(ecbp, to_nuc), ecbp_data);
}

std::vector<int> pyne::ecbp_child(int parent){
  return data_access<int, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp, to_nuc), ecbp_data);
}

std::vector<std::pair<double, double> > pyne::ecbp_xrays(int parent) {
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
    offsetof(ecbp, k_conv_e), ecbp_data);
  std::vector<double> l_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
    offsetof(ecbp, l_conv_e), ecbp_data);
  std::vector<int> children = data_access<int, ecbp>(parent, 0.0, DBL_MAX,
    offsetof(ecbp, to_nuc), ecbp_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
    decay_beta_branch_ratios(parent);
  for(int i = 0; i < k_list.size(); ++i){
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
          k_list[i]*decay_br[j].first, l_list[i]*decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  return result;
}

///////////////////////////
/// Combined decay data ///
///////////////////////////

std::vector<std::pair<double, double> > pyne::gammas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<std::pair<double, double> > energies = gamma_energy(parent_state_id);
  std::vector<std::pair<double, double> > intensities =
    gamma_photon_intensity(parent_state_id);
  std::vector<std::pair<int, int> > children = gamma_from_to(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<std::pair<double, double> > decay_br =
    decay_photon_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i].first) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i].first,
          decay_c*intensities[i].first*decay_br[j].first));
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::alphas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<double> energies = alpha_energy(parent_state_id);
  std::vector<double> intensities = alpha_intensity(parent_state_id);
  std::vector<int> children = alpha_child(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<double> decay_br = decay_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i],
          decay_c*decay_br[j]*intensities[i]));
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::betas(int parent_state_id) {
  std::vector<std::pair<double, double> > result;
  double decay_c = decay_const(parent_state_id);
  std::vector<double> energies = beta_average_energy(parent_state_id);
  std::vector<double> intensities = beta_intensity(parent_state_id);
  std::vector<int> children = beta_child(parent_state_id);
  std::vector<int> decay_children = decay_data_children(parent_state_id);
  std::vector<std::pair<double, double> > decay_br =
    decay_beta_branch_ratios(parent_state_id);
  for (int i = 0; i < children.size(); ++i) {
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        result.push_back(std::make_pair(energies[i],
          decay_c*decay_br[j].first*intensities[i]));
        break;
      }
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::xrays(int parent) {
  double decay_c = decay_const(parent);
  std::vector<std::pair<double, double> > result;
  std::vector<std::pair<double, double> > temp;
  std::vector<double> k_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
    offsetof(ecbp, k_conv_e), ecbp_data);
  std::vector<double> l_list = data_access<double, ecbp>(parent, 0.0, DBL_MAX,
    offsetof(ecbp, l_conv_e), ecbp_data);
  std::vector<int> children = data_access<int, ecbp>(parent, 0.0, DBL_MAX,
                     offsetof(ecbp, to_nuc), ecbp_data);
  std::vector<int> decay_children = decay_data_children(parent);
  std::vector<std::pair<double, double> > decay_br =
    decay_beta_branch_ratios(parent);
  for(int i = 0; i < k_list.size(); ++i){
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(children[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(children[i]),
          k_list[i]*decay_br[j].first, l_list[i]*decay_br[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }
  std::vector<double> gk_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, k_conv_e), gamma_data);
  std::vector<double> gl_list = data_access<double, gamma>(parent, 0.0, DBL_MAX,
    offsetof(gamma, l_conv_e), gamma_data);
  std::vector<int> gchildren = data_access<int, gamma>(parent, 0.0,
    DBL_MAX, offsetof(gamma, from_nuc), gamma_data);
  std::vector<std::pair<double, double> > decay_nrbr =
    decay_photon_branch_ratios(parent);
  for(int i = 0; i < gk_list.size(); ++i){
    for (int j = 0; j < decay_children.size(); ++j) {
      if (nucname::zzzaaa(gchildren[i]) == nucname::zzzaaa(decay_children[j])) {
        temp = calculate_xray_data(nucname::znum(gchildren[i]),
          gk_list[i]*decay_nrbr[j].first, gl_list[i]*decay_nrbr[j].first);
        for (int k = 0; k < temp.size(); ++k) {
          if (!isnan(temp[k].second) && !isnan(temp[k].first)) {
            int found = 0;
            for (int l = 0; l < result.size(); ++l) {
              if (temp[k].first == result[l].first) {
                result[l].second = result[l].second + temp[k].second;
                found = 1;
                break;
              }
            }
            if (found == 0) {
              result.push_back(temp[k]);
            }
          }
        }
        break;
      }
    }
  }

  for(int i = 0; i < result.size(); ++i)
    result[i].second = result[i].second * decay_c;
  return result;
}

//////////////////////////////////////////
//////////// simple xs data //////////////
//////////////////////////////////////////

typedef struct simple_xs {
  int nuc;
  double sigma_t;
  double sigma_s;
  double sigma_e;
  double sigma_i;
  double sigma_a;
  double sigma_gamma;
  double sigma_f;
  double sigma_alpha;
  double sigma_proton;
  double sigma_deut;
  double sigma_trit;
  double sigma_2n;
  double sigma_3n;
  double sigma_4n;
} simple_xs;

std::map<std::string, std::map<int, std::map<int, double> > > pyne::simple_xs_map;

// loads the simple cross section data for the specified energy band from
// the nuc_data.h5 file into memory.
static void _load_simple_xs_map(std::string energy) {
  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  using pyne::rxname::id;
  std::map<unsigned int, size_t> rxns;
  rxns[id("tot")] = offsetof(simple_xs, sigma_t);
  rxns[id("scat")] = offsetof(simple_xs, sigma_s);
  rxns[id("elas")] = offsetof(simple_xs, sigma_e);
  rxns[id("inel")] = offsetof(simple_xs, sigma_i);
  rxns[id("abs")] = offsetof(simple_xs, sigma_a);
  rxns[id("gamma")] = offsetof(simple_xs, sigma_gamma);
  rxns[id("fiss")] = offsetof(simple_xs, sigma_f);
  rxns[id("alpha")] = offsetof(simple_xs, sigma_alpha);
  rxns[id("proton")] = offsetof(simple_xs, sigma_proton);
  rxns[id("deut")] = offsetof(simple_xs, sigma_deut);
  rxns[id("trit")] = offsetof(simple_xs, sigma_trit);
  rxns[id("z_2n")] = offsetof(simple_xs, sigma_2n);
  rxns[id("z_3n")] = offsetof(simple_xs, sigma_3n);
  rxns[id("z_4n")] = offsetof(simple_xs, sigma_4n);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(simple_xs));
  H5Tinsert(desc, "nuc",   HOFFSET(simple_xs, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "sigma_t",  HOFFSET(simple_xs, sigma_t),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_s", HOFFSET(simple_xs, sigma_s), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_e", HOFFSET(simple_xs, sigma_e), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_i", HOFFSET(simple_xs, sigma_i), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_a", HOFFSET(simple_xs, sigma_a), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_gamma", HOFFSET(simple_xs, sigma_gamma), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_f", HOFFSET(simple_xs, sigma_f), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_alpha", HOFFSET(simple_xs, sigma_alpha), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_proton", HOFFSET(simple_xs, sigma_proton), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_deut", HOFFSET(simple_xs, sigma_deut), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_trit", HOFFSET(simple_xs, sigma_trit), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_2n", HOFFSET(simple_xs, sigma_2n), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_3n", HOFFSET(simple_xs, sigma_3n), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "sigma_4n", HOFFSET(simple_xs, sigma_4n), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // build path to prober simple xs table
  std::string path = "/neutron/simple_xs/" + energy;

  // Open the data set
  hid_t simple_xs_set = H5Dopen2(nuc_data_h5, path.c_str(), H5P_DEFAULT);
  hid_t simple_xs_space = H5Dget_space(simple_xs_set);
  int n = H5Sget_simple_extent_npoints(simple_xs_space);

  // Read in the data
  simple_xs* array = new simple_xs[n];
  H5Dread(simple_xs_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

  // close the nuc_data library, before doing anything stupid
  H5Dclose(simple_xs_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the map
  for(int i = 0; i < n; i++) {
    std::map<unsigned int, size_t>::iterator it;
    for (it = rxns.begin(); it != rxns.end(); ++it) {
      double xs = *(double*)((char*)&array[i] + it->second);
      pyne::simple_xs_map[energy][array[i].nuc][it->first] = xs;
    }
  }
  delete[] array;
}

double pyne::simple_xs(int nuc, int rx_id, std::string energy) {
  std::set<std::string> energies;
  energies.insert("thermal");
  energies.insert("thermal_maxwell_ave");
  energies.insert("resonance_integral");
  energies.insert("fourteen_MeV");
  energies.insert("fission_spectrum_ave");

  if (energies.count(energy) == 0) {
    throw InvalidSimpleXS("Energy '" + energy +
        "' is not a valid simple_xs group");
  } else if (simple_xs_map.count(energy) == 0) {
    _load_simple_xs_map(energy);
  }

  if (simple_xs_map[energy].count(nuc) == 0) {
    throw InvalidSimpleXS(rxname::name(rx_id) +
        " is not a valid simple_xs nuclide");
  } else if (simple_xs_map[energy][nuc].count(rx_id) == 0) {
    throw InvalidSimpleXS(rxname::name(rx_id) +
        " is not a valid simple_xs reaction");
  }

  return simple_xs_map[energy][nuc][rx_id];
}

double pyne::simple_xs(int nuc, std::string rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
double pyne::simple_xs(std::string nuc, int rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
double pyne::simple_xs(std::string nuc, std::string rx, std::string energy) {
  return pyne::simple_xs(nucname::id(nuc), rxname::id(rx), energy);
}
