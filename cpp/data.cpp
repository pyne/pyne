// Implements basic nuclear data functions.
#ifndef PYNE_IS_AMALGAMATED
#include "data.h"
#endif


//
// Math Helpers
//

const double pyne::pi = 3.14159265359;
const double pyne::N_A = 6.0221415e+23;
const double pyne::barns_per_cm2 = 1e24;
const double pyne::cm2_per_barn = 1e-24;
const double pyne::sec_per_day = 24.0 * 3600.0;

/********************************/
/*** data_checksums Functions ***/
/********************************/

std::map<std::string, std::string> pyne::get_data_checksums()
{
    std::map<std::string, std::string> temp_map;
    // Initialization of dataset hashes
    temp_map["/atomic_decay"]="09bf73252629077785e20b3532fde8b3";
    temp_map["/atomic_mass"]="10edfdc662e35bdfab91beb89285efff";
    temp_map["/material_library"]="8b10864378fbd88538434679acf908cc";
    temp_map["/neutron/eaf_xs"]="29622c636c4a3a46802207b934f9516c";
    temp_map["/neutron/scattering_lengths"]="a24d391cc9dc0fc146392740bb97ead4";
    temp_map["/neutron/simple_xs"]="3d6e086977783dcdf07e5c6b0c2416be";
    
    return temp_map;
};

std::map<std::string, std::string> pyne::data_checksums = pyne::get_data_checksums();

/*****************************/
/*** atomic_mass Functions ***/
/*****************************/
std::map<int, double> pyne::atomic_mass_map = std::map<int, double>();

void pyne::_load_atomic_mass_map()
{
  // Loads the important parts of atomic_wight table into atomic_mass_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic_mass_struct));
  H5Tinsert(desc, "nuc",   HOFFSET(atomic_mass_struct, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "mass",  HOFFSET(atomic_mass_struct, mass),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "error", HOFFSET(atomic_mass_struct, error), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "abund", HOFFSET(atomic_mass_struct, abund), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atomic_mass_set = H5Dopen2(nuc_data_h5, "/atomic_mass", H5P_DEFAULT);
  hid_t atomic_mass_space = H5Dget_space(atomic_mass_set);
  int atomic_mass_length = H5Sget_simple_extent_npoints(atomic_mass_space);

  // Read in the data
  atomic_mass_struct * atomic_mass_array = new atomic_mass_struct[atomic_mass_length];
  H5Dread(atomic_mass_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_mass_array);

  // close the nuc_data library, before doing anythng stupid
  H5Dclose(atomic_mass_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the map
  for(int n = 0; n < atomic_mass_length; n++){
    atomic_mass_map[atomic_mass_array[n].nuc] = atomic_mass_array[n].mass;
    natural_abund_map[atomic_mass_array[n].nuc] = atomic_mass_array[n].abund;
  }

  delete[] atomic_mass_array;
};


double pyne::atomic_mass(int nuc)
{
  // Find the nuclide's mass in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = atomic_mass_map.find(nuc);
  nuc_end = atomic_mass_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (atomic_mass_map.empty())
  {
    // Don't fail if we can't load the library
    try
    {
      _load_atomic_mass_map();
      return atomic_mass(nuc);
    }
    catch(...){};
  };

  double aw;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state mass...not strictly true, but good guess.
  if (0 < nucid%10000)
  {
    aw = atomic_mass((nucid/10000)*10000);
    atomic_mass_map[nuc] = aw;
    return aw;
  };

  // Finally, if none of these work, 
  // take a best guess based on the 
  // aaa number.
  aw = (double) ((nucid/10000)%1000);
  atomic_mass_map[nuc] = aw;
  return aw;
};


double pyne::atomic_mass(char * nuc)
{
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
};


double pyne::atomic_mass(std::string nuc)
{
  int nuc_zz = nucname::id(nuc);
  return atomic_mass(nuc_zz);
};


/*******************************/
/*** natural_abund functions ***/
/*******************************/

std::map<int, double> pyne::natural_abund_map = std::map<int, double>();

double pyne::natural_abund(int nuc)
{
  // Find the nuclide's natural abundance
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = natural_abund_map.find(nuc);
  nuc_end = natural_abund_map.end();

  // First check if we already have the nuc mass in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (natural_abund_map.empty())
  {
    // Don't fail if we can't load the library
    try
    {
      _load_atomic_mass_map();
      return natural_abund(nuc);
    }
    catch(...){};
  };

  double na;
  int nucid = nucname::id(nuc);

  // If in an excited state, return the ground
  // state abundance...not strictly true, but good guess.
  if (0 < nucid%10000)
  {
    na = natural_abund((nucid/10000)*10000);
    atomic_mass_map[nuc] = na;
    return na;
  };

  // Finally, if none of these work, 
  // take a best guess based on the 
  // aaa number.
  na = 0.0;
  natural_abund_map[nuc] = na;
  return na;
};


double pyne::natural_abund(char * nuc)
{
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
};


double pyne::natural_abund(std::string nuc)
{
  int nuc_zz = nucname::id(nuc);
  return natural_abund(nuc_zz);
};




/***********************************/
/*** scattering length functions ***/
/***********************************/
std::map<int, extra_types::complex_t> pyne::b_coherent_map = std::map<int, extra_types::complex_t>();
std::map<int, extra_types::complex_t> pyne::b_incoherent_map = std::map<int, extra_types::complex_t>();
std::map<int, double> pyne::b_map = std::map<int, double>();


void pyne::_load_scattering_lengths()
{
  // Loads the important parts of atomic_wight table into atomic_mass_map
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(scattering_lengths_struct));
  status = H5Tinsert(desc, "nuc", HOFFSET(scattering_lengths_struct, nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "b_coherent", HOFFSET(scattering_lengths_struct, b_coherent), 
                      h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "b_incoherent", HOFFSET(scattering_lengths_struct, b_incoherent), 
                      h5wrap::PYTABLES_COMPLEX128);
  status = H5Tinsert(desc, "xs_coherent", HOFFSET(scattering_lengths_struct, xs_coherent), 
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs_incoherent", HOFFSET(scattering_lengths_struct, xs_incoherent), 
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "xs", HOFFSET(scattering_lengths_struct, xs), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t scat_len_set = H5Dopen2(nuc_data_h5, "/neutron/scattering_lengths", H5P_DEFAULT);
  hid_t scat_len_space = H5Dget_space(scat_len_set);
  int scat_len_length = H5Sget_simple_extent_npoints(scat_len_space);

  // Read in the data
  scattering_lengths_struct * scat_len_array = new scattering_lengths_struct[scat_len_length];
  status = H5Dread(scat_len_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, scat_len_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(scat_len_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for(int n = 0; n < scat_len_length; n++)
  {
    b_coherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_coherent;
    b_incoherent_map[scat_len_array[n].nuc] = scat_len_array[n].b_incoherent;
  };

  delete[] scat_len_array;
};



//
// Coherent functions 
//


extra_types::complex_t pyne::b_coherent(int nuc)
{
  // Find the nuclide's bound scattering length in cm
  std::map<int, extra_types::complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_coherent_map.find(nuc);
  nuc_end = b_coherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (b_coherent_map.empty())
  {
    _load_scattering_lengths();
    return b_coherent(nuc);
  };

  extra_types::complex_t bc;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (anum == nucname::anum((*nuc_iter).first))
    {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    };
    nuc_iter++;
  };

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (znum == nucname::znum((*nuc_iter).first))
    {
      bc = (*nuc_iter).second;
      b_coherent_map[nuc] = bc;
      return bc;
    };
    nuc_iter++;
  };

  // Finally, if none of these work, 
  // just return zero...
  bc.re = 0.0;
  bc.im = 0.0;
  b_coherent_map[nuc] = bc;
  return bc;
};


extra_types::complex_t pyne::b_coherent(char * nuc)
{
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
};


extra_types::complex_t pyne::b_coherent(std::string nuc)
{
  int nuc_zz = nucname::id(nuc);
  return b_coherent(nuc_zz);
};



//
// Incoherent functions 
//


extra_types::complex_t pyne::b_incoherent(int nuc)
{
  // Find the nuclide's bound inchoherent scattering length in cm
  std::map<int, extra_types::complex_t>::iterator nuc_iter, nuc_end;

  nuc_iter = b_incoherent_map.find(nuc);
  nuc_end = b_incoherent_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (b_incoherent_map.empty())
  {
    _load_scattering_lengths();
    return b_incoherent(nuc);
  };

  extra_types::complex_t bi;
  int nucid = nucname::id(nuc);
  int znum = nucname::znum(nucid);
  int anum = nucname::anum(nucid);

  // Try to find a nuclide with matching A-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (anum == nucname::anum((*nuc_iter).first))
    {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    };
    nuc_iter++;
  };

  // Try to find a nuclide with matching Z-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (znum == nucname::znum((*nuc_iter).first))
    {
      bi = (*nuc_iter).second;
      b_incoherent_map[nuc] = bi;
      return bi;
    };
    nuc_iter++;
  };

  // Finally, if none of these work, 
  // just return zero...
  bi.re = 0.0;
  bi.im = 0.0;
  b_incoherent_map[nuc] = bi;
  return bi;
};


extra_types::complex_t pyne::b_incoherent(char * nuc)
{
  return b_incoherent(nucname::id(nuc));
};


extra_types::complex_t pyne::b_incoherent(std::string nuc)
{
  return b_incoherent(nucname::id(nuc));
};



//
// b functions
//

double pyne::b(int nuc)
{
  // Find the nuclide's bound scattering length in cm
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = b_map.find(nuc);
  nuc_end = b_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, calculate the value from coherent and incoherent lengths
  extra_types::complex_t bc = b_coherent(nuc);
  extra_types::complex_t bi = b_incoherent(nuc);

  double b_val = sqrt(bc.re*bc.re + bc.im*bc.im + bi.re*bi.re + bi.im*bi.im);

  return b_val;
};


double pyne::b(char * nuc)
{
  int nucid = nucname::id(nuc);
  return b(nucid);
};


double pyne::b(std::string nuc)
{
  int nucid = nucname::id(nuc);
  return b(nucid);
};



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
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(wimsdfpy_struct));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(wimsdfpy_struct, from_nuc), 
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(wimsdfpy_struct, to_nuc), 
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yields", HOFFSET(wimsdfpy_struct, yields),
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
  wimsdfpy_struct * wimsdfpy_array = new wimsdfpy_struct[wimsdfpy_length];
  status = H5Dread(wimsdfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, wimsdfpy_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(wimsdfpy_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  for(int n=0; n < wimsdfpy_length; n++) {
    wimsdfpy_data[std::make_pair(wimsdfpy_array[n].from_nuc, 
      wimsdfpy_array[n].to_nuc)] = wimsdfpy_array[n].yields;
  };

  delete[] wimsdfpy_array;
};


std::map<std::pair<int, int>, pyne::ndsfpysub_struct> pyne::ndsfpy_data = \
  std::map<std::pair<int, int>, pyne::ndsfpysub_struct>();

void pyne::_load_ndsfpy() {
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ndsfpy_struct));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ndsfpy_struct, from_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ndsfpy_struct, to_nuc),
                     H5T_NATIVE_INT);
  status = H5Tinsert(desc, "yield_thermal", HOFFSET(ndsfpy_struct, yield_thermal),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_thermal_err", HOFFSET(ndsfpy_struct, yield_thermal_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast", HOFFSET(ndsfpy_struct, yield_fast),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_fast_err", HOFFSET(ndsfpy_struct, yield_fast_err),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV", HOFFSET(ndsfpy_struct, yield_14MeV),
                     H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "yield_14MeV_err", HOFFSET(ndsfpy_struct, yield_14MeV_err),
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
  ndsfpy_struct * ndsfpy_array = new ndsfpy_struct[ndsfpy_length];
  status = H5Dread(ndsfpy_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ndsfpy_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(ndsfpy_set);
  status = H5Fclose(nuc_data_h5);

  ndsfpysub_struct ndsfpysub_temp;

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
  };



  delete[] ndsfpy_array;
};

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
    std::map<std::pair<int, int>, ndsfpysub_struct>::iterator fpy_iter, fpy_end;
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
};

double pyne::fpyield(int from_nuc, int to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc), 
                                     nucname::id(to_nuc)), source, get_error);
};

double pyne::fpyield(char * from_nuc, char * to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc), 
                                     nucname::id(to_nuc)), source, get_error);
};

double pyne::fpyield(std::string from_nuc, std::string to_nuc, int source, bool get_error) {
  return fpyield(std::pair<int, int>(nucname::id(from_nuc), 
                                     nucname::id(to_nuc)), source, get_error);
};



/***********************/
/*** decay functions ***/
/***********************/
std::map<int, double> pyne::half_life_map = std::map<int, double>();
std::map<int, double> pyne::decay_const_map = std::map<int, double>();
std::map<std::pair<int, int>, double> pyne::branch_ratio_map = \
                                              std::map<std::pair<int, int>, double>();
std::map<int, double> pyne::state_energy_map = std::map<int, double>();
std::map<int, std::set<int> > pyne::decay_children_map = \
                                                      std::map<int, std::set<int> >();


void pyne::_load_half_life_decay()
{
  // Loads the important parts of half_life_decay table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(half_life_decay_struct));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(half_life_decay_struct, from_nuc),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "level", HOFFSET(half_life_decay_struct, level), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(half_life_decay_struct, to_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "half_life", HOFFSET(half_life_decay_struct, half_life),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "decay_const", HOFFSET(half_life_decay_struct, decay_const),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(half_life_decay_struct, branch_ratio),
                      H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atom_dec_set = H5Dopen2(nuc_data_h5, "/decay/half_life", H5P_DEFAULT);
  hid_t atom_dec_space = H5Dget_space(atom_dec_set);
  int atom_dec_length = H5Sget_simple_extent_npoints(atom_dec_space);

  // Read in the data
  half_life_decay_struct * atom_dec_array = new half_life_decay_struct[atom_dec_length];
  status = H5Dread(atom_dec_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, atom_dec_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(atom_dec_set);
  status = H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the maps
  // giving precednece to ground state values or those seen first.
  int from_nuc, to_nuc;
  double level;
  std::pair<int, int> from_to;
  for(int n = 0; n < atom_dec_length; n++)
  {
    from_nuc = atom_dec_array[n].from_nuc;
    level = atom_dec_array[n].level;
    to_nuc = atom_dec_array[n].to_nuc;
    from_to = std::pair<int, int>(from_nuc, to_nuc);

    if (0 == half_life_map.count(from_nuc) || 0.0 == level)
      half_life_map[from_nuc] = atom_dec_array[n].half_life;

    if (0 == decay_const_map.count(from_nuc) || 0.0 == level)
      decay_const_map[from_nuc] = atom_dec_array[n].decay_const;

    if (0 == branch_ratio_map.count(from_to) || 0.0 == level)
      branch_ratio_map[from_to] = atom_dec_array[n].branch_ratio;

    state_energy_map[from_nuc] = level;

    if (0.0 != atom_dec_array[n].decay_const)
      decay_children_map[from_nuc].insert(to_nuc);
  };

  delete[] atom_dec_array;
};


//
// Decay constant data
//

double pyne::decay_const(int nuc)
{
  // Find the nuclide's decay constant in 1/s
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = decay_const_map.find(nuc);
  nuc_end = decay_const_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (decay_const_map.empty())
  {
    _load_half_life_decay();
    return decay_const(nuc);
  };

  // Finally, if none of these work, 
  // assume the value is stable
  double dc = 0.0;
  decay_const_map[nuc] = dc;
  return dc;
};


double pyne::decay_const(char * nuc)
{
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
};


double pyne::decay_const(std::string nuc)
{
  int nuc_zz = nucname::id(nuc);
  return decay_const(nuc_zz);
};




//
// Branch ratio data
//

double pyne::branch_ratio(std::pair<int, int> from_to)
{
  // Find the parent/child pair branch ratio as a fraction
  std::map<std::pair<int, int>, double>::iterator br_iter, br_end;

  br_iter = branch_ratio_map.find(from_to);
  br_end = branch_ratio_map.end();

  // First check if we already have the pair in the map
  if (br_iter != br_end)
    return (*br_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (branch_ratio_map.empty())
  {
    _load_half_life_decay();
    return branch_ratio(from_to);
  };

  // Finally, if none of these work, 
  // assume the value is stable
  double br = 0.0;
  branch_ratio_map[from_to] = br;
  return br;
};


double pyne::branch_ratio(int from_nuc, int to_nuc)
{
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc), 
                                          nucname::id(to_nuc)));
};

double pyne::branch_ratio(char * from_nuc, char * to_nuc)
{
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc), 
                                          nucname::id(to_nuc)));
};

double pyne::branch_ratio(std::string from_nuc, std::string to_nuc)
{
  return branch_ratio(std::pair<int, int>(nucname::id(from_nuc), 
                                          nucname::id(to_nuc)));
};

//
// Decay children data
//

std::set<int> pyne::decay_children(int nuc)
{
  // Find the nuclide's decay constant in 1/s
  std::map<int, std::set<int> >::iterator nuc_iter, nuc_end;

  nuc_iter = decay_children_map.find(nuc);
  nuc_end = decay_children_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end) {
    return (*nuc_iter).second;
  };

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (decay_children_map.empty())
  {
    _load_half_life_decay();
    return decay_children(nuc);
  };

  // Finally, if none of these work, 
  // assume the value is stable
  std::set<int> dc = std::set<int>();
  decay_children_map[nuc] = dc;
  return dc;
};


std::set<int> pyne::decay_children(char * nuc)
{
  return decay_children(nucname::id(nuc));
};


std::set<int> pyne::decay_children(std::string nuc)
{
  return decay_children(nucname::id(nuc));
};

std::map<int, pyne::level_struct> pyne::level_data = \
  std::map<int, pyne::level_struct>();

void pyne::_load_level_data()
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
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(level_struct));
  status = H5Tinsert(desc, "nuc_id", HOFFSET(level_struct, nuc_id),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "level", HOFFSET(level_struct, level), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "half_life", HOFFSET(level_struct, half_life),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "metastable", HOFFSET(level_struct, metastable),
                      H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t level_set = H5Dopen2(nuc_data_h5, "/decay/level_list", H5P_DEFAULT);
  hid_t level_space = H5Dget_space(level_set);
  int level_length = H5Sget_simple_extent_npoints(level_space);

  // Read in the data
  level_struct * level_array = new level_struct[level_length];
  status = H5Dread(level_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, level_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(level_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < level_length; ++i) {
    level_data[level_array[i].nuc_id] = level_array[i];
  }

}

int pyne::metastable_id(int nuc, int m) {
  if (m==0) return 0;
  int nostate = (nuc / 1000) * 1000;
  if (level_data.empty()) {
    _load_level_data();
  }

  std::map<int, pyne::level_struct>::iterator nuc_lower, nuc_upper;

  nuc_lower = level_data.lower_bound(nostate);
  nuc_upper = level_data.upper_bound(nostate+9999);
  for (std::map<int, pyne::level_struct>::iterator it=nuc_lower; it!=nuc_upper; ++it) {
    if (it->second.metastable == m)
        return it->second.nuc_id;
  }

  return 0;
}

int pyne::metastable_id(int nuc) {
  return metastable_id(nuc, 1);
}

//
// Excitation state energy data
//

double pyne::state_energy(int nuc)
{
  // Find the nuclide's state energy in MeV
  double energy = 0.0;

    // Find the nuclide's state energy in MeV
  std::map<int, pyne::level_struct>::iterator nuc_iter, nuc_end;

  nuc_iter = level_data.find(nuc);
  nuc_end = level_data.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
     return nuc_iter->second.level;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (level_data.empty())
  {
    _load_level_data();
    state_energy(nuc);
  };
  return 0.0;
};


double pyne::state_energy(char * nuc)
{
  return state_energy(nucname::id(nuc));
};


double pyne::state_energy(std::string nuc)
{
  return state_energy(nucname::id(nuc));
};

//
// Half-life data
//

double pyne::half_life(int nuc)
{

    // Find the nuclide's state energy in MeV
  std::map<int, pyne::level_struct>::iterator nuc_iter, nuc_end;

  nuc_iter = level_data.find(nuc);
  nuc_end = level_data.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
     return nuc_iter->second.half_life;

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (level_data.empty())
  {
    _load_level_data();
    return half_life(nuc);
  };
  return 0.0;
};


double pyne::half_life(char * nuc)
{
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
};


double pyne::half_life(std::string nuc)
{
  int nuc_zz = nucname::id(nuc);
  return half_life(nuc_zz);
};

std::map<std::pair<int, int>, pyne::decay_struct> pyne::decay_data = \
  std::map<std::pair<int, int>, pyne::decay_struct>();

void pyne::_load_decay_data()
{

  // Loads the decay table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(decay_struct));
  status = H5Tinsert(desc, "parent", HOFFSET(decay_struct, parent),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "daughter", HOFFSET(decay_struct, daughter),
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "decay", HOFFSET(decay_struct, decay), H5T_STRING);
  status = H5Tinsert(desc, "half_life", HOFFSET(decay_struct, half_life), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "half_life_error", HOFFSET(decay_struct, half_life_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(decay_struct, branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio", HOFFSET(decay_struct, photon_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_branch_ratio_error", HOFFSET(decay_struct, photon_branch_ratio_error), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio", HOFFSET(decay_struct, beta_branch_ratio), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_branch_ratio_error", HOFFSET(decay_struct, beta_branch_ratio_error), H5T_NATIVE_DOUBLE);
  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t decay_set = H5Dopen2(nuc_data_h5, "/decay/decays", H5P_DEFAULT);
  hid_t decay_space = H5Dget_space(decay_set);
  int decay_length = H5Sget_simple_extent_npoints(decay_space);

  // Read in the data
  decay_struct * decay_array = new decay_struct[decay_length];
  status = H5Dread(decay_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, decay_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(decay_set);
  status = H5Fclose(nuc_data_h5);

  for (int i = 0; i < decay_length; ++i) {
    decay_data[std::make_pair(decay_array[i].parent, decay_array[i].daughter)] = decay_array[i];
  }

}


template<typename T> void pyne::decay_data_access(std::pair<int, int> from_to, T &a, size_t valoffset) {
  // Find the nuclide's state energy in MeV
  std::map<std::pair<int, int> , pyne::decay_struct>::iterator nuc_iter, nuc_end;

  nuc_iter = decay_data.find(from_to);
  nuc_end = decay_data.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    a = T ((char*)&(nuc_iter->second) + valoffset);

  // Next, fill up the map with values from the
  // nuc_data.h5, if the map is empty.
  if (decay_data.empty())
  {
    _load_decay_data();
    decay_data_access(from_to, a, valoffset);
  };

}


void pyne::decay_data_byparentchild(std::pair<int, int> from_to, decay_struct *data){

  if (decay_data.empty()) {
    _load_decay_data();
    decay_data_byparentchild(from_to, data);
  }

  data = (decay_struct*) malloc(sizeof(decay_struct));

  std::map<std::pair<int, int> , pyne::decay_struct>::iterator nuc_iter, nuc_end;

  nuc_iter = decay_data.find(from_to);
  nuc_end = decay_data.end();

  if (nuc_iter != nuc_end)
    data[0] = nuc_iter->second;
}


std::vector<pyne::gamma_struct> pyne::gamma_data;

void pyne::_load_gamma_data()
{

  // Loads the gamma table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(gamma_struct));
  status = H5Tinsert(desc, "energy", HOFFSET(gamma_struct, energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "energy_err", HOFFSET(gamma_struct, energy_err),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity", HOFFSET(gamma_struct, photon_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "photon_intensity_err", HOFFSET(gamma_struct, photon_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity", HOFFSET(gamma_struct, conv_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "conv_intensity_err", HOFFSET(gamma_struct, conv_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity", HOFFSET(gamma_struct, total_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "total_intensity_err", HOFFSET(gamma_struct, total_intensity_err), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(gamma_struct, from_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(gamma_struct, to_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "parent_nuc", HOFFSET(gamma_struct, parent_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(gamma_struct, k_conv_e), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(gamma_struct, l_conv_e), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(gamma_struct, m_conv_e), H5T_NATIVE_DOUBLE);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t gamma_set = H5Dopen2(nuc_data_h5, "/decay/gammas", H5P_DEFAULT);
  hid_t gamma_space = H5Dget_space(gamma_set);
  int gamma_length = H5Sget_simple_extent_npoints(gamma_space);

  // Read in the data
  gamma_struct * gamma_array = new gamma_struct[gamma_length];
  status = H5Dread(gamma_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, gamma_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(gamma_set);
  status = H5Fclose(nuc_data_h5);
  gamma_data.resize(gamma_length);

  for (int i = 0; i < gamma_length; ++i) {
    gamma_data[i] = gamma_array[i];
  }

}


int pyne::gamma_data_byparent(int nuc, gamma_struct *data){
  int count = 0;
  if (gamma_data.empty()) {
    _load_gamma_data();
    count = gamma_data_byparent(nuc, data);
    return count;
  }
  for (int i = 0; i < gamma_data.size(); ++i) {
    if (gamma_data[i].parent_nuc == nuc) ++count;
  }
  data = (gamma_struct*) malloc(sizeof(gamma_struct)*count);
  int n = 0;
  for (int i = 0; i < gamma_data.size(); ++i) {
    if (gamma_data[i].parent_nuc == nuc) {
      data[n] = gamma_data[i];
      ++n;
    }
  }
  return count;
}

int pyne::gamma_data_byen(double en, double pm, gamma_struct *data){
  int count = 0;
  if (gamma_data.empty()) {
    _load_gamma_data();
    count = gamma_data_byen(en, pm, data);
    return count;
  }
  for (int i = 0; i < gamma_data.size(); ++i) {
    if ((en + pm > gamma_data[i].energy) && (gamma_data[i].energy > en - pm)) ++count;
  }
  data = (gamma_struct*) malloc(sizeof(gamma_struct)*count);
  int n = 0;
  for (int i = 0; i < gamma_data.size(); ++i) {
    if ((en + pm > gamma_data[i].energy) && (gamma_data[i].energy > en - pm)) {
      data[n] = gamma_data[i];
      ++n;
    }
  }
  return count;
}


std::vector<pyne::alpha_struct> pyne::alpha_data;

void pyne::_load_alpha_data()
{

  // Loads the alpha table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(alpha_struct));
  status = H5Tinsert(desc, "energy", HOFFSET(alpha_struct, energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(alpha_struct, intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(alpha_struct, from_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(alpha_struct, to_nuc), H5T_NATIVE_INT);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t alpha_set = H5Dopen2(nuc_data_h5, "/decay/alphas", H5P_DEFAULT);
  hid_t alpha_space = H5Dget_space(alpha_set);
  int alpha_length = H5Sget_simple_extent_npoints(alpha_space);

  // Read in the data
  alpha_struct * alpha_array = new alpha_struct[alpha_length];
  status = H5Dread(alpha_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, alpha_array);

  // close the nuc_data library, before doing anything stupid
  status = H5Dclose(alpha_set);
  status = H5Fclose(nuc_data_h5);
  alpha_data.resize(alpha_length);

  for (int i = 0; i < alpha_length; ++i) {
    alpha_data[i] = alpha_array[i];
  }

}

int pyne::alpha_data_byparent(int nuc, alpha_struct *data){
  int count = 0;
  if (alpha_data.empty()) {
    _load_alpha_data();
    count = alpha_data_byparent(nuc, data);
    return count;
  }
  for (int i = 0; i < alpha_data.size(); ++i) {
    if (alpha_data[i].from_nuc == nuc) ++count;
  }
  data = (alpha_struct*) malloc(sizeof(alpha_struct)*count);
  int n = 0;
  for (int i = 0; i < alpha_data.size(); ++i) {
    if (alpha_data[i].from_nuc == nuc) {
      data[n] = alpha_data[i];
      ++n;
    }
  }
  return count;
}

int pyne::alpha_data_byen(double en, double pm, alpha_struct *data){
  int count = 0;
  if (alpha_data.empty()) {
    _load_alpha_data();
    count = alpha_data_byen(en, pm, data);
    return count;
  }
  for (int i = 0; i < alpha_data.size(); ++i) {
    if ((en + pm > alpha_data[i].energy) && (alpha_data[i].energy > en - pm)) ++count;
  }
  data = (alpha_struct*) malloc(sizeof(alpha_struct)*count);
  int n = 0;
  for (int i = 0; i < alpha_data.size(); ++i) {
    if ((en + pm > alpha_data[i].energy) && (alpha_data[i].energy > en - pm)) {
      data[n] = alpha_data[i];
      ++n;
    }
  }
  return count;
}


std::vector<pyne::beta_struct> pyne::beta_data;

void pyne::_load_beta_data()
{

  // Loads the beta table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(beta_struct));
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(beta_struct, endpoint_energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(beta_struct, avg_energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "intensity", HOFFSET(beta_struct, intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(beta_struct, from_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(beta_struct, to_nuc), H5T_NATIVE_INT);


  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t beta_set = H5Dopen2(nuc_data_h5, "/decay/betas", H5P_DEFAULT);
  hid_t beta_space = H5Dget_space(beta_set);
  int beta_length = H5Sget_simple_extent_npoints(beta_space);

  // Read in the data
  beta_struct * beta_array = new beta_struct[beta_length];
  status = H5Dread(beta_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, beta_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(beta_set);
  status = H5Fclose(nuc_data_h5);
  beta_data.resize(beta_length);

  for (int i = 0; i < beta_length; ++i) {
    beta_data[i] = beta_array[i];
  }

}


int pyne::beta_data_byparent(int nuc, beta_struct *data){
  int count = 0;
  if (beta_data.empty()) {
    _load_beta_data();
    count = beta_data_byparent(nuc, data);
    return count;
  }
  for (int i = 0; i < beta_data.size(); ++i) {
    if (beta_data[i].from_nuc == nuc) ++count;
  }
  data = (beta_struct*) malloc(sizeof(beta_struct)*count);
  int n = 0;
  for (int i = 0; i < beta_data.size(); ++i) {
    if (beta_data[i].from_nuc == nuc) {
      data[n] = beta_data[i];
      ++n;
    }
  }
  return count;
}


std::vector<pyne::ecbp_struct> pyne::ecbp_data;

void pyne::_load_ecbp_data()
{

  // Loads the ecbp table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(ecbp_struct));
  status = H5Tinsert(desc, "endpoint_energy", HOFFSET(ecbp_struct, endpoint_energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "avg_energy", HOFFSET(ecbp_struct, avg_energy),
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "beta_plus_intensity", HOFFSET(ecbp_struct, beta_plus_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "ec_intensity", HOFFSET(ecbp_struct, ec_intensity), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "from_nuc", HOFFSET(ecbp_struct, from_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(ecbp_struct, to_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "k_conv_e", HOFFSET(gamma_struct, k_conv_e), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "l_conv_e", HOFFSET(gamma_struct, l_conv_e), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "m_conv_e", HOFFSET(gamma_struct, m_conv_e), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t ecbp_set = H5Dopen2(nuc_data_h5, "/decay/ecbp", H5P_DEFAULT);
  hid_t ecbp_space = H5Dget_space(ecbp_set);
  int ecbp_length = H5Sget_simple_extent_npoints(ecbp_space);

  // Read in the data
  ecbp_struct * ecbp_array = new ecbp_struct[ecbp_length];
  status = H5Dread(ecbp_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecbp_array);

  // close the nuc_data library, before doing anythng stupid
  status = H5Dclose(ecbp_set);
  status = H5Fclose(nuc_data_h5);
  ecbp_data.resize(ecbp_length);

  for (int i = 0; i < ecbp_length; ++i) {
    ecbp_data[i] = ecbp_array[i];
  }

}


int pyne::ecbp_data_byparent(int nuc, ecbp_struct *data){
  int count = 0;
  if (ecbp_data.empty()) {
    _load_ecbp_data();
    count = ecbp_data_byparent(nuc, data);
    return count;
  }
  for (int i = 0; i < ecbp_data.size(); ++i) {
    if (ecbp_data[i].from_nuc == nuc) ++count;
  }
  data = (ecbp_struct*) malloc(sizeof(ecbp_struct)*count);
  int n = 0;
  for (int i = 0; i < ecbp_data.size(); ++i) {
    if (ecbp_data[i].from_nuc == nuc) {
      data[n] = ecbp_data[i];
      ++n;
    }
  }
  return count;
}
