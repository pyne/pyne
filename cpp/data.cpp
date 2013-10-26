// Implements basic nuclear data functions.

#include "data.h"


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
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic_weight_struct));
  H5Tinsert(desc, "nuc",   HOFFSET(atomic_weight_struct, nuc),   H5T_NATIVE_INT);
  H5Tinsert(desc, "mass",  HOFFSET(atomic_weight_struct, mass),  H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "error", HOFFSET(atomic_weight_struct, error), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "abund", HOFFSET(atomic_weight_struct, abund), H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atomic_weight_set = H5Dopen2(nuc_data_h5, "/atomic_weight", H5P_DEFAULT);
  hid_t atomic_weight_space = H5Dget_space(atomic_weight_set);
  int atomic_weight_length = H5Sget_simple_extent_npoints(atomic_weight_space);

  // Read in the data
  atomic_weight_struct * atomic_weight_array = new atomic_weight_struct[atomic_weight_length];
  H5Dread(atomic_weight_set, desc, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_weight_array);

  // close the nuc_data library, before doing anythng stupid
  H5Dclose(atomic_weight_set);
  H5Fclose(nuc_data_h5);

  // Ok now that we have the array of stucts, put it in the map
  for(int n = 0; n < atomic_weight_length; n++){
    atomic_mass_map[atomic_weight_array[n].nuc] = atomic_weight_array[n].mass;
    natural_abund_map[atomic_weight_array[n].nuc] = atomic_weight_array[n].abund;
  }
};


double pyne::atomic_mass(int nuc)
{
  // Find the nuclide's weight in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = atomic_mass_map.find(nuc);
  nuc_end = atomic_mass_map.end();

  // First check if we already have the nuc weight in the map
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
  // state weight...not strictly true, but good guess.
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

  // First check if we already have the nuc weight in the map
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









/******************************/
/*** atomic decay functions ***/
/******************************/
std::map<int, double> pyne::half_life_map = std::map<int, double>();
std::map<int, double> pyne::decay_const_map = std::map<int, double>();
std::map<std::pair<int, int>, double> pyne::branch_ratio_map = \
                                              std::map<std::pair<int, int>, double>();
std::map<int, double> pyne::state_energy_map = std::map<int, double>();
std::map<int, std::set<int> > pyne::decay_children_map = \
                                                      std::map<int, std::set<int> >();


void pyne::_load_atomic_decay()
{
  // Loads the important parts of atomic_decay table into memory
  herr_t status;

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool ish5 = H5Fis_hdf5(pyne::NUC_DATA_PATH.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  hid_t desc = H5Tcreate(H5T_COMPOUND, sizeof(atomic_decay_struct));
  status = H5Tinsert(desc, "from_nuc", HOFFSET(atomic_decay_struct, from_nuc), 
                      H5T_NATIVE_INT);
  status = H5Tinsert(desc, "level", HOFFSET(atomic_decay_struct, level), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "to_nuc", HOFFSET(atomic_decay_struct, to_nuc), H5T_NATIVE_INT);
  status = H5Tinsert(desc, "half_life", HOFFSET(atomic_decay_struct, half_life), 
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "decay_const", HOFFSET(atomic_decay_struct, decay_const), 
                      H5T_NATIVE_DOUBLE);
  status = H5Tinsert(desc, "branch_ratio", HOFFSET(atomic_decay_struct, branch_ratio), 
                      H5T_NATIVE_DOUBLE);

  // Open the HDF5 file
  hid_t nuc_data_h5 = H5Fopen(pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open the data set
  hid_t atom_dec_set = H5Dopen2(nuc_data_h5, "/atomic_decay", H5P_DEFAULT);
  hid_t atom_dec_space = H5Dget_space(atom_dec_set);
  int atom_dec_length = H5Sget_simple_extent_npoints(atom_dec_space);

  // Read in the data
  atomic_decay_struct * atom_dec_array = new atomic_decay_struct[atom_dec_length];
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
};


//
// Half-life data
//

double pyne::half_life(int nuc)
{
  // Find the nuclide's half life in s
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = half_life_map.find(nuc);
  nuc_end = half_life_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (half_life_map.empty())
  {
    _load_atomic_decay();
    return half_life(nuc);
  };

  // Finally, if none of these work, 
  // assume the value is stable
  double hl = 1.0 / 0.0;
  half_life_map[nuc] = hl;
  return hl;
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
    _load_atomic_decay();
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
    _load_atomic_decay();
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
// Excitation state energy data
//

double pyne::state_energy(int nuc)
{
  // Find the nuclide's state energy in MeV
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = state_energy_map.find(nuc);
  nuc_end = state_energy_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (state_energy_map.empty())
  {
    _load_atomic_decay();
    return state_energy(nuc);
  };

  // Finally, if none of these work, 
  // assume the value is stable
  double se = 0.0;
  state_energy_map[nuc] = se;
  return se;
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
// Decay children data
//

std::set<int> pyne::decay_children(int nuc)
{
  // Find the nuclide's decay constant in 1/s
  std::map<int, std::set<int> >::iterator nuc_iter, nuc_end;

  nuc_iter = decay_children_map.find(nuc);
  nuc_end = decay_children_map.end();

  // First check if we already have the nuc in the map
  if (nuc_iter != nuc_end){
    return (*nuc_iter).second;
  };

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (decay_children_map.empty())
  {
    _load_atomic_decay();
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
