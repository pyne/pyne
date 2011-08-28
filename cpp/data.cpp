// Implements basic nuclear data functions.

#include "data.h"


/****************************/
/*** nuc_weight Functions ***/
/****************************/
std::map<int, double> pyne::nuc_weight_map = std::map<int, double>();

void pyne::_load_nuc_weight_map()
{
  // Loads the importnat parts of atomic_wight table into nuc_weight_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool isH5 = H5::H5File::isHdf5(pyne::NUC_DATA_PATH);
  if (!isH5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);

  // Get the HDF5 compound type (table) description
  H5::CompType atomic_weight_desc(sizeof(atomic_weight_struct));
  atomic_weight_desc.insertMember("nuc_name", HOFFSET(atomic_weight_struct, nuc_name), H5::StrType(0, 6));
  atomic_weight_desc.insertMember("nuc_zz",   HOFFSET(atomic_weight_struct, nuc_zz),   H5::PredType::NATIVE_INT);
  atomic_weight_desc.insertMember("mass",     HOFFSET(atomic_weight_struct, mass),     H5::PredType::NATIVE_DOUBLE);
  atomic_weight_desc.insertMember("error",    HOFFSET(atomic_weight_struct, error),    H5::PredType::NATIVE_DOUBLE);
  atomic_weight_desc.insertMember("abund",    HOFFSET(atomic_weight_struct, abund),    H5::PredType::NATIVE_DOUBLE);

  // Open the HDF5 file
  H5::H5File nuc_data_h5 (pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY);

  // Open the data set
  H5::DataSet atomic_weight_set = nuc_data_h5.openDataSet("/atomic_weight");
  H5::DataSpace atomic_weight_space = atomic_weight_set.getSpace();
  int atomic_weight_length = atomic_weight_space.getSimpleExtentNpoints();

  // Read in the data
  atomic_weight_struct * atomic_weight_array = new atomic_weight_struct[atomic_weight_length];
  atomic_weight_set.read(atomic_weight_array, atomic_weight_desc);

  // close the nuc_data library, before doing anythng stupid
  nuc_data_h5.close();

  // Ok now that we have the array of stucts, put it in the map
  for(int n = 0; n < atomic_weight_length; n++)
    nuc_weight_map[atomic_weight_array[n].nuc_zz] = atomic_weight_array[n].mass;
};


double pyne::nuc_weight(int nuc)
{
  // Find the nuclide;s weight in AMU
  std::map<int, double>::iterator nuc_iter, nuc_end;

  nuc_iter = nuc_weight_map.find(nuc);
  nuc_end = nuc_weight_map.end();

  // First check if we already have the nuc weight in the map
  if (nuc_iter != nuc_end)
    return (*nuc_iter).second;

  // Next, fill up the map with values from the 
  // nuc_data.h5, if the map is empty.
  if (nuc_weight_map.empty())
  {
    // Don't fail if we can't load the library
    try
    {
      _load_nuc_weight_map();
      return nuc_weight(nuc);
    }
    catch(...){};
  };

  double aw;
  int nuc_zz = nucname::zzaaam(nuc);

  // If in an excited state, return the ground
  // state weight...not strictly true, but good guess.
  if (0 < nuc_zz%10)
  {
    aw = nuc_weight((nuc_zz/10)*10);
    nuc_weight_map[nuc] = aw;
    return aw;
  };

  // Finally, if none of these work, 
  // take a best guess based on the 
  // aaa number.
  aw = (double) ((nuc_zz/10)%1000);
  nuc_weight_map[nuc] = aw;
  return aw;
};


double pyne::nuc_weight(char * nuc)
{
  int nuc_zz = nucname::zzaaam(nuc);
  return nuc_weight(nuc_zz);
};


double pyne::nuc_weight(std::string nuc)
{
  int nuc_zz = nucname::zzaaam(nuc);
  return nuc_weight(nuc_zz);
};








/***********************************/
/*** scattering length functions ***/
/***********************************/
std::map<int, extra_types::complex_t> pyne::b_coherent_map = std::map<int, extra_types::complex_t>();
std::map<int, extra_types::complex_t> pyne::b_incoherent_map = std::map<int, extra_types::complex_t>();
std::map<int, double> pyne::b_map = std::map<int, double>();


void pyne::_load_scattering_lengths()
{
  // Loads the importnat parts of atomic_wight table into nuc_weight_map

  //Check to see if the file is in HDF5 format.
  if (!pyne::file_exists(pyne::NUC_DATA_PATH))
    throw pyne::FileNotFound(pyne::NUC_DATA_PATH);

  bool isH5 = H5::H5File::isHdf5(pyne::NUC_DATA_PATH);
  if (!isH5)
    throw h5wrap::FileNotHDF5(pyne::NUC_DATA_PATH);


  // Get the HDF5 compound type (table) description
  H5::CompType scat_len_desc(sizeof(scattering_lengths_struct));
  scat_len_desc.insertMember("nuc_name", HOFFSET(scattering_lengths_struct, nuc_name), H5::StrType(0, 6));
  scat_len_desc.insertMember("nuc_zz",   HOFFSET(scattering_lengths_struct, nuc_zz),   H5::PredType::NATIVE_INT);
  scat_len_desc.insertMember("b_coherent", HOFFSET(scattering_lengths_struct, b_coherent), h5wrap::PYTABLES_COMPLEX128);
  scat_len_desc.insertMember("b_incoherent", HOFFSET(scattering_lengths_struct, b_incoherent), h5wrap::PYTABLES_COMPLEX128);
  scat_len_desc.insertMember("xs_coherent", HOFFSET(scattering_lengths_struct, xs_coherent), H5::PredType::NATIVE_DOUBLE);
  scat_len_desc.insertMember("xs_incoherent", HOFFSET(scattering_lengths_struct, xs_incoherent), H5::PredType::NATIVE_DOUBLE);
  scat_len_desc.insertMember("xs", HOFFSET(scattering_lengths_struct, xs), H5::PredType::NATIVE_DOUBLE);

  // Open the HDF5 file
  H5::H5File nuc_data_h5 (pyne::NUC_DATA_PATH.c_str(), H5F_ACC_RDONLY);

  // Open the data set
  H5::DataSet scat_len_set = nuc_data_h5.openDataSet("/neutron/scattering_lengths");
  H5::DataSpace scat_len_space = scat_len_set.getSpace();
  int scat_len_length = scat_len_space.getSimpleExtentNpoints();

  // Read in the data
  scattering_lengths_struct * scat_len_array = new scattering_lengths_struct[scat_len_length];
  scat_len_set.read(scat_len_array, scat_len_desc);

  // close the nuc_data library, before doing anythng stupid
  nuc_data_h5.close();

  // Ok now that we have the array of stucts, put it in the maps
  extra_types::complex_t bc, bi;
  for(int n = 0; n < scat_len_length; n++)
  {
    b_coherent_map[scat_len_array[n].nuc_zz] = scat_len_array[n].b_coherent;
    b_incoherent_map[scat_len_array[n].nuc_zz] = scat_len_array[n].b_incoherent;
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
  int nuc_zz = nucname::zzaaam(nuc);
  int znum = nuc_zz/10000;
  int anum = (nuc_zz/10)%1000;

  // Try to find a nuclide with matching A-number
  nuc_iter = b_coherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (anum == (((*nuc_iter).first)/10)%1000)
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
    if (znum == ((*nuc_iter).first)/10000)
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
  int nuc_zz = nucname::zzaaam(nuc);
  return b_coherent(nuc_zz);
};


extra_types::complex_t pyne::b_coherent(std::string nuc)
{
  int nuc_zz = nucname::zzaaam(nuc);
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
  int nuc_zz = nucname::zzaaam(nuc);
  int znum = nuc_zz/10000;
  int anum = (nuc_zz/10)%1000;

  // Try to find a nuclide with matching A-number
  nuc_iter = b_incoherent_map.begin();
  while (nuc_iter != nuc_end)
  {
    if (anum == (((*nuc_iter).first)/10)%1000)
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
    if (znum == ((*nuc_iter).first)/10000)
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
  int nuc_zz = nucname::zzaaam(nuc);
  return b_incoherent(nuc_zz);
};


extra_types::complex_t pyne::b_incoherent(std::string nuc)
{
  int nuc_zz = nucname::zzaaam(nuc);
  return b_incoherent(nuc_zz);
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
  int nuc_zz = nucname::zzaaam(nuc);
  return b(nuc_zz);
};


double pyne::b(std::string nuc)
{
  int nuc_zz = nucname::zzaaam(nuc);
  return b(nuc_zz);
};
