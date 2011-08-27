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

  // Ok now that we have the array of stucts, put it in a the map
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
