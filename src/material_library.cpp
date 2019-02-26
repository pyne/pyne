#include <unistd.h>
#include <iostream>


#ifndef PYNE_IS_AMALGAMATED
  #include "transmuters.h"
  #include "material.h"
#endif

// Empty Constructor
MaterialLibrary::MaterialLibrary() {
};

// Default constructor
MaterialLibrary::MaterialLibrary(string filei, string datapath) {

  if (!check_file_exists(file)) {
    std::cerr << "The file " << file << " does not exist or is read protected" << std::endl;
    exit(1);
  }
  if (!hdf5_path_exists(file, datapath)) {
    std::cerr << "The datapath " << datapath << " in " << file << " empty." << std::endl;
    exit(1);
  }

  // load materials
  material_library = load_materials(full_filepath);

};


// Destructor
MaterialLibrary::~MaterialLibrary() {
};

// convert convert a filename into path+filename (for pyne)
std::string MaterialLibrary::get_full_filepath(char* filename) {
  std::string file(filename);
  return MaterialLibrary::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string MaterialLibrary::get_full_filepath(std::string filename) {
  // remove all extra whitespace
  filename.erase(std::remove(filename.begin(), filename.end(), ' '), filename.end());
  // use stdlib call
  const char* full_filepath = realpath(filename.c_str(), NULL);
  return std::string(full_filepath);
}

// see if file exists
bool MaterialLibrary::check_file_exists(std::string filename) {
  // from http://stackoverflow.com/questions/12774207/
  // fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  std::ifstream infile(filename.c_str());
  return infile.good();
}

// loads all materials into map

std::map<std::string, pyne::Material> MaterialLibrary::load_pyne_materials(std::string filename, std::string datapath) {
  std::map<std::string, pyne::Material> library; // material library

  const char* data_path = datapath.c_str();

  if (!hdf5_path_exists(filename, data_path))
    return library;

  num_materials = get_length_of_table(filename, datapath);

  for (int i = 0 ; i < num_materials ; i++) {
    pyne::Material mat; // from file
    mat.from_hdf5(filename, datapath, i);
    // renumber material number by position in the library
    mat.metadata["mat_number"] = i + 1;
    library[mat.metadata["name"].asString()] = mat;
  }

  return library;
}


// see if path exists before we go on
bool MaterialLibrary::hdf5_path_exists(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath.c_str());

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return datapath_exists;
}

int MaterialLibrary::get_length_of_table(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);

  hsize_t arr_dims[1];
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return arr_dims[0];
}
