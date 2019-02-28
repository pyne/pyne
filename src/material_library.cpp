#include <unistd.h>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

namespace pyne {
// Empty Constructor
MaterialLibrary::MaterialLibrary(){};

// Default constructor
MaterialLibrary::MaterialLibrary(const std::string& file,
                                 const std::string& datapath) {
  if (!check_file_exists(file)) {
    throw std::runtime_error("File " + file +
                             " not found or no read permission");
    if (!hdf5_path_exists(file, datapath)) {
      throw std::runtime_error("The datapath, " + datapath + ", in "
                               << file + " is empty.");
    }

    // load materials
    from_hdf5(file);
  };

  // Destructor
  MaterialLibrary::~MaterialLibrary(){};

  // convert convert a filename into path+filename (for pyne)
  std::string MaterialLibrary::get_full_filepath(char* filename) {
    std::string file(filename);
    return MaterialLibrary::get_full_filepath(file);
  }

  // convert convert a filename into path+filename (for pyne)
  const std::string& MaterialLibrary::get_full_filepath(
      const std::string& filename) {
    // remove all extra whitespace
    filename.erase(std::remove(filename.begin(), filename.end(), ' '),
                   filename.end());
    // use stdlib call
    const char* full_filepath = realpath(filename.c_str(), NULL);
    return std::string(full_filepath);
  }

  // see if file exists
  bool MaterialLibrary::check_file_exists(const std::string& filename) {
    // from http://stackoverflow.com/questions/12774207/
    // fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
    std::ifstream infile(filename.c_str());
    return infile.good();
  }

  // Append Material to the library from hdf5 file
  void MaterialLibrary::from_hdf5(const std::string& filename,
                                  const std::string& datapath, int protocol) {
    const char* data_path = datapath.c_str();

    if (!hdf5_path_exists(filename, data_path)) return;

    int num_materials = get_length_of_table(filename, datapath);

    for (int i = 0; i < num_materials; i++) {
      pyne::Material mat;  // from file
      // read new mat
      mat.from_hdf5(filename, datapath, i, protocol);

      append_to_nuclist(mat);
      material_library[mat.metadata["name"].asString()] = mat;
      matlist.insert(mat.metadata["name"].asString());
    }
  }

  void MaterialLibrary::add_material(pyne::Material mat) {
    std::string mat_name = mat.metadata["name"].asString();
    auto ins = material_library.insert(std::make_pair(mat_name, mat));

    if (ins.second) {
      append_to_nuclist(mat);
      matlist.insert(mat.metadata["name"].asString());
    }
  }

  void MaterialLibrary::del_material(pyne::Material mat) {
    std::string mat_name = mat.metadata["name"].asString();
    del_material(mat_name);
  }

  void MaterialLibrary::del_material(const std::string& mat_name) {
    material_library.erase(mat_name);
    matlist.erase(mat_name);
  }

  pyne::Material MaterialLibrary::get_material(const std::string& mat_name) {
    auto it = material_library.find(mat_name);
    if (it != material_library.end()) {
      return it->second;
    } else {
      return pyne::Material();
    }
  }

  void MaterialLibrary::write_hdf5(const std::string& filename,
                                   const std::string& datapath,
                                   const std::string& nucpath, int chunksize) {
    // A large part of this is inspired/taken from by material.cpp...
    // Turn off annoying HDF5 errors
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    // Set file access properties so it closes cleanly
    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
    // Create new/open datafile.
    hid_t db;
    if (pyne::file_exists(filename)) {
      bool ish5 = H5Fis_hdf5(filename.c_str());
      if (!ish5) throw h5wrap::FileNotHDF5(filename);
      db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
    } else
      db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

    pyne::mat_map material_library_tmp = material_library;
    pyne::nuc_set nuclist_tmp = nuclist;
    // We are using material_library and nuclist copy to allow adding exisitng
    // nuclide/materials from the file

    std::vector<int> nuclides;
    nuclides.clear();
    std::copy(nuclist_tmp.begin(), nuclist_tmp.end(),
              inserter(nuclides, nuclides.begin()));

    int nuclide_size;
    hsize_t nuclide_dims[1];
    nuclide_size = nuclides.size();

    // Create the data if it doesn't exist
    int nuclide_data[nuclide_size];
    int n = 0;
    for (auto nuclide : nuclides) {
      nuclide_data[n] = (*nuclide);
      n++;
    }
    nuclide_dims[0] = nuclide_size;

    // Write Nuclist in the hdf5 file
    // Not sure if thsi will work with an existing nuclist in the file
    hid_t nuclide_space = H5Screate_simple(1, nuclide_dims, NULL);
    hid_t nuclide_set =
        H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuclide_space,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(nuclide_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             nuclide_data);
    H5Fflush(db, H5F_SCOPE_GLOBAL);
    // Close out the Dataset
    H5Tclose(db);

    // Write the Materials in the file
    for (auto material : material_library_tmp) {
      (material->second).write_hdf5(filename, datapath, nucpath);
    }
  }

  void MaterialLibrary::append_to_nuclist(pyne::Material mat) {
    pyne::comp_map mat_comp = mat.comp;
    for (auto nuclide : mat_comp) {
      nuclist.insert(nuclide->first);
    }
  }

  // see if path exists before we go on
  bool MaterialLibrary::hdf5_path_exists(const std::string& filename,
                                         const std::string& datapath) {
    // Turn off annoying HDF5 errors
    herr_t status;
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    // Set file access properties so it closes cleanly
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

    hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

    bool datapath_exists = h5wrap::path_exists(db, datapath.c_str());

    status = H5Eclear(H5E_DEFAULT);

    // Close the database
    status = H5Fclose(db);

    return datapath_exists;
  }

  int MaterialLibrary::get_length_of_table(const std::string& filename,
                                           const std::string& datapath) {
    // Turn off annoying HDF5 errors
    herr_t status;
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    // Set file access properties so it closes cleanly
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
}  // end of pyne namespace
