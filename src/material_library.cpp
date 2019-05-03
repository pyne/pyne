#include <unistd.h>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

// Empty Constructor
pyne::MaterialLibrary::MaterialLibrary(){
  last_mat_number = 0;
};

// Default constructor
pyne::MaterialLibrary::MaterialLibrary(const std::string& file,
                                       const std::string& datapath) {
  if (!pyne::file_exists(file)) {
    throw std::runtime_error("File " + file +
                             " not found or no read permission");
  }
  if (!hdf5_path_exists(file, datapath)) {
    throw std::runtime_error("The datapath, " + datapath + ", in " + file +
                             " is empty.");
  }
  last_mat_number = 0;
  // load materials
  from_hdf5(file, datapath);
};

// Destructor
pyne::MaterialLibrary::~MaterialLibrary(){

 // for (auto element : material_library) {
 //   delete element.second;
 //   material_library.erase(element.first);
 // }


};

void pyne::MaterialLibrary::from_hdf5(char* fname, char* dpath, char* npath, int protocol) {
  std::string filename(fname);
  std::string datapath(dpath);
  std::string nucpath(npath);
  from_hdf5(filename, datapath, npath, protocol);
}
// Append Material to the library from hdf5 file
void pyne::MaterialLibrary::from_hdf5(const std::string& filename,
                                      const std::string& datapath,
                                      const std::string& nucpath,
                                      int protocol) {
  if (!hdf5_path_exists(filename, datapath)) return;

  int file_num_materials = get_length_of_table(filename, datapath);
  int library_length = material_library.size();
  for (int i = 0; i < file_num_materials; i++) {
    pyne::Material* mat = new pyne::Material();  // from file
    // read new mat
    if( nucpath == ""){
      mat->from_hdf5(filename, datapath, i, protocol);
    }  
    else {
       mat->from_hdf5(filename, datapath, nucpath, i, protocol);
    }
    // if exists, get the material name from metadata make one instead
    std::string mat_name;
    std::pair<pyne::matname_set::iterator, bool> mat_insert;
    if (mat->metadata.isMember("name")) {
      mat_name = mat->metadata["name"].asString();
      mat_insert = matlist.insert(mat->metadata["name"].asString());
    } else {
      // form a mat name as 'mX'
      int mat_number = library_length + i;
      mat_name = "m" + std::to_string(mat_number);
      mat->metadata["name"] = mat_name;
      mat_insert = matlist.insert(mat->metadata["name"].asString());

      // if 'mX' name is already in the Library rename the material
      while (!mat_insert.second) {
        mat_number += file_num_materials;
        mat_name = "m" + std::to_string(library_length + i);
        mat->metadata["name"] = mat_name;
        mat_insert = matlist.insert(mat->metadata["name"].asString());
      }
    }
    if (mat_insert.second) {
      append_to_nuclist(*mat);
      material_library[mat->metadata["name"].asString()] = mat;
      material_library[mat->metadata["name"].asString()]->metadata["mat_number"] = last_mat_number++;
    }
  }
}

void pyne::MaterialLibrary::merge(pyne::MaterialLibrary mat_lib) {
  pyne::matname_set mats_to_add = mat_lib.get_matlist();
  for (auto it = mats_to_add.begin(); it != mats_to_add.end(); it++) {
    pyne::Material mat = Material(mat_lib.get_material(*it));
    (*this).add_material(mat);
  }
}

void pyne::MaterialLibrary::load_json(Json::Value json) {
  
  Json::Value::Members keys = json.getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  mat_map mat_lib_tmp; // required to maintain order in mat number
  for (; ikey != ikey_end; ++ikey) {
    pyne::Material mat = pyne::Material();
    mat.load_json(json[*ikey]);
    mat_lib_tmp[mat.metadata["name"].asString()] = new Material(mat);
  }
  for (auto mat_element : mat_lib_tmp) {
    (*this).add_material(*mat_element.second);
  }
}

Json::Value pyne::MaterialLibrary::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);
  for( auto mat : material_library)
    json[mat.first] = (mat.second)->dump_json();
  return json;
}

void pyne::MaterialLibrary::from_json(char* fname) {
  std::string filename(fname);
  from_json(filename);

}
void pyne::MaterialLibrary::from_json(const std::string& filename) {
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f (filename.c_str(), std::ios::in | std::ios::binary);
  f.seekg(0, std::ios::end);
  s.resize(f.tellg());
  f.seekg(0, std::ios::beg);
  f.read(&s[0], s.size());
  f.close();
  Json::Reader reader;
  Json::Value json;
  reader.parse(s, json);
  load_json(json);
}

void pyne::MaterialLibrary::write_json(char * filename) {
  std::string fname (filename);
  write_json(fname);
}

void pyne::MaterialLibrary::write_json(const std::string& filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}

void pyne::MaterialLibrary::add_material(pyne::Material mat) {
  std::string mat_name = mat.metadata["name"].asString();
  pyne::Material* mat_p = new pyne::Material(mat);
  auto inst = material_library.insert(std::make_pair(mat_name, mat_p));

  if (inst.second) {
    append_to_nuclist(*mat_p);
    matlist.insert(mat_p->metadata["name"].asString());
    material_library[mat_p->metadata["name"].asString()]->metadata["mat_number"] = last_mat_number++;
  }
}

void pyne::MaterialLibrary::add_material(pyne::Material mat, char* mat_name) {
  std::string material_name(mat_name);
  add_material(mat, material_name);
}

void pyne::MaterialLibrary::add_material(pyne::Material mat,
                                         const std::string& mat_name) {
  mat.metadata["name"] = mat_name;
  add_material(mat);
}

void pyne::MaterialLibrary::del_material(pyne::Material mat) {
  if (mat.metadata.isMember("name")) {
    std::string mat_name = mat.metadata["name"].asString();
    del_material(mat_name);
  }
}

void pyne::MaterialLibrary::del_material(const std::string& mat_name) {
  material_library.erase(mat_name);
  matlist.erase(mat_name);
}

pyne::Material pyne::MaterialLibrary::get_material(
    const std::string& mat_name) {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return *(it->second);
  } else {
    return pyne::Material();
  }
}

void pyne::MaterialLibrary::write_hdf5(char* fname, char* dpath, char* npath){
  std::string filename(fname);
  std::string datapath(dpath);
  std::string nucpath(npath);
  write_hdf5(filename, datapath, nucpath);
}

void pyne::MaterialLibrary::write_hdf5(const std::string& filename,
                                       const std::string& datapath,
                                       const std::string& nucpath) {
  // A large part of this is inspiwrite_hdf5red/taken from by material.cpp...
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
    nuclide_data[n] = (nuclide);
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
  for (auto element : material_library_tmp) {
    (element.second)->write_hdf5(filename, datapath, nucpath);
  }
}

void pyne::MaterialLibrary::append_to_nuclist(pyne::Material mat) {
  pyne::comp_map mat_comp = mat.comp;
  for (auto nuclide : mat_comp) {
    nuclist.insert(nuclide.first);
  }
}

// see if path exists before we go on
bool pyne::MaterialLibrary::hdf5_path_exists(const std::string& filename,
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

int pyne::MaterialLibrary::get_length_of_table(const std::string& filename,
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
