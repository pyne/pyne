#include <unistd.h>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

// Empty Constructor
pyne::MaterialLibrary::MaterialLibrary(){};

// Default constructor
pyne::MaterialLibrary::MaterialLibrary(const std::string& file,
                                       const std::string& datapath,
                                       const std::string& nucpath) {
  if (!pyne::file_exists(file)) {
    throw std::runtime_error("File " + file +
                             " not found or no read permission");
  }
  if (!hdf5_path_exists(file, datapath)) {
    throw std::runtime_error("The datapath, " + datapath + ", in " + file +
                             " is empty.");
  }
  // load materials
  from_hdf5(file, datapath, nucpath);
};

void pyne::MaterialLibrary::from_hdf5(char* fname, char* dpath, char* npath,
                                      int protocol) {
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
  if (!hdf5_path_exists(filename, datapath)){
    throw std::runtime_error("The datapath, " + datapath + ", in " + filename +
                             " is empty.");
  }

  int file_num_materials = get_length_of_table(filename, datapath);
  int library_length = material_library.size();
  for (int i = 0; i < file_num_materials; i++) {
    pyne::Material mat = pyne::Material();  
    // Get material from the hdf5 file
    if (nucpath == "") {
      mat.from_hdf5(filename, datapath, i, protocol);
    } else {
      mat.from_hdf5(filename, datapath, nucpath, i, protocol);
    }
    (*this).add_material(mat);
  }
}

void pyne::MaterialLibrary::merge(pyne::MaterialLibrary mat_lib) {
  pyne::matname_set mats_to_add = mat_lib.get_keylist();
  for (auto it = mats_to_add.begin(); it != mats_to_add.end(); it++) {
    pyne::Material mat = Material(mat_lib.get_material(*it));
    (*this).add_material(*it, mat);
  }
}

void pyne::MaterialLibrary::load_json(Json::Value json) {
  Json::Value::Members keys = json.getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  for (; ikey != ikey_end; ++ikey) {
    pyne::Material mat = pyne::Material();
    mat.load_json(json[*ikey]);
    (*this).add_material(*ikey, Material(mat));
  }
}

Json::Value pyne::MaterialLibrary::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);

  for (auto name : name_order) {
    json[name] = material_library[name]->dump_json();
  }
  return json;
}

void pyne::MaterialLibrary::from_json(char* fname) {
  std::string filename(fname);
  from_json(filename);
}
void pyne::MaterialLibrary::from_json(const std::string& filename) {
  if (!pyne::file_exists(filename)) throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
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

void pyne::MaterialLibrary::write_json(char* filename) {
  std::string fname(filename);
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
  // if exists, get the material name from metadata make one instead
  std::string mat_name;
  int mat_number = -1;
  std::set<int>::iterator mat_numb_it;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue ){
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
    }
  } 
  pyne::matname_set::iterator key_it;
  if (mat.metadata.isMember("name")) {
    
    if (Json::intValue <= mat.metadata["name"].type() &&
        mat.metadata["name"].type() <= Json::realValue ){
      mat_name = std::to_string(mat.metadata["name"].asInt());
      mat.metadata["name"] = mat_name; 
    } else {
      mat_name = mat.metadata["name"].asString();
    }
    
  } else {
    if (mat_number == -1) {
      mat_number = keylist.size(); //set a temp mat_number to form mat name 
    }
    mat_name = "_" + std::to_string(mat_number);
    mat.metadata["name"] = mat_name;
  }
  
  add_material(mat_name, mat);  
}

void pyne::MaterialLibrary::add_material(char* key, pyne::Material mat) {
  std::string key_name(key);
  add_material(key_name, mat);
}

void pyne::MaterialLibrary::add_material(const std::string& key, pyne::Material mat) {
  std::string mat_name;
  int mat_number = -1;
  std::set<int>::iterator mat_numb_it;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue ){
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
      mat.metadata["mat_number"] = mat_number;
    }
    mat_numb_it = mat_number_set.find(mat_number);
    if (mat_numb_it != mat_number_set.end()) {
      std::string msg = "Material number ";
      msg += mat_number;
      msg += " is already in the library.";
      warning(msg);
    }
  }  
  if ( !mat.metadata.isMember("name")) {
    mat.metadata["name"] = key;
  } 

  append_to_nuclist(mat);
  if(mat_number > 0)
    mat_number_set.insert(mat_number);
  std::pair<std::set<std::string>::iterator, bool> key_insert;
  key_insert = keylist.insert(key);
  // if the key was not present in the MaterialLibrary add it, otherwise overwrite the existing material
  if( key_insert.second == true)
    name_order.push_back(key);
  material_library[key] = new Material(mat);
}

void pyne::MaterialLibrary::del_material(pyne::Material mat) {
  if (mat.metadata.isMember("name")) {
    std::string mat_name = mat.metadata["name"].asString();
    del_material(mat_name);
  }
}

void pyne::MaterialLibrary::del_material(const std::string& key) {
  material_library.erase(key);
  keylist.erase(key);
  for (auto name = name_order.begin(); name != name_order.end(); name++) {
    if (*name == key) {
      name_order.erase(name);
      return;
    }
  }
}

pyne::Material pyne::MaterialLibrary::get_material(
    const std::string& mat_name) const {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return *(it->second);
  } else {
    return pyne::Material();
  }
}

pyne::Material* pyne::MaterialLibrary::get_material_ptr(
    const std::string& mat_name) {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return it->second;
  } else {
    return new pyne::Material();
  }
}

void pyne::MaterialLibrary::replace(int num, pyne::Material mat) {
  if (!mat.metadata.isMember("name")) {
    mat.metadata["name"] = name_order[num];
  } else if (mat.metadata["name"] != name_order[num]) {
    material_library.erase(name_order[num]);
    name_order[num] = mat.metadata["name"].asString();
  }
  material_library[name_order[num]] = new pyne::Material(mat);
}

void pyne::MaterialLibrary::write_hdf5(char* fname, char* dpath, char* npath) {
  std::string filename(fname);
  std::string datapath(dpath);
  std::string nucpath(npath);
  write_hdf5(filename, datapath, nucpath);
}

void pyne::MaterialLibrary::write_hdf5(const std::string& filename,
                                       const std::string& datapath,
                                       const std::string& nucpath) {
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

  std::vector<int> nuclides;
  nuclides.clear();
  std::copy(nuclist.begin(), nuclist.end(),
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
  H5Fclose(db);

  // Write the Materials in the file
  for (auto name : name_order) {
    material_library[name]->write_hdf5(filename, datapath, nucpath);
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
