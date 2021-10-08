#ifndef _WIN32
#include <unistd.h>
#endif
#include <fstream>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

// Empty Constructor
pyne::MaterialLibrary::MaterialLibrary(){};

// Default constructor when loading from HDF5 file
pyne::MaterialLibrary::MaterialLibrary(const std::string& filename,
                                       const std::string& datapath) {
  // load materials

  // Check that the file is there
  if (!pyne::file_exists(filename)) throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename, datapath);
  else if (datapath != "") {
    std::string wrng =
        filename + " is not an hdf5 file but a datapath was prodided!";
    pyne::warning(wrng);
  }
  from_json(filename);
};

// Append Material to the library from hdf5 file
void pyne::MaterialLibrary::from_hdf5(const std::string& filename,
                                      const std::string& datapath) {
  if (!pyne::file_exists(filename)) {
    throw std::runtime_error("File " + filename +
                             " not found or no read permission");
  }
  std::string nucpath;
  std::string full_datapath = datapath;

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  // Get datapath status in the hdf5 file
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  hid_t matlib_group_id = db;
  H5O_info_t object_info;
  status = H5Oget_info_by_name(db, datapath.c_str(), &object_info, H5P_DEFAULT);

  // if datapath exist and is a dataset -> old hdf5 mat_lib format
  // else if datapath does not exist -> new hdf5 mat lib format
  // else don't know what to do with it: fail !
  if (status == 0 && object_info.type == H5O_TYPE_DATASET) {
    // Open dataset to retrieve nucclide list location form nucpath attribute
    hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    pyne::detect_nuclidelist(data_set, nucpath);
    H5Dclose(data_set);
  } else if (status != 0) {
    full_datapath = "/material_library" + datapath + "/composition";
    nucpath = "/material_library" + datapath + "/nuclidelist";
  } else {
    throw std::runtime_error(datapath + " exist but is not a dataset.");
  }

  int file_num_materials = get_length_of_table(db, full_datapath);

  for (int i = 0; i < file_num_materials; i++) {
    pyne::Material mat = pyne::Material();

    mat._load_comp_protocol1(db, full_datapath, nucpath, i);
    add_material(mat);
  }
  H5Fclose(db);
}

void pyne::MaterialLibrary::merge(const pyne::MaterialLibrary& mat_lib) {
  pyne::matname_set mats_to_add = mat_lib.get_keylist();
  for (auto it = mats_to_add.begin(); it != mats_to_add.end(); it++) {
    pyne::Material mat = Material(mat_lib.get_material(*it));
    (*this).add_material(*it, mat);
  }
}

void pyne::MaterialLibrary::merge(pyne::MaterialLibrary* mat_lib) {
  merge(*mat_lib);
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

  for (auto name : keylist) {
    json[name] = material_library[name]->dump_json();
  }
  return json;
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

void pyne::MaterialLibrary::write_json(const std::string& filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}

void pyne::MaterialLibrary::write_openmc(const std::string& filename) const
{
  std::ofstream f(filename.c_str());
  // write header
  f << "<?xml version=\"1.0\"?>\n";
  f << "<materials>\n";
  // write materials
  for (auto& mat : material_library) {
    f << mat.second->openmc("atom");
  }
  // write footer
  f << "</materials>";

  // close file
  f.close();
}

int pyne::MaterialLibrary::ensure_material_number(pyne::Material& mat) const {
  int mat_number = -1;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue) {
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
      mat.metadata["mat_number"] = mat_number;
    }
    std::set<int>::iterator mat_numb_it;
    mat_numb_it = mat_number_set.find(mat_number);
    if (mat_numb_it != mat_number_set.end()) {
      std::string msg = "Material number ";
      msg += std::to_string(mat_number);
      msg += " is already in the library.";
      warning(msg);
    }
  }
  if (mat_number == -1) {
    if (!mat_number_set.empty())
      mat_number = *mat_number_set.rbegin() + 1;
    else {
      mat_number = 1;
    }
    mat.metadata["mat_number"] = mat_number;
  }
  return mat_number;
}

std::string pyne::MaterialLibrary::ensure_material_name_and_number(
    pyne::Material& mat) const {
  std::string mat_name = "";
  int mat_number = ensure_material_number(mat);

  if (mat.metadata.isMember("name")) {
    if (Json::intValue <= mat.metadata["name"].type() &&
        mat.metadata["name"].type() <= Json::realValue) {
      mat_name = std::to_string(mat.metadata["name"].asInt());
      mat.metadata["name"] = mat_name;
    } else {
      mat_name = mat.metadata["name"].asString();
    }
  } else {
    mat_name = "_" + std::to_string(mat_number);
    mat.metadata["name"] = mat_name;
  }
  return mat_name;
}

void pyne::MaterialLibrary::add_material(pyne::Material mat) {
  // if exists, get the material name from metadata make one instead
  std::string mat_name = ensure_material_name_and_number(mat);
  add_material(mat_name, mat);
}

void pyne::MaterialLibrary::add_material(const std::string& key,
                                         const pyne::Material& mat) {
  pyne::shr_mat_ptr new_mat = std::make_shared<pyne::Material>(mat);

  int mat_number = ensure_material_number(*new_mat);
  if (!new_mat->metadata.isMember("name")) {
    new_mat->metadata["name"] = key;
  }
  append_to_nuclist(*new_mat);
  if (mat_number > 0) mat_number_set.insert(mat_number);
  keylist.insert(key);
  material_library[key] = new_mat;
}

void pyne::MaterialLibrary::del_material(const std::string& key) {
  material_library.erase(key);
  keylist.erase(key);
}

pyne::Material pyne::MaterialLibrary::get_material(
    const std::string& mat_name) const {
  return *(this->get_material_ptr(mat_name));
}

pyne::shr_mat_ptr pyne::MaterialLibrary::get_material_ptr(
    const std::string& mat_name) const {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return it->second;
  } else {
    return std::make_shared<pyne::Material>(pyne::Material());
  }
}

void pyne::MaterialLibrary::write_hdf5(const std::string& filename,
                                       const std::string& datapath,
                                       bool h5_overwrite) const {
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  // Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  hid_t matlib_grp_id;
  hid_t data_id;
  hid_t db;
  if (!pyne::file_exists(filename)) {
    // Create the file
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  } else {
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  }
  // Check if root_path exist and what type it is
  std::string root_path = "/material_library";

  H5O_info_t object_info;
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check if root_path exist and what type it is
  status =
      H5Oget_info_by_name(db, root_path.c_str(), &object_info, H5P_DEFAULT);
  if (status == 0) {
    // "/material_library" not a group: fail!
    if (object_info.type != H5O_TYPE_GROUP) {
      throw std::runtime_error(
          "Non-group/non-dataset object /material_library already exists in "
          "the file. Can't write the Material");
    } else {  // "/material_library" is a group get its hid
      matlib_grp_id = H5Gopen2(db, root_path.c_str(), H5P_DEFAULT);
    }
  } else {  // "/material_library" does not exist -> create it !
    matlib_grp_id = H5Gcreate2(db, root_path.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT);
  }

  // Group "/root_path/datapath" does exist throw an error or remove it
  std::string full_path = root_path + datapath;
  if (h5wrap::path_exists(db, full_path)) {
    if (h5_overwrite) {
      H5Ldelete(db, full_path.c_str(), H5P_DEFAULT);
    } else {
      throw std::runtime_error(
          "datapath" + full_path +
          "already exist, use \"h5_overwrite=true\" option to "
          "overwrite existing datapath");
    }
  }

  // create "/root_path/datapath" Group
  data_id =
      H5Gcreate2(db, full_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::string compath = "composition";
  std::string nucpath = "nuclidelist";
  write_hdf5_nucpath(data_id, nucpath);
  std::vector<int> nuc_list;
  nuc_list.assign(nuclist.begin(), nuclist.end());
  for (auto mat : material_library) {
    mat.second->write_hdf5_datapath(data_id, compath, -0.0,
                                    DEFAULT_MAT_CHUNKSIZE, nuc_list);
  }

  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Gclose(data_id);
  H5Gclose(matlib_grp_id);
  H5Fclose(db);
}

void pyne::MaterialLibrary::write_hdf5_nucpath(hid_t db,
                                               std::string nucpath) const {
  int nuc_size;
  nuc_size = nuclist.size();
  std::vector<int> nuc_data;
  nuc_data.assign(nuclist.begin(), nuclist.end());

  hsize_t nuc_dims[1];
  nuc_dims[0] = nuc_size;
  hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
  hid_t nuc_set = H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuc_space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           nuc_data.data());

  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Sclose(nuc_space);
  H5Dclose(nuc_set);
}

void pyne::MaterialLibrary::append_to_nuclist(const pyne::Material& mat) {
  pyne::comp_map mat_comp = mat.comp;
  for (auto nuclide : mat_comp) {
    nuclist.insert(nuclide.first);
  }
}

int pyne::MaterialLibrary::get_length_of_table(
    hid_t db, const std::string& datapath) const {
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);
  // failure to read the dapaspace return 0 for the size
  if (arr_space < 0) return 0;

  hsize_t arr_dims[1];
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  herr_t status = H5Eclear(H5E_DEFAULT);
  return arr_dims[0];
}
