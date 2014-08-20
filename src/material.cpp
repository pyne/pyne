// Material.cpp
// The very central Material class
// -- Anthony Scopatz

#include <string>
#include <vector>
#include <iomanip>  // std::setprecision
#include <math.h>   // modf

#ifndef PYNE_IS_AMALGAMATED
#include "material.h"
#include "nucname.h"
#endif

// h5wrap template
template double h5wrap::get_array_index(hid_t, int, hid_t);



/***************************/
/*** Protected Functions ***/
/***************************/

double pyne::Material::get_comp_sum() {
  // Sums the weights in the composition dictionary
  double sum = 0.0;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    sum = sum + i->second;
  }
  return sum;
};



void pyne::Material::norm_comp() {
  double sum = get_comp_sum();
  if (sum != 1.0 && sum != 0.0) {
    for (comp_iter i = comp.begin(); i != comp.end(); i++)
      i->second = i->second / sum;
  }

  if (mass < 0.0)
    mass = sum;
}






void pyne::Material::_load_comp_protocol0(hid_t db, std::string datapath, int row) {
  hid_t matgroup = H5Gopen2(db, datapath.c_str(), H5P_DEFAULT);
  hid_t nucset;
  double nucvalue;
  ssize_t nuckeylen;
  std::string nuckey;

  // get the number of members in the material group
  H5G_info_t group_info;
  H5Gget_info(matgroup, &group_info);
  hsize_t matG = group_info.nlinks;

  // Iterate over datasets in the group.
  for (int matg = 0; matg < matG; matg++) {
    nuckeylen = 1 + H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                        NULL, 0, H5P_DEFAULT);
    char * nkey = new char[nuckeylen];
    nuckeylen = H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                    nkey, nuckeylen, H5P_DEFAULT);
    nuckey = nkey;
    nucset = H5Dopen2(matgroup, nkey, H5P_DEFAULT);
    nucvalue = h5wrap::get_array_index<double>(nucset, row);

    if (nuckey == "Mass" || nuckey == "MASS" || nuckey == "mass")
      mass = nucvalue;
    else
      comp[pyne::nucname::id(nuckey)] = nucvalue;

    H5Dclose(nucset);
    delete[] nkey;
  };

  // Set meta data
  atoms_per_molecule = -1.0;
};



void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath, int row) {
  std::string nucpath;
  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  hsize_t data_offset[1] = {row};
  if (row < 0) {
    // Handle negative row indices
    hid_t data_space = H5Dget_space(data_set);
    hsize_t data_dims[1];
    H5Sget_simple_extent_dims(data_space, data_dims, NULL);
    data_offset[0] += data_dims[0];
  };

  // Grab the nucpath
  hid_t nuc_attr = H5Aopen(data_set, "nucpath", H5P_DEFAULT);
  H5A_info_t nuc_info;
  H5Aget_info(nuc_attr, &nuc_info);
  hsize_t nuc_attr_len = nuc_info.data_size;
  hid_t str_attr = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_attr, nuc_attr_len);
  char * nucpathbuf = new char [nuc_attr_len];
  H5Aread(nuc_attr, str_attr, nucpathbuf);
  nucpath = std::string(nucpathbuf, nuc_attr_len);
  delete[] nucpathbuf;

  // Grab the nuclides
  std::vector<int> nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
  int nuc_size = nuclides.size();
  hsize_t nuc_dims[1] = {nuc_size};

  // Get the data hyperslab
  hid_t data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

  // Get memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, NULL);

  // Get material type
  size_t material_data_size = sizeof(pyne::material_data) + sizeof(double)*nuc_size;
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule", HOFFSET(pyne::material_data, atoms_per_mol),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_data * mat_data = new material_data [material_data_size];

  // Finally, get data and put in on this instance
  H5Dread(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  mass = (*mat_data).mass;
  density = (*mat_data).density;
  atoms_per_molecule = (*mat_data).atoms_per_mol;
  for (int i = 0; i < nuc_size; i++)
    comp[nuclides[i]] = (double) (*mat_data).comp[i];

  delete[] mat_data;
  H5Tclose(str_attr);

  //
  // Get metadata from associated dataset, if available
  //
  std::string attrpath = datapath + "_metadata";
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (!attrpath_exists)
    return;

  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;
  hvl_t attrdata [1];

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // Get the metadata from the file
  metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, NULL);
  H5Dread(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT, attrdata);

  // convert to in-memory JSON
  Json::Reader reader;
  reader.parse((char *) attrdata[0].p, (char *) attrdata[0].p+attrdata[0].len, metadata, false);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);
};





void pyne::Material::from_hdf5(char * filename, char * datapath, int row, int protocol) {
  std::string fname (filename);
  std::string dpath (datapath);
  from_hdf5(fname, dpath, row, protocol);
};



void pyne::Material::from_hdf5(std::string filename, std::string datapath, int row, int protocol) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(filename);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (!datapath_exists)
    throw h5wrap::PathNotFound(filename, datapath);

  // Clear current content
  comp.clear();

  // Load via various protocols
  if (protocol == 0)
    _load_comp_protocol0(db, datapath, row);
  else if (protocol == 1)
    _load_comp_protocol1(db, datapath, row);
  else
    throw pyne::MaterialProtocolError();

  // Close the database
  status = H5Fclose(db);

  // Renormalize the composition, just to be safe.
  norm_comp();
};





void pyne::Material::write_hdf5(char * filename, char * datapath, char * nucpath, float row, int chunksize) {
  std::string fname (filename);
  std::string groupname (datapath);
  std::string nuclist (nucpath);
  write_hdf5(fname, groupname, nuclist, row, chunksize);
};



void pyne::Material::write_hdf5(std::string filename, std::string datapath,
                                std::string nucpath, float row, int chunksize) {
  int row_num = (int) row;

  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
  // Create new/open datafile.
  hid_t db;
  if (pyne::file_exists(filename)) {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  }
  else
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  //
  // Read in nuclist if available, write it out if not
  //
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  std::vector<int> nuclides;
  int nuc_size;
  hsize_t nuc_dims[1];

  if (nucpath_exists) {
    nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
    nuc_size = nuclides.size();
    nuc_dims[0] = nuc_size;
  } else {
    nuclides = std::vector<int>();
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
      nuclides.push_back(i->first);
    nuc_size = nuclides.size();

    // Create the data if it doesn't exist
    int nuc_data [nuc_size];
    for (int n = 0; n != nuc_size; n++)
      nuc_data[n] = nuclides[n];
    nuc_dims[0] = nuc_size;
    hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
    hid_t nuc_set = H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuc_space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nuc_data);
    H5Fflush(db, H5F_SCOPE_GLOBAL);
  };


  //
  // Write out the data itself to the file
  //
  hid_t data_set, data_space, data_hyperslab;
  int data_rank = 1;
  hsize_t data_dims[1] = {1};
  hsize_t data_max_dims[1] = {H5S_UNLIMITED};
  hsize_t data_offset[1] = {0};

  size_t material_data_size = sizeof(pyne::material_data) + sizeof(double)*nuc_size;
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule", HOFFSET(pyne::material_data, atoms_per_mol),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp),
            comp_values_array_type);

  material_data * mat_data  = new material_data[material_data_size];
  (*mat_data).mass = mass;
  (*mat_data).density = density;
  (*mat_data).atoms_per_mol = atoms_per_molecule;
  for (int n = 0; n != nuc_size; n++) {
    if (0 < comp.count(nuclides[n]))
      (*mat_data).comp[n] = comp[nuclides[n]];
    else
      (*mat_data).comp[n] = 0.0;
  };

  // get / make the data set
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (datapath_exists) {
    data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    data_space = H5Dget_space(data_set);
    data_rank = H5Sget_simple_extent_dims(data_space, data_dims, data_max_dims);

    // Determine the row size.
    if (std::signbit(row))
      row_num = data_dims[0] + row;  // careful, row is negative

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(data_set, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    // Get full space
    data_space = H5Screate_simple(1, data_dims, data_max_dims);

    // Make data set properties to enable chunking
    hid_t data_set_params = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] ={chunksize};
    H5Pset_chunk(data_set_params, 1, chunk_dims);
    H5Pset_deflate(data_set_params, 1);

    material_data * data_fill_value  = new material_data[material_data_size];
    (*data_fill_value).mass = -1.0;
    (*data_fill_value).density= -1.0;
    (*data_fill_value).atoms_per_mol = -1.0;
    for (int n = 0; n != nuc_size; n++)
      (*data_fill_value).comp[n] = 0.0;
    H5Pset_fill_value(data_set_params, desc, &data_fill_value);

    // Create the data set
    data_set = H5Dcreate2(db, datapath.c_str(), desc, data_space, H5P_DEFAULT,
                            data_set_params, H5P_DEFAULT);
    H5Dset_extent(data_set, data_dims);

    // Add attribute pointing to nuc path
    hid_t nuc_attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(nuc_attr_type, nucpath.length());
    hid_t nuc_attr_space = H5Screate(H5S_SCALAR);
    hid_t nuc_attr = H5Acreate2(data_set, "nucpath", nuc_attr_type, nuc_attr_space,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(nuc_attr, nuc_attr_type, nucpath.c_str());
    H5Fflush(db, H5F_SCOPE_GLOBAL);

    // Remember to de-allocate
    delete[] data_fill_value;
  };

  // Get the data hyperslab
  data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

  // Get a memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, data_max_dims);

  // Write the row...
  H5Dwrite(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  // Close out the Dataset
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(data_set);
  H5Sclose(data_space);
  H5Tclose(desc);

  //
  // Write out the metadata to the file
  //
  std::string attrpath = datapath + "_metadata";
  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // get / make the data set
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (attrpath_exists) {
    metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
    metadatapace = H5Dget_space(metadataet);
    attrrank = H5Sget_simple_extent_dims(metadatapace, data_dims, data_max_dims);

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(metadataet, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    hid_t metadataetparams;
    hsize_t attrchunkdims [1];

    // Make data set properties to enable chunking
    metadataetparams = H5Pcreate(H5P_DATASET_CREATE);
    attrchunkdims[0] = chunksize;
    H5Pset_chunk(metadataetparams, 1, attrchunkdims);
    H5Pset_deflate(metadataetparams, 1);

    hvl_t attrfillvalue [1];
    attrfillvalue[0].len = 3;
    attrfillvalue[0].p = (char *) "{}\n";
    H5Pset_fill_value(metadataetparams, attrtype, &attrfillvalue);

    // make dataset
    metadatapace = H5Screate_simple(1, data_dims, data_max_dims);
    metadataet = H5Dcreate2(db, attrpath.c_str(), attrtype, metadatapace,
                         H5P_DEFAULT, metadataetparams, H5P_DEFAULT);
    H5Dset_extent(metadataet, data_dims);
    H5Pclose(metadataetparams);
  };

  // set the attr string
  hvl_t attrdata [1];
  Json::FastWriter writer;
  std::string metadatatr = writer.write(metadata);
  attrdata[0].p = (char *) metadatatr.c_str();
  attrdata[0].len = metadatatr.length();

  // write the attr
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, data_max_dims);
  H5Dwrite(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT, attrdata);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);
  // Remember the milk!
  // ...by which I mean to deallocate
  delete[] mat_data;
};

std::string pyne::Material::mcnp(std::string frac_type) {
  //////////////////// Begin card creation ///////////////////////
  std::ostringstream oss;
  // 'name'
  if (metadata.isMember("name")) {
    oss << "C name: " << metadata["name"].asString() << std::endl;
  }
  // 'density'
  if (density != -1.0) {
     std::stringstream ds;
     ds << std::setprecision(1) << std::fixed << "C density = " << density << std::endl;
     oss << ds.str();
  }
  // 'source'
  if (metadata.isMember("source")) {
     oss << "C source: " << metadata["source"].asString() << std::endl;
  }
  // Metadata comments
  if (metadata.isMember("comments")) {
    std::string comment_string = "comments: " + metadata["comments"].asString();
    // Include as is if short enough
    if (comment_string.length() <= 77) {
      oss << "C " << comment_string << std::endl;
    }
    else { // otherwise create a remainder string and iterate/update it
      oss << "C " << comment_string.substr(0,77) << std::endl;
      std::string remainder_string = comment_string.substr(77);
      while (remainder_string.length() > 77) {
        oss << "C " << remainder_string.substr(0,77) << std::endl;
        remainder_string.erase(0,77);
      }
      if (remainder_string.length() > 0) {
        oss << "C " << remainder_string << std::endl;
      }
    }
  }

  // Metadata mat_num
  oss << "m";
  if (metadata.isMember("mat_number")) {
    int mat_num = metadata["mat_number"].asInt();
    oss << mat_num << std::endl;
  } else {
    oss << "?" << std::endl;
  }

  // Set up atom or mass frac map
  std::map<int, double> fracs;
  std::string frac_sign;

  if ("atom" == frac_type) {
    fracs = to_atom_frac();
    frac_sign = "";
  } else {
    fracs = comp;
    frac_sign = "-";
  }

  // iterate through frac map
  // This is an awkward pre-C++11 way to put an int to a string
  std::stringstream ss;
  std::string nucmcnp;
  std::string table_item;
  for(pyne::comp_iter i = fracs.begin(); i != fracs.end(); ++i) {
    // Clear first
    ss.str(std::string() );
    ss.clear();
    ss << pyne::nucname::mcnp(i->first );
    nucmcnp = ss.str();

    int mcnp_id;
    mcnp_id = pyne::nucname::mcnp(i->first );
    // Spaces are important for tests
    table_item = metadata["table_ids"][nucmcnp].asString();
    if ( !table_item.empty() ) {
      oss << "     " << mcnp_id << "." << table_item << " ";
    } else {
      oss << "     " << mcnp_id << " ";
    }
    // The int needs a little formatting
    std::stringstream fs;
    fs << std::setprecision(4) << std::scientific << frac_sign << i->second \
       << std::endl;
    oss << fs.str();
  }

  return oss.str();
}

//---------------------------------------------------------------------------//
// Create a set out of the static string array. 
std::set<std::string> FLUKA_builtin(pyne::flukaMatStrings, 
                                    pyne::flukaMatStrings+pyne::FLUKA_MAT_NUM);

//---------------------------------------------------------------------------//
// notBuiltin
//---------------------------------------------------------------------------//
// Convenience function
// This is written as a negative because that is what we care about
bool pyne::Material::notBuiltin(std::string fluka_name)
{
  return (FLUKA_builtin.find(fluka_name) == FLUKA_builtin.end());
}


//---------------------------------------------------------------------------//
// fluka
//---------------------------------------------------------------------------//
// Main external call
std::string pyne::Material::fluka(int id )
{
  std::stringstream rs;

  // Element, one nucid
  if (comp.size() == 1) {
    rs << write_material(id);
  } else if (comp.size() > 1 ) {
  // Compound
    rs << write_compound(id);
  } else {
    rs << "There is no nuclide information in the Material Object" << std::endl;
  }
  return rs.str();
}

//---------------------------------------------------------------------------//
// write_material
//---------------------------------------------------------------------------//
//
// Requirement:  the material upon which this function is called has
//               exactly one nucid component, i.e. it is elemental
// Do not assume fluka_name is defined in the metadata.  This function
// may be called from a user-defined material, i.e. on that is not 
// read out of a UW^2-tagged geometry file, and thus does not have
// certain metadata.
// Result:        The MATERIAL line is returned.  
//                The fluka_name metadata is set
//                The true, final material index metadata is set
#define FAKE_ZNUM 999
std::string pyne::Material::write_material(int id)
{
  std::stringstream ms;
  std::string fluka_name; 
  int za_id = FAKE_ZNUM;

  pyne::nucname::zzname_t zfd = pyne::nucname::get_zz_fluka();
   
  if (metadata.isMember("fluka_name")) {
    fluka_name = metadata["fluka_name"].asString();
  } else {  // Should be elemental
    if (comp.size() > 1 ) {
      std::cerr << "Error: this mix is a compound, there should be a fluka_name defined."
                << std::endl;
      // ToDo:  Perhaps should throw an exception here
      return ms.str();
    }

    // The nucid of the first component is the only nucid
    int nucid = comp.begin()->first;
    fluka_name = zfd[nucid];
    metadata["fluka_name"] = fluka_name;
  }

  if (notBuiltin(fluka_name)) {  
    int nucid = comp.begin()->first;
    std::stringstream ss;
    ss << id;
    metadata["fluka_material_index"] = ss.str();
    ms << material_component(id, nucid, fluka_name);
  }
  // could be empty
  return ms.str();
}
    
//---------------------------------------------------------------------------//
// material_component
//---------------------------------------------------------------------------//
// Material has only one component, 
// Density is either object density or it is ignored ==> use object density
std::string pyne::Material::material_component(int fid, int nucid, std::string fluka_name)
{
  int znum;
  if (FAKE_ZNUM != nucid) {
    znum = pyne::nucname::znum(nucid); 
  } else {
    znum = nucid;
  }

  double atomic_mass;
  if (0 != pyne::NUC_DATA_PATH.length() ) { 
    // for compounds (i.e., unrecognized nucids), this will be 0
    atomic_mass = pyne::atomic_mass(nucid);
  } else {
    atomic_mass = -1; 
  }  

  return material_line(znum, atomic_mass, fid, fluka_name);
}

//---------------------------------------------------------------------------//
// material_line
//---------------------------------------------------------------------------//
// Given all the info, return the Material string
std::string pyne::Material::material_line (int znum, double atomic_mass, 
                                           int fid, std::string fluka_name)
{

  std::stringstream ls;

  if (metadata.isMember("comments") ) {
     std::string comment = metadata["comments"].asString();
     ls << "* " << comment;
     ls << std::endl;
  }
  ls << std::setw(10) << std::left << "MATERIAL";
  ls << std::setprecision(0) << std::fixed << std::showpoint <<
        std::setw(10) << std::right << (float)znum; 

  // Control the maximum precision:  this will print however many digits 
  // there are after the decimal, up to a maximum of six
  ls.unsetf(std::ios::showpoint);
  ls.unsetf(std::ios::floatfield);
  ls.precision(6);
  ls << std::setw(10) << std::right << atomic_mass;

  // Density formatting depends on whether it is a true float
  // Note this is the current object density, and may or may not be meaningful
  double intpart;
  modf (density, &intpart);
  if (density == intpart) {
    // There is no part after the decimal
    ls << std::setprecision(0) << std::fixed << std::showpoint 
       << std::setw(10) << std::right << density;
  } else {
    // Use the previous field's setting (atomic mass: max precision = 6)
    ls << std::setw(10) << std::right << density;
  }

  ls << std::setprecision(0) << std::fixed << std::showpoint <<
        std::setw(10) << std::right << (float)fid;
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::left << fluka_name << std::endl;

  return ls.str();
}
//---------------------------------------------------------------------------//
// write_compound
//---------------------------------------------------------------------------//
// Returns
// a) MATERIAL lines for those components that need it
// b) MATERIAL line for compound
// c) COMPOUND lines
std::string pyne::Material::write_compound(int id)
{
  std::stringstream ss;
  std::map<double, std::string> frac_name_map;
  std::string compound_string = "";
  std::vector<std::string> material_names;

  pyne::nucname::zzname_t zfd = pyne::nucname::get_zz_fluka();

  // The nucid doesn't make sense for a compound
  int znum = 999;
  double atomic_mass = 999.;
  // This better be true
  std::string compound_name;
  if (metadata.isMember("fluka_name")) {
    compound_name = metadata["fluka_name"].asString();
  } else {
    // ToDo:  what to do?  
    std::cerr << "Error:  metadata \"fluka_name\" expected." << std::endl;
    compound_name = "NotFound";
  }  
  ss << material_line(znum, atomic_mass, id, compound_name);
  
  int count = 1;
  for (comp_iter nuc = comp.begin(); nuc != comp.end(); ++nuc) {
    if ( count%2 == 0 ) {
	ss << std::setw(10) << std::right << zfd[nuc->first];
	ss << std::setw(10) << std::right << nuc->second;
    } else if ( count%3 == 0 ) {
	ss << std::setw(10) << std::right << zfd[nuc->first];
	ss << std::setw(10) << std::right << nuc->second;
	ss << std::setw(10) << std::left << compound_name << std::endl;
	count = 0;
    } else  {
	ss << std::setw(10) << std::left << "COMPOUND";
	ss << std::setw(10) << std::right << zfd[nuc->first];
	ss << std::setw(10) << std::right << nuc->second;
    }
    count++;
  }
  ss<< std::endl;
  return ss.str();
}

void pyne::Material::from_text(char * filename) {
  std::string fname (filename);
  from_text(fname);
};


void pyne::Material::from_text(std::string filename) {
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // New filestream
  std::ifstream f;
  f.open(filename.c_str());

  // Read in
  comp.clear();
  std::string keystr, valstr;

  while ( !f.eof() ) {
    f >> keystr;

    if (0 == keystr.length())
      continue;

    if (keystr == "Mass"){
      f >> valstr;
      mass = pyne::to_dbl(valstr);
    } else if (keystr == "Density") {
      f >> valstr;
      density = pyne::to_dbl(valstr);
    } else if (keystr == "APerM") {
      f >> valstr;
      atoms_per_molecule = pyne::to_dbl(valstr);
    } else if (pyne::nucname::isnuclide(keystr)) {
      f >> valstr;
       comp[pyne::nucname::id(keystr)] = pyne::to_dbl(valstr);
    } else {
      getline(f, valstr);
      valstr= valstr.substr(0, valstr.length()-1);
      metadata[keystr]= valstr;
      continue;
   }
   };

   f.close();
   norm_comp();
};



void pyne::Material::write_text(char * filename) {
  std::string fname (filename);
  write_text(fname);
};


void pyne::Material::write_text(std::string filename) {
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);

  Json::Reader reader;
  std::vector<std::string> obj = metadata.getMemberNames();

  if (0 <= mass)
    f << "Mass    " << mass << "\n";

  if (0 <= density)
    f << "Density "  << density << "\n";

  if (0 <= atoms_per_molecule)
    f << "APerM   " << atoms_per_molecule << "\n";

  for (int i=0; i < metadata.size(); i=i+2){
    f <<metadata.get(obj.at(i), "") << metadata.get(obj.at(i+1), "");
  }

  std::string nuc_name;
  for(pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    nuc_name = pyne::nucname::name( i->first ) + "  ";
    while (nuc_name.length() < 8)
      nuc_name += " ";
    f << nuc_name << i->second << "\n";
  };

  f.close();
};


void pyne::Material::load_json(Json::Value json) {
  Json::Value::Members keys = json["comp"].getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  comp.clear();
  for (; ikey != ikey_end; ++ikey)
    comp[nucname::id(*ikey)] = json["comp"][*ikey].asDouble();
  norm_comp();
  mass = json["mass"].asDouble();
  density = json["density"].asDouble();
  atoms_per_molecule = json["atoms_per_molecule"].asDouble();
  metadata = json["metadata"];
};


Json::Value pyne::Material::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);
  Json::Value jcomp = Json::Value(Json::objectValue);
  json["mass"] = mass;
  json["density"] = density;
  json["atoms_per_molecule"] = atoms_per_molecule;
  json["metadata"] = metadata;
  for(comp_iter i = comp.begin(); i != comp.end(); i++)
    jcomp[nucname::name(i->first)] = (i->second);
  json["comp"] = jcomp;
  return json;
};


void pyne::Material::from_json(char * filename) {
  std::string fname (filename);
  from_json(fname);
};

void pyne::Material::from_json(std::string filename) {
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
};


void pyne::Material::write_json(char * filename) {
  std::string fname (filename);
  write_json(fname);
};

void pyne::Material::write_json(std::string filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
};


/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/

pyne::Material::Material() {
  // Empty Material constructor
  mass = -1.0;
  density = -1.0;
  atoms_per_molecule = -1.0;
  metadata = Json::Value(Json::objectValue);
}


pyne::Material::Material(pyne::comp_map cm, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;
  if (!comp.empty())
    norm_comp();
};



pyne::Material::Material(char * filename, double m, double d, double apm,
                         Json::Value attributes) {
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  std::string fname (filename);
  if (!pyne::file_exists(fname))
    throw pyne::FileNotFound(fname);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(fname.c_str());
  if (ish5)
    from_hdf5(fname);
  else
    from_text(fname);
};


pyne::Material::Material(std::string filename, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename);
  else
    from_text(filename);
};


pyne::Material::~Material() {
};



/*--- Method definitions ---*/


std::ostream& operator<<(std::ostream& os, pyne::Material mat) {
  //print the Mass Stream to stdout
  os << "\tMass: " << mat.mass << "\n";
  os << "\t---------\n";
  for(pyne::comp_iter i = mat.comp.begin(); i != mat.comp.end(); i++)
  {
    os << "\t" << pyne::nucname::name( i->first ) << "\t" << i->second << "\n";
  };
  return os;
};

// Note this refines << for an inheritor of std::ostream.
std::ostringstream& operator<<(std::ostringstream& os, pyne::Material mat) {
  return os;
}

void pyne::Material::normalize () {
  // normalizes the mass
  mass = 1.0;
};


pyne::comp_map pyne::Material::mult_by_mass() {
  // bypass calculation if already normalized.
  if (mass == 1.0)
    return comp;

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    cm[i->first] = (i->second) * mass;
  };
  return cm;
};



double pyne::Material::molecular_mass(double apm) {
  // Calculate the atomic weight of the Material
  double inverseA = 0.0;

  for (pyne::comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++)
    inverseA += (nuc->second) / pyne::atomic_mass(nuc->first);

  if (inverseA == 0.0)
    return inverseA;

  // select the atoms per mol
  double atsperm = 1.0; // default to 1.0
  if (0.0 <= apm) {
    atsperm = apm;            // take the function argument, if valid
    if (atoms_per_molecule < 0.0)
      atoms_per_molecule = apm;     // Store the function argument on class, if class has no value
  } else if (0.0 <= atoms_per_molecule)
    atsperm = atoms_per_molecule;  // select the class's value

  return atsperm / inverseA;
};


pyne::Material pyne::Material::expand_elements() {
  // Expands the natural elements of a material and returns a new material note
  // that this implementation relies on the fact that maps of ints are stored in
  // a sorted manner in C++.
  int n, nabund, znuc, zabund;
  comp_map newcomp;
  std::map<int, double>::iterator abund_itr, abund_end;
  if (pyne::natural_abund_map.empty())
    pyne::_load_atomic_mass_map();
  abund_itr = pyne::natural_abund_map.begin();
  abund_end = pyne::natural_abund_map.end();
  zabund = nucname::znum((*abund_itr).first);
  for (comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++) {
    if(abund_itr == abund_end)
      newcomp.insert(*nuc);
    else if(0 == nucname::anum((*nuc).first)) {
      n = (*nuc).first;
      znuc = nucname::znum(n);
      if (znuc < zabund) {
        newcomp.insert(*nuc);
        continue;
      }
      while(zabund <= znuc) {
        nabund = (*abund_itr).first;
        if (zabund == znuc && 0 != nucname::anum(nabund) && 0.0 != (*abund_itr).second)
          newcomp[nabund] = (*abund_itr).second * (*nuc).second * \
                            atomic_mass_map[nabund] / atomic_mass_map[n];
        else if (n == nabund && 0.0 == (*abund_itr).second)
          newcomp.insert(*nuc);
        abund_itr++;
        if (abund_itr == abund_end) {
          zabund = INT_MAX;
          break;
        }
        zabund = nucname::znum(nabund);
      };
    } else
      newcomp.insert(*nuc);
  };
  return Material(newcomp, mass, density, atoms_per_molecule, metadata);
};

pyne::Material pyne::Material::collapse_elements(std::set<int> exception_ids) {
  ////////////////////////////////////////////////////////////////////////
  // Assumptions
  //    - list passed in is formed of znum-anum of Fluka-named isotopes,
  //      since we want to preserve the full nucid of any such material in
  //      the problem
  // Algorithm
  // for each component listed in this material that has a nonzero frac or 
  //    weight amount, look at its 'stripped' nucid, that is, the last four 
  //    places replaced by zeros.
  //    if it's on the exception list, copy the component
  //    else it is to be collapsed 
  //       => add its frac to the component of the znum
  //  
  // * When from_hdf5 reads from a file the comp iterator will produce a 
  //   hit for EVERY nucid in EVERY material in the file.  Only the nucids
  //   belonging to the CURRENT material have a nonzero fraction/mass amount
  ///////////////////////////////////////////////////////////////////////// 
  pyne::comp_map cm;
  
  for (pyne::comp_iter ptr = comp.begin(); ptr != comp.end(); ptr++) {
      if (0 < ptr->second) {
        // There is a nonzero amount of this nucid in the current material, 
    // check if znum and anum are in the exception list, 
        int cur_stripped_id = nucname::znum(ptr->first)*10000000 
                        + nucname::anum(ptr->first)*10000;
        if (0 < exception_ids.count(cur_stripped_id)) {
          // The znum/anum combination identify the current material as a 
      // fluka-named exception list => copy, don't collapse
          cm[ptr->first] = (ptr->second) * mass;
        } else {
          // Not on exception list => add frac to id-component
          int znum_id = nucname::id(nucname::znum(ptr->first));
          cm[znum_id] += (ptr->second) * mass;
        }
      }
  }
  // Copy 
  pyne::Material collapsed = pyne::Material(cm, mass, density, 
                                            atoms_per_molecule, metadata);
  return collapsed;
}

// Wrapped version for calling from python
pyne::Material pyne::Material::collapse_elements(int** int_ptr_arry ) {
    std::set<int> nucvec;
    // Set first pointer to first int pointed to by arg
    int *int_ptr = *int_ptr_arry;
    while (int_ptr != NULL)
    {
      nucvec.insert(*int_ptr);
      int_ptr++;
    }
    return collapse_elements(nucvec);
}

double pyne::Material::mass_density(double num_dens, double apm) {
  if (0.0 <= num_dens) {
    double mw = molecular_mass(apm);
    density = num_dens * mw / pyne::N_A / atoms_per_molecule;
  };
  return density;
};


double pyne::Material::number_density(double mass_dens, double apm) {
  if (0 <= mass_dens)
    density = mass_dens;
  double mw = molecular_mass(apm);
  double num_dens = density * pyne::N_A * atoms_per_molecule / mw;
  return num_dens;
};


/*--- Stub-Stream Computation ---*/

pyne::Material pyne::Material::sub_mat(std::set<int> nucset) {
  // Grabs a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ( 0 < nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, -1);
};



pyne::Material pyne::Material::sub_mat(std::set<std::string> nucset) {
  // Grabs a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  };

  return sub_mat(iset);
};



pyne::Material pyne::Material::set_mat (std::set<int> nucset, double value) {
  // Sets a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;

  // Add non-set components
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  // Add set component
  for (std::set<int>::iterator nuc = nucset.begin(); nuc != nucset.end(); nuc++)
    cm[*nuc] = value;

  return pyne::Material(cm, -1, -1);
};



pyne::Material pyne::Material::set_mat(std::set<std::string> nucset, double value) {
  // Sets a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  };

  return set_mat(iset, value);
};




pyne::Material pyne::Material::del_mat(std::set<int> nucset) {
  // Removes a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    // Only add to new comp if not in nucset
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, -1);
};



pyne::Material pyne::Material::del_mat (std::set<std::string> nucset) {
  // Removes a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  };

  return del_mat(iset);
};






pyne::Material pyne::Material::sub_range(int lower, int upper) {
  // Grabs a sub-material from this mat based on a range of integers.
  if (upper < lower)
  {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  };

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1,-1);
};



pyne::Material pyne::Material::set_range(int lower, int upper, double value) {
// Sets a sub-material from this mat based on a range of integers.
if (upper < lower) {
int temp_upper = upper;
upper = lower;
lower = temp_upper;
};

pyne::comp_map cm;
for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
  if ((lower <= (i->first)) && ((i->first) < upper))
    cm[i->first] = value;
  else
    cm[i->first] = (i->second) * mass;
};

  return pyne::Material(cm, -1,-1);
};



pyne::Material pyne::Material::del_range(int lower, int upper) {
  // Removes a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  };

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((upper <= (i->first)) || ((i->first) < lower))
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, -1);
};










pyne::Material pyne::Material::sub_elem(int elem) {
  // Returns a material of the element that is a submaterial of this one.
  return sub_range(elem, elem + 10000000);
};



pyne::Material pyne::Material::sub_lan() {
  // Returns a material of Lanthanides that is a sub-material of this one.
  return sub_range(570000000, 720000000);
};



pyne::Material pyne::Material::sub_act() {
  //Returns a material of Actindes that is a sub-material of this one.
  return sub_range(890000000, 1040000000);
};


pyne::Material pyne::Material::sub_tru() {
  // Returns a material of Transuranics that is a sub-material of this one.
  return sub_range(930000000, INT_MAX);
};



pyne::Material pyne::Material::sub_ma() {
  // Returns a material of Minor Actinides that is a sub-material of this one.
  return sub_range(930000000, 1040000000).del_range(940000000, 950000000);
};



pyne::Material pyne::Material::sub_fp() {
  // Returns a material of Fission Products that is a sub-material of this one.
  return sub_range(0, 890000000);
};




/*--- Atom Frac Functions ---*/

std::map<int, double> pyne::Material::to_atom_frac() {
  // Returns an atom fraction map from this material's composition
  // the material's molecular mass
  double mat_mw = molecular_mass();

  std::map<int, double> atom_fracs = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_fracs[ci->first] = (ci->second) * mat_mw / pyne::atomic_mass(ci->first);

  return atom_fracs;
};


void pyne::Material::from_atom_frac(std::map<int, double> atom_fracs) {
  // atom frac must be of the form {nuc: af}, eg, water
  //  80160: 1.0
  //  10010: 2.0

  // clear existing components
  comp.clear();
  atoms_per_molecule = 0.0;

  for (std::map<int, double>::iterator afi = atom_fracs.begin(); afi != atom_fracs.end(); afi++) {
    comp[afi->first] = (afi->second) * pyne::atomic_mass(afi->first);
    atoms_per_molecule += (afi->second);
  };

  norm_comp();
};



std::vector<std::pair<double, double> > pyne::Material::gammas() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
        state_id = nucname::id_to_state_id(ci->first);
    else
        state_id = ci->first;

    std::vector<std::pair<double, double> > raw_gammas = pyne::gammas(state_id);
    for (int i = 0; i < raw_gammas.size(); ++i) {
      result.push_back(std::make_pair(raw_gammas[i].first,
        atom_fracs[ci->first]*raw_gammas[i].second));
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::Material::xrays() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
        state_id = nucname::id_to_state_id(ci->first);
    else
        state_id = ci->first;

    std::vector<std::pair<double, double> > raw_xrays = pyne::xrays(state_id);
    for (int i = 0; i < raw_xrays.size(); ++i) {
      result.push_back(std::make_pair(raw_xrays[i].first,
        atom_fracs[ci->first]*raw_xrays[i].second));
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::Material::photons(bool norm) {
  std::vector<std::pair<double, double> >  txray = this->xrays();
  std::vector<std::pair<double, double> >  tgammas = this->gammas();
  for (int i = 0; i < txray.size(); ++i)
    tgammas.push_back(txray[i]);
  if (norm)
    tgammas = normalize_radioactivity(tgammas);
  return tgammas;
}

std::vector<std::pair<double, double> > pyne::Material::normalize_radioactivity(
std::vector<std::pair<double, double> > unnormed) {
  std::vector<std::pair<double, double> > normed;
  double sum = 0.0;
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second))
      sum = sum + unnormed[i].second;
  }
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second)) {
      normed.push_back(std::make_pair(unnormed[i].first,
        (unnormed[i].second)/sum));
    }
  }
  return normed;
}


pyne::Material pyne::Material::operator+ (double y) {
  // Overloads x + y
  return pyne::Material(comp, mass + y, density);
};



pyne::Material pyne::Material::operator+ (Material y) {
  // Overloads x + y
  pyne::comp_map cm;
  pyne::comp_map xwgt = mult_by_mass();
  pyne::comp_map ywgt = y.mult_by_mass();

  for (pyne::comp_iter i = xwgt.begin(); i != xwgt.end(); i++) {
    if ( 0 < ywgt.count(i->first) )
      cm[i->first] = xwgt[i->first] + ywgt[i->first];
    else
      cm[i->first] = xwgt[i->first];
  };

  for (pyne::comp_iter i = ywgt.begin(); i != ywgt.end(); i++) {
    if ( 0 == cm.count(i->first) )
      cm[i->first] = ywgt[i->first];
  };

  return pyne::Material(cm, -1, -1);
};



pyne::Material pyne::Material::operator* (double y) {
  // Overloads x * y
  return pyne::Material(comp, mass * y, density);
};



pyne::Material pyne::Material::operator/ (double y) {
  // Overloads x / y
  return pyne::Material(comp, mass / y, density );
}
