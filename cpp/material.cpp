// Material.cpp
// The very central Material class
// -- Anthony Scopatz

#include "material.h"


// h5wrap template
template double h5wrap::get_array_index(hid_t, int, hid_t);



/***************************/
/*** Protected Functions ***/
/***************************/

double pyne::Material::get_comp_sum()
{
  // Sums the weights in the composition dictionary
  double sum = 0.0;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    sum = sum + i->second;
  }
  return sum;
};



void pyne::Material::norm_comp()
{
  double sum = get_comp_sum();
  if (sum != 1.0 && sum != 0.0)
  {
    for (comp_iter i = comp.begin(); i != comp.end(); i++)
      i->second = i->second / sum;
  }

  if (mass < 0.0)
    mass = sum;
}






void pyne::Material::_load_comp_protocol0(hid_t db, std::string datapath, int row)
{
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
  for (int matg = 0; matg < matG; matg++)
  {
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
      comp[pyne::nucname::zzaaam(nuckey)] = nucvalue;

    H5Dclose(nucset);
    delete[] nkey;
  };

  // Set meta data
  name = datapath.substr(datapath.rfind("/")+1, datapath.length());
  atoms_per_mol = -1.0;
};



void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath, int row)
{
  std::string nucpath;
  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  hsize_t data_offset[1] = {row};
  if (row < 0)
  {
    // Handle negative row indices
    hid_t data_space = H5Dget_space(data_set);
    hsize_t data_dims[1];
    H5Sget_simple_extent_dims(data_space, data_dims, NULL);
    data_offset[0] += data_dims[0];
  };

  // Grab the nucpath
  hid_t nuc_attr = H5Aopen(data_set, "nucpath", H5P_DEFAULT);
  hsize_t nuc_attr_len = H5Aget_storage_size(nuc_attr) / sizeof(char);
  hid_t str_attr = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_attr, nuc_attr_len);
  char * nucpathbuf = new char [nuc_attr_len];
  H5Aread(nuc_attr, str_attr, nucpathbuf);
  nucpath = nucpathbuf;
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
  size_t material_struct_size = sizeof(pyne::material_struct) + sizeof(double)*nuc_size;
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_struct_size);
  hid_t str20 = H5Tcopy(H5T_C_S1);
  H5Tset_size(str20, 20);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "name", HOFFSET(pyne::material_struct, name), str20);
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_struct, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_mol", HOFFSET(pyne::material_struct, atoms_per_mol), 
              H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_struct, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_struct * mat_data = new material_struct [material_struct_size];

  // Finally, get data and put in on this instance
  H5Dread(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  name = std::string((*mat_data).name);
  mass = (*mat_data).mass;
  atoms_per_mol = (*mat_data).atoms_per_mol;
  for (int i = 0; i < nuc_size; i++)
    comp[nuclides[i]] = (double) (*mat_data).comp[i];

  delete[] mat_data;
  H5Tclose(str_attr);
  H5Tclose(str20);

  //
  // Get attrs from associated dataset, if available
  //
  std::string attrpath = datapath + "_attrs";
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (!attrpath_exists)
    return;

  hid_t attrspace, attrtype, attrset, attrslab, attrmemspace;
  int attrrank; 
  hvl_t attrdata [1];

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // Get the attrs from the file
  attrset = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
  attrslab = H5Dget_space(attrset);
  H5Sselect_hyperslab(attrslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, NULL);
  H5Dread(attrset, attrtype, attrmemspace, attrslab, H5P_DEFAULT, attrdata);

  // convert to in-memory JSON
  Json::Reader reader;
  reader.parse((char *) attrdata[0].p, (char *) attrdata[0].p+attrdata[0].len, attrs, false);
  
  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(attrset);
  H5Sclose(attrspace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);
};





void pyne::Material::from_hdf5(char * fchar, char * dchar, int row, int protocol)
{
  std::string filename (fchar);
  std::string datapath (dchar);
  from_hdf5(filename, datapath, row, protocol);  
};



void pyne::Material::from_hdf5(std::string filename, std::string datapath, int row, int protocol)
{
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

  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

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

  // Renomalize the composition, just to be safe.
  norm_comp();
};





void pyne::Material::write_hdf5(char * fchar, char * gchar, char * nchar, float row, int chunksize)
{
  std::string filename (fchar);
  std::string groupname (gchar);
  std::string nuclist (nchar);
  write_hdf5(filename, groupname, nuclist, row, chunksize);  
};



void pyne::Material::write_hdf5(std::string filename, std::string datapath, std::string nucpath, 
                                float row, int chunksize)
{
  int row_num = (int) row;

  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Create new/open datafile.
  hid_t db;
  if (pyne::file_exists(filename))
  {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  }
  else
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //
  // Read in nuclist if available, write it out if not
  //
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  std::vector<int> nuclides;
  int nuc_size;
  hsize_t nuc_dims[1];

  if (nucpath_exists)
  {
    nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
    nuc_size = nuclides.size();
    nuc_dims[0] = nuc_size;
  }
  else
  {
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

  size_t material_struct_size = sizeof(pyne::material_struct) + sizeof(double)*nuc_size;
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_struct_size);
  hid_t str20 = H5Tcopy(H5T_C_S1);
  H5Tset_size(str20, 20);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "name", HOFFSET(pyne::material_struct, name), str20);
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_struct, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_mol", HOFFSET(pyne::material_struct, atoms_per_mol), 
                      H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_struct, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_struct * mat_data  = new material_struct[material_struct_size];
  int name_len = name.length();
  for (int i=0; i < 20; i++)
  {
    if (i < name_len)
      (*mat_data).name[i] = name[i];
    else
      (*mat_data).name[i] = NULL;
  };
  (*mat_data).mass = mass;
  (*mat_data).atoms_per_mol = atoms_per_mol;
  for (int n = 0; n != nuc_size; n++)
  {
    if (0 < comp.count(nuclides[n]))
      (*mat_data).comp[n] = comp[nuclides[n]];
    else
      (*mat_data).comp[n] = 0.0;
  };

  // get / make the data set
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (datapath_exists)
  {
    data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    data_space = H5Dget_space(data_set);
    data_rank = H5Sget_simple_extent_dims(data_space, data_dims, data_max_dims);

    // Determine the row size.
    if (std::signbit(row))
      row_num = data_dims[0] + row;  // careful, row is negative

    if (data_dims[0] <= row_num)
    {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(data_set, data_dims);
    }
    else if (data_dims[0] < 0)
      throw h5wrap::HDF5BoundsError();

    data_offset[0] = row_num;
  }
  else
  {
    // Get full space
    data_space = H5Screate_simple(1, data_dims, data_max_dims);

    // Make data set properties to enable chunking
    hid_t data_set_params = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] ={chunksize}; 
    H5Pset_chunk(data_set_params, 1, chunk_dims);

    material_struct * data_fill_value  = new material_struct[material_struct_size];
    for (int i=0; i < 20; i++)
      (*data_fill_value).name[i] = NULL;
    (*data_fill_value).mass = -1.0;
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
  H5Tclose(str20);
  H5Tclose(desc);

  //
  // Write out the attrs to the file
  //
  std::string attrpath = datapath + "_attrs";
  hid_t attrspace, attrtype, attrset, attrslab, attrmemspace;
  int attrrank; 

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // get / make the data set
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (attrpath_exists)
  {
    attrset = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
    attrspace = H5Dget_space(attrset);
    attrrank = H5Sget_simple_extent_dims(attrspace, data_dims, data_max_dims);

    if (data_dims[0] <= row_num)
    {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(attrset, data_dims);
    }
    else if (data_dims[0] < 0)
      throw h5wrap::HDF5BoundsError();

    data_offset[0] = row_num;
  }
  else
  {
    hid_t attrsetparams;
    hsize_t attrchunkdims [1];

    // Make data set properties to enable chunking
    attrsetparams = H5Pcreate(H5P_DATASET_CREATE);
    attrchunkdims[0] = chunksize; 
    H5Pset_chunk(attrsetparams, 1, attrchunkdims);

    hvl_t attrfillvalue [1];
    attrfillvalue[0].len = 3;
    attrfillvalue[0].p = (char *) "{}\n";
    H5Pset_fill_value(attrsetparams, attrtype, &attrfillvalue);

    // make dataset
    attrspace = H5Screate_simple(1, data_dims, data_max_dims);
    attrset = H5Dcreate2(db, attrpath.c_str(), attrtype, attrspace, 
                         H5P_DEFAULT, attrsetparams, H5P_DEFAULT);
    H5Dset_extent(attrset, data_dims);
  };

  // set the attr string
  hvl_t attrdata [1];
  Json::FastWriter writer;
  std::string attrstr = writer.write(attrs);
  attrdata[0].p = (char *) attrstr.c_str();
  attrdata[0].len = attrstr.length();

  // write the attr
  attrslab = H5Dget_space(attrset);
  H5Sselect_hyperslab(attrslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, data_max_dims);
  H5Dwrite(attrset, attrtype, attrmemspace, attrslab, H5P_DEFAULT, attrdata);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(attrset);
  H5Sclose(attrspace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);

  // Remember the milk!  
  // ...by which I mean to deallocate
  delete[] mat_data;
};











void pyne::Material::from_text(char * fchar)
{
  std::string filename (fchar);
  from_text(filename);
};


void pyne::Material::from_text(std::string filename)
{
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // New filestream
  std::ifstream f;
  f.open(filename.c_str());

  // Read in
  comp.clear();
  std::string keystr, valstr;

  while ( !f.eof() )
  {
    f >> keystr;

    if (0 == keystr.length())
      continue;
    else 
      f >> valstr;

    if (keystr == "Name")
      name = valstr;
    else if (keystr == "Mass")
      mass = pyne::to_dbl(valstr);
    else if (keystr == "APerM")
      atoms_per_mol = pyne::to_dbl(valstr);
    else
      comp[pyne::nucname::zzaaam(keystr)] = pyne::to_dbl(valstr);
  };

  f.close();
  norm_comp();
};



void pyne::Material::write_text(char * fchar)
{
  std::string filename (fchar);
  write_text(filename);
};


void pyne::Material::write_text (std::string filename)
{
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);

  if (0 < name.length())
    f << "Name    " << name << "\n";

  if (0 <= mass)
    f << "Mass    " << mass << "\n";

  if (0 <= atoms_per_mol)
    f << "APerM   " << atoms_per_mol << "\n";

  std::string nuc_name;
  for(pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    nuc_name = pyne::nucname::name( i->first ) + "  ";
    while (nuc_name.length() < 8)
      nuc_name += " ";
    f << nuc_name << i->second << "\n";
  };

  f.close();
};



/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/

pyne::Material::Material()
{
  // Empty Material constructor
  mass = -1.0;
  name = std::string("");
  atoms_per_mol = -1.0;
  attrs = Json::Value(Json::objectValue);
}


pyne::Material::Material(pyne::comp_map cm, double m, std::string s, double apm,
                         Json::Value attributes)
{
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  name = s;
  atoms_per_mol = apm;
  attrs = attributes;
  if (!comp.empty()) 
    norm_comp();
};



pyne::Material::Material(char * fchar, double m, std::string s, double apm,
                         Json::Value attributes)
{
  // Initializes the mass stream based on an isotopic composition file with a (char *) name.
  mass = m;
  name = s;
  atoms_per_mol = apm;
  attrs = attributes;

  // Check that the file is there
  std::string filename (fchar);
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename);
  else
    from_text(filename);
};


pyne::Material::Material(std::string filename, double m, std::string s, double apm,
                         Json::Value attributes)
{
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  name = s;
  atoms_per_mol = apm;
  attrs = attributes;

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


pyne::Material::~Material()
{
};



/*--- Method definitions ---*/


std::ostream& operator<<(std::ostream& os, pyne::Material mat)
{
  //print the Mass Stream to stdout
  os << "Material: " << mat.name << "\n";
  os << "\tMass: " << mat.mass << "\n";
  os << "\t---------\n";
  for(pyne::comp_iter i = mat.comp.begin(); i != mat.comp.end(); i++)
  {
    os << "\t" << pyne::nucname::name( i->first ) << "\t" << i->second << "\n";
  };
  return os;
};


void pyne::Material::normalize ()
{
  // normalizes the mass
  mass = 1.0;
};


pyne::comp_map pyne::Material::mult_by_mass()
{
  // bypass calculation if already normalized.
  if (mass == 1.0)
    return comp;
    
  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    cm[i->first] = (i->second) * mass;
  };
  return cm;
};



double pyne::Material::molecular_weight(double apm)
{
  // Calculate the atomic weight of the Material
  double inverseA = 0.0;
  for (pyne::comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++)
    inverseA += (nuc->second) / pyne::atomic_mass(nuc->first);

  if (inverseA == 0.0)
    return inverseA;

  // select the atoms per mol
  double atsperm = 1.0; // default to 1.0
  if (0.0 <= apm)
  {
    atsperm = apm;            // take the function argument, if valid
    if (atoms_per_mol < 0.0)
      atoms_per_mol = apm;     // Store the function argument on class, if class has no value
  }
  else if (0.0 <= atoms_per_mol)
    atsperm = atoms_per_mol;  // select the class's value

  return atsperm / inverseA;
};







/*--- Stub-Stream Computation ---*/

pyne::Material pyne::Material::sub_mat (std::set<int> nucset,  std::string n)
{
  // Grabs a sub-material from this mat based on a set of integers.
  // Integers can either be of zzaaam form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    if ( 0 < nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, n);
};



pyne::Material pyne::Material::sub_mat (std::set<std::string> nucset,  std::string n)
{
  // Grabs a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++)
  {
    iset.insert(pyne::nucname::zzaaam(*i));
  };

  return sub_mat(iset, n);
};



pyne::Material pyne::Material::set_mat (std::set<int> nucset, double value, std::string n)
{
  // Sets a sub-material from this mat based on a set of integers.
  // Integers can either be of zzaaam form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  
  // Add non-set components
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  // Add set component
  for (std::set<int>::iterator nuc = nucset.begin(); nuc != nucset.end(); nuc++)
    cm[*nuc] = value;
  
  return pyne::Material(cm, -1, n);
};



pyne::Material pyne::Material::set_mat (std::set<std::string> nucset, double value, std::string n)
{
  // Sets a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++)
  {
    iset.insert(pyne::nucname::zzaaam(*i));
  };

  return set_mat(iset, value, n);
};




pyne::Material pyne::Material::del_mat (std::set<int> nucset,  std::string n)
{
  // Removes a sub-material from this mat based on a set of integers.
  // Integers can either be of zzaaam form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    // Only add to new comp if not in nucset
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, n);
};



pyne::Material pyne::Material::del_mat (std::set<std::string> nucset,  std::string n)
{
  // Removes a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++)
  {
    iset.insert(pyne::nucname::zzaaam(*i));
  };

  return del_mat(iset, n);
};






pyne::Material pyne::Material::sub_range(int lower, int upper, std::string n)
{
  // Grabs a sub-material from this mat based on a range of integers.
  if (upper < lower)
  {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  };

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, n);
};



pyne::Material pyne::Material::set_range(int lower, int upper, double value, std::string n)
{
  // Sets a sub-material from this mat based on a range of integers.
  if (upper < lower)
  {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  };

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = value;
    else
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, n);
};



pyne::Material pyne::Material::del_range(int lower, int upper, std::string n)
{
  // Removes a sub-material from this mat based on a range of integers.
  if (upper < lower)
  {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  };

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
  {
    if ((upper <= (i->first)) || ((i->first) < lower))
      cm[i->first] = (i->second) * mass;
  };

  return pyne::Material(cm, -1, n);
};










pyne::Material pyne::Material::sub_u (std::string n)
{
  // Returns a material of Uranium that is a submaterial of this one.
  return sub_range(920000, 930000, n);
};



pyne::Material pyne::Material::sub_pu (std::string n)
{
  // Returns a material of Plutonium that is a sub-material of this one.
  return sub_range(940000, 950000, n);
};



pyne::Material pyne::Material::sub_lan (std::string n)
{
  // Returns a material of Lanthanides that is a sub-material of this one.
  return sub_range(570000, 720000, n);
};



pyne::Material pyne::Material::sub_act (std::string n)
{
  //Returns a material of Actindes that is a sub-material of this one.
  return sub_range(890000, 1040000, n);
};


pyne::Material pyne::Material::sub_tru (std::string n)
{
  // Returns a material of Transuranics that is a sub-material of this one.
  return sub_range(930000, 10000000, n);
};



pyne::Material pyne::Material::sub_ma (std::string n)
{
  // Returns a material of Minor Actinides that is a sub-material of this one.
  return sub_range(930000, 1040000).del_range(940000, 950000, n);
};



pyne::Material pyne::Material::sub_fp (std::string n)
{
  // Returns a material of Fission Products that is a sub-material of this one.
  return sub_range(0, 890000, n);
};




/*--- Atom Frac Functions ---*/

std::map<int, double> pyne::Material::to_atom_frac()
{
  // Returns an atom fraction map from this material's composition

  // the material's molecular weight
  double mat_mw = molecular_weight();

  std::map<int, double> atom_fracs = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_fracs[ci->first] = (ci->second) * mat_mw / pyne::atomic_mass(ci->first);

  return atom_fracs;
};


void pyne::Material::from_atom_frac(std::map<int, double> atom_fracs)
{
  // atom frac must be of the form {nuc: af}, eg, water
  //  80160: 1.0
  //  10010: 2.0

  // clear existing components
  comp.clear();
  atoms_per_mol = 0.0;

  for (std::map<int, double>::iterator afi = atom_fracs.begin(); afi != atom_fracs.end(); afi++)
  {
    comp[afi->first] = (afi->second) * pyne::atomic_mass(afi->first);
    atoms_per_mol += (afi->second);
  };

  norm_comp();
};





/*--- Overloaded Operators ---*/

pyne::Material pyne::Material::operator+ (double y)
{
  // Overloads x + y
  return pyne::Material(comp, mass + y, name);
};



pyne::Material pyne::Material::operator+ (Material y)
{
  // Overloads x + y
  pyne::comp_map cm;
  pyne::comp_map xwgt = mult_by_mass();
  pyne::comp_map ywgt = y.mult_by_mass();

  for (pyne::comp_iter i = xwgt.begin(); i != xwgt.end(); i++)
  {
    if ( 0 < ywgt.count(i->first) )
      cm[i->first] = xwgt[i->first] + ywgt[i->first];
    else
      cm[i->first] = xwgt[i->first];
  };
    
  for (pyne::comp_iter i = ywgt.begin(); i != ywgt.end(); i++)
  {
    if ( 0 == cm.count(i->first) )
      cm[i->first] = ywgt[i->first];			
  };

  return pyne::Material(cm, -1, "");
};



pyne::Material pyne::Material::operator* (double y)
{
  // Overloads x * y
  return pyne::Material(comp, mass * y, name);
};



pyne::Material pyne::Material::operator/ (double y)
{
  // Overloads x / y
  return pyne::Material(comp, mass / y, name);
}

