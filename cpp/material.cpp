// Material.cpp
// The very central Material class
// -- Anthony Scopatz

#include "material.h"


// h5wrap template
template double h5wrap::get_array_index(H5::DataSet *, int, H5::DataType);



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






void pyne::Material::_load_comp_protocol0(H5::H5File * db, std::string datapath, int row)
{
  H5::Group matgroup = (*db).openGroup(datapath);
  H5::DataSet nucset;

  double nucvalue;
  hsize_t matG = matgroup.getNumObjs();

  // Iterate over datasets in the group.
  for (int matg = 0; matg < matG; matg++)
  {
    std::string nuckey = matgroup.getObjnameByIdx(matg);
    nucset = matgroup.openDataSet(nuckey);
    nucvalue = h5wrap::get_array_index<double>(&nucset, row);

    if (nuckey == "Mass" || nuckey == "MASS" || nuckey == "mass")
      mass = nucvalue;
    else
      comp[pyne::nucname::zzaaam(nuckey)] = nucvalue;

    nucset.close();
  };

  // Set meta data
  name = datapath.substr(datapath.rfind("/")+1, datapath.length());
  atoms_per_mol = -1.0;
};



void pyne::Material::_load_comp_protocol1(H5::H5File * db, std::string datapath, int row)
{
  H5::DataSet data_set = (*db).openDataSet(datapath);

  hsize_t data_offset[1] = {row};
  if (row < 0)
  {
    // Handle negative row indecies
    H5::DataSpace data_space = data_set.getSpace();
    hsize_t data_dims[1];
    int data_rank = data_space.getSimpleExtentDims(data_dims);
    data_offset[0] += data_dims[0];
  };

  // Grab the nucpath
  void * npath;
  H5::Attribute nuc_attr = data_set.openAttribute("nucpath");
  hsize_t nuc_attr_len = nuc_attr.getStorageSize() / sizeof(char);
  H5::StrType nuc_attr_type(0, nuc_attr_len);
  nuc_attr.read(nuc_attr_type, npath);
  std::string nucpath ((char *) npath, ((char *) npath)+nuc_attr_len);

  // Grab the nuclides
  std::vector<int> nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5::PredType::NATIVE_INT);
  int nuc_size = nuclides.size();
  hsize_t nuc_dims[1] = {nuc_size};

  // Get the data hyperslab
  H5::DataSpace data_hyperslab = data_set.getSpace();
  hsize_t data_count[1] = {1};
  data_hyperslab.selectHyperslab(H5S_SELECT_SET, data_count, data_offset);

  // Get memory space for writing
  H5::DataSpace mem_space (1, data_count);

  // Get material type
  size_t material_struct_size = sizeof(pyne::material_struct) + sizeof(double)*(nuc_size);
  H5::CompType data_desc(material_struct_size);
  H5::ArrayType comp_values_array_type (H5::PredType::NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  data_desc.insertMember("name", HOFFSET(pyne::material_struct, name), H5::StrType(0, 20));
  data_desc.insertMember("mass", HOFFSET(pyne::material_struct, mass), H5::PredType::NATIVE_DOUBLE);
  data_desc.insertMember("atoms_per_mol", HOFFSET(pyne::material_struct, atoms_per_mol), H5::PredType::NATIVE_DOUBLE);
  data_desc.insertMember("comp", HOFFSET(pyne::material_struct, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_struct * mat_data = (material_struct *) malloc(material_struct_size);

  // Finally, get data and put in on this instance
  data_set.read(mat_data, data_desc, mem_space, data_hyperslab);

  name = std::string((*mat_data).name);
  mass = (*mat_data).mass;
  atoms_per_mol = (*mat_data).atoms_per_mol;
  for (int i = 0; i < nuc_size; i++)
    comp[nuclides[i]] = (double) (*mat_data).comp[i];

  free(mat_data);
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
  //H5::Exception::dontPrint();

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (!isH5)
    throw h5wrap::FileNotHDF5(filename);

  // Open the database
  H5::H5File db (filename, H5F_ACC_RDONLY);

  bool datapath_exists = h5wrap::path_exists(&db, datapath);
  if (!datapath_exists)
    throw h5wrap::PathNotFound(filename, datapath);

  // Clear current content
  comp.clear();

  // Load via various protocols
  if (protocol == 0)
    _load_comp_protocol0(&db, datapath, row);
  else if (protocol == 1)
    _load_comp_protocol1(&db, datapath, row);
  else
    throw pyne::MaterialProtocolError();

  // Close the database
  db.close();

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



void pyne::Material::write_hdf5(std::string filename, std::string datapath, std::string nucpath, float row, int chunksize)
{
  // Turn off annoying HDF5 errors
  H5::Exception::dontPrint();

  // Create new/open datafile.
  H5::H5File db;
  if (pyne::file_exists(filename))
  {
    bool isH5 = H5::H5File::isHdf5(filename);
    if (!isH5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5::H5File(filename, H5F_ACC_RDWR);
  }
  else
    db = H5::H5File(filename, H5F_ACC_TRUNC);

  //
  // Read in nuclist if available, write it out if not
  //
  bool nucpath_exists = h5wrap::path_exists(&db, nucpath);
  std::vector<int> nuclides;
  int nuc_size;
  hsize_t nuc_dims[1];
  
  if (nucpath_exists)
  {
    nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(&db, nucpath, H5::PredType::NATIVE_INT);
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
    H5::DataSpace nuc_space(1, nuc_dims);
    H5::DataSet nuc_set = db.createDataSet(nucpath, H5::PredType::NATIVE_INT, nuc_space);
    nuc_set.write(nuc_data, H5::PredType::NATIVE_INT);
    db.flush(H5F_SCOPE_GLOBAL);
  };



  //
  // Write out to the file
  //
  H5::DataSet data_set;
  H5::DataSpace data_space, data_hyperslab;
  int data_rank = 1;
  hsize_t data_dims[1] = {1};
  hsize_t data_max_dims[1] = {H5S_UNLIMITED};
  hsize_t data_offset[1] = {0};

  size_t material_struct_size = sizeof(pyne::material_struct) + sizeof(double)*(nuc_size);
  H5::CompType data_desc(material_struct_size);
  H5::ArrayType comp_values_array_type (H5::PredType::NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  data_desc.insertMember("name", HOFFSET(pyne::material_struct, name), H5::StrType(0, 20));
  data_desc.insertMember("mass", HOFFSET(pyne::material_struct, mass), H5::PredType::NATIVE_DOUBLE);
  data_desc.insertMember("atoms_per_mol", HOFFSET(pyne::material_struct, atoms_per_mol), H5::PredType::NATIVE_DOUBLE);
  data_desc.insertMember("comp", HOFFSET(pyne::material_struct, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_struct * mat_data  = (material_struct *) malloc(material_struct_size);
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
  bool datapath_exists = h5wrap::path_exists(&db, datapath);
  if (datapath_exists)
  {
    data_set = db.openDataSet(datapath);
    data_space = data_set.getSpace();
    data_rank = data_space.getSimpleExtentDims(data_dims, data_max_dims);

    // Determine the row size.
    int row_num = (int) row;

    if (signbit(row))
      row_num = data_dims[0] + row;  // careful, row is negative

    if (data_dims[0] <= row_num)
    {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      data_set.extend(data_dims);
    }
    else if (data_dims[0] < 0)
      throw h5wrap::HDF5BoundsError();

    data_offset[0] = row_num;
  }
  else
  {
    // Get full space
    data_space = H5::DataSpace(1, data_dims, data_max_dims);

    // Make data set properties to enable chunking
    H5::DSetCreatPropList data_set_params;
    hsize_t chunk_dims[1] ={chunksize}; 
    data_set_params.setChunk(1, chunk_dims);

    material_struct * data_fill_value  = (material_struct *) malloc(material_struct_size);
    for (int i=0; i < 20; i++)
      (*data_fill_value).name[i] = NULL;
    (*data_fill_value).mass = -1.0;
    (*data_fill_value).atoms_per_mol = -1.0;
    for (int n = 0; n != nuc_size; n++)
      (*data_fill_value).comp[n] = 0.0;
    data_set_params.setFillValue(data_desc, &data_fill_value);

    // Create the data set
    data_set = db.createDataSet(datapath, data_desc, data_space, data_set_params);
    data_set.extend(data_dims);

    // Add attribute pointing to nuc path
    H5::StrType nuc_attr_type(0, nucpath.length());
    H5::DataSpace nuc_attr_space(H5S_SCALAR);
    H5::Attribute nuc_attr = data_set.createAttribute("nucpath", nuc_attr_type, nuc_attr_space);
    nuc_attr.write(nuc_attr_type, nucpath);

    // Remember to de-allocate
    free(data_fill_value);
  };

  // Get the data hyperslab
  data_hyperslab = data_set.getSpace();
  hsize_t data_count[1] = {1};
  data_hyperslab.selectHyperslab(H5S_SELECT_SET, data_count, data_offset);

  // Get a memory space for writing
  H5::DataSpace mem_space (1, data_count, data_max_dims);

  // Write the row...
  data_set.write(mat_data, data_desc, mem_space, data_hyperslab);

  // Close out the HDF5 file
  db.close();

  // Remember the milk!  
  // ...by which I mean to deallocate
  free(mat_data);
};











void pyne::Material::from_text(char * fchar)
{
  std::string filename (fchar);
  from_text(filename);
};


void pyne::Material::from_text (std::string filename)
{
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // New filestream
  std::ifstream f;
  f.open(filename.c_str());

  // Read in
  while ( !f.eof() )
  {
    std::string isostr, wgtstr;
    f >> isostr;
    if (0 < isostr.length())
    {
      f >> wgtstr;
      comp[pyne::nucname::zzaaam(isostr)] = pyne::to_dbl(wgtstr);
    };
  };

  f.close();
  norm_comp();
};



/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/

pyne::Material::Material()
{
  // Empty Material constructor
  mass = -1.0;
  name = "";
  atoms_per_mol = -1.0;
}


pyne::Material::Material(pyne::comp_map cm, double m, std::string s, double apm)
{
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  name = s;
  atoms_per_mol = apm;
  norm_comp();
};



pyne::Material::Material(char * fchar, double m, std::string s, double apm)
{
  // Initializes the mass stream based on an isotopic composition file with a (char *) name.
  mass = m;
  name = s;
  atoms_per_mol = apm;

  // Check that the file is there
  std::string filename (fchar);
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (isH5)
    from_hdf5(filename);
  else
    from_text(filename);
};


pyne::Material::Material(std::string filename, double m, std::string s, double apm)
{
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  name = s;
  atoms_per_mol = apm;

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (isH5)
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
    inverseA += (nuc->second) / pyne::nuc_weight(nuc->first);

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
    atom_fracs[ci->first] = (ci->second) * mat_mw / pyne::nuc_weight(ci->first);

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
    comp[afi->first] = (afi->second) * pyne::nuc_weight(afi->first);
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

