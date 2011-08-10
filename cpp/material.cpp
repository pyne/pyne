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



void pyne::Material::load_from_hdf5(char * fchar, char * gchar, int row)
{
  std::string filename (fchar);
  std::string groupname (gchar);
  load_from_hdf5(filename, groupname, row);  
};



void pyne::Material::load_from_hdf5(std::string filename, std::string groupname, int row)
{
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (!isH5)
    throw h5wrap::FileNotHDF5(filename);

  H5::Exception::dontPrint();

  H5::H5File matfile (filename, H5F_ACC_RDONLY);
  H5::Group matgroup;

  try
  {
    matgroup = matfile.openGroup(groupname);
  }
  catch (H5::Exception fgerror)
  {
    throw h5wrap::GroupNotFound(filename, groupname);
  }    

  // Clear current content
  comp.clear();

  // Iterate over elements of the group.
  H5::DataSet nucset;
  double nucvalue;
  hsize_t matG = matgroup.getNumObjs();
  for (int matg = 0; matg < matG; matg++)
  {
    std::string nuckey = matgroup.getObjnameByIdx(matg);
    nucset = matgroup.openDataSet(nuckey);
    nucvalue = h5wrap::get_array_index<double>(&nucset, row);

    if (nuckey == "Mass" || nuckey == "MASS" || nuckey == "mass")
      mass = nucvalue;
    else
      comp[nucname::zzaaam(nuckey)] = nucvalue;

    nucset.close();
  };

  // FIXME: Set the material name here. (based on groupname)

  matfile.close();

  norm_comp();
};



void pyne::Material::load_from_text (char * fchar)
{
  std::string filename (fchar);
  load_from_text(filename);
};


void pyne::Material::load_from_text (std::string filename)
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
      comp[nucname::zzaaam(isostr)] = pyne::to_dbl(wgtstr);
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
}


pyne::Material::Material(pyne::comp_map cm, double m, std::string s)
{
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  name = s;
  norm_comp();
};



pyne::Material::Material(char * fchar, double m, std::string s)
{
  // Initializes the mass stream based on an isotopic composition file with a (char *) name.
  mass = m;
  name = s;

  // Check that the file is there
  std::string filename (fchar);
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (isH5)
    load_from_hdf5(filename);
  else
    load_from_text(filename);
};


pyne::Material::Material(std::string filename, double m, std::string s)
{
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  name = s;

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(filename);
  if (isH5)
    load_from_hdf5(filename);
  else
    load_from_text(filename);
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
    os << "\t" << nucname::name( i->first ) << "\t" << i->second << "\n";
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



double pyne::Material::molecular_weight()
{
  // Calculate the atomic weight of the Material
  double inverseA = 0.0;
  for (pyne::comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++)
    inverseA += (nuc->second) / nucname::nuc_weight(nuc->first);

  if (inverseA == 0.0)
    return inverseA;

  return 1.0 / inverseA;
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
    iset.insert(nucname::zzaaam(*i));
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
    iset.insert(nucname::zzaaam(*i));
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
    iset.insert(nucname::zzaaam(*i));
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

