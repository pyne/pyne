// Tally.cpp
// Central Tally Class
// -- Andrew Davis

#include <iomanip>
#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
  #include "particle.h"
  #include "tally.h"
#endif

enum entity_type_enum {VOLUME, SURFACE, MESH}; // Enumeration for entity types
enum tally_type_enum {FLUX, CURRENT};  // Enumeration for tally types

const std::string tally_type_enum2string[] = {"Flux", "Current"};
const std::string entity_type_enum2string[] = {"Volume", "Surface", "Mesh"};
const std::string geometry_type_enum2string[] = {"Cartesian", "Cylinder"};

/***************************/
/*** Protected Functions ***/
/***************************/

// there are no protected functions currently
// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Tally::Tally() {
  // Empty Tally Constructor
  tally_type = "";
  entity_id = -1;
  entity_type = "";
  entity_name = "";
  tally_name = "";
  entity_size = -1.0;
  normalization = 1.0;
}

// Default constructor
pyne::Tally::Tally(std::string type, std::string part_name, 
       int ent, std::string ent_type, 
       std::string ent_name, std::string tal_name,
       double size, double norm ) {

  tally_type = type;
  particle_names.push_back(pyne::particle::name(part_name));
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string type, std::vector<std::string> part_names, 
       int ent, std::string ent_type, 
       std::string ent_name, std::string tal_name,
       double size, double norm ) {
  
  tally_type = type;
  particle_names = part_names;
  for (int i = 0; i < particle_names.size(); i++){
    particle_names[i] = pyne::particle::name(particle_names[i]);
  }
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string part_name, std::string ent_geom,
                   std::vector<double> orgn,
                   std::vector<double> mesh_i, std::vector<double> mesh_j, std::vector<double> mesh_k,
                   std::vector<int> ints_i, std::vector<int> ints_j, std::vector<int> ints_k,
                   std::vector<double> e_bounds_, std::vector<int> e_ints_,
                   std::vector<double> axs_, std::vector<double> vec_,
                   std::string tal_name, double norm) {
  // Empty Tally Constructor
  entity_type = "Mesh";
  entity_name = "";
  particle_names.push_back(pyne::particle::name(part_name));
  entity_geometry = ent_geom;
  tally_name = tal_name;
  entity_size = -1;

  origin = orgn;
  vec = vec_;
  axs = axs_;
  meshes[0] = mesh_i;
  meshes[1] = mesh_j;
  meshes[2] = mesh_k;
  ints[0] = ints_i;
  ints[1] = ints_j;
  ints[2] = ints_k;
  e_bounds = e_bounds_;
  e_ints = e_ints_;
  normalization = norm;
}

// Destructor
pyne::Tally::~Tally() {}

/*--- Method definitions ---*/
//
void pyne::Tally::from_hdf5(char * filename, char *datapath, int row) {
  std::string fname(filename);
  std::string dpath(datapath);
  from_hdf5(fname,dpath,row);
}

//
void pyne::Tally::from_hdf5(std::string filename, std::string datapath, 
          int row) { 
  // line of data to acces
  int data_row = row;

  // check for file existence
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());
  if (!is_h5)
    throw h5wrap::FileNotHDF5(filename);

  // Open file and dataset.
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);

  // Get dataspace and allocate memory for read buffer.
  hid_t space = H5Dget_space(dset);
  int rank  = H5Sget_simple_extent_ndims(space);
  hsize_t dims[1]; // for length of dataset 

  // get the length of the dataset
  int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
  
  // determine if chunked
  hid_t prop = H5Dget_create_plist(dset);
  
  hsize_t chunk_dimsr[1];
  int rank_chunk;

  if(H5D_CHUNKED == H5Pget_layout(prop))
    rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);
  
  // allocate memory for data from file
  tally_struct* read_data = new tally_struct[dims[0]];

  // if row number is larger than data set only give last element
  if ( row >= dims[0] )
    data_row = dims[0]-1;
    
  // Create variable-length string datatype.
  hid_t strtype = H5Tcopy(H5T_C_S1);
  int status  = H5Tset_size(strtype, H5T_VARIABLE);
  
  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();
  
  // Create the compound datatype for the file
  hid_t filetype = create_filetype();
  
  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);

  // unpack the data and set values
  entity_id = read_data[data_row].entity_id;
  entity_type = entity_type_enum2string[read_data[data_row].entity_type];
  tally_type = tally_type_enum2string[read_data[data_row].tally_type];
  particle_names = pyne::split_string(read_data[data_row].particle_name, ",");
  tally_name = std::string(read_data[data_row].tally_name);
  entity_name = std::string(read_data[data_row].entity_name);
  entity_size = read_data[data_row].entity_size;
  normalization = read_data[data_row].normalization;

  // close the data sets
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(filetype);
  status = H5Fclose(file);

  // tidy up
  delete[] read_data;
 
}

// Dummy Wrapper around C Style Functions
void pyne::Tally::write_hdf5(char * filename, char * datapath) {
  std::string fname(filename);
  std::string groupname(datapath);
  write_hdf5(fname,groupname);
}

// create filetype
hid_t pyne::Tally::create_filetype() {
  herr_t status;  // iostatus

  // create string type
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);

  hid_t filetype = H5Tcreate(H5T_COMPOUND, 8 + 8 + 8 + 
           (3*sizeof(hvl_t)) + 8 + 8);
  status = H5Tinsert(filetype, "entity_id", 0, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "entity_type", 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "tally_type", 8 + 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "particle_name", 8 + 8 + 8, strtype);
  status = H5Tinsert(filetype, "entity_name", 8 + 8 + 8 + 
         sizeof(hvl_t), strtype);
  status = H5Tinsert(filetype, "tally_name", 8 + 8 + 8 + 
         (2*sizeof(hvl_t)) , strtype);
  status = H5Tinsert(filetype, "entity_size", 8 + 8 + 8 + 
         (3*sizeof(hvl_t)), H5T_IEEE_F64BE);
  status = H5Tinsert(filetype, "normalization", 8 + 8 + 8 + 
         (3*sizeof(hvl_t)) + 8, H5T_IEEE_F64BE);
  return filetype;
}

// create memory type 
hid_t pyne::Tally::create_memtype() {
  // iostatus
  herr_t status;

  //
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);
  // Create the compound datatype for memory.
  hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(tally_struct));
  status = H5Tinsert(memtype, "entity_id",
         HOFFSET(tally_struct, entity_id), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "entity_type",
         HOFFSET(tally_struct, entity_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "tally_type",
         HOFFSET(tally_struct, tally_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "particle_name",
         HOFFSET(tally_struct, particle_name), strtype);
  status = H5Tinsert(memtype, "entity_name",HOFFSET(tally_struct, entity_name),
         strtype);
  status = H5Tinsert(memtype, "tally_name",HOFFSET(tally_struct, tally_name),
         strtype);
  status = H5Tinsert(memtype, "entity_size",
         HOFFSET(tally_struct, entity_size), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(memtype, "normalization",
         HOFFSET(tally_struct, normalization), H5T_NATIVE_DOUBLE);
  return memtype;
}

hid_t pyne::Tally::create_dataspace(hid_t file, std::string datapath) {
    // enable chunking 
    hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
    // set chunk size
    hsize_t chunk_dimensions[1]={1};
    herr_t status = H5Pset_chunk(prop, 1, chunk_dimensions);
    
    // allow varaible length strings
    hid_t strtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(strtype, H5T_VARIABLE);

    // Create the compound datatype for memory.
    hid_t memtype = create_memtype();
    
    // Create the compound datatype for the file
    hid_t filetype = create_filetype();
    
    // max dims unlimted
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    // only ever let 1 tally object be added
    hsize_t dims[1] = {1}; 
    // Create dataspace.  Setting maximum size to NULL sets the maximum
    hid_t space = H5Screate_simple(1, dims, max_dims);

    // Create the dataset and write the compound data to it.
    return H5Dcreate2(file, datapath.c_str(), filetype, space, H5P_DEFAULT, prop,
           H5P_DEFAULT);
}

// Appends Tally object to dataset if file & datapath already exists
// if file exists & data path doesnt creates new datapath, 
// otherwise creates new file
void pyne::Tally::write_hdf5(std::string filename, std::string datapath) {

  // turn of annoying hdf5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  tally_struct tally_data[1]; // storage for the tally to add

  // setup the data to write
  tally_data[0].entity_id = entity_id;
  // entity type
  if (entity_type.find("Volume") != std::string::npos)
    tally_data[0].entity_type = VOLUME;
  else if (entity_type.find("Surface") != std::string::npos)
    tally_data[0].entity_type = SURFACE;

  // tally kind
  if (tally_type.find("Flux") != std::string::npos)
    tally_data[0].tally_type = FLUX;
  else if (tally_type.find("Current") != std::string::npos)
    tally_data[0].tally_type = CURRENT;

  // unpack from class to struct array
  std::string particle_name = pyne::join_to_string(particle_names, ",");
  tally_data[0].entity_id = entity_id;
  tally_data[0].entity_name = entity_name.c_str();
  tally_data[0].particle_name = particle_name.c_str();
  tally_data[0].tally_name = tally_name.c_str();
  tally_data[0].entity_size = entity_size;
  tally_data[0].normalization = normalization;

  
  // check for file existence
  bool is_exist = pyne::file_exists(filename);
  // create new file

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());

  if (is_exist && !is_h5)
    throw h5wrap::FileNotHDF5(filename);

  if (!is_exist ) { // is a new file        
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
         H5P_DEFAULT);
    // create a dataspace
    hid_t dset = create_dataspace(file, datapath);

    hid_t memtype = create_memtype();
  
    herr_t status; // iostatus

    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

    // close the data sets
    status = H5Dclose(dset);
    //    status = H5Sclose(space);
    //    status = H5Tclose(filetype);
    status = H5Fclose(file);
  
  }  else if ( is_exist && is_h5 ) {// already exists and is an hdf file
    // then we append the data to the end
    herr_t data_status; // iostatus

    // Open file and dataset.
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // see if path exists
    data_status = H5Gget_objinfo (file, datapath.c_str(), 0, NULL);

    hid_t dset;
    // if fails neet to create dataset
    // still need to check that the datapath exists
    if (data_status != 0) // doesnt exist
      {
  dset = create_dataspace(file,datapath.c_str());
  hid_t memtype = create_memtype();
  herr_t status; // iostatus

  status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

  // close the data sets
  status = H5Dclose(dset);
  status = H5Fclose(file);
      } else {
      
        dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);
      
  // Get dataspace and allocate memory for read buffer.
  hid_t space = H5Dget_space(dset);
  int rank  = H5Sget_simple_extent_ndims(space);
  hsize_t dims[1]; // for length of dataset 
  
  // get the length of the dataset
  int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
  
  // determine if chunked
  hid_t prop = H5Dget_create_plist(dset);

  hsize_t chunk_dimsr[1];
  int rank_chunk;
  if (H5D_CHUNKED == H5Pget_layout(prop))
    rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);
  
  // allocate memory for data from file
  tally_struct* read_data = new tally_struct[dims[0]];
  
  // Create variable-length string datatype.
  hid_t strtype = H5Tcopy(H5T_C_S1);
  int status  = H5Tset_size(strtype, H5T_VARIABLE);
  
  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();
  
  // Create the compound datatype for the file
  hid_t filetype = create_filetype();
  
  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);
      
  // resize dims
  dims[0] += 1;
  
  // Extend the dataset
  status = H5Dextend(dset,dims);
  hid_t filespace = H5Dget_space(dset);
  // calculate the existing offset
  hsize_t offset[1] = {dims[0] - 1};  
  
  // select hyerslab
  hsize_t new_length[1] = {1};
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,offset , NULL,
             new_length, NULL);
  
  // create dataspace for new data
  space = H5Screate_simple(1,new_length, NULL);
  
  // Write the dataset to memory
  status = H5Dwrite(dset, memtype, space, filespace, H5P_DEFAULT, tally_data);
  
  // tidy up
  status = H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, read_data);
  delete[] read_data;
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(memtype);
  status = H5Tclose(strtype);
  status = H5Fclose(file);
    }
  }
}

std::ostream& operator<<(std::ostream& os, pyne::Tally tal) {
  //print the Tally to ostream
  os << "\t---------\n";
  os << "\t Tallying " << pyne::join_to_string(tal.particle_names, ", ") << " " << tal.tally_type << std::endl;
  os << "\t in/on " << tal.entity_type << " " << tal.entity_id << std::endl;
  return os;
}

// Sets string to valid mcnp formatted tally
// Takes mcnp version as arg, like 5 or 6
std::string pyne::Tally::mcnp(int tally_index, std::string mcnp_version,
                              std::string out) {
  std::stringstream output;  // output stream
  std::string particle_token = "";

  if (particle_names.size() == 0) {
      particle_token = "?";
  }

  // particle token
  for (int i = 0; i < particle_names.size(); i++ ){
    if (i > 0)
      particle_token += ",";
    if (mcnp_version.find("mcnp5") != std::string::npos)
      particle_token += pyne::particle::mcnp(particle_names[i]);
    else if (mcnp_version.find("mcnp6") != std::string::npos)
      particle_token += pyne::particle::mcnp6(particle_names[i]);
    else
      particle_token += "?";
  }
  // print out comment line
  output << "C " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(6);

  if (normalization != 1.0) output << std::scientific;
  int tally_id = 0;

  // neednt check entity type
  if (entity_type.find("Surface") != std::string::npos) {
    if (tally_type.find("Current") != std::string::npos) {
      tally_id = 1;
    } else if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 2;
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Volume") != std::string::npos) {
    if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 4;
    } else if (tally_type.find("Current") != std::string::npos) {
      // makes no sense in mcnp
      return "";
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Mesh") != std::string::npos) {
    output << form_mcnp_meshtally(tally_index, particle_token, entity_geometry,
                                  axs, vec, origin, meshes, ints, e_bounds,
                                  e_ints, out);

  } else {
    std::cout << "tally/entity combination makes no sense for MCNP"
              << std::endl;
  }
  // print sd card if area/volume specified
  return output.str();
}

template <typename T>
bool pyne::Tally::is_zero(T vect) {
  int size = sizeof(vect) / sizeof(vect[0]);
  bool result = true;
  for (int i = 0; i < size; i++) result &= (vect[i] == 0);
  return result;
}


// Form the tally line as function of its properties
std::string pyne::Tally::form_mcnp_tally(int tally_index, int type,
                                         std::string particle_token,
                                         int entity_id, double entity_size,
                                         double normalization) {
  std::stringstream tally_stream;  // tally stream
  tally_stream << std::setiosflags(std::ios::fixed) << std::setprecision(6);
  if (normalization != 1.0) tally_stream << std::scientific;

  tally_stream << "F" << tally_index << type << ":" << particle_token << " "
               << entity_id << std::endl;

  if (entity_size > 0.0)
    tally_stream << "SD" << tally_index << type << " " << entity_size
                 << std::endl;

  if (normalization != 1.0)
    tally_stream << "FM" << tally_index << type << " " << normalization
                 << std::endl;

  return tally_stream.str();
}


// Form the mesh tally line as function of its properties
std::string pyne::Tally::form_mcnp_meshtally(
    int tally_index, std::string particle_token, std::string entity_geometry,
    std::vector<double> axs, std::vector<double> vec,
    std::vector<double> origin, std::vector<double> meshes[3],
    std::vector<int> ints[3], std::vector<double> e_bounds,
    std::vector<int> e_ints, std::string out) {
  std::stringstream mtally_stream;
  // indentation block
  std::string indent_block = "          ";

  mtally_stream << "FMESH" << tally_index << "4:" << particle_token << " ";
  mtally_stream << "GEOM=";

  if (entity_geometry.find("Cartesian") != std::string::npos) {
    mtally_stream << "XYZ" << std::endl;
    mtally_stream << indent_block;
  } else if (entity_geometry.find("Cylinder") != std::string::npos) {
    mtally_stream << "CYL" << std::endl;
    if (!is_zero(axs)) {
      mtally_stream << indent_block << "AXS=" << pyne::join_to_string(axs) << std::endl;
    }
    if (!is_zero(vec)) {
      mtally_stream << indent_block << "VEC=" << pyne::join_to_string(vec) << std::endl;
    }
    mtally_stream << indent_block;
  }

  mtally_stream << "ORIGIN=" << pyne::join_to_string(origin) << std::endl;
  std::string dir_name[3] = {"I", "J", "K"};

  for (int j = 0; j < 3; j++) {
    mtally_stream << indent_block;
    mtally_stream << dir_name[j] << "MESH=" << pyne::join_to_string(meshes[j]);
    if (ints[j].size() > 0) {
      mtally_stream << " " << dir_name[j]
                    << "INTS=" << pyne::join_to_string(ints[j]);
    }
    mtally_stream << std::endl;
  }
  if (e_bounds.size() > 0) {
    mtally_stream << indent_block << "EMESH=" << pyne::join_to_string(e_bounds);
  }
  mtally_stream << std::endl;
  if (e_ints.size() > 0) {
    mtally_stream << indent_block << "EINTS=" << pyne::join_to_string(e_ints);
  }
  if (out.size() > 0) {
    mtally_stream << std::endl << indent_block << "OUT=" << out;
  }
  return mtally_stream.str();
}


// Produces valid fluka tally
std::string pyne::Tally::fluka(std::string unit_number) {
  std::stringstream output; // output stream
  
  if (particle_names.size() != 1) {
    std::cout << "Fluka multi-particles tally have not been yet implemented!" << std:: endl;
    exit(1);
  }
  std::string particle_name = particle_names[0];

  // check entity type
  if (entity_type.find("Volume") != std::string::npos) {
    // ok
  }  else if (entity_type.find("Surface") != std::string::npos) {
      std::cout << "Surface tally not valid in FLUKA" << std::endl;
  } else {
      std::cout << "Unknown entity type" << std::endl;
  }

  output << "* " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(1);
  // check tally type
  if (tally_type.find("Flux") != std::string::npos) {
      output << std::setw(10) << std::left  << "USRTRACK";
      output << std::setw(10) << std::right << "     1.0";
      output << std::setw(10) << std::right 
             << pyne::particle::fluka(particle_name);
      output << std::setw(10) << std::right << unit_number;
      output << std::setw(10) << std::right << entity_name;
      if(entity_size > 0.0) {
        output << std::scientific;
        output << std::setprecision(4);
        output << std::setw(10) << std::right << entity_size;
      }
      else
        output << std::setw(10) << std::right << 1.0;

      output << std::setw(10) << std::right << "   1000."; // number of eints
      tally_name.resize(8);
      output << std::setw(8) << std::left 
             << tally_name; // may need to make sure less than 10 chars
      output << std::endl;
      output << std::setw(10) << std::left  << "USRTRACK";
      output << std::setw(10) << std::right << "   10.E1";
      output << std::setw(10) << std::right << "   1.E-3";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(8) << std::left << "       &";      
      // end of usrtrack
  } else if (tally_type.find("Current") != std::string::npos) {
      output << std::setw(10) << std::left  << "USRBDX  ";    
      output << std::setw(10) << std::right << "   110.0";
      output << std::setw(10) << std::right 
             << pyne::particle::fluka(particle_name);
      output << std::setw(10) << std::right << unit_number;
      output << std::setw(10) << std::right << entity_name; // upstream
      output << std::setw(10) << std::right << entity_name; // downstream
      if ( entity_size > 0.0 )
        output << std::setw(10) << std::right << entity_size; // area
      else
        output << std::setw(10) << std::right << 1.0;

      tally_name.resize(8);
      output << std::setw(8) << std::right 
             << tally_name; // may need to make sure less than 10 chars
      output << std::endl;
      output << std::setw(10) << std::left  << "USRBDX  ";    
      output << std::setw(10) << std::right << "  10.0E1";
      output << std::setw(10) << std::right << "     0.0";
      output << std::setw(10) << std::right << "  1000.0"; // number of ints
      output << std::setw(10) << std::right << "12.56637"; // 4pi
      output << std::setw(10) << std::right << "     0.0";
      output << std::setw(10) << std::right 
             << "   240.0"; // number of angular ints
      output << std::setw(8) << std::left << "       &";      
      // end of usrbdx
  } else {
    std::cout << "Unknown tally type" << std::endl;
  }
  return output.str();
}
