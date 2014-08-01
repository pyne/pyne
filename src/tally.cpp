// Tally.cpp
// Central Tally Class
// -- Andrew Davis

#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
  #include "tally.h"
#endif

enum entity_type {VOLUME,SURFACE}; // Enumeration for entity types
enum tally_type  {FLUX,CURRENT};   // Enumeration for tally types


/***************************/
/*** Protected Functions ***/
/***************************/

/// there are no protected functions currently
/// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Tally::Tally()
{
  // Empty Tally Constructor
  tally_type = "";
  particle_name = "";
  entity_id = -1;
  entity_type = "";
  entity_name = "";
}

pyne::Tally::Tally(std::string type, std::string part_name, 
		   int ent, std::string ent_type, 
		   std::string ent_name)
{
  // Empty Tally Constructor
  tally_type = type;
  particle_name = part_name;
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
}

// Destructor
pyne::Tally::~Tally()
{
};

/*--- Method definitions ---*/

//
void pyne::Tally::from_hdf5(char * filename, char *datapath, int row) 
{
  std::string fname(filename);
  std::string dpath(datapath);
  from_hdf5(fname,dpath,row);
}

//
void pyne::Tally::from_hdf5(std::string filename, std::string datapath, int row) 
{
  
  // check for file existence
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());
  if (!is_h5)
    throw h5wrap::FileNotHDF5(filename);
 
}

// Dummy Wrapper around C Style Functions
void pyne::Tally::write_hdf5( char * filename, char * datapath) {
  std::string fname(filename);
  std::string groupname(datapath);
  write_hdf5(fname,groupname);
}

// Appends Tally object to dataset if file & datapath already exists
// if file exists & data path doesnt creates new datapath, 
// otherwise creates new file
void pyne::Tally::write_hdf5(std::string filename, std::string datapath) {

  tally_t tally_data;

  // setup the data to write
  tally_data.entity_id = entity_id;
  // entity type
  if (tally_type.find("Volume") != std::string::npos)
    tally_data.entity_type = VOLUME;
  else if (tally_type.find("Surface") != std::string::npos)
    tally_data.entity_type = SURFACE;

  // tally kind
  if (tally_type.find("Flux") != std::string::npos)
    tally_data.entity_type = FLUX;
  else if (tally_type.find("Current") != std::string::npos)
    tally_data.entity_type = CURRENT;

  // entity id
  tally_data.entity_id = entity_id;
  // entity_name
  tally_data.entity_name = entity_name.c_str();
  // particle name
  tally_data.particle_name = particle_name.c_str();

  
  // check for file existence
  bool is_exist = pyne::file_exists(filename);
    // create new file

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());
  if ( is_exist && !is_h5)
    throw h5wrap::FileNotHDF5(filename);

  if ( !is_exist ) { // is a new file        
    hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // enable chunking 
    hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
    // set chunk size
    hsize_t chunk_dimensions[1]={2};
    herr_t status = H5Pset_chunk(prop,1,chunk_dims);
    
    // allow varaible length strings
    hid_t strtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (strtype, H5T_VARIABLE);

    // Create the compound datatype for memory.
    hid_t memtype = H5Tcreate (H5T_COMPOUND, sizeof (tally_t));
    status = H5Tinsert (memtype, "Entity ID",
		HOFFSET (tally_t, entity_id), H5T_NATIVE_INT);
    status = H5Tinsert (memtype, "Entity Type",
                HOFFSET (tally_t, entity_type), H5T_NATIVE_INT);
    status = H5Tinsert (memtype, "Tally Name", HOFFSET (tally_t, tally_name),
                strtype);
    status = H5Tinsert (memtype, "Entity Size",
		HOFFSET (tally_t, entity_size), H5T_NATIVE_DOUBLE);

    // Create the compound datatype for the file.
    hid_t filetype = H5Tcreate (H5T_COMPOUND, 8 + 8 + sizeof (hvl_t) + 8);
    status = H5Tinsert (filetype, "Entity ID", 0, H5T_STD_I64BE);
    status = H5Tinsert (filetype, "Entity Type", 8, H5T_STD_I64BE);
    status = H5Tinsert (filetype, "Tally Name", 8+8, strtype);
    status = H5Tinsert (filetype, "Entity Size", 8+8 + sizeof (hvl_t),
			H5T_IEEE_F64BE);

    // only ever let 1 tally object be added
    hsize_t dims[1] = {1};  
    // Create dataspace.  Setting maximum size to NULL sets the maximum
    hid_t space = H5Screate_simple (1, dims, NULL);

    // Create the dataset and write the compound data to it.
    hid_t dset = H5Dcreate2 (file, datapath.c_str(), filetype, space, H5P_DEFAULT, prop,
			     H5P_DEFAULT);
    status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);


    // close the data sets
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Fclose (file);
  
  }  else if ( is_exist && is_h5 ) {// already exists and is an hdf file
     // then we append the data to the end
 
    // Open file and dataset.
    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t dset = H5Dopen2 (file, datapath.c_str(), H5P_DEFAULT);

    // Get dataspace and allocate memory for read buffer.
    hid_t space = H5Dget_space (dset);
    hsize_t dims[1]; // for length of dataset 

    // get the length of the dataset
    int ndims = H5Sget_simple_extent_dims (space, dims, NULL);

    // determine if chunked
    hid_t prop = H5Dget_create_plist(dset);
    if(H5D_CHUNKED == H5Pget_layout(prop))
      rank_chunk = H5Pget_chunk(prop,rank,chunk_dimsr);

    // allocate memory for data from file
    read_data = new tally_t[dims[0]];

    // Create variable-length string datatype.
    strtype = H5Tcopy (H5T_C_S1);
    status  = H5Tset_size (strtype, H5T_VARIABLE);

    // Create the compound datatype for memory.
    memtype = H5Tcreate (H5T_COMPOUND, sizeof (tally_t));
    status = H5Tinsert (memtype, "Entity ID",
                HOFFSET (tally_t, entity_id), H5T_NATIVE_INT);
    status = H5Tinsert (memtype, "Entity Type",
                HOFFSET (tally_t, entity_type), H5T_NATIVE_INT);
    status = H5Tinsert (memtype, "Tally Name", HOFFSET (tally_t, tally_name),
                strtype);
    status = H5Tinsert (memtype, "Entity Size",
		HOFFSET (tally_t, entity_size), H5T_NATIVE_DOUBLE);

    /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */
    filetype = H5Tcreate (H5T_COMPOUND, 8 + 8 + sizeof (hvl_t) + 8);
    status = H5Tinsert (filetype, "Entity ID", 0, H5T_STD_I64BE);
    status = H5Tinsert (filetype, "Entity Type", 8, H5T_STD_I64BE);
    status = H5Tinsert (filetype, "Tally Name", 8+8, strtype);
    status = H5Tinsert (filetype, "Entity Size", 8+8 + sizeof (hvl_t),
			H5T_IEEE_F64BE);

    // Read the data.
    status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);

    // resize dims
    dims[0] += 1;

    // Extend the dataset
    status = H5Dextend(dset,dims);
    hid_t filespace = H5Dget_space(dset);
    // calculate the existing offset
    hsize_t     offset[1] = {dims[0] - 1};  

    // select hyerslab
    hsize_t new_length[1] = {1};
    status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET,offset , NULL,
                                  new_length, NULL);

    // create dataspace for new data
    space = H5Screate_simple(1,new_length, NULL);

    // Write the dataset to memory
    status = H5Dwrite (dset, memtype, space, filespace, H5P_DEFAULT, tally_data);
    std::cout << "status = " << status << std::endl;

    // tidy up
    status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, read_data);
    delete[] read_data;
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (memtype);
    status = H5Tclose (strtype);
    status = H5Fclose (file);
   
  }
}

std::ostream& operator<<(std::ostream& os, pyne::Tally tal)
{
  //print the Tally to ostream
  os << "\t---------\n";
  os << "\t Tallying " << tal.particle_name << " " << tal.tally_type << "\n";
  os << "\t in/on " << tal.entity_type << " " << tal.entity_id << "\n";
  return os;
};
