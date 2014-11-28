#include <iostream>

#define OUR_NAME_TAG_SIZE 128

#include "MBCore.hpp"
#include "MBTagConventions.hpp"

#define CHKERR(rval,msg)  if (MB_SUCCESS != rval) { std::cerr << msg << std::endl; return rval;}
#define CHKERR1(rval,msg,data)  if (MB_SUCCESS != rval) { std::cerr << msg << data << std::endl; return rval;}

MBInterface *MBI() 
{
    static MBCore instance;
    return &instance;
}

void help_msg()
{
  std::cerr <<
    "Usage: dagmc_add_density <in_filename> <volume_id> <density> <out_filename>" << std::endl << 
    std::endl <<
    "\tin_filename   is the H5M file to load" << std::endl <<
    "\tvolume_id     is the volume for which the density will be set (integer)" << std::endl <<
    "\tdensity       is the value of the density which will be set (double)" << std::endl <<
    "\tout_filename  is the H5M file to write the modified geometry" << std::endl;
}

// create a meshset of all volumes to narrow GLOBAL_ID search
MBErrorCode create_set_of_vols(MBEntityHandle& set_of_vols)
{
  MBErrorCode rval;

  // get the handle for the GEOM_DIMENSION tag
  MBTag geom_tag;
  rval = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag);
  CHKERR(rval,"Failed to find GEOM_DIMENSION tag.");

  // use the geom dim tag to the meshsets for all entities of dimension=3
  int geom_dim = 3;
  MBRange vols;
  const void* const geom_dim_search[] = {&geom_dim};
  rval = MBI()->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag,
                                             geom_dim_search, 1, vols);
  CHKERR1(rval,"Failed to find entity sets of dimension ",geom_dim);
    
  // generate a new meshset
  rval = MBI()->create_meshset( MESHSET_SET, set_of_vols);
  CHKERR(rval,"Failed to create meshset.");
  
  // insert volumes into this meshset
  rval = MBI()->add_entities(set_of_vols,vols);
  CHKERR(rval,"Failed to add volume entities to set of volumes.");

  return MB_SUCCESS;
}

MBErrorCode find_vol_with_id(MBEntityHandle set_of_vols, int find_vol_id, MBEntityHandle& vol)
{
  MBErrorCode rval;

  // get handle for querying the GLOBAL_ID
  MBTag id_tag;
  rval = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag);
  CHKERR(rval,"Failed to find GLOBAL_ID tag.");

  // Query the GLOBAL_ID tag to find the meshset for the volume
  const void* const vol_id_search[] = {&find_vol_id};
  MBRange vols;
  rval = MBI()->get_entities_by_type_and_tag( set_of_vols, MBENTITYSET, &id_tag,
                                              vol_id_search, 1, vols);
  CHKERR(rval,"Failed to get the requested volume");

  // there should only be one volume with a given ID
  if (vols.size() > 1)
    {
      std::cout << "Found multiple volumes with ID: " << find_vol_id << std::endl
                << "Using first volume." << std::endl;
    }

  vol = *vols.begin();

  return MB_SUCCESS;

}

int main(int argc, char **argv) 
{
  MBErrorCode rval;
  
  if (argc < 5)
    {
      help_msg();
      return 1;
    }

  // Load a file specified on the command line
  char* in_filename = argv[1];
  // identify a volume ID on the command line
  int find_vol_id = atoi(argv[2]);
  // specify a density on the command line
  double density = atof(argv[3]);
  // specify an output filename on the command line
  char* out_filename = argv[4];

  // Load H5M file
  rval = MBI()->load_file(in_filename);
  CHKERR1(rval,"Failed to load file: ",in_filename);

  // create a meshset containing all the volume meshsets
  // this will narrow the search for GLOBAL_ID below
  MBEntityHandle set_of_vols;
  rval = create_set_of_vols(set_of_vols);

  // Find the volume in question
  MBEntityHandle vol;
  rval = find_vol_with_id(set_of_vols,find_vol_id,vol);

  // Create a new MOAB Tag
  MBTag density_tag;
  char density_tag_name[] = "density";
  rval = MBI()->tag_get_handle( density_tag_name, 1, MB_TYPE_DOUBLE, density_tag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHKERR(rval,"Failed to create density tag!");

  // tag the volume with the density
  rval = MBI()->tag_set_data(density_tag, &vol, 1, &density);
  CHKERR(rval,"Failed to set the density on the requested volume.");

  double check_density;
  rval = MBI()->tag_get_data(density_tag, &vol, 1, &check_density);
  CHKERR(rval,"Failed to get the density.");
  std::cout << "Retrieved density on volume id " << find_vol_id << " = " << check_density << std::endl;

  // create new string tag
  MBTag name_tag;
  char name_tag_name [OUR_NAME_TAG_SIZE] = "name\0";

  rval = MBI()->tag_get_handle( name_tag_name, OUR_NAME_TAG_SIZE, MB_TYPE_OPAQUE, name_tag,
				moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHKERR(rval,"Failed to create name tag");

  // set the value
  char* name_tag_value = "My name is bob\0";
  rval = MBI()->tag_set_data(name_tag,&vol,1,&name_tag_value);
  CHKERR(rval,"Failed to set the entity name");

  // check it can be retrieved
  char* name_tag_return;
  rval = MBI()->tag_get_data(name_tag, &vol, 1, &name_tag_return);
  CHKERR(rval,"Failed to get the density.");
  std::cout << "Retrieved density on volume id " << find_vol_id << " = " << name_tag_return << std::endl;


  // save the mesh
  rval = MBI()->write_mesh(out_filename);
  CHKERR1(rval,"Failed to write file: ",out_filename);

  return MB_SUCCESS;  
}


