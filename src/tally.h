/// \brief The tally class and helper functions.
///
/// The tally class is in essesence a structure containing attributes
/// related to tallies.

#ifndef PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
#define PYNE_IQ4M73STINHJDPRV6KWUZZXOYE

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#ifndef PYNE_IS_AMALGAMATED
  #include "h5wrap.h"
  #include "utils.h"
#endif


namespace pyne
{
  // simplifying implementation...
  using std::vector;
  // default empty parameters for constructors
  const vector<double> null_v_dbl;
  const vector<int> null_v_int;
  
  class Tally 
  {
  public:
    std::map<std::string, std::string> rx2fluka;
    std::map<std::string, std::string> rx2mcnp5;
    std::map<std::string, std::string> rx2mcnp6;

    /// Tally Constructors
    Tally (); /// empty constructor

    /// Constructor for regular tally
    /// \param type the type of tally (flux or current)
    /// \param particle_name the name of the particle type
    /// \param entity the entity id of the tally (eg. surface index,
    ///          volume number)
    /// \param entity_type (volume or surface)
    /// \param entity_name string identifying the entity
    /// \param tally_name string identifying the tally
    /// \param entity_size the physical size of the tally volume
    /// \param normalization the number required to normalize your tally
    Tally(std::string type, std::string particle_name, int entity,
          std::string entity_type, std::string entity_name,
          std::string tally_name = "", double entity_size = 0.0, 
          double normalization = 1.0);
    
    /// Constructor for regular tally
    /// \param type the type of tally (flux or current)
    /// \param particles_name vector of the name of the particle types
    /// \param entity the entity id of the tally (eg. surface index,
    ///          volume number)
    /// \param entity_type (volume or surface)
    /// \param entity_name string identifying the entity
    /// \param tally_name string identifying the tally
    /// \param entity_size the physical size of the tally volume
    /// \param normalization the number required to normalize your tally
    Tally(std::string type, std::vector<std::string> part_names, int entity,
          std::string entity_type, std::string entity_name,
          std::string tally_name = "", double entity_size = 0.0, 
          double normalization = 1.0);
 
    /// Constructor for mesh tally
    /// \param particle_name the name of the particle type
    /// \param entity_geom (xyz, cylinder)
    /// \param orgn double[3] origin of the meshtally
    /// \param i_mesh, j_mesh, k_mesh vector<double> coordinate of the coarse meshes
    /// \param i_ints, j_ints, k_ints vector<int> coordinate of the coarse meshes
    /// \param e_bounds vector<double values of the coearse mesh points in energy in MeV
    /// \param e_ints vector<int> number of fine mesh points within each corresponding
    ///        coarse mesh in energy
    /// \param axs vector<double> vector giving the direction of the axis of the
    ///        cylindrical mesh
    /// \param vec vector<double> vector defining along with axs the plane for Theta=0
    /// \param tally_name string identifying the tally
    /// \param norm the number required to normalize your tally
    Tally(std::string particle_name, std::string entity_geom, vector<double> orgn,
    vector<double> i_mesh, vector<double> j_mesh, vector<double> k_mesh,
    vector<int> i_ints = null_v_int, 
    vector<int> j_ints = null_v_int, 
    vector<int> k_ints = null_v_int, 
    vector<double> e_bounds = null_v_dbl, vector<int> e_ints = null_v_int,
    vector<double> axs =null_v_dbl, vector<double> vec = null_v_dbl,
    std::string tally_name = "", double norm = 1.0);

    ~Tally();  /// default destructor
 
    // Create hdf5 datatable for tallies
    hid_t create_dataspace(hid_t file, std::string datapath);

    // create filetype of data
    hid_t create_filetype();

    // create memory type for tally struct
    hid_t create_memtype();

    /// Dummy read method wrapper around c style strings
    /// \param filename the filename of the file to read from
    /// \param datapath _name the name of the region where tallies
    ///          are stored
    /// \param row  the array index of data to access
    void from_hdf5(char* filename, char* datapath, int row = -1);

    /// Main read tally method
    /// \param filename the filename of the file to read from
    /// \param datapath _name the name of the region where tallies
    ///          are stored
    /// \param row  the array index of data to access
    void from_hdf5(std::string filename, std::string datapath, int row = -1);

    /// Dummy write method wrapper around c style strings
    /// \param filename the filename of the file to write to
    /// \param datapath _name the name of the region where tallies
    ///          are to be stored
    void write_hdf5(char* filename, char* datapath);

    /// Main write tally method
    /// \param filename the filename of the file to write to
    /// \param datapath _name the name of the region where tallies
    ///          are to be stored
    void write_hdf5(std::string filename, std::string datapath);

    // mcnp tally
    std::string mcnp(int tally_index = 1, std::string mcnp_version = "mcnp5",
                     std::string out = "");

    template<typename T> bool is_zero(T vect);

    // Form the mcnp tally line as function of its properties
    std::string form_mcnp_tally(int tally_index, int type,
                                std::string particle_token, int entity_id,
                                double entity_size, double normalization);
    // Form the mesh tally line as function of its properties
    std::string form_mcnp_meshtally(
        int tally_index, std::string particle_token, std::string entity_geometry,
        std::vector<double> axs, std::vector<double> vec,
        std::vector<double> origin, std::vector<double> meshes[3],
        std::vector<int> ints[3], std::vector<double> e_bounds,
        std::vector<int> e_ints, std::string out);
    
    // fluka tally
    std::string fluka(std::string unit_number = "-21");

    /// fundamental Tally variables
    /// Volume/Surface Tally variables
    std::string entity_type;   ///< the type of entity (volume,surface)
    std::string entity_name;   ///< the name of the entity (optional)
    std::vector<std::string> particle_names; ///< particle names string
    std::string tally_type;    ///< type of tally flux or current
    std::string tally_name;    ///< name of the tally 
    int entity_id;             ///< id number of the entity being tallied upon    
    double entity_size;        ///< the physical size of the entity 
    double normalization;      ///< the tally normalization 

    /// Mesh Tally variables
    std::string entity_geometry;   ///< Shape of the Mesh, XYZ or Cylindrical
    std::vector<double> origin;    ///< Origin of the Mesh Tally
    std::vector<double> vec;       ///< defining with axs the place for Theta=0 (Cyl Mesh)
    std::vector<double> axs;       ///< dircetion of the axis of the cylindrical Mesh
    std::vector<double> meshes[3]; ///< mesh interval coordinates
    std::vector<int> ints[3];      ///< Bin per mesh interval
    std::vector<double> e_bounds;  ///< Energy Mesh
    std::vector<int> e_ints;       ///< Bin per energy
  };

  /// Converts a Tally to a string stream representation.
  std::ostream& operator<< (std::ostream& os, Tally tally);


  /// A stuct for reprensenting fundemental data in a tally
  /// Maybe Useful for HDF5 representations.
  /// following scoptaz's lead here
  typedef struct tally_struct {
    int entity_id;
    int entity_type;
    int tally_type;
    const char* particle_name;
    const char* entity_name;
    const char* tally_name;
    double entity_size;
    double normalization;
  } tally_struct;

// End pyne namespace
}

#endif
