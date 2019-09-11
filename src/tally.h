/// \brief The tally class and helper functions.
///
/// The tally class is in essesence a structure containing attributes
/// related to tallies.

#ifndef PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
#define PYNE_IQ4M73STINHJDPRV6KWUZZXOYE

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>

#ifndef PYNE_IS_AMALGAMATED
  #include "h5wrap.h"
  #include "utils.h"
#endif


namespace pyne
{
  class Tally
  {
  public:
    std::map<std::string, std::string> rx2fluka;
    std::map<std::string, std::string> rx2mcnp5;
    std::map<std::string, std::string> rx2mcnp6;

    /// Tally Constructors
    Tally (); /// empty constructor

    /// Constructor from passed in vars
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

    ~Tally (); /// default destructor


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
    void from_hdf5(char * filename, char *datapath, int row = -1);

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
    void write_hdf5( char * filename, char * datapath);

    /// Main write tally method
    /// \param filename the filename of the file to write to
    /// \param datapath _name the name of the region where tallies
    ///          are to be stored
    void write_hdf5(std::string filename, std::string datapath);

    // mcnp tally
    std::string mcnp(int tally_index = 1, std::string mcnp_version = "mcnp5" );

    // Form the mcnp tally line as function of its properties
    std::string form_mcnp_tally(int tally_index, int type,
                                    std::string particle_token, int entity_id,
                                    double entity_size, double normalization);

    // fluka tally
    std::string fluka(std::string unit_number = "-21");

    /// fundamental tally variables
    std::string entity_type; ///< the type of entity (volume,surface)
    std::string entity_name; ///< the name of the entity (optional)
    std::string particle_name; ///< particle name string
    std::string tally_type; ///< type of tally flux or current
    std::string tally_name; ///< name of the tally
    int entity_id; ///< id number of the entity being tallied upon
    double entity_size; ///< the physical size of the entity
    double normalization; ///< the tally normalization
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
    const char * particle_name;
    const char * entity_name;
    const char * tally_name;
    double entity_size;
    double normalization;
  } tally_struct;

// End pyne namespace
}

#endif
