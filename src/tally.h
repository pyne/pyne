#ifndef PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
#define PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
/// \file tally.h
/// \author Andrew Davis (davisa\@engr.wisc.edu)
///
/// \brief The tally class and helper functions
/// 
/// The tally class is in essesence a structure containing attributes
/// related to tallies

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
    // Tally Constructors
    Tally (); ///< empty constructor

    /// Constructor from passed in vars
    /// \param type the type of tally (flux or current)
    /// \param particle_name the name of the particle type
    /// \param entity the entity id of the tally (eg. surface index, 
    ///          volume number)
    /// \param entity_type (volume or surface)
    /// \param entity_name string identifying the entity
    Tally(std::string type, std::string particle_name, int entity,
  	std::string entity_type, std::string entity_name);
   
    ~Tally (); ///< default destructor

    // Dummy read method wrapper around c style strings
    /// \param filename the filename of the file to read from
    /// \param datapath _name the name of the region where tallies 
    ///          are stored
    void from_hdf5(char * filename, char *datapath, int row = -1);

    // Main read tally method
    /// \param filename the filename of the file to read from
    /// \param datapath _name the name of the region where tallies 
    ///          are stored
    void from_hdf5(std::string filename, std::string datapath, int row = -1);

    // Dummy write method wrapper around c style strings
    /// \param filename the filename of the file to write to
    /// \param datapath _name the name of the region where tallies 
    ///          are to be stored
    void write_hdf5( char * filename, char * datapath);

    // Main write tally method
    /// \param filename the filename of the file to write to
    /// \param datapath _name the name of the region where tallies 
    ///          are to be stored
    void write_hdf5(std::string filename, std::string datapath);

    // fundamental tally variables
    std::string entity_type; // the type of entity (volume,surface)
    std::string entity_name; // the name of the entity (optional)
    std::string particle_name; // particle name string
    std::string tally_type; // type of tally flux or current
    int entity_id; // id number of the entity being tallied upon    
  };

  /// Converts a Tally to a string stream representation.
  std::ostream& operator<< (std::ostream& os, Tally tally);


  /// A stuct for reprensenting fundemental data in a tally
  /// Maybe Useful for HDF5 representations.
  /// following scoptaz's lead here
  typedef struct tally_struct {
    std::string tally_type;
    std::string particle_name;
    int entity_id;
    std::string entity_type;
    std::string entity_name;
  } tally_struct;
  
// End pyne namespace
};

#endif
