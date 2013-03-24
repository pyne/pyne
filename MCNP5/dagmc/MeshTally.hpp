// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESHTALLY_H
#define DAGMC_MESHTALLY_H

#include <string>
#include <map>
#include <vector>

#include <cassert>

#include "moab/Range.hpp"

//===========================================================================//
/**
 * \struct MeshTallyInput
 * \brief Input data needed to set up a MeshTally
 */
//===========================================================================//
struct MeshTallyInput
{
    /// Typedef for map that stores optional MeshTally input parameters
    typedef std::multimap<std::string, std::string> TallyOptions;

    /// User-specified ID for this MeshTally
    int tally_id;

    // TODO temporary, see if we can remove this from the MeshTallyInput struct later on
    /// The index in the fortran 'fm' array, also used to index arrays in meshtal_funcs.cpp
    int fmesh_index;

    // TODO change pointer to Fortran memory into vector in C++ memory  
    /// Energy bin boundaries defined for all mesh tally points
    //std::vector<double> energy_bin_bounds;
    const double* energy_bin_bounds;

    // TODO temporary, see if this is even needed based on size of vector
    /// The length of energy_bin_boundaries
    int num_ebin_bounds;

    /// If true, add an extra energy bin to tally all energy levels
    bool total_energy_bin;

    /// Optional input parameters requested by user for this MeshTally
    TallyOptions options;
};

// forward declaration
namespace moab{
  class Interface;
}

//===========================================================================//
/**
 * MOAB-based mesh tally class
 */
//===========================================================================//
class MeshTally {

protected:


  MeshTally( const MeshTallyInput& input ):
    fmesh(input)
  {
    ebins = fmesh.num_ebin_bounds;
    if( !fmesh.total_energy_bin ){
      ebins = fmesh.num_ebin_bounds - 1;
    }
    assert( ebins > 0 );
  }

public:
  virtual ~MeshTally(){}
  
  /**
   * Print / write results to the AMeshTally's output file.
   * @param sp_norm The number of source particles, as reported from within mcnp's fortran code.
   * @param fmesh_fact Multiplication factor from fmesh card.  
   */
  virtual void print( double sp_norm, double fmesh_fact) = 0;

  /**
   * Updates tally information when the history of a particle ends.
   */
  virtual void end_history() = 0;

  /**
   * get_tally_data(), get_error_data(), get_scratch_data() : 
   * These functions get pointer to tally data arrays, with length as output parameter.
   * They are used to load and reload tally data by the runtpe and MPI functions.
   * A subclass need not use these arrays to store its data, but unless it does,
   * runtpe and MPI features will not work.  
   */

  virtual double* get_tally_data( int& length ){
    length = tally_data.size();
    return &(tally_data[0]);
  }

  virtual double* get_error_data( int& length ){
    length = error_data.size();
    return &(error_data[0]);
  }

  virtual double* get_scratch_data( int& length ){
    length = temp_tally_data.size();
    return &(temp_tally_data[0]);
  }

  virtual void zero_tally_data( ){
    std::fill( tally_data.begin(), tally_data.end(), 0 );
    std::fill( error_data.begin(), error_data.end(), 0 );
    std::fill( temp_tally_data.begin(), temp_tally_data.end(), 0 );
  }

  int get_fmesh_index() { return fmesh.fmesh_index; }

protected:
  /**
   * Resize data storage arrays to hold a given number of points.
   * Arrays will be resized to the given size * the number of energy bins
   */
  void resize_data_arrays( unsigned int size ){
    tally_data.resize( size * ebins, 0 );
    error_data.resize( size * ebins, 0 );
    temp_tally_data.resize( size * ebins, 0 );
  }

  unsigned int ent_idx( moab::EntityHandle eh ){
    
    unsigned int ret = tally_ents.index( eh );
    assert( ret < tally_ents.size() );
    return ret;

  }

  double& data_ref( std::vector<double>& data, moab::EntityHandle eh, unsigned ebin = 0){
    assert( ebin < ebins );
    return data[ ent_idx(eh)*ebins + ebin ];
  }

  moab::ErrorCode setup_tags( moab::Interface* mbi, const char* prefix="" );

  /// User's MCNP input parameters for this mesh tally
  MeshTallyInput fmesh;

  /// data arrays
  std::vector<double> tally_data, error_data, temp_tally_data;

  /// actual number of energy bins implemented in the data arrays
  unsigned ebins; 

  /// entities on which to compute tally
  moab::Range tally_ents;

  /// Tag arrays
  std::vector< moab::Tag > tally_tags, error_tags; 

};

#endif // DAGMC_MESHTALLY_H

// end of MCNP5/dagmc/MeshTally.hpp
