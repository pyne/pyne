// MCNP5/dagmc/Tally.cpp
#include <iostream>
#include "Tally.hpp"
#include "TrackLengthMeshTally.hpp"
#include "KDEMeshTally.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Tally::Tally(const TallyInput& input) 
     : input_data(input), total_energy_bin(true)
{    
       // Determine the total number of energy bins requested
       num_energy_bins = input_data.energy_bin_bounds.size();

       assert(num_energy_bins >= 2);
       if (num_energy_bins == 2)
       {
          --num_energy_bins;
          total_energy_bin = false;
       }
}

//---------------------------------------------------------------------------//
// create_tally(..)  Factory Method
//---------------------------------------------------------------------------//

    /**
     *  \brief Factory method for creation of Tally Observers
     *
     *  3D Flux Tally Types:
     *
     *    Other      | Estimator type | Generic Type     Name         Status
     *                                 Mesh, Cell, surf 
     *   ------        --------------   -------------   ----------    -----------
     *  Unstructured | Track Length   | Mesh Tally   || unstr_track   implemented
     *  KDE          | Integral Track | Mesh Tally   || kde_track     KD's Thesis
     *  KDE          | SubTrack       | Mesh Tally   || kde_subtrack  implemented
     *  KDE          | Collision      | Mesh Tally   || kde_coll      implemented
     *               | Collision      | Cell         || coll_cell     testing, not implemented 
     *               | Track Length   | Cell         || track_cell    testing, not implemented 
     *
    */
Tally *Tally::create_tally(const TallyInput& input)
{
      Tally *newTally = NULL;
          
      if (input.tally_type == "unstr_track")
      {
         newTally = new moab::TrackLengthMeshTally(input);
      }
      else if (input.tally_type == "kde_track")
      {
         KDEMeshTally::Estimator estimator = KDEMeshTally::INTEGRAL_TRACK;
         newTally = new KDEMeshTally(input, estimator); 
      }
      else if (input.tally_type == "kde_subtrack")
      {
         KDEMeshTally::Estimator estimator = KDEMeshTally::SUB_TRACK;
         newTally = new KDEMeshTally(input, estimator); 
      }
      else if (input.tally_type == "kde_coll")
      {
         // This line is not necessary because COLLISION is the default estimator.
         KDEMeshTally::Estimator estimator = KDEMeshTally::COLLISION;
         newTally = new KDEMeshTally(input, estimator); 
      }
      else if (input.tally_type == "coll_cell")
      {
         // newTally = new CollCellTally(..)
         std::cout << "Warning: " << input.tally_type << " is not implemented." << std::endl;
      }
      else if (input.tally_type == "track_cell")
      {
         // newTally = new TrackCellTally(..)
         std::cout << "Warning: " << input.tally_type << " is not implemented." << std::endl;
      }
      else 
      {
         std::cout << "Warning: " << input.tally_type << " is not a valid tally type." << std::endl;
      }
         
     return newTally;
}

// end of MCNP5/dagmc/Tally.cpp
