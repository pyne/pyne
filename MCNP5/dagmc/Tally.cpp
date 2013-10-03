// MCNP5/dagmc/Tally.cpp

#include <iostream>

#include "Tally.hpp"
#include "TrackLengthMeshTally.hpp"
#include "KDEMeshTally.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Tally::Tally(const TallyInput& input) 
    : input_data(input), data(NULL)
{
    bool total_energy_bin = true;
    // determine total number of energy bins requested
    int num_energy_bins = input_data.energy_bin_bounds.size();
    assert(num_energy_bins >= 2);

    // turn off total energy bin if only one bin exists
    if (num_energy_bins == 2)
    {
        --num_energy_bins;
       total_energy_bin = false;
    }

    data = new TallyData(num_energy_bins, total_energy_bin); 

}
//---------------------------------------------------------------------------//
// FACTORY METHOD
//
// 3D Flux Tally Types:
//    Other      | Estimator type | Generic Type     Name         Status
//                                 Mesh, Cell, Surf 
//   ------        --------------   -------------   ----------    -----------
//  Unstructured | Track Length   | Mesh Tally   || unstr_track   implemented
//  KDE          | Integral Track | Mesh Tally   || kde_track     KD's Thesis
//  KDE          | SubTrack       | Mesh Tally   || kde_subtrack  implemented
//  KDE          | Collision      | Mesh Tally   || kde_coll      implemented
//               | Collision      | Cell         || coll_cell     testing, not implemented 
//               | Track Length   | Cell         || track_cell    testing, not implemented 
//---------------------------------------------------------------------------//
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
        // This line is not necessary because COLLISION is default estimator.
        KDEMeshTally::Estimator estimator = KDEMeshTally::COLLISION;
        newTally = new KDEMeshTally(input, estimator); 
    }
    else if (input.tally_type == "coll_cell")
    {
        // newTally = new CollCellTally(..)
        std::cout << "Warning: " << input.tally_type
                  << " is not implemented." << std::endl;
    }
    else if (input.tally_type == "track_cell")
    {
        // newTally = new TrackCellTally(..)
        std::cout << "Warning: " << input.tally_type
                  << " is not implemented." << std::endl;
    }
    else 
    {
        std::cout << "Warning: " << input.tally_type
                  << " is not a valid tally type." << std::endl;
    }
         
    return newTally;
}
//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
void Tally::end_history()
{
    data->my_end_history();
}
//---------------------------------------------------------------------------//
void Tally::add_score_to_tally(unsigned int tally_point, double score, unsigned int ebin)
{
    data->my_add_score_to_tally(tally_point, score, ebin);
}

// end of MCNP5/dagmc/Tally.cpp
