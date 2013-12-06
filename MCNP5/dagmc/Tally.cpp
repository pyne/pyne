// MCNP5/dagmc/Tally.cpp

#include <cassert>
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
// DESTRUCTOR
//---------------------------------------------------------------------------//
Tally::~Tally()
{
   delete data;
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
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void Tally::end_history()
{
    data->end_history();
}
//---------------------------------------------------------------------------//
// PROTECTED INTERFACE
//---------------------------------------------------------------------------//
unsigned int Tally::get_energy_bin(double energy)
{
    unsigned int bin_no = 0;

    // Test:  print out the energy and all the boundaries
    std::cout << energy << std::endl;
    for (std::vector<double>::iterator bdry = input_data.energy_bin_bounds.begin(); bdry != input_data.energy_bin_bounds.end(); ++bdry)
    {
        std::cout << *bdry << ' ';
    }
    std::cout << '\n';

    std::vector<double>::iterator curIt = input_data.energy_bin_bounds.begin();
    double prev = *curIt;
    // Move the iterator to point to the second boundary
    ++curIt;
    for (unsigned int bin = 0; curIt != input_data.energy_bin_bounds.end(); ++curIt, ++bin)
    {
       if (prev <= energy && energy < *curIt)
       {
          bin_no = bin;
          break;
       }
       prev = *curIt;
    }

    return bin_no;
}
// end of MCNP5/dagmc/Tally.cpp
