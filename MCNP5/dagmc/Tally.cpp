// MCNP5/dagmc/Tally.cpp

#include <cassert>
#include <iostream>
#include <cmath>

#include "Tally.hpp"
#include "TrackLengthMeshTally.hpp"
#include "KDEMeshTally.hpp"
#include "CellTally.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Tally::Tally(const TallyInput& input) 
    : input_data(input), data(NULL)
{
    assert(input_data.energy_bin_bounds.size() > 1);

    // This is a placeholder for a future option to set t.e.b. false via the
    // TallyInput 
    bool total_energy_bin = true;

    unsigned int num_energy_bins = input_data.energy_bin_bounds.size() - 1;

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
//               | Track Length   | Cell         || cell_track    implemented 
//               | Collision      | Cell         || cell_coll     implemented 
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
    else if (input.tally_type == "cell_track")
    {
        newTally = new CellTally(input, TallyEvent::TRACK);
    }
    else if (input.tally_type == "cell_coll")
    {
        newTally = new CellTally(input, TallyEvent::COLLISION);
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
const TallyData& Tally::getTallyData()
{
    return *data;
}
//---------------------------------------------------------------------------//
std::string Tally::get_tally_type()
{
    return input_data.tally_type;
}
//---------------------------------------------------------------------------//
// PROTECTED INTERFACE
//---------------------------------------------------------------------------//
bool Tally::get_energy_bin(double energy, unsigned int& ebin)
{
    bool bin_exists = false;

    if (energy_in_bounds(energy))
    {
        // in bounds, energy bin index must exist
        bin_exists = true;

        if (data->get_num_energy_bins() == 1)
        {
            ebin = 0;
        }
        else  // in bounds and more than one energy bin
        {
            unsigned int max_ebound = input_data.energy_bin_bounds.size() - 1;

            // Pre-load ebin with maximum bin as default
            ebin =  max_ebound - 1;

            // find ebin if not maximum bin
	    for (unsigned int i=0; i < max_ebound; ++i)
            {
                if (input_data.energy_bin_bounds.at(i) <= energy &&
                    energy < input_data.energy_bin_bounds.at(i+1))
                {
                    ebin = i;
		    break;
                }
            }  // end for
        }  // end else in bounds and >1 energy bin
    }  // end if in bounds

    return bin_exists;
}
//---------------------------------------------------------------------------//
bool Tally::energy_in_bounds(double energy)
{
    unsigned int max_ebound = input_data.energy_bin_bounds.size() - 1;

    return !(energy < input_data.energy_bin_bounds.at(0) ||
             energy > input_data.energy_bin_bounds.at(max_ebound));
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/Tally.cpp
