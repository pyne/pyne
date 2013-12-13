// MCNP5/dagmc/CellTally.cpp

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <utility>

#include "CellTally.hpp" 

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
CellTally::CellTally(const TallyInput& input, TallyEvent::EventType eventType)
    : Tally(input),
      cell_id(1),
      cell_volume(1.0),
      expected_type(eventType)
{
    // Set up CellTally member variables from TallyInput
    parse_tally_options();

    // Initialize the data arrays to store one tally point
    data->resize_data_arrays(1);
}
//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from Tally.hpp
//---------------------------------------------------------------------------//
void CellTally::compute_score(const TallyEvent& event)
{
    // Return if current cell or particle energy is incompatible with CellTally
    unsigned int ebin = 0;

    if (event.current_cell != cell_id ||
        !get_energy_bin(event.particle_energy, ebin))
    {
        return;
    }

    // Compute score based on event type and add it to this CellTally
    double event_score = event.get_score_multiplier(input_data.multiplier_id);

    if (event.type == TallyEvent::TRACK && event.type == expected_type)
    {
        event_score *= event.track_length;         
    }
    else if (event.type == TallyEvent::COLLISION && event.type == expected_type)
    {
        event_score /= event.total_cross_section;
    }
    else // NONE, return from this method
    {
        return;
    }

    int tally_index = 0;
    data->add_score_to_tally(tally_index, event_score, ebin);
}
//---------------------------------------------------------------------------//
void CellTally::write_data(double num_histories)
{
    std::cout << "Writing data for CellTally " << input_data.tally_id
              << ": " << std::endl;

    std::cout << "cell id = " << cell_id << std::endl;
    std::cout << "type = " <<          
	     (expected_type == TallyEvent::COLLISION ? "collision " :     
              expected_type == TallyEvent::TRACK     ? "track "     : "none");
    std::cout << std::endl;

    std::cout << "volume = " << cell_volume << std::endl << std::endl; 

    // Get data for CellTally and print final results to std::cout
    unsigned int point_index = 0;
    unsigned int num_bins = data->get_num_energy_bins();

    for (unsigned int i = 0; i < num_bins; ++i)
    {
        if (data->has_total_energy_bin() && (i == num_bins - 1))
        {
            std::cout << "Total Energy Bin: " << std::endl;
        }
        else
        {
            std::cout << "Energy bin (" << input_data.energy_bin_bounds.at(i)
                      << ", " << input_data.energy_bin_bounds.at(i+1) << "):\n";
        }

        std::pair <double, double> tally_data = data->get_data(point_index, i);
        double tally = tally_data.first;
        double error = tally_data.second;

        // compute relative error for the tally result
        double rel_error = 0.0;

        if (error != 0.0)
        {
            rel_error = sqrt(error / (tally * tally) - 1.0 / num_histories);
        }

        // normalize mesh tally result by the number of source particles
        tally /= (num_histories * cell_volume);

        std::cout << "    tally = " << tally << std::endl;
        std::cout << "    error = " << rel_error << std::endl;
        std::cout << std::endl;
    }
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
void CellTally::parse_tally_options()
{
    const TallyInput::TallyOptions& options = input_data.options;  
    TallyInput::TallyOptions::const_iterator it;

    for (it = options.begin(); it != options.end(); ++it)
    {
        std::string key = it->first;
        std::string value = it->second;

        // process tally option according to key
        if (key == "cell")
        {
            char* end;
            cell_id = strtol(value.c_str(), &end, 10);

            if (value.c_str() == end)
            {
                std::cerr << "Warning: '" << value << "' is an invalid value"
                          << " for the cell id" << std::endl;
                cell_id = 1;
                std::cerr << "    setting cell id to " << cell_id
                      << " for a CellTally option." << std::endl;   
            }
        }
        else if (key == "volume")
        {
            char* end;
            cell_volume = strtod(value.c_str(), &end);

            if (value.c_str() == end)
            {
               std::cerr << "Warning: '" << value << "' is an invalid value"
                         << " for the cell volume" << std::endl;
            }
        }
        else // invalid tally option
        {
            std::cerr << "Warning: input data for cell tally "
                      << input_data.tally_id
                      << " has unknown key '" << key << "'" << std::endl;
        }
    }
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/CellTally.cpp
