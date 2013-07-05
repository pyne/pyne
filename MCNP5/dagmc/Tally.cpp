// MCNP5/dagmc/Tally.cpp
#include "Tally.hpp"

    Tally::Tally(int id, const TallyInput& input)
    {    
       tally_id   = id; 
       input_data = input;

       // Determine the total number of energy bins requested
       num_energy_bins = input_data.energy_bin_bounds.size();

       if(!input_data.total_energy_bin)
       {
          --num_energy_bins;
       }

       assert(num_energy_bins > 0);
    }

    /**
     *  \brief Factory method for creation of Tally Observers
     * 
     *
    */
    Tally *Tally::create_tally(int id, const TallyInput& input)
    {
       return NULL;
    }

// end of MCNP5/dagmc/Tally.cpp
