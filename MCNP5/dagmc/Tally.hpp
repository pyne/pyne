// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP

#include <string>
#include <map>
#include <vector>
#include <cassert>

//===========================================================================//
/**
 * \struct TallyInput
 * \brief Input data needed to set up a tally (could be mesh or other)
 *  May add time bins, angle bins etc.
 */
//===========================================================================//
struct TallyInput
{
    /// Typedef for map that stores optional tally input parameters
    typedef std::multimap<std::string, std::string> TallyOptions;

    /// Energy bin boundaries defined for all tally points
    std::vector<double> energy_bin_bounds;

    /// If true, add an extra energy bin to tally all energy levels
    bool total_energy_bin;

    /// Optional input parameters requested by user
    TallyOptions options;
};


//===========================================================================//
/**
 * \class Tally
 * \brief Defines a single tally event
 *
 * Tally is a pure virtual class that is the Observer class of the Observer pattern.
 * It is part of a Monte Carlo particle transport simulation.  
 *
 */
//===========================================================================//
class Tally
{
  public:
    
    virtual ~Tally(){}

    /**
     * \brief Defines update
     *
     * Pass Particle state data
     */
    virtual void update() = 0;
 
    virtual void end_history() = 0;

    virtual void write_data() = 0;

    /**
     *  \brief Factory method for creation of Tally Observers
     * 
     *
    */
    static Tally *create_tally(int id, const TallyInput& input);

    unsigned int tally_id;

  protected:

    /// Input data defined by user for this tally
    TallyInput input_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

    Tally(int id, const TallyInput input)
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
};


#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
