// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP

#include <string>
#include <map>
#include <vector>

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

    /// User-specified ID for this tally
    unsigned int tally_id;

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
    static Tally *create_tally(const TallyInput& input);
};

#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
