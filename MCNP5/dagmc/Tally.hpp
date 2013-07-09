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
    virtual void update(const ParticleState& state) = 0;
/*
    virtual void update_track(const ParticleState& state) = 0;
    virtual void update_collision(const ParticleState& state) = 0;
*/ 
    virtual void end_history() = 0;

    virtual void write_data(double num_particles, double multiplier = 1.0) = 0;

    /**
     *  \brief Factory method for creation of Tally Observers
     * 
     *
    */
    static Tally *create_tally(int id, const TallyInput& input);

    unsigned int tally_id;

  protected:

    Tally(int id, const TallyInput& input);

    /// Input data defined by user for this tally
    TallyInput input_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

};


#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
