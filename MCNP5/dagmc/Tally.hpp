// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP

#include <string>
#include <map>
#include <vector>
#include <cassert>

// Forward declare because it's only referenced here
struct TallyEvent;

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

    /// Optional input parameters requested by user
    TallyOptions options;

    /// Type of tally to create as a concrete class
    std::string tally_type;
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
  protected:
   /**
     * \brief Constructor
     * \param id    mesh 
     * \param input user-defined input parameters for this tally
     */
    Tally(int id, const TallyInput& input);

  public:
    
    virtual ~Tally(){}

    /**
     * \brief Defines update
     *
     * Pass Particle state data
     */
    virtual void compute_score(const TallyEvent& event) = 0;

    virtual void end_history() = 0;

    virtual void write_data(double num_histories);

    /**
     *  \brief Factory method for creation of Tally Observers
     * 
    */
    static Tally *create_tally(int id, const TallyInput& input);


  protected:

    /// Input data defined by user for this tally
    TallyInput input_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

    bool total_energy_bin;

    unsigned int tally_id;
};


#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
