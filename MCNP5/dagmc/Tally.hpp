// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP

#include <cassert>
#include <map>
#include <string>
#include <vector>

// Forward declare because it's only referenced here
struct TallyEvent;

//===========================================================================//
/**
 * \struct TallyInput
 * \brief Data structure for storing input data needed to set up a Tally
 *
 * Data stored in this struct is set by TallyManager and read by the different
 * Tally implementations.  Required variables for all tallies are tally_id,
 * tally_type and energy_bin_bounds.  The TallyOptions multimap is optional
 * and can be used to store all key-value pairs that define options that are
 * specific to each Tally implementation.
 */
//===========================================================================//
struct TallyInput
{
    /// Unique ID for a single Tally instance
    unsigned int tally_id;

    /// Type of Tally to create as a concrete class
    std::string tally_type;

    /// Energy bin boundaries defined for all tally points
    std::vector<double> energy_bin_bounds;

    /// Typedef for map that stores optional Tally input parameters
    typedef std::multimap<std::string, std::string> TallyOptions;

    /// Optional input parameters requested by user
    TallyOptions options;
};

//===========================================================================//
/**
 * \class Tally
 * \brief Defines an abstract Tally interface
 *
 * Tally is an abstract Base class that defines the variables and methods
 * needed to implement a generic Tally for use in Monte Carlo particle
 * codes.  All Derived classes must implement the following
 *
 *     1) compute_score(const TallyEvent& event)
 *     2) end_history()
 *     3) write_data(double num_histories)
 *
 * These three update methods are called by TallyManager (Observable) for all
 * Tally objects (Observers) that have been added.
 */
//===========================================================================//
class Tally
{
  protected:
    /**
     * \brief Constructor
     * \param input user-defined input parameters for this Tally
     */
    Tally(const TallyInput& input);

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~Tally(){}

    // >>> FACTORY METHOD

    /**
     * \brief Creates a Tally based on the given TallyInput
     * \param input user-defined input parameters for this Tally
     * \return pointer to the new Tally that was created
     *
     * NOTE: if an invalid tally_type is requested, a NULL pointer is returned
     */
    static Tally *create_tally(const TallyInput& input);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Computes scores for this Tally based on the given TallyEvent
     * \param event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event) = 0;

    /**
     * \brief Updates Tally when a particle history ends
     */
    virtual void end_history() = 0;

    /**
     * \brief Write results to an output file for this Tally
     * \param num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results to an output file for this Tally, typically normalized by
     * the number of particle histories that were tracked during the Monte
     * Carlo simulation.
     */
    virtual void write_data(double num_histories) = 0;

  protected:
    /// Input data defined by user for this tally
    TallyInput input_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

    /// Set to true by default; determines if total energy bin is included
    bool total_energy_bin;
};

#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
