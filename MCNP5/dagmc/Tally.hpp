// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP

#include <map>
#include <string>
#include <vector>

#include "TallyData.hpp"

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
 * specific to each Tally implementation.  The multiplier_id is also optional
 * and is used to indicate what energy-dependent multiplier is to be used
 * with the Tally.
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
 
    /// Support a single multiplier for each tally; this id refers to an
    /// index in the TallyEvent::multipliers vector
    int multiplier_id;
};

//===========================================================================//
/**
 * \class Tally
 * \brief Defines an abstract Tally interface
 *
 * Tally is an abstract Base class that defines the variables and methods
 * needed to implement a generic Tally for use in Monte Carlo particle
 * codes.  All Derived classes must implement the following methods
 *
 *     1) compute_score(const TallyEvent& event)
 *     2) end_history()
 *     3) write_data(double num_histories)
 *
 * These three update methods are called by TallyManager (Observable) for all
 * Tally objects (Observers) that have been added.
 *
 * Note that Tally provides a default end_history() method that should be
 * sufficient for most Tally objects that use the TallyData structure for
 * storing their data.  If a different data structure is used, or alternative
 * behavior is desired, then Derived classes can override this method.
 */
//===========================================================================//
class Tally
{
  protected:
    /**
     * \brief Constructor
     * \param[in] input user-defined input parameters for this Tally
     */
    explicit Tally(const TallyInput& input);

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~Tally();

    // >>> FACTORY METHOD

    /**
     * \brief Creates a Tally based on the given TallyInput
     * \param[in] input user-defined input parameters for this Tally
     * \return pointer to the new Tally that was created
     *
     * NOTE: if an invalid tally_type is requested, a NULL pointer is returned
     */
    static Tally *create_tally(const TallyInput& input);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Computes scores for this Tally based on the given TallyEvent
     * \param[in] event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event) = 0;

    /**
     * \brief Updates Tally when a particle history ends
     */
    virtual void end_history();

    /**
     * \brief Write results for this Tally
     * \param[in] num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results to std::out or an output file for this Tally, typically
     * normalized by the number of particle histories that were tracked during
     * the Monte Carlo simulation.
     */
    virtual void write_data(double num_histories) = 0;

    /**
     * \brief Provide access to data for testing
     *
     * Provide read-only access to the protected TallyData *data
     */
    const TallyData& getTallyData();

  protected:
    /// Input data defined by user for this tally
    TallyInput input_data;

    /// All of the tally data for this tally
    TallyData *data;

    /**
     * \brief Get the bin index for the current energy
     * \param[in] energy the current particle energy
     * \param[out] ebin the energy bin index corresponding to the energy
     * \return true if energy bin is found; false otherwise
     *
     * Calls energy_in_bounds first to check if energy is valid for this Tally.
     */
    bool get_energy_bin(double energy, unsigned int& ebin);

    /// The purpose of this is to allow TallyManager to use the data
    friend class TallyManager;

  private:
    /**
     * \brief Check that the energy is not outside the allowed range
     * \param[in] energy the current particle energy
     * \return true if energy is within the allowed range; false otherwise
     */
    bool energy_in_bounds(double energy);
};

#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
