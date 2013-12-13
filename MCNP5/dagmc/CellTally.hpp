// MCNP5/dagmc/CellTally.hpp

#ifndef DAGMC_CELL_TALLY_HPP
#define DAGMC_CELL_TALLY_HPP

#include "Tally.hpp"
#include "TallyEvent.hpp"

//===========================================================================//
/**
 * \class CellTally
 * \brief Defines a simple CellTally interface
 *
 * CellTally is a concrete class derived from Tally that provides a simple
 * cell-based tally implementation.  Two types of cell tallies can be created
 *
 *     1) TallyEvent::COLLISION cell tally, based on particle collisions
 *     2) TallyEvent::TRACK cell tally, based on full particle tracks
 *
 * This class was designed to be used primarily for unit testing purposes, but
 * can also be used as a fully-functioning cell tally within a Monte Carlo
 * particle transport simulation.  Both energy bins and tally multipliers are
 * supported.
 *
 * ==========
 * TallyInput
 * ==========
 *
 * The TallyInput struct needed to construct a CellTally object is defined
 * in Tally.hpp and is set through the TallyManager when a Tally is created.
 * Options that are currently available for CellTally objects include
 *
 * 1) "cell"="value"
 * -----------------
 * Sets the cell ID to the given value, which should represent the index of an
 * actual geometric cell.  Currently only one cell can be tallied per CellTally,
 * and the default value is 1.
 *
 * 2) "volume"="value"
 * -------------------
 * Sets the volume for the cell ID that will be used to normalize the final
 * tally result.  Note that this quantity is not computed by the CellTally,
 * and the default value is 1.0.
 */
//===========================================================================//
class CellTally : public Tally
{
  public:
    /**
     * \brief Constructor
     * \param[in] input user-defined input parameters for this CellTally
     * \param[in] eventType the type of event that is to be tallied
     */
    CellTally(const TallyInput& input, TallyEvent::EventType eventType);

    /**
     * \brief Destructor
     */
    virtual ~CellTally(){}

    // >>> PUBLIC INTERFACE

    /**
     * \brief Computes scores for this CellTally based on the given TallyEvent
     * \param[in] event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event);

    /**
     * \brief Write results for this CellTally
     * \param[in] num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results to std::out for this CellTally, normalized by both the
     * number of particle histories that were tracked and the volume of the
     * cell for which the results were computed.
     */
    virtual void write_data(double num_histories);

  private:
    // ID for the geometric cell to be tallied
    int cell_id;

    // Volume for the geometric cell to be tallied
    double cell_volume;

    // Event type used by this CellTally(COLLISION or TRACK, not both)
    TallyEvent::EventType expected_type;

    /**
     * \brief Parse the TallyInput options for this CellTally
     */
    void parse_tally_options();
};

#endif // DAGMC_CELL_TALLY_HPP

// end of MCNP5/dagmc/CellTally.hpp
