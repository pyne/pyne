// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP

#include "moab/CartVect.hpp"

//===========================================================================//
/**
 * \struct TallyEvent
 * \brief Data structure for storing event data to be tallied
 *
 * Data stored in this struct is set by TallyManager and read by the different
 * Tally implementations.  The variables that will be used by each Tally depend
 * on the type of event that is defined.  All events set type, particle_energy
 * and particle_weight. Collision events add the current_cell, position
 * (i.e. collision point) and total_cross_section.  Track events add the
 * current_cell, position (i.e. start of track), direction and track_length.
 */
//===========================================================================//
struct TallyEvent
{
    /**
     * \brief Constructor
     */
    TallyEvent();

    /**
     * \brief Defines type of tally event
     *
     *     0) NONE indicates no event has been set yet
     *     1) COLLISION indicates a collision event has been set
     *     2) TRACK indicates a track-based event has been set
     */
    enum EventType {NONE = 0, COLLISION = 1, TRACK = 2};

    EventType type;
    
    /// Geometric cell in which the event occurred
    int current_cell;

    /// Total length of track segment
    double track_length;

    /// Position of particle (x, y, z)
    moab::CartVect position;

    /// Direction in which particle is traveling (u, v, w)
    moab::CartVect direction;

    /// Total macroscopic cross section for cell in which collision occurred
    double total_cross_section;

    /// Energy and weight of particle when event occurred
    double particle_energy;
    double particle_weight;

    // >>> PUBLIC INTERFACE

    /**
     * \brief Defines this tally event as a track-based event
     */
    void set_track_event(double track_length,
                         double x, double y, double z,
                         double u, double v, double w,
                         double particle_energy,
                         double particle_weight);

    /**
     * \brief Defines this tally event as a collision event
     */
    void set_collision_event(double total_cross_section,
                             double x, double y, double z,
                             double particle_energy,
                             double particle_weight);

    /**
     * \brief Sets tally multiplier for the event
     * \param value energy-dependent value for the tally multiplier
     */
    void set_tally_multiplier(double value);

    /**
     * \brief get_tally_multiplier()
     * \return current value of the tally multiplier for this event
     */
    double get_tally_multiplier() const;

    /**
     * \brief Gets the total weighting factor for this event
     * \return product of tally multiplier and particle weight
     */
    double get_weighting_factor() const;

  private:

    /// Energy-dependent multiplier for this event
    double tally_multiplier;
};

#endif // DAGMC_TALLY_EVENT_HPP

// end of MCNP5/dagmc/TallyEvent.hpp
