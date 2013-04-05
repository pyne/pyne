// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLYEVENT_H
#define DAGMC_TALLYEVENT_H

#include "moab/CartVect.hpp"

//===========================================================================//
/**
 * \class TallyEvent
 * \brief Defines a single tally event
 *
 * TallyEvent is a class that represents a tally event that can be triggered
 * as part of a Monte Carlo particle transport simulation.  Both collision
 * and track-based events can be created.  Note that once an event has been
 * created, it cannot be changed.
 *
 * In addition to setting the required variables for each type of tally event,
 * there is also an optional function available to set the tally multiplier.
 * If the tally multiplier is not set, then it will default to 1.
 */
//===========================================================================//
class TallyEvent
{
  public:
    /**
     * \brief Constructor
     * \param energy, weight properties of particle triggering tally event
     */
    TallyEvent(double energy, double weight);

    /**
     * \brief Sets required variables for a track-based event
     * \param track_length the length of the track segment
     * \param start_point the starting location of the track (xo, yo, zo)
     * \param direction the direction the particle is traveling (uo, vo, wo)
     */
    void set_track_event(double track_length,
                         moab::CartVect& start_point,
                         moab::CartVect& direction);

    /**
     * \brief Sets required variables for a collision event
     * \param total_cross_section the total macroscopic cross section
     * \param collision_point the (x, y, z) location of the collision
     */
    void set_collision_event(double total_cross_section,
                             moab::CartVect& collision_point);

    /**
     * \brief Sets tally multiplier for the event
     * \param value energy-dependent value for the tally multiplier
     */
    void set_tally_multiplier(double value);

  private:
    /// MeshTally implementations that can access private TallyEvent data
    friend class KDEMeshTally;
    //friend class TrackLengthMeshTally;

    /**
     * \brief Defines type of tally event
     *
     *     1) NONE indicates no event has been set yet
     *     2) COLLISION indicates a collision event has been set
     *     3) TRACK indicates a track-based event has been set
     */
    enum EventType {NONE = 1, COLLISION = 2, TRACK = 3};

    /// Type of tally event this event represents
    EventType event_type;

    /// Physical quantity needed to compute the score for this event
    double event_value;

    /// Energy and weight of particle when event occurred
    double particle_energy;
    double particle_weight;

    /// Energy-dependent multiplier for this event
    double tally_multiplier;

    /// Position and direction of particle when event occurred
    moab::CartVect position;
    moab::CartVect direction;
};

#endif // DAGMC_TALLYEVENT_H

// end of MCNP5/dagmc/TallyEvent.hpp
