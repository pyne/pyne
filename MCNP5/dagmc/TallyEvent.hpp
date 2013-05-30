// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP

#include <utility>

#include "moab/CartVect.hpp"

//===========================================================================//
/**
 * \struct TrackData
 * \brief Stores physics data needed to construct a track-based tally event
 */
//===========================================================================//
struct TrackData
{
    /// Total length of track segment
    double track_length;

    /// Starting location of track segment (xo, yo, zo)
    moab::CartVect start_point;

    /// Direction in which particle is traveling (uo, vo, wo)
    moab::CartVect direction;
};

//===========================================================================//
/**
 * \struct CollisionData
 * \brief Stores physics data needed to construct a collision tally event
 */
//===========================================================================//
struct CollisionData
{
    /// Total macroscopic cross section for cell in which collision occurred
    double total_cross_section;

    /// Location of the collision (x, y, z)
    moab::CartVect collision_point;
};

//===========================================================================//
/**
 * \class TallyEvent
 * \brief Defines a single tally event
 *
 * TallyEvent is a class that represents a tally event that can be triggered
 * as part of a Monte Carlo particle transport simulation.  Both collision
 * and track-based events can be created.  Note that once an event type has
 * been set it cannot be changed.
 *
 * In addition to setting the required variables for each type of tally event,
 * there is also an optional method available to set the tally multiplier.
 * If the tally multiplier is not set, then it will default to 1.
 */
//===========================================================================//
class TallyEvent
{
  public:
    /**
     * \brief Defines type of tally event
     *
     *     1) NONE indicates no event has been set yet
     *     2) COLLISION indicates a collision event has been set
     *     3) TRACK indicates a track-based event has been set
     */
    enum EventType {NONE = 1, COLLISION = 2, TRACK = 3};

    /**
     * \brief Constructor
     * \param energy, weight properties of particle triggering event
     */
    TallyEvent(double energy, double weight);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Defines this tally event as a track-based event
     * \param data the physics data needed to set a track-based event
     */
    void set_track_event(const TrackData& data);

    /**
     * \brief Defines this tally event as a collision event
     * \brief data the physics data needed to set a collision event
     */
    void set_collision_event(const CollisionData& data);

    /**
     * \brief Sets tally multiplier for the event
     * \param value energy-dependent value for the tally multiplier
     */
    void set_tally_multiplier(double value);

    // >>> TALLY EVENT DATA ACCESS METHODS

    /**
     * \brief Gets particle energy and weight for this tally event
     * \return the particle data in form of pair<energy, weight>
     */
    std::pair<double, double> get_particle_data() const;

    /**
     * \brief get_track_data(), get_collision_data()
     * \param data the struct to which the tally event data will be copied
     * \return true if data was copied successfully, false otherwise
     */
    bool get_track_data(TrackData& data) const;
    bool get_collision_data(CollisionData& data) const;

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

    /**
     * \brief get_event_type()
     * \return type of tally event this event represents
     */
    TallyEvent::EventType get_event_type() const;

  private:

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

    /// Type of tally event this event represents
    EventType event_type;
};

#endif // DAGMC_TALLY_EVENT_HPP

// end of MCNP5/dagmc/TallyEvent.hpp
