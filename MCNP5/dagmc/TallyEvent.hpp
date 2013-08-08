// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP
#include "moab/CartVect.hpp"
 

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
 *
 * TODO: This is a temporary class structure that will be changed once we
 * merge with the observer pattern implementation in meshtally_reftest branch.
 * Taking this first step will make it easier to merge with KDEMeshTally since
 * I am making some substantial changes to it as a part of the boundary
 * correction method.
 */
//===========================================================================//
class TallyEvent
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

    /// Total length of track segment
    double track_length;

    /// Position of particle (x, y, z)
    moab::CartVect position;

    /// Direction in which particle is traveling (u, v, w)
    moab::CartVect direction;

    /// Collision: Total macroscopic cross section for cell in
    /// which collision occurred
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
