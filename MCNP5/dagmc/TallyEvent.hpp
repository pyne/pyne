// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP

// ToDo:  Clean up;  This file is being modified for Issue #90
#include <map>
#include "Tally.hpp"

#include "moab/CartVect.hpp"

struct ParticleState
{
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
    double energy;
    double weight;
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
 * 
 * This class is a Subject (Observable) class that keeps a list of Observers
 * to be Notified.  The Observers are added via the method attach(..)
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
     * Does not compile without this, due to event_type setting
     * in implementation, maybe
     */
    TallyEvent();
    
    // Keep a record of the Observers
    std::map <int, Tally *> observers; 

    // >>> PUBLIC INTERFACE

    // Add a new Tally
    void addTally(int, Tally *);

    // Remove a Tally
    void removeTally(int map_index, Tally *obs);

    // Call action on every tally
    void update_tallies();

    /**
     * \brief Sets tally multiplier for the event
     * \param value energy-dependent value for the tally multiplier
     */
    void set_tally_multiplier(double value);

    // >>> TALLY EVENT DATA ACCESS METHODS

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

    /**
     * \brief fill the ParticleState
     * \param position interpreted differently for TrackLength and Collision 
     */
    void set_event(const moab::CartVect& position, 
                           const moab::CartVect& direction,
                           double track_length, double total_cross_section, 
                           double particle_energy, double particle_weight, 
                           EventType type);

  private:

    /// Energy-dependent multiplier for this event
    double tally_multiplier;

    /// Type of tally event this event represents
    EventType event_type;

    // Store particle state for the event
    ParticleState particle;
};


#endif // DAGMC_TALLYEVENT_H

// end of MCNP5/dagmc/TallyEvent.hpp
