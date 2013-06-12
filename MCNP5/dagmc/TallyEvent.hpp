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
    
    // Keep a record of the Observers
    // 3.  Coupled only to the base Observer class
    std::map <int, Tally *> observers; 

    // >>> PUBLIC INTERFACE

    // Add a new Tally
    void addTally(int map_index, Tally *obs)
    {
        observers.insert(map_index, obs);   
    }

    // Remove a Tally
    void removeTally(int map_index, Tally *obs)
    {
        std::map<int, Tally *>::iterator it;	
 	it = observers.find(map_index);
	observers.erase(it);
    }

    void update_tallies()
    {
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->update();
       }
    }

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

  private:

    /// Energy-dependent multiplier for this event
    double tally_multiplier;

    /// Type of tally event this event represents
    EventType event_type;

#endif // DAGMC_TALLYEVENT_H

// end of MCNP5/dagmc/TallyEvent.hpp
