// MCNP5/dagmc/TallyManager.hpp

#ifndef DAGMC_TALLY_MANAGER_HPP
#define DAGMC_TALLY_MANAGER_HPP

#include <iostream>
#include "Tally.hpp"
#include "TallyEvent.hpp"

//===========================================================================//
/**
 * \class TallyManager
 * \brief Defines a single tally event
 *
 * ToDo:  Update this doc
 * TallyManager is a class that represents a tally event that can be triggered
 * as part of a Monte Carlo particle transport simulation.  Both collision
 * and track-based events can be created.  Note that once an event type has
 * been set it cannot be changed.
 *
 * This class is a Subject (Observable) class that keeps a list of Observers
 * to be Notified.  The Observers are added via the method attach(..)
 */
//===========================================================================//
class TallyManager
{
  public:
    /**
     * \brief Constructor
     * Does not compile without this, due to event_type setting
     * in implementation, maybe
     */
    TallyManager();
    

    // >>> PUBLIC INTERFACE
    // Create and a new Tally
    // ToDo: turn into doc
    // ToDo:  Do we need this to be multimap?
    /// User-specified ID for this tally
    /// Type of tally to create
    /// Options  map that stores optional tally input parameters
    /// Energy bin boundaries defined for all tally points
    /// If true, add an extra energy bin to tally all energy levels
    Tally *createTally(unsigned int tally_id,
                   std::string  tally_type,
                   std::multimap<std::string, std::string>& options, 
                   const std::vector<double>& energy_bin_bounds);

    // Add a Tally to the observer list
    void addTally(int tally_id, Tally *obs);


    // Create a Tally and add it to the observer list    
    void addNewTally(unsigned int tally_id,
                   std::string tally_type,
                   std::multimap<std::string, std::string>& options, 
                   const std::vector<double>& energy_bin_bounds);

    // Remove a Tally
    void removeTally(int tally_id);

    // Call action on every tally
    void update_tallies();

    // Call end_history() on every tally
    void end_history();

    // Call write_data() on every tally
    void write_data(double num_histories);

    /**
     * \brief fill the ParticleState
     * \param position interpreted differently for TrackLength and Collision 
     * Requires one of track_length and total_cross_section to be nonzero
     * ToDo:  may not need u, v, w in signature 
     */
    void set_event(double x, double y, double z,
                   double u, double v, double w,                           
                   double particle_energy, double particle_weight,
                   double track_length = 0.0, double total_cross_section = 0.0); 
                           

    /**
     *  \brief Reset TallyEvent data
     *
     *  Set eventType to NONE and clear the particle data
     */
    void clear_last_event();

  private:

    // Keep a record of the Observers
    std::map <int, Tally*> observers; 

    // Store particle state for the event
    TallyEvent event;
};

#endif // DAGMC_TALLY_MANAGER_HPP

// end of MCNP5/dagmc/TallyManager.hpp
