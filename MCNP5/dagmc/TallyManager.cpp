// MCNP5/dagmc/TallyManager.cpp

#include <cstdlib>
#include "TallyManager.hpp"
#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyManager::TallyManager() 
{
    event.type = TallyEvent::NONE;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
// Create and add a new Tally
void TallyManager::addNewTally(unsigned int tally_id,
                   std::string tally_type,
                   std::multimap<std::string, std::string>& options, 
                   const std::vector<double>& energy_bin_bounds)
{
	Tally *newTally = createTally(tally_id, tally_type, options,  energy_bin_bounds);
        observers.insert(std::pair<int, Tally *>(tally_id, newTally));   
}

// Remove a Tally - Observer pattern best practise
void TallyManager::removeTally(int tally_id)
{
        std::map<int, Tally *>::iterator it;	
 	it = observers.find(tally_id);
        if (it != observers.end())
        {
	   observers.erase(it);
        }
        else
        {
           std::cerr << "Warning: Tally " << tally_id << " does not exist and cannot be removed. " << std::endl;
        }
}

////////////////////////////////////////////////////////////////////
void TallyManager::end_history()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->end_history();
       }
}

////////////////////////////////////////////////////////////////////
void TallyManager::write_data(double num_histories)
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->write_data(num_histories);
       }
}

//---------------------------------------------------------------------------//
void TallyManager::clear_last_event()
{
    event.type = TallyEvent::NONE;
    event.position  = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.0, 0.0, 0.0);
    event.particle_energy     = 0.0;
    event.particle_weight     = 0.0;
    event.track_length        = 0.0;
    event.total_cross_section = 0.0;
    event.current_cell        = 0;
}
    /**
     * \brief fill the TallyEvent for a collision event
     */
    bool TallyManager::set_collision_event(double x, double y, double z,
                   double particle_energy, double particle_weight,
                   double total_cross_section, int cell_id)
    { 
       if (total_cross_section < 0.0)
       {
          std::cerr << "Warning: total_cross_section, " << total_cross_section << ", cannot be less than zero." << std::endl;
          return false;;
       }
       return set_event(TallyEvent::COLLISION, 
                 x, y, z, 0.0, 0.0, 0.0,
                 particle_energy, particle_weight, 
                 0.0, total_cross_section, 
                 cell_id); 
    } 

    /**
     * \brief fill the TallyEvent for a track event
     */
    bool TallyManager::set_track_event(double x, double y, double z,
                   double u, double v, double w,                           
                   double particle_energy, double particle_weight,
                   double track_length, int cell_id)
    { 
       if (track_length < 0.0)
       {
          std::cerr << "Warning: track_length, " << track_length << ", cannot be less than zero." << std::endl;
          return false;
       }
       return set_event(TallyEvent::TRACK, 
                 x, y, z, u, v, w,
                 particle_energy, particle_weight,
                 track_length, 0.0,
                 cell_id);
    } 

////////////////////////////////////////////////////////////////////
void TallyManager::update_tallies()
{
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;
        tally->compute_score(event);
    }
    clear_last_event();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
// Create a new Tally with the implementation that calls this method
Tally *TallyManager::createTally(unsigned int tally_id,
                   std::string  tally_type,
                   std::multimap<std::string, 
                   std::string>& options, 
                   const std::vector<double>& energy_bin_bounds)
{
        TallyInput input; 

        // Set up the input structure from the passed parameters
        input.tally_id = tally_id;
        input.options  = options;
        input.energy_bin_bounds = energy_bin_bounds;
        input.tally_type        = tally_type;
        
        return Tally::create_tally(input);
}

////////////////////////////////////////////////////////////////////
bool TallyManager::set_event(TallyEvent::EventType type,  
                           double x, double y, double z, 
                           double u, double v, double w,                           
                           double particle_energy, double particle_weight, 
                           double track_length, double total_cross_section,
                           int cell_id) 
{
    // Test whether an error condition has occurred for this event
    bool errflag = false;

    /// Set the particle state object
    event.position  = moab::CartVect(x, y, z);
    event.direction = moab::CartVect(u, v, w);
    // This should already be normalized
    event.direction.normalize();
    
    if (particle_energy < 0.0)
    {
        std::cerr << "Warning: particle_energy, " << particle_energy << ", cannot be less than zero." << std::endl;
        errflag = true; 
    }
    else
    {
        event.particle_energy = particle_energy;
    }

    if (particle_energy < 0.0)
    {
        std::cerr << "Warning: particle_weight, " << particle_weight << ", cannot be less than zero." << std::endl;
        errflag = true; 
    }
    else
    {
        event.particle_weight = particle_weight;
    }

    event.track_length        = track_length;
    event.total_cross_section = total_cross_section;
    event.current_cell        = cell_id;
    
    if (type != TallyEvent::NONE)
    {
 	event.type = type;
        errflag = true; 
    }
    else
    {
        std::cerr << "Warning: Cannot set a tally event of type NONE." << std::endl;
    }
    if (errflag)
    {
        clear_last_event();
    }    
    bool event_is_set = !errflag;
    return event_is_set;

    // Only update (that is, notify observers) if all is good
    // update_tallies();
}

// end of MCNP5/dagmc/TallyEvent.cpp
