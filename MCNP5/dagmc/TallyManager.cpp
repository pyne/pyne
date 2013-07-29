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
// Create a new Tally with the implementation that calls this method
Tally *createTally(std::string  tally_type,
                   std::multimap<std::string, std::string>& options, 
                   const std::vector<double>& energy_bin_bounds)
{
	Tally *ret;
        TallyInput input; 

        // Set up the input structure from the passed parameters
        input.options  = options;
        input.energy_bin_bounds = energy_bin_bounds;
        input.tally_type        = tally_type;
        
        ret = Tally::create_tally(input);
        return ret;
}

// Add a Tally  
void TallyManager::addTally(int tally_id, Tally *obs)
{
        observers.insert(std::pair<int, Tally *>(tally_id, obs));   
}

// Add a newly created Tally
void TallyManager::addNewTally(unsigned int tally_id,
                   std::string tally_type,
                   std::multimap<std::string, std::string>& options, 
                   const std::vector<double>& energy_bin_bounds)
{
	Tally *newTally = createTally(tally_id, tally_type, options,  energy_bin_bounds);
        addTally(tally_id, newTally);
}

// Remove a Tally - Observer pattern best practise
void TallyManager::removeTally(int tally_id)
{
        std::map<int, Tally *>::iterator it;	
 	it = observers.find(tally_id);
	observers.erase(it);
}

////////////////////////////////////////////////////////////////////
// UPDATE
void TallyManager::update_tallies()
{
    if (event.type == TallyEvent::NONE)
    {
        std::cerr << "Error:  No event type has been defined." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;
        tally->compute_score(event);
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

void TallyManager::write_data(double num_histories)
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->write_data(num_histories);
       }
}

void TallyManager::set_event(double x, double y, double z, 
                           double u, double v, double w,                           
                           double particle_energy, double particle_weight, 
                           double track_length, double total_cross_section,
                           int cell_id) 
{
    /// Set the particle state object
    event.position            = moab::CartVect(x, y, z);
    event.direction           = moab::CartVect(u, v, w);
    event.particle_energy     = particle_energy;
    event.particle_weight     = particle_weight;
    event.track_length        = track_length;
    event.total_cross_section = total_cross_section;
    event.current_cell        = cell_id;
 
    // If more event types are needed this should become a nested if statement
    // event.type = track_length > 0.0 ? TallyEvent::TRACK : (total_cross_section > 0.0 ? 
    //                                   TallyEvent::COLLISION : 
    //                                   TallyEvent::NONE);
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
     
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyEvent.cpp
