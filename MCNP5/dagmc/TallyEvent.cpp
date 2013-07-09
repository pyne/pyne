// MCNP5/dagmc/TallyEvent.cpp

#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyEvent::TallyEvent(): event_type(NONE){}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
// Create a new Tally with the implementation that calls this method
// ToDo:  Create a signature with an input file instead of the multimap arg
Tally *createTally(std::multimap<std::string, std::string>& options, 
                   unsigned int tally_id,
                   const std::vector<double>& energy_bin_bounds,
                   bool total_energy_bin)
{
	Tally *ret;
        TallyInput input; 

        input.options  = options;
        input.energy_bin_bounds = energy_bin_bounds;
        input.total_energy_bin  = total_energy_bin; 
        
        ret = Tally::create_tally(tally_id, input);
        return ret;
}

// Add a Tally  
void TallyEvent::addTally(int tally_id, Tally *obs)
{
        observers.insert(std::pair<int, Tally *>(tally_id, obs));   
}

// Add a newly created Tally
void TallyEvent::addNewTally(std::multimap<std::string, std::string>& options, 
                   unsigned int tally_id,
                   const std::vector<double>& energy_bin_bounds,
                   bool total_energy_bin)

{
	Tally *newTally = createTally(options, tally_id, energy_bin_bounds, total_energy_bin);
        addTally(tally_id, newTally);
}

// Remove a Tally - Observer pattern best practise
void TallyEvent::removeTally(int tally_id)
{
        std::map<int, Tally *>::iterator it;	
 	it = observers.find(tally_id);
	observers.erase(it);
}

////////////////////////////////////////////////////////////////////
// UPDATE
void TallyEvent::update_tallies()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->update();
       }
}
/*
void TallyEvent::update_track_tallies()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->update_track();
       }
}
void TallyEvent::update_collision_tallies()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->update_collision_tallies();
       }
}
*/
////////////////////////////////////////////////////////////////////
void TallyEvent::end_history()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->end_history();
       }
}

void TallyEvent::write_data()
{
       std::map<int, Tally*>::iterator map_it;
       for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
       {
           Tally *tally = map_it->second;
	   tally->write_data();
       }
}

void TallyEvent::set_event(double x, double y, double z, 
                           double u, double v, double w,                           
                           double particle_energy, double particle_weight, 
                           double track_length, double total_cross_section) 
{
    if (track_length == 0.0 && total_cross_section == 0.0)
    {
       std::cerr << "Error:  No event type has been defined." << std::endl;
       return;
    }

    /// Set the particle state object
    particle.position            = moab::CartVect(x, y, z);
    // ToDo:  Direction is set for all event_types, but not used for collision. 
    particle.direction           = moab::CartVect(u, v, w);
    particle.energy              = particle_energy;
    particle.weight              = particle_weight;
    particle.track_length        = track_length;
    particle.total_cross_section = total_cross_section;
 
    // If more event types are needed this should become a nested if statement
    event_type = track_length > 0.0 ? TRACK : (total_cross_section > 0.0 ? COLLISION : NONE);
}

//---------------------------------------------------------------------------//
void TallyEvent::clear_last_event()
{
    event_type = NONE;
    particle.position  = moab::CartVect(0.0, 0.0, 0.0);
    particle.direction = moab::CartVect(0.0, 0.0, 0.0);
    particle.energy              = 0.0;
    particle.weight              = 0.0;
    particle.track_length        = 0.0;
    particle.total_cross_section = 0.0;
     
}
//---------------------------------------------------------------------------//
void TallyEvent::set_tally_multiplier(double value)
{
    tally_multiplier = value;
    return;
}
//---------------------------------------------------------------------------//
// TALLY EVENT ACCESS METHODS
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
double TallyEvent::get_tally_multiplier() const
{
    return tally_multiplier;
}
//---------------------------------------------------------------------------//
double TallyEvent::get_weighting_factor() const
{
    return tally_multiplier * particle.weight;
} 
//---------------------------------------------------------------------------//
TallyEvent::EventType TallyEvent::get_event_type() const
{
    return event_type;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyEvent.cpp
