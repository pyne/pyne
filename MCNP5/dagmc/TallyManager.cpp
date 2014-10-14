// MCNP5/dagmc/TallyManager.cpp

#include <cstdlib>
#include <iostream>

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
void TallyManager::addNewTally(unsigned int tally_id,
                               std::string tally_type,
                               unsigned int particle,
                               const std::vector<double>& energy_bin_bounds,
                               const std::multimap<std::string, std::string>& options)
{
    Tally *newTally = createTally(tally_id, tally_type, particle,
                                  energy_bin_bounds, options);

    if (newTally != NULL)
    {
        observers.insert(std::pair<int, Tally*>(tally_id, newTally));   
    }
    else
    {
        std::cerr << "Warning: Tally will be ignored." << std::endl;
    }
}
//---------------------------------------------------------------------------//
void TallyManager::addNewMultiplier(unsigned int multiplier_id)
{
    // pad multipliers vector up to a size one greater than the multiplier_id
    // NOTE: this would not be needed if we use an unordered map over a vector
    while (event.multipliers.size() <= multiplier_id)
    {
        event.multipliers.push_back(1.0);
    }
}
//---------------------------------------------------------------------------//
void TallyManager::addMultiplierToTally(unsigned int multiplier_id,
                                        unsigned int tally_id)
{
    std::map<int, Tally *>::iterator it;	
    it = observers.find(tally_id);

    if (event.multipliers.size() > multiplier_id && it != observers.end())
    {
        Tally *tally = it->second;  
        tally->input_data.multiplier_id = multiplier_id;
    }
    else
    {
        std::cerr << "Warning: Cannot set multiplier id for Tally " << tally_id
                  << ".  Tally and/or multiplier are/is invalid." << std::endl;
    }
}
//---------------------------------------------------------------------------//
void TallyManager::updateMultiplier(unsigned int multiplier_id, double value)
{
    if (event.multipliers.size() > multiplier_id)
    {
        event.multipliers.at(multiplier_id) = value; 
    }
}
//---------------------------------------------------------------------------//
unsigned int TallyManager::numTallies()
{
    return observers.size();
}
//---------------------------------------------------------------------------//
void TallyManager::removeTally(unsigned int tally_id)
{
    std::map<int, Tally *>::iterator it;	
    it = observers.find(tally_id);

    if (it != observers.end())
    {
        // release memory allocated to Tally and remove it from the map
        delete it->second;
        observers.erase(it);
    }
    else
    {
        std::cerr << "Warning: Tally " << tally_id
                  << " does not exist and cannot be removed. " << std::endl;
    }
}
//---------------------------------------------------------------------------//
bool TallyManager::setCollisionEvent(unsigned int particle, 
                                     double x, double y, double z,
                                     double particle_energy, double particle_weight,
                                     double total_cross_section, int cell_id)
{ 
    if (total_cross_section < 0.0)
    {
        std::cerr << "Warning: total_cross_section, " << total_cross_section
                  << ", cannot be less than zero." << std::endl;
        return false;;
    }

    return setEvent(TallyEvent::COLLISION, particle,
                     x, y, z, 0.0, 0.0, 0.0,
                     particle_energy, particle_weight, 
                     0.0, total_cross_section, 
                     cell_id);
} 
//---------------------------------------------------------------------------//
bool TallyManager::setTrackEvent(unsigned int particle,
                                 double x, double y, double z,
                                 double u, double v, double w,                           
                                 double particle_energy, double particle_weight,
                                 double track_length, int cell_id)
{ 
    if (track_length < 0.0)
    {
        std::cerr << "Warning: track_length, " << track_length
                  << ", cannot be less than zero." << std::endl;
        return false;
    }

    return setEvent(TallyEvent::TRACK, particle,
                     x, y, z, u, v, w,
                     particle_energy, particle_weight,
                     track_length, 0.0,
                     cell_id);
}
//---------------------------------------------------------------------------//
void TallyManager::clearLastEvent()
{
    event.type = TallyEvent::NONE;
    event.particle  = 0;
    event.position  = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.0, 0.0, 0.0);
    event.particle_energy     = 0.0;
    event.particle_weight     = 0.0;
    event.track_length        = 0.0;
    event.total_cross_section = 0.0;
    event.current_cell        = 0;
}
//---------------------------------------------------------------------------//
// Note: the event is set just before updateTallies is called
void TallyManager::updateTallies()
{
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;

        // skip events involving particles not expected by the tally
        if (tally->input_data.particle == event.particle)
        { 
           tally->compute_score(event);
        }
    }
    clearLastEvent();
}
//---------------------------------------------------------------------------//
void TallyManager::endHistory()
{
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;
	tally->end_history();
    }
}
//---------------------------------------------------------------------------//
void TallyManager::writeData(double num_histories)
{
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;
	tally->write_data(num_histories);
    }
}
//---------------------------------------------------------------------------//
// TALLY DATA ACCESS METHODS
//---------------------------------------------------------------------------//
// TODO: These will only work if TallyData is used to store all data.
// Future addition could add similar functions to the Tally interface so that
// each implementation can choose how to store its data.
double* TallyManager::getTallyData(int tally_id, int& length)
{
    std::map<int, Tally *>::iterator it;	
    it = observers.find(tally_id);

    if (it != observers.end())
    {
        Tally *tally = it->second;
        return tally->data->get_tally_data(length);;
    }
    else
    {
        std::cerr << "Warning: Tally " << tally_id
                  << " does not exist and cannot be accessed for tally data. " << std::endl;
        return NULL;
    }
}
//---------------------------------------------------------------------------//
double* TallyManager::getErrorData(int tally_id, int& length)
{
    std::map<int, Tally *>::iterator it;	
    it = observers.find(tally_id);

    if (it != observers.end())
    {
        Tally *tally = it->second;
        return tally->data->get_error_data(length);;
    }
    else
    {
        std::cerr << "Warning: Tally " << tally_id
                  << " does not exist and cannot be accessed for error data. " << std::endl;
        return NULL;
    }
}
//---------------------------------------------------------------------------//
double* TallyManager::getScratchData(int tally_id, int& length)
{
    std::map<int, Tally *>::iterator it;	
    it = observers.find(tally_id);

    if (it != observers.end())
    {
        Tally *tally = it->second;
        return tally->data->get_scratch_data(length);;
    }
    else
    {
        std::cerr << "Warning: Tally " << tally_id
                  << " does not exist and cannot be accessed for scratch data. " << std::endl;
        return NULL;
    }
}
//---------------------------------------------------------------------------//
void TallyManager::zeroAllTallyData()
{
    std::map<int, Tally*>::iterator map_it;
    for (map_it = observers.begin(); map_it != observers.end(); ++map_it)
    {
        Tally *tally = map_it->second;
        tally->data->zero_tally_data();
    }
    clearLastEvent();
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
Tally *TallyManager::createTally(unsigned int tally_id,
                                 std::string  tally_type,
                                 unsigned int particle,
                                 const std::vector<double>& energy_bin_bounds,
                                 const std::multimap<std::string, std::string>& options)
{
    TallyInput input; 

    // Check that there are at least 2 energy bin boundaries
    if (energy_bin_bounds.size() < 2)
    {
        std::cerr << "Warning: energy bin boundaries for Tally " << tally_id
                  << " are invalid." << std::endl;
        return NULL;
    }

    switch (particle)
    {
      case 1:
        input.particle = TallyInput::NEUTRON;
        break;
 
      case 2:
        input.particle = TallyInput::PHOTON;
        break;

      case 3:
        input.particle = TallyInput::ELECTRON;
        break;

      default:
        input.particle = TallyInput::NEUTRON;
    }

    // Set up the input structure from the passed parameters
    input.tally_id          = tally_id;
    input.tally_type        = tally_type;
    input.energy_bin_bounds = energy_bin_bounds;
    input.options           = options;
    input.multiplier_id     = -1;      // Turn off multipliers by default
   
    return Tally::create_tally(input);
}
//---------------------------------------------------------------------------//
bool TallyManager::setEvent(TallyEvent::EventType type, unsigned int particle,
                           double x, double y, double z, 
                           double u, double v, double w,                           
                           double particle_energy, double particle_weight, 
                           double track_length, double total_cross_section,
                           int cell_id) 
{
    // Test whether an error condition has occurred for this event
    bool errflag = false;

    // Set the particle state object
    event.particle = particle;
    event.position  = moab::CartVect(x, y, z);
    event.direction = moab::CartVect(u, v, w);
    // This should already be normalized
    event.direction.normalize();
    
    if (particle_energy < 0.0)
    {
        std::cerr << "Warning: particle_energy, " << particle_energy
                  << ", cannot be less than zero." << std::endl;
        errflag = true; 
    }
    else
    {
        event.particle_energy = particle_energy;
    }

    if (particle_weight < 0.0)
    {
        std::cerr << "Warning: particle_weight, " << particle_weight
                  << ", cannot be less than zero." << std::endl;
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
    }
    else
    {
        errflag = true; 
        std::cerr << "Warning: Cannot set a tally event of type NONE." << std::endl;
    }
    if (errflag)
    {
        clearLastEvent();
    }    
    bool event_is_set = !errflag;
    return event_is_set;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyManager.cpp
