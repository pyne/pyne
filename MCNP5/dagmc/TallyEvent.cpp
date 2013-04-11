// MCNP5/dagmc/TallyEvent.cpp

#include <cstdlib>
#include <iostream>

#include "moab/CartVect.hpp"

#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyEvent::TallyEvent()
    : particle_energy(0), particle_weight(1), tally_multiplier(1),
      position(moab::CartVect(0, 0, 0)), direction(moab::CartVect(0, 0, 0))
{
    // set this tally event to begin with no type
    event_type = NONE;
    event_value = 0;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void TallyEvent::set_track_event(const TrackData& data)
{
    if (event_type != NONE)
    {
        std::cerr << "Error: Tally event has already been set" << std::endl;
        exit(EXIT_FAILURE);
    }

    // set variables for track-based event
    event_value = data.track_length;
    position = data.start_point;
    direction = data.direction;
    particle_energy = data.particle_energy;
    particle_weight = data.particle_weight;

    // set event type
    event_type = TRACK;
  
    return;
}
//---------------------------------------------------------------------------//
void TallyEvent::set_collision_event(const CollisionData& data)
{
    if (event_type != NONE)
    {
        std::cerr << "Error: Tally event has already been set" << std::endl;
        exit(EXIT_FAILURE);
    }

    // set variables for collision event
    event_value = data.total_cross_section;
    position = data.collision_point;
    particle_energy = data.particle_energy;
    particle_weight = data.particle_weight;

    // set event type
    event_type = COLLISION;
  
    return;
}
//---------------------------------------------------------------------------//
void TallyEvent::set_tally_multiplier(double value)
{
    tally_multiplier = value;
    return;
}
//---------------------------------------------------------------------------//
// TALLY EVENT ACCESS FUNCTIONS
//---------------------------------------------------------------------------//
bool TallyEvent::get_track_data(TrackData& data) const
{
    if (event_type == TRACK)
    {
        // copy variables for track-based event
        data.track_length = event_value;
        data.start_point = position;
        data.direction = direction;
        data.particle_energy = particle_energy;
        data.particle_weight = particle_weight;

        return true;
    }
    else
    {
        // Not a valid track-based event
        return false;
    }
}
//---------------------------------------------------------------------------//
bool TallyEvent::get_collision_data(CollisionData& data) const
{
    if (event_type == COLLISION)
    {
        // copy variables for collision event
        data.total_cross_section = event_value;
        data.collision_point = position;
        data.particle_energy = particle_energy;
        data.particle_weight = particle_weight;

        return true;
    }
    else
    {
        // Not a valid collision event
        return false;
    }
}
//---------------------------------------------------------------------------//
double TallyEvent::get_tally_multiplier() const
{
    return tally_multiplier;
}
//---------------------------------------------------------------------------//
TallyEvent::EventType TallyEvent::get_event_type() const
{
    return event_type;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyEvent.cpp
