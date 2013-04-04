// MCNP5/dagmc/TallyEvent.cpp

#include <iostream>

#include "moab/CartVect.hpp"

#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyEvent::TallyEvent(double energy, double weight)
    : particle_energy(energy), particle_weight(weight), tally_multiplier(1)
{
    // set this tally event to begin with no type
    event_type = NONE;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void TallyEvent::set_track_event(double track_length,
                                 moab::CartVect& start_point,
                                 moab::CartVect& direction)
{
    if (event_type != NONE)
    {
        std::cout << "Error: Tally event has already been set" << std::endl;
        return;
    }

    // set variables for track-based event
    event_value = track_length;
    position = start_point;
    this->direction = direction;

    // set event type
    event_type = TRACK;
  
    return;
}
//---------------------------------------------------------------------------//
void TallyEvent::set_collision_event(double total_cross_section,
                                     moab::CartVect& collision_point)
{
    if (event_type != NONE)
    {
        std::cout << "Error: Tally event has already been set" << std::endl;
        return;
    }

    // set variables for collision event
    event_value = total_cross_section;
    position = collision_point;
    direction = moab::CartVect(0, 0, 0);  // not used for collision events

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

// end of MCNP5/dagmc/TallyEvent.cpp
