// MCNP5/dagmc/KDENeighborhood.cpp

#include <cstdlib>
#include <iostream>

#include "KDENeighborhood.hpp"
#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
KDENeighborhood::KDENeighborhood(const TallyEvent& event,
                                 const moab::CartVect& bandwidth)
{
    // define the neighborhood region for this tally event
    TallyEvent::EventType type = event.get_event_type();

    if (type == TallyEvent::COLLISION)
    {
        // set neighborhood region for collision event
        CollisionData data;
        event.get_collision_data(data);
        set_neighborhood(data.collision_point, bandwidth);
    }
    else if (type == TallyEvent::TRACK)
    {
        // set neighborhood region for the track-based event
        TrackData data;
        event.get_track_data(data);
        set_neighborhood(data.track_length,
                         data.start_point,
                         data.direction,
                         bandwidth);
    }
    else
    {
        // Neighborhood region does not exist
        std::cerr << "\nError: Could not define neighborhood for tally event";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void KDENeighborhood::get_calculation_points()
{
// TODO implement this function when points_in_box exists here
}
//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
void KDENeighborhood::set_neighborhood(const moab::CartVect& collision_point,
                                       const moab::CartVect& bandwidth)
{
    for (int i = 0; i < 3; ++i)
    {
        min_corner[i] = collision_point[i] - bandwidth[i];
        max_corner[i] = collision_point[i] + bandwidth[i];
    }

    // maximum radius is not used for collision events
    radius = 0;
}
//---------------------------------------------------------------------------//
void KDENeighborhood::set_neighborhood(double track_length,
                                       const moab::CartVect& start_point,
                                       const moab::CartVect& direction,
                                       const moab::CartVect& bandwidth)
{
    for (int i = 0; i < 3; ++i)
    {
        // default case where coordinate of direction vector is zero
        min_corner[i] = start_point[i] - bandwidth[i];
        max_corner[i] = start_point[i] + bandwidth[i];

        // adjust for direction being positive or negative
        if (direction[i] > 0)
        {
            max_corner[i] += track_length * direction[i];
        }
        else if (direction[i] < 0)
        {
            min_corner[i] += track_length * direction[i];
        }
    }

    // set maximum radius around the track to hx^2 + hy^2 + hz^2
    radius = bandwidth.length();
}                              
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDENeighborhood.cpp
