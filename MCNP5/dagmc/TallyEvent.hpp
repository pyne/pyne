// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP
#include "moab/CartVect.hpp"
 

struct TallyEvent
{
    enum EventType {NONE = 1, COLLISION = 2, TRACK = 3};

    EventType type;

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
    double particle_energy;
    double particle_weight;
};

#endif // DAGMC_TALLY_EVENT_HPP

// end of MCNP5/dagmc/TallyEvent.hpp
