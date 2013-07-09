// MCNP5/dagmc/ParticleState.hpp

#ifndef DAGMC_PARTICLE_STATE_HPP
#define DAGMC_PARTICLE_STATE_HPP
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

#endif // DAGMC_PARTICLE_STATE_HPP
