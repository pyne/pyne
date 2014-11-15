// MCNP5/dagmc/TallyEvent.hpp

#ifndef DAGMC_TALLY_EVENT_HPP
#define DAGMC_TALLY_EVENT_HPP

#include "moab/CartVect.hpp"
#include <vector>

//===========================================================================//
/**
 * \struct TallyEvent
 * \brief Data structure for storing event data to be tallied
 *
 * Data stored in this struct is set by TallyManager and read by the different
 * Tally implementations.  The variables that will be used by each Tally depend
 * on the type of event that is defined.  All events set type, particle,
 * particle_energy and particle_weight. Collision events add the current_cell,
 * position (i.e. collision point) and total_cross_section.  Track events add
 * the current_cell, position (i.e. start of track), direction and track_length.
 *
 * An optional tally multipliers vector is also stored in TallyEvent.  Each
 * Tally can set a multiplier_id in TallyInput through the TallyManager that
 * corresponds to an index in this vector.  This multiplier_id can then be
 * used with get_score_multiplier() to return the product of the appropriate
 * multiplier and the particle_weight.
 */
//===========================================================================//
struct TallyEvent
{
    /**
     * \brief Defines type of tally event
     *
     *     0) NONE indicates no event has been set yet
     *     1) COLLISION indicates a collision event has been set
     *     2) TRACK indicates a track-based event has been set
     */
    enum EventType {NONE = 0, COLLISION = 1, TRACK = 2};

    EventType type;
 
    /// Type of particle being tallied: NEUTRON = 1, PHOTON = 2, ELECTRON = 3.
    unsigned int particle;
    
    /// Geometric cell in which the event occurred
    int current_cell;

    /// Total length of track segment
    double track_length;

    /// Position of particle (x, y, z)
    moab::CartVect position;

    /// Direction in which particle is traveling (u, v, w)
    moab::CartVect direction;

    /// Total macroscopic cross section for cell in which collision occurred
    double total_cross_section;

    /// Energy and weight of particle when event occurred
    double particle_energy;
    double particle_weight;

    /// Energy-dependent tally multipliers: variable with each event
    std::vector<double> multipliers; 

    /**
     * \brief returns multiplier * particle_weight for the current tally event
     * \param[in] multiplier_index the index of the multipliers vector to access
     * \return the score multiplier
     *
     * Note that if the multiplier_index is invalid or the tally does not use
     * multipliers, then only the particle weight will be returned.
     */
    double get_score_multiplier(int multiplier_index) const
    {
        int size = multipliers.size();

        if (multiplier_index <= -1 || multiplier_index >= size)
        {
            return particle_weight;
        }

        return multipliers.at(multiplier_index) * particle_weight; 
    }
};

#endif // DAGMC_TALLY_EVENT_HPP

// end of MCNP5/dagmc/TallyEvent.hpp
