// MCNP5/dagmc/KDENeighborhood.hpp

#ifndef DAGMC_KDE_NEIGHBORHOOD_H
#define DAGMC_KDE_NEIGHBORHOOD_H

#include "moab/CartVect.hpp"

#include "TallyEvent.hpp"

/**
 * \class KDENeighborhood
 * \brief Defines a neighborhood region for a KDE mesh tally
 */
class KDENeighborhood
{
  public:
    /**
     * \brief Constructor
     * \param event the tally event for which the neighborhood is desired
     * \param bandwidth the bandwidth vector (hx, hy, hz)
     */
    KDENeighborhood(const TallyEvent& event, const moab::CartVect& bandwidth);

    // >>> PUBLIC INTERFACE

    /**
     * \brief TODO define function
     */
    void get_calculation_points();

  private:
    /// Minimum and maximum corner of a rectangular neighborhood region
    double min_corner[3];
    double max_corner[3];

    /// Radius of a cylindrical neighborhood region
    double radius;

    // >>> PRIVATE FUNCTIONS

    /**
     * \brief Sets the neighborhood region for a collision event
     * \param collision_point the location of the collision (x, y, z)
     * \param bandwidth the bandwidth vector (hx, hy, hz)
     */
    void set_neighborhood(const moab::CartVect& collision_point,
                          const moab::CartVect& bandwidth);

    /**
     * \brief Sets the neighborhood region for a track-based event
     * \param track_length the total length of the track segment
     * \param start_point the starting location of the particle (xo, yo, zo)
     * \param direction the direction the particle is traveling (uo, vo, wo)
     * \param bandwidth the bandwidth vector (hx, hy, hz)
     */
    void set_neighborhood(double track_length,
                          const moab::CartVect& start_point,
                          const moab::CartVect& direction,
                          const moab::CartVect& bandwidth);
};

#endif // DAGMC_KDE_NEIGHBORHOOD_H

// end of MCNP5/dagmc/KDENeighborhood.hpp
