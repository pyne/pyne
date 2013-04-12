// MCNP5/dagmc/KDENeighborhood.hpp

#ifndef DAGMC_KDE_NEIGHBORHOOD_H
#define DAGMC_KDE_NEIGHBORHOOD_H

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"

#include "TallyEvent.hpp"

/**
 * \class KDENeighborhood
 * \brief Defines a neighborhood region for a KDE mesh tally event
 *
 * KDENeighborhood is a class that defines an approximation of the
 * neighborhood region for either a collision or track-based event that is
 * to be tallied as part of a Monte Carlo particle transport simulation.  This
 * neighborhood region is formally defined as the region in space for which the
 * kernel function produces a non-trivial result for any mesh node that exists
 * inside that region.  This set of mesh nodes is also known as the set of
 * calculation points for the KDE mesh tally.
 *
 * Once a KDENeighborhood has been created, the calculation points for the
 * corresponding tally event can be obtained using the get_points() function.
 * Note that this may still include some mesh nodes that produce trivial tally
 * contributions because KDENeighborhood is only an approximation of the true
 * neighborhood region.
 */
class KDENeighborhood
{
  public:
    /**
     * \brief Constructor
     * \param event the tally event for which the neighborhood is desired
     * \param bandwidth the bandwidth vector (hx, hy, hz)
     * \param tree the KD-Tree containing all mesh nodes in the input mesh
     * \param tree_root the root of the KD-Tree
     */
    KDENeighborhood(const TallyEvent& event,
                    const moab::CartVect& bandwidth,
                    moab::AdaptiveKDTree& tree,
                    moab::EntityHandle& tree_root);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Gets the calculation points for this neighborhood region
     * \param points the vector to which the calculation points are copied
     * \return the MOAB ErrorCode value
     *
     * The neighborhood region is currently assumed to be rectangular.  This
     * function will therefore return all of the calculation points that exist
     * within the boxed region defined by min_corner and max_corner.
     */
    moab::ErrorCode get_points(std::vector<moab::EntityHandle>& points);

  private:
    /// Minimum and maximum corner of a rectangular neighborhood region
    double min_corner[3];
    double max_corner[3];

    /// Radius of a cylindrical neighborhood region
    double radius;

    /// KD-Tree containing all mesh nodes in the input mesh
    moab::AdaptiveKDTree* tree;
    moab::EntityHandle tree_root;

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

    /**
     * \brief Finds the vertices that exist inside a rectangular region
     * \param points stores unique set of vertices that exist inside box
     * \return the MOAB ErrorCode value
     */
    moab::ErrorCode points_in_box(std::vector<moab::EntityHandle>& points);
};

#endif // DAGMC_KDE_NEIGHBORHOOD_H

// end of MCNP5/dagmc/KDENeighborhood.hpp
