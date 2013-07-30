// MCNP5/dagmc/KDENeighborhood.cpp

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"

#include "KDENeighborhood.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
KDENeighborhood::KDENeighborhood(const TallyEvent& event,
                                 const moab::CartVect& bandwidth,
                                 moab::AdaptiveKDTree& kd_tree,
                                 moab::EntityHandle& kd_tree_root)
    : event(event), kd_tree(&kd_tree), kd_tree_root(kd_tree_root)
{
    // define and set the neighborhood region for this tally event
    if (event.type == TallyEvent::COLLISION)
    {
        set_neighborhood(event.position, bandwidth);
    }
    else if (event.type == TallyEvent::TRACK)
    {
        set_neighborhood(event.track_length,
                         event.position,
                         event.direction,
                         bandwidth);
    }
    else
    {
        // neighborhood region does not exist
        std::cerr << "\nError: Could not define neighborhood for tally event";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
std::set<moab::EntityHandle> KDENeighborhood::get_points() const
{
    // find all the points that exist within a rectangular neighborhood region
    return points_in_box();
}
//---------------------------------------------------------------------------//
bool KDENeighborhood::check_point_in_region(const moab::CartVect& coords) const
{
    bool in_region = false;

    // check point is in the rectangular neighborhood region
    for (int i = 0; i < 3; ++i)
    {
        // account for boundary cases first
        double min_diff = fabs(coords[i] - min_corner[i]);
        double max_diff = fabs(coords[i] - max_corner[i]);

        if (min_diff < 1e-12 || max_diff < 1e-12 ||
            coords[i] > min_corner[i] && coords[i] < max_corner[i])
        {
            in_region = true;
        }
        else // point is not in the region
        {
            return false;
        }  
    }

    return in_region;
}
//---------------------------------------------------------------------------//
bool KDENeighborhood::point_within_max_radius(const moab::CartVect& point) const
{
    // process track-based tally event only
    if (event.type == TallyEvent::TRACK)
    {
        // create a vector from starting position to point being tested
        moab::CartVect temp;

        for (int i = 0; i < 3; ++i)
        {
            temp[i] = point[i] - event.position[i];
        }

        // compute perpendicular distance from point being tested to line
        // defined by track segment using the cross-product method
        double distance_to_track = (event.direction * temp).length();

        // return true if distance is less than radius of cylindrical region
        if (distance_to_track < radius) return true;
    }

    // otherwise return false
    return false;
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
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
    radius = 0.0;
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

    // set maximum radius around the track to sqrt(hx^2 + hy^2 + hz^2)
    radius = bandwidth.length();
}                              
//---------------------------------------------------------------------------//
std::set<moab::EntityHandle> KDENeighborhood::points_in_box() const
{
    assert(kd_tree != NULL);

    std::set<moab::EntityHandle> points;

    // determine the center point of the box
    double box_center[3];

    for (int i = 0; i < 3; ++i)
    {
        box_center[i] = 0.5 * (max_corner[i] + min_corner[i]);
    }

    // set radius equal to distance from center to max corner of the box
    moab::CartVect center_to_max_corner(max_corner);
    moab::CartVect center(box_center);
    center_to_max_corner -= center;
    double radius = center_to_max_corner.length();

    // find all leaves of the kd-tree within the given radius
    std::vector<moab::EntityHandle> leaves;
    moab::ErrorCode rval = moab::MB_SUCCESS;

    rval = kd_tree->leaves_within_distance(kd_tree_root, box_center, radius, leaves);
    assert(rval == moab::MB_SUCCESS);

    // obtain the set of unique points in the box
    std::vector<moab::EntityHandle>::iterator i; 
    moab::Interface* mb = kd_tree->moab();
    moab::CartVect coords(0.0, 0.0, 0.0);
    moab::Range leaf_points;
    moab::Range::iterator j;

    for (i = leaves.begin(); i != leaves.end(); ++i)
    {
        leaf_points.clear();
        rval = mb->get_entities_by_type(*i, moab::MBVERTEX, leaf_points);
        assert(rval == moab::MB_SUCCESS);
  
        // iterate through the points in each leaf  
        for (j = leaf_points.begin(); j != leaf_points.end(); ++j)
        {
            // get coordinates of the point
            moab::EntityHandle point = *j;
            rval = mb->get_coords(&point, 1, coords.array());
            assert(rval == moab::MB_SUCCESS);

            // add the point to the set if it is in the box
            if (check_point_in_region(coords))
            {
                points.insert(point);
            }
        }
    }

    return points;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDENeighborhood.cpp
