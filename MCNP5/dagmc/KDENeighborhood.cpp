// MCNP5/dagmc/KDENeighborhood.cpp

#include <cassert>
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
                                 moab::AdaptiveKDTree& tree,
                                 moab::EntityHandle& tree_root)
{
    // copy KD-Tree to this neighborhood for the given tally event
    this->tree = &tree;
    this->tree_root = tree_root;

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
std::vector<moab::EntityHandle> KDENeighborhood::get_points() const
{
    // find all the points that exist within a rectangular neighborhood region
    std::vector<moab::EntityHandle> points;
    moab::ErrorCode rval = moab::MB_SUCCESS;

    rval = points_in_box(points);
    assert(rval == moab::MB_SUCCESS);

    return points;
}
//---------------------------------------------------------------------------//
bool KDENeighborhood::point_within_max_radius(const TallyEvent& event,
                                              const moab::CartVect& point) const
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
moab::ErrorCode
    KDENeighborhood::points_in_box(std::vector<moab::EntityHandle>& points) const
{
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

    // find all leaves of the tree within the given radius
    std::vector<moab::EntityHandle> leaves;
    moab::ErrorCode rval = moab::MB_SUCCESS;

    (*tree).leaves_within_distance(tree_root, box_center, radius, leaves);
  
    if (moab::MB_SUCCESS != rval) return rval;

    // obtain the set of unique points in the box 
    std::vector<moab::EntityHandle>::iterator i; 
    moab::Interface* mb = tree->moab();
    moab::Range leaf_points;
    moab::Range::iterator j;
    moab::EntityHandle point;
    double coords[3];
  
    // iterate through the leaves
    for (i = leaves.begin(); i != leaves.end(); ++i)
    {
        leaf_points.clear();
        rval = mb->get_entities_by_type(*i, moab::MBVERTEX, leaf_points);
    
        if (moab::MB_SUCCESS != rval) return rval;
  
        // iterate through the points in each leaf  
        for (j = leaf_points.begin(); j != leaf_points.end(); ++j)
        {
            point = *j;
            rval = mb->get_coords(&point, 1, coords);
    
            if (moab::MB_SUCCESS != rval) return rval;

            // check point is in the box
            bool in_box = true;
            int k = 0;

            do
            {
                if (coords[k] < min_corner[k] || coords[k] > max_corner[k])
                {
                    in_box = false;
                }

                ++k;
            }
            while (true && k < 3);
      
            // add the point to the set if it is in the box
            if (in_box)
            {
                points.push_back(point);
            }
        }
    }
  
    // remove duplicates from points vector
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
  
    return moab::MB_SUCCESS;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDENeighborhood.cpp
