// MCNP5/dagmc/tools/boundary_correction/boundary.cpp

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <DagMC.hpp>
#include <moab/CartVect.hpp>
#include <moab/Core.hpp>
#include <moab/Interface.hpp>

// set up MOAB and DAGMC instances
moab::Interface* mbi = new moab::Core();
moab::DagMC* dagmc = moab::DagMC::instance();

// >>> HELPER METHODS

// loads mesh data from mesh_input_file into MOAB instance
void load_moab_instance(char* mesh_input_file)
{
    moab::ErrorCode mb_error = mbi->load_file(mesh_input_file);

    if (mb_error != moab::MB_SUCCESS)
    {
        std::cerr << "Error: Could not load mesh data from input file: "
                  << mesh_input_file << std::endl;
        exit(EXIT_FAILURE);
    }
}

// loads geometry data from geometry_input_file into DAGMC instance
void load_dagmc_instance(char* geometry_input_file)
{
    moab::ErrorCode dagmc_error = dagmc->load_file(geometry_input_file);

    if (dagmc_error != moab::MB_SUCCESS)
    {
        std::cerr << "Error: DAGMC failed to load the geometry input file: "
                  << geometry_input_file << std::endl;
        exit(EXIT_FAILURE);
    }

    // initialize OBB tree and implicit complement for loaded geometry
    dagmc_error = dagmc->init_OBBTree();

    if (dagmc_error != moab::MB_SUCCESS)
    {
        std::cerr << "Error: DAGMC failed to initialize the loaded geometry\n";
        exit(EXIT_FAILURE);
    }

    // define and parse valid metadata properties
    std::vector<std::string> dagmc_keywords;
    dagmc_keywords.push_back("spec.reflect");
    dagmc_keywords.push_back("graveyard");
    dagmc_error = dagmc->parse_properties(dagmc_keywords);

    if (dagmc_error != moab::MB_SUCCESS)
    {
        std::cerr << "Error: DAGMC failed to parse metadata properties\n";
        exit(EXIT_FAILURE);
    }
}

// returns volume in which coords are located, based on direction
moab::EntityHandle find_volume_ID(const moab::CartVect& coords,
                                  const moab::CartVect& direction)
{
    // get number of 3D volumes from DAGMC instance
    int num_volumes = dagmc->num_entities(3);

    // iterate through volumes to find the one in which coords is located
    moab::EntityHandle volume_ID = 0;
    int inside_volume = 0;
    int volume_index = 1;
    
    while (inside_volume == 0 && volume_index <= num_volumes)
    {
        volume_ID = dagmc->entity_by_index(3, volume_index);
        moab::ErrorCode dagmc_error = dagmc->point_in_volume(volume_ID,
                                                             coords.array(),
                                                             inside_volume,
                                                             direction.array());

        assert(dagmc_error == moab::MB_SUCCESS);

        ++volume_index;
    }

    if (inside_volume == 0)
    {
        std::cerr << "Warning: coords (" << coords[0] << "," << coords[1] << ","
                  << coords[2] << ") do not exist in any DAGMC volume.\n";
        return 0;
    }

    return volume_ID;
}

// returns distance to boundary from coords for given direction and volume_ID
// will return -1.0 if distance is > max_distance or no intersection exists
double distance_to_boundary(double max_distance,
                            const moab::CartVect& coords,
                            const moab::CartVect& direction,
                            const moab::EntityHandle& volume_ID)
{
    assert(max_distance > 0.0);

    // create RayHistory to store surface crossing history
    moab::DagMC::RayHistory surface_history;

    // call ray fire until a reflecting or vacuum boundary is found
    moab::CartVect position = coords;
    moab::EntityHandle current_volume = volume_ID;
    moab::EntityHandle next_surface = 0;
    bool found_boundary = false;
    double distance = 0.0;
    double next_surface_distance = 0.0;

    do
    {
        // get distance to next boundary; ray_fire will only return next
        // surface if less than or equal to max_distance
        moab::ErrorCode dagmc_error = moab::MB_SUCCESS;

        dagmc_error = dagmc->ray_fire(current_volume,
                                      position.array(),
                                      direction.array(),
                                      next_surface,
                                      next_surface_distance,
                                      &surface_history,
                                      max_distance);

        assert(dagmc_error == moab::MB_SUCCESS);

        // return if no intersection is found less than max_distance
        if (next_surface == 0) return -1.0;

        // otherwise update distance to boundary and current position
        distance += next_surface_distance;
        position += next_surface_distance * direction;

        // update the current volume to the next volume
        dagmc_error = dagmc->next_vol(next_surface, current_volume, current_volume);
        assert(dagmc_error == moab::MB_SUCCESS);

        if (dagmc->has_prop(next_surface, "spec.reflect") ||
            dagmc->has_prop(current_volume, "graveyard"))
        {
            found_boundary = true;
        }
    }
    while (found_boundary == false);

    return distance;
}

// updates volume ID if on a boundary to the next volume ID
bool get_adjacent_volume(const moab::CartVect& coords,
                         const moab::CartVect& direction,
                         moab::EntityHandle& volume_ID)
{
    moab::ErrorCode dagmc_error = moab::MB_SUCCESS;
    moab::EntityHandle next_surface = 0;
    double next_surface_distance = 0.0;    

    // check distance to next surface in current direction    
    dagmc_error = dagmc->ray_fire(volume_ID,
                                  coords.array(),
                                  direction.array(),
                                  next_surface,
                                  next_surface_distance);

    assert(dagmc_error == moab::MB_SUCCESS);

    if (next_surface != 0 && next_surface_distance < 1e-10)
    {
        dagmc_error = dagmc->next_vol(next_surface, volume_ID, volume_ID);
        assert(dagmc_error == moab::MB_SUCCESS);
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * Determines if node is less than or equal to one bandwidth from any external
 * boundary in the DAGMC geometry.  If true, the boundary array will be set to
 * 0 for a lower boundary OR 1 for an upper boundary.  The distance array will
 * reflect how far that boundary is from the node.  If false, then the boundary
 * and distance arrays may return meaningless results.
 */
bool node_near_boundary(const moab::EntityHandle& mesh_node,
                        const moab::CartVect& bandwidth,
                        int boundary[3],
                        double distance[3])
{
    // get coordinates of the mesh node
    moab::CartVect coords;
    moab::ErrorCode mb_error = mbi->get_coords(&mesh_node, 1, coords.array());
    assert(mb_error == moab::MB_SUCCESS);

    // determine volume in DAGMC geometry in which mesh node is located
    moab::CartVect direction(1.0, 0.0, 0.0);
    moab::EntityHandle volume_ID = find_volume_ID(coords, direction);

    if (dagmc->has_prop(volume_ID, "graveyard"))
    {
        bool volume_change = false;

        // check if coords are on inner graveyard boundary and update volume_ID
        for (int i = 0; i < 3; ++i)
        {
            if (!volume_change)
            {
                direction[i] = 1.0;
                volume_change = get_adjacent_volume(coords, direction, volume_ID);
            }
            
            if (!volume_change)
            {
                direction[i] = -1.0;
                volume_change = get_adjacent_volume(coords, direction, volume_ID);
            }

            direction[i] = 0.0;
        }

        if (!volume_change)
        {
            std::cerr << "Warning: node " << coords
                      << " was found inside graveyard" << std::endl;
            return false;
        }
    }

    if (volume_ID == 0) return false;

    // iterate through all directions x = 0, y = 1, and z = 2
    bool boundary_point = false;

    for (int i = 0; i < 3; ++i)
    {
        // compute distance to upper boundary
        direction[i] = 1.0;
        distance[i] = distance_to_boundary(bandwidth[i],
                                           coords,
                                           direction,
                                           volume_ID);

        // test if mesh node is NOT within one bandwidth of upper boundary
        if (distance[i] < 0)
        {
            // compute distance to lower boundary
            direction[i] = -1.0;
            distance[i] = distance_to_boundary(bandwidth[i],
                                               coords,
                                               direction,
                                               volume_ID);

            // test if mesh node is NOT within one bandwidth of lower boundary
            if (distance[i] < 0)
            {
                boundary[i] = -1;
            }
            else // mesh node is within one bandwidth of lower boundary
            {
                boundary_point = true;
                boundary[i] = 0;
            }
        }
        else // mesh node is within one bandwidth of upper boundary
        {
            boundary_point = true;
            boundary[i] = 1;
        }

        // reset direction for next iteration
        direction[i] = 0.0;
    }

    return boundary_point;
}

// converts a C-string value into a valid bandwidth value
double parse_bandwidth_value(char* value)
{
    char* end;
    double bandwidth_value = strtod(value, &end);

    if (value != end && bandwidth_value > 0.0)
    {
        return bandwidth_value;
    }
    else
    {
        std::cerr << "Warning: bandwidth = " << value << " is invalid\n";
        std::cerr << "    using default value of 0.1" << std::endl;
        return 0.1;
    }
}

// >>> MAIN METHOD

int main(int argc, char* argv[])
{
    if (argc != 6)
    {
        std::cerr << "usage: " << argv[0] << " <input mesh> <input geometry> "
                  << "<hx hy hz>" << std::endl;
        exit( EXIT_FAILURE );
    }

    // load the mesh data into the MOAB instance
    load_moab_instance(argv[1]);

    // load the geometry data into the DAGMC instance
    load_dagmc_instance(argv[2]);

    // read bandwidth vector from command line
    moab::CartVect bandwidth;
    
    for (int i = 0; i < 3; ++i)
    {
        bandwidth[i] = parse_bandwidth_value(argv[i+3]);
    }

    std::cout << "Using a bandwidth value of " << bandwidth
              << " to determine boundary points" << std::endl;
    
    // create new BOUNDARY and DISTANCE_TO_BOUNDARY tags
    int tag_size = 3;
    const std::vector<int> default_boundary(3, -1);
    const moab::CartVect default_distance(-1.0);
    moab::ErrorCode mb_error = moab::MB_SUCCESS;
    moab::Tag boundary_tag, distance_tag;

    mb_error = mbi->tag_get_handle("BOUNDARY", tag_size,
                                   moab::MB_TYPE_INTEGER,
                                   boundary_tag,
                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
                                   &default_boundary[0]);

    assert(mb_error == moab::MB_SUCCESS);

    mb_error = mbi->tag_get_handle("DISTANCE_TO_BOUNDARY", tag_size,
                                   moab::MB_TYPE_DOUBLE,
                                   distance_tag,
                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
                                   default_distance.array());
                                   
    assert(mb_error == moab::MB_SUCCESS);

    // obtain all mesh nodes loaded into the MOAB instance
    // TODO consider using Skinner to reduce set of mesh nodes to check
    moab::EntityHandle mesh_set = 0;
    std::vector<moab::EntityHandle> mesh_nodes;

    mb_error = mbi->get_entities_by_type(mesh_set, moab::MBVERTEX, mesh_nodes);
    assert(mb_error == moab::MB_SUCCESS);

    // iterate through all mesh nodes
    std::vector<moab::EntityHandle>::iterator it;
    int node_count = 0;

    for(it = mesh_nodes.begin(); it != mesh_nodes.end(); ++it)
    {
        moab::EntityHandle node = *it;
        int boundary[3];
        double distance[3];

        // determine if mesh node is near any external boundaries
        bool boundary_point = node_near_boundary(node, bandwidth, boundary, distance);

        // if boundary point, set BOUNDARY and DISTANCE_TO_BOUNDARY tag data
        if (boundary_point)
        {
            mb_error = mbi->tag_set_data(boundary_tag, &node, 1, boundary);
            assert(mb_error == moab::MB_SUCCESS);

            mb_error = mbi->tag_set_data(distance_tag, &node, 1, distance);
            assert(mb_error == moab::MB_SUCCESS);
            ++node_count;
        }
    }

    std::cout << node_count << " of " << mesh_nodes.size()
              << " were tagged as boundary points" << std::endl;

    // write mesh to output file with new tag data
    moab::Tag output_tags[2] = {boundary_tag, distance_tag};

    mb_error = mbi->write_file("output_mesh.vtk",
                               NULL, NULL, NULL, 0,
                               output_tags, 2);

    // Memory management
    delete mbi;
    // DAGMC has no destructor

    return 0;
}

// end of MCNP5/dagmc/tools/boundary_correction/boundary.cpp
