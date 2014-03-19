// MCNP5/dagmc/KDEMeshTally.cpp

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "moab/Core.hpp"
#include "moab/Range.hpp"

#include "KDEMeshTally.hpp"

// initialize static variables
bool KDEMeshTally::seed_is_set = false;
const char* const KDEMeshTally::kde_estimator_names[] = {"collision",
                                                         "integral-track",
                                                         "sub-track"};

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
KDEMeshTally::KDEMeshTally(const TallyInput& input,
                           KDEMeshTally::Estimator type)
    : MeshTally(input),
      estimator(type),
      bandwidth(moab::CartVect(0.01, 0.01, 0.01)),
      kernel(NULL),
      use_kd_tree(true),
      region(NULL),
      use_boundary_correction(false),
      num_subtracks(3),
      quadrature(NULL),
      mbi(new moab::Core())
{
    std::cout << "Creating KDE " << kde_estimator_names[estimator]
              << " mesh tally " << input.tally_id << std::endl;

    std::cout << "    for input mesh: " << input_filename
              << ", output file: " << output_filename << std::endl;

    // set up KDEMeshTally member variables from TallyInput
    parse_tally_options();
    assert(kernel != NULL);
 
    std::cout << "    using " << kernel->get_kernel_name()
              << " kernel and bandwidth " << bandwidth << std::endl;

    // turn off boundary correction for higher-order kernels
    if (use_boundary_correction == true && kernel->get_order() > 2)
    {
        std::cerr << "Warning: boundary correction is not yet supported "
                  << "for higher-order kernels" << std::endl;
        use_boundary_correction = false;
    }

    if (estimator == INTEGRAL_TRACK)
    {
        // set up quadrature rule for the integral_track estimator
        // NOTE: this will only work correctly for polynomial kernel functions
        int num_points = 3 * kernel->get_min_quadrature(0) - 2;
        std::cout << "    using " << num_points << "-pt quadrature scheme\n";
        quadrature = new Quadrature(num_points);
    }
    else if (estimator == SUB_TRACK)
    {
        std::cout << "    splitting full tracks into "
                  << num_subtracks << " sub-tracks" << std::endl;

        // set seed value if not already set by another instance
        if (!seed_is_set)
        {
            srand(time(NULL));
            seed_is_set = true;
        }
    }

    // initialize MeshTally member variables representing the mesh data
    moab::ErrorCode rval = initialize_mesh_data();

    if (rval != moab::MB_SUCCESS)
    {
        std::cerr << "Error: Could not load mesh data for KDE mesh tally "
                  << input.tally_id << " from input file "
                  << input_filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // initialize running variance variables
    max_collisions = false;
    num_collisions = 0;
    mean = moab::CartVect(0, 0, 0);
    variance = moab::CartVect(0, 0, 0);
}
//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
KDEMeshTally::~KDEMeshTally()
{
    delete kernel;
    delete region;
    delete quadrature;
    delete mbi;
}
//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from Tally.hpp
//---------------------------------------------------------------------------//
void KDEMeshTally::compute_score(const TallyEvent& event)
{
    double weight = event.get_score_multiplier(input_data.multiplier_id);

    // set up tally event based on KDE mesh tally type
    std::vector<moab::CartVect> subtrack_points;

    if (event.type == TallyEvent::TRACK && estimator != COLLISION)
    {
        if (estimator == SUB_TRACK)
        {
            // multiply weight by track length and set up sub-track points
            weight *= event.track_length;
            subtrack_points = choose_points(num_subtracks, event);
        }
    }
    else if (event.type == TallyEvent::COLLISION && estimator == COLLISION)
    {
        // divide weight by cross section and update optimal bandwidth
        weight /= event.total_cross_section;
        update_variance(event.position);
    }
    else // NONE, return from this method
    {
 	return;
    }

    unsigned int ebin;
    if (!get_energy_bin(event.particle_energy, ebin))
    {  
        return;
    }

    // update the neighborhood region and find all of the calculations points
    region->update_neighborhood(event, bandwidth);
    std::set<moab::EntityHandle> calculation_points = region->get_points();

    // iterate through calculation points and compute their final scores
    std::set<moab::EntityHandle>::iterator i;
    CalculationPoint X;

    for (i = calculation_points.begin(); i != calculation_points.end(); ++i)
    {
        // get coordinates of this point
        moab::EntityHandle point = *i;
        moab::ErrorCode rval = mbi->get_coords(&point, 1, X.coords);
	if(rval != moab::MB_SUCCESS)
	  {
	    std::cout << "Failed to get coordinates" << std::endl;
	    exit(1);
	  }
        assert(rval == moab::MB_SUCCESS);

        // get tag data for this point if user requested boundary correction
        if (use_boundary_correction)
        {
            rval = mbi->tag_get_data(boundary_tag, &point, 1, X.boundary_data);
            assert(rval == moab::MB_SUCCESS);

            rval = mbi->tag_get_data(distance_tag, &point, 1, X.distance_data);
            assert(rval == moab::MB_SUCCESS);
        }

        // compute the final contribution to the tally for this point
        double score = 0.0;

        if (estimator == INTEGRAL_TRACK)
        {
            score = integral_track_score(X, event);
        }
        else if (estimator == SUB_TRACK)
        {
            score = subtrack_score(X, subtrack_points);
        }
        else // estimator == COLLISION
        {
            score = evaluate_kernel(X, event.position);
        }

        add_score_to_mesh_tally(point, weight, score,  ebin);
    }  // end calculation_points iteration
}
//---------------------------------------------------------------------------//
void KDEMeshTally::write_data(double num_histories)
{
    // display the optimal bandwidth if it was computed
    if (estimator == COLLISION)
    {
        std::cout << std::endl << "optimal bandwidth for " << num_collisions
                  << " collisions is: " << get_optimal_bandwidth() << std::endl;
    }

    // tag tally and relative error results to the mesh for each tally point
    moab::ErrorCode rval = moab::MB_SUCCESS;
    if(rval != moab::MB_SUCCESS)
      {
	std::cout << "This can't fail" << std::endl;
	exit(1);
      }
    moab::Range::iterator i;

    for (i = tally_points.begin(); i != tally_points.end(); ++i)
    {
        moab::EntityHandle point = *i;
        unsigned int point_index = get_entity_index(point);

        for (unsigned int j = 0; j < data->get_num_energy_bins(); ++ j)
        {
            std::pair <double,double> tally_data = data->get_data(point_index,j);
            double tally = tally_data.first;
            double error = tally_data.second;

            // compute relative error for the tally result
            double rel_error = 0.0;

            if (error != 0.0)
            {
                rel_error = sqrt(error / (tally * tally) - 1.0 / num_histories);
            }

            // normalize mesh tally result by the number of source particles
            tally /= num_histories;

            // set tally and error tag values for this entity
            rval = mbi->tag_set_data(tally_tags[j], &point, 1, &tally);

            assert(moab::MB_SUCCESS == rval);
      
            rval = mbi->tag_set_data(error_tags[j], &point, 1, &rel_error);

            assert(moab::MB_SUCCESS == rval);
        }
    }

    // create a global tag to store the bandwidth value
    moab::Tag bandwidth_tag;
    rval = mbi->tag_get_handle("BANDWIDTH_TAG", 3,
                               moab::MB_TYPE_DOUBLE,
                               bandwidth_tag,
                               moab::MB_TAG_MESH|moab::MB_TAG_CREAT);

    assert(moab::MB_SUCCESS == rval);

    // add bandwidth tag to the root set
    moab::EntityHandle bandwidth_set = mbi->get_root_set();
    rval = mbi->tag_set_data(bandwidth_tag, &bandwidth_set, 1, &bandwidth);

    assert(moab::MB_SUCCESS == rval);

    // define list of tags to include and write mesh to output file
    std::vector<moab::Tag> output_tags = tally_tags;
    output_tags.insert(output_tags.end(), error_tags.begin(), error_tags.end());
    output_tags.push_back(bandwidth_tag);

    rval = mbi->write_file(output_filename.c_str(),
                           NULL, NULL,
                           &tally_mesh_set, 1,
                           &(output_tags[0]),
                           output_tags.size());

    assert(moab::MB_SUCCESS == rval);
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
void KDEMeshTally::set_bandwidth_value(const std::string& key,
                                       const std::string& value,
                                       unsigned int i)
{
    // extract the bandwidth value for the given key
    char* end;
    double bandwidth_value = strtod(value.c_str(), &end);

    // set the bandwidth value for the given index
    if (value.c_str() != end && bandwidth_value > 0.0)
    {
        bandwidth[i] = bandwidth_value;
    }
    else
    {
        std::cerr << "Warning: invalid bandwidth value for " << key
                  << " = " << value << std::endl;
        std::cerr << "    using default value " << key << " = "
                  << bandwidth[i] << std::endl;
    }
}
//---------------------------------------------------------------------------//
void KDEMeshTally::parse_tally_options()
{
    const TallyInput::TallyOptions& options = input_data.options;  
    TallyInput::TallyOptions::const_iterator it;
    std::string kernel_type = "epanechnikov";
    int kernel_order = 2;

    for (it = options.begin(); it != options.end(); ++it)
    {
        std::string key = it->first;
        std::string value = it->second;

        // process tally option according to key
        if      (key == "hx") set_bandwidth_value(key, value, 0);
        else if (key == "hy") set_bandwidth_value(key, value, 1);
        else if (key == "hz") set_bandwidth_value(key, value, 2);
        else if (key == "kernel")
        {
            kernel_type = value;
        }
        else if (key == "order")
        {
            char* end;
            kernel_order = strtol(value.c_str(), &end, 10);

            if (value.c_str() == end || kernel_order % 2 != 0)
            {
                std::cerr << "Warning: '" << value << "' is an invalid value"
                          << " for the kernel order" << std::endl;
                std::cerr << "    using default value " << key << " = 2\n";
                kernel_order = 2;
            }
        }
        else if (key == "neighborhood" && value == "off")
        {
            std::cout << "    using neighborhood-search: " << value << std::endl;
            use_kd_tree = false;
        }
        else if (key == "boundary" && value == "default")
        {
            std::cout << "    using boundary correction: " << value << std::endl;
            use_boundary_correction = true;
        }
        else if (key == "seed" && estimator == SUB_TRACK)
        {
            // override random number seed if requested by user
            unsigned long int seed = strtol(value.c_str(), NULL, 10);
            srand(seed);
            seed_is_set = true;
            std::cout << "    setting random seed to " << seed
                      << " for choosing sub-track points" << std::endl;
        }
        else if (key == "subtracks" && estimator == SUB_TRACK)
        {
            char* end;
            int subtracks = strtol(value.c_str(), &end, 10);

            if (value.c_str() == end || subtracks <= 0)
            {
                std::cerr << "Warning: '" << value << "' is an invalid value"
                          << " for the number of subtracks" << std::endl;
                std::cerr << "    using default value " << key << " = 3\n";
                subtracks = 3;
            }
            
            num_subtracks = subtracks;
        }
        else // invalid tally option
        {
            std::cerr << "Warning: input data for KDE mesh tally "
                      << input_data.tally_id
                      << " has unknown key '" << key << "'" << std::endl;
        }
    }

    // create a kernel function to use for this KDEMeshTally
    kernel = KDEKernel::createKernel(kernel_type, kernel_order);

    if (kernel == NULL)
    {
        std::cerr << "Warning: invalid kernel type '" << kernel_type << "'\n";
        kernel = KDEKernel::createKernel("epanechnikov", kernel_order);
    }
}
//---------------------------------------------------------------------------//
moab::ErrorCode KDEMeshTally::initialize_mesh_data()
{
    // load the MOAB mesh data from the input file for this KDE mesh tally
    moab::EntityHandle moab_mesh_set;
    moab::ErrorCode rval = load_moab_mesh(mbi, moab_mesh_set);

    if (rval != moab::MB_SUCCESS) return rval;

    // get all of the mesh nodes from the MOAB mesh set
    moab::Range mesh_nodes;
    rval = mbi->get_entities_by_type(moab_mesh_set, moab::MBVERTEX, mesh_nodes);

    if (rval != moab::MB_SUCCESS) return rval;

    // initialize MeshTally::tally_points to include all mesh nodes
    set_tally_points(mesh_nodes);

    // set up the KDE neighborhood region from the mesh nodes
    region = new KDENeighborhood(mbi, mesh_nodes, use_kd_tree);

    // reduce the loaded MOAB mesh set to include only 3D elements
    moab::Range mesh_cells;
    rval = reduce_meshset_to_3D(mbi, moab_mesh_set, mesh_cells);  

    if (rval != moab::MB_SUCCESS) return rval;

    // initialize MeshTally::tally_mesh_set and set up tags for energy bins
    tally_mesh_set = moab_mesh_set;
    rval = setup_tags(mbi, "KDE_");

    if (rval != moab::MB_SUCCESS) return rval;

    // if requested, set up tags for boundary correction method
    if (use_boundary_correction)
    {
        moab::ErrorCode tag1 = mbi->tag_get_handle("BOUNDARY", 3,
                                                   moab::MB_TYPE_INTEGER,
                                                   boundary_tag);

        moab::ErrorCode tag2 = mbi->tag_get_handle("DISTANCE_TO_BOUNDARY", 3,
                                                   moab::MB_TYPE_DOUBLE,
                                                   distance_tag);

        if (tag1 == moab::MB_TAG_NOT_FOUND || tag2 == moab::MB_TAG_NOT_FOUND)
        {
            std::cerr << "Warning: no valid boundary tags were found\n"
                      << "    ignoring request for boundary correction\n";
            use_boundary_correction = false;
        }
        else if (tag1 != moab::MB_SUCCESS || tag2 != moab::MB_SUCCESS)
        {
            return moab::MB_FAILURE;
        }
    }

    return moab::MB_SUCCESS; 
}
//---------------------------------------------------------------------------//
void KDEMeshTally::update_variance(const moab::CartVect& collision_point)
{
    if (num_collisions != LLONG_MAX)
    {
        ++num_collisions;

        // compute new values for the mean and variance
        if (num_collisions == 1)
        {
            mean = collision_point;
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                // get difference between point and previous mean
                double value = collision_point[i] - mean[i];

                // update mean and variance variables
                mean[i] += value / num_collisions;
                variance[i] += value * (collision_point[i] - mean[i]);
            }
        }
    }
    else if (!max_collisions)
    {
        std::cerr << "Warning: number of collisions exceeds maximum\n"
                  << "    optimal bandwidth will be based on "
                  << num_collisions << " collisions" << std::endl;

        max_collisions = true;
    }
}
//---------------------------------------------------------------------------//
moab::CartVect KDEMeshTally::get_optimal_bandwidth() const
{
    double stdev = 0.0;
    moab::CartVect optimal_bandwidth;

    for (int i = 0; i < 3; ++i)
    {
        stdev = sqrt(variance[i] / (num_collisions - 1));
        optimal_bandwidth[i] = 0.968625 * stdev * pow(num_collisions, -1.0/7.0);
    }

    return optimal_bandwidth;
}
//---------------------------------------------------------------------------//
// KDE ESTIMATOR METHODS
//---------------------------------------------------------------------------//
double KDEMeshTally::PathKernel::evaluate(double s) const
{
    moab::CartVect observation = event.position + s * event.direction;  
    return tally.evaluate_kernel(X, observation);
}
//---------------------------------------------------------------------------//
double KDEMeshTally::evaluate_kernel(const CalculationPoint& X,
                                     const moab::CartVect& observation) const
{
    // evaluate the 3D kernel function
    double kernel_value = 1.0;

    for (int i = 0; i < 3; ++i)
    {
        double u = (X.coords[i] - observation[i]) / bandwidth[i];

        if (use_boundary_correction && X.boundary_data[i] != -1)
        {
            // use boundary kernel for this direction
            kernel_value *= kernel->evaluate(u, bandwidth[i],
                                             X.distance_data[i],
                                             X.boundary_data[i]) / bandwidth[i];
        }
        else // use standard kernel for this direction
        {
            kernel_value *= kernel->evaluate(u) / bandwidth[i];
        }
    }

    return kernel_value;
}                    
//---------------------------------------------------------------------------//
double KDEMeshTally::integral_track_score(const CalculationPoint& X,
                                          const TallyEvent& event) const
{
    // determine the limits of integration
    std::pair<double, double> limits;  
    bool valid_limits = set_integral_limits(event,
                                            moab::CartVect(X.coords),
                                            limits);

    // compute value of the integral only if valid limits exist
    if (valid_limits)
    {
        // construct a PathKernel and return value of its integral
        PathKernel path_kernel(*this, event, X);
        return quadrature->integrate(limits.first, limits.second, path_kernel);
    }
    else // integration limits are not valid so no score is computed
    {
        return 0.0;
    }
}
//---------------------------------------------------------------------------//
bool KDEMeshTally::set_integral_limits(const TallyEvent& event,
                                       const moab::CartVect& coords,
                                       std::pair<double, double>& limits) const
{
    bool valid_limits = false;

    // set initial integral limits to the full track length (default values)
    limits = std::make_pair(0, event.track_length);

    // check limits against the valid path length interval for each dimension
    for (int i = 0; i < 3; ++i)
    {
        double path_min = limits.first;
        double path_max = limits.second;

        // compute valid path length interval Si = [path_min, path_max]
        if (event.direction[i] > 0)
        {
            path_min = coords[i] - event.position[i] - bandwidth[i];
            path_min /= event.direction[i];

            path_max = coords[i] - event.position[i] + bandwidth[i];
            path_max /= event.direction[i];
        }
        else if (event.direction[i] < 0)
        {
            path_min = coords[i] - event.position[i] + bandwidth[i];
            path_min /= event.direction[i];

            path_max = coords[i] - event.position[i] - bandwidth[i];
            path_max /= event.direction[i];
        }

        // set lower limit to highest minimum
        if (path_min > limits.first)
        {
            limits.first = path_min;
        }

        // set upper limit to lowest maximum
        if (path_max < limits.second)
        {
            limits.second = path_max;
        }
    }
  
    // limits are only valid if upper limit is greater than lower limit
    if (limits.first < limits.second)
    {
        valid_limits = true;
    }

    return valid_limits;
}
//---------------------------------------------------------------------------//
double KDEMeshTally::subtrack_score(const CalculationPoint& X,
                                    const std::vector<moab::CartVect>& points) const
{
    // iterate through the sub-track points
    std::vector<moab::CartVect>::const_iterator i;
    double score = 0.0;

    if (!points.empty())
    {
        for (i = points.begin(); i != points.end(); ++i)
        {
            // add kernel contribution for sub-track point to sum
            score += evaluate_kernel(X, *i);
        }

        // normalize by the total number of sub-track points
        score /= points.size();
    }

    return score;
}
//---------------------------------------------------------------------------//
std::vector<moab::CartVect> KDEMeshTally::choose_points(unsigned int p,
                                                        const TallyEvent& event) const
{
    // make sure the number of sub-tracks is valid
    assert(p > 0);

    // compute sub-track length, assumed to be equal for all sub-tracks
    double sub_track_length = event.track_length / p;

    // set the starting point to the beginning of the track segment
    moab::CartVect start_point = event.position;

    // choose a random position along each sub-track
    std::vector<moab::CartVect> random_points;

    for (unsigned int i = 0; i < p; ++i)
    {
        double path_length = rand() * sub_track_length / RAND_MAX;
        
        // add the coordinates of the corresponding point
        random_points.push_back(start_point + path_length * event.direction);

        // shift starting point to the next sub-track
        start_point += sub_track_length * event.direction;
    }
 
    return random_points;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDEMeshTally.cpp
