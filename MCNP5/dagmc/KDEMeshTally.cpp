// KDEMeshTally.cpp

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"

#include "KDECollision.hpp"
#include "KDEKernel.hpp"
#include "KDETrack.hpp"
#include "KDEMeshTally.hpp"
#include "meshtal_funcs.h"

moab::CartVect KDEMeshTally::default_bandwidth( 1, 1, 1 );


//-----------------------------------------------------------------------------
static double parse_bandwidth_param( const std::string& key,
                                     const std::string& val, 
                                     double default_value = 1 )
{
  char* end;
  double ret = strtod( val.c_str(), &end );
  if( val.c_str() == end || ret <= 0 ){
    // parsing failed or value is invalid
    std::cerr << "Error parsing bandwidth param " << key << ": '" << val << "' is invalid" << std::endl;
    std::cerr << "      Using default value " << key << " = " << default_value << std::endl;
    ret = default_value;
  }
  return ret;
}
//-----------------------------------------------------------------------------
static KDEKernel::KernelType set_kernel_type( const std::string& key,
                                              const std::string& val )
{

  KDEKernel::KernelType k;

  if ( val == "epanechnikov" || val == "e" )
    k = KDEKernel::EPANECHNIKOV;
  else if ( val == "uniform" || val == "u" )
    k = KDEKernel::UNIFORM;
  else {

    k = KDEKernel::EPANECHNIKOV;
    std::cerr << "\nWarning: " << val << " is not a valid kernel function" << std::endl;
        std::cerr << "      Using default " << key << " = epanechnikov\n" << std::endl;

  }

  return k;

}
//-----------------------------------------------------------------------------
KDEMeshTally* KDEMeshTally::setup( const MeshTallyInput& fmesh, 
                                   moab::Interface* mbi, TallyType type )
{

  bool use_dagmc_mesh = true; // true if mesh data should be pulled from DagMC object
  std::string input_filename, output_filename;
  moab::CartVect bandwidth = default_bandwidth;
  KDEKernel::KernelType kernel = KDEKernel::EPANECHNIKOV;
  unsigned int subtracks = 0;

  const MeshTallyInput::TallyOptions& params = fmesh.options;
  int id = fmesh.tally_id;

  for( MeshTallyInput::TallyOptions::const_iterator i = params.begin(); i != params.end(); ++i )
  {
    std::string key = (*i).first, val = (*i).second;

    if( key == "inp" ){ input_filename = val; use_dagmc_mesh = false; }
    else if( key == "out" ) output_filename = val;
    else if( key == "hx" ) bandwidth[0] = parse_bandwidth_param( key, val );
    else if( key == "hy" ) bandwidth[1] = parse_bandwidth_param( key, val );
    else if( key == "hz" ) bandwidth[2] = parse_bandwidth_param( key, val );
    else if( key == "kernel" ) kernel = set_kernel_type( key, val );
    else if( key == "subtracks" && type != COLLISION ) { 
      subtracks = atoi( val.c_str() );
      if( subtracks == 0 ) {
        std::cerr << "\nWarning: number of subtracks requested is invalid" << std::endl;
        std::cerr << "      Using default value " << key << " = 3\n" << std::endl;
        subtracks = 3;
      }
    }
    else{
      std::cerr << "Warning: KDE tally's FC" << id << " card has unknown key '" << key << "'" << std::endl;
    }
  }
  
  std::stringstream strbuf;
  if( output_filename.length() == 0 ){
    strbuf << "meshtal" << id << ".h5m";
    strbuf >> output_filename;
    strbuf.clear();
  }

  std::cout << "Creating KDE ";

  if ( type == COLLISION ) std::cout << "collision ";
  else if ( type == TRACKLENGTH ) std::cout << "track ";
  else std::cout << "subtrack ";
 
  std::cout << "fmesh" << id 
            << ", input: " << (use_dagmc_mesh ? "(pre-loaded DagMC data)" : input_filename.c_str())  
            << ", output: " << output_filename << std::endl;
  
  std::cout << "   using the " << KDEKernel::kernel_names[kernel] << " kernel";
  std::cout << " with bandwidth = " << bandwidth;

  if ( subtracks != 0 )
    std::cout << "\n   and splitting tracks into " << subtracks << " subtracks" << std::endl;
  else
    std::cout << std::endl;

  moab::EntityHandle moab_set = NULL ; // TODO: this should be queried from DagMC
  if( !use_dagmc_mesh ){
   
    moab::ErrorCode rval;
    rval = mbi->create_meshset( moab::MESHSET_SET, moab_set );
    assert( rval == moab::MB_SUCCESS );

    rval = mbi->load_file( input_filename.c_str(), &moab_set );

    if( rval != moab::MB_SUCCESS ){
      std::cerr << "Error: could not load KDE tally mesh file " << input_filename << std::endl;
      exit( EXIT_FAILURE );
    }

  }

  KDEMeshTally *kde = new KDEMeshTally( fmesh, mbi, moab_set, bandwidth,
                                        type, kernel, subtracks );
  kde->output_filename = output_filename;

  // create a tally set that contains only the 3D mesh cells (i.e. hexes/tets)
  moab::ErrorCode rval;
  rval = mbi->create_meshset( moab::MESHSET_SET, kde->tally_set );
  assert( rval == moab::MB_SUCCESS );

  moab::Range mesh_cells;
  rval = mbi->get_entities_by_dimension( moab_set, 3, mesh_cells );
  assert( rval == moab::MB_SUCCESS );
   
  rval = mbi->add_entities( kde->tally_set, mesh_cells );
  assert( rval == moab::MB_SUCCESS );

  return kde;

}
//-----------------------------------------------------------------------------
KDEMeshTally::KDEMeshTally( const MeshTallyInput& settings,
                            moab::Interface* moabMesh,
                            moab::EntityHandle moabSet,
                            moab::CartVect bandwidth,
                            KDEMeshTally::TallyType type,
                            KDEKernel::KernelType k,
                            unsigned int numSubtracks )
: MeshTally(settings), mb( moabMesh ), bandwidth( bandwidth ),
  kde_tally( type ), kernel( new KDEKernel(k) )
{

  build_tree(moabSet); 

  // initialize running variance variables
  max_collisions = false;
  numCollisions = 0;
  Mn = moab::CartVect( 0, 0, 0 );
  Sn = moab::CartVect( 0, 0, 0 );

  // check numSubtracks is a valid parameter for SUBTRACK tallies
  if ( numSubtracks == 0 && kde_tally == SUBTRACK ) {

    std::cerr << "\nWarning: number of subtracks must be non-zero for KDE subtrack tallies" << std::endl;
    std::cerr << "      Using default value subtracks = 3\n" << std::endl;
    subtracks = 3;

  }
  else
    subtracks = numSubtracks;
  
}
//-----------------------------------------------------------------------------
KDEMeshTally::~KDEMeshTally()
{

  delete tree;
  delete kernel;
  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::tally_collision( const KDEWeightParam & param,
                                    const moab::CartVect & collision_loc, 
                                    int ebin )
{

  // make sure cross section parameter has been set properly
  if ( param.total_xs == NULL ) {

    std::cerr << "\nError: cross section parameter set to NULL" << std::endl;
    exit( EXIT_FAILURE );

  } 
 
  // make a KDECollision object to represent a single collision location
  KDECollision collision( collision_loc, bandwidth, kernel );

  // get valid neighborhood dimensions for non-zero contributions to tally
  double min[3];
  double max[3];

  collision.get_neighborhood( min, max );

  // find all the calculation points within the valid neighborhood
  std::vector<moab::EntityHandle> calc_points;
  moab::ErrorCode rval = points_in_box( *tree, tree_root, min, max, calc_points );
  assert( moab::MB_SUCCESS == rval );  

  // get the tally weighting factor for this collision
  double weight = get_score_weight( param );

  // compute the contribution for all calculation points in this neighborhood
  std::vector<moab::EntityHandle>::iterator i;
  moab::EntityHandle point;
  double contribution = 0;
  double coords[3];

  for ( i = calc_points.begin() ; i != calc_points.end() ; ++i ) {

    point = *i;
    rval = mb->get_coords( &point, 1, coords );
    assert( moab::MB_SUCCESS == rval );  

    contribution = weight * collision.compute_contribution( coords );
    
    // add results to the mesh tally for the current history
    add_score_to_tally( point, contribution, ebin );

  }
    
  // add collision to the running variance used in computing optimal bandwidth
  update_bandwidth_variance( collision_loc );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::tally_track( const KDEWeightParam & param,
                                const moab::CartVect & start_point,
                                const moab::CartVect & direction, 
                                int ebin )
{

  // make sure track length parameter has been set properly
  if ( param.tracklength == NULL ) {

    std::cerr << "\nError: track length parameter set to NULL" << std::endl;
    exit( EXIT_FAILURE );

  } 

  // set the number of subtracks to be tallied only for SUBTRACK tallies
  unsigned int tally_subtracks = 0;

  if ( kde_tally == SUBTRACK )
    tally_subtracks = subtracks;

  // make a KDETrack object to represent a single track segment
  KDETrack track( start_point, direction, bandwidth, *(param.tracklength),
                  kernel, tally_subtracks );
    
  // get valid neighborhood dimensions for non-zero contributions to tally
  double min[3];
  double max[3];
  double radius = 0;

  track.get_neighborhood( min, max, radius );

  // find all the calculation points within the valid neighborhood
  std::vector<moab::EntityHandle> calc_points;
  moab::ErrorCode rval = points_in_box( *tree, tree_root, min, max, calc_points );
  assert( moab::MB_SUCCESS == rval );  

  // get the tally weighting factor for this track
  double weight = get_score_weight( param );

  // compute the contribution for all calculation points in this neighborhood
  std::vector<moab::EntityHandle>::iterator i;
  moab::EntityHandle point;
  double contribution = 0;
  double coords[3];

  for ( i = calc_points.begin() ; i != calc_points.end() ; ++i ) {

    point = *i;
    rval = mb->get_coords( &point, 1, coords );
    assert( moab::MB_SUCCESS == rval );  

    contribution = weight * track.compute_contribution( coords );

    // add results to the mesh tally for the current history
    add_score_to_tally( point, contribution, ebin );
    
  }
    
  // TODO need to determine how to compute an optimal bandwidth with a track
  // This could be done using the subtracks/choose_points feature for both
  // SUBTRACK and TRACKLENGTH tallies.  For SUBTRACK tallies, simply would
  // need access to the subtrack_points vector in KDETrack.  For TRACKLENGTH
  // tallies, would need to use the choose_points function to get the points.
  // update_bandwidth_variance( collision_loc );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::end_history()
{
  
  std::set<moab::EntityHandle>::iterator i;
 
  // add results from current history to the tally for each calculation point
  for ( i = visited_this_history.begin() ; i != visited_this_history.end() ; ++i ) {
    
    for( unsigned int j = 0; j < ebins; ++j ){
      double& history = data_ref( temp_tally_data, *i, j );
      double& tally =   data_ref( tally_data, *i, j );
      double& error =   data_ref( error_data, *i, j );
      
      tally += history;
      error += ( history * history );
      
      // reset temp_tally for the next particle history
      history = 0;
      
    }
  }
  visited_this_history.clear();

}
//-----------------------------------------------------------------------------
void KDEMeshTally::print( double sp_norm, double fmesh_fact )
{

  // tags tally/error results to the nodes and writes mesh to output file
  write_results( sp_norm, fmesh_fact );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::write_results( double sp_norm, double fmesh_fact )
{

  double tally = 0;
  double error = 0, rel_err = 0;
  
  moab::ErrorCode rval = moab::MB_SUCCESS;

  // print the optimal bandwidth if it was computed
  if ( kde_tally == COLLISION ) {
  
    std::cout << std::endl << "optimal bandwidth for " << numCollisions;
    std::cout  << " collisions is: " << get_optimal_bandwidth() << std::endl;

  }

  // tag tally and relative error results to the mesh for each entity
  moab::Range::iterator i;
  
  for ( i = tally_ents.begin() ; i != tally_ents.end() ; ++i ) {

    moab::EntityHandle point = *i;

    for ( unsigned int j = 0; j < ebins; ++ j){

      tally = data_ref( tally_data, point, j);
      error = data_ref( error_data, point, j );
      
      // compute relative error for the tally
      // Use 0 as the rel_err value if nothing has been computed for this tally point;
      // this reflects MCNP's approach to avoiding a divide-by-zero situation.
      rel_err = 0; 
      if( error != 0 ){
        rel_err = sqrt( error / ( tally * tally ) - 1.0 / sp_norm );
      }
      
      // normalizing mesh tally results by the number of source particles
      tally /= sp_norm;
      
      // applying the fmesh multiplication FACTOR to the mesh tally results
      tally *= fmesh_fact;
      
      // set tally and error tag values for this entity
      rval = mb->tag_set_data( tally_tags[j], &point, 1, &tally );
      assert( moab::MB_SUCCESS == rval );
      
      rval = mb->tag_set_data( error_tags[j], &point, 1, &rel_err );
      assert( moab::MB_SUCCESS == rval ); 
      
    } 
    
  }

  // create a global tag to store the bandwidth value
  moab::Tag bandwidth_tag;
  rval = mb->tag_get_handle( "BANDWIDTH_TAG", 3, moab::MB_TYPE_DOUBLE, bandwidth_tag,
                             moab::MB_TAG_MESH|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == rval );

  // add bandwidth tag to the root set
  moab::EntityHandle bandwidth_set = mb->get_root_set();
  rval = mb->tag_set_data( bandwidth_tag, &bandwidth_set, 1, &bandwidth );
  assert( moab::MB_SUCCESS == rval );

  // define list of tags to include and write mesh to output file
  std::vector<moab::Tag> output_tags = tally_tags;
  output_tags.insert( output_tags.end(), error_tags.begin(), error_tags.end() );
  output_tags.push_back( bandwidth_tag );  

  rval = mb->write_file( output_filename.c_str(),
                         NULL, NULL, &tally_set, 1, &(output_tags[0]), output_tags.size() );
  assert( moab::MB_SUCCESS == rval );
  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::build_tree( moab::EntityHandle meshset )
{

  // Obtain all of the calculation points in the mesh and store into tally_ents
  moab::EntityType type = moab::MBVERTEX;
  
  moab::ErrorCode rval = mb->get_entities_by_type( meshset, type, tally_ents );
  assert( moab::MB_SUCCESS == rval );  

  resize_data_arrays( tally_ents.size() );
    
  // Measure the number of divisions in the moab::Range used to represent the tally points
  // If there are many divisions (for some rather arbitrary definition of "many"), print
  // a warning about performance compromise
  int psize = tally_ents.psize();
  std::cout << "   Tally range has psize: " << psize << std::endl;
  if( psize > 4 ){
    std::cerr << "Warning: large tally range psize " << psize 
              << ", may reduce performance." << std::endl;
  }
  // Build a KD-tree using all of the calculation points in the mesh
  moab::AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 10;
  settings.maxTreeDepth = 30;

  tree = new moab::AdaptiveKDTree( mb );
  rval = tree->build_tree( tally_ents, tree_root, &settings );
  assert( moab::MB_SUCCESS == rval );  

  rval = setup_tags( mb, "KDE_" );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::update_bandwidth_variance( const moab::CartVect & collision_loc )
{
 
  if ( numCollisions != LLONG_MAX ) {

    ++numCollisions;
    
    // obtain previous value for the mean
    moab::CartVect Mn_prev = Mn;
  
    // compute new values for the mean and variance
    if ( numCollisions == 1 )
      Mn = collision_loc;
    else {

      Mn += ( collision_loc - Mn_prev ) / numCollisions;
    
      for ( int i = 0 ; i < 3 ; ++i )
        Sn[i] += ( collision_loc[i] - Mn_prev[i] ) * ( collision_loc[i] - Mn[i] );
    
    }

  }
  else if ( !max_collisions ) {
  
    std::cerr << "/nWarning: number of collisions exceeds maximum/n"
              << "  optimal bandwidth will be based on " << numCollisions
              << " collisions./n/n";

    max_collisions = true;

  }

}
//-----------------------------------------------------------------------------
moab::CartVect KDEMeshTally::get_optimal_bandwidth()
{
  
  double stdev = 0;
  moab::CartVect optimal_bandwidth;
  
  for ( int i = 0 ; i < 3 ; ++i ) {

    stdev = sqrt( Sn[i] / ( numCollisions - 1 ) );
    optimal_bandwidth[i] = 0.968625 * stdev * pow( numCollisions, -1.0/7.0 );

  }

  return optimal_bandwidth;

}
//-----------------------------------------------------------------------------
double KDEMeshTally::get_score_weight( const KDEWeightParam & param )
{

  double score_weight = 1;
  double d = 1;  // track length set to 1 so it does not change weight
  
  if ( kde_tally == SUBTRACK )
    d = *(param.tracklength);  // use actual track length for SUBTRACK tallies
  
  // determine the tally weighting factor from MCNP
  mcnp_weight_calculation( param.fmesh_index, param.energy, param.particle_wgt,
                           &d, &score_weight );
  
  // divide weight by total cross section for COLLISION tallies only
  if ( kde_tally == COLLISION )
    score_weight /= *(param.total_xs);

  return score_weight;

}
//-----------------------------------------------------------------------------
void KDEMeshTally::add_score_to_tally( moab::EntityHandle mesh_point,
                                       double score,
                                       int ebin )
{

  data_ref( temp_tally_data, mesh_point, ebin ) += score;

  // tally the total energy bin if requested
  if ( fmesh.total_energy_bin )
    data_ref( temp_tally_data, mesh_point, (ebins-1) ) += score;

  visited_this_history.insert( mesh_point );

}
//-----------------------------------------------------------------------------
moab::ErrorCode
  KDEMeshTally::points_in_box( moab::AdaptiveKDTree & tree,
                               moab::EntityHandle tree_root,
                               const double box_min[3],
                               const double box_max[3],
                               std::vector<moab::EntityHandle> & points )
{
  
  // determine the center point of the box
  double box_center[3];
  
  for ( int i = 0 ; i < 3 ; ++i )
    box_center[i] = 0.5 * ( box_max[i] + box_min[i] );  
  
  // set radius equal to distance from the center to the max corner of box
  moab::CartVect max_corner( box_max );
  moab::CartVect center( box_center );
  max_corner -= center;  
  double radius = max_corner.length();
  
  // find all leaves of the tree within the given radius
  std::vector<moab::EntityHandle> leaves;
  moab::ErrorCode rval =
    tree.leaves_within_distance( tree_root, box_center, radius, leaves );
  
  if ( moab::MB_SUCCESS != rval )
    return rval;

  // obtain the set of unique points in the box 
  std::vector<moab::EntityHandle>::iterator i; 
  moab::Interface* mb = tree.moab();
  moab::Range leaf_points;
  moab::Range::iterator j;
  moab::EntityHandle point;
  double coords[3];
  
  // iterate through the leaves
  for ( i = leaves.begin() ; i != leaves.end() ; ++i) {

    leaf_points.clear();
    rval = mb->get_entities_by_type( *i, moab::MBVERTEX, leaf_points );
    
    if ( moab::MB_SUCCESS != rval )
      return rval;
  
    // iterate through the points in each leaf  
    for ( j = leaf_points.begin() ; j != leaf_points.end() ; ++j ) {

      point = *j;
      rval = mb->get_coords( &point, 1, coords );
    
      if ( moab::MB_SUCCESS != rval )
        return rval;

      // check point is in the box
      bool in_box = true;
      int k = 0;

      do {

        if ( coords[k] < box_min[k] || coords[k] > box_max[k] )
          in_box = false;

        ++k;

      }
      while ( true && k < 3 );
      
      // add the point to the set if it is in the box
      if ( in_box )
        points.push_back( point );

    }

  }
  
  // remove duplicates from points vector
  std::sort( points.begin(), points.end() );
  points.erase( std::unique( points.begin(), points.end() ), points.end() );
  
  return moab::MB_SUCCESS;

}
//-----------------------------------------------------------------------------
