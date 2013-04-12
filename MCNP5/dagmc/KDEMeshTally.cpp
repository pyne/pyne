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
#include "KDENeighborhood.hpp"
#include "KDETrack.hpp"
#include "KDEMeshTally.hpp"
#include "TallyEvent.hpp"
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
void KDEMeshTally::compute_score(const TallyEvent& event, int ebin)
{
  // make sure tally event has been set
  TallyEvent::EventType type = event.get_event_type();

  if (type == TallyEvent::NONE)
  {
    std::cerr << "\nError: Tally event type has not been set" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (type == TallyEvent::TRACK)
  {
    tally_track(event, ebin);
  }
  else if (type == TallyEvent::COLLISION)
  {
    tally_collision(event, ebin);
  }
}
//-----------------------------------------------------------------------------
void KDEMeshTally::tally_collision(const TallyEvent& event, int ebin)
{
  // make sure tally event is a collision event
  if (event.get_event_type() != TallyEvent::COLLISION)
  {
    std::cerr << "\nError: Tally event is not a collision event" << std::endl;
    exit(EXIT_FAILURE);
  } 
 
  // make a KDECollision object to represent a single collision location
  CollisionData data;
  bool set_data = event.get_collision_data(data);

  if (!set_data)
  {
    std::cerr << "\nError: Invalid collision event data" << std::endl;
    exit(EXIT_FAILURE);
  }

  KDECollision collision( data.collision_point, bandwidth, kernel );

  // create the neighborhood region for this tally event
  KDENeighborhood region(event, bandwidth, *tree, tree_root);

  // find all the calculation points within this neighborhood region
  std::vector<moab::EntityHandle> calculation_points;
  moab::ErrorCode rval = region.get_points(calculation_points);
  assert(moab::MB_SUCCESS == rval);

  // get the tally weighting factor for this collision
  double weight = event.get_tally_multiplier() * data.particle_weight;

  // divide by the total cross section
  weight /= data.total_cross_section;

  // compute the contribution for all calculation points in this neighborhood
  std::vector<moab::EntityHandle>::iterator i;
  moab::EntityHandle point;
  double contribution = 0;
  double coords[3];

  for (i = calculation_points.begin(); i != calculation_points.end(); ++i)
  {

    point = *i;
    rval = mb->get_coords( &point, 1, coords );
    assert( moab::MB_SUCCESS == rval );  

    contribution = weight * collision.compute_contribution( coords );
    
    // add results to the mesh tally for the current history
    add_score_to_tally( point, contribution, ebin );

  }
    
  // add collision to the running variance used in computing optimal bandwidth
  update_bandwidth_variance(data.collision_point);

}
//-----------------------------------------------------------------------------
void KDEMeshTally::tally_track(const TallyEvent& event, int ebin)
{
  // make sure tally event is a track-based event
  if (event.get_event_type() != TallyEvent::TRACK)
  {
    std::cerr << "\nError: Tally event is not a track-based event" << std::endl;
    exit(EXIT_FAILURE);
  } 

  // set the number of subtracks to be tallied only for SUBTRACK tallies
  unsigned int tally_subtracks = 0;

  if ( kde_tally == SUBTRACK )
    tally_subtracks = subtracks;

  // make a KDETrack object to represent a single track segment
  TrackData data;
  bool set_data = event.get_track_data(data);

  if (!set_data)
  {
    std::cerr << "\nError: Invalid track-based event data" << std::endl;
    exit(EXIT_FAILURE);
  }

  KDETrack track(data.start_point, data.direction, bandwidth, data.track_length,
                 kernel, tally_subtracks);

  // create the neighborhood region for this tally event
  KDENeighborhood region(event, bandwidth, *tree, tree_root);

  // find all the calculation points that exist in this neighborhood region
  std::vector<moab::EntityHandle> calculation_points;
  moab::ErrorCode rval = region.get_points(calculation_points);
  assert(moab::MB_SUCCESS == rval);

  // get the tally weighting factor for this track
  double weight = event.get_tally_multiplier() * data.particle_weight;

  // if SUBTRACK mesh tally, multiply weight by the total track length
  if (kde_tally == SUBTRACK)
    weight *= data.track_length;

  // compute the contribution for all calculation points in this neighborhood
  std::vector<moab::EntityHandle>::iterator i;
  moab::EntityHandle point;
  double contribution = 0;
  double coords[3];

  for (i = calculation_points.begin(); i != calculation_points.end(); ++i)
  {

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
    
    for( unsigned int j = 0; j < num_energy_bins; ++j ){
      double& history = get_data( temp_tally_data, *i, j );
      double& tally =   get_data( tally_data, *i, j );
      double& error =   get_data( error_data, *i, j );
      
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
  
  for ( i = tally_points.begin() ; i != tally_points.end() ; ++i ) {

    moab::EntityHandle point = *i;

    for ( unsigned int j = 0; j < num_energy_bins; ++ j){

      tally = get_data( tally_data, point, j);
      error = get_data( error_data, point, j );
      
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

  // Obtain all of the calculation points in the mesh and store into tally_points
  moab::EntityType type = moab::MBVERTEX;
  
  moab::ErrorCode rval = mb->get_entities_by_type( meshset, type, tally_points );
  assert( moab::MB_SUCCESS == rval );  

  resize_data_arrays( tally_points.size() );
    
  // Measure the number of divisions in the moab::Range used to represent the tally points
  // If there are many divisions (for some rather arbitrary definition of "many"), print
  // a warning about performance compromise
  int psize = tally_points.psize();
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
  rval = tree->build_tree( tally_points, tree_root, &settings );
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
void KDEMeshTally::add_score_to_tally( moab::EntityHandle mesh_point,
                                       double score,
                                       int ebin )
{

  get_data( temp_tally_data, mesh_point, ebin ) += score;

  // tally the total energy bin if requested
  if ( input_data.total_energy_bin )
    get_data( temp_tally_data, mesh_point, (num_energy_bins-1) ) += score;

  visited_this_history.insert( mesh_point );

}
//-----------------------------------------------------------------------------
