#include <iostream>
#include <sstream>
#include <cmath>
#include <set>

#include "moab/Range.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/Skinner.hpp"
#include "moab/CN.hpp"
  

/* Two macros are available: 
 * MESHTAL_DEBUG: produce much debugging output, with histories of particular particle tracks.
 *
 * MESHTAL_FORCE_ASSERTS: force use of assert(), regardless of compile settings.  Useful for 
 *     testing at scale; asserts do not lead to substantial slowdown at run time.
 *     MESHTAL_DEBUG implies MESHTAL_FORCE_ASSERTS.
 */

#ifdef MESHTAL_DEBUG
#define MESHTAL_FORCE_ASSERTS
#endif 

#ifdef MESHTAL_FORCE_ASSERTS
#undef NDEBUG 
#endif

#include <cassert>

// the header file has at least one assert, so keep this include below the macro checks
#include "TallyEvent.hpp"
#include "TrackLengthMeshTally.hpp"
#include "meshtal_funcs.h"


// tolerance for ray-triangle intersection tests
// (note: this paramater is ignored by GeomUtil, so don't bother trying to tune it)
#define TRIANGLE_INTERSECTION_TOL 1e-6

// If the following is defined, use OBB trees for ray-triangle intersections,
// otherwise use KD tree
//#define USE_OBB_TREE_RAY_TRACING

/* Tetrahedron volume code taken from MOAB/tools/measure.cpp */
inline static double tet_volume( const moab::CartVect& v0,
                                 const moab::CartVect& v1,
                                 const moab::CartVect& v2,
                                 const moab::CartVect& v3 )
{
  return 1./6. * ( ((v1 - v0) * (v2 - v0)) % (v3 - v0) );
}



// Adapted from MOAB's convert.cpp
// Parse list of integer ranges, e.g. "1,2,5-10,12"
static
bool parse_int_list( const char* string, std::set<int>& results )
{
  bool okay = true;
  char* mystr = strdup( string );
  for (const char* ptr = strtok(mystr, ", \t"); ptr; ptr = strtok(0,", \t"))
  {
    char* endptr;
    long val = strtol( ptr, &endptr, 0 );
    if (endptr == ptr) {
      std::cerr << "Not an integer: \"" << ptr << '"' << std::endl;
      okay = false;
      break;
    }
    
    long val2 = val;
    if (*endptr == '-') {
      const char* sptr = endptr+1;
      val2 = strtol( sptr, &endptr, 0 );
      if (endptr == sptr) {
        std::cerr << "Not an integer: \"" << sptr << '"' << std::endl;
        okay = false;
        break;
      }
      if (val2 < val) {
        std::cerr << "Invalid id range: \"" << ptr << '"' << std::endl;
        okay = false;
        break;
      }
    }
    
    if (*endptr) {
      okay = false;
      break;
    }
    
    for (; val <= val2; ++val)
      results.insert( (int)val );

  }
  
  free( mystr );
  return okay;    
}

extern "C"{
extern int namchg_( int*, int* );
}

static 
bool map_conformal_names( std::set<int>& input, std::set<int>& output ){
  
  for( std::set<int>::iterator i = input.begin(); i!=input.end(); ++i){
    int x, y, one = 1;
    x = *i;
    y = namchg_( &one, &x );
#ifdef MESHTAL_DEBUG
    std::cerr << "namchg mapped cell " << *i << " to name " << y << std::endl;
#endif
    if( y == 0 ){
        std::cerr << " conformality cell " << *i << " does not exist." << std::endl;
        return false;
    }
    output.insert( y );
  }
  return true;
}

namespace moab { 



TrackLengthMeshTally* TrackLengthMeshTally::setup( const MeshTallyInput& fmesh_params, 
                                                   Interface* mbi, const int* current_mcnp_cell )
{

  const MeshTallyInput::TallyOptions& fc_params = fmesh_params.options;

  std::string input_filename, output_filename, tag_name;
  std::vector<std::string> tag_values;

  bool convex_flag = false, conformal_flag = false;
  std::set<int> conf_tmp, *conformality = NULL;
  bool conformal_surface_source = false;

  int id = fmesh_params.tally_id;

  for( MeshTallyInput::TallyOptions::const_iterator i = fc_params.begin();
       i != fc_params.end(); ++i )
  {
    std::string key = (*i).first, val = (*i).second;
    if( key == "inp" ) input_filename = val;
    else if( key == "out" ) output_filename = val;
    else if( key == "tag" ) tag_name = val;
    else if( key == "tagval" ) tag_values.push_back(val);
    else if( key == "convex" && (val == "t" || val == "true" ) ) convex_flag = true; 
    else if( key == "conf_surf_src" && (val == "t" || val == "true" ) ) conformal_surface_source = true;
    else if( key == "conformal" ) { 

      if( !conformality ) conformality = new std::set<int>();
      if( !parse_int_list( val.c_str(), *conformality ) ){
        std::cerr << "Error: FC" << id << " card has bad conformality value '" << val << "'" << std::endl;
        exit(EXIT_FAILURE);
      }
      if( !map_conformal_names( *conformality, conf_tmp ) ){
        std::cerr << "Error: a conformal cell does not exist in the problem!" << std::endl;
        exit(EXIT_FAILURE);
      }

    }
    else{
      std::cerr << "Warning: FC" << id << " card has unknown key '" << key << "'" << std::endl;
    }
  }
  
  std::stringstream strbuf;
  if( input_filename.length() == 0 ){
    strbuf << "fmesh" << id << ".h5m";
    strbuf >> input_filename;
    strbuf.clear();
  }
  if( output_filename.length() == 0 ){
    strbuf << "meshtal" << id << ".h5m";
    strbuf >> output_filename;
    strbuf.clear();
  }

  std::cout << "Creating dagmc fmesh" << id 
            << ", input: " << input_filename 
            << ", output: " << output_filename << std::endl;

  if( tag_name != "" ){
    std::cout << "  using tag name='" << tag_name << "'";
    if( tag_values.size() > 0 ){
      std::cout <<", and tag values= " << std::endl;
      for( unsigned int i = 0; i < tag_values.size(); ++i ){
        std::cout << "    '" << tag_values[i] << "'" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  if( convex_flag ){ 
    std::cout << "  user asserts that this tally mesh has convex geometry." << std::endl;
  }
  
  if( conformality ){
    std::cout << "  conformal to cells " << std::flush;
    for( std::set<int>::iterator i = conformality->begin(); i!=conformality->end(); ){
      std::cout << *i; 
      if( ++i != conformality->end() ) std::cout << ", ";
    }
    std::cout << std::endl; 
    *conformality = conf_tmp;
  }

  moab::TrackLengthMeshTally *mt = new moab::TrackLengthMeshTally( fmesh_params, mbi, output_filename );
  
  if( convex_flag && conformal_flag ){
    std::cerr << "Warning: FC" << id << " specifies both conformal and convex logic; using conformal logic." << std::endl;
  }
  if( conformality ){
    mt->mcnp_current_cell = current_mcnp_cell;
    mt->conformality = conformality;
    mt->conformal_surface_source = conformal_surface_source;
  }
  else mt->set_convex_flag( convex_flag );

  moab::ErrorCode rval;
  rval = mt->load_mesh( input_filename, tag_name, tag_values );
  if( rval != moab::MB_SUCCESS ){
    std::cerr << "** DAGMC TrackLengthMeshTally creation failed!" << std::endl;
    std::cerr << "** Tally " << id << " failed to initialize." << std::endl;
    exit( EXIT_FAILURE );
    
  }
  
  return mt;

}

  
TrackLengthMeshTally::TrackLengthMeshTally( const MeshTallyInput& fmesh, Interface* mb_p, 
                                            const std::string& output_filename_p ) :
  MeshTally( fmesh ),
  mb( mb_p ), output_filename(output_filename_p),  
  obb_tool( new OrientedBoxTreeTool(mb_p) ),
  last_visited_tet( 0 ), 
  convex( false ), conformality( NULL ), conformal_surface_source( false ),
  mcnp_current_cell( NULL ), last_cell( -1 ), num_negative_tracks(0)
{}
  
TrackLengthMeshTally::~TrackLengthMeshTally(){
  delete obb_tool;
}

/**
 * Load the given file as an input mesh
 * Member variable tally_set will contain the mesh contents
 */
ErrorCode TrackLengthMeshTally::load_mesh( const std::string& input_filename, 
                                           std::string tag_name, std::vector<std::string>& tag_values ){
  
  ErrorCode rval;

  EntityHandle loaded_file_set;
  rval = mb->create_meshset( MESHSET_SET, loaded_file_set );
  assert( rval == MB_SUCCESS );
  rval = mb->load_file( input_filename.c_str(), &loaded_file_set );
  if( rval != MB_SUCCESS ){
    std::cerr << "Warning: could not load tally mesh " << input_filename << std::endl;
    return rval;
  }


  rval = mb->create_meshset( MESHSET_SET, tally_set );
  assert( rval == MB_SUCCESS );
  
 

  if( tag_name.length() > 0 ){

    std::cout << "  User-specified tag to load:  " << tag_name << std::endl;

    /* Until there is more certainty about the type and parameters of the tag the user specified,
       use MB_TAG_ANY to get access to any tag with the given name */
    Tag user_spec_tag;
    rval = mb->tag_get_handle( tag_name.c_str(), 0, MB_TYPE_OPAQUE, user_spec_tag, MB_TAG_ANY );
    assert( rval == MB_SUCCESS );
    
    int user_spec_tag_length = 0;
    rval = mb->tag_get_bytes( user_spec_tag, user_spec_tag_length );
    assert( rval == MB_SUCCESS );

    std::cout << "  user tag length: " << user_spec_tag_length << " bytes" << std::endl;

    Range user_sets;
    rval = mb->get_entities_by_type_and_tag( loaded_file_set, MBENTITYSET, &user_spec_tag, NULL, 1, user_sets );
    assert( rval == MB_SUCCESS );

    std::cout << "  Found " << user_sets.size() << " sets with this tag." << std::endl;

    for( Range::iterator i = user_sets.begin(); i!=user_sets.end(); ++i){
      EntityHandle s = *i;
      char* name = new char[ user_spec_tag_length + 1];
      
      rval = mb->tag_get_data( user_spec_tag, &s, 1, name );
      assert( rval == MB_SUCCESS );

      // if user specified no tag value, list the available ones for informational purposes
      if( tag_values.size() == 0 ){
        std::cout << "    available tag value: " << name << std::endl; 
      }
      
      if( std::find( tag_values.begin(), tag_values.end(),std::string(name) ) != tag_values.end() ){
        std::cout << "  Successfully found a set with tag value " << name << std::endl;
        rval = mb->unite_meshset( tally_set, s );
        assert( rval == MB_SUCCESS );
      }
      delete[] name;
    }

  }
  else{ // no user-specified tag filter

    rval = mb->unite_meshset( tally_set, loaded_file_set );
    assert( rval == MB_SUCCESS );
  }

  rval = setup_tags( mb );
  assert( rval == MB_SUCCESS );

  int num_tets;
  rval = mb->get_number_entities_by_dimension( tally_set, 3, num_tets );
  assert( rval == MB_SUCCESS );
  std::cerr << "  There are " << num_tets << " tetrahedrons in this tally mesh." << std::endl;


  tet_baryc_data.resize( num_tets );  
  resize_data_arrays( num_tets );

  Range all_tets;
  rval = mb->get_entities_by_dimension( tally_set, 3, all_tets );
  assert( rval == MB_SUCCESS );
  assert( all_tets.size() == (unsigned)num_tets );

  // restruct tally_set to contain only 3-dimensional elements
  rval = mb->clear_meshset( &tally_set, 1  );
  assert( rval == MB_SUCCESS );
  rval = mb->add_entities( tally_set, all_tets );
  assert( rval == MB_SUCCESS );

  tally_points = all_tets;

  // Measure the number of divisions in the moab::Range used to represent the tally tets
  // If there are many divisions (for some rather arbitrary definition of "many"), print
  // a warning about performance compromise
  int psize = tally_points.psize();
  std::cout << "  Tally range has psize: " << psize << std::endl;
  if( psize > 4 ){
    std::cerr << "Warning: large tally range psize " << psize << ", may reduce performance." << std::endl;
  }

  /**
   * Iterate over all tets and compute barycentric matrices 
   */
  for( Range::const_iterator i=all_tets.begin(); i!=all_tets.end(); ++i){
    EntityHandle tet = *i;

    const EntityHandle* verts;
    int num_verts;
    rval = mb->get_connectivity( tet, verts, num_verts );
    assert( rval == MB_SUCCESS );
    
    if( num_verts != 4 ){
      std::cerr << "Error: DAGMC TrackLengthMeshTally cannot handle non-tetrahedral meshes yet," << std::endl;
      std::cerr << "       but your mesh has at least one cell with " << num_verts << " vertices." << std::endl;
      return MB_NOT_IMPLEMENTED;
    }
    
    CartVect p[4];
    rval = mb->get_coords( verts, 4, p[0].array() );
    assert( rval == MB_SUCCESS );

    Matrix3 a( p[1]-p[0], p[2]-p[0], p[3]-p[0] );
    a = a.transpose().inverse();
    tet_baryc_data.at( get_entity_index(tet) ) = a;

  }

  // prepare to build KD tree and OBB tree
  Range all_tris;
  Skinner skinner(mb);
  skinner.find_skin( all_tets, 2, all_tris );
  std::cout << "  Tally mesh skin has " << all_tris.size() << " triangles." << std::endl;


#ifdef USE_OBB_TREE_RAY_TRACING
  std::cout << " Building OBB tree of size " << all_tris.size() << "... " << std::flush;
  
  rval = obb_tool->build( all_tris, obbtree_root );
  assert( rval == MB_SUCCESS );
  std::cout << "done." << std::endl;

#else

  // put tris with tets to be rolled into KD tree
  all_tets.merge( all_tris );

#endif

  // build KD tree of all tetrahedra
  std::cout << "  Building KD tree of size " << all_tets.size() << "... " << std::flush;
  kdtree = new AdaptiveKDTree( mb );
  AdaptiveKDTree::Settings kdtree_settings;
  //kdtree_settings.candidatePlaneSet = MBAdaptiveKDTree::VERTEX_SAMPLE;
  kdtree->build_tree( all_tets, kdtree_root, &kdtree_settings );
  std::cout << "done." << std::endl << std::endl;;

  return MB_SUCCESS;

}

  void TrackLengthMeshTally::print( double sp_norm, double mult_fact ){
    write_results( sp_norm, mult_fact );
  }

/**
 * Write out the mesh with tally and error tags attached
 * @param sp_norm Tally for each tetrahedron will be divided by sp_norm*volume of tet.
 * @param mult_fact Tally for each tetrahedron will be multiplied by this parameter.
 * @param override_output_filename If non-NULL, use this filename for output instead
 *                                 of this->output_filename
 */
ErrorCode TrackLengthMeshTally::write_results( double sp_norm, double mult_fact, 
                                    const std::string* override_output_filename )
{

  std::string filename = override_output_filename ? *override_output_filename : output_filename;

  ErrorCode rval;

  Range all_tets;
  rval = mb->get_entities_by_dimension( tally_set, 3, all_tets );
  assert( rval == MB_SUCCESS );


  for( Range::const_iterator i=all_tets.begin(); i!=all_tets.end(); ++i){
    EntityHandle t = *i;

    CartVect v[4];

    std::vector<EntityHandle> vtx;    
    mb->get_connectivity( &t, 1, vtx );
    assert( vtx.size() == 4);

    int k = 0;
    for( std::vector<EntityHandle>::iterator j = vtx.begin(); j!=vtx.end(); ++j){
      EntityHandle vertex = *j;
      mb->get_coords( &vertex, 1, v[k++].array() );
    }

    double volume = tet_volume( v[0], v[1], v[2], v[3] );

    for( unsigned j = 0; j < num_energy_bins; ++j ){

      double tally = get_data( tally_data, t, j );
      double error = get_data( error_data, t, j );
      double score = (tally / (volume*sp_norm)) * mult_fact;
      
      rval = mb->tag_set_data( tally_tags[j], &t, 1, &score );
      assert( rval == MB_SUCCESS );

      // Use 0 as the error output value if nothing has been computed for this mesh cell;
      // this reflects MCNP's approach to avoiding a divide-by-zero situation.
      double rel_err = 0;
      if( error != 0 ){
        rel_err = sqrt( (error / (tally*tally)) - (1./sp_norm) );
      }        

      rval = mb->tag_set_data( error_tags[j], &t, 1, &rel_err );
      assert( rval == MB_SUCCESS );

    }
  }

  std::vector<Tag> output_tags = tally_tags;
  output_tags.insert( output_tags.end(), error_tags.begin(), error_tags.end() );

  rval = mb->write_file( filename.c_str(), NULL, NULL, &tally_set, 1, &(output_tags[0]), output_tags.size() );
  assert (rval == MB_SUCCESS );
 
  if ( num_negative_tracks != 0 )
  {
    std::cout << std::endl;
    std::cout << num_negative_tracks << " negative tracks occurred during the simulation." << std::endl;
    std::cout << "These tracks were not included in the final tally results." << std::endl << std::endl;
  }

  return MB_SUCCESS;
}


/**
 * Return true if the point falls inside tet.  Assumes tet is part of this TrackLengthMeshTally.
 */
bool TrackLengthMeshTally::point_in_tet( const CartVect& point, const EntityHandle* tet )
{ 
  ErrorCode rval;
                    
  const EntityHandle* verts;
  int num_verts;
  rval = mb->get_connectivity( *tet, verts, num_verts );
  assert( rval == MB_SUCCESS );
  
  CartVect p0;
  rval = mb->get_coords( verts, 1, p0.array() );
  assert( rval == MB_SUCCESS );

  Matrix3& Ainverse = tet_baryc_data[ get_entity_index(*tet) ];

  CartVect bary = (Ainverse) * (point-p0);
  
  bool in_tet = ( bary[0]>= 0 && bary[1] >= 0 && bary[2] >= 0 &&
                  bary[0]+bary[1]+bary[2] <= 1. );

  return in_tet;
}


/**
 * Add a score to a given mesh cell
 */
void TrackLengthMeshTally::add_score_to_mesh_cell( EntityHandle mesh_cell, double score, int ebin ){
  
  visited_this_history.insert( mesh_cell );

  get_data( temp_tally_data, mesh_cell, ebin ) += score;

  if( input_data.total_energy_bin ){
    get_data( temp_tally_data, mesh_cell, (num_energy_bins-1) ) += score;
  }
}

/**
 * Finish adding a set of scores for a particular monte carlo particle track
 */
void TrackLengthMeshTally::end_history () {

  for( std::set< EntityHandle >::iterator i=visited_this_history.begin(); i!=visited_this_history.end(); ++i){

    for( unsigned j = 0; j < num_energy_bins; ++j ){
      double& d =     get_data(temp_tally_data, *i, j );
      double& tally = get_data(tally_data, *i, j );
      double& error = get_data(error_data, *i, j );
      
      tally += d;
      error += (d * d);
      d = 0.0;
    }

  }

  visited_this_history.clear();

  if( conformality ){ last_cell = -1; } 
}

/**
 * Given an EntityHandle representing a mesh skin triangle, return the tetrahedron it belongs to
 * as well as the handles of its vertices.  Calling this function asserts that the triangle 
 * "belongs" to exaclty one tetrahedron.  This should always be a valid assumption, since this
 * code only creates explicit triangle entities for skin triangles.
 * 
 * @param triangle The input triangle
 * @param tetrahedron (output) The parent tetrahedron of triangle
 * @param vertices (output) The vertices of triangle
 *
 */
void
TrackLengthMeshTally::get_skin_triangle_adjacencies( EntityHandle triangle, 
                                                     EntityHandle& tetrahedron, EntityHandle vertices[3] ){

  ErrorCode rval;
  const EntityHandle* tri_conn;
  int num_verts; 
  rval = mb->get_connectivity( triangle, tri_conn, num_verts );
  assert( rval == MB_SUCCESS );
  assert( num_verts == 3 );
  memcpy( vertices, tri_conn, sizeof(EntityHandle)*3 );
  
  Range tri_sides; 
  rval = mb->get_adjacencies( &triangle, 1, 3, false, tri_sides );
  assert( rval == MB_SUCCESS );
  assert( tri_sides.size() == 1 );
  tetrahedron = tri_sides[0];

}

/**
 * Given a ray that starts outside the mesh, find the next tetrahedron (if any) 
 * that will be intersected by the ray.
 * @param start Beginning of ray
 * @param vec Unit vector of ray
 * @param length Length of ray
 * @param first_tri If return value is non-zero, will contain the skin triangle that 
 *                  ray intersects.
 * @param first_t Value of t along the ray where first_tri is intersected.
 * @param last_crossing A skin triangle to ignore for this ray fire.  Used if a ray is
 *                      being fired away from the surface of the mesh.
 * @return The tetrahedron intersected by this ray fire, or zero if none.
 */
EntityHandle 
TrackLengthMeshTally::find_next_tet_by_ray_fire(const CartVect& start, const CartVect& vec, double length, 
                                      EntityHandle first_tri[3], double& first_t, 
                                      EntityHandle last_crossing)
{
  ErrorCode rval;
  EntityHandle first_tet = 0;

  std::vector< EntityHandle > triangles;
  std::vector< double > distances;

#ifdef USE_OBB_TREE_RAY_TRACING
  double vlen = length;
  rval = obb_tool->ray_intersect_triangles( distances, triangles,
                                            obbtree_root, TRIANGLE_INTERSECTION_TOL, 
                                            start.array(), vec.array(), &vlen );
#else
  rval = kdtree->ray_intersect_triangles( kdtree_root, TRIANGLE_INTERSECTION_TOL, 
                                          vec.array(), start.array(),
                                          triangles, distances, 0, length );  
#endif

  assert( rval == MB_SUCCESS );
  
  if( distances.size() ){
    // ray goes into mesh
    int min_idx = 0;
    for( unsigned i = 1; i < distances.size(); ++i ){
      if( distances[i] < distances[min_idx] &&
          triangles[i] != last_crossing )
        { min_idx = i; }  
    }
    
    if( distances[min_idx] < length && triangles[min_idx] != last_crossing ){ 

#if MESHTAL_DEBUG
      std::cout << "    ray fire success, d=" << distances[min_idx] << std::endl;
#endif

      EntityHandle first_tri_eh = triangles[min_idx];
      first_t   = distances[min_idx];

      get_skin_triangle_adjacencies( first_tri_eh, first_tet, first_tri );
      
    }
  }

  return first_tet;
}

/**
 * Find the first tetrahedron intersected by a track segment, using conformal logic.
 * This is invoked when a surface crossing has just occurred in a conformal mesh tally,
 * allowing us to assume that one of the surface triangles of this mesh has just been crossed.
 * Thus the search for the nearest tetrahedron is simply a search for the surface triangle
 * nearest to the starting point of the ray.
 *
 * This function should always successfully return values, but it will print a warning if 
 * the nearest available result is suspicious.
 *
 * @param start Beginning of track 
 * @param first_tri (output) The vertices of the skin triangle at which track begins
 * @return The first tetrahedron along the ray.
 */
EntityHandle
TrackLengthMeshTally::get_starting_tet_conformal(const CartVect& start, EntityHandle first_tri[3])
{
  ErrorCode rval;
  CartVect tri_start;
  EntityHandle first_tri_eh, first_tet=0;
  
  rval = kdtree->closest_triangle( kdtree_root, start.array(), tri_start.array(), first_tri_eh );
  assert( rval == MB_SUCCESS );
  
  get_skin_triangle_adjacencies( first_tri_eh, first_tet, first_tri );
  
  double distance = (start-tri_start).length();

#ifdef MESHTAL_DEBUG
  std::cout << "  start: conformal mesh entry, dist =" << distance << std::endl;    
#endif

  if( distance > 1e-6 ){
    std::cerr << "Warning: suspicious epsilon distance in conformal tally: " << distance << std::endl;
  }
  
  return first_tet;

}

/**
 * Find the first tetrahedron intersected by the given track segment.
 * @param start Beginning of ray
 * @param vec Unit vector of ray
 * @param length Length of ray
 * @param first_tri If return value is non-zero, and a ray fire was performed, 
 *                  will contain the skin triangle that ray intersects.
 * @param first_t Value of t along the ray where first_tri is intersected.
 * @param conformal_begin_track If true, track is known to enter mesh at a conformal 
 *                              MCNP cell at t=0.  Use simplified logic in this case.
 * @return The first tetrahedron along the ray, or zero if none found.
 */

EntityHandle 
TrackLengthMeshTally::get_starting_tet(const CartVect& start, const CartVect& vec, double length, 
                                       EntityHandle first_tri[3], double& first_t )
{
  ErrorCode rval;
  EntityHandle first_tet = 0; 

  
  // First check last visited tetrahedron
  if( last_visited_tet && this->point_in_tet( start, &last_visited_tet ) ){
    first_tet = last_visited_tet; 

#ifdef MESHTAL_DEBUG
    std::cout << "   start: cached first tet" << std::endl;
#endif

  }
  else{ 
    // Check to see if starting point begins inside a tet
    AdaptiveKDTreeIter tree_iter;
    rval = kdtree->leaf_containing_point( kdtree_root, start.array(), tree_iter );
    if( rval == MB_SUCCESS ){
      
      EntityHandle leaf = tree_iter.handle();
      Range candidate_tets;
      rval = mb->get_entities_by_dimension( leaf, 3, candidate_tets, false );
      assert( rval == MB_SUCCESS );
      
      for( Range::const_iterator i = candidate_tets.begin(); i!=candidate_tets.end(); ++i){
        if( this->point_in_tet( start, &(*i) )){
          first_tet = *i;
        }
      }
    }
#ifdef MESHTAL_DEBUG
    if( first_tet ){
      std::cout << "   start: pt in volume" << std::endl;
    }
#endif 

  }
  
  if( first_tet == 0 ){
    // Starting point is outside of the mesh: fire the ray to see if it enters mesh.
    first_tet = this->find_next_tet_by_ray_fire( start, vec, length, first_tri, first_t );

#ifdef MESHTAL_DEBUG
    if( first_tet ){
      std::cout << "   start: ray fire"  << std::endl;
    }
    if( conformality ){
      std::cout << "    WARNING: conformality enabled, but a start was selected by ray fire!" << std::endl;
    }
#endif

  }
  
  return first_tet;
}

static inline bool tris_eq( const EntityHandle *t1, const EntityHandle *t2 ){
  int ignored1, ignored2;
  return CN::ConnectivityMatch( t1, t2, 3, ignored1, ignored2 );
}

void TrackLengthMeshTally::compute_score(const TallyEvent& event, int ebin)
{
  // make sure tally event is a track-based event
  if (event.get_event_type() != TallyEvent::TRACK)
  {
    std::cerr << "\nError: Tally event is not a track-based event" << std::endl;
    exit(EXIT_FAILURE);
  } 

  ErrorCode rval;

  bool conformal_begin_track = false;

  if( conformality ){

    bool cell_change = (last_cell != *mcnp_current_cell);
    bool newparticle = (last_cell == -1);
    
    conformal_begin_track = cell_change;
    // new particles should only use conformal crossing logic if a conformal surface source was declared
    if( newparticle )
        conformal_begin_track = conformal_surface_source; 

#ifdef MESHTAL_DEBUG
    if( newparticle ){ std::cout << "Started new particle in cell " << *mcnp_current_cell << std::endl; } 
    else if( cell_change ){ 
        std::cout << "Crossed surface from " << last_cell 
                  << " into " << *mcnp_current_cell<< std::endl; }
#endif

    last_cell = *mcnp_current_cell;

    // if the new cell is not part of this tally, return immediately
    if( conformality->find( *mcnp_current_cell ) == conformality->end() ) {
      return;
    }
  }

  // Get the physics data for this tally event
  TrackData data;
  bool set_data = event.get_track_data(data);

  if (!set_data)
  {
    std::cerr << "\nError: Invalid track-based event data" << std::endl;
    exit(EXIT_FAILURE);
  }

  EntityHandle last_crossed_tri[3] = {0,0,0};
  double last_t = 0;
  EntityHandle first_tet;

  if( conformal_begin_track ){
    first_tet = get_starting_tet_conformal(data.start_point, last_crossed_tri);
  }
  else{
    first_tet = get_starting_tet(data.start_point, data.direction, data.track_length, last_crossed_tri, last_t);
  }

  if( first_tet == 0 ){
    // this ray never touches the tally mesh
    return;
  }

  EntityHandle next_tet = first_tet;
  //double first_t = last_t; //useful as a debugging var
  int tet_count = 0;

  while( next_tet != 0 ){

    EntityHandle tet = next_tet; // the tetrahedron being currently handled
    tet_count++;

    const EntityHandle* tet_verts; 
    int num_tet_verts; 
    rval = mb->get_connectivity( tet, tet_verts, num_tet_verts );
    assert( rval == MB_SUCCESS );
    assert( num_tet_verts == 4 );

    bool found_crossing = false;


    for( int i = 0; i < 4 && !found_crossing; ++i ){

      EntityHandle tri[3]; 
      int three;

      CN::SubEntityConn( tet_verts, MBTET, 2, i, tri, three );
      assert( three == 3 );

      if( tris_eq( tri, last_crossed_tri ) ) continue;

      CartVect tri_pts[3];
      rval = mb->get_coords( tri, 3, tri_pts[0].array() );
      assert( rval == MB_SUCCESS );

      double t;
      if( GeomUtil::ray_tri_intersect( tri_pts, data.start_point, data.direction, TRIANGLE_INTERSECTION_TOL, t ) ){

        double track_length;

        if( t >= data.track_length ){
          // track ends in this tetrahedron
          track_length = data.track_length - last_t;
          next_tet = 0;
          this->last_visited_tet = tet;

#ifdef MESHTAL_DEBUG
          std::cout << " track ends in mesh" << std::endl;
#endif

        }
        else{ 
          // t < length, track proceeds to next tetrahedron
          track_length = t - last_t;
          last_t = t;
          
          Range tri_sides; 
          rval = mb->get_adjacencies( tri, 3, 3, false, tri_sides );
          assert( rval == MB_SUCCESS );
          assert( tri_sides.size() <= 2 );
          assert( tri_sides[0] == tet || tri_sides[1] == tet );
          assert( tri_sides.size() == 1 || (tri_sides[0] != tet || tri_sides[1] != tet) );

          if( tri_sides.size() == 1 ){
            // the mesh ends here
            CartVect crossing = data.start_point + (data.direction*t);
#ifdef MESHTAL_DEBUG
            std::cout << "  Ray departs from mesh at t = " << t << ", "  << crossing  << std::endl;
#endif
      

            if( convex || conformality ) {
              // input file made assertion that this mesh tally is convex, 
              // or conformality assures that no single track will reenter tally mesh
            next_tet = 0;
            }
            else{
              Range last_tri_eh;
              rval = mb->get_adjacencies( tri, 3, 2, false, last_tri_eh );
              assert( rval == MB_SUCCESS );
              assert( last_tri_eh.size() == 1 );
              next_tet = find_next_tet_by_ray_fire( crossing, data.direction, data.track_length - t, last_crossed_tri, last_t, last_tri_eh[0] );
            }

#ifdef MESHTAL_DEBUG
            if( next_tet ){
              std::cout << "  Ray enters mesh again at t = " << last_t << std::endl;
            }
            else{
              std::cout << " track ends." << std::endl;
            }
#endif
            last_t += t;
            this->last_visited_tet = 0;
          }
          else{ 
#ifdef MESHTAL_DEBUG
            std::cout << "  track proceeds, t = " << last_t << std::endl;
#endif
            next_tet = (tri_sides[0]==tet) ? tri_sides[1] : tri_sides[0];
            memcpy( last_crossed_tri, tri, sizeof(last_crossed_tri) );
          }

        }
        // If track length is negative do not add score to the tally and end track early
        if( track_length < 0 ){
          std::cerr << "Warning-- negative track length.   Dubious continuation.." << std::endl;
          ++num_negative_tracks;
          return;
        }

        double weight = event.get_weighting_factor();
        double score = weight * track_length;

        add_score_to_mesh_cell( tet, score, ebin );
        found_crossing = true;
      }

    }
    if( !found_crossing ){
      std::cerr << "ERROR: No triangle crossings found.  Ending a meshtal track early." << std::endl;
      return;
    }
    assert( found_crossing );

  }

  return;

}

} // namespace moab
