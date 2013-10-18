// MCNP5/dagmc/TrackLengthMeshTally.cpp

#include <iostream>
#include <sstream>
#include <cmath> 
#include <set>

#include "moab/Core.hpp"
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
#include "TrackLengthMeshTally.hpp"

// tolerance for ray-triangle intersection tests
// (note: this paramater is ignored by GeomUtil, so don't bother trying to tune it)
#define TRIANGLE_INTERSECTION_TOL 1e-6

// Stores tag name and values expected in input mesh
std::string tag_name;
std::vector<std::string> tag_values;

// used to store the intersection data
struct ray_data {
  double intersect;
  moab::EntityHandle triangle;
};

//---------------------------------------------------------------------------//
// MISCELLANEOUS FILE SCOPE METHODS
//---------------------------------------------------------------------------/
/* sorting function */
// used to sort the ray triangle intersection data
inline static bool compare(const ray_data &a, const ray_data &b)
{
    return a.intersect < b.intersect;
}

// Adapted from MOAB's convert.cpp
// Parse list of integer ranges, e.g. "1,2,5-10,12"
static bool parse_int_list(const char* string, std::set<int>& results)
{
  bool okay = true;
  char* mystr = strdup( string );
  for (const char* ptr = strtok(mystr, ", \t"); ptr; ptr = strtok(0,", \t"))
  {
    char* endptr;
    long val = strtol( ptr, &endptr, 0 );
    if (endptr == ptr)
    {
      std::cerr << "Not an integer: \"" << ptr << '"' << std::endl;
      okay = false;
      break;
    }
    
    long val2 = val;
    if (*endptr == '-')
    {
      const char* sptr = endptr+1;
      val2 = strtol( sptr, &endptr, 0 );
      if (endptr == sptr) 
      {
        std::cerr << "Not an integer: \"" << sptr << '"' << std::endl;
        okay = false;
        break;
      }
      if (val2 < val)
      {
        std::cerr << "Invalid id range: \"" << ptr << '"' << std::endl;
        okay = false;
        break;
      }
    }
    
    if (*endptr) 
    {
      okay = false;
      break;
    }
    
    for (; val <= val2; ++val)
      results.insert( (int)val );
  }
  
  free( mystr );
  return okay;    
}

/* Tetrahedron volume code taken from MOAB/tools/measure.cpp */
inline static double tet_volume(const moab::CartVect& v0,
                                const moab::CartVect& v1,
                                const moab::CartVect& v2,
                                const moab::CartVect& v3)
{
  return 1./6. * ( ((v1 - v0) * (v2 - v0)) % (v3 - v0) );
}
//---------------------------------------------------------------------------//
static inline bool tris_eq(const moab::EntityHandle *t1,
                           const moab::EntityHandle *t2)
{
  int ignored1, ignored2;
  return moab::CN::ConnectivityMatch( t1, t2, 3, ignored1, ignored2 );
}

namespace moab { 
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TrackLengthMeshTally::TrackLengthMeshTally(const TallyInput& input)
    : MeshTally(input),
      mb (new moab::Core()),  
      obb_tool(new OrientedBoxTreeTool(mb)),
      last_visited_tet(0),
      last_cell(-1),
      convex(false),
      conformal_surface_source(false),
      num_negative_tracks(0)
{
   std::cout << "Creating dagmc mesh tally" << input.tally_id 
            << ", input: " << input_filename 
            << ", output: " << output_filename << std::endl;

   parse_tally_options();
   set_tally_meshset();

   // reduce the loaded MOAB mesh set to include only 3D elements
   Range all_tets;
   ErrorCode rval = reduce_meshset_to_3D(mb, tally_mesh_set, all_tets);  
   if(rval != MB_SUCCESS)
     {
       std::cout << "Failed to reduce meshset to 3d" << std::endl;
       exit(1);
     }
   assert (rval == MB_SUCCESS);

   // initialize MeshTally::tally_points to include all mesh cells
   set_tally_points(all_tets);
   // build the kdtree
   // Does not change all_tets
   rval = compute_barycentric_data(all_tets);
   assert (rval == MB_SUCCESS);
  
   // Perform tasks
   rval = setup_tags( mb );
   assert (rval == MB_SUCCESS);

   build_trees(all_tets);

}
//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
TrackLengthMeshTally::~TrackLengthMeshTally()
{
  delete mb;
  delete obb_tool;
}



//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from Tally.hpp
//---------------------------------------------------------------------------//
void TrackLengthMeshTally::compute_score(const TallyEvent& event)
{
  EntityHandle tet; // tet 

  // If it's not the type we want leave immediately
  if (event.type != TallyEvent::TRACK) return;
  
  std::vector<double> intersections;     // vector of distance to triangular facet intersections
  std::vector< EntityHandle > triangles; // vector of entityhandles that belong to the triangles hit
  //  std::vector<ray_data> hit_information; // vector of reformated ray triangle intersections

  // get all ray-triangle intersections along the ray 
  ErrorCode rval = get_all_intersections(event.position,event.direction,event.track_length,triangles,intersections);
  if (rval != MB_SUCCESS )
    {
      std::cout << "we have a problem finding intersections" << std::endl;
      exit(1);
    }

  if( intersections.size() == 0 )
    // ray is so short it either does not intersect a triangular face, or it inside the mesh
    // but can't reach
    {
       tet = TrackLengthMeshTally::point_in_which_tet(event.position);
      // if tet value is greater than 0 then in a tet, otherwise not
      if( tet == 0 )
	{
	  return;
	}
      else
	{
	  // determine tracklength to return
	  double weight = event.particle_weight;
	  double score = weight * event.track_length;
	  
	  // ToDo:  fix fake ebin
	  int ebin = 0;
	  unsigned int tet_index = get_entity_index(tet);
	  data->add_score_to_tally(tet_index, score, ebin);
	  //	  found_crossing = true;
	  return;
	}
    }

  // sort the intersection data
  sort_intersection_data(intersections,triangles);
  // compute the tracklengths
  compute_tracklengths(event,intersections,triangles);


  return;
}


//---------------------------------------------------------------------------//
// This may not need to be overridden, depending on whether conformality
void TrackLengthMeshTally::end_history () 
{
  MeshTally::end_history();
  if( !conformality.empty() ){ last_cell = -1; } 
}

//---------------------------------------------------------------------------//
void TrackLengthMeshTally::write_data(double num_histories)
{
  ErrorCode rval;

  Range all_tets;
  rval = mb->get_entities_by_dimension( tally_mesh_set, 3, all_tets );
  if(rval != MB_SUCCESS )
    {
      std::cout << "Failed to get 3d entities" << std::endl;
      exit(1);
    }
  assert( rval == MB_SUCCESS );

  for( Range::const_iterator i=all_tets.begin(); i!=all_tets.end(); ++i){
    EntityHandle t = *i;

    CartVect v[4];

    std::vector<EntityHandle> vtx;    
    mb->get_connectivity( &t, 1, vtx );
    assert( vtx.size() == 4);

    int k = 0;
    for( std::vector<EntityHandle>::iterator j = vtx.begin(); j!=vtx.end(); ++j)
    {
      EntityHandle vertex = *j;
      mb->get_coords( &vertex, 1, v[k++].array() );
    }

    double volume = tet_volume( v[0], v[1], v[2], v[3] );
    unsigned int tet_index = get_entity_index(t);

    for( unsigned j = 0; j < data->get_num_energy_bins(); ++j )
    {
      std::pair <double,double> tally_data = data->get_data(tet_index,j);
      double tally = tally_data.first;
      double error = tally_data.second;
      // double tally = get_data( tally_data, t, j );
      // double error = get_data( error_data, t, j );
      double score = (tally / (volume*num_histories));
      
      rval = mb->tag_set_data( tally_tags[j], &t, 1, &score );
      assert( rval == MB_SUCCESS );

      // Use 0 as the error output value if nothing has been computed for this mesh cell;
      // this reflects MCNP's approach to avoiding a divide-by-zero situation.
      double rel_err = 0;
      if( error != 0 )
      {
        rel_err = sqrt( (error / (tally*tally)) - (1./num_histories) );
      }        

      rval = mb->tag_set_data( error_tags[j], &t, 1, &rel_err );
      assert( rval == MB_SUCCESS );
    }
  }

  std::vector<Tag> output_tags = tally_tags;
  output_tags.insert( output_tags.end(), error_tags.begin(), error_tags.end() );

  rval = mb->write_file( output_filename.c_str(), NULL, NULL, &tally_mesh_set, 1, &(output_tags[0]), output_tags.size() );
  assert (rval == MB_SUCCESS );
 
}
//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
void TrackLengthMeshTally::parse_tally_options()
{
  const TallyInput::TallyOptions& options = input_data.options;
  TallyInput::TallyOptions::const_iterator it;

  for(it = options.begin(); it != options.end(); ++it )
  {
    std::string key = (*it).first, val = (*it).second;
    if( key == "tag" ) tag_name = val;
    else if( key == "tagval" ) tag_values.push_back(val);
    else if( key == "convex" && (val == "t" || val == "true" ) ) convex = true; 
    else if( key == "conf_surf_src" && (val == "t" || val == "true" ) ) conformal_surface_source = true;
    else if( key == "conformal" ) 
    { 
      // Since the options are a multimap, the conformal tag could (illogically) occur more than once
      if (conformality.empty())
      {
         if( !parse_int_list( val.c_str(), conformality ) )
         {
           std::cerr << "Error: Tally " << input_data.tally_id << " input has bad conformality value '" << val << "'" << std::endl;
           exit(EXIT_FAILURE);
         }
      }
    }
    else
    {
      std::cerr << "Warning: Tally " << input_data.tally_id << " input has unknown key '" << key << "'" << std::endl;
    }
  }
  if( tag_name != "" )
  {
    std::cout << "  using tag name='" << tag_name << "'";
    if( tag_values.size() > 0 )
    {
      std::cout <<", and tag values= " << std::endl;
      for( unsigned int i = 0; i < tag_values.size(); ++i )
      {
        std::cout << "    '" << tag_values[i] << "'" << std::endl;
      }
    }
    std::cout << std::endl;
  }
}  
//---------------------------------------------------------------------------//
void TrackLengthMeshTally::set_tally_meshset()
{
  // load the MOAB mesh data from the input file for this mesh tally
  moab::EntityHandle loaded_file_set;
  moab::ErrorCode rval = load_moab_mesh(mb, loaded_file_set);
  if(rval != MB_SUCCESS)
    {
      std::cout << "Failed to load moab mesh" << std::endl;
      exit(1);
    }
  assert( rval == MB_SUCCESS ); 

  rval = mb->create_meshset( MESHSET_SET, tally_mesh_set );
  assert( rval == MB_SUCCESS ); 

  if( tag_name.length() > 0 )
  {
    std::cout << "  User-specified tag to load:  " << tag_name << std::endl;

    // Until there is more certainty about the type and parameters of the tag the user specified,
    //   use MB_TAG_ANY to get access to any tag with the given name 
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

    for( Range::iterator i = user_sets.begin(); i!=user_sets.end(); ++i)
    {
      EntityHandle s = *i;
      char* name = new char[ user_spec_tag_length + 1];
      
      rval = mb->tag_get_data( user_spec_tag, &s, 1, name );
      assert( rval == MB_SUCCESS );

      // if user specified no tag value, list the available ones for informational purposes
      if( tag_values.size() == 0 )
      {
        std::cout << "    available tag value: " << name << std::endl; 
      }
      
      if( std::find( tag_values.begin(), tag_values.end(),std::string(name) ) != tag_values.end() )
      {
        std::cout << "  Successfully found a set with tag value " << name << std::endl;
        rval = mb->unite_meshset( tally_mesh_set, s );
        assert( rval == MB_SUCCESS );
      }
      delete[] name;
    }
  }
  else
  { // no user-specified tag filter
    rval = mb->unite_meshset( tally_mesh_set, loaded_file_set );
    if(rval != MB_SUCCESS)
      {
	std::cout << "Failed to unite meshset" << std::endl;
	exit(1);
      }
    assert (rval == MB_SUCCESS);
  }
} 
//---------------------------------------------------------------------------//
ErrorCode TrackLengthMeshTally::compute_barycentric_data(const Range& all_tets)
{
  ErrorCode rval;

  // Iterate over all tets and compute barycentric matrices 
  int num_tets = all_tets.size();
  std::cerr << "  There are " << num_tets << " tetrahedrons in this tally mesh." << std::endl;

  if (num_tets != 0)
  {
     tet_baryc_data.resize (num_tets);  
  }

  for( Range::const_iterator i=all_tets.begin(); i!=all_tets.end(); ++i)
  {
    EntityHandle tet = *i;

    const EntityHandle* verts;
    int num_verts;
    rval = mb->get_connectivity (tet, verts, num_verts);
    if(rval != MB_SUCCESS)
      {
	std::cout << "Failed to get connectivity information" << std::endl;
	exit(1);
      }
    assert( rval == MB_SUCCESS );
    
    if( num_verts != 4 )
    {
      std::cerr << "Error: DAGMC TrackLengthMeshTally cannot handle non-tetrahedral meshes yet," << std::endl;
      std::cerr << "       but your mesh has at least one cell with " << num_verts << " vertices." << std::endl;
      return MB_NOT_IMPLEMENTED;
    }
    
    CartVect p[4];
    rval = mb->get_coords (verts, 4, p[0].array());
    if(rval != MB_SUCCESS)
      {
	std::cout << "Failed to get coordinate data" << std::endl;
	exit(1);
      }
    assert( rval == MB_SUCCESS );

    Matrix3 a( p[1]-p[0], p[2]-p[0], p[3]-p[0] );
    a = a.transpose().inverse();
    tet_baryc_data.at( get_entity_index(tet) ) = a;
  }
  return MB_SUCCESS;
}
//---------------------------------------------------------------------------//
void TrackLengthMeshTally::build_trees (Range& all_tets)
{
  // prepare to build KD tree and OBB tree
  Range all_tris;
  // get the triangles that belong to the mesh
  Range new_triangles = TrackLengthMeshTally::get_adjacency_info(all_tets);

  std::cout << "  Tally mesh has " << new_triangles.size() << " triangles." << std::flush;

  // put tris with tets to be rolled into KD tree
  all_tets.merge( new_triangles );

  // build KD tree of all tetrahedra and triangles
  std::cout << "  Building KD tree of size " << all_tets.size() << "... " << std::flush;
  kdtree = new AdaptiveKDTree( mb );
  kdtree->build_tree( all_tets, kdtree_root );
  std::cout << "done." << std::endl << std::endl;;
}
//---------------------------------------------------------------------------//
bool TrackLengthMeshTally::point_in_tet(const CartVect& point,
                                        const EntityHandle* tet)
{ 
  ErrorCode rval;
                    
  const EntityHandle* verts;
  int num_verts;
  rval = mb->get_connectivity( *tet, verts, num_verts );
  if(rval != MB_SUCCESS)
    {
      std::cout << "Failed to get connectivity information" << std::endl;
      exit(1);
    }
  assert( rval == MB_SUCCESS );
  
  CartVect p0;
  rval = mb->get_coords( verts, 1, p0.array() );
  if(rval != MB_SUCCESS)
    {
      std::cout << "Failed to get coordinate information" << std::endl;
      exit(1);
    }
  assert( rval == MB_SUCCESS );

  Matrix3& Ainverse = tet_baryc_data[ get_entity_index(*tet) ];

  CartVect bary = (Ainverse) * (point-p0);
  
  bool in_tet = ( bary[0]>= 0 && bary[1] >= 0 && bary[2] >= 0 &&
                  bary[0]+bary[1]+bary[2] <= 1. );

  return in_tet;
}
//---------------------------------------------------------------------------//
/*
 * return the list of intersections
 */
ErrorCode TrackLengthMeshTally::get_all_intersections(CartVect position, CartVect direction, double track_length, 
				std::vector<EntityHandle> &triangles,std::vector<double> &intersections)
  {
    
    ErrorCode result = kdtree->ray_intersect_triangles( kdtree_root, TRIANGLE_INTERSECTION_TOL, 
					    direction.array(), position.array(), triangles,
					    intersections, 0, track_length );
    return result;
  }


/*
 * loop through all tets to find which one we are in
 */
EntityHandle TrackLengthMeshTally::point_in_which_tet (CartVect point)
{
  ErrorCode rval;
  AdaptiveKDTreeIter tree_iter;
  
  // Check to see if starting point begins inside a tet
  rval = kdtree->leaf_containing_point( kdtree_root, point.array(), tree_iter );
  if( rval == MB_SUCCESS )
    {
      EntityHandle leaf = tree_iter.handle();
      Range candidate_tets;
      rval = mb->get_entities_by_dimension( leaf, 3, candidate_tets, false );
      assert( rval == MB_SUCCESS );
      for( Range::const_iterator i = candidate_tets.begin(); i!=candidate_tets.end(); ++i)
	{
	  if( point_in_tet( point,&(*i)) )
	    {
	      return *i;
	    }
	}
    }
  return 0;
}

/*
 * Returns the tet in which the remaining length to score ends in
 */
EntityHandle TrackLengthMeshTally::remainder( CartVect start, CartVect dir, double distance, double left_over)
{
  CartVect pos_check = start+(dir*(distance+left_over));
  return TrackLengthMeshTally::point_in_which_tet (pos_check);
}

/* 
 * wrapper function to return the triangles that belong to each tet
 * including no duplicate entities for shared triangles
 */
Range TrackLengthMeshTally::get_adjacency_info(Range input_handles)
{
  Range adjacencies,test;
  Range::iterator inh;
  ErrorCode rval;
  int dimension;
  
  // get the dimensionality of the input set
  dimension = mb->dimension_from_handle(input_handles[0]);
  
  // generate all edges for these tets
  rval = mb->get_adjacencies(input_handles, dimension-1, true, 
			     adjacencies, Interface::UNION); 
  
  if( rval!=MB_SUCCESS )
    std::cout << "ERROR could not determine adjacancy data" << std::endl;
  
  return adjacencies;
}

// function to sort the ray-mesh intersection data
void TrackLengthMeshTally::sort_intersection_data(std::vector<double> &intersections,
						  std::vector<EntityHandle> &triangles)
{
  std::vector<ray_data> hit_information; // vector of reformated ray triangle intersections)

  // copy the data from intersections and triangles to temporary array of structs for the
  // purpose of sorting the intersection data
  for ( unsigned int i = 0 ; i < intersections.size() ; i++ ) 
    {
      ray_data data;
      data.intersect=intersections[i];
      data.triangle=triangles[i];
      hit_information.push_back(data);
    }
  // at some point the sort will be done by moab rather than by us
  std::sort(hit_information.begin(), hit_information.end(), compare); // sort the intersections

  // copy the sorted data back into array
  for ( unsigned int i = 0 ; i < intersections.size() ; i++) 
    {
      intersections[i]=hit_information[i].intersect;
      triangles[i]=hit_information[i].triangle;
    }

  return;
}

// determine the score for the given tet
void TrackLengthMeshTally::determine_score(const TallyEvent event, double tracklength, const EntityHandle tet)
{
    double weight = event.particle_weight;
    double score = weight * tracklength;
    
    // ToDo:  fix fake ebin
    int ebin = 0;
    unsigned int tet_index = get_entity_index(tet);
    data->add_score_to_tally(tet_index, score, ebin);
    return;
}

// function to compute the track lengths
void TrackLengthMeshTally::compute_tracklengths(const TallyEvent event, 
						const std::vector<double> intersections,
						const std::vector<EntityHandle> triangles)
{
  double track_length; // track_length to add to the tet
  CartVect hit_p; //position on the triangular face of the hit
  std::vector<CartVect> hit_point; // array of all hit points
  CartVect tet_centroid; // centroid position between intersect point
  EntityHandle tet;
  hit_point.push_back(event.position); // add the origin of the ray to the point to the list

  EntityHandle next_tet = 0;
  // loop over all intersections
  for (unsigned int i = 0 ; i < intersections.size() ; i++) 
    {
      // hit point
      //      std::cout << "pos " << event.position << std::endl;
      hit_p = (event.direction*intersections[i]) + event.position;
      hit_point.push_back(hit_p); // add to list of hit points
      tet_centroid = ((hit_point[i+1]-hit_point[i])/2.0)+hit_point[i]; // centre of the tet
      //      std::cout << "centroid " << tet_centroid << std::endl;
      // determine the tet that the point belongs to
      tet = TrackLengthMeshTally::point_in_which_tet(tet_centroid);
      //      std::cout << tet << std::endl;
      if ( tet > 0 )
	{
	  if ( i != 0 ) // detrmine the track_length
	    track_length = intersections[i]-intersections[i-1];
	  else
	    track_length = intersections[i];

	  //  std::cout << tet << " " << track_length << " " << intersections[i] << std::endl;

	  if (track_length < 0.0 )
	    {
	      std::cout << "!!! Negative Track Length !!!" << std::endl;
	      std::cout << track_length << " " << intersections[i] << " " << intersections[i-1] << std::endl;
	      std::cout << tet << " " << next_tet << std::endl;
	    }

	  TrackLengthMeshTally::determine_score(event,track_length,tet);
	  /*
	  double weight = event.particle_weight;
	  double score = weight * track_length;
	  
	  // ToDo:  fix fake ebin
	  int ebin = 0;
	  unsigned int tet_index = get_entity_index(tet);
	  data->add_score_to_tally(tet_index, score, ebin);
	  */
	}
    }

  // it is possible that there is some tracklength left to allocate at the end, where the ray ends in the middle of a tet
  // or the ray could end in free space
  if ( intersections[intersections.size()-1] < event.track_length )
    {
      track_length = event.track_length-intersections[intersections.size()-1];
      tet = TrackLengthMeshTally::remainder(event.position,event.direction,
					   intersections[intersections.size()-1],
					   track_length);
      if (track_length < 0.0 )
	{
	  std::cout << "Negative Track Length!!" << std::endl;
	  std::cout << track_length << " " << intersections[intersections.size()-1] << " " <<  event.track_length << std::endl;
	  std::cout << tet << " " << next_tet << std::endl;
	}
      
      if ( tet > 0 ) 
	{
	  TrackLengthMeshTally::determine_score(event,track_length,tet);
	  /*
	  double weight = event.particle_weight;
	  double score = weight * track_length;
	  
	  // ToDo:  fix fake ebin
	  int ebin = 0;
	  unsigned int tet_index = get_entity_index(tet);
	  data->add_score_to_tally(tet_index, score, ebin);
	  */
	}
    }

}

//---------------------------------------------------------------------------//

} // end namespace moab

// end of MCNP5/dagmc/TrackLengthMeshTally.cpp
