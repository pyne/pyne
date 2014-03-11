// FLUDag/src/moab_utils.cpp

//----------------------------------*-C++ -*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/moab_utils.cpp
 * \author Julie Zachman 
 * \date   Thu Dec 19 2013 
 * \brief  Functions called for parse_flu and others
 * \note   After DagMC.cpp
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//
#include "moab_utils.hpp"

#include "moab/Core.hpp"
#include "MBCartVect.hpp"


#include <iostream>  // cerr
#include <iomanip>
#include <sstream>
#include <set>
#include <cstring>
#include <string>

// globals

#include <numeric>
#include <fstream>

using namespace moab;

  int MCRecordInfo::num_entities( int dimension )
  {
    assert(0 <= dimension && 3 >= dimension);

    return entHandles[dimension].size() - 1;
  }

/* SECTION III */

EntityHandle MCRecordInfo::entity_by_id( Interface* mbi, int dimension, int id )
{
  assert(0 <= dimension && 3 >= dimension);
  const Tag tags[] = { idTag, geomTag };
  const void* const vals[] = { &id, &dimension };
  ErrorCode rval;

  Range results;
  rval = mbi->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, results );

  if ( MB_SUCCESS != rval )
      return 0;

  if ( results.empty() ){
    // old versions of dagmc did not set tags correctly on the implicit complement 'volume',
    // causing it to not be found by the call above.  This check allows this function to work
    // correctly, even on reloaded files from older versions.
    if( dimension == 3 && get_entity_id(mbi, impl_compl_handle) == id )
      return impl_compl_handle;
    else
      return 0;
  }

  return results.front();
}


int MCRecordInfo::get_entity_id(Interface* mbi, EntityHandle this_ent) 
{
  int id = 0;
  ErrorCode result = mbi->tag_get_data(idTag, &this_ent, 1, &id);
  if (MB_TAG_NOT_FOUND == result)
    id = mbi->id_from_handle(this_ent);
    
  return id;
}


int MCRecordInfo::id_by_index( moab::Interface* mbi, int dimension, int index )
{
  EntityHandle h = entity_by_index( dimension, index );
  if (!h)
    return 0;

  int result = 0;
  mbi->tag_get_data( idTag, &h, 1, &result );
  return result;
}


EntityHandle MCRecordInfo::entity_by_index( int dimension, int index )
{
  assert(2 <= dimension && 3 >= dimension && (unsigned) index < entHandles[dimension].size());
  return entHandles[dimension][index];
}


ErrorCode MCRecordInfo::get_impl_compl(Interface* mbi)
{
  Range entities;
  const void* const tagdata[] = {implComplName};
  ErrorCode rval = mbi->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                         &nameTag, tagdata, 1,
                                                         entities );
  // query error
  if (MB_SUCCESS != rval) {
    std::cerr << "Unable to query for implicit complement." << std::endl;
    return rval;
  }

  // found too many
  if (entities.size() > 1) {
    std::cerr << "Too many implicit complement sets." << std::endl;
    return MB_MULTIPLE_ENTITIES_FOUND;
  }

  // found none
  if (entities.empty()) {
    rval = mbi->create_meshset(MESHSET_SET,impl_compl_handle);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to create mesh set for implicit complement." << std::endl;
      return rval;
    }
      // tag this entity with name for implicit complement
    rval = mbi->tag_set_data(nameTag,&impl_compl_handle,1,&implComplName);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to tag new entity as implicit complement." << std::endl;
    }

    return rval;

  } else {
    // found a single implicit complement
    impl_compl_handle = entities.front();
    return MB_SUCCESS;
  }
}

//---------------------------------------------------------------------------//
// make_property_string
//---------------------------------------------------------------------------//
// For a given volume, find all properties associated with it, and any and all 
//     values associated with each property
// Copied and modified from obb_analysis.cpp
/*
static std::string MCRecordInfo::make_property_string (Interface* mbi, EntityHandle eh, std::vector<std::string> &properties)
{
  ErrorCode ret;
  std::string propstring;
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( has_prop(mbi, eh, *p) )
     {
        std::vector<std::string> vals;
        // ret = DAG->prop_values(eh, *p, vals);
        ret = prop_values(mbi, eh, *p, vals);
        // CHECKERR(*DAG, ret);
        propstring += *p;
        if (vals.size() == 1)
        {
 	   propstring += "=";
           propstring += vals[0];
        }
        else if (vals.size() > 1)
        {
 	   // this property has multiple values; list within brackets
           propstring += "=[";
	   for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
           {
	       propstring += *i;
               propstring += ",";
           }
           // replace the last trailing comma with a close bracket
           propstring[ propstring.length() -1 ] = ']';
        }
        propstring += ", ";
     }
  }
  if (propstring.length())
  {
     propstring.resize( propstring.length() - 2); // drop trailing comma
  }
  return propstring;
}
*/
/************************* Modified from DagMC *****************************/
/* SECTION V: Metadata handling */

// ToDo:  Consider using StringSplit
void MCRecordInfo::tokenize( const std::string& str,
                      std::vector<std::string>& tokens,
                      const char* delimiters ) 
                      // const char* delimiters ) const
{
  std::string::size_type last = str.find_first_not_of( delimiters, 0 );
  std::string::size_type pos  = str.find_first_of( delimiters, last );
  if ( std::string::npos == pos )
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back( str.substr( last, pos - last ) );
      last = str.find_first_not_of( delimiters, pos );
      pos  = str.find_first_of( delimiters, last );
      if(std::string::npos == pos)
        pos = str.size();
    }
}

ErrorCode MCRecordInfo::get_group_name( moab::Interface* mbi, EntityHandle group_set, std::string& name )
{
  ErrorCode rval;
  const void* v = NULL;
  int ignored;
  rval = mbi->tag_get_by_ptr(name_tag(), &group_set, 1, &v, &ignored);
  if( MB_SUCCESS != rval ) return rval;
  name = static_cast<const char*>(v);
  return MB_SUCCESS;
}

ErrorCode MCRecordInfo::append_packed_string( moab::Interface* mbi,  Tag tag, EntityHandle eh,
                                       std::string& new_string )
{
    // When properties have multiple values, the values are tagged in a single character array
    // with the different values separated by null characters
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = mbi->tag_get_by_ptr( tag, &eh, 1, &p, &len );
  if( rval == MB_TAG_NOT_FOUND ){
    // This is the first entry, and can be set directly
    p = new_string.c_str();
    return mbi->tag_clear_data( tag, &eh, 1, p, new_string.length()+1);
  }
  else if( rval != MB_SUCCESS ) return rval;
  else{
    str = static_cast<const char*>(p);
  }

  // append a new value for the property to the existing property string
  unsigned int tail_len = new_string.length() + 1;
  char* new_packed_string = new char[ len + tail_len ];
  memcpy( new_packed_string, str, len );
  memcpy( new_packed_string + len, new_string.c_str(), tail_len );

  int new_len = len + tail_len;
  p = new_packed_string;
  rval = mbi->tag_set_by_ptr( tag, &eh, 1, &p, &new_len );
  delete[] new_packed_string;
  return rval;
}


ErrorCode MCRecordInfo::unpack_packed_string( moab::Interface* mbi, Tag tag, EntityHandle eh,
                                       std::vector< std::string >& values )
{
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = mbi->tag_get_by_ptr( tag, &eh, 1, &p, &len );
  if( rval != MB_SUCCESS ) return rval;
  str = static_cast<const char*>(p);
  int idx = 0;
  while( idx < len ){
    std::string item(str + idx);
    values.push_back( item );
    idx += item.length() + 1;
  }
  return MB_SUCCESS;
}

ErrorCode MCRecordInfo::parse_properties( Interface* mbi,  const std::vector<std::string>& keywords,
                              const std::map<std::string, std::string>& keyword_synonyms )
{
  ErrorCode rval;

  // master keyword map, mapping user-set words in cubit to canonical property names
  std::map< std::string, std::string > keyword_map( keyword_synonyms );

  for( std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i )
  {
    keyword_map[*i] = *i;
  }

  // the set of all canonical property names
  std::set< std::string > prop_names;
  for( prop_map::iterator i = keyword_map.begin();
       i != keyword_map.end(); ++i )
  {
    prop_names.insert((*i).second);
  }

  // set up DagMC's property tags based on what's been requested
  for( std::set<std::string>::iterator i = prop_names.begin();
       i != prop_names.end(); ++i )
  {
    std::string tagname("DAGMCPROP_");
    tagname += (*i);

    Tag new_tag;
    rval = mbi->tag_get_handle( tagname.c_str(), 0, MB_TYPE_OPAQUE, new_tag,
                                MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_CREAT );
    if( MB_SUCCESS != rval ) return rval;
    property_tagmap[(*i)] = new_tag;
  }

  // now that the keywords and tags are ready, iterate over all the actual geometry groups
  for( std::vector<EntityHandle>::iterator grp=group_handles().begin();
       grp != group_handles().end(); ++grp )
  {

    prop_map properties;
    rval = parse_group_name(mbi, *grp, properties );
    if( rval == MB_TAG_NOT_FOUND ) continue;
    else if( rval != MB_SUCCESS ) return rval;

    Range grp_sets;
    rval = mbi->get_entities_by_type( *grp, MBENTITYSET, grp_sets);
    if( MB_SUCCESS != rval ) return rval;
    if( grp_sets.size() == 0 ) continue;

    for( prop_map::iterator i = properties.begin();
         i != properties.end(); ++i )
    {
      std::string groupkey = (*i).first;
      std::string groupval = (*i).second;

      if( property_tagmap.find( groupkey ) != property_tagmap.end() ){
        Tag proptag = property_tagmap[groupkey];
        const unsigned int groupsize = grp_sets.size();
        for( unsigned int j = 0; j < groupsize; ++j){
            rval = append_packed_string( mbi, proptag, grp_sets[j], groupval );
        }
      }
    }
  }
  return MB_SUCCESS;
}


ErrorCode MCRecordInfo::prop_value( moab::Interface* mbi, EntityHandle eh, const std::string& prop, std::string& value )
{
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return MB_TAG_NOT_FOUND;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = mbi->tag_get_by_ptr( proptag, &eh, 1, &data, &ignored );
  if( rval != MB_SUCCESS ) return rval;
  value = static_cast<const char*>(data);
  return MB_SUCCESS;
}

ErrorCode MCRecordInfo::prop_values( moab::Interface* mbi, EntityHandle eh, const std::string& prop,
                              std::vector< std::string >& values )
{

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  return unpack_packed_string(mbi, proptag, eh, values );

}

bool MCRecordInfo::has_prop( moab::Interface* mbi, EntityHandle eh, const std::string& prop )
{
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if( it == property_tagmap.end() ){
      return false;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = mbi->tag_get_by_ptr( proptag, &eh, 1, &data, &ignored );
  return ( rval == MB_SUCCESS );
}

ErrorCode MCRecordInfo::detect_available_props(moab::Interface* mbi, std::vector<std::string>& keywords_list )
{
  ErrorCode rval;
  std::set< std::string > keywords;
  for( std::vector<EntityHandle>::const_iterator grp=group_handles().begin();
       grp != group_handles().end(); ++grp )
  {
    std::map< std::string, std::string > properties;
    rval = parse_group_name(mbi, *grp, properties );
    if( rval == MB_TAG_NOT_FOUND ) continue;
    else if( rval != MB_SUCCESS ) return rval;

    for( prop_map::iterator i = properties.begin();
         i != properties.end(); ++i )
    {
      keywords.insert( (*i).first );
    }
  }
  keywords_list.assign( keywords.begin(), keywords.end() );
  return MB_SUCCESS;
}

ErrorCode MCRecordInfo::parse_group_name(moab::Interface*mbi, EntityHandle group_set, prop_map& result )
{
  ErrorCode rval;
  std::string group_name;
  rval = get_group_name(mbi, group_set, group_name );
  if( rval != MB_SUCCESS ) return rval;

  std::vector< std::string > group_tokens;
  tokenize( group_name, group_tokens, "_" );

  // iterate over all the keyword positions 
  // keywords are even indices, their values (optional) are odd indices
  for( unsigned int i = 0; i < group_tokens.size(); i += 2 ){
    std::string groupkey = group_tokens[i];
    std::string groupval;
    if( i < group_tokens.size() - 1 )
      groupval = group_tokens[i+1];
    result[groupkey] = groupval;
  }
  return MB_SUCCESS;
}


// calculate volume of polyhedron
ErrorCode MCRecordInfo::measure_volume( moab::Interface* mbi, EntityHandle volume, double& result )
{
  ErrorCode rval;
  std::vector<EntityHandle> surfaces, surf_volumes;
  result = 0.0;
  
   // don't try to calculate volume of implicit complement
  if (volume == impl_compl_handle) {
    result = 1.0;
    return MB_SUCCESS;
  }

    // get surfaces from volume
  rval = mbi->get_child_meshsets( volume, surfaces );
  if (MB_SUCCESS != rval) return rval;
  
    // get surface senses
  std::vector<int> senses( surfaces.size() );
  rval = surface_sense( mbi, volume, surfaces.size(), &surfaces[0], &senses[0] );
  if (MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }
  
  for (unsigned i = 0; i < surfaces.size(); ++i) {
      // skip non-manifold surfaces
    if (!senses[i])
      continue;
    
      // get triangles in surface
    Range triangles;
    rval = mbi->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (MB_SUCCESS != rval) 
      return rval;
    if (!triangles.all_of_type(MBTRI)) {
      std::cout << "WARNING: Surface " << get_entity_id(mbi, surfaces[i])
                << " contains non-triangle elements. Volume calculation may be incorrect." 
                << std::endl;
      triangles.clear();
      rval = mbi->get_entities_by_type( surfaces[i], MBTRI, triangles );
      if (MB_SUCCESS != rval) return rval;
    }
    
      // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const EntityHandle *conn;
    int len;
    CartVect coords[3];
    for (Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = mbi->get_connectivity( *j, conn, len, true );
      if (MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = mbi->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;
    
      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }
  
  result /= 6.0;
  return MB_SUCCESS;
}


// get sense of surface(s) wrt volume
ErrorCode MCRecordInfo::surface_sense( Interface* mbi, EntityHandle volume, 
                           int num_surfaces,
                           const EntityHandle* surfaces,
                           int* senses_out )
{

  /* The sense tags do not reference the implicit complement handle.
     All surfaces that interact with the implicit complement should have
     a null handle in the direction of the implicit complement. */
  //if (volume == impl_compl_handle)
  //  volume = (EntityHandle) 0;

  std::vector<EntityHandle> surf_volumes( 2*num_surfaces );
  ErrorCode rval = mbi->tag_get_data( sense_tag(), surfaces, num_surfaces, &surf_volumes[0] );
  if (MB_SUCCESS != rval)  return rval;
  
  const EntityHandle* end = surfaces + num_surfaces;
  std::vector<EntityHandle>::const_iterator surf_vols = surf_volumes.begin();
  while (surfaces != end) {
    EntityHandle forward = *surf_vols; ++surf_vols;
    EntityHandle reverse = *surf_vols; ++surf_vols;
    if (volume == forward) 
      *senses_out = (volume != reverse); // zero if both, otherwise 1
    else if (volume == reverse)
      *senses_out = -1;
    else 
      return MB_ENTITY_NOT_FOUND;
    
    ++surfaces;
    ++senses_out;
  }
  
  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
/*
ErrorCode DagMC::surface_sense( EntityHandle volume, 
                                  EntityHandle surface,
                                  int& sense_out )
{
  // The sense tags do not reference the implicit complement handle.
  // All surfaces that interact with the implicit complement should have
  // a null handle in the direction of the implicit complement. 
  //if (volume == impl_compl_handle)
  //  volume = (EntityHandle) 0;

    // get sense of surfaces wrt volumes
  EntityHandle surf_volumes[2];
  ErrorCode rval = MBI->tag_get_data( sense_tag(), &surface, 1, surf_volumes );
  if (MB_SUCCESS != rval)  return rval;
  
  if (surf_volumes[0] == volume)
    sense_out = (surf_volumes[1] != volume); // zero if both, otherwise 1
  else if (surf_volumes[1] == volume) 
    sense_out = -1;
  else
    return MB_ENTITY_NOT_FOUND;
  
  return MB_SUCCESS;
}
*/
// }   // namespace moab
