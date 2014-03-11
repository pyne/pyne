// FLUDag/src/moab_utils.hpp

#ifndef FLUDAG_SRC_MOAB_UTILS_HPP
#define FLUDAG_SRC_MOAB_UTILS_HPP

#include <string>
#include <list>
#include <map>
#include <cassert>

#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"

using namespace moab;
//===========================================================================//
/**
 * \class DagUtil
 * \brief Defines a utility class similar to, and a subset of, DAGMC
 *
 *  This class provides utility calls that are based only on MOAB so that 
 *  different physics codes may use them to prepare records for Monte Carlo 
 *  runs, e.g. FLUKA or Geant4
 * 
 *  This is for extracting the information from the geometry file in order to
 *  assist the user of a variety of Monte Carlo codes in preparing the input
 *  records for the actual physics run.
 */
//===========================================================================//
class MCRecordInfo
{
public:

  /** Calculate the volume contained in a 'volume' */
  ErrorCode measure_volume( Interface* mbi, EntityHandle volume, double& result );

  /** Get the sense of surfaces wrt a volume.  Sense values are:
   *  {-1 -> reversed, 0 -> both, 1 -> forward}
   */
  ErrorCode surface_sense( Interface* mbi, EntityHandle volume, 
                           int num_surfaces,
                           const EntityHandle* surfaces,
                           int* senses_out );


  /* Local utilities also defined in DagMC, but referring to MOAB tagging only */
  /** map from dimension & base-1 ordinal index to EntityHandle */
  EntityHandle entity_by_index( int dimension, int index );
  /** map from dimension & base-1 ordinal index to global ID */
  int id_by_index(Interface* mbi,  int dimension, int index );
  /** map from dimension & global ID to EntityHandle */
  EntityHandle entity_by_id( Interface* mbi, int dimension, int id );
  /** PPHW: Missing dim & global ID ==> base-1 ordinal index */
  /** map from EntityHandle to base-1 ordinal index */
  int index_by_handle( EntityHandle handle );
  /** map from EntityHandle to global ID */
  int get_entity_id(Interface* mbi, EntityHandle this_ent);


  /** test for pre-existing implicit complement definition, or return a new one */
  ErrorCode get_impl_compl(Interface* mbi);
  
  /*
   * Prepare a descriptive string that creates the properties of the volume whose index is index
   */
  std::string mat_property_string (Interface* mbi, int index, std::vector<std::string> &properties);

  // Also defined in DagMC as a DAG function
  // metadata


  /** Add a string value to a property tag for a given entity */
  ErrorCode append_packed_string( Interface* mbi, Tag, EntityHandle, std::string& );
  /** Convert a property tag's value on a handle to a list of strings */
  ErrorCode unpack_packed_string( Interface* mbi, Tag tag, EntityHandle eh,
                                  std::vector< std::string >& values );

  // a common type within the property and group name functions
  typedef std::map<std::string, std::string> prop_map;

 /* SECTION V: Metadata handling */
  /** Detect all the property keywords that appear in the loaded geometry
   *
   * @param keywords_out The result list of keywords.  This list could be
   *        validly passed to parse_properties().
   */
  ErrorCode detect_available_props( Interface* mbi, std::vector<std::string>& keywords_out );

  /** Parse properties from group names per metadata syntax standard
   * 
   * @param keywords A list of keywords to parse.  These are considered the canonical
   *                 names of the properties, and constitute the valid inputs to 
   *                 has_prop() and prop_value().
   * @param synonyms An optional mapping of synonym keywords to canonical keywords. 
   *                 This allows more than one group name keyword to take on the same
   *                 meaning
   *                 e.g. if synonyms["rest.of.world"] = "graveyard", then volumes
   *                 in the "rest.of.world" group will behave as if they were in a
   *                 group named "graveyard".
   */
  
  // an empty synonym map to provide as a default argument to parse_properties()
  static const std::map<std::string,std::string> no_synonyms;
  ErrorCode parse_properties( Interface* mbi,  const std::vector<std::string>& keywords,
                              const std::map<std::string, std::string>& keyword_synonyms=no_synonyms );

  /** Get the value of a property on a volume or surface
   *
   * @param eh The entity handle to get a property value on
   * @param prop The canonical property name
   * @param value Output parameter, the value of the property.  If no value was
   *              set on the handle, this will be the empty string.
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from 
   *         MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_value( Interface* mbi, EntityHandle eh, const std::string& prop, std::string& value );

  /** Get the value of a property on a volume or surface
   *
   * @param eh The entity handle to get a property value on
   * @param prop The canonical property name
   * @param values Output parameter, the values of the property will be appended to this list.  If no value was
   *               set on the handle, no entries will be added.
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from 
   *         MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_values(Interface* mbi,  EntityHandle eh, const std::string& prop,
                        std::vector< std::string >& value );

  /** tokenize the metadata stored in group names - basically borroed from ReadCGM.cpp */
  // ToDo:  Consider using StringSplit
  void tokenize( const std::string& str,
                 std::vector<std::string>& tokens,
                 const char* delimiters );

  bool has_prop(Interface* mbi,  EntityHandle eh, const std::string& prop );

  /** Store the name of a group in a string */
  ErrorCode get_group_name(Interface* mbi, EntityHandle group_set, std::string& name );

  /** Parse a group name into a set of key:value pairs */
  ErrorCode parse_group_name(Interface* mbi, EntityHandle group_set, prop_map& result );

  int num_entities( int dimension );

private:

  // a map from the canonical property names to the tags representing them
  std::map<std::string, Tag> property_tagmap;

  // DagMC.hpp private
  std::vector<EntityHandle> entHandles[5];
  std::vector<EntityHandle>& group_handles() {return entHandles[4];}

  // Tag obbTag, geomTag, idTag, nameTag, senseTag, facetingTolTag;
  Tag geomTag, idTag, nameTag, senseTag;
  char implComplName[NAME_TAG_SIZE];
  EntityHandle impl_compl_handle;
  /** get the tag for the "name" of a surface == global ID */
  Tag name_tag() {return nameTag;}
  Tag sense_tag() {return senseTag;}
  bool debug;

}; // class MCRecordInfo

#endif  // FLUDAG_SRC_MOAB_UTILS_HPP
