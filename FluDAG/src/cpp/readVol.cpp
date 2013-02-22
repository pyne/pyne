// Test code that calls DagMC's readVol
#include "DagMC.hpp"
#include <iostream>
#include <stdlib.h>
#include <string>
#include "moab/Types.hpp"

using namespace moab;

#define DAG DagMC::instance()

static std::string make_property_string( DagMC& , EntityHandle , std::vector<std::string> & );

ErrorCode readVol(char *fileptr)
{
  // See if this works
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols in " << fileptr << " is " << num_vols << std::endl;

  ErrorCode code;
  double volume_measure;
  
  // use iterator on Range vols to get EntityHandle vols
/*
  DagMC& dagmc = *DagMC::instance();
  Interface& moab = *dagmc.moab_instance();

  Tag dim_tag = dagmc.geom_tag();
  Range vols;
  const int three = 3;
  const void* ptr = &three;
  code = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, &ptr, 1, vols);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "____iterator_____" << std::endl;
  std::cout << "\tnumber of volumes is " << vols.size() << std::endl;
  Range::iterator iter = vols.begin();
  
  for (unsigned i = 0; i<vols.size(); ++i, ++iter)
  {
      code = dagmc.measure_volume(*iter, volume_measure);
      std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
  }
*/

 // Get all properties (will depend on the input file) and list them.
  std::vector<std::string> detected;
  std::string spaces4 = "    ";
  DAG->detect_available_props(detected);
  DAG->parse_properties (detected);
  std::cout << detected.size() << " metadata properties detected:" << std::endl;
  for( std::vector<std::string>::iterator kwdp = detected.begin();
          kwdp != detected.end(); ++kwdp )
  {
       std::string keyword = *kwdp;
       std::cout << spaces4 << keyword << std::endl;
  }
// Answer for test.h5m:
///////////////////////

// 5 metadata properties detected:
//    graveyard
//    mat
//    rho
//    surf.flux
//    tally
/////////////////

  // Loop through 3d entities.  In model_complete.h5m there are 90 vols
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "___entity_by_id___" << std::endl;
  std::cout << "\tnum_vols is " << num_vols << std::endl;
  std::cout << "Graveyard list: " << std::flush;
  EntityHandle entity = NULL;
  std::string gystr = "graveyard";
  for (unsigned i = 0; i<num_vols; i++)
  {
      entity = DAG->entity_by_id(3, i);
      std::string propstring;
  //  Test code: no longer needed
  //  code = DAG->measure_volume(entity, volume_measure);
  //  std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
      // make a property string; after obb_analysis.cpp:make_property_string(dag, eh, vector)
      if (DAG->has_prop(entity, gystr))
      {
         std::vector< std::string> vals;
         code = DAG->prop_values (entity, gystr, vals);
         propstring += gystr; 
         if (vals.size() == 1)
         {
            propstring += "=";
            propstring += vals[0];
         }
         else if (vals.size() > 1)
         {
            propstring += "=[";
            for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
            {
                propstring += *i;
                propstring += ", ";
            }
           // replace the last trailing comma with a close bracket
           propstring[propstring.length()-1] = ']';
         } // end if vals.size() > 1
         propstring += ", "; 
      } // end if DAG->has_prop
      std::ostream& out = std::cout;
      out << "\nVolume " << i  << " " << std::flush;
      if (propstring.length() )
      {
         propstring.resize (propstring.length() - 2); // drop trailing comma
      // end make_property_string fragment
//         out << "\nVolume " << i  << " " << std::flush;
         if (DAG->is_implicit_complement(entity) ) out << "Properties: " << propstring << std::endl;
         out << "(" << propstring << ")";
      }
      else
      {
         out << "  this volume is not a graveyard ";
      }
  }  // end this volume
  std::cout << std::endl;
}

ErrorCode readGraveyard()
{
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols is " << num_vols << std::endl;
  std::cout << "Graveyard list: " << std::endl;
  EntityHandle entity = NULL;
  // std::string gystr = "graveyard";
  ErrorCode ret;
  std::map<std::string, std::string> keyword_synonyms;

  std::vector< std::string > keywords;
  ret = DAG->detect_available_props( keywords );
  // If you don't do this props from make_property_string will be empty.
  ret = DAG->parse_properties( keywords);

   std::cout << keywords.size() << " metadata properties detected:" << std::endl;
   for( std::vector<std::string>::iterator i = keywords.begin(); i != keywords.end(); ++i )
   {
       std::cout << "    " << (*i) << std::endl;
   }
  // Loop through 3d entities.  In model_complete.h5m there are 90 vols
  for (unsigned i = 0; i<num_vols; i++)
  {
      entity = DAG->entity_by_index(3, i);
      std::string props = make_property_string(*DAG, entity, keywords);
      if (props.length()) std::cout << "Parsed props: " << props << std::endl; 
      if (DAG->has_prop(entity, "graveyard"))
      {
	 std::cout << "        graveyard" << std::endl;
      }
  }
}


static std::string make_property_string( DagMC& dag, EntityHandle eh, std::vector<std::string> &properties )
{
  ErrorCode ret;
  std::string propstring;
  for( std::vector<std::string>::iterator p = properties.begin();
    p != properties.end(); ++p )
  {
    if( dag.has_prop( eh, *p ) ){
      std::vector< std::string> vals;
      ret = dag.prop_values( eh, *p, vals );
      propstring += *p;
      if( vals.size() == 1 ){
        propstring += "=";
        propstring += vals[0];
      }
      else if( vals.size() > 1 ){
        // this property has multiple values, list within brackets
        propstring += "=[";
        for( std::vector<std::string>::iterator i = vals.begin();
             i != vals.end(); ++i )
        {
            propstring += *i;
            propstring += ",";
        }
        // replace the last trailing comma with a close braket
        propstring[ propstring.length()-1 ] = ']';
      }
      propstring += ", ";
    }
  }
  if( propstring.length() ){
    propstring.resize( propstring.length() - 2 ); // drop trailing comma
  }
  return propstring;
}

