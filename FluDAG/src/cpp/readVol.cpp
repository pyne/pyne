<<<<<<< HEAD
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   storage/DAGMC/FluDAG/src/cpp/readVol.cpp
 * \author Julie Zachman 
 * \date   Fri Mar 8 2013 
 * \brief  Read properties from an h5m file
 * \note   Includes a main
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//

#include "fluka_funcs.h"
#include "dagmc_utils.hpp"

#include "moab/Interface.hpp"
=======
// Test code that calls DagMC's readVol
>>>>>>> svalinn/master
#include "DagMC.hpp"
#include <iostream>
#include <stdlib.h>
#include <string>
#include "moab/Types.hpp"

using namespace moab;

#define DAG DagMC::instance()

static std::string make_property_string( DagMC& , EntityHandle , std::vector<std::string> & );

<<<<<<< HEAD
//---------------------------------------------------------------------------//
// readVol
//---------------------------------------------------------------------------//
/// Highest level cal; reads number of volumes
// 
// sample metadata properties detected:
//    graveyard
//    mat
//    rho
//    surf.flux
//    tally
/////////////////
ErrorCode readVol(char *fileptr, const std::string* override_output_filename=NULL )
{
  // This code fragment is from TrackLengthMeshTally:write_results has no purpose other than to test 
  // setting the arg NULL, so that I can implement the same signature in both extensions of MeshTally
  // Case 1:  DONE no second arg to readVol:  Override = "myfilename"
  // Case 2:  DONE second arg to readVol is null:  Override = "myfilename"
  // Case 3:  DONE second arg to readVol is not null:  Override is same as second arg to readvol
  // std::string output_filename = "myfilename";
  // std::string override = override_output_filename ? *override_output_filename : output_filename;
  // std::cout << "Output_filename = " << output_filename << std::endl;
  // std::cout << "Override = " <<  override << std::endl;
  // Note:  see main() in this file for the call
  
=======
ErrorCode readVol(char *fileptr)
{
  // See if this works
>>>>>>> svalinn/master
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols in " << fileptr << " is " << num_vols << std::endl;

  ErrorCode code;
<<<<<<< HEAD
  // double volume_measure;

  // Get all properties (will depend on the input file) and list them.
=======
  double volume_measure;
  

 // Get all properties (will depend on the input file) and list them.
>>>>>>> svalinn/master
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
<<<<<<< HEAD

  // Loop through 3d entities.  
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "___entity_by_id___" << std::endl;

  std::cout << "Graveyard list: " << std::flush;
  EntityHandle entity = NULL;
  std::string gystr = "graveyard";
  std::vector<int> graveyardList;
=======
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
>>>>>>> svalinn/master
  for (unsigned i = 0; i<num_vols; i++)
  {
      entity = DAG->entity_by_id(3, i);
      std::string propstring;
<<<<<<< HEAD
=======
  //  Test code: no longer needed
>>>>>>> svalinn/master
  //  code = DAG->measure_volume(entity, volume_measure);
  //  std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
      // make a property string; after obb_analysis.cpp:make_property_string(dag, eh, vector)
      if (DAG->has_prop(entity, gystr))
      {
<<<<<<< HEAD
	 graveyardList.push_back(i);
=======
>>>>>>> svalinn/master
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
<<<<<<< HEAD
      else if (DAG->has_prop(entity, "MAT"))
      {
      }
=======
>>>>>>> svalinn/master
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
<<<<<<< HEAD
=======
      else
      {
         out << "  this volume is not a graveyard ";
      }
>>>>>>> svalinn/master
  }  // end this volume
  std::cout << std::endl;
}

<<<<<<< HEAD
//---------------------------------------------------------------------------//
// readGraveyard
//---------------------------------------------------------------------------//
// 
=======
>>>>>>> svalinn/master
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
<<<<<<< HEAD
  // If you don't call parse_properties, props from make_property_string will be empty.
=======
  // If you don't do this props from make_property_string will be empty.
>>>>>>> svalinn/master
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


<<<<<<< HEAD
//---------------------------------------------------------------------------//
// makePropertyString
//---------------------------------------------------------------------------//
/// 
=======
>>>>>>> svalinn/master
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
<<<<<<< HEAD
      CHECKERR(dag,ret);
=======
>>>>>>> svalinn/master
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
<<<<<<< HEAD
//---------------------------------------------------------------------------//
// main
//---------------------------------------------------------------------------//
// 
int main(int argc, char* argv[]) 
{
  ErrorCode code;
  if (argc < 2)
  {
     std::cerr << "Usage:  " << argv[0] << " filename.h5m" << std::endl;
     exit(0);
  }
  
 
  int max_pbl = 1;
  bool flukarun = false;
  char *fileptr = argv[1];
  cpp_dagmcinit(fileptr, 0, max_pbl,flukarun); 

  // Test code
  // std::string test;
  // if (argc > 2)
  // {
  //  test = std::string(argv[2]);
  // }
  // code = readVol(argv[1], &test);

  code = readVol(argv[1]);
  return 0;
}
=======

>>>>>>> svalinn/master
