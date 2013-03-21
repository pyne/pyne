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
#include "DagMC.hpp"
#include <fstream>
#include "moab/Types.hpp"
#include "moab/Core.hpp"
#include "DagWrapUtils.hh"

using namespace moab;

#define DAG DagMC::instance()

static std::vector<EntityHandle> null_list;

static std::string make_property_string( DagMC& , EntityHandle , std::vector<std::string> & );
ErrorCode get_prop_values (EntityHandle entity, std::string property, std::string& propstring, 
                           std::vector<EntityHandle> entityList);
std::string makeRegionName(int );
ErrorCode mcnp5_property_names( Interface* MBI );


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
 
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols in " << fileptr << " is " << num_vols << std::endl;

  ErrorCode code;
  // double volume_measure;

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


  std::cout << "Property list: " << std::flush;
  EntityHandle entity = NULL;
  std::string gystr = "graveyard";
  std::vector<EntityHandle> graveyardList;
  std::string matstr = "MAT";
  std::vector<EntityHandle> matVols;
  std::string usrtrackstr = "USRTRACK";
  std::vector<EntityHandle> usrtrackVols;

  std::string propstring;
  // Loop through vols looking for properties of interest.
  for (unsigned i = 0; i<num_vols; i++)
  {
      entity = DAG->entity_by_id(3, i);
      std::string propstring;
  //  code = DAG->measure_volume(entity, volume_measure);
  //  std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
      // Append to a property string for selected properties
      get_prop_values(entity, matstr,             propstring, matVols);
      get_prop_values(entity, usrtrackstr,        propstring, usrtrackVols);
      // get_prop_values(entity, CATEGORY_TAG_NAME,  propstring, null_list);
      // get_prop_values(entity, NAME_TAG_NAME,      propstring, null_list);
      get_prop_values(entity, gystr,              propstring, graveyardList);

      Core mbi;
      mcnp5_property_names( &mbi );

      std::ostream& out = std::cout;
      out << "\nVolume " << i  << " " << std::flush;
      if (propstring.length() )
      {
         propstring.resize (propstring.length() - 2); // drop trailing comma
         if (DAG->is_implicit_complement(entity) ) 
         {
            out << "Properties: " << propstring << std::endl;
         }
         out << "(" << propstring << ")";
      }
  }  // end this volume
  
  std::cout << std::endl;
}

//---------------------------------------------------------------------------//
// get_prop_values()  
//---------------------------------------------------------------------------//
/// For the given property, append its values to the string.
/// @param entity     - the volume to check
/// @param property   - the given property
/// @param propstring - string to append to
/// @param entityList - fills with those volumes that have the given property. 
///                     If the statically-defined empty null_list was passed, don't do anything.
ErrorCode get_prop_values (EntityHandle entity, std::string property, std::string& propstring, 
                           std::vector<EntityHandle> entityList)
{
  ErrorCode code;
  std::vector<std::string> vals;
  if (DAG->has_prop(entity, property))
  {
	
     if (&entityList != &null_list)
     {
        entityList.push_back(entity);
     }
     code = DAG->prop_values (entity, property, vals);
     propstring += property; 
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
  return code;
}



//---------------------------------------------------------------------------//
// get_signed_volume(..)
//---------------------------------------------------------------------------//
/// From dagmc_preproc
/// 1. fill category tag with tag names
ErrorCode mcnp5_property_names( Interface* MBI )
{
  ErrorCode rval;
  Tag category_tag;
  rval = MBI->tag_get_handle( CATEGORY_TAG_NAME, 32, MB_TYPE_OPAQUE, category_tag );
  if (MB_SUCCESS != rval)
    return rval;
  
  // group_val[]:  Set an array of pointers consisting of one ptr to string "Group"
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category+CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};

  // groups:  Range filled with entities that have tags for group.  I think
  Range groups;
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &category_tag,
                                           group_val, 1, groups);
  if (MB_SUCCESS != rval)
    return rval;

  Tag name_tag;
  rval = MBI->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE, name_tag );
  if (MB_SUCCESS != rval)
    return rval;
  // std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  // Run through every entity in "groups"
  for( Range::iterator i = groups.begin(); i != groups.end(); ++i )
  {
    EntityHandle grp = *i;
    const void* p;
    int ignored;
    // Get the name_tag data for the current group
    rval = MBI->tag_get_by_ptr( name_tag, &grp, 1, &p, &ignored );
    if( MB_SUCCESS != rval ) return rval;

    // string modname:  convert the name_tag data (in &p from previous call)
    const char* grpname = static_cast<const char*>(p);
    std::string modname(grpname);
    std::cout << modname << std::endl;
    /*
    size_t idx;
    if( modname.find("tally_") == 0 ){
        std::string arg = modname.substr(6);
        // replace up to 2 underscores
        int count = 0;
        while( count < 2 &&  (idx = arg.find_first_of("_")) != arg.npos )
        {
          count ++;
          arg[idx] = '.';
        }
        modname = modname.substr(0,6) + arg;
    }
    else if( (idx = modname.find("imp_")) != modname.npos ){
        modname[idx+3] = '.';
    }

    if( modname != grpname ){
      std::cout << "Group name " << grpname << " changed to " << modname << std::endl;
      p = static_cast<const void*>(modname.c_str());
      int length = NAME_TAG_SIZE;
      rval = MBI->tag_set_by_ptr( name_tag, &grp, 1, &p, &length);
      if( MB_SUCCESS != rval ) return rval;
    }
*/
  } // End iteration through groups
  return MB_SUCCESS;
}

// readProperties()  -- not called
//---------------------------------------------------------------------------//
// 
/*
ErrorCode readProperties()
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
  // If you don't call parse_properties, props from make_property_string will be empty.
  ret = DAG->parse_properties( keywords);

   std::cout << keywords.size() << " metadata properties detected:" << std::endl;
   // for divgroup.h5m, there are 6 metadata properties detected.  They are:
   // 1.0, FLUX, H20, MAT, STUFFING, USRTRACK
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
      // Optional check to see if the current volume has property "graveyard"
      if (DAG->has_prop(entity, "graveyard"))
      {
	 std::cout << "        graveyard" << std::endl;
      }
      std::string propValue;
      std::string prop = "MAT";
      if (DAG->has_prop(entity, "MAT"))
      {
         ret = prop_value(entity,  "MAT", propValue);
         std::cout << "MAT value is " << propValue << std::endl;
      }
  }
}
*/

//---------------------------------------------------------------------------/
// obbstat_write
//---------------------------------------------------------------------------//
/// This is copied and modified from obb_analysis.cpp in the moab/tools/dagmc dir. 
/// Changes:
/// - static DAG used instead of passed in
/// - std::ostream& out -> std::cout, moved to the body
/*
ErrorCode obbstat_write( std::vector<int> &volumes,
                         std::vector<std::string> &properties )
{

  std::ostream out = std::cout;
  ErrorCode ret = MB_SUCCESS;
  OrientedBoxTreeTool& obbtool = *dag.obb_tree();

  // can assume that volume numbers are valid.
  for( std::vector<int>::iterator i = volumes.begin(); i!=volumes.end(); ++i){
    EntityHandle vol_root;
    EntityHandle vol = dag.entity_by_id(3,*i);
    CHECKERR(dag,ret);

    if( vol == 0 ){
      std::cerr << "ERROR: volume " << *i << " has no entity." << std::endl;
      continue;
    }

    ret = dag.get_root( vol, vol_root );
    CHECKERR(dag,ret);

    out << "\nVolume " << *i << " " << std::flush;

    if( dag.is_implicit_complement(vol) ) out << "(implicit complement) ";
    out << std::endl;

    std::string propstring = make_property_string( dag, vol, properties );
    if( propstring.length() ) out << "Properties: " << propstring << std::endl;

    // get all surfaces in volume
    Range surfs;
    ret = dag.moab_instance()->get_child_meshsets( vol, surfs );
    CHECKERR(dag,ret);

    out << "   with " << surfs.size() << " surfaces" << std::endl;

    TriStats ts( dag.moab_instance(), &obbtool, vol_root );
    ret = obbtool.preorder_traverse( vol_root, ts );
    CHECKERR(dag,ret);
    ts.write_results( out );
 
if( verbose ){
      out << "Surface list: " << std::flush;
      for( Range::iterator j = surfs.begin(); j!=surfs.end(); ++j){
        out << dag.get_entity_id(*j);
        std::string props = make_property_string( dag, *j, properties );
        if( props.length() ) out << "(" << props << ")";
        if( j+1 != surfs.end() ) out << ",";
      }
      out << std::endl;
      ret = obbtool.stats( vol_root, out );
      CHECKERR(dag,ret);
    }

    out << "\n    ------------ " << std::endl;

  }

  return ret;
}
*/
//---------------------------------------------------------------------------/
// makePropertyString
//---------------------------------------------------------------------------//
/// 
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
      CHECKERR(dag,ret);
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
//////////////////////////////////////////////////////////////////////////////
//////////////////////////  Create Mat File Section //////////////////////////
//////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------/
// createFlukaMatFile()
//---------------------------------------------------------------------------//
/// 
void createFlukaMatFile() 
{
  
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG createFlukaMatFile()" << std::endl;
  std::cout << "================== FILEWR =================" << std::endl;
#endif 


  std::ofstream vos("test1.inp");
  PrintEntityRegionNames(vos);
  vos.close();

}

//---------------------------------------------------------------------------/
// PrintEntityRegionNames()
//---------------------------------------------------------------------------//
/// 
void PrintEntityRegionNames(std::ostream& os)
{
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG PrintEntityNames()" << std::endl;
#endif
  PrintHeader(os, "VOLUMES, and Names");
  std::string Vname;
  unsigned int numVol = DAG->num_entities(3);
  for(unsigned int l=0; l < numVol; l++) 
  {
    int iRegion = l+1;
    Vname = makeRegionName(iRegion);
    //Print index and region name in some fixed format
    writeRegionLine(os, iRegion, Vname);
  }
  int iRegion = numVol + 1;
  Vname = "BLACKHOL";
  writeRegionLine(os, iRegion, Vname);

#ifdef DAGGEOMETRY_DEBUG
  std::cout << "<== FluDAG PrintEntityRegionNames()" << std::endl;
#endif
}

//---------------------------------------------------------------------------/
// makeRegionName
//---------------------------------------------------------------------------//
/// 
std::string makeRegionName(int l)
{

  std::string VVname;
  std::string Vname;
  EntityHandle entity = NULL;
  entity = DAG->entity_by_id(3, l);
  // Don't add 1
  int iRegion = l;
  std::cout << iRegion << " index,name " << std::endl;
  char vname[8];
  sprintf(vname,"%-8u",iRegion);
  Vname.replace(0,8,vname);
  std::cout<<iRegion<<" vname" << vname <<" Vname " << Vname<<std::endl;
  unsigned int found=Vname.rfind(" ",7);
  // take out blanks
  while (found<Vname.size()-1)
  {
      found=Vname.rfind(" ",7);
      if (found<Vname.size()-1)
      {
	  std::string temp=Vname.substr(found+1);
	  Vname.replace(found, temp.size(), temp);
	  Vname.erase(Vname.size()-1);
	  std::cout << Vname << std::endl;
      } 
  }

    unsigned int nameSize=Vname.size();
    if (nameSize > 8 ) {VVname=Vname.substr(0,7);}
    else {VVname=Vname;}
    std::cout << "VVname "<<VVname<< std::endl;
    // ToDo:  Re-implement the maps if we need to check that the name is unique and there
    // isn't another way.
    // check that name is unique, if not add numbering
/*
    unsigned int matrep = 1;
    unsigned int ii=VVname.size()+3;
    unsigned int  newSize = ii < 8 ? ii : 8;
    bool old=true;
    char smatrep[3];
       while (old)
       {
	 old=false;
         for ( DRegionIterator i=dRegionNameMap.begin(); i!=dRegionNameMap.end();i++)
         {
           sprintf(smatrep,"%03u",matrep);
           Vname=(*i).second;
           if (Vname==VVname)
           {
   	      old=true;
   	      if (VVname.size()<=5)
              {
   	          VVname+=smatrep;
              }
   	      else
              {
   	          VVname.resize(newSize);
   	          VVname.replace(VVname.size()-3, 3,smatrep);
              }
	   matrep++;
           }
         }
       }
      std::cout<< "VVname "<< VVname<< std::endl;
*/
      return VVname;
}
///////			End makeRegionName
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------/
// writeRegionLine
//---------------------------------------------------------------------------//
/// Write one line of the *.inp file
/// Convenience function that is re-used
void writeRegionLine(std::ostream& os, int iRegion, std::string name)
{
    os.setf(std::ios::left, std::ios::adjustfield);
    os << setw10 << iRegion;
    // ToDo:  replace "entname" if possible
    // os << std::setw(20) << entity->GetName() << std::setw(20) << "";
    os << std::setw(20) << "entname"  << std::setw(20) << "";
    os << std::setw(5) << name << std::setw(5) << "";
    //If volume is a replica... print some more stuff
    os << std::endl;
}
///////			End writeRegionLine
/////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------/
// PrintHeader
//---------------------------------------------------------------------------//
/// 
std::ostream& PrintHeader(std::ostream& os, const char* title) 
{
  os << "*\n" << "*\n" << "*\n";
  os << "*********************  " << title << " *********************\n"
     << "*\n";
  os << "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7..."
     << std::endl;
  os << "*" << std::endl;

  return os;
}
// End PrintHeader
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
///////////////// End Create Mat File Section ////////////////////////
/////////////////////////////////////////////////////////////////////

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
