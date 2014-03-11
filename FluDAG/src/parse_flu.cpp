// FLUDag/src/parse_flu.cpp

//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/fluka_funcs.cpp
 * \author Julie Zachman 
 * \date   Mon Mar 22 2013 
 * \brief  Functions called by fluka
 * \note   After mcnp_funcs
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//

#include "parse_flu.hpp"
#include "moab_utils.hpp"
#include "UnitNumberManager.hpp"
#include "chkerr.hpp"

#include "moab/Core.hpp"

#include <iomanip>
#include <sstream>
#include <set>
#include <cstring>
#include <string>

// #ifdef CUBIT_LIBS_PRESENT
// #include <fenv.h>
// #endif

// globals

#include <fstream>
// #include <numeric>

#define DEBUG 1
/* These 37 strings are predefined FLUKA materials. Any ASSIGNMAt of unique 
 * materials not on this list requires a MATERIAL card. */
std::string flukaMatStrings[] = {"BLCKHOLE", "VACUUM", "HYDROGEN",
"HELIUM", "BERYLLIU", "CARBON", "NITROGEN", "OXYGEN", "MAGNESIU",      
"ALUMINUM", "IRON", "COPPER", "SILVER", "SILICON", "GOLD", "MERCURY",  
"LEAD", "TANTALUM", "SODIUM", "ARGON", "CALCIUM", "TIN", "TUNGSTEN",   
"TITANIUM", "NICKEL", "WATER", "POLYSTYR", "PLASCINT", "PMMA",         
"BONECOMP", "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON",  
"POLYETHY", "AIR"};

int NUM_FLUKA_MATS = 37;

/* Create a set out of the hardcoded string array. */
std::set<std::string> FLUKA_mat_set(flukaMatStrings, flukaMatStrings+NUM_FLUKA_MATS); 

/* Maximum character-length of a cubit-named material property */
int MAX_MATERIAL_NAME_SIZE = 32;

std::multimap<std::string, unsigned int> scoring_vol_map;

MCRecordInfo       recordInfo = MCRecordInfo();
UnitNumberManager unit_no_mgr = UnitNumberManager();

// Empty synonym map for MCRecordInfo::parse_properties()
const std::map<std::string, std::string> MCRecordInfo::no_synonyms;

bool debug = false; //true ;

//---------------------------------------------------------------------------//
// fludagwrite_assignma
//---------------------------------------------------------------------------//
/// Called from main() 
//  This function 
// 	o writes out a simple numerical material assignment to the named argument file
//  Example usage:  parseMain dagmc.html
//  Outputs
//           mat.inp  contains MATERIAL and ASSIGNMAt records for the input geometry.
//                    The MATERIAL is gotten by parsing the Cubit volume name on underscores.  
//                    The string after "M_" is considered to be the material for that volume.
//                    There are no MATERIAL cards for the materials in the FLUKA_mat_set list
//                    For the remaining materials, there is one MATERIAL card apiece (no dups)
//                    User-named (not predefined) materials are TRUNCATED to 8 chars.
//                    User-named material id's start at 25 and increment by 1 for each MATERIAL card
//           index-id.txt  Map of FluDAG volume index vs Cubit volume ids, for info only.
//  The geometry information is contained in the Interface* object, mbi, which holds the meshed h5m
//  geometry extracted during a previous step.
//  The name of the (currently hardcoded) output file is "mat.inp"
//  The graveyard is assumed to be the last region.
void fludagwrite_assignma(Interface* mbi, std::string filename_to_write)  // file with cell/surface cards
{
  int num_vols = recordInfo.num_entities(3);
  // std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols is " << num_vols << std::endl;
  ErrorCode ret;
  EntityHandle entity = 0;
  int id;

  std::vector< std::string > keywords;
  ret = recordInfo.detect_available_props( mbi, keywords );

  // parse data from geometry so that property can be found
  ret = recordInfo.parse_properties( mbi, keywords );
  if (MB_SUCCESS != ret) 
  {
    std::cerr << "Failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // jcz debug: lists DEN, M, NEUTRON, S, USRTRACK
  std::vector< std::string >::iterator vit;
  for (vit=keywords.begin(); vit!=keywords.end(); ++vit)
  {
    std::cout << *vit << std::endl;
  }
  // Preprocessing loop:  make a string, "props",  for each vol entity
  // This loop could be removed if the user doesn't care to see terminal output
  std::cout << "Property list: " << std::endl;
  for (unsigned i=1; i<=num_vols; i++)
  {
     
      std::string props = mat_property_string(mbi, i, keywords);
      id = recordInfo.id_by_index(mbi, 3, i);
      if (props.length()) 
      {
         std::cout << "Vol " << i << ", id=" << id << ": parsed props: " << props << std::endl; 
      }
      else
      {
         std::cout << "Vol " << i << ", id=" << id << " has no props: " <<  std::endl; 
      }
  }

  std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";

  // Open an outputstring for mat.inp
  std::ostringstream ostr;
  std::ostringstream graveyard_str;
  std::ostringstream impl_compl_str;

  // open an outputstring for the M_  (ASSIGNMAt) portions
  std::ostringstream A_filestr;

  // Open an outputstring for index-id table and put a header in it
  std::ostringstream idstr;
  idstr << std::setw(5) <<  "Index" ;
  idstr << std::setw(5) <<  "   Id" << std::endl;

  // Prepare a list to contain unique materials not in Fluka's list
  std::list<std::string> uniqueMatList;

  // Loop through 3d entities (vols).  
  std::vector<std::string> vals;
  std::string material_trunc;
  char buffer[MAX_MATERIAL_NAME_SIZE];
  for (unsigned int i=1; i<=num_vols ; i++)
  {  
      vals.clear();
      entity = recordInfo.entity_by_index(3, i);

      // Create the id-index string for this vol
      addToIDIndexFile(mbi, i, idstr);

      // Create the mat.inp string for this vol
      if (recordInfo.has_prop(mbi, entity, "graveyard"))
      {
	 graveyard_str << std::setw(10) << std::left  << "ASSIGNMAt";
	 graveyard_str << std::setw(10) << std::right << "BLCKHOLE";
	 graveyard_str << std::setw(10) << std::right << i << std::endl;
      }
      if (recordInfo.has_prop(mbi, entity, "M"))
      {
         recordInfo.prop_values(mbi, entity, "M", vals);

         process_Mi(mbi, A_filestr, entity, uniqueMatList, i);
      } // end processing of "M_" property
      if (recordInfo.has_prop(mbi, entity, "S"))
      {
         process_Si(mbi, entity, i);
      } // end processing of "S_" property
  }  // end of volume processing loop

  std::ostringstream S_filestr;
  std::ostringstream r_filestr;
  std::ostringstream ut_filestr;
  std::ostringstream uc_filestr;
  std::ostringstream ub_filestr;
  std::ostringstream uy_filestr;

  // Print out the scoring by volume map
  std::cout << "All scoring.particle and volumes" << std::endl;
  unsigned int counter = 0;
  std::map<std::string, unsigned int>::iterator uit;
  
  // Go through the map that was created while going through the volumes
  for (uit = scoring_vol_map.begin(); uit != scoring_vol_map.end(); ++uit)
  {
     counter++;
     // std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
     std::cout << counter << ". " << uit->first << " => " << uit->second << std::endl;

     // Get the volume id of the current volume, whose scoring info we are pulling out
     int iVol = uit->second;

     // Create a vector of '.'-delimited strings from the original "key" part of the 
     // group name (now the "value" part of the map).  
     // The vector may be of size 1, 2, or 3, depending on the score type;
     // No more than the first two values are used in the unit map.
     std::vector<std::string> score_words = StringSplit(uit->first,".");
     std::string score_name;
     char particle[3];
     
     if (score_words.size() == 0)
     {
        std::cout << "The score group for volume " << iVol << " is empty." << std::endl;
        continue;
     }
     if (score_words.size() > 0)
     {
        score_name = score_words[0];
     }
     if (score_words.size() > 1)
     {

        std::size_t len = score_words[1].copy(particle,2);
        particle[len] = '\0';
     }
    
     char strDetName[10];
     float fortran_unit;
     double measurement;

     // The RESNUCLEI section:  use to_string when c11 comes
     if (score_name.compare("RESNUCLEI") == 0)
     {
        if (r_filestr.tellp() == 0)
        {
           r_filestr << "* RESNUCLEI scoring requests." << std::endl;
           r_filestr << header << std::endl;
        }
        // jcz ToDo;  Error return.  How serious an error is -1?  Is returning warranted?
        measurement = measurementOfVol(mbi, iVol);
        fortran_unit = get_score_particle_unit(score_words);
        sprintf (strDetName, "%s%d","RES_", iVol);

        r_filestr << std::setw(10) << std::left << "RESNUCLEI";
        r_filestr << std::setw(20) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
	r_filestr << std::setw(30) << std::right << std::fixed  << std::setprecision(1) << (float)iVol;
        r_filestr << std::setw(9)  << std::right << measurement << " ";
        r_filestr << std::setw(10) << std::left <<  strDetName;

        r_filestr << std::endl;
     }
     // The USRTRACK and USRCOLL section 
     else if (score_name.compare("USRTRACK") == 0 || score_name.compare("USRCOLL") == 0)
     {
        // These tallies need both a score and a particle name
	if (score_words.size() >= 2)
        {
           fortran_unit = get_score_particle_unit(score_words);
           if (score_name.compare("USRTRACK") == 0)
           {
              if (ut_filestr.tellp() == 0)
              {
                 ut_filestr << "* USRTRACK scoring requests." << std::endl;
                 ut_filestr << header << std::endl;
              }
              sprintf (strDetName, "%s_%s_%d","TR",  particle, iVol);
              basic_score(mbi, ut_filestr, score_words, iVol, fortran_unit, measurement, strDetName);
           }
 	   else  // USRCOLL case
           {
              if (uc_filestr.tellp() == 0)
              {
                 uc_filestr << "* USRCOLL scoring requests." << std::endl;
                 uc_filestr << header << std::endl;
              }
              sprintf (strDetName, "%s_%s_%d","CO",  particle, iVol);
              basic_score(mbi, uc_filestr, score_words, iVol, fortran_unit, measurement, strDetName);
           }
        }
        else   // Error:  no particle was added to the group name
        {
           std::cerr << "Error: the " << score_name << " score does not include a particle. " 
                     << "Please label the group with particle" << std::endl;
        }
     } 
     // The USRBDX section
     else if (score_name.compare("USRBDX") == 0 || score_name.compare("USRYIELD") == 0)
     {
        // These tallies need a score, particle name, and a FROM volume
        if (score_words.size() >= 3)
        {
           fortran_unit = get_score_particle_unit(score_words);
           if (score_name.compare("USRBDX") == 0)     // USRBDX case 
           {
              if (ub_filestr.tellp() == 0)
              {
                 ub_filestr << "* USRBDX scoring requests." << std::endl;
                 ub_filestr << header << std::endl;
              }
              two_vol_score(ub_filestr, score_words, iVol, fortran_unit, "BDX");
           }
           else   // USRYIELD case
           {
              if (uy_filestr.tellp() == 0)
              {
                 uy_filestr << "* USRYIELD scoring requests." << std::endl;
                 uy_filestr << header << std::endl;
              }
              two_vol_score(uy_filestr, score_words, iVol, fortran_unit, "YIELD");
           }
        }
        else   // Error:  missing particle or volume
        {
           std::cerr << "Error: the " << score_name << " score is missing a particle or a volume. " 
                     << "Please label the group with both particle and FROM volume" << std::endl;
        }
     } 
  }

  // Collect all the scoring records into one stream
  S_filestr <<  r_filestr.str();
  S_filestr << ut_filestr.str();
  S_filestr << uc_filestr.str();
  S_filestr << ub_filestr.str();
  S_filestr << uy_filestr.str();
  // Optional: send all the scores to the screen
  std::cout <<  S_filestr.str();

  // Finish the ostr with the implicit complement card
  std::string implicit_comp_comment = "* The next volume is the implicit complement";
  impl_compl_str << implicit_comp_comment << std::endl;
  impl_compl_str << std::setw(10) << std::left  << "ASSIGNMAt";
  impl_compl_str << std::setw(10) << std::right << "VACUUM";
  impl_compl_str << std::setw(10) << std::right << num_vols << std::endl;

  // Prepare the MATERIAL cards as a string stream using a list
  // of materials that has no duplicates
  uniqueMatList.sort();
  uniqueMatList.unique();
  std::ostringstream MAT_filestr;
  processUniqueMaterials(MAT_filestr, uniqueMatList, header);
  // Print the final list
  if (debug)
  {
     std::list<std::string>::iterator it; 
     for (it=uniqueMatList.begin(); it!=uniqueMatList.end(); ++it)
     {
        std::cout << *it << std::endl;
     }
  
     // Show the output string just created
     std::cout << ostr.str();
  }

  // Prepare a file of the given name; put a header in it;  This file contains the FLUKA input cards 
  // (records) the user has requested via geometry tags.
  std::ofstream records_filestr( filename_to_write.c_str());
  records_filestr << header << std::endl;

  // Put all the filestr parts together
  records_filestr << MAT_filestr.str();    // The material list is created separately
  records_filestr << graveyard_str.str();  // the graveyard
  records_filestr << A_filestr.str();      // ASIGNMAt statements
  records_filestr << impl_compl_str.str(); // implicit complement
  records_filestr << S_filestr.str();      // Detector tallies
  records_filestr.close();

  std::cout << "Writing input file = " << filename_to_write << std::endl; 

  // Prepare an output file named "index_id.txt" for idstr
   writeToFileNamed(idstr, "index_id.txt");
}
// End fludagwrite_assignma

//---------------------------------------------------------------------------//
// measurementOfVol
//---------------------------------------------------------------------------//
// Convenience function to call dag measure_volume(..)
double measurementOfVol(moab::Interface* mbi, int iVol)
{
   ErrorCode ret;

   // Calculate some things we'll need that are based on the volume id 
   EntityHandle handle = recordInfo.entity_by_id(mbi, 3, iVol);

   double measurement;
   ret = recordInfo.measure_volume(mbi, handle, measurement);
   if (MB_SUCCESS != ret) 
   {
      std::cerr << "Failed to get the measured volume of region " <<  iVol <<  std::endl;
      measurement = -1.0;
   }
   return measurement; 
}
//---------------------------------------------------------------------------//
// usrtrackRecord
//---------------------------------------------------------------------------//
// Specialized record-preparer for USRTRACK defaults
/*
void usrtrackRecord(std::ostringstream& ostr, 
                 std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, std::string score_prefix)
{
     // Prepare the detector name to go at the end of the line
     char strDetName[10];
     std::string subname (score_prefix + "_" + score_words[2] + "_");
     char *cstr = new char [subname.length() + 1];
     std::strcpy (cstr, subname.c_str());
     sprintf (strDetName, "%s%d", cstr, iVol);

     // We are guaranteed there are three values in score_words
     ostr << std::setw(10) << std::left << score_words[0];           
     ostr << std::setw(10) << std::right << "-1.0";
     ostr << std::setw(10) << std::right << score_words[1];           
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << score_words[2];
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;

     ostr << std::setw(10) << std::right << " ";
     ostr << std::setw(10) << std::left <<  strDetName << std::endl;

     sdum_endline(ostr, score_words[0]);
     if (score_words.size() > 3)  // hmm, there is a particle piece, but also more
     {
           std::cerr << "Error:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
*/
//---------------------------------------------------------------------------//
// two_vol_score
//---------------------------------------------------------------------------//
// Some records are very similar, differing only in name
void two_vol_score(std::ostringstream& ostr, 
                 std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, std::string score_prefix)
{
     // Prepare the detector name to go at the end of the line
     char strDetName[10];
     std::string subname (score_prefix + "_" + score_words[2] + "_");
     char *cstr = new char [subname.length() + 1];
     std::strcpy (cstr, subname.c_str());
     sprintf (strDetName, "%s%d", cstr, iVol);

     // We are guaranteed there are three values in score_words
     ostr << std::setw(10) << std::left << score_words[0];           
     ostr << std::setw(10) << std::right << "-1.0";
     ostr << std::setw(10) << std::right << score_words[1];           
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << score_words[2];
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;

     ostr << std::setw(10) << std::right << " ";
     ostr << std::setw(10) << std::left <<  strDetName << std::endl;

     sdum_endline(ostr, score_words[0]);
     if (score_words.size() > 3)  // hmm, there is a particle piece, but also more
     {
           std::cerr << "Error:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
//---------------------------------------------------------------------------//
// basic_score
//---------------------------------------------------------------------------//
// Some records are very similar, differing only in name
// This handls USRCOLL and USRTRACK score requests
void basic_score(moab::Interface* mbi, std::ostringstream& ostr, 
                 std::vector <std::string> score_words, 
                 int iVol, float fortran_unit, double measurement,
                 std::string name)
{
     measurement = measurementOfVol(mbi, iVol);
     if (score_words.size() >= 2)  // all is good
     {
           ostr << std::setw(10) << std::left << score_words[0];           
           ostr << std::setw(10) << std::right << "-1.0";
           ostr << std::setw(10) << std::right << score_words[1];           
           ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
	   ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;
           ostr << std::setw(10) << std::right << measurement << " ";
           // Default number of energy bins
           ostr << std::setw(8)  << std::right << "1.";

           // ostr << std::setw(10) << std::left <<  name << std::endl;;
           ostr << " " << name << std::endl;;
	   sdum_endline(ostr, score_words[0], true);
     }
     if (score_words.size() > 2)  // hmm, there is a particle word, but also more
     {
           std::cerr << "Warning:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
//---------------------------------------------------------------------------//
// sdum_endline
//---------------------------------------------------------------------------//
// Create a standard continuation line, putting the '&' on the 71st space
// The line is tacked on to the output stream
void sdum_endline(std::ostringstream& ostr, std::string score, bool isBasic)
{
     ostr << std::setw(10) << std::left << score;           
     if (isBasic)  // This path is for using the defaults common to USRTRACK and USRYIELD
     {
        // Default maximum energy
        ostr << std::setw(10) << std::right << "1E20";           
        // Default minimum energy
        ostr << std::setw(10) << std::right << "1.0";           
        ostr << std::setw(41) << std::right << "&";
     }
     else
     {
        ostr << std::setw(61) << std::right << "&";
     }
     ostr << std::endl;
}

//---------------------------------------------------------------------------//
// process_Si
//---------------------------------------------------------------------------//
// Process group names that have S as the key and track.[p'le] as the value
// Examples:  USRTRACK.PROTON, RESNUCLEI, USRCOLL.NEUTRON
// A multimap is used to store every score request for every volume,
void process_Si(Interface* mbi, EntityHandle entity, unsigned int vol_id)
{
    ErrorCode ret;
    std::vector<std::string> vals;

    // We only get here if has_prop(... "S" ...) is true
    ret = recordInfo.prop_values(mbi, entity, "S", vals);
    if (MB_SUCCESS != ret) 
    {
       std::cerr << "Failed to get S_ properties" <<  std::endl;
       return;
    }
    for (int i=0; i<vals.size(); i++) 
    { 
        scoring_vol_map.insert(std::pair<std::string, unsigned int>(vals[i], vol_id)); 
    }
}

//---------------------------------------------------------------------------//
// get_score_particle_mapname
//---------------------------------------------------------------------------//
/// Convert the first two values of the group name (in the case there are two) to a 
//  standard name for the fortran unit number getter
std::string get_score_particle_mapname(std::string score_name, std::string particle_name)
{
    std::string keyword (score_name + "." + particle_name);
    char *cstr = new char [keyword.length() + 1];
    std::strcpy (cstr, keyword.c_str());
    return cstr;
}

//---------------------------------------------------------------------------//
// get_score_particle_unit
//---------------------------------------------------------------------------//
/// Parse the vector of strings from the scoring group name values and determine
//  the correct fortran unit number for them.  
//  There can be 1, 2, or more words in teh score_words vectore. 
//  If 1, it is used.  If 2, both are used.  If more than 2, only the first two are used.
// This function relies on the UnitNumberManager class, whose key method returns
// an int, however we always need a float unit number for the fluka cards, so
// a float is returned. 
//  If the vector is empty, this function returns -1
float get_score_particle_unit(std::vector<std::string> score_words)
{
    std::string mapname;
    if (score_words.size() == 1)
    {
        mapname = score_words[0];
    }
    else if (score_words.size() > 1)
    {
        mapname = get_score_particle_mapname(score_words[0], score_words[1]);
    }
    return (float)unit_no_mgr.getUnitNumber(mapname);
}
//---------------------------------------------------------------------------//
// processUniqueMaterials
//---------------------------------------------------------------------------//
// Convenience method to create MATERIAL cards
void processUniqueMaterials(std::ostringstream& ostr, 
                            std::list<std::string> uniqueList, 
                            std::string header)
{
  // Prepare the MATERIAL cards as a string stream
  if (uniqueList.size() != 0)
  {
     int matID = 25;
     std::list<std::string>::iterator it; 
     for (it=uniqueList.begin(); it!=uniqueList.end(); ++it)
     {
        ostr << std::setw(10) << std::left << "MATERIAL";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << ++matID;
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::left << *it << std::endl;
     }
  }
  if (uniqueList.size() !=0)
  {
     ostr << header << std::endl;
  }
  return;
}
//---------------------------------------------------------------------------//
// process_MI
//---------------------------------------------------------------------------//
// Add a record to ostr of the form "ASSIGNMA ....."
void process_Mi(Interface* mbi, std::ostringstream& ostr, 
                EntityHandle entity, 
                std::list<std::string> &matList, unsigned i)
{
    ErrorCode ret;
    std::vector<std::string> vals;
    char buffer[MAX_MATERIAL_NAME_SIZE];
    std::string material_trunc;

    ret = recordInfo.prop_values(mbi, entity, "M", vals);
    if (MB_SUCCESS != ret) 
    {
       std::cerr << "Failed to get M_ properties" <<  std::endl;
       return;
    }

    if (vals.size() >= 1)
    {
       // Make a copy of string in vals[0]; full string needs to be compared to
       // FLUKA materials list; copy is for potential truncation
       std::strcpy(buffer, vals[0].c_str());
       material_trunc = std::string(buffer);
            
       if (vals[0].size() > 8)
       {
           material_trunc.resize(8);
       }

       if (FLUKA_mat_set.find(vals[0]) == FLUKA_mat_set.end())
       {
          // current material is not in the pre-existing FLUKA material list
          matList.push_back(material_trunc); 
          std::cout << "Adding material " << material_trunc << " to the MATERIAL card list" << std::endl;
       }
    }
    else
    {
         material_trunc = "moreThanOne";
    }
     ostr << std::setw(10) << std::left  << "ASSIGNMAt";
     ostr << std::setw(10) << std::right << material_trunc;
     ostr << std::setw(10) << std::right << i << std::endl;
}
//---------------------------------------------------------------------------//
// addToIDIndexMap(int i, 
//---------------------------------------------------------------------------//
// Convenience method to connect the geometry id to the ith volume 
void addToIDIndexFile(moab::Interface* mbi, int i, std::ostringstream &idstr)
{
      idstr << std::setw(5) << std::right << i;
      idstr << std::setw(5) << std::right << recordInfo.id_by_index(mbi, 3,i) << std::endl;
}
//---------------------------------------------------------------------------//
// writeStringToFile
//---------------------------------------------------------------------------//
// Convenience method to write a prepared stream to a file
void writeToFileNamed(std::ostringstream& file_contents, std::string filename)
{
     std::ofstream file_stream(filename.c_str());
     file_stream << file_contents.str();
     file_stream.close(); 
}

//---------------------------------------------------------------------------//
// mat_property_string
//---------------------------------------------------------------------------//
// For a given volume, find the values of all properties named "MAT".   
// Create a string with these properites
// Modified from make_property_string
// This function helps with debugging, but is not germane to writing cards.
std::string mat_property_string (moab::Interface* mbi, int index, std::vector<std::string> &properties)
{
  ErrorCode ret;
  std::string propstring;
  EntityHandle entity = recordInfo.entity_by_index(3,index);
  int id = recordInfo.id_by_index(mbi, 3, index);
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( recordInfo.has_prop(mbi, entity, *p) )
     {
        std::vector<std::string> vals;
        ret = recordInfo.prop_values(mbi, entity, *p, vals);
        CHECKERR(mbi, ret);
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

// } // namespace moab
