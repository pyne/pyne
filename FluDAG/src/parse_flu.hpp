#ifndef FLUDAG_SRC_PARSE_FLU_HPP
#define FLUDAG_SRC_PARSE_FLU_HPP

#include <string>
#include <list>
#include <map>
#include <cassert>

#include "moab/Interface.hpp"
#include "moab/Types.hpp"

using namespace moab;

  /*
   * Convenience function to open a file for writing
   */
  void writeToFileNamed(std::ostringstream& oss, std::string index_id_filename);
  void addToIDIndexFile(Interface* mbi, int i, std::ostringstream& idstr);
  
  //---------------------------------------------------------------------------//
  // basic_score
  //---------------------------------------------------------------------------//
  /// This function encapsulates very similar records for USRTRACK and USRCOLL
/*
  void basic_score(std::ostringstream& ostr, 
                 std::string score_words, std::string particle_words, 
                 float fVol, float fortran_unit, double measurement,
                 std::string name);
*/
  void basic_score(Interface* mbi, std::ostringstream& ostr, 
                 std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, double measurement,
                 std::string name);
  //---------------------------------------------------------------------------//
  // two_vol_score
  //---------------------------------------------------------------------------//
  /// This function encapsulates very similar records for USRBDX and USRYIELD
  void two_vol_score(std::ostringstream& ostr, std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, std::string name);

  /* Convenience function for tokenizing.  
   * Ref: http://stackoverflow.com/questions/10051679/c-tokenize-string 
   */
  std::vector<std::string> inline StringSplit(const std::string &source, 
                                              const char *delimiter = " ", 
					      bool keepEmpty = false)
  {
        std::vector<std::string> results;
        size_t prev = 0;
	size_t next = 0;		

        // Reminder: string::npos ==> what 'next' will be if nothing found
	while ( (next = source.find_first_of(delimiter, prev)) != std::string::npos)
	{
 	    if (keepEmpty || (next - prev != 0))
	    {
		results.push_back(source.substr(prev, next-prev));
	    }		
	    prev = next + 1;
	}
	
	if (prev < source.size())
	{
 	    results.push_back(source.substr(prev));
	}
	
	return results;
  }

  /*
   * Process a material assignment for the ith volume
   */
  void process_Mi(Interface* mbi, std::ostringstream& ostr, 
                  EntityHandle entity, std::list<std::string> &matList, unsigned i);
  /*
   * Process all scoring requests for the ith volume.
   * 	- create a scoring-request vs volume-id multimap (numerous keys per value)
   *    - create a unique-scoring-request vs fortran-unit-number map.  
   *      Note 1: the unit number will be negative so that unformatted data will be produced 
   *      Note 2: the unit number is translated from an int to a floating point in the card
   */
  void process_Si(Interface* mbi, EntityHandle entity, unsigned i);

  void processUniqueMaterials(std::ostringstream& ostr, std::list<std::string> uniqueList,std::string header);

  void sdum_endline(std::ostringstream& ostr, std::string score, bool isBasic=false);
  std::string get_score_particle_mapname(std::string score_name, std::string particle_name);
  float get_score_particle_unit(std::vector<std::string> score_words);
  double measurementOfVol(Interface* mbi, int iVol);
  /*
   * Prepare a descriptive string that creates the properties of the volume whose index is index
   */
  std::string mat_property_string (moab::Interface* mbi, int index, std::vector<std::string> &properties);

#endif  // FLUDAG_SRC_PARSE_FLU
