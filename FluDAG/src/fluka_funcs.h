#ifndef DAGMC_MCNP_IFACE_H
#define DAGMC_MCNP_IFACE_H

#include <string>
#include <vector>

#include "moab/Types.hpp"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"

#ifdef __cplusplus
extern "C" {
#endif


/* 
 * This signature is intended to be called from a c++ function.  
 * It assumes the caller does not need to know the dagmc_version or the moab_version.
 * The original dagmcinit_, which is to be called by a Fortran method, has been reworked 
 * to set up variables and then call this version.
 * 
 * 14 Feb 2013  jcz added boolean running_with_fluka, which prepends "../" to the file name
 */
  void cpp_dagmcinit(const std::string infile, 
                int parallel_file_mode, // parallel read mode
                int max_pbl);

  void slow_check(double pos[3], const double dir[3], int &oldReg);
  // check we are where we say we are
  MBEntityHandle check_reg(MBEntityHandle volume, double point[3], double dir[3]); 
  /*
   * Convenience function to open a file for writing
   */
  void writeToFileNamed(std::ostringstream& oss, std::string index_id_filename);
  void addToIDIndexFile(int i, std::ostringstream& idstr);
  /*
   * Process a material assignment for the ith volume
   */
  void process_Mi(std::ostringstream& ostr, MBEntityHandle entity, std::list<std::string> &matList, unsigned i);
  /*
   * Process all scoring requests for the ith volume.
   * 	- create a scoring-request vs volume-id multimap (numerous keys per value)
   *    - create a unique-scoring-request vs fortran-unit-number map.  
   *      Note 1: the unit number will be negative so that unformatted data will be produced 
   *      Note 2: the unit number is translated from an int to a floating point in the card
   */
  void process_Si(MBEntityHandle entity, unsigned i);
  void processUniqueMaterials(std::ostringstream& ostr, std::list<std::string> uniqueList,std::string header);
  void sdum_endline(std::ostringstream& ostr, std::string score, bool isBasic=false);
  std::string get_score_particle_mapname(std::string score_name, std::string particle_name);
  float get_score_particle_unit(std::vector<std::string> score_words);
  double measurementOfVol(int iVol);
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
  void basic_score(std::ostringstream& ostr, 
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

  /* get the sense of a region with respect to the global next_surf,
   * which is set by a call to rayfire
  */
  int getSense(int region);
  /*
   * Prepare a descriptive string that creates the properties of the volume whose index is index
   */
  std::string mat_property_string (int index, std::vector<std::string> &properties);
  /*
   * Write records to a *.inp file for each volume tagged with material info and scoring requests.
   */
  void fludagwrite_assignma(std::string matfile);  

  void dagmc_version_(double* dagmcVersion);
#ifdef __cplusplus
} // extern "C"
#endif

#endif /* DAGMC_MCNP_IFACE_H */
