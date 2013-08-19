#include "meshtal_funcs.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <map>
#include <cstring>

#include "TallyManager.hpp"

// For tally multipliers
void mcnp_weight_calculation( int* index, double* erg, double* wgt, 
                              double* dist, double* score_result )
{
    FMESH_FUNC(dagmc_mesh_score)( index, erg, wgt, dist, score_result );
}

/*
ToDo:  This has been moved from TrackLengthMeshTally to here:  it is mcnp-related
ToDo:  Get this working in meshtal_funcs
       Probably:  Converts the name read in from the input file 
                  to the MCNP index
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
*/


/********************************************************************
 * Initialization and setup functions
 ********************************************************************/ 

TallyManager tallyManager = TallyManager();

// static bool initialized = false;

/** 
 * Called at least once from fmesh_mod on program initialization;
 * in runtpe or MPI modes may be called multiple times.  
 */
/*
void dagmc_fmesh_initialize_( const int* mcnp_icl ){

  if( initialized ) return;

  //std::cerr << "Executed DAGMC fmesh initialize." << std::endl;

  current_mcnp_cell = mcnp_icl;

  initialized = true;
}
*/

/**
 * Convert the contents of an FC card to an fmesh_params multimap<string,string>
 * @param fc_content The FC card's comment content as a string
 * @param fmesh_params The output data and comments, as a multimap
 * @param fcid The tally ID of the FC card
 * @return true on success, or false if the input has serious enough formatting problems
 *         to make parameter parsing impossible.
 */
static void parse_fc_card( std::string& fc_content, std::multimap<std::string, std::string>& fmesh_params, int fcid )
{
  // convert '=' chars to spaces 
  size_t found;
   
  found = fc_content.find_first_of('=');
  while (found!= fc_content.npos )
  {
    fc_content[found] = ' ';
    found = fc_content.find_first_of('=',found+1);
  }

  std::stringstream tokenizer(fc_content);

  // skip tokens until 'dagmc' found
  bool found_dagmc = false;
  while( tokenizer ){
    std::string dagmc; 
    tokenizer >> dagmc;
    if( dagmc == "dagmc" ){
      found_dagmc = true;
      break;
    }
  }

  if( !found_dagmc )
  {
    std::cerr << "Error: FC" << fcid << " card is incorrectly formatted" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string last_key;
  while(tokenizer){
    std::string token;
    tokenizer >> token;

    if( token == "" ) continue;
    if( token == "-dagmc") break; // stop parsing if -dagmc encountered

    if( last_key == "" ){ last_key = token; }
    else{ 
      fmesh_params.insert(std::make_pair(last_key,token));
      last_key = "";
    }
  }

  if( last_key != "" ){
    std::cerr << "Warning: FC" << fcid << " card has unused key '" << last_key << "'" << std::endl;
  }
}

// Convenience methods
// Copy comment string out of fortran's data structure, and get it into comment_str
std::string copyComments(char* fort_comment, int* n_comment_lines)
{
    std::string comment_str; 

    const unsigned int fort_line_len = 75;
    unsigned int comment_len = fort_line_len * *n_comment_lines;

    // Need to turn it into a c-style string first
    char* c_comment = new char[(comment_len+1)];
    
    memcpy(c_comment,fort_comment,comment_len);
    c_comment[comment_len]='\0';
    
    comment_str = c_comment;
    delete[] c_comment;
    return comment_str;
 }


void dagmc_fmesh_setup_mesh_( int* /*ipt*/, int* id, 
                              double* energy_mesh, int* n_energy_mesh, int* tot_energy_bin, 
                              char* fort_comment, int* n_comment_lines, int* is_collision_tally )
{

  std::cout << "Mesh tally " << *id << " has these " << *n_energy_mesh << " energy bins: " << std::endl;
  for( int i = 0; i < *n_energy_mesh; ++i )
  {
    std::cout << "     " << energy_mesh[i] << std::endl;
  }
  // ToDo:  We may want to do an action if there is not a total energy bin
  std::cout << "tot bin: " << (*tot_energy_bin ? "yes" : "no") << std::endl;
  
  if( *n_comment_lines <= 0 )
  {
    std::cerr << "FMESH" << *id << " has geom=dag without matching FC card" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string comment_str = copyComments(fort_comment, n_comment_lines);

  // Copy emesh bin boundaries from MCNP (includes 0 MeV)
  std::vector<double> emesh_boundaries;
  for( int i = 0; i < *n_energy_mesh; ++i )
  {
    emesh_boundaries.push_back(energy_mesh[i]);
  }

  // Parse FC card and create input data for MeshTally
  std::multimap<std::string, std::string> fc_settings;
  parse_fc_card( comment_str, fc_settings, *id );

  // determine the user-specified tally type
  std::string type = "unstr_track";
  if( fc_settings.find("type") != fc_settings.end() )
  {
    type = (*fc_settings.find("type")).second;
    if( fc_settings.count("type") > 1 )
    {
      std::cerr << "Warning: FC" << *id << " has multiple 'type' keywords, using " << type << std::endl;
    }
    // remove the type keywords
    fc_settings.erase("type"); 
  }
     
  // Set whether the tally type is a collision tally 
  if (type.find("coll") != std::string::npos)
  {
       *is_collision_tally = true;
  }
  else
  {
       *is_collision_tally = false;
  } 

  tallyManager.addNewTally(*id, type, fc_settings, emesh_boundaries);
}


/********************************************************************
 * Runtape and MPI calls
 ********************************************************************/ 

/**
 * Get a fortran pointer to the tally array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */
/*
void dagmc_fmesh_get_tally_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
 
  data = all_tallies[*fmesh_index]->get_tally_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}
*/

/**
 * Get a fortran pointer to the error array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */

/*
void dagmc_fmesh_get_error_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
 
  data = all_tallies[*fmesh_index]->get_error_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}
*/
/**
 * Get a fortran pointer to the scratch array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */
/*
void dagmc_fmesh_get_scratch_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
  
  data = all_tallies[*fmesh_index]->get_scratch_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}
*/
/**
 * Set the tally and error arrays of the specified mesh tally to all zeros.
 * Called when an MPI subtask has just sent all its tally and error values
 * back to the master task.
 */
/*
void dagmc_fmesh_clear_data_( int* fmesh_index ){

  all_tallies[*fmesh_index]->zero_tally_data( );
}
*/

/**
 * Add the values in this mesh's scratch array to its tally array.
 * Called when merging together values from MPI subtasks at the master task.
 */
/*
void dagmc_fmesh_add_scratch_to_tally_( int* fmesh_index ){
  double* data, *scratch;
  int length, scratchlength;

  data = all_tallies[*fmesh_index]->get_tally_data( length );
  scratch = all_tallies[*fmesh_index]->get_scratch_data( scratchlength );
  
  assert( scratchlength >= length );

  for( int i = 0; i < length; ++i ){
    data[i] += scratch[i];
  }
}
*/
/**
 * Add the values in this mesh's scratch array to its error array.
 * Called when merging together values from MPI subtasks at the master task.
 */
/*
void dagmc_fmesh_add_scratch_to_error_( int* fmesh_index ){
  double* data, *scratch;
  int length, scratchlength;

  data = all_tallies[*fmesh_index]->get_error_data( length );
  scratch = all_tallies[*fmesh_index]->get_scratch_data( scratchlength );
  
  assert( scratchlength >= length );

  for( int i = 0; i < length; ++i ){
    data[i] += scratch[i];
  }
}
*/
/********************************************************************
 * Routine calls from fmesh_mod: track length reports and print commands
 ********************************************************************/ 

/**
 * Called from fortran when a source particle history ends
 */
void dagmc_fmesh_end_history_()
{
  tallyManager.end_history();

#ifdef MESHTAL_DEBUG
  std::cout << "* History ends *" << std::endl;
#endif
}

/**
 * Called from fortran to score a particular track length on fmesh with fmesh_index 
 * This is a tracklength tally
 * @param fmesh_index index of mesh tally
 * @param x, y, z - particle location
 * @param u, v, w - particle direction
 * @param erg particle energy
 * @param wgt particle weight
 * @param d   track length
 * @param icl MCNP global variable, current cell id 
 */
void dagmc_fmesh_score_( double *x, double *y, double *z,
                         double *u, double *v, double *w, 
                         double *erg,double *wgt, 
                         double *d, int *icl )
{
  tallyManager.set_track_event(*x, *y, *z, *u, *v, *w, *erg, *wgt, *d, *icl);
  tallyManager.update_tallies();

#ifdef MESHTAL_DEBUG
    std::cout << "meshtal particle loc: " << *x << ", " << *y << ", " << *z << std::endl;
    std::cout << "meshtal particle dir: " << *u << ", " << *v << ", " << *w << std::endl;
    std::cout << "meshtal track length: " << *d << std::endl;
#endif
}

/**
 * Called from fortan to instruct tallies to print data to the appropriate file
 * @param sp_norm "Source Particle Normalization" - the number of source particles so far
 */
void dagmc_fmesh_print_( double* sp_norm)
{
  tallyManager.write_data(*sp_norm);
}

/**  
 *   Obtains the collision position (x,y,z), the particle weight (wgt), the
 *   total macroscopic cross section of the current cell (ple), and the
 *   particle energy (erg) from MCNP for use in the KDE collision tally.
 *
 *   called from hstory.F90
 */
void dagmc_collision_score_( double* x,   double* y, double* z, 
                             double* erg, double* wgt,
                             double* ple, int* icl )
{
  tallyManager.set_collision_event(*x, *y, *z, *erg, *wgt, *ple, *icl);
  tallyManager.update_tallies();
}
