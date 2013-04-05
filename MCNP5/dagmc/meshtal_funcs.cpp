#include "meshtal_funcs.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>

#include "moab/Core.hpp"
#include "MeshTally.hpp"
#include "TallyEvent.hpp"
#include "TrackLengthMeshTally.hpp"

using moab::TrackLengthMeshTally;

#include "KDEMeshTally.hpp"

/********************************************************************
 * File statics defining global meshtal data
 ********************************************************************/ 


/// number of source particle histories - used for debugging only
static int history_count = 0; 

/// number of calls to dagmc_fmesh_score - used for debugging only
static int score_count = 0;

/// The following four lists are indexed by the fmesh_index values
/// used to index fmesh tallies in the fortran code.  Entries in 
/// these lists may be NULL.

/// List of all track length tallies that have been created.
/// These tallies respond to fmesh mesh_score calls
static std::vector< TrackLengthMeshTally* > tracklen_tallies;

/// List of all KDE collision tallies that have been created.
/// These tallies respond to dagmc_kde_tally calls
static std::vector< KDEMeshTally* > kde_coll_tallies;

/// List of all KDE track/subtrack tallies that have been created.
/// These tallies respond to fmesh mesh_score calls
static std::vector< KDEMeshTally* > kde_track_tallies;

/// List of all dagmc-typed fmesh tallies: the entries are duplicates
/// of the above lists with the same indexing.
/// All our tallies are stored here, and respond to fmesh_print commands.
static std::vector< MeshTally* > all_tallies;

/// pointer to MCNP5's current cell ID variable (icl) 
static const int* current_mcnp_cell;

void mcnp_weight_calculation( int* index, double* erg, double* wgt, 
                              double* dist, double* score_result )
{
    FMESH_FUNC(dagmc_mesh_score)( index, erg, wgt, dist, score_result );
}



/********************************************************************
 * Initialization and setup functions
 ********************************************************************/ 


static bool initialized = false;

/** 
 * Called at least once from fmesh_mod on program initialization;
 * in runtpe or MPI modes may be called multiple times.  
 */
void dagmc_fmesh_initialize_( const int* mcnp_icl ){

  if( initialized ) return;

  //std::cerr << "Executed DAGMC fmesh initialize." << std::endl;

  current_mcnp_cell = mcnp_icl;

  initialized = true;
}

/**
 * Convert the contents of an FC card to an fmesh_params_t (i.e. a multimap<string,string>)
 * @param fc_content The FC card's comment content as a string
 * @param results The output data
 * @param fcid The tally ID of the FC card
 * @return true on success, or false if the input has serious enough formatting problems
 *         to make parameter parsing impossible.
 */
static bool parse_fc_card( std::string& fc_content, MeshTallyInput::TallyOptions& results, int fcid ){

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

  if( !found_dagmc ){
    std::cerr << "Error: FC" << fcid << " card is incorrectly formatted" << std::endl;
    return false;
  }

  std::string last_key;
  while(tokenizer){
    std::string token;
    tokenizer >> token;

    if( token == "" ) continue;
    if( token == "-dagmc") break; // stop parsing if -dagmc encountered

    if( last_key == "" ){ last_key = token; }
    else{ 
      results.insert(std::make_pair(last_key,token));
      last_key = "";
    }

  }

  if( last_key != "" ){
    std::cerr << "Warning: FC" << fcid << " card has unused key '" << last_key << "'" << std::endl;
  }

  return true;

}



void dagmc_fmesh_setup_mesh_( int* /*ipt*/, int* id, int* fmesh_index, 
                              double* energy_mesh, int* n_energy_mesh, int* tot_energy_bin, 
                              char* fort_comment, int* n_comment_lines, int* is_collision_tally  )
{

  std::cerr << "Mesh tally " << *id << " has these " << *n_energy_mesh << " energy bins: " << std::endl;
  for( int i = 0; i < *n_energy_mesh; ++i ){
    std::cerr << "     " << energy_mesh[i] << std::endl;
  }
  std::cerr << "tot bin: " << (*tot_energy_bin ? "yes" : "no") << std::endl;
  
  if( *n_comment_lines <= 0 ){
    std::cerr << "FMESH" << *id << " has geom=dag without matching FC card" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string comment_str; 

  
  {
    // Copy comment string out of fortran's data structure, and get it into comment_str
    // Need to turn it into a c-style string first

    const unsigned int fort_line_len = 75;
    unsigned int comment_len = fort_line_len * *n_comment_lines;
    char* c_comment = new char[(comment_len+1)];
    
    memcpy(c_comment,fort_comment,comment_len);
    c_comment[comment_len]='\0';
    
    comment_str = c_comment;
    delete[] c_comment;

  }

  // Copy emesh bin boundaries from MCNP (includes 0 MeV)
  std::vector<double> emesh_boundaries;

  for( int i = 0; i < *n_energy_mesh; ++i ){
    emesh_boundaries.push_back(energy_mesh[i]);
  }
  
  MeshTallyInput fmesh_settings;
  fmesh_settings.tally_id = *id;
  fmesh_settings.energy_bin_bounds = emesh_boundaries;
  fmesh_settings.total_energy_bin = (*tot_energy_bin == 1);

  MeshTallyInput::TallyOptions& fc_settings = fmesh_settings.options;

  bool success = parse_fc_card( comment_str, fc_settings, *id );
  if( !success ){
    exit(EXIT_FAILURE);
  }

  // pad all tally lists with nulls up to a max of (*fmesh_index)
  while( all_tallies.size() <= (unsigned)(*fmesh_index) ){
    tracklen_tallies.push_back(NULL);
    kde_coll_tallies.push_back(NULL);
    kde_track_tallies.push_back(NULL);
    all_tallies.push_back(NULL);
    
  }

  // determine the user-specified tally type
  std::string type = "tracklen";
  *is_collision_tally = 0;

  if( fc_settings.find("type") != fc_settings.end() ){

    type = (*fc_settings.find("type")).second;
    if( fc_settings.count("type") > 1 ){
      std::cerr << "Warning: FC" << *id << " has multiple 'type' keywords, using " << type << std::endl;
    }
    
    // remove the type keywords
    fc_settings.erase("type"); 
  }
  
  moab::Interface* mbi = new moab::Core();

  MeshTally *new_tally;

  if( type == "tracklen" ){

    TrackLengthMeshTally* t = TrackLengthMeshTally::setup( fmesh_settings, mbi, current_mcnp_cell );
    new_tally = tracklen_tallies[*fmesh_index] = t;

  }
  else if( type == "kde_coll" ){
  
    *is_collision_tally = 1; 

    KDEMeshTally* kde = KDEMeshTally::setup( fmesh_settings, mbi );
    new_tally = kde_coll_tallies[*fmesh_index] = kde;

  }
  else if( type == "kde_track" || type == "kde_subtrack" ){

    KDEMeshTally::TallyType kde_type = KDEMeshTally::TRACKLENGTH;

    if ( type == "kde_subtrack" )
      kde_type = KDEMeshTally::SUBTRACK;

    KDEMeshTally* kde = KDEMeshTally::setup( fmesh_settings, mbi, kde_type );
    new_tally = kde_track_tallies[*fmesh_index] = kde;

  }
  else{
    std::cerr << "FC" << *id << " error: cannot make mesh tally of type " << type << std::endl;
    exit( EXIT_FAILURE );
  }

  all_tallies[*fmesh_index] = new_tally;
  
}

/********************************************************************
 * Runtape and MPI calls
 ********************************************************************/ 

/**
 * Get a fortran pointer to the tally array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */
void dagmc_fmesh_get_tally_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
 
  data = all_tallies[*fmesh_index]->get_tally_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}

/**
 * Get a fortran pointer to the error array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */

void dagmc_fmesh_get_error_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
 
  data = all_tallies[*fmesh_index]->get_error_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}

/**
 * Get a fortran pointer to the scratch array for the specified mesh tally.
 * Called when this data needs to be written or read from a runtpe file or 
 * an MPI stream.
 */
void dagmc_fmesh_get_scratch_data_( int* fmesh_index, void* fortran_data_pointer ){
  double* data; 
  int length;
  
  data = all_tallies[*fmesh_index]->get_scratch_data( length );
  FMESH_FUNC( dagmc_make_fortran_pointer )( fortran_data_pointer, data, &length );
}

/**
 * Set the tally and error arrays of the specified mesh tally to all zeros.
 * Called when an MPI subtask has just sent all its tally and error values
 * back to the master task.
 */
void dagmc_fmesh_clear_data_( int* fmesh_index ){

  all_tallies[*fmesh_index]->zero_tally_data( );

}

/**
 * Add the values in this mesh's scratch array to its tally array.
 * Called when merging together values from MPI subtasks at the master task.
 */
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

/**
 * Add the values in this mesh's scratch array to its error array.
 * Called when merging together values from MPI subtasks at the master task.
 */
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

/********************************************************************
 * Routine calls from fmesh_mod: track length reports and print commands
 ********************************************************************/ 

/**
 * Called from fortran when a source particle history ends
 */
void dagmc_fmesh_end_history_(){
  history_count += 1;

  for( std::vector< MeshTally* >::iterator i = all_tallies.begin();
       i != all_tallies.end(); ++i )
  {
    if( *i ){
      (*i)->end_history( );
    }
  }

#ifdef MESHTAL_DEBUG
  std::cout << "* History ends *" << std::endl;
#endif

}

/**
 * Called from fortran to score a particular track length on fmesh with fmesh_index 
 * @param fmesh_index index of mesh tally
 * @param x, y, z - particle location
 * @param u, v, w - particle direction
 * @param erg particle energy
 * @param wgt particle weight
 * @param d track length
 */
void dagmc_fmesh_score_( int *fmesh_index, double *x, double *y, double *z,
                         double *u, double *v, double *w, double *erg,double *wgt,double *d, int* ien )
{

  if( tracklen_tallies[*fmesh_index] ){
    score_count += 1;
    
    moab::CartVect ray = moab::CartVect(*x,*y,*z);
    moab::CartVect vec = moab::CartVect(*u,*v,*w);
    
    MCNPTrackParam param( fmesh_index, erg, wgt );
    
#ifdef MESHTAL_DEBUG
    std::cout << "meshtal particle: " << ray << " " << vec << " " << *d << std::endl;
#endif
    
    tracklen_tallies[*fmesh_index]->add_track_segment( ray, vec, *d, (*ien)-1, &param );

  }
  else if ( kde_track_tallies[*fmesh_index] ) {

    moab::CartVect loc( *x, *y, *z );
    moab::CartVect dir( *u, *v, *w );

    TallyEvent event(*erg, *wgt);
    event.set_track_event(*d, loc, dir);
    
    kde_track_tallies[*fmesh_index]->tally_track(event, (*ien)-1); 

  }

}

/**
 * Called from fortan to instruct a particular tally to print its data to the appropriate file
 * @param fmesh_index The mesh to be printed right now
 * @param sp_norm "Source Particle Normalization" - the number of source particles so far
 * @param fmesh_fact Multiplication factor for this tally, as recorded in fmesh_mod
 */
void dagmc_fmesh_print_( int* fmesh_index, double* sp_norm, double* fmesh_fact ){

  if( all_tallies[*fmesh_index] ){
    all_tallies[*fmesh_index]->print( *sp_norm, *fmesh_fact );
  }

}

/**  
 *   Obtains the collision position (x,y,z), the particle weight (wgt), the
 *   total macroscopic cross section of the current cell (ple), and the
 *   particle energy (erg) from MCNP for use in the KDE collision tally.
 *
 *   called from hstory.F90
 */
void dagmc_kde_tally_( double* x, double* y, double* z, double* wgt,
                       double* ple, double* erg )
{
  // counter for determining the fmesh_index of the KDE collision tallies
  int fmesh_index = 0;

  // Record collision on all valid KDE tallies
  for( std::vector<KDEMeshTally*>::iterator i = kde_coll_tallies.begin(); i!=kde_coll_tallies.end(); ++i ){
    if( *i ){

      int ien; // index of the energy bin for this collision 

      // ask Fortran to pick the energy bin for this collision
      FMESH_FUNC( dagmc_mesh_choose_ebin )( &fmesh_index, erg, &ien );

      if( ien == -1 ) continue; // erg falls outside of requested energy bins for this mesh

      ien -= 1; // convert fortran array index to C index

      moab::CartVect collision_loc( *x, *y, *z );

      TallyEvent event(*erg, *wgt);
      event.set_collision_event(*ple, collision_loc);

      (*i)->tally_collision(event, ien);
    }
    ++fmesh_index;
  }
}

