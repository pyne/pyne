// MCNP5/dagmc/meshtal_funcs.cpp

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>

#include "meshtal_funcs.h"
#include "TallyManager.hpp"

// create a tally manager to handle all DAGMC tally actions
TallyManager tallyManager;

//---------------------------------------------------------------------------//
// INITIALIZATION AND SETUP METHODS
//---------------------------------------------------------------------------//
/**
 * \brief Convert contents of FC card to multimap<string, string>
 * \param[in] fc_content the FC card's comment content as a string
 * \param[out] fmesh_params the output data and comments, as a multimap
 * \param[in] fcid the tally ID of the FC card
 * \return true on success; false if input has serious formatting problems
 *         to make parameter parsing impossible
 */
static void parse_fc_card(std::string& fc_content,
                          std::multimap<std::string, std::string>& fmesh_params,
                          int fcid)
{
    // convert '=' chars to spaces 
    size_t found;
    found = fc_content.find_first_of('=');

    while (found!= fc_content.npos)
    {
        fc_content[found] = ' ';
        found = fc_content.find_first_of('=',found+1);
    }

    std::stringstream tokenizer(fc_content);

    // skip tokens until 'dagmc' found
    bool found_dagmc = false;

    while(tokenizer)
    {
        std::string dagmc; 
        tokenizer >> dagmc;
        if(dagmc == "dagmc")
        {
            found_dagmc = true;
            break;
        }
    }

    if(!found_dagmc)
    {
        std::cerr << "Error: FC" << fcid << " card is incorrectly formatted" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string last_key;

    while(tokenizer)
    {
        std::string token;
        tokenizer >> token;

        if( token == "" ) continue;
        if( token == "-dagmc") break; // stop parsing if -dagmc encountered

        if( last_key == "" ){ last_key = token; }
        else
        { 
            fmesh_params.insert(std::make_pair(last_key,token));
            last_key = "";
        }
    }

    if( last_key != "" )
    {
        std::cerr << "Warning: FC" << fcid << " card has unused key '" 
                  << last_key << "'" << std::endl;
    }
}
//---------------------------------------------------------------------------//
/**
 * \brief Copy comment string from fortran's data structure
 * \param[in] fort_comment the comment string to copy
 * \param[in] n_comment_lines the number of comment lines
 * \return a std::string version of the original comment string
 */
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
//---------------------------------------------------------------------------//
/**
 * \brief Sets up a DAGMC mesh tally in Fortran
 * \param[in] ipt the type of particle; currently unused
 * \param[in] id the unique ID for the tally defined by FMESH
 * \param[in] energy_mesh the energy bin boundaries
 * \param[in] n_energy_mesh the number of energy bin boundaries
 * \param[in] tot_energy_bin determines if total energy bin was requested
 * \param[in] fort_comment the FC card comment
 * \param[in] n_comment_lines the number of comment lines
 * \param[out] is_collision_tally indicates that tally uses collision estimator
 */
void dagmc_fmesh_setup_mesh_(int* /*ipt*/, int* id, int* fmesh_idx,
                             double* energy_mesh, int* n_energy_mesh,
                             int* tot_energy_bin, 
                             char* fort_comment, int* n_comment_lines,
                             int* is_collision_tally)
{
    std::cout << "Mesh tally " << *id << " has these " << *n_energy_mesh
              << " energy bins: " << std::endl;

    for(int i = 0; i < *n_energy_mesh; ++i)
    {
        std::cout << "     " << energy_mesh[i] << std::endl;
    }

    // TODO: Total energy bin is currently always on unless one bin is used
    std::cout << "tot bin: " << (*tot_energy_bin ? "yes" : "no") << std::endl;

    if(*n_comment_lines <= 0)
    {
        std::cerr << "FMESH" << *id
                  << " has geom=dag without matching FC card" << std::endl;
        exit(EXIT_FAILURE);
    }
  
    // Copy emesh bin boundaries from MCNP (includes 0.0 MeV)
    std::vector<double> emesh_boundaries;
    for(int i = 0; i < *n_energy_mesh; ++i)
    {
        emesh_boundaries.push_back(energy_mesh[i]);
    }

    // Parse FC card and create input data for MeshTally
    std::multimap<std::string, std::string> fc_settings;
    std::string comment_str = copyComments(fort_comment, n_comment_lines);
    parse_fc_card(comment_str, fc_settings, *id);

    // determine the user-specified tally type
    std::string type = "unstr_track";

    if(fc_settings.find("type") != fc_settings.end())
    {
        type = (*fc_settings.find("type")).second;

        if(fc_settings.count("type") > 1)
        {
            std::cerr << "Warning: FC" << *id
                      << " has multiple 'type' keywords, using " << type << std::endl;
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

    tallyManager.addNewTally(*id, type, emesh_boundaries, fc_settings);

    // Add tally multiplier, if it exists  
    if (*fmesh_idx != -1)
    {
       // Create a zero-based version of fmesh_idx to use as multiplier_id
       int multiplier_id = *fmesh_idx - 1;

       tallyManager.addNewMultiplier(multiplier_id);
       tallyManager.addMultiplierToTally(multiplier_id, *id);
    }
}
//---------------------------------------------------------------------------//
// RUNTAPE AND MPI METHODS
//---------------------------------------------------------------------------//
/**
 * \brief Get fortran pointer to the tally data for the given tally
 * \param[in] tally_id the unique ID of the tally
 * \param[out] fortran_data_pointer pointer to the tally data
 *
 * Called when the tally data needs to be written or read from a runtpe file
 * or an MPI stream.
 */
void dagmc_fmesh_get_tally_data_(int* tally_id, void* fortran_data_pointer)
{
    double* data;
    int length;

    data = tallyManager.getTallyData(*tally_id, length);
    FMESH_FUNC(dagmc_make_fortran_pointer)(fortran_data_pointer, data, &length);
}
//---------------------------------------------------------------------------//
/**
 * \brief Get fortran pointer to the error data for the given tally
 * \param[in] tally_id the unique ID of the tally
 * \param[out] fortran_data_pointer the pointer to the error data
 *
 * Called when the error data needs to be written or read from a runtpe file
 * or an MPI stream.
 */
void dagmc_fmesh_get_error_data_(int* tally_id, void* fortran_data_pointer)
{
    double* data;
    int length;

    data = tallyManager.getErrorData(*tally_id, length);
    FMESH_FUNC(dagmc_make_fortran_pointer)(fortran_data_pointer, data, &length);
}
//---------------------------------------------------------------------------//
/**
 * \brief Get fortran pointer to the scratch data for the given tally
 * \param[in] tally_id the unique ID of the tally
 * \param[out] fortran_data_pointer the pointer to the scratch data
 *
 * Called when the scratch data needs to be written or read from a runtpe file
 * or an MPI stream.
 */
void dagmc_fmesh_get_scratch_data_(int* tally_id, void* fortran_data_pointer)
{
    double* data;
    int length;

    data = tallyManager.getScratchData(*tally_id, length);
    FMESH_FUNC(dagmc_make_fortran_pointer)(fortran_data_pointer, data, &length);
}
//---------------------------------------------------------------------------//
/**
 * \brief Reset all data for all tallies to zeroes
 *
 * Called when an MPI subtask has just sent all its tally and error values
 * back to the master task.
 */
void dagmc_fmesh_clear_data_()
{
    tallyManager.zeroAllTallyData();
}
//---------------------------------------------------------------------------//
/**
 * \brief Add values in the given tally's scratch data to its tally data
 * \param[in] tally_id the unique ID of the tally
 *
 * Called when merging together values from MPI subtasks at the master task.
 */
void dagmc_fmesh_add_scratch_to_tally_(int* tally_id)
{
    double* data;
    double* scratch;
    int length, scratchlength;

    data = tallyManager.getTallyData(*tally_id, length);
    scratch = tallyManager.getScratchData(*tally_id, scratchlength);

    assert(scratchlength >= length);

    for(int i = 0; i < length; ++i)
    {
        data[i] += scratch[i];
    }
}
//---------------------------------------------------------------------------//
/**
 * \brief Add values in the given tally's scratch data to its error data
 * \param[in] tally_id the unique ID of the tally
 *
 * Called when merging together values from MPI subtasks at the master task.
 */
void dagmc_fmesh_add_scratch_to_error_(int* tally_id)
{
    double* error_data;
    double* scratch;
    int length, scratchlength;

    error_data = tallyManager.getErrorData(*tally_id, length);
    scratch = tallyManager.getScratchData(*tally_id, scratchlength);

    assert(scratchlength >= length);

    for(int i = 0; i < length; ++i)
    {
        error_data[i] += scratch[i];
    }
}
//---------------------------------------------------------------------------//
// ROUTINE FMESH METHODS
//---------------------------------------------------------------------------//
/**
 * \brief Called from fortran when a particle history ends
 */
void dagmc_fmesh_end_history_()
{
    tallyManager.endHistory();

#ifdef MESHTAL_DEBUG
    std::cout << "* History ends *" << std::endl;
#endif
}
//---------------------------------------------------------------------------//
/**
 * \brief Called from fortran to score a track event
 * \param[in] x, y, z the position of the particle
 * \param[in] u, v, w the direction of the particle
 * \param[in] erg the energy of the particle
 * \param[in] wgt the weight of the particle
 * \param[in] d the track length
 * \param[in] icl the current cell ID (MCNP global variable)
 * 
 * This function is called once per track event.
 */
void dagmc_fmesh_score_(double *x, double *y, double *z,
                        double *u, double *v, double *w, 
                        double *erg,double *wgt, 
                        double *d, int *icl)
{
#ifdef MESHTAL_DEBUG
    std::cout << "particle loc: " << *x << ", " << *y << ", " << *z << std::endl;
    std::cout << "particle dir: " << *u << ", " << *v << ", " << *w << std::endl;
    std::cout << "track length: " << *d << std::endl;
#endif

    tallyManager.setTrackEvent(*x, *y, *z, *u, *v, *w, *erg, *wgt, *d, *icl);
    tallyManager.updateTallies();
}
//---------------------------------------------------------------------------//
/**
 * \brief Called from fortran to instruct tallies to print data to file
 * \param[in] sp_norm "Source Particle Normalization", number of source particles
 */
void dagmc_fmesh_print_(double* sp_norm)
{
    tallyManager.writeData(*sp_norm);
}
//---------------------------------------------------------------------------//
/**
 * \brief Called from hstory.F90 to score a collision event
 * \param[in] x, y, z the position of the particle
 * \param[in] erg the energy of the particle
 * \param[in] wgt the weight of the particle
 * \param[in] ple the total macroscopic cross section of the current cell
 * \param[in] icl the current cell ID (MCNP global variable)
 * 
 * This function is called once per collision event.
 */
void dagmc_collision_score_(double* x,   double* y, double* z, 
                            double* erg, double* wgt,
                            double* ple, int* icl)
{
    tallyManager.setCollisionEvent(*x, *y, *z, *erg, *wgt, *ple, *icl);
    tallyManager.updateTallies();
}
//---------------------------------------------------------------------------//
/**
 * \brief Called from fmesh_mod.F90 to update tally multipliers
 * \param[in] fmesh_idx the fmesh index for multiplier to be updated
 * \param[in] value the value of the multiplier
 */
void dagmc_update_multiplier_(int* fmesh_idx, double* value)
{
   tallyManager.updateMultiplier(*fmesh_idx-1, *value);
}

//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/meshtal_funcs.cpp
