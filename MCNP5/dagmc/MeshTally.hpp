// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESHTALLY_H
#define DAGMC_MESHTALLY_H

#include <cassert>
#include <map>
#include <string>
#include <vector>

#include "moab/Range.hpp"

//===========================================================================//
/**
 * \struct MeshTallyInput
 * \brief Input data needed to set up a mesh tally
 */
//===========================================================================//
struct MeshTallyInput
{
    /// Typedef for map that stores optional mesh tally input parameters
    typedef std::multimap<std::string, std::string> TallyOptions;

    /// User-specified ID for this mesh tally
    int tally_id;

    /// Energy bin boundaries defined for all mesh tally points
    std::vector<double> energy_bin_bounds;

    /// If true, add an extra energy bin to tally all energy levels
    bool total_energy_bin;

    /// Optional input parameters requested by user
    TallyOptions options;
};

// forward declaration
namespace moab{
  class Interface;
}

//===========================================================================//
/**
 * \class MeshTally
 * \brief Defines a basic mesh tally interface
 *
 * MeshTally is a Base class that defines the variables and methods that are
 * typically needed to implement a mesh tally for use in Monte Carlo particle
 * transport codes.  Some basic functionality is already included but the
 * following functions must be implemented in all Derived classes
 * 
 *     1) get_score_weight
 *     2) add_score_to_tally
 *     3) end_history
 *     4) write_results
 *
 * Note that three arrays are defined for storing mesh tally data
 *
 *    1) tally_data: stores sum of scores for all particle histories
 *    2) error_data: stores data needed to determine error in tally results
 *    3) temp_tally_data: stores sum of scores for a single history
 */
//===========================================================================//
class MeshTally {

  protected:
    /**
     * \brief Constructor
     * \param input user-defined input parameters for this mesh tally
     */
    explicit MeshTally(const MeshTallyInput& input);

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~MeshTally(){}

    // >>> PUBLIC INTERFACE

    /**
     * \brief Computes weighting factor for individual tally scores
     */
    //virtual void get_score_weight() = 0;

    /**
     * \brief Adds score to the mesh tally
     */ 
    //virtual void add_score_to_tally() = 0;

    /**
     * \brief Updates tally information when the history of a particle ends
     */
    virtual void end_history() = 0;

    /**
     * \brief Write tally and error results to the mesh tally's output file
     */ 
    //virtual void write_results() = 0;

    // TODO remove print function once we change the Derived class implementations
    /**
     * Print / write results to the AMeshTally's output file.
     * @param sp_norm The number of source particles, as reported from within mcnp's fortran code.
     * @param fmesh_fact Multiplication factor from fmesh card.  
     */
    virtual void print(double sp_norm, double fmesh_fact) = 0;

    // >>> TALLY DATA ACCESS FUNCTIONS

    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param length output parameter containing size of data array
     * \return pointer to the tally data array
     */
    virtual double* get_tally_data(int& length);
    virtual double* get_error_data(int& length);
    virtual double* get_scratch_data(int& length);

    /**
     * /brief Resets all of the tally data arrays to zero
     */
    virtual void zero_tally_data();

  protected:
    /// Input data defined by user for this mesh tally
    MeshTallyInput input_data;

    /// Data array for storing sum of scores for all particle histories
    std::vector<double> tally_data;

    /// Data array for determining error in tally results
    std::vector<double> error_data;

    /// Data array for storing sum of scores for a single history
    std::vector<double> temp_tally_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins; 

    // >>> PROTECTED FUNCTIONS

    /**
     * /brief Resize data arrays to hold a given number of points
     * /param size defines the new size of the data array
     *
     * Arrays will be resized to the given size * the number of energy bins
     */
    void resize_data_arrays(unsigned int size);

    // >>> MOAB-BASED DATA/FUNCTIONS TODO still need to remove this MOAB dependency

    unsigned int ent_idx(moab::EntityHandle eh)
    {
        unsigned int ret = tally_ents.index(eh);
        assert(ret < tally_ents.size());
        return ret;
    }

    double& data_ref(std::vector<double>& data,
                     moab::EntityHandle eh,
                     unsigned ebin = 0)
    {
        assert(ebin < num_energy_bins);
        return data[ ent_idx(eh) * num_energy_bins + ebin ];
    }

    moab::ErrorCode setup_tags(moab::Interface* mbi, const char* prefix="");

    /// entities on which to compute tally
    moab::Range tally_ents;

    /// Tag arrays
    std::vector<moab::Tag> tally_tags, error_tags; 
};

#endif // DAGMC_MESHTALLY_H

// end of MCNP5/dagmc/MeshTally.hpp
