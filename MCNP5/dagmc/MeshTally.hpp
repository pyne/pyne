// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESHTALLY_H
#define DAGMC_MESHTALLY_H

#include <cassert>
#include <map>
#include <string>
#include <vector>

#include "moab/Range.hpp"

#include "TallyEvent.hpp"

// forward declaration
namespace moab {
  class Interface;
}

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
 *     2) compute_score
 *     3) add_score_to_tally
 *     4) end_history
 *     5) write_results
 *
 * Note that three arrays are available for storing mesh tally data
 *
 *    1) tally_data: stores sum of scores for all particle histories
 *    2) error_data: stores data needed to determine error in tally results
 *    3) temp_tally_data: stores sum of scores for a single history
 *
 * Each element in these three data arrays represents one tally point and
 * one energy bin.  They are ordered first by tally point, and then by energy
 * bin.  Derived classes can easily access/modify individual elements using
 * the get_data function.
 */
//===========================================================================//
class MeshTally
{

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
     * \brief Computes mesh tally scores for the given tally event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin index representing energy bin
     * TODO remove ebin as parameter since this can be computed from energy?
     */
    virtual void compute_score(const TallyEvent& event, int ebin) = 0;

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
     * \return pointer to the data array
     */
    virtual double* get_tally_data(int& length);
    virtual double* get_error_data(int& length);
    virtual double* get_scratch_data(int& length);

    /**
     * \brief Resets all of the mesh tally data arrays to zero
     */
    virtual void zero_tally_data();

  protected:
    /// Input data defined by user for this mesh tally
    MeshTallyInput input_data;

    /// Set of tally points (cells, nodes, etc) for this mesh tally
    moab::Range tally_points;

    /// Tag arrays for storing energy bin labels
    std::vector<moab::Tag> tally_tags, error_tags;

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
     * \brief Resize data arrays to hold all of the mesh tally data
     * \param num_tally_points number of tally points included in mesh tally
     *
     * Arrays will be resized to the given number of tally points multiplied
     * by the number of energy bins.
     */
    void resize_data_arrays(unsigned int num_tally_points);

    /**
     * \brief Determines entity index corresponding to tally point
     * \param tally_point entity handle representing tally point
     * \return entity index for given tally point 
     */
    unsigned int get_entity_index(moab::EntityHandle tally_point);

    /**
     * \brief Determines location of element in data array
     * \param data array containing element to be accessed
     * \param tally_point entity handle representing tally point
     * \param energy_bin index representing energy bin
     * \return reference to element in data array
     *
     * Enables direct access to the mesh tally data for the given tally point
     * and energy bin.
     */
    double& get_data(std::vector<double>& data,
                     moab::EntityHandle tally_point,
                     unsigned energy_bin = 0);

    /**
     * \brief Sets up tally value and error labels for all energy bins
     * \param mbi the MOAB interface for this mesh tally
     * \param prefix additional string to be added before each label
     * \return the MOAB ErrorCode value
     *
     * Note that labels are stored as MOAB tag handles in the tally_tags
     * and error_tags arrays.
     */
    moab::ErrorCode setup_tags(moab::Interface* mbi, const char* prefix="");
};

#endif // DAGMC_MESHTALLY_H

// end of MCNP5/dagmc/MeshTally.hpp
