// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESH_TALLY_HPP
#define DAGMC_MESH_TALLY_HPP

#include <map>
#include <string>
#include <vector>

#include "moab/Range.hpp"
#include "Tally.hpp"

// forward declaration
namespace moab {
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
 * following methods must be implemented in all Derived classes
 * 
 *     1) compute_score
 *     2) end_history
 *     3) print
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
 * the get_data() method.
 */
//===========================================================================//
class MeshTally : public Tally
{

  protected:
    /**
     * \brief Constructor
     * \param input user-defined input parameters for this mesh tally
     */
    MeshTally(const TallyInput& input);

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~MeshTally(){}

    // >>> PUBLIC INTERFACE

    // ToDo:  These comments, perhaps in modified form, will go with 
    //        the implementation of these methods in derived classes.
    /**
     * \brief Computes mesh tally scores for the given tally event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin index representing energy bin
     * TODO remove ebin as parameter since this can be computed from energy?
     */
    // jcz note: compare to update() in base class
    // virtual void compute_score(const TallyEvent& event, int ebin) = 0;

    /**
     * \brief Updates tally information when a particle history ends
     */
    virtual void end_history();

    // >>> TALLY DATA ACCESS METHODS

    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param length output parameter containing size of data array
     * \return pointer to the data array
     *
     * These three methods can be useful for implementing MeshTally
     * functionality in parallel.
     */
    virtual double* get_tally_data(int& length);
    virtual double* get_error_data(int& length);
    virtual double* get_scratch_data(int& length);

    /**
     * \brief Resets all of the mesh tally data arrays to zero
     */
    virtual void zero_tally_data();

  protected:
    /// Name of file to which the final tally results will be written
    std::string output_filename;
    /// Name of the file that contains the mesh description separate from the geometry
    std::string input_filename;

    /// Entity handle for the MOAB mesh data used for this mesh tally
    moab::EntityHandle tally_mesh_set;

    /// Set of tally points (cells, nodes, etc) for this mesh tally
    moab::Range tally_points;

    /// Tag arrays for storing energy bin labels
    std::vector<moab::Tag> tally_tags, error_tags;

    // >>> TALLY DATA ARRAYS

    /// Data array for storing sum of scores for all particle histories
    std::vector<double> tally_data;

    /// Data array for determining error in tally results
    std::vector<double> error_data;

    /// Data array for storing sum of scores for a single history
    std::vector<double> temp_tally_data;

    /// Entity handles updated in current history; cleared by end_history()
    std::set<moab::EntityHandle> visited_this_history;

    // >>> PROTECTED METHODS

    /** 
     * \brief Add score to the tally for the given tally point
     * \param tally_point entity handle representing tally point
     * \param score the contribution to add to the tally
     * \param ebin the energy bin to which the score will be added
     */
    virtual void add_score_to_tally(moab::EntityHandle tally_point,
                                    double score,
                                    int ebin);

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
     * \brief Loads the MOAB mesh data from the input file for this mesh tally
     * \param mbi the MOAB interface for this mesh tally
     * \param mesh_set entity handle for the mesh set that will be created
     * \return the MOAB ErrorCode value
     */
    moab::ErrorCode load_moab_mesh(moab::Interface* mbi,
                                   moab::EntityHandle& mesh_set);

    /**
     * \brief Defines the set of tally points to use for this mesh tally
     * \param mesh_elements the set of mesh elements to use as tally points
     *
     * Note that this method calls resize_data_arrays() to set the tally
     * data arrays for the given number of tally points.
     */
    void set_tally_points(const moab::Range& mesh_elements);

    /**
     * \brief Reduces a MOAB mesh set to include only its 3D elements
     * \param mbi the MOAB interface for this mesh tally
     * \param mesh_set entity handle for the mesh set that will be reduced
     * \param mesh_elements stores 3D elements that were added to the mesh set
     * \return the MOAB ErrorCode value
     *
     * NOTE: this method will overwrite the mesh set
     */
    moab::ErrorCode reduce_meshset_to_3D(moab::Interface* mbi,
                                         moab::EntityHandle& mesh_set,
                                         moab::Range& mesh_elements);

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

#endif // DAGMC_MESHTALLY_HPP

// end of MCNP5/dagmc/MeshTally.hpp
