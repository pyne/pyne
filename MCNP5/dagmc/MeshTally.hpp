// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESH_TALLY_HPP
#define DAGMC_MESH_TALLY_HPP

#include <map>
#include <set>
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
 * \brief Defines an abstract MeshTally interface
 *
 * MeshTally is an abstract class derived from Tally that defines all the
 * variables and methods needed to implement a generic MeshTally for use in
 * Monte Carlo particle transport codes.  All Derived classes must implement
 * the following methods from Tally
 *
 *     1) compute_score(const TallyEvent& event)
 *     2) write_data(double num_histories)
 *
 * A simple version of the end_history() method is already included, so this
 * does not need to be implemented unless additional or alternative actions
 * are required by a Derived class.
 *
 * Note that three arrays are available for storing the mesh tally data
 *
 *    1) tally_data: stores sum of scores for all particle histories
 *    2) error_data: stores data needed to determine error in tally results
 *    3) temp_tally_data: stores sum of scores for a single history
 *
 * Each element in these three data arrays represents one tally point and
 * one energy bin.  They are ordered first by tally point, and then by energy
 * bin.  Derived classes can easily access/modify individual elements using
 * the protected get_data() method.
 *
 * ==================
 * Input/Output Files
 * ==================
 *
 * All MeshTally objects are REQUIRED to include "inp"="input_filename" as
 * a TallyOption in the TallyInput struct defined in Tally.hpp and set through
 * the TallyManager.  This input file contains all of the mesh data that is
 * needed to compute the mesh tally scores.  It must be created in a file format
 * that is supported by the Mesh-Oriented Database (MOAB), which includes both
 * H5M and VTK options.  Source code and more information on MOAB can be found
 * at https://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB
 *
 * In addition to the "inp" key, all MeshTally objects can also include an
 * optional "out"="output_filename" key-value pair.  If the "out" key is not
 * included as a TallyOption, then the default is a H5M file format named
 * meshtal<tally_id>.h5m.  To write to a different file format that is
 * supported by MOAB, simply add the desired extension to the output filename
 * (i.e. "out"="filename.vtk" will write results to the VTK format).
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

    /**
     * \brief Updates MeshTally when a particle history ends
     */
    virtual void end_history();

    /**
     * \brief Write tally and error results to the mesh tally's output file
     * \param num_particles the number of source particles tracked
     * \param multiplier an optional constant multiplication factor
     */
    virtual void print(double num_particles, double multiplier = 1.0) = 0;

    // >>> TALLY DATA ACCESS METHODS

    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param length output parameter containing size of data array
     * \return pointer to the data array
     *
     * These three methods can be useful for implementing MeshTally
     * functionality in parallel.
     */
    // virtual double* get_tally_data(int& length);
    // virtual double* get_error_data(int& length);
    // virtual double* get_scratch_data(int& length);

    /**
     * \brief Resets all of the mesh tally data arrays to zero
     */
    // virtual void zero_tally_data();

  protected:
    /// Name of file to which the final tally results will be written
    std::string output_filename;

    /// Name of file that contains the mesh description (separate from geometry)
    std::string input_filename;

    /// Entity handle for the MOAB mesh data used for this mesh tally
    moab::EntityHandle tally_mesh_set;

    /// Set of tally points (cells, nodes, etc) for this mesh tally
    moab::Range tally_points;

    /// Tag arrays for storing energy bin labels
    std::vector<moab::Tag> tally_tags, error_tags;

    // >>> TALLY DATA ARRAYS

    /// Data array for storing sum of scores for all particle histories
    // std::vector<double> tally_data;

    /// Data array for determining error in tally results
    // std::vector<double> error_data;

    /// Data array for storing sum of scores for a single history
    // std::vector<double> temp_tally_data;

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
 //   void resize_data_arrays(unsigned int num_tally_points);

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
  //  double& get_data(std::vector<double>& data,
   //                  moab::EntityHandle tally_point,
    //                 unsigned energy_bin = 0);

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
