// MCNP5/dagmc/MeshTally.hpp

#ifndef DAGMC_MESH_TALLY_HPP
#define DAGMC_MESH_TALLY_HPP

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
 * In general, the default end_history() method that is implemented in Tally
 * will also be sufficient for most MeshTally objects that use the TallyData
 * structure for storing their data.  However, it can be extended or overridden
 * as needed by Derived classes.
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
     * \param[in] input user-defined input parameters for this mesh tally
     */
    explicit MeshTally(const TallyInput& input);

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~MeshTally(){}

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

    // >>> PROTECTED METHODS

    /**
     * \brief Determines entity index corresponding to tally point
     * \param[in] tally_point entity handle representing tally point
     * \return entity index for given tally point
     */
    unsigned int get_entity_index(moab::EntityHandle tally_point);

    /**
     * \brief Loads the MOAB mesh data from the input file for this mesh tally
     * \param[in] mbi the MOAB interface for this mesh tally
     * \param[out] mesh_set entity handle for the mesh set that will be created
     * \return the MOAB ErrorCode value
     */
    moab::ErrorCode load_moab_mesh(moab::Interface* mbi,
                                   moab::EntityHandle& mesh_set);

    /**
     * \brief Defines the set of tally points to use for this mesh tally
     * \param[in] mesh_elements the set of mesh elements to use as tally points
     *
     * Note that this method calls resize_data_arrays() to set the tally
     * data arrays for the given number of tally points.
     */
    void set_tally_points(const moab::Range& mesh_elements);

    /**
     * \brief Reduces a MOAB mesh set to include only its 3D elements
     * \param[in] mbi the MOAB interface for this mesh tally
     * \param[in, out] mesh_set entity handle for the mesh set that will be reduced
     * \param[out] mesh_elements stores 3D elements that were added to the mesh set
     * \return the MOAB ErrorCode value
     *
     * NOTE: this method will overwrite the mesh set
     */
    moab::ErrorCode reduce_meshset_to_3D(moab::Interface* mbi,
                                         moab::EntityHandle& mesh_set,
                                         moab::Range& mesh_elements);

    /**
     * \brief Sets up tally value and error labels for all energy bins
     * \param[in] mbi the MOAB interface for this mesh tally
     * \param[in] prefix additional string to be added before each label
     * \return the MOAB ErrorCode value
     *
     * Note that labels are stored as MOAB tag handles in the tally_tags
     * and error_tags arrays.
     */
    moab::ErrorCode setup_tags(moab::Interface* mbi, const char* prefix="");

    /**
     * \brief Adds weight * score to the mesh tally for the tally point
     * \param[in] tally_point entity handle representing tally point
     * \param[in] weight the multiplier value for the score to be tallied
     * \param[in] score the score that is to be tallied
     * \param[in] ebin the energy bin index corresponding to the energy
     *
     * The weight and ebin can be obtained using Tally::get_score_multiplier
     * and Tally::get_energy_bin respectively.
     */
    void add_score_to_mesh_tally(const moab::EntityHandle& tally_point, 
                                 double weight, double score, unsigned int ebin);
};

#endif // DAGMC_MESHTALLY_HPP

// end of MCNP5/dagmc/MeshTally.hpp
