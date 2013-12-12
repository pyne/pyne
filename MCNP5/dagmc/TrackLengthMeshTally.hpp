// MCNP5/dagmc/TrackLengthMeshTally.hpp

#ifndef DAGMC_TRACK_LENGTH_MESH_TALLY_HPP
#define DAGMC_TRACK_LENGTH_MESH_TALLY_HPP

#include <string>
#include <cassert>
#include <set>

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"

#include "Matrix3.hpp"
#include "MeshTally.hpp"
#include "TallyEvent.hpp"

namespace moab{

/* Forward Declarations */
class AdaptiveKDTree;
class OrientedBoxTreeTool;

//===========================================================================//
/**
 * \class TrackLengthMeshTally
 * \brief Represents an unstructured mesh tally based on particle tracks
 *
 * TrackLengthMeshTally is a concrete class derived from MeshTally that
 * implements the Tally interface to tally particle tracks on an unstructured
 * mesh as part of a Monte Carlo particle transport simulation.  If a
 * TrackLengthMeshTally object receives a TallyEvent type that is not
 * TallyEvent::TRACK, then no scores are computed.
 *
 * ==========
 * TallyInput
 * ==========
 *
 * The TallyInput struct needed to construct a TrackLengthMeshTally object is
 * defined in Tally.hpp and is set through the TallyManager when a Tally is
 * created.  Options that are currently available for TrackLengthMeshTally
 * objects include
 *
 * 1) "inp"="input_filename", "out"="output_filename"
 * --------------------------------------------------
 * These two options are processed through the MeshTally constructor.  The
 * "inp" key is REQUIRED for TrackLengthMeshTally objects, whereas the "out"
 * key is optional.  See MeshTally.hpp for more information.
 *
 * 2) "convex"="t/f", "convex"="true/fale"
 * ----------------------------------------------
 * If the convex option is set to true, the user is asserting that the input
 * mesh has convex geometry.  Single particle tracks that leave a convex mesh
 * are no longer tracked as they cannot re-enter the mesh.  The default value
 * for this option is false.
 *
 * 3) "conformal"="values"
 * -----------------------
 * If the input mesh is exactly conformal to the geometry, then conformal
 * logic should be used for scoring particle tracks across mesh cells or
 * errors may occur. This option allows the user to identify what geometric
 * cells are conformal to the input mesh.  These values can be defined as a
 * list of cells separated by commas (i.e. 1, 2, 5), and/or a range of cells
 * (i.e. 2-5).  Like the convex case, particle tracks that leave a conformal
 * mesh are no longer tracked.
 *
 * 4) "conf_surf_src"="t/f", "conf_surf_src"="true/false"
 * ------------------------------------------------------
 * If the conformal surface source option is set to true, then the user is
 * asserting that a surface source has been defined that is conformal to a
 * geometric surface.  This means that new particles that are born from
 * this source will use conformal logic to determine the first mesh cell
 * that is crossed.  The default value for this option is false.
 *
 * 5) "tag"="name", "tagval"="value"
 * ---------------------------------
 * The "tag" and "tagval" options can be used to define a subset of the mesh
 * elements in the input file that will be tallied.  These tags are defined
 * on the input mesh itself using the MOAB tagging feature.  Note that "tag"
 * name can only be set once, whereas multiple "tagval" values can be added.
 * This option is only used during setup to define the set of tally points.
 */
//===========================================================================//
class TrackLengthMeshTally : public MeshTally
{ 
  public:
    /**
     * \brief Constructor
     * \param[in] input user-defined input parameters for this TrackLengthMeshTally
     */
    explicit TrackLengthMeshTally(const TallyInput& input);

    /**
     * \brief Virtual destructor
     */
    virtual ~TrackLengthMeshTally();

    // >>> DERIVED PUBLIC INTERFACE from Tally.hpp

    /**
     * \brief Computes scores for this TrackLengthMeshTally based on the given TallyEvent
     * \param[in] event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event);

    /**
     * \brief Updates TrackLengthMeshTally when a particle history ends
     *
     * Calls MeshTally::end_history() and if input mesh is conformal sets the
     * last_cell to -1 to indicate that a new particle will be born.
     */
    virtual void end_history();

    /**
     * \brief Write results to the output file for this TrackLengthMeshTally
     * \param[in] num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results for all of the mesh cells defined as tally points to the
     * output_filename set for this TrackLengthMeshTally.  These values are
     * normalized by both the number of particle histories that were tracked
     * and the volume of the mesh cell for which the results were computed.
     */
    virtual void write_data(double num_histories);

  protected:
    /// Copy constructor and operator= methods are not implemented
    TrackLengthMeshTally(const TrackLengthMeshTally& obj);
    TrackLengthMeshTally& operator=(const TrackLengthMeshTally& obj);

  protected:
    // MOAB instance that stores all of the mesh data
    moab::Interface* mb;

    // KD-Tree used with unstructured meshes to get starting tetrahedron
    moab::AdaptiveKDTree* kdtree;  
    moab::EntityHandle kdtree_root;

    // Oriented Box Tree variables
    OrientedBoxTreeTool* obb_tool;
    EntityHandle obbtree_root;

    // Variables needed to keep track of mesh cells visited
    EntityHandle last_visited_tet;
    int last_cell;
    
    // Optional convex mesh and conformal surface source flags
    bool convex;
    bool conformal_surface_source;

    // If not empty, user has asserted mesh tally geometry
    // conforms to the cells identified in this set
    std::set<int> conformality;

    // Stores barycentric data for tetrahedrons
    std::vector<Matrix3> tet_baryc_data;

    // >>> PROTECTED METHODS

    /**
     * \brief Parse the TallyInput options for this TrackLengthMeshTally
     */
    void parse_tally_options();

    /**
     * \brief Loads and sets the mesh data used by this TrackLengthMeshTally
     *
     * Mesh data is loaded from the input_filename provided in the TallyInput
     * options.  If no tag name and values are found, the set of tally points
     * will be all of the mesh cells defined in input_filename.  Otherwise,
     * the set of tally points stored in MeshTally::tally_mesh_set will only
     * contain the mesh cells that have the given tag name and tag values.
     */
    void set_tally_meshset();

    /**
     * \brief Computes the barycentric matrices for all tetrahedrons
     * \param[in] all_tets the set of tets extracted from the input mesh
     * \return the MOAB ErrorCode value
     */
    ErrorCode compute_barycentric_data (const Range& all_tets);

    /**
     * \brief Constructs the KD and OBB trees from the mesh data
     * \param[in, out] all_tets the set of tets extracted from the input mesh
     *
     * Also adds the set of skin triangles to all_tets.
     */
    void build_trees (Range& all_tets);

    /**
     * \brief returns all triangle intersections and their distances form the point position
     * \param[in] position cart vect of the origin
     * \param[in] direction cart vect unit vector of direction
     * \param[in] track_length distance to intersect upto 
     * \param[out] triangles returns the handles to the triangles we hit
     * \param[out] intersections returns the distances to the triangles we hit
     * \returns the error code should call the kdtree fail
     */
    ErrorCode get_all_intersections(const CartVect& position, const CartVect& direction, double track_length, 
				    std::vector<EntityHandle> &triangles,std::vector<double> &intersections);

    /**
     * \brief Checks if the given point is inside the given tet
     * \param[in] point the coordinates of the point to test
     * \param[in] tet the EntityHandle of the tet
     * \return true if the point falls inside tet; false otherwise
     *
     * Assumes that tet is part of this TrackLengthMeshTally.
     */                 
    bool point_in_tet(const CartVect& point, const EntityHandle* tet);

  /**
   * \brief loop through all tets to find which tet, the point belong to
   * \param [in] point point to test
   * \return entity handle to the tet which the point belongs to
   */
  EntityHandle point_in_which_tet (const CartVect& point);

  /**
   * \brief return the tet_element in which the ray ends
   * \param[in] start the begining of the ray
   * \param[in] dir unit direction vector of the ray
   * \param[in] distance to the last intersecction
   * \param[in] left_over remaining track_length
   * \return the tet which the end point of the ray belongs to, zero if none found
   */
  EntityHandle remainder(const CartVect& start, const CartVect& dir, double distance, double left_over);

  /** 
   * \brief return the MBRange of triangles that belong in the tet mesh there are no repliacted surfaces in this list
   * \param[in] all_tets the MBRange of tets in the problem
   * \return range of triangles that make up faces in the mesh
   */
  Range get_adjacency_info(const Range& all_tets);

  /** 
   * \brief return the sorted (by distance) list of triangles and intersections
   * \param[in, out] vector<double> intersections list of all the intersections
   * \param[in, out] vector<EntityHandle> intersections list of the triangle entity handles that correspond to the intersections
   * \return void
   */
  void sort_intersection_data(std::vector<double> &intersections, std::vector<EntityHandle> &triangles);

  /** 
   * \brief return the tracklengths of the ray in each tet
   * \param[in] event the tally event, direction, position, track_length, etc
   * \param[in] ebin the energy bin index corresponding to the energy
   * \param[in] weight the multiplier value for the score to be tallied
   * \param[in] vector<double> intersections list of all the intersections
   * \param[in] vector<EntityHandle> triangles list of the triangle entity handles that correspond to the intersections
   * \return void
   */
  void compute_tracklengths(const TallyEvent& event, 
                            unsigned int ebin, double weight,
   			    const std::vector<double>& intersections,
			    const std::vector<EntityHandle>& triangles);

};

} // end namespace moab

#endif // DAGMC_TRACK_LENGTH_MESH_TALLY_HPP

// end of MCNP5/dagmc/TrackLengthMeshTally.hpp
