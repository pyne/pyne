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
     * \param input user-defined input parameters for this KDEMeshTally
     */
    explicit TrackLengthMeshTally(const TallyInput& input);

    /**
     * \brief Virtual destructor
     */
    virtual ~TrackLengthMeshTally();

    // >>> DERIVED PUBLIC INTERFACE from Tally.hpp

    /**
     * \brief Computes scores for this TrackLengthMeshTally based on the given TallyEvent
     * \param event the parameters needed to compute the scores
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
     * \param num_histories the number of particle histories tracked
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

    // Counts how many negative tracks occurred that were not tallied
    int num_negative_tracks;

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
     * \param all_tets the set of tets extracted from the input mesh
     * \return the MOAB ErrorCode value
     */
    ErrorCode compute_barycentric_data (const Range& all_tets);

    /**
     * \brief Constructs the KD and OBB trees from the mesh data
     * \param all_tets the set of tets extracted from the input mesh
     *
     * Also adds the set of skin triangles to all_tets.
     */
    void build_trees (Range& all_tets);


    /**
     * \brief returns all triangle intersections and their distances form the point position
     * \param position cart vect of the origin
     * \param direction cart vect unit vector of direction
     * \param track_length distance to intersect upto 
     * \param triangles returns the handles to the triangles we hit
     * \param intersections returns the distances to the triangles we hit
     * \returns the error code should call the kdtree fail
     * Also adds the set of skin triangles to all_tets.
     */
    ErrorCode get_all_intersections(CartVect position, CartVect direction, double track_length, 
				    std::vector<EntityHandle> &triangles,std::vector<double> &intersections);

    /**
     * \brief Checks if the given point is inside the given tet
     * \param point the coordinates of the point to test
     * \param tet the EntityHandle of the tet
     * \return true if the point falls inside tet; false otherwise
     *
     * Assumes that tet is part of this TrackLengthMeshTally.
     */                 
    bool point_in_tet(const CartVect& point, const EntityHandle* tet);

    /**
     * \brief Get parent tetrahedron and vertices from given triangle
     * \param triangle the input triangle
     * \param tetrahedron (output) the parent tetrahedron of triangle
     * \param vertices (output) the verticles of triangle
     *
     * Given an EntityHandle representing a mesh skin triangle, return the
     * tetrahedron it belongs to as well as the handles of its vertices.
     * Calling this function asserts that the triangle "belongs" to exactly
     * one tetrahedron.  This should always be a valid assumption, since
     * this code only creates explicit triangle entities for skin triangles.
     */
    void get_skin_triangle_adjacencies(EntityHandle triangle, 
                                       EntityHandle& tetrahedron,
                                       EntityHandle vertices[3]);

    /**
     * \brief Finds the next tetrahedron (if any) intersected by a ray
     * \param start the beginning of ray
     * \param vec the unit vector of ray
     * \param length the length of ray
     * \param first_tri if return value is non-zero, will contain the skin
     *                  triangle that ray intersects
     * \param first_t value of t along the ray where first_tri is intersected
     * \param last_crossing a skin triangle to ignore for this ray fire (Used
     *                      if a ray is being fired away from surface of mesh)
     * \return tetrahedron intersected by this ray fire, or zero if none
     *
     * Given a ray that starts outside the mesh, find the next tetrahedron
     * (if any) that will be intersected by the ray.
     */
    EntityHandle find_next_tet_by_ray_fire(const CartVect& start,
                                           const CartVect& vec, double length, 
                                           EntityHandle first_tri[3], double& first_t, 
                                           EntityHandle last_crossing = 0);

    /**
     * \brief Find first tet intersected by a track, using conformal logic
     * \param start the beginning of track
     * \param first_tri (output) vertices of skin triangle at which track begins
     * \return first tetrahedron intersected
     *
     * Find the first tetrahedron intersected by a track segment, using
     * conformal logic.  This is invoked when a surface crossing has just
     * occurred in a conformal mesh tally, allowing us to assume that one of
     * the surface triangles of this mesh has just been crossed.  Thus the
     * search for the nearest tetrahedron is simply a search for the surface
     * triangle nearest to the starting point of the ray.
     *
     * This function should always successfully return values, but it will
     * print a warning if the nearest available result is suspicious.
     */
    EntityHandle get_starting_tet_conformal(const CartVect& start,
                                            EntityHandle first_tri[3]);

    /**
     * \brief Find first tet intersected by the given track segment
     * \param start the beginning of track
     * \param vec the unit vector of track
     * \param length the track length
     * \param first_tri if return value is non-zero, and ray fire was performed,
     *                  will contain the skin triangle that ray intersects
     * \param first_t value of t along the ray where first_tri intersects
     * \return the first tetrahedron along the ray, or zero if none found
     */ 
    EntityHandle get_starting_tet(const CartVect& start,
                                  const CartVect& vec, double length, 
                                  EntityHandle first_tri[3], double& first_t);

  /**
   * \brief loop through all tets to find which tet, the point belong to
   * \param point point to test
   * \return entity handle to the tet which the point belongs to
   */
  EntityHandle point_in_which_tet (CartVect point);

  /**
   * \brief return the tet_element in which the ray ends
   * \param start the begining of the ray
   * \param dir unit direction vector of the ray
   * \param distance to the last intersecction
   * \param left_over remaining track_length
   * \return the tet which the end point of the ray belongs to, zero if none found
   */
  EntityHandle remainder( CartVect start, CartVect dir, double distance, double left_over);

  /** 
   * \brief return the MBRange of triangles that belong in the tet mesh there are no repliacted surfaces in this list
   * \param all_tets the MBRange of tets in the problem
   * \return range of triangles that make up faces in the mesh
   */
  Range get_adjacency_info(Range all_tets);

  /** 
   * \brief return the sorted (by distance) list of triangles and intersections
   * \param vector<double> intersections list of all the intersections
   * \param vector<EntityHandle> intersections list of the triangle entity handles that correspond to the intersections
   * \return void
   */
  void sort_intersection_data(std::vector<double> &intersections, std::vector<EntityHandle> &triangles);

  /** 
   * \brief return the tracklengths of the ray in each tet
   * \param vector<double> intersections list of all the intersections
   * \param vector<EntityHandle> intersections list of the triangle entity handles that correspond to the intersections
   * \return void
   */
  void compute_tracklengths(const TallyEvent event,
			    const std::vector<double> intersections,
			    const std::vector<EntityHandle> triangles);

  /** 
   * \brief determine the score to add to the results
   * \param tally event the tally event, direction, position etc
   * \param double tracklength, the physical tracklength of the current history
   * \param EntityHandle tet, the MOAB entityhandle of the current tet
   * \return void
   */
  void determine_score(const TallyEvent event, double tracklength, EntityHandle tet);

};

} // end namespace moab

#endif // DAGMC_TRACK_LENGTH_MESH_TALLY_HPP

// end of MCNP5/dagmc/TrackLengthMeshTally.hpp
