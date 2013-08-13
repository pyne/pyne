// MCNP5/dagmc/KDEMeshTally.hpp

#ifndef DAGMC_KDE_MESH_TALLY_HPP
#define DAGMC_KDE_MESH_TALLY_HPP

#include <utility>
#include <vector>

#include "moab/CartVect.hpp"

#include "KDEKernel.hpp"
#include "MeshTally.hpp"
#include "TallyEvent.hpp"
#include "Quadrature.hpp"

// forward declarations
namespace moab {
  class AdaptiveKDTree;
  class Interface;
}

/**
 * \class KDEMeshTally
 * \brief Represents a mesh tally based on a kernel density estimator approach
 *
 * KDEMeshTally is a Derived class that implements the MeshTally interface to
 * enable the use of the kernel density estimator as a mesh tally option in a
 * Monte Carlo particle transport simulation.  Three different types of KDE
 * mesh tallies can be created
 *
 *     1) KDE COLLISION mesh tally, based only on particle collisions
 *     2) KDE INTEGRAL_TRACK mesh tally, based on full particle tracks
 *     3) KDE SUB_TRACK mesh tally, based on approximated particle tracks
 *
 * The primary difference between these three options is what estimator is used
 * to compute the tally scores for the collisions or particle tracks that are
 * received from the Monte Carlo particle transport code.  Note that this cannot
 * be changed once a KDEMeshTally object has been created.
 *
 * For more information on using kernel density estimator methods for mesh tally
 * purposes, please refer to the following two references
 *
 *     1) K. L. Dunn and P. P. H. Wilson, "Kernel Density Estimators for Monte
 *        Carlo Tallies on Unstructured Meshes," Transactions of the American
 *        Nuclear Society, 107, pp. 490-493 (2012)
 *
 *     2) K. L. Dunn and P. P. H. Wilson, "Monte Carlo Mesh Tallies based on a
 *        Kernel Density Estimator Approach using Integrated Particle Tracks,"
 *        Proceedings of the International Conference on Mathematics and
 *        Computational Methods (M&C 2013), Sun Valley, Idaho, May 5-9 (2013)
 *
 * ==============
 * TallyInput
 * ==============
 *
 * In addition to choosing the type of estimator, the user must also prepare a
 * TallyInput object that will be used to set all of the other required and
 * optional input parameters.  This is a struct defined in MeshTally.hpp that
 * must include the tally ID number, the energy bin boundaries, whether or not
 * a total energy bin is required, the name of the mesh input file, and a
 * std::map with all other tally options as key-value pairs.  Tally options
 * that are currently available include the following
 *
 *     1) "out": defines the name of the KDE mesh tally output file
 *     2) "hx", "hy", "hz": sets the value of the bandwidth vector components
 *     3) "kernel": defines the kernel function to be used
 *     4) "seed": overrides random number seed for determining sub-track points
 *     5) "subtracks": sets number of sub-tracks to use
 *
 * The last two options are only valid for KDE sub-track mesh tallies.  Note
 * also that if an invalid "kernel" is requested, or this key is not included
 * in the MeshTallyInput options, then the default 2nd-order "epanechnikov"
 * kernel function is used.  Other options include "uniform", "biweight",
 * or "triweight".
 *
 * TODO need to add higher-order kernels/boundary kernel as available keys in
 * tally options as well
 *
 * ===============
 * Mesh Input File
 * ===============
 *
 * The mesh input file that is added to MeshTallyInput must be created in a file
 * format that is supported by the Mesh-Oriented Database (MOAB), which includes
 * both H5M and VTK options.  Source code and more information on MOAB can be
 * found at https://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB
 *
 * ============================
 * KDE Mesh Tally Functionality
 * ============================
 *
 * After a KDEMeshTally object has been created, there are a few simple steps
 * that are needed to ensure that the tally will work correctly.  The first step
 * is to use the compute_score() method for every tally event that is reported
 * by the Monte Carlo particle transport code being used.  These tally events
 * can be either collision or track-based, as the compute_score() method will
 * use the type of estimator that has been associated with the KDEMeshTally
 * object for computing the corresponding tally scores.  These scores are
 * accumulated in the MeshTally::temp_tally_data array.
 *
 * The second step is to use the end_history() method once all tally scores for
 * a single particle history has been completed.  This step takes the current
 * sum of scores that is in MeshTally::temp_tally_data and adds it to the
 * MeshTally::tally_data array.  It is also important for calculating the
 * relative standard error of the final tally results.
 *
 * The third and final step is to use the write_data() method to write all of the
 * tally results and their corresponding relative standard errors to the output
 * file.  If the "out" key is not included as a tally option, then the default
 * case is a H5M file format named meshtal<tally_id>.h5m.  To write to another
 * file format that is supported by MOAB, simply add the "out" key and use a
 * name with the desired extension.  For example, to write to a VTK file format,
 * add "out" = "file_name.vtk" to the MeshTallyInput options.
 */
class KDEMeshTally : public MeshTally
{
  public:
    /**
     * \brief Defines type of estimator used in computing KDE mesh tally scores
     *
     *     0) COLLISION estimator is used for KDE collision tallies
     *     1) INTEGRAL_TRACK estimator is used for KDE integral-track tallies
     *     2) SUB_TRACK estimator is used for KDE sub-track tallies
     */
    enum Estimator {COLLISION = 0, INTEGRAL_TRACK = 1, SUB_TRACK = 2};
    static const char* const kde_estimator_names[];

    /**
     * \brief Constructor
     * \param input user-defined input parameters for this KDE mesh tally
     * \param type the type of estimator to use with this KDE mesh tally
     */
    KDEMeshTally(const TallyInput& input, Estimator type = COLLISION);

    /**
     * \brief Virtual destructor
     */
    virtual ~KDEMeshTally();

    // >>> DERIVED PUBLIC INTERFACE from MeshTally.hpp

    /**
     * \brief Computes mesh tally scores for the given tally event
     * \param event the parameters needed to compute the mesh tally scores
     */
    virtual void compute_score(const TallyEvent& event, int ebin);

    /**
     * \brief Write tally and error results to the mesh tally's output file
     * \param num_histories the number of source particles tracked
     * \param multiplier an optional constant multiplication factor
     *
     * The write_data() method writes the current tally and relative standard error
     * results to the output file defined for this mesh tally, normalized by
     * the number of source particles that were tracked during the Monte Carlo
     * simulation.
     *
     * If the multiplier parameter is provided, then all final tally results
     * will be multiplied by this value.  Note that this multiplier is not the
     * same as the standard tally multiplier, which is typically applied to
     * individual scores instead.
     */
    virtual void write_data(double num_histories);

  private:
    /// Copy constructor and operator= methods are not implemented
    KDEMeshTally(const KDEMeshTally& obj);
    KDEMeshTally& operator=(const KDEMeshTally& obj);

  private:
    /// Type of estimator used by this KDE mesh tally
    Estimator estimator;

    /// Bandwidth vector (hx, hy, hz) used to compute KDE mesh tally scores
    moab::CartVect bandwidth;

    /// Kernel function used to compute KDE mesh tally scores
    KDEKernel* kernel;

    /// Variables used if boundary correction method is requested by user
    bool use_boundary_correction;
    moab::Tag boundary_tag;
    moab::Tag distance_tag;

    /// Number of sub-tracks used to compute KDE sub-track mesh tally scores
    unsigned int num_subtracks;

    /// Quadrature used to compute KDE integral-track mesh tally scores
    Quadrature* quadrature;

    /// MOAB instance that stores all of the mesh data
    moab::Interface* mbi;

    /// KD-Tree used with unstructured meshes to locate calculation points
    moab::AdaptiveKDTree* kd_tree;  
    moab::EntityHandle kd_tree_root;

    /// Running variance variables for computing optimal bandwidth at runtime
    bool max_collisions;
    long long int num_collisions;       
    moab::CartVect mean;
    moab::CartVect variance;

    /// If true, another instance already set the random number generator seed
    static bool seed_is_set;

    // >>> PRIVATE METHODS

    /**
     * \brief Sets bandwidth[i] = value for the given key
     * \param key the name of the bandwidth vector component
     * \param value the value of the bandwidth vector component
     * \param i the index in the bandwidth vector associated with the key
     */
    void set_bandwidth_value(const std::string& key,
                             const std::string& value,
                             unsigned int i);

    /**
     * \brief Parse the TallyInput options for this KDE mesh tally
     */
    void parse_tally_options();

    /**
     * \brief Initializes MeshTally member variables representing the mesh data
     * \return the MOAB ErrorCode value
     *
     * MeshTally variables set by this method include tally_points and
     * tally_mesh_set.  The tally_points will be defined as the set of mesh
     * nodes, whereas tally_mesh_set stores the set of all 3D mesh elements.
     * This method also calls MeshTally::setup_tags() to set the tag names for
     * the energy bins, and initializes tag handles for boundary correction.
     */
    moab::ErrorCode initialize_mesh_data();

    /**
     * \brief Adds the collision point to the running variance formula
     * \param collision_point the coordinates of the collision point
     *
     * The update_variance() method updates mean and variance variables with
     * the coordinates of the new collision point, which can then be used by
     * get_optimal_bandwidth() to compute the optimal bandwidth vector.
     */
    void update_variance(const moab::CartVect& collision_point);

    /**
     * \brief Computes the optimal bandwidth vector
     * \return h_optimal = (hx_optimal, hy_optimal, hz_optimal)
     *
     * The get_optimal_bandwidth() method determines what the bandwidth should
     * be to achieve good results using the Kernel Density Estimator.  For the
     * 3D case it is defined globally by the following formula
     *
     *     h_optimal[i] = 0.968625 * stdev[i] * N^(-1.0/7.0)
     *
     * where stdev[i] is the standard deviation of the ith component of the
     * collision point locations, and N is the number of collisions that have
     * occurred.
     *
     * NOTE: stdev is calculated during runtime by calling update_variance()
     * whenever a collision occurs.
     */
    moab::CartVect get_optimal_bandwidth() const;
  
    // >>> KDE ESTIMATOR METHODS

    // Declare test fixtures as friends for testing KDE estimator methods
    friend class KDEIntegralTrackTest;
    friend class KDESubtrackTest;
    friend class KDECollisionTest;

    // Defines common data needed for computing score for a calculation point
    struct CalculationPoint
    {
        double coords[3];
        int boundary_data[3];
        double distance_data[3];
    };

    /**
     * \class PathKernel
     * \brief Defines the path length kernel function K(X, s)
     */
    class PathKernel : public Function
    {
      public:
        /**
         * \brief Constructor
         * \param tally the KDEMeshTally object using this path kernel
         * \param event the tally event containing the track segment data
         * \param X the calculation point
         */
        PathKernel(const KDEMeshTally& tally,
                   const TallyEvent& event,
                   const CalculationPoint& X)
            : tally(tally), event(event), X(X) {}

        /**
         * \brief Evaluates the path length kernel function
         * \param s the path length for which to evaluate this function
         * \return K(X, s)
         */
        double evaluate(double s) const;

      private:
        const KDEMeshTally& tally;
        const TallyEvent& event;
        const CalculationPoint& X;
    };

    /**
     * \brief Computes value of the 3D kernel function K(x, y, z)
     * \param X the calculation point
     * \param observation the random observation point (Xi, Yi, Zi)
     * \return K(x, y, z)
     *
     * Evaluates the general 3D kernel function for the given calculation and
     * observation points.  Uses the bandwidth vector stored in this KDEMeshTally
     * and adjusts for the boundary correction as needed.
     */
    double evaluate_kernel(const CalculationPoint& X,
                           const moab::CartVect& observation) const;                    

    /**
     * \brief Computes tally score based on the integral-track estimator
     * \param X the calculation point
     * \param event the tally event containing the track segment data
     * \return the tally score for the calculation point
     *
     * The integral_track_score() method computes the integral of the 3D
     * path-length dependent kernel function K(X, s) with respect to path-
     * length s for the given calculation point X, using the limits of
     * integration as determined by the set_integral_limits() method.
     */
    double integral_track_score(const CalculationPoint& X,
                                const TallyEvent& event) const;

    /**
     * \brief Determines integration limits for the integral-track estimator
     * \param event the tally event containing the track segment data
     * \param coords the (x, y, z) coordinates of the calculation point X
     * \param limits stores integration limits in form of pair<lower, upper>
     * \return true if valid limits exist, false otherwise
     *
     * The set_integral_limits() method determines the integration limits for
     * the integral of the 3D path-length dependent kernel function K(X, s)
     * with respect to path-length s for the given calculation point X.
     *
     * NOTE: if this method returns false, that means there are no valid
     * integration limits within the range from 0 to the total track length.
     * This essentially means that the resulting integrand over this range
     * would have been zero and the calculation point can be ignored.
     */
    bool set_integral_limits(const TallyEvent& event,
                             const moab::CartVect& coords,
                             std::pair<double, double>& limits) const;

    /**
     * \brief Computes tally score based on the sub-track estimator
     * \param X the calculation point
     * \param points the set of sub-track points needed for computing score
     * \return the tally score for the calculation point
     *
     * The subtrack_score() method computes the average 3D kernel contribution
     * for the number of sub-tracks requested for the given calculation point X,
     * using the randomly chosen points along the track that were computed by
     * the choose_points() method.
     */
    double subtrack_score(const CalculationPoint& X,
                          const std::vector<moab::CartVect>& points) const;

    /**
     * \brief Chooses p random points along a track segment
     * \param p the number of random points requested
     * \param event the tally event containing the track segment data
     * \return the vector of p random points
     *
     * The choose_points() method sub-divides the track segment into p
     * sub-tracks of equal length and randomly chooses the coordinates of
     * one point from each sub-track.
     */
    std::vector<moab::CartVect> choose_points(unsigned int p,
                                              const TallyEvent& event) const;
};

#endif // DAGMC_KDE_MESH_TALLY_HPP

// end of MCNP5/dagmc/KDEMeshTally.hpp
