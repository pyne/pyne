// MCNP5/dagmc/KDEMeshTally.hpp

#ifndef DAGMC_KDE_MESH_TALLY_HPP
#define DAGMC_KDE_MESH_TALLY_HPP

#include <utility>
#include <vector>

#include "moab/CartVect.hpp"

#include "KDEKernel.hpp"
#include "KDENeighborhood.hpp"
#include "MeshTally.hpp"
#include "Quadrature.hpp"

// forward declarations
namespace moab {
  class Interface;
}

//===========================================================================//
/**
 * \class KDEMeshTally
 * \brief Represents a mesh tally based on a kernel density estimator
 *
 * KDEMeshTally is a concrete class derived from MeshTally that implements
 * the Tally interface to enable the use of the kernel density estimator as a
 * mesh tally option in a Monte Carlo particle transport simulation.  Three
 * different types of KDE mesh tallies can be created
 *
 *     1) KDE COLLISION mesh tally, based only on particle collisions
 *     2) KDE INTEGRAL_TRACK mesh tally, based on full particle tracks
 *     3) KDE SUB_TRACK mesh tally, based on approximated particle tracks
 *
 * The primary difference between these three options is what estimator is used
 * to compute the tally scores for the TallyEvent::COLLISION or TallyEvent::TRACK
 * events that are set through the TallyManager.  If a KDEMeshTally object
 * receives an invalid TallyEvent type for its estimator type, then no scores
 * are computed.
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
 * ==========
 * TallyInput
 * ==========
 *
 * The TallyInput struct needed to construct a KDEMeshTally object is defined
 * in Tally.hpp and is set through the TallyManager when a Tally is created.
 * Options that are currently available for KDEMeshTally objects include
 *
 * 1) "inp"="input_filename", "out"="output_filename"
 * --------------------------------------------------
 * These two options are processed through the MeshTally constructor.  The
 * "inp" key is REQUIRED for KDEMeshTally objects, whereas the "out" key is
 * optional.  See MeshTally.hpp for more information.
 *
 * 2) "hx"="value", "hy"="value", "hz"="value"
 * -------------------------------------------
 * Sets the value of the bandwidth vector components.  Note that these values
 * must all be positive and non-zero.  The default value is (0.01, 0.01, 0.01).
 *
 * 3) "kernel"="type", "order"="value"
 * -----------------------------------
 * Defines the type and order of the kernel function to be used.  Valid types
 * include "epanechnikov", "uniform", "biweight", or "triweight".  Note that
 * only symmetric versions of these kernels are supported, which means that
 * "order" must be a multiple of 2.  If "kernel" is omitted or invalid, then
 * the default is "epanechnikov".  Similarly, if "order" is omitted or invalid,
 * then the default is 2nd-order.
 *
 * 4) "neighborhood"="off"
 * -----------------------
 * Turns off the neighborhood-search and computes scores for all calculation
 * points.  This key will also be used to request different neighborhood-search
 * methods in the future.  The default is the kd-tree method.
 *
 * 5) "boundary"="default"
 * -----------------------
 * Indicates that boundary correction is needed for tally points within one
 * bandwidth of an external boundary.  The "default" method uses the boundary
 * kernel approach, and is currently the only option available.  Note that
 * this feature will only work properly for 2nd-order kernels.
 *
 * 6) "seed"="value", "subtracks"="value"
 * --------------------------------------
 * These two options are only available for KDE sub-track mesh tallies.  The
 * "seed" option overrides the random number seed value that is used for
 * determining sub-track points.  The "subtracks" option sets the number
 * of sub-tracks to use for computing scores.
 */
//===========================================================================//
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
     * \param[in] input user-defined input parameters for this KDEMeshTally
     * \param[in] type the type of estimator to use with this KDEMeshTally
     */
    KDEMeshTally(const TallyInput& input, Estimator type = COLLISION);

    /**
     * \brief Virtual destructor
     */
    virtual ~KDEMeshTally();

    // >>> DERIVED PUBLIC INTERFACE from Tally.hpp

    /**
     * \brief Computes scores for this KDEMeshTally based on the given TallyEvent
     * \param[in] event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event);

    /**
     * \brief Write results to the output file for this KDEMeshTally
     * \param[in] num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results for all of the mesh nodes to the output_filename set for
     * this KDEMeshTally.  These values are normalized only by the number of
     * particle histories that were tracked.
     */
    virtual void write_data(double num_histories);

  private:
    // Copy constructor and operator= methods are not implemented
    KDEMeshTally(const KDEMeshTally& obj);
    KDEMeshTally& operator=(const KDEMeshTally& obj);

  private:
    // Type of estimator used by this KDE mesh tally
    Estimator estimator;

    // Bandwidth vector (hx, hy, hz) used to compute KDE mesh tally scores
    moab::CartVect bandwidth;

    // Kernel function used to compute KDE mesh tally scores
    KDEKernel* kernel;

    // Defines neighborhood region for computing scores
    bool use_kd_tree;
    KDENeighborhood* region;

    // Variables used if boundary correction method is requested by user
    bool use_boundary_correction;
    moab::Tag boundary_tag;
    moab::Tag distance_tag;

    // Number of sub-tracks used to compute KDE sub-track mesh tally scores
    unsigned int num_subtracks;

    // Quadrature used to compute KDE integral-track mesh tally scores
    Quadrature* quadrature;

    // MOAB instance that stores all of the mesh data
    moab::Interface* mbi;

    // Running variance variables for computing optimal bandwidth at runtime
    bool max_collisions;
    long long int num_collisions;       
    moab::CartVect mean;
    moab::CartVect variance;

    // If true, another instance already set the random number generator seed
    static bool seed_is_set;

    // >>> PRIVATE METHODS

    /**
     * \brief Sets bandwidth[i] = value for the given key
     * \param[in] key the name of the bandwidth vector component
     * \param[in] value the value of the bandwidth vector component
     * \param[in] i the index in the bandwidth vector associated with the key
     */
    void set_bandwidth_value(const std::string& key,
                             const std::string& value,
                             unsigned int i);

    /**
     * \brief Parse the TallyInput options for this KDEMeshTally
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
     * \param[in] collision_point the coordinates of the collision point
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

    // Declare test fixture as friend for testing KDE estimator methods
    friend class KDEMeshTallyTest;

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
         * \param[in] tally the KDEMeshTally object using this path kernel
         * \param[in] event the tally event containing the track segment data
         * \param[in] X the calculation point
         */
        PathKernel(const KDEMeshTally& tally,
                   const TallyEvent& event,
                   const CalculationPoint& X)
            : tally(tally), event(event), X(X) {}

        /**
         * \brief Evaluates the path length kernel function
         * \param[in] s the path length for which to evaluate this function
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
     * \param[in] X the calculation point
     * \param[in] observation the random observation point (Xi, Yi, Zi)
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
     * \param[in] X the calculation point
     * \param[in] event the tally event containing the track segment data
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
     * \param[in] event the tally event containing the track segment data
     * \param[in] coords the (x, y, z) coordinates of the calculation point X
     * \param[out] limits stores integration limits in form of pair<lower, upper>
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
     * \param[in] X the calculation point
     * \param[in] points the set of sub-track points needed for computing score
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
     * \param[in] p the number of random points requested
     * \param[in] event the tally event containing the track segment data
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
