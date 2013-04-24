// KDEMeshTally.hpp

#ifndef KDEMESHTALLY_H
#define KDEMESHTALLY_H

#include <iosfwd>
#include <set>
#include <string>
#include <utility>

#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"

#include "KDEKernel.hpp"
#include "MeshTally.hpp"
#include "TallyEvent.hpp"

// forward declarations
class AdaptiveKDTree;

/** 
 * A class that represents a mesh tally based on a kernel density estimator
 * approach for use in a Monte Carlo particle transport simulation.
 */
class KDEMeshTally : public MeshTally {

  public:

     static moab::CartVect default_bandwidth;  // bandwidth choice for KDE tallies

    /**
     * \brief Defines type of estimator used in computing KDE mesh tally scores
     *
     *     1) COLLISION tallies use the KDE collision estimator
     *     2) INTEGRAL_TRACK tallies use the KDE integral-track estimator
     *     3) SUB_TRACK tallies use the KDE sub-track estimator
     */
    enum TallyType { COLLISION = 1, INTEGRAL_TRACK = 2, SUB_TRACK = 3 };

    /**
     * Allocate and return the specified tally object
     * @param settings The FC card parameters
     * @param mbi The associated MOAB interface
     * @param type The tally type
     * @return A newly allocated object, ready to receive collisions or tracks 
     */
    static KDEMeshTally* setup( const MeshTallyInput& settings, 
                                moab::Interface* mbi,
                                TallyType type = COLLISION );

    /**
     * Constructs a mesh tally object based on a kernel density estimator
     * approach.  The tally type cannot be changed once a KDEMeshTally object
     * has been created.
     *
     * NOTE: the parameter "numSubtracks" is required for SUB_TRACK tallies,
     * but it is optional for INTEGRAL_TRACK tallies for computing an optimal
     * bandwidth.  It has no meaning for COLLISION tallies.
     *
     * @param settings the FC card parameters 
     * @param moabMesh the MOAB instance containing the mesh information
     * @param moabSet the MOAB set of entities used in the tally
     * @param bandwidth the set of bandwidth values (hx, hy, hz)
     * @param type (optional) the type of tally to be used
     * @param k (optional) the KernelType function to be used in the computation
     * @param numSubtracks (optional) the number of subtracks to be used
     */
    KDEMeshTally( const MeshTallyInput& settings,
                  moab::Interface* moabMesh,
                  moab::EntityHandle moabSet,
                  moab::CartVect bandwidth,
                  TallyType type = COLLISION,
                  KDEKernel::KernelType k = KDEKernel::EPANECHNIKOV,
                  unsigned int numSubtracks = 0 );

    /**
     * Destructor.
     */
    ~KDEMeshTally();

    /**
     * \brief Computes mesh tally scores for the given tally event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin index representing energy bin
     * TODO remove ebin as a parameter
     */
    void compute_score(const TallyEvent& event, int ebin);

    /** 
     * Implement MeshTally end_history interface 
     */ 
    virtual void end_history();

    /** 
     * Implement MeshTally print interface 
     */ 
    virtual void print( double sp_norm, double fmesh_fact );

    /**
     * Write mesh with tags for current tally and error results attached.
     *
     * @param sp_norm the number of source particles
     * @param mult_fact the multiplication factor from the FMESH card
     */
    void write_results( double sp_norm, double fmesh_fact );

  private:

    moab::Interface* mb;           // the MOAB instance
    
    moab::CartVect bandwidth;      // the set of bandwidth values (hx, hy, hz)
    TallyType kde_tally;           // specifies type of KDE tally 
    KDEKernel* kernel;             // kernel function tally is based on
    unsigned int num_subtracks;    // number of subtracks used in tally       
    
    moab::EntityHandle tally_set;  // the MOAB set of entities used in tally

    moab::AdaptiveKDTree* tree;    
    moab::EntityHandle tree_root;
    
    // The entityhandles updated in the current particle history; cleared by end_history()
    std::set<moab::EntityHandle> visited_this_history;

    /**
     * Variables needed to compute a running variance for calculating the
     * optimal bandwidth during runtime.
     */
    bool max_collisions;
    long long int numCollisions;       
    moab::CartVect Mn;  // mean
    moab::CartVect Sn;  // variance

    // subtrack estimator parameters
    static bool seed_is_set;
        
    /**
     * Builds a KD-Tree from all the mesh node entities in the specified meshset,
     * or if the set is 0, from all mesh nodes in MOAB.
     */
    void build_tree( moab::EntityHandle meshset );

    /**
     * Adds the given collision coordinates to the running variance formula
     * used to compute the optimal bandwidth.
     */
    void update_variance(const moab::CartVect& collision_point);

    /**
     * Computes and returns the optimal bandwidth to achieve the best results
     * in the kernel density approximation.
     */
    moab::CartVect get_optimal_bandwidth();
   
    /**
     * Adds a tally contribution to one calculation point on the mesh.
     *
     * @param mesh_point the calculation point on the mesh
     * @param score the tally contribution to add
     * @param ebin the energy bin for the tally contribution
     */
    void add_score_to_tally( moab::EntityHandle mesh_point,
                             double score,
                             int ebin );

    /// unimplemented 
    KDEMeshTally( const KDEMeshTally & obj );
    /// unimplemented
    KDEMeshTally & operator=( const KDEMeshTally & obj );

    // >>> KDE ESTIMATOR FUNCTIONS

    /**
     * \brief Computes tally score based on the integral-track estimator
     * \param data the physics data for the track segment being tallied
     * \param X the (x, y, z) coordinates of the calculation point
     * \return the tally score for the calculation point
     *
     * The integral_track_score() function computes the integral of the
     * 3D path-length dependent kernel function K(X, s) with respect to path-
     * length s for the given calculation point X, using the limits of
     * integration as determined by the set_integral_limits() function.
     */
    double integral_track_score(const TrackData& data,
                                const moab::CartVect& X);

    /**
     * \brief Determines integration limits for the integral-track estimator
     * \param data the physics data for the track segment being tallied
     * \param X the (x, y, z) coordinates of the calculation point
     * \param limits stores integration limits in form of pair<lower, upper>
     * \return true if valid limits exist, false otherwise
     *
     * The set_integral_limits() function determines the integration limits
     * for the integral of the 3D path-length dependent kernel function K(X, s)
     * with respect to path-length s for the given calculation point X.
     *
     * NOTE: if this function returns false, that means there are no valid
     * integration limits within the range from 0 to the total track length.
     * This essentially means that the resulting integrand over this range
     * would have been zero and the calculation point can be ignored.
     */
    bool set_integral_limits(const TrackData& data,
                             const moab::CartVect& X,
                             std::pair<double, double>& limits);

    /**
     * \brief Computes tally score based on the sub-track estimator
     * \param points the set of sub-track points needed for computing score
     * \param X the (x, y, z) coordinates of the calculation point
     * \return the tally score for the calculation point
     *
     * The subtrack_score() function computes the average 3D kernel
     * contribution for the number of sub-tracks requested for the given
     * calculation point X, using the randomly chosen points along the
     * track that were computed by the choose_points function.
     */
    double subtrack_score(const std::vector<moab::CartVect>& points,
                          const moab::CartVect& X);

    /**
     * \brief Chooses p random points along a track segment
     * \param data the physics data for the track segment being tallied
     * \param p the number of random points requested
     * \return the vector of p random points
     *
     * The choose_points() function sub-divides the track segment into p
     * sub-tracks of equal length and randomly chooses the coordinates of
     * one point from each sub-track.
     */
    std::vector<moab::CartVect> choose_points(const TrackData& data, int p);

    /**
     * \brief Computes tally score based on the collision estimator
     * \param data the physics data for the collision being tallied
     * \param X the (x, y, z) coordinates of the calculation point
     * \return the tally score for the calculation point
     *
     * The collision_score() function computes the 3D kernel contribution
     * for the given calculation point X.
     */
    double collision_score(const CollisionData& data, const moab::CartVect& X);
};

#endif
