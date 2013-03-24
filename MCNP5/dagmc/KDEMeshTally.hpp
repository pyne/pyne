// KDEMeshTally.hpp

#ifndef KDEMESHTALLY_H
#define KDEMESHTALLY_H

#include <iosfwd>
#include <set>
#include <string>

#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"

#include "KDEKernel.hpp"
#include "MeshTally.hpp"

// forward declarations
class AdaptiveKDTree;

/**
 * A struct that contains all the parameters needed to compute the tally
 * weighting factor of a KDEMeshTally score.
 *
 * @param fmesh_index the index of the KDE mesh tally
 * @param particle_wgt the current weight of the particle
 * @param energy the current energy of the particle
 * @param tracklength (optional) the track length the particle has traveled
 * @param total_xs (optional) the total macroscopic cross section
 */
struct KDEWeightParam {

  // parameters needed for all KDE tally types
  int* fmesh_index;
  double* particle_wgt;
  double* energy;

  // optional parameters, default to NULL
  double* tracklength;  // only needed for TRACKLENGTH/SUBTRACK tallies
  double* total_xs;     // only needed for COLLISION tallies

  KDEWeightParam( int* index, double* wgt, double* erg ):
    fmesh_index( index ), particle_wgt( wgt ), energy( erg ),
    tracklength( NULL ), total_xs( NULL ) {}

};

/** 
 * A class that represents a mesh tally based on a kernel density estimator
 * approach for use in a Monte Carlo particle transport simulation.
 */
class KDEMeshTally : public MeshTally {

  public:

     static moab::CartVect default_bandwidth;  // bandwidth choice for KDE tallies

    /**
     * An enumerative type that specifies the type of tally being used by this
     * KDEMeshTally object.
     *
     *  - COLLISION tallies are based on the KDE collision estimator
     *  - TRACKLENGTH tallies are based on the KDE integral track estimator
     *  - SUBTRACK tallies are based on the KDE subtrack estimator
     */
    enum TallyType { COLLISION = 0, TRACKLENGTH = 1, SUBTRACK = 2 };

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
     * NOTE: the parameter "numSubtracks" is required for SUBTRACK tallies,
     * but it is optional for TRACKLENGTH tallies for computing an optimal
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
     * Computes the contribution for the given collision location and adds the
     * results to the current mesh tally.  Note that the total cross section
     * parameter in the KDEWeightParam object must be set to a non-NULL value.
     *
     * @param param the parameters needed to compute the tally weighting factor
     * @param collision_loc a point where a collision has occurred
     * @param ebin the energy bin to tally this collision into (calculated by fortran)
     */
    void tally_collision( const KDEWeightParam & param,
                          const moab::CartVect & collision_loc, 
                          int ebin );

    /**
     * Computes the contribution for the given track segment and adds the
     * results to the current mesh tally.  Note that the track length parameter
     * in the KDEWeightParam object must be set to a non-NULL value.
     *
     * @param param the parameters needed to compute the tally weighting factor
     * @param start_point the starting location of the track (xo, yo, zo)
     * @param direction the direction the particle is traveling (uo, vo, wo)
     * @param ebin the energy bin to tally this track into (calculated by fortran)
     */
    void tally_track( const KDEWeightParam & param,
                      const moab::CartVect & start_point,
                      const moab::CartVect & direction, 
                      int ebin );

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
    unsigned int subtracks;        // number of subtracks used in tally       
    
    std::string output_filename;
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
    
    /**
     * Builds a KD-Tree from all the mesh node entities in the specified meshset,
     * or if the set is 0, from all mesh nodes in MOAB.
     */
    void build_tree( moab::EntityHandle meshset );

    /**
     * Adds the given collision coordinates to the running variance formula
     * used to compute the optimal bandwidth.
     */
    void update_bandwidth_variance( const moab::CartVect & collision_loc );

    /**
     * Computes and returns the optimal bandwidth to achieve the best results
     * in the kernel density approximation.
     */
    moab::CartVect get_optimal_bandwidth();

    /**
     * Determines the tally weighting factor for a collision or track
     * contribution.  This weight includes the particle weight and any tally
     * multipliers that were used for both tally types.
     * 
     * @param param the parameters needed to compute the tally weighting factor
     * @return the tally weighting factor for a collision or track
     */
    double get_score_weight( const KDEWeightParam & param );

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

    /**
     * Gets a unique set of all of the vertices that lie within the given box.
     *
     * @param tree a KD-Tree containing all vertices
     * @param tree_root the root of the KD-Tree
     * @param box_min the minimum corner of the box
     * @param box_max the maximum corner of the box
     * @param points  the unique set of vertices in the box
     * @return the MOAB ErrorCode value
     */
    moab::ErrorCode points_in_box( moab::AdaptiveKDTree & tree,
                                   moab::EntityHandle tree_root,
                                   const double box_min[3],
                                   const double box_max[3],
                                   std::vector<moab::EntityHandle> & points );

    /// unimplemented 
    KDEMeshTally( const KDEMeshTally & obj );
    /// unimplemented
    KDEMeshTally & operator=( const KDEMeshTally & obj );

};

#endif
