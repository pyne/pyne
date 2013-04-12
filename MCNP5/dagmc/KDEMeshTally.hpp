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
     * \brief Computes mesh tally scores for the given tally event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin index representing energy bin
     */
    void compute_score(const TallyEvent& event, int ebin);

    // TODO combine tally_collision and tally_track into one function
    /**
     * \brief Computes mesh tally scores for the given collision event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin the energy bin to tally this collision into (calculated by fortran)
     * TODO remove ebin as a parameter
     */
    void tally_collision(const TallyEvent& event, int ebin);

    /**
     * \brief Computes mesh tally scores for the given track-based event
     * \param event the parameters needed to compute the mesh tally scores
     * \param ebin the energy bin to tally this track into (calculated by fortran)
     * TODO remove ebin as a parameter
     */
    void tally_track(const TallyEvent& event, int ebin);

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

};

#endif
