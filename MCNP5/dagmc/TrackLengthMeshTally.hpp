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


class TrackLengthMeshTally : public MeshTally
{ 

public:
  TrackLengthMeshTally( const TallyInput& input );

  /**
   * \brief Computes mesh tally scores for the given tally event
   * \param event the parameters needed to compute the mesh tally scores
   * \param ebin index representing energy bin
   */
  virtual void compute_score(const TallyEvent& event);
  
  virtual void end_history ();
 
  virtual void write_data( double num_histories);
    
  ~TrackLengthMeshTally();


protected:

  ErrorCode compute_barycentric_data (const Range& all_tets);
  void build_trees (Range& all_tets);
                       
  bool point_in_tet( const CartVect& point, const EntityHandle* tet );

 
  void get_skin_triangle_adjacencies( EntityHandle triangle, 
                                      EntityHandle& tetrahedron, EntityHandle vertices[3] );


  EntityHandle find_next_tet_by_ray_fire(const CartVect& start, const CartVect& vec, double length, 
                                         EntityHandle first_tri[3], double& first_t, 
                                         EntityHandle last_crossing = 0);

  EntityHandle get_starting_tet(const CartVect& start, const CartVect& vec, double length, 
                                EntityHandle first_tri[3], double& first_t);
    

  EntityHandle get_starting_tet_conformal(const CartVect& start, EntityHandle first_tri[3]);


  Interface* mb;

  AdaptiveKDTree* kdtree;
  EntityHandle kdtree_root;

  OrientedBoxTreeTool* obb_tool;
  EntityHandle obbtree_root;

  EntityHandle last_visited_tet;

  std::vector< Matrix3 > tet_baryc_data;

  bool convex; // true if user asserts this mesh tally has convex geometry

  // if not empty, user has asserted mesh tally geometry
  // conforms to the cells identified in this set
  std::set<int> conformality; 
  bool conformal_surface_source;

  int last_cell;
    

private:
  void parse_tally_options();
  void set_tally_meshset();
  TrackLengthMeshTally& operator=( const TrackLengthMeshTally& mt ); // unimplemented
  TrackLengthMeshTally( const TrackLengthMeshTally& mt ); // unimplemented

  int num_negative_tracks;
};

} // end namespace moab


#endif /* DAGMC_TRACK_LENGTH_MESH_TALLY_HPP */
