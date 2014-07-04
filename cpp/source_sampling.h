/// \file source_sampling.h
/// \author Elliott Biondo (biondo\@wisc.edu)
///
/// \brief Mesh-based Monte Carlo source sampling.
/// 
/// The Sampler class is used for Monte Carlo source sampling from mesh-based
/// sources.  The source density distribution and optional biased source density
/// distribution are defined on a MOAB mesh. Upon instantiation, a Sampler  
/// object reads this mesh and creates an alias table for randomly sampling
/// particle birth parameters. The particle_birth member function is supplied 
/// with 6 pseudo-random numbers and returns the position, energy, and weight 
/// of a particle upon birth. 
/// There are three sampling modes: analog, uniform, and user-speficied
/// In analog sampling, no source biasing is used and birth weights
/// are all 1. In uniform sampling, all phase space is sampled evenly and weights
/// are adjusted accordingly. In user-speficied mode, a supplied biased source
/// density distribution is used for sampling and particle weights are adjusted
/// accordingly. The biased source density distribution must have the same number 
/// of energy groups as the unbiased distribution. Alternatively, it may have 
/// exactly 1 energy group, in which case this value is used for all energies
/// within a particular mesh volume element.
 
#ifndef PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY
#define PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY

#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept> 
#include <sstream>
#include <string>

#include "moab/Range.hpp"
#include "MBCore.hpp"
#include "measure.h"
#include "MBCartVect.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<double> vect_d;
typedef std::string str;

/// Stores 4 connected points in a mesh volume element
struct edge_points {
  MBCartVect o_point;
  MBCartVect x_vec;
  MBCartVect y_vec;
  MBCartVect z_vec;
};

/// A data structure for O(1) source sampling
class AliasTable
{
public:
  /// Constructor
  /// \param p A normalized probability distribution function
  AliasTable(vect_d p);
  /// Samples the alias table
  /// \param rand1 A random number in range [0, 1].
  /// \param rand2 A random number in range [0, 1].
  int sample_pdf(double ran1, double ran2);
  ~AliasTable(){};
private:
  std::vector<double> prob;
  std::vector<int> alias;
  int n;
};

/// Problem modes
enum Mode {USER, ANALOG, UNIFORM};

/// Mesh based Monte Carlo source sampling.
class Sampler
{
public:
  /// Constuctor for analog and uniform sampling
  /// \param filename The path to the MOAB mesh (.h5m) file
  /// \param src_tag_name The name of the tag that describes the unbiased 
  ///                     source density distribution.
  /// \param e_bounds The energy boundaries, note there are N + 1 energy
  ///                 bounds for N energy groups
  /// \param uniform If false, analog sampling is used. If true, uniform
  ///                sampling is used.
  Sampler(str filename, str src_tag_name, vect_d e_bounds, bool uniform);
  /// Constuctor for analog and uniform sampling
  /// \param filename The path to the MOAB mesh (.h5m) file
  /// \param src_tag_name The name of the tag with the unbiased source density
  ///                     distribution.
  /// \param e_bounds The energy boundaries, note there are N + 1 energy
  ///                 bounds for N energy groups
  /// \param bias_tag_name The name of the tag describing the biased source
  ///                       density distribution. The tag must have the same
  ///                       number of energy groups as <src_tag_name> or 1.
  ///                       If 1 (i.e. spacial biasing only), all energy groups
  ///                       within a mesh volume element are sampled equally.
  Sampler(str filename, str src_tag_name, vect_d e_bounds, str bias_tag_name);
  /// Samples particle birth parameters
  /// \param rands Six pseudo-random numbers in range [0, 1].
  /// \return A vector containing the x position, y, position, z, point, energy
  ///         and weight of a particle (in that order).
  vect_d particle_birth(vect_d rands);
  ~Sampler() {
    delete _mesh;
    delete _at;
  };

// member variables
private:
  // problem parameters
  str _filename; ///< MOAB mesh file path
  str _src_tag_name; ///< Unbiased source density distribution
  str _bias_tag_name; ///< Biased source density distribution
  vect_d _e_bounds;  ///< Energy boundaries
  int _num_e_groups; ///< Number of groups in tag \a _src_tag_name
  int _num_bias_groups; ///< Number of groups tag \a _bias_tag_name
  Mode _mode; ///< Problem mode: analog, uniform, user
  // mesh
  MBInterface *_mesh; ///< MOAB mesh
  int _num_ves; ///< Number of mesh volume elements on \a mesh.
  MBEntityType _ve_type; ///< Type of mesh volume: MBTET or MBHEX
  int _verts_per_ve; ///< Number of verticles per mesh volume element
  // sampling
  std::vector<edge_points> _all_edge_points; ///< Four connected points on a VE.
  std::vector<double> _biased_weights; ///< Birth weights for biased sampling.
  AliasTable* _at; ///< Alias table used for sampling.

// member functions
private:
  // instantiation
  void setup();
  void get_mesh_geom_data(MBRange ves, vect_d &volumes);
  void get_mesh_tag_data(MBRange ves, const vect_d volumes);
  // select birth parameters
  MBCartVect get_xyz(int ve_idx, vect_d rands);
  double get_e(int e_idx, double rand);
  double get_w(int pdf_idx);
  // helper functions
  void normalize_pdf(vect_d & pdf);
  int get_num_groups(MBTag tag);
  vect_d get_bias_pdf(MBRange ves, vect_d volumes);
};

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY
