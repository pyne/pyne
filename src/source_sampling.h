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
/// are all 1. In uniform sampling, the position of the particle (but not the 
/// energy) is sampled uniformly and weights are adjusted accordingly. In 
/// user-speficied mode, a supplied biased source density distribution is used 
/// for sampling and particle weights are adjusted accordingly. The biased 
/// source density distribution must have the same number of energy groups as 
/// the unbiased distribution. Alternatively, it may have exactly 1 energy
/// group, in which case only spatial biasing is done, and energies are sampled
/// in analog.
 
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
#include "moab/Core.hpp"
#include "measure.h"
#include "moab/CartVect.hpp"

#ifdef __cplusplus
extern "C" {
#endif

namespace pyne {

  /// MCNP interface for source sampling setup
  /// \param mode The sampling mode: 1 = analog, 2 = uniform, 3 = user-specified
  void sampling_setup_(int* mode);
  /// MCNP interface to sample particle birth parameters after sampling setup
  /// \param rands Six pseudo-random numbers supplied from the Fortran side.
  /// \param x The sampled x position returned by this function
  /// \param y The sampled y position returned by this function
  /// \param z The sampled z position returned by this function
  /// \param e The sampled energy returned by this function
  /// \param w The sampled y statistical weight returned by this function
  void particle_birth_(double* rands,
                            double* x,
                            double* y,
                            double* z,
                            double* e,
                            double* w);
  /// Helper function for MCNP interface that reads energy boudaries from a file
  /// \param e_bounds_file A file containing the energy group boundaries.
  std::vector<double> read_e_bounds(std::string e_bounds_file);

  /// Stores 4 connected points in a mesh volume element
  struct edge_points {
    moab::CartVect o_point;
    moab::CartVect x_vec;
    moab::CartVect y_vec;
    moab::CartVect z_vec;
  };
  
  /// A data structure for O(1) source sampling
  class AliasTable {
  public:
    /// Constructor
    /// \param p A normalized probability distribution function
    AliasTable(std::vector<double> p);
    /// Samples the alias table
    /// \param rand1 A random number in range [0, 1].
    /// \param rand2 A random number in range [0, 1].
    int sample_pdf(double rand1, double rand2);
    ~AliasTable(){};
    int n; /// Number of bins in the PDF.
    std::vector<double> prob; /// Probabilities.
    std::vector<int> alias; /// Alias probabilities.
  };
  
  /// Problem modes
  enum Mode {USER, ANALOG, UNIFORM};
  
  /// Mesh based Monte Carlo source sampling.
  class Sampler {
  public:
    /// Constuctor for analog and uniform sampling
    /// \param filename The path to the MOAB mesh (.h5m) file
    /// \param src_tag_name The name of the tag that describes the unbiased 
    ///                     source density distribution.
    /// \param e_bounds The energy boundaries, note there are N + 1 energy
    ///                 bounds for N energy groups
    /// \param uniform If false, analog sampling is used. If true, uniform
    ///                sampling is used.
    Sampler(std::string filename, 
            std::string src_tag_name, 
            std::vector<double> e_bounds, 
            bool uniform);
    /// Constuctor for analog and uniform sampling
    /// \param filename The path to the MOAB mesh (.h5m) file
    /// \param src_tag_name The name of the tag with the unbiased source density
    ///                     distribution.
    /// \param e_bounds The energy boundaries, note there are N + 1 energy
    ///                 bounds for N energy groups
    /// \param bias_tag_name The name of the tag describing the biased source
    ///                       density distribution. The tag must have the same
    ///                       number of energy groups as <src_tag_name> or 1.
    ///                       If 1 (i.e. spatial biasing only), all energy groups
    ///                       within a mesh volume element are sampled equally.
    Sampler(std::string filename, 
            std::string src_tag_name, 
            std::vector<double> e_bounds, 
            std::string bias_tag_name);
    /// Samples particle birth parameters
    /// \param rands Six pseudo-random numbers in range [0, 1].
    /// \return A vector containing the x position, y, position, z, point, energy
    ///         and weight of a particle (in that order).
    std::vector<double> particle_birth(std::vector<double> rands);
    ~Sampler() {
      delete mesh;
      delete at;
    };
  
  // member variables
  private:
    // problem parameters
    std::string filename; ///< MOAB mesh file path
    std::string src_tag_name; ///< Unbiased source density distribution
    std::string bias_tag_name; ///< Biased source density distribution
    std::vector<double> e_bounds;  ///< Energy boundaries
    int num_e_groups; ///< Number of groups in tag \a _src_tag_name
    int num_bias_groups; ///< Number of groups tag \a _bias_tag_name
    Mode mode; ///< Problem mode: analog, uniform, user
    // mesh
    moab::Interface* mesh; ///< MOAB mesh
    int num_ves; ///< Number of mesh volume elements on \a mesh.
    moab::EntityType ve_type; ///< Type of mesh volume: moab::TET or moab::HEX
    int verts_per_ve; ///< Number of verticles per mesh volume element
    // sampling
    std::vector<edge_points> all_edge_points; ///< Four connected points on a VE.
    std::vector<double> biased_weights; ///< Birth weights for biased sampling.
    AliasTable* at; ///< Alias table used for sampling.
  
  // member functions
  private:
    // instantiation
    void setup();
    void mesh_geom_data(moab::Range ves, std::vector<double> &volumes);
    void mesh_tag_data(moab::Range ves, const std::vector<double> volumes);
    // select birth parameters
    moab::CartVect sample_xyz(int ve_idx, std::vector<double> rands);
    double sample_e(int e_idx, double rand);
    double sample_w(int pdf_idx);
    // helper functions
    void normalize_pdf(std::vector<double> & pdf);
    int num_groups(moab::Tag tag);
    std::vector<double> read_bias_pdf(moab::Range ves, std::vector<double> volumes, 
                                      std::vector<double> pdf);
  };
} //end namespace pyne

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY
