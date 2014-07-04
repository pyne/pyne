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
/// of a particle upon birth.  There are three sampling modes: analog, uniform, 
/// and user. In analog sampling, no source biasing is used and birth weights
/// are all 1. In uniform sampling, all phase space is sampled evenly and weights
/// are adjusted accordingly. In user mode, a user-supplied biased source
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

void get_e_bounds_data(std::string e_bounds_file);
typedef std::vector<double> vect_d;
typedef std::string str;

struct edge_points {
  MBCartVect o_point;
  MBCartVect x_vec;
  MBCartVect y_vec;
  MBCartVect z_vec;
};

class AliasTable
{
private:
  std::vector<double> prob;
  std::vector<int> alias;
  int n;
public:
  int sample_pdf(double ran1, double ran2);
  AliasTable(vect_d p);
  ~AliasTable(){};
};

enum Mode {USER, ANALOG, UNIFORM};

class Sampler
{
public:
  // Sampler constructors
  // \param filename the filename of the mesh file
  Sampler(str filename, str src_tag_name, vect_d e_bounds, bool uniform);
  Sampler(str filename, str src_tag_name, vect_d e_bounds, str bias_tag_name);
  ~Sampler() {
    delete _mesh;
    delete _at;
  };
  vect_d particle_birth(vect_d rands);

private:
  // member variables
  vect_d _e_bounds;
  Mode _mode;
  MBInterface *_mesh;
  str _filename;
  str _src_tag_name;
  str _bias_tag_name;
  int _num_e_groups;
  int _num_bias_groups;
  int _num_ves;
  MBEntityType _ve_type;
  int _verts_per_ve;
  std::vector<edge_points> _all_edge_points;
  AliasTable* _at;
  std::vector<double> _biased_weights;

private:
  // functions and classes
  void setup();
  void get_mesh_geom_data(MBRange ves, vect_d &volumes);
  void get_mesh_tag_data(MBRange ves, const vect_d volumes);
  MBCartVect get_xyz(int ve_idx, vect_d rands);
  double get_e(int e_idx, double rand);
  double get_w(int pdf_idx);
  void normalize_pdf(vect_d & pdf);
  int get_num_groups(MBTag tag);
  vect_d get_bias_pdf(MBRange ves, vect_d volumes);
};

#ifdef __cplusplus
} //  extern "C"
#endif

#endif
