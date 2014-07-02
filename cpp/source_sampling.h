#ifndef PYNE_AIQ4M73S
#define PYNE_AIQ4M73S
/// \file source_sampling.h
/// \author Elliott Biondo (biondo\@wisc.edu)
///
/// \brief The tally class and helper functions
/// 
/// The tally class is in essesence a structure containing attributes
/// related to tallies

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
#include "measure.hpp"
#include "MBCartVect.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void get_e_bounds_data(std::string e_bounds_file);
typedef std::vector<double> vect_d;
typedef std::string str;


struct Sample {
    vect_d xyz; /// \param type
    double e;
    double w;
};

struct sample_vects{
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
  void setup();
  ~Sampler(){
     delete mesh;
     delete at;
     };
  Sample particle_birth(vect_d rands);


private:
  vect_d e_bounds;
  Mode mode;
  MBInterface *mesh;
  str filename;
  str src_tag_name;
  str bias_tag_name;
  int num_e_groups;
  MBEntityType ve_type;
  sample_vects vects;
  int verts_per_ve;
  std::vector<sample_vects> cart_sampler;
  AliasTable* at;
  std::vector<double> biased_weights;

private:
  // functions and classes
  void get_mesh_geom_data(MBRange ves, vect_d &volumes);
  void get_mesh_tag_data(MBRange ves, vect_d volumes);
  MBCartVect get_xyz(int ve_idx, vect_d rands);
  double get_e(int e_idx, double rand);
  double get_w(int pdf_idx);
};

#ifdef __cplusplus
} //  extern "C"
#endif

#endif
