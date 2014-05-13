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
#define MBI moab_instance()
#define SI Sampling::instance()


#ifdef __cplusplus
extern "C" {
#endif
// Fortran interface
//void Fsampling_setup_(char* file_name, char* src_tag_name, char* e_bound_tag_name, bool analog);
//void Fsampling_setup_(char* file_name, char* src_tag_name, char* e_bounds_file_name, bool analog, char* bias_tag_name);
//void Fparticle_birth_(double* rands, double &x, double &y, double &z, double &e, double &w);

void mcnp_sampling_setup_(bool* analog);
void fsampling_setup_(char*, char*, char*, bool*);
void fsampling_setup2_(char*, char*, char*, bool*, char*);
void fparticle_birth_(double*, double*, double *, double *, double *, double *);

class Sampling
{
public:
  static Sampling *instance(MBInterface *mb_impl = NULL);
  MBInterface* moab_instance() {return mbImpl;}
  ~Sampling();
  void sampling_setup(char* file_name, char* src_tag_name, char* e_bound_tag_name, bool analog);
  void sampling_setup(char* file_name, char* src_tag_name, char* e_bounds_file_name, bool analog, char* bias_tag_name);
  void particle_birth(double* rands, double* x, double* y, double* z, double* e, double* w);


private:
  // functions and classes
  void get_mesh_geom_data(MBRange ves, std::vector<double> &volumes);
  void get_mesh_tag_data(MBRange ves, std::vector<double>volumes);
  void get_e_bounds_data(char* e_bounds_file);
  void get_xyz(int ve_idx, double* rands, double* x, double *y, double *z);
  void get_e(int e_idx, double rand, double* e);
  void get_w(int pdf_idx, double* w);
  class AliasTable
  {
  private:
    std::vector<double> prob;
    std::vector<int> alias;
    int n;
  public:
    int sample_pdf(double ran1, double ran2);
    AliasTable(std::vector<double> p);
  };
  //  member variables
public:
  bool analog;
  bool uniform;
  char* src_tag_name;
  char* bias_tag_name;
  int num_e_groups;
  MBEntityType ve_type;
  std::vector<double> e_bounds;
  int verts_per_ve;
  struct vector_points{
    MBCartVect o_point;
    MBCartVect x_vec;
    MBCartVect y_vec;
    MBCartVect z_vec;
  };
  std::vector<vector_points> cart_sampler;
  AliasTable* at;
  std::vector<double> biased_weights;

private:
  Sampling(MBInterface *mb_impl);
  static void create_instance(MBInterface *mb_impl = NULL);
  static Sampling *instance_;
  MBInterface *mbImpl;
};


inline Sampling *Sampling::instance(MBInterface *mb_impl)
{
  if (NULL == instance_) create_instance(mb_impl);
  return instance_;
}

#ifdef __cplusplus
} //  extern "C"
#endif

