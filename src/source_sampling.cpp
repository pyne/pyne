#ifndef PYNE_IS_AMALGAMATED
#include "source_sampling.h"
#endif

// Global sampler instance
static pyne::Sampler* sampler = NULL;

// Fortran API
void pyne::sampling_setup_(int* mode) {
  if (sampler == NULL) {
    std::string filename ("source.h5m");
    std::string src_tag_name ("source_density");
    std::string e_bounds_file ("e_bounds");
    std::vector<double> e_bounds = read_e_bounds(e_bounds_file);
    if (*mode == 0) {
      sampler = new pyne::Sampler(filename, src_tag_name, e_bounds, false);
    } else if (*mode == 1) {
      sampler = new pyne::Sampler(filename, src_tag_name, e_bounds, true);
    } else if (*mode == 2) {
      std::string bias_tag_name ("biased_source_density");
      sampler = new pyne::Sampler(filename, src_tag_name, e_bounds, bias_tag_name);
    }
  }
}

void pyne::particle_birth_(double* rands,
                           double* x,
                           double* y,
                           double* z,
                           double* e,
                           double* w) {
    std::vector<double> rands2(rands, rands + 6);
    std::vector<double> samp = sampler->particle_birth(rands2);
    *x = samp[0];
    *y = samp[1];
    *z = samp[2];
    *e = samp[3];
    *w = samp[4];
}

std::vector<double> pyne::read_e_bounds(std::string e_bounds_file){
  std::vector<double> e_bounds;
  std::ifstream inputFile(e_bounds_file.c_str());
  double value;
  if (inputFile) {
    while (inputFile >> value)
      e_bounds.push_back(value);
  }
  return e_bounds;
}


// C++ API
pyne::Sampler::Sampler(std::string filename, 
                 std::string src_tag_name, 
                 std::vector<double> e_bounds, 
                 bool uniform)
  : filename(filename), src_tag_name(src_tag_name), e_bounds(e_bounds) {
  mode = (uniform) ? UNIFORM : ANALOG;
  setup();
}

pyne::Sampler::Sampler(std::string filename,  
                 std::string src_tag_name, 
                 std::vector<double> e_bounds, 
                 std::string bias_tag_name)
  : filename(filename), 
    src_tag_name(src_tag_name), 
    e_bounds(e_bounds), 
    bias_tag_name(bias_tag_name) {
  mode = USER;
  setup();
}

std::vector<double> pyne::Sampler::particle_birth(std::vector<double> rands) {
  // select mesh volume and energy group
  //
  int pdf_idx =at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx/num_e_groups;
  int e_idx = pdf_idx % num_e_groups;

  // Sample uniformly within the selected mesh volume element and energy
  // group.
  std::vector<double> samp;
  std::vector<double> xyz_rands;
  xyz_rands.push_back(rands[2]);
  xyz_rands.push_back(rands[3]);
  xyz_rands.push_back(rands[4]);
  moab::CartVect pos = sample_xyz(ve_idx, xyz_rands);
  samp.push_back(pos[0]); 
  samp.push_back(pos[1]); 
  samp.push_back(pos[2]); 
  samp.push_back(sample_e(e_idx, rands[5]));
  samp.push_back(sample_w(pdf_idx));
  return samp;
}

void pyne::Sampler::setup() {
  moab::ErrorCode rval;
  moab::EntityHandle loaded_file_set;
  // Create MOAB instance
  mesh = new moab::Core();
  rval = mesh->create_meshset(moab::MESHSET_SET, loaded_file_set);
  rval = mesh->load_file(filename.c_str(), &loaded_file_set);
  if (rval != moab::MB_SUCCESS)
    throw std::invalid_argument("Could not load mesh file.");

  // Get mesh volume elemebts 
  moab::Range ves;
  rval = mesh->get_entities_by_dimension(loaded_file_set, 3, ves);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem entities of dimension 3");
  num_ves = ves.size();
  int num_hex, num_tet;
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBHEX, num_hex);
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBTET, num_tet);
  if (num_hex == num_ves) {
    ve_type = moab::MBHEX;
    verts_per_ve = 8;
  } else if (num_tet == num_ves) {
    ve_type = moab::MBTET;
    verts_per_ve = 4;
  }
  else throw std::invalid_argument("Mesh file must contain only tets or hexes.");

  // Process all the spatial and tag data and create an alias table.
  std::vector<double> volumes(num_ves);
  mesh_geom_data(ves, volumes);
  mesh_tag_data(ves, volumes);
}

void pyne::Sampler::mesh_geom_data(moab::Range ves, std::vector<double> &volumes) {
  // Get connectivity.
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> connect;
  rval = mesh->get_connectivity_by_type(ve_type, connect);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting mesh connectivity.");

  // Grab the coordinates that define 4 connected points within a mesh volume
  // element and setup a data structure to allow uniform sampling with each 
  // mesh volume element.
  double coords[verts_per_ve*3];
  int i;
  for (i=0; i<num_ves; ++i) {
    rval = mesh->get_coords(&connect[verts_per_ve*i], verts_per_ve, &coords[0]);
    if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem vertex coordinates.");
    volumes[i] = measure(ve_type, verts_per_ve, &coords[0]);
    if (ve_type == moab::MBHEX) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[9], coords[10], coords[11]);
      moab::CartVect z(coords[12], coords[13], coords[14]);
      edge_points ep = {o, x-o, y-o, z-o};
      all_edge_points.push_back(ep);
    } else if (ve_type == moab::MBTET) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[6], coords[7], coords[8]);
      moab::CartVect z(coords[9], coords[10], coords[11]);
      edge_points ep = {o, x-o, y-o, z-o};
      all_edge_points.push_back(ep);
    }
  }
}

void pyne::Sampler::mesh_tag_data(moab::Range ves, 
                                  const std::vector<double> volumes) {
  moab::ErrorCode rval;
  moab::Tag src_tag;
  rval = mesh->tag_get_handle(src_tag_name.c_str(),
                              moab::MB_TAG_VARLEN, 
                              moab::MB_TYPE_DOUBLE, 
                              src_tag);
  // THIS rval FAILS because we do not know number of energy groups a priori.
  // That's okay. That's what the next line is all about:
  num_e_groups = num_groups(src_tag);
  std::vector<double> pdf(num_ves*num_e_groups); 
  rval = mesh->tag_get_data(src_tag, ves, &pdf[0]);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting source tag data.");

  // Multiply the source densities by the VE volumes
  int i, j;
  for (i=0; i<num_ves; ++i) {
    for (j=0; j<num_e_groups; ++j) {
       pdf[i*num_e_groups + j] *= volumes[i];
    }
  }
  normalize_pdf(pdf);

  // Setup alias table based off PDF or biased PDF
  if (mode == ANALOG) {
    at = new AliasTable(pdf);
  } else {
    std::vector<double> bias_pdf = read_bias_pdf(ves, volumes, pdf);
    normalize_pdf(bias_pdf);
    //  Create alias table based off biased pdf and calculate birth weights.
    biased_weights.resize(num_ves*num_e_groups);
    for (i=0; i<num_ves*num_e_groups; ++i) {
      biased_weights[i] = pdf[i]/bias_pdf[i];
    }
    at = new AliasTable(bias_pdf);
  }
}

std::vector<double> pyne::Sampler::read_bias_pdf(moab::Range ves, 
                                                 std::vector<double> volumes,
                                                 std::vector<double> pdf) {
    std::vector<double> bias_pdf(num_ves*num_e_groups);
    int i, j;
    moab::ErrorCode rval;
    if (mode == UNIFORM) {
      // Uniform sampling: uniform in space, analog in energy. Biased PDF is
      // found by normalizing the total photon emission density to 1 in each
      // mesh volume element and multiplying by the volume of the element.
      double q_in_group;
      for (i=0; i<num_ves; ++i) {
        q_in_group = 0;
        for (j=0; j<num_e_groups; ++j){
          q_in_group += pdf[i*num_e_groups + j];
        }
        if (q_in_group > 0) {
          for (j=0; j<num_e_groups; ++j) {
            bias_pdf[i*num_e_groups + j] =
                volumes[i]*pdf[i*num_e_groups + j]/q_in_group;
          }
        } else {
          for (j=0; j<num_e_groups; ++j) {
            bias_pdf[i*num_e_groups + j] = 0.0;
          }
        }
      }
    } else if (mode == USER) {
      // Get the biased PDF from the mesh
      moab::Tag bias_tag;
      rval = mesh->tag_get_handle(bias_tag_name.c_str(), 
                                  moab::MB_TAG_VARLEN, 
                                  moab::MB_TYPE_DOUBLE, 
                                  bias_tag);
      num_bias_groups = num_groups(bias_tag);

      if (num_bias_groups == num_e_groups) {
        rval = mesh->tag_get_data(bias_tag, ves, &bias_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        for (i=0; i<num_ves; ++i) {
          for (j=0; j<num_e_groups; ++j)
             bias_pdf[i*num_e_groups + j] *=  volumes[i];
        }
      } else if (num_bias_groups == 1) {
        // Spatial biasing only: the supplied bias PDF values are applied
        // to all energy groups within a mesh volume element, which are
        // sampled in analog.
        std::vector<double> spatial_pdf(num_ves); 
        rval = mesh->tag_get_data(bias_tag, ves, &spatial_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        double q_in_group;
        for (i=0; i<num_ves; ++i) {
          q_in_group = 0;
          for (j=0; j<num_e_groups; ++j){
            q_in_group += pdf[i*num_e_groups + j];
          }
          if (q_in_group > 0){
            for (j=0; j<num_e_groups; ++j){
              bias_pdf[i*num_e_groups + j] = 
                spatial_pdf[i]*volumes[i]*pdf[i*num_e_groups + j]/q_in_group;
            }
          } else {
            for (j=0; j<num_e_groups; ++j)
              bias_pdf[i*num_e_groups + j] =  0;
          }
        }
      } else {
        throw std::length_error("Length of bias tag must equal length of the"
                                "  source tag, or 1.");
      }
    }
return bias_pdf;
}

moab::CartVect pyne::Sampler::sample_xyz(int ve_idx, std::vector<double> rands) {
  double s = rands[0];
  double t = rands[1];
  double u = rands[2];

  // Transform s, t, u to uniformly sample a tetrahedron. See:
  // C. Rocchini and P. Cignoni, “Generating Random Points in a Tetrahedron,” 
  //  Journal of Graphics Tools, 5, 200–202 (2001).
  if (ve_type == moab::MBTET) {
    if (s + t > 1) {
      s = 1.0 - s;
      t = 1.0 - t;
    }
    if (s + t + u > 1) {
      if (t + u > 1) {
        double old_t = t;
        t = 1.0 - u;
        u = 1.0 - s - old_t;
      }else if (t + u <= 1) {
        double old_s = s;
        s = 1.0 - t - u;
        u = old_s + t + u - 1;
      }
    }
  }

 return s*all_edge_points[ve_idx].x_vec + \
        t*all_edge_points[ve_idx].y_vec + \
        u*all_edge_points[ve_idx].z_vec + \
          all_edge_points[ve_idx].o_point;
}

double pyne::Sampler::sample_e(int e_idx, double rand) {
   double e_min = e_bounds[e_idx];
   double e_max = e_bounds[e_idx + 1];
   return rand * (e_max - e_min) + e_min;
}

double pyne::Sampler::sample_w(int pdf_idx) {
  return (mode == ANALOG) ? 1.0 : biased_weights[pdf_idx];
}

void pyne::Sampler::normalize_pdf(std::vector<double> & pdf) {
  double sum = 0;
  int i;
  for (i=0; i<num_ves*num_e_groups; ++i)
    sum += pdf[i];
  for (i=0; i<num_ves*num_e_groups; ++i)
    pdf[i] /= sum;
}

int pyne::Sampler::num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mesh->tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}


// Random-number sampling using the Walker-Vose alias method,
// Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
// M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
// A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)

pyne::AliasTable::AliasTable(std::vector<double> p) {
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  std::vector<double> small(n);
  std::vector<double> large(n);
  int i, a, g;

  for (i=0; i<n; ++i) 
    p[i] *= n;

  // Set separate index lists for small and large probabilities:
  int n_s = 0;
  int n_l = 0;
  for (i=n-1; i>=0; --i) {
    // at variance from Schwarz, we revert the index order
    if (p[i] < 1)
      small[n_s++] = i;
    else
      large[n_l++] = i;
  }

  // Work through index lists
  while(n_s && n_l) {
    a = small[--n_s]; // Schwarz's l
    g = large[--n_l]; // Schwarz's g
    prob[a] = p[a];
    alias[a] = g;
    p[g] = p[g] + p[a] - 1;
    if (p[g] < 1)
      small[n_s++] = g;
    else
      large[n_l++] = g;
  }

  while(n_l)
    prob[large[--n_l]] = 1;

  while(n_s)
    // can only happen through numeric instability
    prob[small[--n_s] ] = 1;
}

int pyne::AliasTable::sample_pdf(double rand1, double rand2) {
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}
