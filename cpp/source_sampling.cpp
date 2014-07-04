#ifndef PYNE_IS_AMALGAMATED
#include "source_sampling.h"
#endif

Sampler::Sampler(str filename, str src_tag_name, vect_d e_bounds, bool uniform)
  : _filename(filename), _src_tag_name(src_tag_name), _e_bounds(e_bounds) {
  _mode = (uniform) ? UNIFORM : ANALOG;
  setup();
}

Sampler::Sampler(str filename,  
                 str src_tag_name, 
                 vect_d e_bounds, 
                 str bias_tag_name)
  : _filename(filename), 
    _src_tag_name(src_tag_name), 
    _e_bounds(e_bounds), 
    _bias_tag_name(bias_tag_name) {
  _mode = USER;
  setup();
}

vect_d Sampler::particle_birth(vect_d rands) {
  // select mesh volume and energy group
  int pdf_idx = _at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx/_num_e_groups;
  int e_idx = pdf_idx % _num_e_groups;

  // Sample uniformly within the selected mesh volume elemenet and energy
  // group.
  vect_d samp;
  vect_d xyz_rands;
  xyz_rands.push_back(rands[2]);
  xyz_rands.push_back(rands[3]);
  xyz_rands.push_back(rands[4]);
  MBCartVect pos = get_xyz(ve_idx, xyz_rands);
  samp.push_back(pos[0]); 
  samp.push_back(pos[1]); 
  samp.push_back(pos[2]); 
  samp.push_back(get_e(e_idx, rands[5]));
  samp.push_back(get_w(pdf_idx));
  return samp;
}

void Sampler::setup() {
  MBErrorCode rval;
  MBEntityHandle loaded_file_set;
  // Create MOAB instance
  _mesh = new MBCore();
  rval = _mesh->create_meshset(MESHSET_SET, loaded_file_set);
  rval = _mesh->load_file(_filename.c_str(), &loaded_file_set);
  if (rval != moab::MB_SUCCESS)
    throw std::invalid_argument("Could not load mesh file.");

  // Get mesh volume elemebts 
  MBRange ves;
  rval = _mesh->get_entities_by_dimension(loaded_file_set, 3, ves);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem entities of dimension 3");
  _num_ves = ves.size();
  int num_hex, num_tet;
  rval = _mesh->get_number_entities_by_type(loaded_file_set, MBHEX, num_hex);
  rval = _mesh->get_number_entities_by_type(loaded_file_set, MBTET, num_tet);
  if (num_hex == _num_ves) {
    _ve_type = MBHEX;
    _verts_per_ve = 8;
  } else if (num_tet == _num_ves) {
    _ve_type = MBTET;
    _verts_per_ve = 4;
  }
  else throw std::invalid_argument("Mesh file must contain only tets or hexes.");

  // Process all the spacial and tag data and create an alias table.
  vect_d volumes(_num_ves);
  get_mesh_geom_data(ves, volumes);
  get_mesh_tag_data(ves, volumes);
}

void Sampler::get_mesh_geom_data(MBRange ves, vect_d &volumes) {
  // Get connectivity.
  MBErrorCode rval;
  std::vector<MBEntityHandle> connect;
  rval = _mesh->get_connectivity_by_type(_ve_type, connect);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting mesh connectivity.");

  // Grab the coordinates that define 4 connected points within a mesh volume
  // element and setup a data structure to allow uniform sampling with each 
  // mesh volume element.
  double coords[_verts_per_ve*3];
  int i;
  for (i=0; i<_num_ves; ++i) {
    rval=_mesh->get_coords(&connect[_verts_per_ve*i], _verts_per_ve, &coords[0]);
    if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem vertex coordinates.");
    volumes[i] = measure(_ve_type, _verts_per_ve, &coords[0]);
    if (_ve_type == MBHEX) {
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[9], coords[10], coords[11]);
      MBCartVect z(coords[12], coords[13], coords[14]);
      edge_points ep = {o, x-o, y-o, z-o};
      _all_edge_points.push_back(ep);
    } else if (_ve_type == MBTET) {
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[6], coords[7], coords[8]);
      MBCartVect z(coords[9], coords[10], coords[11]);
      edge_points ep = {o, x-o, y-o, z-o};
      _all_edge_points.push_back(ep);
    }
  }
}

void Sampler::get_mesh_tag_data(MBRange ves, const vect_d volumes) {
  MBErrorCode rval;
  MBTag src_tag;
  rval = _mesh->tag_get_handle(_src_tag_name.c_str(),
                              moab::MB_TAG_VARLEN, 
                              MB_TYPE_DOUBLE, 
                              src_tag);
  // THIS rval FAILS because we do not know number of energy groups a priori.
  // That's okay. That's what the next line is all about:
  _num_e_groups = get_num_groups(src_tag);
  vect_d pdf(_num_ves*_num_e_groups); 
  rval = _mesh->tag_get_data(src_tag, ves, &pdf[0]);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting source tag data.");

  // Multiply the source densities by the VE volumes
  int i, j;
  for (i=0; i<_num_ves; ++i) {
    for (j=0; j<_num_e_groups; ++j) {
       pdf[i*_num_e_groups + j] *=  volumes[i];
    }
  }
  normalize_pdf(pdf);

  // Setup alias table based off PDF or biased PDF
  if (_mode == ANALOG) {
    _at = new AliasTable(pdf);
  } else {
    vect_d bias_pdf = get_bias_pdf(ves, volumes);
    normalize_pdf(bias_pdf);
    //  Create alias table based off biased pdf and calculate birth weights.
    _biased_weights.resize(_num_ves*_num_e_groups);
      for (i=0; i<_num_ves*_num_e_groups; ++i) {
        _biased_weights[i] = pdf[i]/bias_pdf[i];
      }
    _at = new AliasTable(bias_pdf);
  }
}

vect_d Sampler::get_bias_pdf(MBRange ves, vect_d volumes) {
    vect_d bias_pdf(_num_ves*_num_e_groups);
    int i, j;
    MBErrorCode rval;
    if (_mode == UNIFORM) {
      // In unform sampling, the biased PDF is just the volume of the mesh
      // volume element
      for (i=0; i<_num_ves; ++i) {
        for (j=0; j<_num_e_groups; ++j)
           bias_pdf[i*_num_e_groups + j] =  volumes[i];
      }
    } else if (_mode == USER) {
      // Get the biased PDF from the mesh
      MBTag bias_tag;
      rval = _mesh->tag_get_handle(_bias_tag_name.c_str(), 
                                  moab::MB_TAG_VARLEN, 
                                  MB_TYPE_DOUBLE, 
                                  bias_tag);
      _num_bias_groups = get_num_groups(bias_tag);

      if (_num_bias_groups == _num_e_groups) {
        rval = _mesh->tag_get_data(bias_tag, ves, &bias_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        for (i=0; i<_num_ves; ++i) {
          for (j=0; j<_num_e_groups; ++j)
             bias_pdf[i*_num_e_groups + j] *=  volumes[i];
        }
      } else if (_num_bias_groups == 1) {
        // Spacial biasing only: the supplied bias PDF values are supplied
        // to all energy groups within a mesh volume element
        vect_d spacial_pdf(_num_ves); 
        rval = _mesh->tag_get_data(bias_tag, ves, &spacial_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        for (i=0; i<_num_ves; ++i) {
          for (j=0; j<_num_e_groups; ++j)
            bias_pdf[i*_num_e_groups + j] =  spacial_pdf[i]*volumes[i];
        }
      } else {
        throw std::length_error("Length of bias tag must equal length of the"
                                "  source tag, or 1.");
      }
    }
return bias_pdf;
}

MBCartVect Sampler::get_xyz(int ve_idx, vect_d rands) {
  double s = rands[0];
  double t = rands[1];
  double u = rands[2];

  // Transform s, t, u to uniformly sample a tetrahedron
  if (_ve_type == MBTET) {
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

 return s*_all_edge_points[ve_idx].x_vec + \
        t*_all_edge_points[ve_idx].y_vec + \
        u*_all_edge_points[ve_idx].z_vec + \
          _all_edge_points[ve_idx].o_point;
}

double Sampler::get_e(int e_idx, double rand) {
   double e_min = _e_bounds[e_idx];
   double e_max = _e_bounds[e_idx + 1];
   return rand * (e_max - e_min) + e_min;
}

double Sampler::get_w(int pdf_idx) {
  return (_mode == ANALOG) ? 1.0 : _biased_weights[pdf_idx];
}

void Sampler::normalize_pdf(vect_d & pdf) {
  double sum = 0;
  int i;
  for (i=0; i<_num_ves*_num_e_groups; ++i)
    sum += pdf[i];
  for (i=0; i<_num_ves*_num_e_groups; ++i)
    pdf[i] /= sum;
}

int Sampler::get_num_groups(MBTag tag) {
  MBErrorCode rval;
  int tag_size;
  rval = _mesh->tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}


// Random-number sampling using the Walker-Vose alias method,
// Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
// M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
// A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)

AliasTable::AliasTable(vect_d p) {
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  vect_d small(n);
  vect_d large(n);
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

int AliasTable::sample_pdf(double rand1, double rand2) {
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}
