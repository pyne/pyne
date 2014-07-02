#ifndef PYNE_IS_AMALGAMATED
#include "source_sampling.h"
#endif


//    throw std::invalid_argument("bias_tag_name should not be specified for analog sampling");

void get_e_bounds_data(str e_bounds_file){
  std::ifstream inputFile(&e_bounds_file[0]);
  // test file open   
  vect_d e_bounds;
  if (inputFile) {        
    double value;
    while (inputFile >> value)
      e_bounds.push_back(value);
  }
}

Sampler::Sampler(str _filename, str _src_tag_name, vect_d _e_bounds, bool _uniform)
  : filename(_filename), src_tag_name(_src_tag_name), e_bounds(_e_bounds) {
  if(_uniform){
   mode = UNIFORM;
  } else {
   mode = ANALOG;
  }
  setup();
}

Sampler::Sampler(str _filename, str _src_tag_name, 
                  vect_d _e_bounds, str _bias_tag_name)
  : filename(_filename), src_tag_name(_src_tag_name), e_bounds(_e_bounds), bias_tag_name(_bias_tag_name) {
  mode = USER;
  setup();
}

void Sampler::setup(){
  MBErrorCode rval;
  MBEntityHandle loaded_file_set;
  mesh = new MBCore();
  rval = mesh->create_meshset(MESHSET_SET, loaded_file_set);
  rval = mesh->load_file(filename.c_str(), &loaded_file_set);
  MBRange ves;
  rval = mesh->get_entities_by_dimension(loaded_file_set, 3, ves);

  int num_hex, num_tet;
  rval = mesh->get_number_entities_by_type(loaded_file_set, MBHEX, num_hex);
  rval = mesh->get_number_entities_by_type(loaded_file_set, MBTET, num_tet);
  if(num_hex == ves.size()){
    ve_type = MBHEX;
    verts_per_ve = 8;
  } else if (num_tet == ves.size()){
    ve_type = MBTET;
    verts_per_ve = 4;
  }
  else throw std::invalid_argument("Mesh file must contain only tets or hexes");

  vect_d volumes(ves.size());
  get_mesh_geom_data(ves, volumes);
  get_mesh_tag_data(ves, volumes);
}

void Sampler::get_mesh_geom_data(MBRange ves, vect_d &volumes){
  MBErrorCode rval;
  std::vector<MBEntityHandle> connect;
  rval = mesh->get_connectivity_by_type(ve_type, connect);
  double coords[verts_per_ve*3];
  int i;
  for(i=0; i<ves.size(); ++i){
    rval=mesh->get_coords(&connect[verts_per_ve*i], verts_per_ve, &coords[0]);
    volumes[i] = measure(ve_type, verts_per_ve, &coords[0]);

    if(ve_type == MBHEX){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[9], coords[10], coords[11]);
      MBCartVect z(coords[12], coords[13], coords[14]);
      sample_vects vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
   }else if (ve_type == MBTET){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[6], coords[7], coords[8]);
      MBCartVect z(coords[9], coords[10], coords[11]);
      sample_vects vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
    }
  }
}

void Sampler::get_mesh_tag_data(MBRange ves, vect_d volumes){
  MBErrorCode rval;
  MBTag src_tag;
  rval = mesh->tag_get_handle(src_tag_name.c_str(),
                              1,
                              // moab::MB_TAG_VARLEN, 
                              MB_TYPE_DOUBLE, 
                              src_tag);
  // THIS ASSERT FAILS because we do not know number of energy groups a priori.
  //assert( rval == MB_SUCCESS );
  int src_tag_size;
  rval = mesh->tag_get_bytes(src_tag, *(&src_tag_size));
  num_e_groups = src_tag_size/sizeof(double);

  vect_d pdf(ves.size()*num_e_groups); 
  rval = mesh->tag_get_data(src_tag, ves, &pdf[0]);

  int i, j;
  for(i=0; i<ves.size(); ++i){
    for(j=0; j<num_e_groups; ++j){
       pdf[i*num_e_groups + j] *=  volumes[i];
    }
  }
  //normalize
  double sum = 0;
  for(i=0; i<ves.size()*num_e_groups; ++i){
    sum += pdf[i];
  }
  for(i=0; i<ves.size()*num_e_groups; ++i){
    pdf[i] /= sum;
  }

  if(mode != ANALOG){
    vect_d bias_pdf(ves.size()*num_e_groups); 
    if(mode == UNIFORM){
      for(i=0; i<ves.size(); ++i){
        for(j=0; j<num_e_groups; ++j)
           bias_pdf[i*num_e_groups + j] =  volumes[i];
      }
    }else if(mode == USER){
      MBTag bias_tag;
      rval = mesh->tag_get_handle(bias_tag_name.c_str(), 
                                  moab::MB_TAG_VARLEN, 
                                  MB_TYPE_DOUBLE, 
                                  bias_tag);
      int bias_tag_size;
      rval = mesh->tag_get_bytes(bias_tag, *(&bias_tag_size));
      int num_bias_groups = bias_tag_size/sizeof(double);

      if (num_bias_groups == num_e_groups){
        rval = mesh->tag_get_data(bias_tag, ves, &bias_pdf[0]);
        for(i=0; i<ves.size(); ++i){
          for(j=0; j<num_e_groups; ++j)
             bias_pdf[i*num_e_groups + j] *=  volumes[i];
        }
      }else if(num_bias_groups == 1){
        vect_d spacial_pdf(ves.size()); 
        rval = mesh->tag_get_data(bias_tag, ves, &spacial_pdf[0]);
        for(i=0; i<ves.size(); ++i){
          bias_pdf[i*num_e_groups + j] *=  spacial_pdf[i]*volumes[i];
        }
      }
      else
        std::cout<<num_bias_groups<<" "<<num_e_groups<<std::endl;
        throw std::length_error("Length of bias tag must equal length of the source tag, or 1.");

      sum = 0;
      for(i=0; i<ves.size()*num_e_groups; ++i){
        sum += bias_pdf[i];
      }
      for(i=0; i<ves.size()*num_bias_groups; ++i){
        bias_pdf[i] /= sum;
      }
    }
    //  Create alias table based off biased pdf and calculate birth weights.
    at = new AliasTable(bias_pdf);
    biased_weights.resize(ves.size()*num_e_groups);
      for(i=0; i<ves.size()*num_e_groups; ++i){
        biased_weights[i] = pdf[i]/bias_pdf[i];
      }
  }else if(mode == ANALOG){
    at = new AliasTable(pdf);
  }
}

vect_d Sampler::particle_birth(vect_d rands){
  // get indices
  int pdf_idx = at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx/num_e_groups;
  int e_idx = pdf_idx % num_e_groups;

  vect_d samp;
  vect_d xyz_rands;
  xyz_rands.push_back(rands[2]);
  xyz_rands.push_back(rands[3]);
  xyz_rands.push_back(rands[4]);
  MBCartVect pos = get_xyz(ve_idx, xyz_rands);

  samp.push_back(pos[0]); 
  samp.push_back(pos[1]); 
  samp.push_back(pos[2]); 
  samp.push_back(get_w(pdf_idx));
  samp.push_back(get_e(e_idx, rands[5]));
  return samp;
}

MBCartVect Sampler::get_xyz(int ve_idx, vect_d rands){
  double s = rands[0];
  double t = rands[1];
  double u = rands[2];

  if (ve_type == MBTET){
    if(s + t > 1){
      s = 1.0 - s;
      t = 1.0 - t;
    }
    if(s + t + u > 1){
      if(t + u > 1){
        double old_t = t;
        t = 1.0 - u;
        u = 1.0 - s - old_t;
      }else if (t + u <= 1){
        double old_s = s;
        s = 1.0 - t - u;
        u = old_s + t + u - 1;
      }
    }
  }

 return s*cart_sampler[ve_idx].x_vec + \
        t*cart_sampler[ve_idx].y_vec + \
        u*cart_sampler[ve_idx].z_vec + \
          cart_sampler[ve_idx].o_point;
}

double Sampler::get_e(int e_idx, double rand){
   double e_min = e_bounds[e_idx];
   double e_max = e_bounds[e_idx + 1];
   return rand * (e_max - e_min) + e_min;
}

double Sampler::get_w(int pdf_idx){
  if(mode == ANALOG){
    return 1.0;
  }else{
    return biased_weights[pdf_idx];
  }
}


////////////////////
// ALIAS TABLE STUFF
////////////////////


AliasTable::AliasTable(vect_d p){
  // Random-number sampling using the Walker-Vose alias method,
  // Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
  // M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
  // A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  vect_d small(n);
  vect_d large(n);
  int i, a, g;

  for(i=0; i<n; ++i) 
    p[i] *= n;

  // Set separate index lists for small and large probabilities:
  int n_s = 0;
  int n_l = 0;
  for(i=n-1; i>=0; --i) {
    // at variance from Schwarz, we revert the index order
    if(p[i] < 1)
      small[n_s++] = i;
    else
      large[n_l++] = i;
  }

  // Work through index lists
  while(n_s && n_l){
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

int AliasTable::sample_pdf(double rand1, double rand2){
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}
