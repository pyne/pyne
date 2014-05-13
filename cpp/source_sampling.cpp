#include "sampling.hpp"
//#define MBI moab_instance()
//#define SI Sampling::instance()

Sampling *Sampling::instance_ = NULL;

/*
 ( FORTRAN API
*/

void mcnp_sampling_setup_(bool* analog){
  SI->sampling_setup((char*)&"tagged_unstr.h5m", (char*)&"phtn_src2", (char*)&"e_bounds_file", *analog);
}

void fsampling_setup_(char* file_name, char* src_tag_name, char* e_bounds_tag_name, bool* analog){
  SI->sampling_setup(file_name, src_tag_name, e_bounds_tag_name, analog);
}

void fsampling_setup2_(char* file_name, char* src_tag_name, char* e_bounds_file_name, bool* analog, char* bias_tag_name){
  SI->sampling_setup(file_name, src_tag_name, e_bounds_file_name, analog, bias_tag_name);
}

void fparticle_birth_(double* rands, double* x, double* y, double* z, double* e, double* w){
  SI->particle_birth(rands, x, y, z, e, w);
}

/*
 * IPA NARTROF
 */

void Sampling::create_instance(MBInterface *mb_impl)
{
  if (NULL == mb_impl) mb_impl = new MBCore();
  instance_ = new Sampling(mb_impl);
}

Sampling::Sampling(MBInterface *mb_impl)
   : mbImpl(mb_impl), uniform(false){}


void Sampling::sampling_setup(char* file_name, char* src_tag_name, char* e_bounds_file, bool _analog){

  analog = _analog;
  if(analog == false)
    uniform = true;

  char* bias_tag_name = NULL;
  sampling_setup(file_name, src_tag_name, e_bounds_file, analog, bias_tag_name);
}


void Sampling::sampling_setup(char* file_name, char* _src_tag_name, char* e_bounds_file, bool analog, char* _bias_tag_name){

  if(analog == true && _bias_tag_name != NULL){
    throw std::invalid_argument("bias_tag_name should not be specified for analog sampling");
  }

  // If this function is not being called from the other overload of this 
  // function, then a bias tag has been specified. 
  src_tag_name = _src_tag_name;
  bias_tag_name = _bias_tag_name;

  MBErrorCode rval;
  MBEntityHandle loaded_file_set;
  rval = MBI->create_meshset(MESHSET_SET, loaded_file_set );
  //assert( rval == MB_SUCCESS );
  rval = MBI->load_file( file_name, &loaded_file_set );
  //assert( rval == MB_SUCCESS );
  MBRange ves;
  rval = MBI->get_entities_by_dimension(loaded_file_set, 3, ves);

  int num_hex, num_tet;
  rval = MBI->get_number_entities_by_type(loaded_file_set, MBHEX, num_hex);
  rval = MBI->get_number_entities_by_type(loaded_file_set, MBTET, num_tet);
  if(num_hex == ves.size()){
    ve_type = MBHEX;
    verts_per_ve = 8;
  } else if (num_tet == ves.size()){
    ve_type = MBTET;
    verts_per_ve = 4;
  }
  else
   throw std::invalid_argument("Mesh file must contain only tets or hexes");

  std::vector<double> volumes(ves.size());
  get_mesh_geom_data(ves, volumes);
  get_mesh_tag_data(ves, volumes);
  get_e_bounds_data(e_bounds_file);
}

void Sampling::get_mesh_geom_data(MBRange ves, std::vector<double> &volumes){
  MBErrorCode rval;
  std::vector<MBEntityHandle> connect;
  rval = MBI->get_connectivity_by_type(ve_type, connect);
  double coords[verts_per_ve*3];
  int i;
  for(i=0; i<ves.size(); ++i){
    rval=MBI->get_coords(&connect[verts_per_ve*i], verts_per_ve, &coords[0]);
    volumes[i] = measure(ve_type, verts_per_ve, &coords[0]);

    if(ve_type == MBHEX){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[9], coords[10], coords[11]);
      MBCartVect z(coords[12], coords[13], coords[14]);
      vector_points vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
   }else if (ve_type == MBTET){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[6], coords[7], coords[8]);
      MBCartVect z(coords[9], coords[10], coords[11]);
      vector_points vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
    }
  }
}

void Sampling::get_mesh_tag_data(MBRange ves, std::vector<double>volumes){
  MBErrorCode rval;
  MBTag src_tag;
  rval = MBI->tag_get_handle(src_tag_name, moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, src_tag);
  // THIS ASSERT FAILS because we do not know number of energy groups a priori.
  //assert( rval == MB_SUCCESS );
  int src_tag_size;
  rval = MBI->tag_get_bytes(src_tag, *(&src_tag_size));
  //assert( rval == MB_SUCCESS );
  num_e_groups = src_tag_size/sizeof(double);

  std::vector<double> pdf(ves.size()*num_e_groups); 
  rval = MBI->tag_get_data(src_tag, ves, &pdf[0]);
  //assert( rval == MB_SUCCESS );

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

  if(analog == false){
    std::vector<double> bias_pdf(ves.size()*num_e_groups); 
    if(uniform == true){
      for(i=0; i<ves.size(); ++i){
        for(j=0; j<num_e_groups; ++j)
           bias_pdf[i*num_e_groups + j] =  volumes[i];
      }
    }else if(uniform == false){
      MBTag bias_tag;
      rval = MBI->tag_get_handle(bias_tag_name, moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, bias_tag);
      int bias_tag_size;
      rval = MBI->tag_get_bytes(bias_tag, *(&bias_tag_size));
      int num_bias_groups = bias_tag_size/sizeof(double);

      if (num_bias_groups == num_e_groups){
        rval = MBI->tag_get_data(bias_tag, ves, &bias_pdf[0]);
        for(i=0; i<ves.size(); ++i){
          for(j=0; j<num_e_groups; ++j)
             bias_pdf[i*num_e_groups + j] *=  volumes[i];
        }
      }else if(num_bias_groups == 1){
        std::vector<double> spacial_pdf(ves.size()); 
        rval = MBI->tag_get_data(bias_tag, ves, &spacial_pdf[0]);
        for(i=0; i<ves.size(); ++i){
          bias_pdf[i*num_e_groups + j] *=  spacial_pdf[i]*volumes[i];
        }
      }
      else
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
  }else if(analog == true){
    at = new AliasTable(pdf);
  }

}

void Sampling::get_e_bounds_data(char* e_bounds_file){
/* E_TAG STUFF
  MBTag e_tag;
  std::cout << e_bounds_file_name << std::endl;
  //rval = MBI->tag_get_handle(e_bounds_file_name, 3, MB_TYPE_DOUBLE, e_tag);
  rval = MBI->tag_get_handle("e_bounds2", 3, MB_TYPE_DOUBLE, e_tag);
  std::cout << "error code: " << rval << std::endl;
  std::cout << "e_tag"<< e_tag<< std::endl;

  MBRange entities;
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &e_tag, NULL, 1, entities);
  //rval = MBI->get_entities_by_type_and_tag(0, MBHEX, &src_tag, NULL, 1, entities);
  std::cout << "error code: " << rval << std::endl;
  std::cout << entities << std::endl; 
  std::cout <<"iterator length: "<<entities.size() << std::endl;

  for (MBRange::const_iterator s = entities.begin(); s != entities.end(); ++s) {
      MBEntityHandle set = *s;
      rval = MBI->tag_get_data( e_tag, &set, 1, &e_bounds[0] );
      std::cout << rval << std::endl;
      std::cout << "hello" << std::endl;
  }

 // std::cout << rval << std::endl;
 // assert(rval == MB_SUCCESS);
 // rval = MBI->tag_get_data(e_tag, entities, 1, &e_bounds[0]);
  std::cout << rval << std::endl;
  //assert(rval == MB_SUCCESS);
 // std::cout << e_bounds[0] << std::endl;
  //std::cout << e_bounds[0] << e_bounds[1] << e_bounds[2] << std::endl;
  */
  std::ifstream inputFile(e_bounds_file);
  // test file open   
  if (inputFile) {        
    double value;
    while (inputFile >> value)
      e_bounds.push_back(value);
  }

/*
    std::ifstream input (e_bounds_file);
    std::string lineData;

    while(getline(input, lineData)){
        double d;
        std::vector<double> row;
        std::stringstream lineStream(lineData);

        while (lineStream >> d)
            row.push_back(d);

        e_bounds.push_back(row);
    }
*/
}


void Sampling::particle_birth(double* rands, double* x, double *y, double *z, double* e, double* w){
  // get indices
  int pdf_idx = at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx/num_e_groups;
  int e_idx = pdf_idx % num_e_groups;
  
  // get x, y, z
  double e_rand = rands[5];
  get_e(e_idx, e_rand, e);

  double xyz_rands[3] = {rands[2], rands[3], rands[4]};
  get_xyz(ve_idx, xyz_rands, x, y, z);
  get_w(pdf_idx, w);
}

void Sampling::get_xyz(int ve_idx, double* rands, double* x, double* y, double* z){

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

  MBCartVect birth_location = s*cart_sampler[ve_idx].x_vec + \
                              t*cart_sampler[ve_idx].y_vec + \
                              u*cart_sampler[ve_idx].z_vec + \
                                cart_sampler[ve_idx].o_point;
  *x = birth_location[0];
  *y = birth_location[1];
  *z = birth_location[2];
}

void Sampling::get_e(int e_idx, double rand, double* e){
   double e_min = e_bounds[e_idx];
   double e_max = e_bounds[e_idx + 1];
   *e = rand * (e_max - e_min) + e_min;
}

void Sampling::get_w(int pdf_idx, double* w){
  if(analog == true){
    *w = 1.0;
  }else{
    *w = biased_weights[pdf_idx];
  }
}

Sampling::AliasTable::AliasTable(std::vector<double> p){
  // Random-number sampling using the Walker-Vose alias method,
  // Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
  // M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
  // A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  std::vector<double> small(n);
  std::vector<double> large(n);
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

int Sampling::AliasTable::sample_pdf(double rand1, double rand2){
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}
