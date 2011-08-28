// Implements basic nuclear data functions.

#ifndef _PYNE_DATA_
#define _PYNE_DATA_
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#include "H5Cpp.h"
#include "h5wrap.h"

#include "extra_types.h"

#include "pyne.h"
#include "nucname.h"

namespace pyne
{
  extern std::string NUC_DATA_PATH;

  /****************************/
  /*** nuc_weight Functions ***/
  /****************************/
  extern std::map<int, double> nuc_weight_map;

  typedef struct atomic_weight_struct {
    char nuc_name[6];
    int nuc_zz;
    double mass;
    double error;
    double abund;
  } atomic_weight_struct;

  void _load_nuc_weight_map();
  double nuc_weight(int);
  double nuc_weight(char *);
  double nuc_weight(std::string);


  /***********************************/
  /*** scattering length functions ***/
  /***********************************/
  extern std::map<int, extra_types::complex_t> b_coherent_map;
  extern std::map<int, extra_types::complex_t> b_incoherent_map;
  extern std::map<int, double> b_map;

  typedef struct scattering_lengths_struct {
    char nuc_name[6];
    int nuc_zz;
    extra_types::complex_t b_coherent;
    extra_types::complex_t b_incoherent;
    double xs_coherent;
    double xs_incoherent;
    double xs;
  } scattering_lengths_struct;

  void _load_scattering_lengths();

  extra_types::complex_t b_coherent(int);
  extra_types::complex_t b_coherent(char *);
  extra_types::complex_t b_coherent(std::string);

  extra_types::complex_t b_incoherent(int);
  extra_types::complex_t b_incoherent(char *);
  extra_types::complex_t b_incoherent(std::string);

  double b(int);
  double b(char *);
  double b(std::string);
}

#endif
