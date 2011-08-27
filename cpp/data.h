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
}

#endif
