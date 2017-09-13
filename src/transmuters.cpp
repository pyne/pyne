extern "C" {
#include "cram.hpp"
}
#include "utils.h"
#include "transmuters.h"


std::map<int, double> pyne::transmuters::cram(std::vector<double>& A,
                                              const std::map<int, double>& n0,
                                              const int order) {
  using std::map;
  using std::vector;
  // Get intial condition vector
  vector<double> b (pyne_cram_transmute_info.n, 0.0);
  map<int, double>::const_iterator it;
  int i = -1;
  for (it = n0.begin(); it != n0.end(); ++it) {
    i = pyne_cram_transmute_nucid_to_i(it->first);
    if (i < 0) {
      continue;
    }
    b[i] = it->second;
  }

  // perform decay
  vector<double> x (pyne_cram_transmute_info.n);
  switch(order) {
    case 6:
      pyne_cram_expm_multiply6(A.data(), b.data(), x.data());
      break;
    case 8:
      pyne_cram_expm_multiply8(A.data(), b.data(), x.data());
      break;
    case 10:
      pyne_cram_expm_multiply10(A.data(), b.data(), x.data());
      break;
    case 12:
      pyne_cram_expm_multiply12(A.data(), b.data(), x.data());
      break;
    case 14:
      pyne_cram_expm_multiply14(A.data(), b.data(), x.data());
      break;
    case 16:
      pyne_cram_expm_multiply16(A.data(), b.data(), x.data());
      break;
    case 18:
      pyne_cram_expm_multiply18(A.data(), b.data(), x.data());
      break;
    default:
      throw pyne::ValueError("Order selected not available for CRAM, please use"
                             " order 6, 8, 10, 12, 14, 16, or 18.");
      break;
  }

  // convert back to map
  map<int, double> n1;
  for (i=0; i < pyne_cram_transmute_info.n; ++i) {
    if (x[i] > 0.0) {
      n1[(pyne_cram_transmute_info.nucids)[i]] = x[i];
    }
  }
  return n1;
}