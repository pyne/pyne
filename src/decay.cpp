// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//              WARNING
// This file has been auto generated
// Do not modify directly. You have
// been warned. This is that warning
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "decay.h"

namespace pyne {
namespace decayers {

std::map<int, double> decay(std::map<int, double> comp, double t) {
  // setup
  using std::map;
  int nuc;
  double out [3] = {};  // init to zero
  map<int, double> outcomp;
  
  // body
  map<int, double>::const_iterator it = comp.begin();
  for (; it != comp.end(); ++it) {
    switch (it->first) {
      case 10010000: {
        out[0] += it->second;
        break;
      } case 10020000: {
        out[1] += it->second;
        break;
      } case 10030000: {
        out[2] += (it->second) * (exp2(-2.572085049840012e-09*t));
        break;
      } default: {
        outcomp.insert(*it);
        break;
      }
    }
  }
  
  // cleanup
  for (int i = 0; i < 3; ++i)
    if (out[i] > 0.0)
      outcomp[all_nucs[i]] = out[i];
  return outcomp;
}

const int all_nucs [3] = {
  10010000, 10020000, 10030000
};

}  // namespace decayers
}  // namespace pyne