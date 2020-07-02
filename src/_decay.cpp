#ifdef PYNE_DECAY_IS_DUMMY
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//              WARNING
// This file has been auto generated
// Do not modify directly. You have
// been warned. This is that warning
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef PYNE_IS_AMALGAMATED
#include "decay.h"
#endif

namespace pyne {
namespace decayers {

void decay_h(double t, std::map<int, double>::const_iterator &it, std::map<int, double> &outcomp, double (&out)[4]) {
  //using std::exp2;
  switch (it->first) {
    case 10010000: {
      out[0] += it->second;
      break;
    } case 10020000: {
      out[1] += it->second;
      break;
    } case 10030000: {
      double b0 = exp2(-2.572085e-09*t);
      out[2] += (it->second) * (b0);
      out[3] += (it->second) * (-1.000000e+00*b0 + 1.0);
      break;
    } default: {
      outcomp.insert(*it);
      break;
    }
  }
}

void decay_he(double t, std::map<int, double>::const_iterator &it, std::map<int, double> &outcomp, double (&out)[4]) {
  //using std::exp2;
  switch (it->first) {
    case 20030000: {
      out[3] += it->second;
      break;
    } default: {
      outcomp.insert(*it);
      break;
    }
  }
}

std::map<int, double> decay(std::map<int, double> comp, double t) {
  // setup
  using std::map;
  int nuc;
  int i = 0;
  double out [4] = {};  // init to zero
  map<int, double> outcomp;

  // body
  map<int, double>::const_iterator it = comp.begin();
  for (; it != comp.end(); ++it) {
    switch (nucname::znum(it->first)) {
      case 1:
        decay_h(t, it, outcomp, out);
        break;
      case 2:
        decay_he(t, it, outcomp, out);
        break;
      default:
        outcomp.insert(*it);
        break;
    }
  }

  // cleanup
  for (i = 0; i < 4; ++i)
    if (out[i] > 0.0)
      outcomp[all_nucs[i]] = out[i];
  return outcomp;
}

const int all_nucs [4] = {
  10010000, 10020000, 10030000, 20030000
};

}  // namespace decayers
}  // namespace pyne

#endif  // PYNE_DECAY_IS_DUMMY
