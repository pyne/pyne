// Enrichment Cascade
#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_cascade.h"
#endif

namespace pyne_enr = pyne::enrichment;

/*********************************/
/*** Enrichment Helper Classes ***/
/*********************************/
pyne_enr::Cascade::Cascade() {
  alpha = 0.0;
  Mstar = 0.0;

  j = 0;
  k = 0;

  N = 0.0;
  M = 0.0;

  x_feed_j = 0.0;
  x_prod_j = 0.0;
  x_tail_j = 0.0;

  mat_feed = pyne::Material();
  mat_prod = pyne::Material();
  mat_tail = pyne::Material();

  l_t_per_feed = 0.0;
  swu_per_feed = 0.0;
  swu_per_prod = 0.0;
}


pyne_enr::Cascade::~Cascade() {
}


void pyne_enr::Cascade::_reset_xjs() {
  // resets the key enriment member variables
  x_feed_j = mat_feed.comp[j];
  x_prod_j = mat_prod.comp[j];
  x_tail_j = mat_tail.comp[j];
}

