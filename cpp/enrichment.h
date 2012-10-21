#if !defined(_PYNE_ENRICHMENT_)
#define _PYNE_ENRICHMENT_

#include "enrichment_symbolic.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

namespace pyne {
namespace enrichment {

  Cascade _fill_default_uranium_cascade();
  extern Cascade default_uranium_cascade;

  double prod_per_feed(double, double, double);
  double tail_per_feed(double, double, double);
  double tail_per_prod(double, double, double);

  double alphastar_i(double, double, double);

  void _recompute_nm(Cascade &, double=1.0E-7);
  void _recompute_prod_tail_mats(Cascade &);
  Cascade _norm_comp_secant(Cascade &, double=1.0E-7, int=100);
  double _deltaU_i_OverG(Cascade &, int);

  Cascade ltot_per_feed(Cascade &, double=1.0E-7, int=100);
  Cascade multicomponent(Cascade &, double=1.0E-7, int=100);

  /******************/
  /*** Exceptions ***/
  /******************/
  class EnrichmentInfiniteLoopError: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Inifinite loop found while calculating enrichment cascade.";
    };
  };

  class EnrichmentIterationLimit: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Iteration limit hit durring enrichment calculation.";
    };
  };

  class EnrichmentIterationNaN: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Iteration has hit a point where some values are not-a-number.";
    };
  };

// end enrichment
};
// end pyne
};

#endif


