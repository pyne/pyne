#if !defined(_PYNE_ENRICHMENT_)
#define _PYNE_ENRICHMENT_

#include "pyne.h"
#include "nucname.h"
#include "data.h"
#include "material.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

namespace pyne {
namespace enrichment {

  class Cascade 
  {
    /** Set of physical parameters used to specify an enrichment cascade. **/

  public:
    // Constructors
    Cascade();
    ~Cascade();

    // Attributes
    double alpha; // stage separation factor
    double Mstar; // mass separation factor

    int j; // Component to enrich (U-235), zzaaam form
    int k; // Component to de-enrich, or strip (U-238), zzaaam form

    double N; // number of enriching stages
    double M; // number of stripping stages

    double x_feed_j; // enrichment of the jth isotope in the feed stream
    double x_prod_j; // enrichment of the jth isotope in the product stream
    double x_tail_j; // enrichment of the jth isotope in the tails stream

    pyne::Material mat_feed; // feed material
    pyne::Material mat_prod; // product material
    pyne::Material mat_feed; // tails material
  };


  Cascade fillUraniumEnrichmentDefaults();
  extern Cascade UraniumEnrichmentDefaults;


  class Enrichment : public FCComp
  {
  // Reprocessing class
  public:
    // Reprocessing Constructors
    Enrichment ();
    Enrichment (std::string);
    Enrichment (Cascade, std::string = "");
    ~Enrichment ();

    // Public data
    double alpha_0;           // specify on init.
    double Mstar_0;           // specify on init.
    double Mstar;             // Current Mstar
    pyne::Material mat_tail;  // Tails Stream

    // key isotopic info
    int j;          // The jth isotope is the key, in zzaaam form, must be in mat_feed.
    int k;          // The kth isotope is the other key to separate j away from.
    double xP_j;    // Product enrichment of jth isotope
    double xW_j;    // Tails/Tails enrichment of the jth isotope

    // Stage info
    double N;       // N Enriching Stages
    double M;       // M Stripping Stages
    double N0;      // initial guess of N-stages
    double M0;      // initial guess of M-stages

    // Flow Rates
    double TotalPerFeed;    // Total flow rate per feed rate.
    double swu_per_feed;      // This is the SWU for 1 kg of Feed material.
    double swu_per_prod;   // This is the SWU for 1 kg of Product material.



    double prod_per_feed(double, double, double);
    double tail_per_feed(double, double, double);
    double tail_per_prod(double, double, double);

    double alphastar_i(double, double, double);
    double Ei (double, double);
    double Si (double, double);
    void FindNM();

    double xP_i(int);
    double xW_i(int);
    void SolveNM();
    void _norm_comp_secant(double=1.0E-7);
    void _norm_comp_other(double=1.0E-7);
    double deltaU_i_OverG(int);
    void ltot_per_feed();
    void multicomponent(double, double=1.0E-7);
  };


  /******************/
  /*** Exceptions ***/
  /******************/
  class EnrichmentInfiniteLoopError: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Inifinite loop found while calculating enrichment cascade!  Breaking...";
    };
  };


  class EnrichmentIterationLimit: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Iteration limit hit durring enrichment calculation!  Breaking...";
    };
  };


  class EnrichmentIterationNaN: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Iteration has hit a point where some values are not-a-number!  Breaking...";
    };
  };

// end enrichment
};
// end pyne
};

#endif


