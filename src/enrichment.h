/// \file enrichment.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Top-level enrichment functionality.

#ifndef PYNE_B3ANNCKDQ5HEJLI33RPZPDNX6A
#define PYNE_B3ANNCKDQ5HEJLI33RPZPDNX6A

#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_symbolic.h"
#endif

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

namespace pyne {
//! Enrichment Component Class and Functions
namespace enrichment {

  /// Greates a cascade instance with default values for a uranium enrichment.
  Cascade _fill_default_uranium_cascade();
  /// a cascade instance with default values for a uranium enrichment.
  extern Cascade default_uranium_cascade;

  /// Computes the feed per product mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double feed_per_prod(double x_feed, double x_prod, double x_tail);
  /// Computes the feed per tails mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double feed_per_tail(double x_feed, double x_prod, double x_tail);
  /// Computes the product per feed mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double prod_per_feed(double x_feed, double x_prod, double x_tail);
  /// Computes the product per tails mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double prod_per_tail(double x_feed, double x_prod, double x_tail);
  /// Computes the tails per feed mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double tail_per_feed(double x_feed, double x_prod, double x_tail);
  /// Computes the tails per product mass ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double tail_per_prod(double x_feed, double x_prod, double x_tail);
  /// Computes the value or separation potential of an assay \a x.
  double value_func(double x);
  /// Computes the swu per feed ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double swu_per_feed(double x_feed, double x_prod, double x_tail);
  /// Computes the swu per product ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double swu_per_prod(double x_feed, double x_prod, double x_tail);
  /// Computes the swu per tails ratio for feed, product, and tails
  /// enrichments \a x_feed, \a x_prod, and \a x_tails.
  double swu_per_tail(double x_feed, double x_prod, double x_tail);

  /// Computes the nuclide-specific stage separation factor from the
  /// overall stage separation factor \a alpha, the key mass \a Mstar,
  /// and the nulide's atomic mass \a M_i.
  double alphastar_i(double alpha, double Mstar, double M_i);

  /// \{ Numeric Solver Functions
  /// \{
  /// These functions enable the numeric cascade solver.

  /// Finds the total flow rate (L) over the feed flow rate (F), the number of
  /// enriching stages (N), and the number of stripping stages (M).
  /// \param orig_casc Original input cascade.
  /// \param tolerance Maximum numerical error allowed in L/F, N, and M.
  /// \param max_iter Maximum number of iterations for to perform.
  /// \return A cascade whose N & M coorespond to the L/F value.
  Cascade solve_numeric(Cascade & orig_casc, double tolerance=1.0E-7,
                                             int max_iter=100);
  /// So,ves for valid stage numbers N &nd M of a casc.
  /// \param casc Cascade instance, modified in-place.
  /// \param tolerance Maximum numerical error allowed in N and M.
  void _recompute_nm(Cascade & casc, double tolerance=1.0E-7);
  /// This function takes a given initial guess number of enriching and stripping
  /// stages for a given composition of fuel with a given jth key component, knowing
  /// the values that are desired in both Product and Tails streams.  Having this it
  /// solves for what the actual N and M stage numbers are and also what the product
  /// and waste streams compositions are.
  /// \param casc Cascade instance, modified in-place.
  void _recompute_prod_tail_mats(Cascade & casc);
  /// This function solves the whole system of equations.  It uses
  /// _recompute_prod_tail_mats() to find the roots for the enriching and stripping
  /// stage numbers.  It then checks to see if the product and waste streams meet
  /// their target enrichments for the jth component like they should.  If they don't
  /// then it trys other values of N and M varied by the Secant method.
  /// \param casc Input cascade.
  /// \param tolerance Maximum numerical error allowed in L/F, N, and M.
  /// \param max_iter Maximum number of iterations for to perform.
  /// \return A cascade whose N & M coorespond to the L/F value.
  Cascade _norm_comp_secant(Cascade & casc, double tolerance=1.0E-7, int max_iter=100);
  /// Solves for a stage separative power relevant to the ith component
  /// per unit of flow G.  This is taken from Equation 31 divided by G
  /// from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
  /// 'On a Criterion Efficiency for Multi-Isotope Mixtures Separation',
  /// Separation Science and Technology, 34:3, 343 - 357"
  /// To link to this article: DOI: 10.1081/SS-100100654
  /// URL: http://dx.doi.org/10.1081/SS-100100654
  /// \param casc Input cascade.
  /// \param i nuclide in id form.
  double _deltaU_i_OverG(Cascade & casc, int i);
  /// \}

  /// \name Multicomponent Functions
  /// \{
  /// Finds a value of Mstar by minimzing the seperative power.
  /// Note that Mstar on \a orig_casc represents an intial guess at what Mstar might
  /// be. This is the final function that actually solves for an optimized M* that
  /// makes the cascade!
  /// \param orig_casc Original input cascade.
  /// \param solver flag for solver to use, may be 'symbolic' or 'numeric'.
  /// \param tolerance Maximum numerical error allowed in L/F, N, and M.
  /// \param max_iter Maximum number of iterations for to perform.
  /// \return A cascade whose N & M coorespond to the L/F value.
  Cascade multicomponent(Cascade & orig_casc, char * solver,
                         double tolerance=1.0E-7, int max_iter=100);
  Cascade multicomponent(Cascade & orig_casc, std::string solver="symbolic",
                         double tolerance=1.0E-7, int max_iter=100);
  /// \}

  /// Custom exception for when an enrichment solver has entered an infinite loop.
  class EnrichmentInfiniteLoopError: public std::exception
  {
    /// Returns a moderately helpful error message.
    virtual const char* what() const throw()
    {
      return "Inifinite loop found while calculating enrichment cascade.";
    }
  };

  /// Custom exception for when an enrichment solver has reached its maximum
  /// number of iterations.
  class EnrichmentIterationLimit: public std::exception
  {
    /// Returns a moderately helpful error message.
    virtual const char* what() const throw()
    {
      return "Iteration limit hit durring enrichment calculation.";
    }
  };

  /// Custom exception for when an enrichment solver iteration has produced a NaN.
  class EnrichmentIterationNaN: public std::exception
  {
    /// Returns a moderately helpful error message.
    virtual const char* what() const throw()
    {
      return "Iteration has hit a point where some values are not-a-number.";
    }
  };

// end enrichment
}
// end pyne
}

#endif
