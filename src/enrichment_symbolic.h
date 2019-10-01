
/// \file enrichment_symbolic.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief A multicomponent enrichment cascade solver using
///     a symbolic solution to the mass flow rate equations.

/*********************************************************/
/***            Symbolic Enrichment Functions          ***/
/*** WARNING: This file is auto-generated.             ***/
/***                  DO NOT MODIFY!!!                 ***/
/*********************************************************/

#ifndef PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM
#define PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM

#include <math.h>

#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_cascade.h"
#endif

namespace pyne {
namespace enrichment {

  /// A multicomponent enrichment cascade solver using 
  /// a symbolic solution to the mass flow rate equations.
  /// \param orig_casc The original state of the cascade.
  /// \return A cascade solved for new N, M, and total flow
  ///         rates.
  Cascade solve_symbolic(Cascade & orig_casc);

// end enrichment
}
// end pyne
}

#endif
