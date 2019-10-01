#ifndef JSON_IS_AMALGAMATION
  #define JSON_IS_AMALGAMATION
#endif

#ifndef PYNE_IS_AMALGAMATED
#include "json-forwards.h"
#include "json.h"
//#include "dagmc_bridge.h"  // figure this out later...
extern "C" {
#include "cram.hpp"
}
#include "transmuters.h"
#include "data.h"
#include "decay.h"
#include "enrichment_cascade.h"
#include "enrichment.h"
#include "enrichment_symbolic.h"
#include "extra_types.h"
#include "h5wrap.h"
#include "material.h"
#include "nucname.h"
#include "rxname.h"
#include "tally.h"
#include "utils.h"
#endif
