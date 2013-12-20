#ifndef DAGMC_MCNP_PARSE_FLU_HPP
#define DAGMC_MCNP_PARSE_FLU_HPP

#include <string>
#include "moab/Core.hpp"

  /*
   * Write records to a *.inp file for each volume tagged with material info and scoring requests.
   */
   void fludagwrite_assignma(moab::Interface* mbi, std::string matfile);  

#endif  /* DAGMC_MCNP_PARSE_FLU_H */
