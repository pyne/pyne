#ifndef DAGMC_UTILS_H
#define DAGMC_UTILS_H

#include <iostream>
#include <stdlib.h>
#include <moab/Interface.hpp>

using namespace moab;
namespace moab
{
    class DagMC;
}

/**
 * Generic halt-on-error error checking
 */
void chkerr( Interface* mbi, ErrorCode code, int line, const char* file );
void chkerr( Interface& mbi, ErrorCode code, int line, const char* file );
void chkerr( DagMC& dag, ErrorCode code, int line, const char* file );

#define CHECKERR(M,C) chkerr(M,C,__LINE__,__FILE__)

#endif /* DAGMC_UTILS_H */

