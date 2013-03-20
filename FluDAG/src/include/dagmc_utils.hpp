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

// features provided by obb_analysis.cpp

// ErrorCode obbvis_create( DagMC& dag, std::vector<int> &volumes, int grid, std::string& filename );
// ErrorCode obbstat_write( DagMC& dag, std::vector<int> &volumes,
//                          std::vector<std::string> &properties, std::ostream& out );


#endif /* DAGMC_UTILS_H */

