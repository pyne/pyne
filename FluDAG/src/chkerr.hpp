#ifndef PARSE_UTILS_HPP
#define PARSE_UTILS_HPP

#include <iostream>
#include <moab/Interface.hpp>

using namespace moab;

void chkerr( Interface* mbi, ErrorCode code, int line, const char* file );
void chkerr( Interface& mbi, ErrorCode code, int line, const char* file );

#define CHECKERR(M,C) chkerr(M,C,__LINE__,__FILE__)

#endif /* PARSE_UTILS_HPP */
