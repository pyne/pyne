//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /testing/iface/test/g1_test.cpp
 * \author Julie C. Zachman
 * \date   Wed Apr 3 
 * \brief  Unit Tests for Issue
 * \note   Test
 */
//---------------------------------------------------------------------------//
// $Id:
//---------------------------------------------------------------------------//
#include "gtest/gtest.h"
#include "fludag_utils.h"
#include "fluka_funcs.h"

std::string goodfile = "../iface/cases/test.h5m";
char *goodfileptr = &goodfile[0];
std::string prefix = "../iface/cases/";
/*
void g1wr(double& pSx, 
          double& pSy, 
          double& pSz, 
          double* pV,
          int& oldReg,         // pass through
          const int& oldLttc,  // ignore
          double& propStep,    // .
          int& nascFlag,       // .
          double& retStep,     // reset in this method
          int& newReg,         // return from callee
          double& saf,         // ignore
          int& newLttc,        // .
          int& LttcFlag,       // . 
          double* sLt,         // .
          int* jrLt)           // .
*/
double pSx;
double pSy;
double pSz;

double pV;
int oldReg;

TEST(g1_Test, MissingFile)
{
   std::string badFile = "MissingFile.h5m";
   EXPECT_FALSE(isFileReadable(badFile));

   EXPECT_TRUE(isFileReadable(goodfile));
}

TEST(Issue22_Test, NonEmptyFile)
{
   std::string s = prefix.append("EmptyFile.h5m");
   char *fileptr = &s[0];;
   EXPECT_FALSE(isFileNonEmpty(fileptr));

   EXPECT_TRUE(isFileNonEmpty(goodfileptr));
}

TEST(Issue22_Test, Combo)
{
   std::string s = "MissingFile.h5m";
   char *fileptr = &s[0];
   EXPECT_FALSE(checkFile(s));

   s = prefix.append("EmptyFile.h5m");
   fileptr = &s[0];
   EXPECT_FALSE(checkFile(fileptr));

   EXPECT_TRUE(checkFile(goodfileptr));
}

TEST(Issue22_Test, MalformedData)
{
   
   std::string s = prefix.append("SomeTextFile.h5m");
   char *fileptr = &s[0];
   std::string fail_msg = "DAGMC failed to read input file: ";

   ASSERT_DEATH({cpp_dagmcinit(fileptr,0,1);},fail_msg); 
}
