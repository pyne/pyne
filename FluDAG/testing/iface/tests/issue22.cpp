//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /tests/issue.cpp
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
#include "dagmc_utils.hpp"

std::string goodfile = "../iface/cases/test.h5m";
char *goodfileptr = &goodfile[0];
std::string prefix = "../iface/cases/";

TEST(Issue22_Test, MissingFile)
{
//    const std::string s;
  //  EXPECT_STREQ(NULL, s.c_str());

   std::string badFile = "MissingFile.h5m";
   // char *fileptr = &s[0];;
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
