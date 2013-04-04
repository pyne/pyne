//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /testing/issue.cpp
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

TEST(Issue22_Test, MissingFile)
{
//    const std::string s;
  //  EXPECT_STREQ(NULL, s.c_str());

   std::string s = "MissingFile.h5m";
   char *fileptr;
   fileptr = &s[0];
   EXPECT_FALSE(checkFile(fileptr));
}

