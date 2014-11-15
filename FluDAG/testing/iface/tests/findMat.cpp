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
// #include "fludag_utils.h"
// #include "fluka_funcs.h"
// #include "dagmc_utils.hpp"

/* These 37 strings are predefined FLUKA materials. Any ASSIGNMAt of unique 
 * materials not on this list requires a MATERIAL card. */
std::string flukaMatStrings[] = {"BLCKHOLE", "VACUUM", "HYDROGEN",
"HELIUM", "BERYLLIU", "CARBON", "NITROGEN", "OXYGEN", "MAGNESIU",
"ALUMINUM", "IRON", "COPPER", "SILVER", "SILICON", "GOLD", "MERCURY",
"LEAD", "TANTALUM", "SODIUM", "ARGON", "CALCIUM", "TIN", "TUNGSTEN",
"TITANIUM", "NICKEL", "WATER", "POLYSTYR", "PLASCINT", "PMMA",
"BONECOMP", "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON",
"POLYETHY", "AIR"};

int NUM_FLUKA_MATS = 37;

/* Create a set out of the hardcoded string array. */
std::set<std::string> FLUKA_mat_set(flukaMatStrings, flukaMatStrings+NUM_FLUKA_MATS);

TEST(FindMat, matchString)
{

    std::set<std::string>::iterator it;
    for (int index = 0; index < NUM_FLUKA_MATS; index++)
    {
       it = FLUKA_mat_set.find(flukaMatStrings[index]);
       std::string findit = *it;
       std::string match = flukaMatStrings[index];
       ASSERT_STREQ(match.c_str(),findit.c_str());
          std::cout << "String at index " << index << " matches " <<  match << std::endl;
    }
}

