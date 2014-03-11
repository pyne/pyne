// FluDAG/src/test/test_UnitNumber.cpp

#include "gtest/gtest.h"
#include "../UnitNumberManager.hpp"
#include <string>        // std::string
#include <cmath>         // std::abs
#include <iostream>      // std::cout
#include <sstream>       // std::ostringstream

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class UnitNumberManagerTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    // Use TEST_F for any test that assigns
    virtual void SetUp()
    {
        manager = NULL;
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete manager;
    }

  protected:
    // data needed for each test
    UnitNumberManager* manager;
};
// non-static data member initializers available c-11
std::string name1 = "Name1";
std::string name2 = "Name2";

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: UnitNumberManagerTest
//---------------------------------------------------------------------------//
// Test the inital value of the UnitNumberManager class memeber
TEST_F(UnitNumberManagerTest, InitializesCorrectly)
{
    manager = new UnitNumberManager();
    EXPECT_EQ(0, manager->getNumUnitsInUse());
}
//---------------------------------------------------------------------------//
// Test the expected error return for an empty name
TEST_F(UnitNumberManagerTest, HandlesErrorInput)
{
    manager = new UnitNumberManager();
    EXPECT_EQ(-1, manager->getUnitNumber(""));
}
//---------------------------------------------------------------------------//
// Test adding new names, fetching an existing one 
TEST_F(UnitNumberManagerTest, AddFetchNames)
{
    manager = new UnitNumberManager();
    EXPECT_EQ(UnitNumberManager::START_UNIT, manager->getUnitNumber(name1));
    EXPECT_EQ(1, manager->getNumUnitsInUse());

    EXPECT_EQ(-22, manager->getUnitNumber(name2));
    EXPECT_EQ(2, manager->getNumUnitsInUse());

    EXPECT_EQ(UnitNumberManager::START_UNIT, manager->getUnitNumber(name1));
    EXPECT_EQ(2, manager->getNumUnitsInUse());
}

//---------------------------------------------------------------------------//
// Test adding new names, fetching an existing one 
TEST_F(UnitNumberManagerTest, MaxUnits)
{
    manager = new UnitNumberManager();
    std::ostringstream oss;
    std::string s;
   
    int expectedNum = std::abs(UnitNumberManager::END_UNIT - UnitNumberManager::START_UNIT) + 1;
    for (int i=1; i<=expectedNum; ++i)
    {
        oss << "name" << i;
        s = oss.str();
        manager->getUnitNumber(s);
        EXPECT_EQ(i, manager->getNumUnitsInUse());
        oss.str("");
    }
    EXPECT_EQ(expectedNum, manager->getNumUnitsInUse());

    EXPECT_EQ(0, manager->getUnitNumber("OneTooMany"));
    // too many requests should not change the number of units in use
    EXPECT_EQ(expectedNum, manager->getNumUnitsInUse());
}
//---------------------------------------------------------------------------//
// end of FluDAG/src/test/test_UnitNumber.cpp
