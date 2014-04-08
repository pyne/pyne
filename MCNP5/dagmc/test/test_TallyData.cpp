// MCNP5/dagmc/test/test_TallyData.cpp

#include <vector>

#include "gtest/gtest.h"

#include "moab/CartVect.hpp"

#include "../TallyData.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class TallyDataTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {

        // 'Normal' case
        tallyData1 = new TallyData(1,false);
        
        // For looking at total_energy_bin
        num_energy_bins = 5;
        tallyData2 = new TallyData(num_energy_bins, true);
    
        tallyData3 = new TallyData(9, false); 
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete tallyData1;
        delete tallyData2;
        delete tallyData3;
    }

  protected:
    // Clarifying global 
    unsigned int num_energy_bins;

   TallyData* tallyData1;
   TallyData* tallyData2;
   TallyData* tallyData3;
};
//---------------------------------------------------------------------------//
// SIMPLE TESTS
//---------------------------------------------------------------------------//
TEST(TallyDataInputTest, TotalEnergyBinLogic)
{
   TallyData tallyData1(1,true);
   TallyData tallyData2(1,false);
   TallyData tallyData3(7,true);
   TallyData tallyData4(8,false);
    
   EXPECT_EQ(1, tallyData1.get_num_energy_bins());
   EXPECT_FALSE(tallyData1.has_total_energy_bin());
   EXPECT_EQ(1, tallyData2.get_num_energy_bins());
   EXPECT_FALSE(tallyData2.has_total_energy_bin());

   EXPECT_EQ(8, tallyData3.get_num_energy_bins());
   EXPECT_TRUE(tallyData3.has_total_energy_bin());

   EXPECT_EQ(8, tallyData4.get_num_energy_bins());
   EXPECT_FALSE(tallyData4.has_total_energy_bin());
}

//---------------------------------------------------------------------------//
TEST(TallyDataDeathTest, ZeroEnergyBins)
{
   EXPECT_EXIT(TallyData(0, true),  ::testing::ExitedWithCode(EXIT_FAILURE), "");
   EXPECT_EXIT(TallyData(0, false), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyDataTest
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
void ResizeDataArraysHelper(TallyData* tallyData, 
                          unsigned int num_tally_points, 
                          unsigned int expected_length)
{
     int length;

     double *tallyArray1, *tallyArray2, *tallyArray3;

     tallyData->resize_data_arrays(num_tally_points);

     tallyArray1 = tallyData->get_tally_data(length);
     EXPECT_EQ(expected_length, length);

     tallyArray2 = tallyData->get_scratch_data(length);
     EXPECT_EQ(expected_length, length);

     tallyArray3 = tallyData->get_error_data(length);
     EXPECT_EQ(expected_length, length);
     // test arrays initialized
     for (int i=0; i<length; i++)
     {
        EXPECT_DOUBLE_EQ(0.0, tallyArray1[i]);
        EXPECT_DOUBLE_EQ(0.0, tallyArray2[i]);
        EXPECT_DOUBLE_EQ(0.0, tallyArray3[i]);
     }
}

//---------------------------------------------------------------------------//
TEST_F(TallyDataTest, ResizeDataArrays)
{
     // tallyData1 has 1 energy bin, with total=false
     ResizeDataArraysHelper(tallyData1, 1, 1);
     ResizeDataArraysHelper(tallyData1, 3, 3);
    
     // tallyData2 has 5 energy bins, with total=true
     ResizeDataArraysHelper(tallyData2, 1, 6);
     ResizeDataArraysHelper(tallyData2, 3, 18);

     // tallyData3 has 9 energy bins, with total=false
     ResizeDataArraysHelper(tallyData3, 1, 9);
     ResizeDataArraysHelper(tallyData3, 3, 27);
}

//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_TallyData.cpp
