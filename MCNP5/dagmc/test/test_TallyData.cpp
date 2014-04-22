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
        tallyData2 = new TallyData(5, true);
    
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

   TallyData* tallyData1;
   TallyData* tallyData2;
   TallyData* tallyData3;
};

//---------------------------------------------------------------------------//
void fillTally(TallyData& tallyData)
{
        double fill_tally[]   = {1.2, -3.4, 5.6, 0.0, 8.9, 7.8};
        double fill_error[]   = {0.2, -0.4, 0.06, 0.0, .089, .078};
        double fill_scratch[] = {10.2, -30.4, 50.6, 0.0, 80.9, 70.8};
        // Set it up
        int length;
        double* tally_data = tallyData.get_tally_data(length); 
        double* error_data = tallyData.get_error_data(length);
        double* scratch_data = tallyData.get_scratch_data(length);
        for (int i=0; i<6; i++)
        { 
            tally_data[i] = fill_tally[i];
            error_data[i] = fill_error[i];
            scratch_data[i] = fill_scratch[i];
        } 
}
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
   ::testing::FLAGS_gtest_death_test_style="threadsafe";
   EXPECT_DEATH(TallyData(0, true),  "");
   EXPECT_DEATH(TallyData(0, false),  "");
}
//---------------------------------------------------------------------------//
TEST(FilledTallyTest, ZeroTallyData)
{
   // Set it up
   TallyData tallyData(2, false);
   tallyData.resize_data_arrays(3);
   ///////////////////////////////////////////////////////////////
   // Fill it with some random data
   fillTally(tallyData);
 
   int length;
   double* tally_ary = tallyData.get_tally_data(length);
   EXPECT_EQ(6, length);

   double* error_ary   = tallyData.get_error_data(length);
   EXPECT_EQ(6, length);

   double* scratch_ary = tallyData.get_scratch_data(length);
   EXPECT_EQ(6, length);
   ///////////////////////////////////////////////////////////////
  
   // Clear the data
   tallyData.zero_tally_data();
   
   // Check the zeroing worked
   for (int i=0; i<length; i++)
   { 
       EXPECT_DOUBLE_EQ(0.0, tally_ary[i]); 
       EXPECT_DOUBLE_EQ(0.0, error_ary[i]); 
       EXPECT_DOUBLE_EQ(0.0, scratch_ary[i]); 
   } 
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyDataTest
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
TEST(FilledTallyTest, GetDataTest)
{
   TallyData tallyData(2, true);
   tallyData.resize_data_arrays(2);
   fillTally(tallyData);

   std::pair<double, double> returnPair;
   // Tally Point 0: two energy bins and energy total
   returnPair = tallyData.get_data(0,0);
   EXPECT_DOUBLE_EQ(1.2,returnPair.first);
   EXPECT_DOUBLE_EQ(0.2,returnPair.second);

   returnPair = tallyData.get_data(0,1);
   EXPECT_DOUBLE_EQ(-3.4,returnPair.first);
   EXPECT_DOUBLE_EQ(-0.4,returnPair.second);

   returnPair = tallyData.get_data(0,2);
   EXPECT_DOUBLE_EQ(5.6,returnPair.first);
   EXPECT_DOUBLE_EQ(0.06,returnPair.second);
   
   // Tally Point 1: two energy bins and energy total
   returnPair = tallyData.get_data(1,0);
   EXPECT_DOUBLE_EQ(0.0,returnPair.first);
   EXPECT_DOUBLE_EQ(0.0,returnPair.second);

   returnPair = tallyData.get_data(1,1);
   EXPECT_DOUBLE_EQ(8.9,returnPair.first);
   EXPECT_DOUBLE_EQ(.089,returnPair.second);

   returnPair = tallyData.get_data(1,2);
   EXPECT_DOUBLE_EQ(7.8,returnPair.first);
   EXPECT_DOUBLE_EQ(0.078,returnPair.second);

}
// end of MCNP5/dagmc/test/test_TallyData.cpp
