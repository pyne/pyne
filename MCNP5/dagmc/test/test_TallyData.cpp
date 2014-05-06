// MCNP5/dagmc/test/test_TallyData.cpp

#include "gtest/gtest.h"

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
        // 'Normal' case with single energy bin
        tallyData1 = new TallyData(1,false);
        
        // Case with multiple energy bins plus total energy bin
        tallyData2 = new TallyData(5, true);

        // Case with multiple energy bins without total energy bin
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
// Helper function to fill initial tally data
void fillTally(TallyData& tallyData)
{
    double fill_tally[]   = {1.2, -3.4, 5.6, 0.0, 8.9, 7.8};
    double fill_error[]   = {0.2, -0.4, 0.06, 0.0, .089, .078};
    double fill_scratch[] = {10.2, -30.4, 50.6, 0.0, 80.9, 70.8};

    // Populate tally data arrays for the tallyData object
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
   // Set up empty tally with 3 tally points and 2 energy bins
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
TEST_F(TallyDataTest, NullEndHistory)
{
      int length;

      // One tally point, one energy bin 
      tallyData1->resize_data_arrays(1); 
      // Three tally points, 5 energy bins, total  
      tallyData2->resize_data_arrays(3);

      // This is a null end history - no scoring has been done
      tallyData1->end_history();
      // test arrays initialized
      double* data    = tallyData1->get_tally_data(length);
      double* error   = tallyData1->get_error_data(length); 
      double* scratch = tallyData1->get_scratch_data(length); 

      for (int i=0; i<length; i++)
      {
        EXPECT_DOUBLE_EQ(0.0, data[i]);
        EXPECT_DOUBLE_EQ(0.0, error[i]);
        EXPECT_DOUBLE_EQ(0.0, scratch[i]);
      }

      tallyData2->end_history(); 
      data    = tallyData2->get_tally_data(length);
      error   = tallyData2->get_error_data(length); 
      scratch = tallyData2->get_scratch_data(length); 

      for (int i=0; i<length; i++)
      {
        EXPECT_DOUBLE_EQ(0.0, data[i]);
        EXPECT_DOUBLE_EQ(0.0, error[i]);
        EXPECT_DOUBLE_EQ(0.0, scratch[i]);
      }
}
//---------------------------------------------------------------------------//
TEST_F(TallyDataTest, EndHistory)
{
      int length;

      // One tally point, one energy bin 
      tallyData1->resize_data_arrays(1); 
      
      tallyData1->add_score_to_tally(0,5.6,0);
      double* tally_data = tallyData1->get_tally_data(length); 
      double* error_data = tallyData1->get_error_data(length);
      double* scratch_data = tallyData1->get_scratch_data(length);
      
      tallyData1->end_history();

      EXPECT_DOUBLE_EQ(5.6, tally_data[0]);
      EXPECT_DOUBLE_EQ(31.36, error_data[0]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[0]);
      //////////////////////////////////////////////////////////    
   
      // Two tally points, 5 energy bins, total  
      tallyData2->resize_data_arrays(2);
      double fill_tally[]   = {1.2, -3.4, 5.6, 0.0, 8.9, 12.3};

      // Set it up - note we don't include the total bin, that 
      // gets done as each score is tallied
      for (int i=0; i<5; i++)
      { 
         tallyData2->add_score_to_tally(0,fill_tally[i], i); 
      } 

      tallyData2->end_history();
      tally_data = tallyData2->get_tally_data(length); 
      error_data = tallyData2->get_error_data(length);
      scratch_data = tallyData2->get_scratch_data(length);

      // The first half of the data arrays will have data.
      double total_value = 0.0;
      for (int i=0; i<5; i++)
      { 
         double tally_value = fill_tally[i];
         double error_value = fill_tally[i]*fill_tally[i];
         total_value += fill_tally[i];
         EXPECT_DOUBLE_EQ(tally_value, tally_data[i]);
         EXPECT_DOUBLE_EQ(error_value, error_data[i]);
         EXPECT_DOUBLE_EQ(0.0, scratch_data[i]);
      }
 
      // Check the total energy_bin separately
      EXPECT_DOUBLE_EQ(fill_tally[5], total_value);

      double total_error = total_value * total_value;
      EXPECT_DOUBLE_EQ(total_value, tally_data[5]);
      EXPECT_DOUBLE_EQ(total_error, error_data[5]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[5]);
     
      // Check the second half of the data arrays
      for (int i=6; i<12; i++)
      { 
         EXPECT_DOUBLE_EQ(0.0, tally_data[i]);
         EXPECT_DOUBLE_EQ(0.0, error_data[i]);
         EXPECT_DOUBLE_EQ(0.0, scratch_data[i]);
      } 
}
//---------------------------------------------------------------------------//
TEST_F(TallyDataTest, AddScoreToTally)
{
      int length; 

      // One tally point, one energy bin 
      tallyData1->resize_data_arrays(1); 
      double* scratch_data = tallyData1->get_scratch_data(length);

      tallyData1->add_score_to_tally(0,5.6,0);
      tallyData1->add_score_to_tally(0,7.8,0);
      EXPECT_DOUBLE_EQ(13.4, scratch_data[0]);
      
      tallyData1->add_score_to_tally(0,0.0,0);
      EXPECT_DOUBLE_EQ(13.4, scratch_data[0]);

      tallyData1->add_score_to_tally(0,-2.1,0);
      EXPECT_DOUBLE_EQ(11.3, scratch_data[0]);
      //////////////////////////////////////////////////////////    
   
      // Two tally points, 5 energy bins, total  
      tallyData2->resize_data_arrays(2);
      scratch_data = tallyData2->get_scratch_data(length);

      tallyData2->add_score_to_tally(1,5.6,0);
      tallyData2->add_score_to_tally(1,7.8,0);
      EXPECT_DOUBLE_EQ(13.4, scratch_data[6]);
      EXPECT_DOUBLE_EQ(13.4, scratch_data[11]);
      
      tallyData2->add_score_to_tally(1,5.6,2);
      tallyData2->add_score_to_tally(1,7.8,2);
      EXPECT_DOUBLE_EQ(13.4, scratch_data[8]);
      EXPECT_DOUBLE_EQ(26.8, scratch_data[11]);
      
      tallyData2->add_score_to_tally(0,1.1,3);
      EXPECT_DOUBLE_EQ(1.1, scratch_data[3]);
      EXPECT_DOUBLE_EQ(1.1, scratch_data[5]);

      tallyData2->add_score_to_tally(0,-2.1,0);
      EXPECT_DOUBLE_EQ(-2.1, scratch_data[0]);
      EXPECT_DOUBLE_EQ(-1.0, scratch_data[5]);
 
      EXPECT_DOUBLE_EQ(0.0, scratch_data[1]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[2]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[4]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[7]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[9]);
      EXPECT_DOUBLE_EQ(0.0, scratch_data[10]);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_TallyData.cpp
