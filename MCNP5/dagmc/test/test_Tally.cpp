// MCNP5/dagmc/test/test_Tally.cpp

#include "gtest/gtest.h"

#include "moab/CartVect.hpp"

#include "../Tally.hpp"
#include "../TallyEvent.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class TallyFactoryTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        tally = NULL;
	// Minimum required member settings
        input.energy_bin_bounds.push_back(0.0);
        input.energy_bin_bounds.push_back(10.0);
	input.tally_id = 1;
	input.options.insert(std::make_pair("inp", "../unstructured_mesh.h5m"));
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete tally;
    }

  protected:
    // data needed for each test
    Tally* tally;
    TallyInput input;
};

//---------------------------------------------------------------------------//
// SIMPLE TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyTest
//---------------------------------------------------------------------------//
// Tests for the factory method
TEST_F(TallyFactoryTest, TallyConstructor1)
{
  input.tally_type = "cell_coll";
  tally = Tally::create_tally(input); 
  const TallyData& data = tally->getTallyData();
  EXPECT_EQ(1, data.get_num_energy_bins());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, TallyConstructor2)
{
  input.tally_type = "cell_coll";
  input.energy_bin_bounds.push_back(12.0);
  input.energy_bin_bounds.push_back(15.0);
  input.energy_bin_bounds.push_back(21.0);
  tally = Tally::create_tally(input); 
  const TallyData& data = tally->getTallyData();
  EXPECT_EQ(5, data.get_num_energy_bins());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateTrackLengthMeshTally)
{
  input.tally_type = "unstr_track";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("unstr_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//

TEST_F(TallyFactoryTest, CreateKDETrackMeshTally)
{
  input.tally_type = "kde_track";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateKDESubtrackMeshTally)
{
  input.tally_type = "kde_subtrack";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_subtrack", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateKDECollisionMeshTally)
{
  input.tally_type = "kde_coll";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_coll", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateCellTrackTally)
{
  input.tally_type = "cell_track";
  input.options.clear();
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("cell_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateCellCollisionTally)
{
  input.tally_type = "cell_coll";
  input.options.clear();
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("cell_coll", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, TallyTypeNotSet)
{
  input.options.clear(); 
  
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally == NULL);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// end of MCNP5/dagmc/test/test_Tally.cpp
