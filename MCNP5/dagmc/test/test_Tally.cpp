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
	// Minimum required member settings
        input.energy_bin_bounds.push_back(0.0);
        input.energy_bin_bounds.push_back(10.0);
	input.tally_id = 1;
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
class TallyEnergyBinTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        input.tally_id = 1;
        input.tally_type = "cell_track";

	// Set a default energy bin case
        input.energy_bin_bounds.push_back(0.0);
        input.energy_bin_bounds.push_back(10.0);
        input.energy_bin_bounds.push_back(20.3);
        input.energy_bin_bounds.push_back(34.5);
        input.energy_bin_bounds.push_back(67.9);

	event.type = TallyEvent::TRACK;
	event.current_cell = 1;
	event.track_length = 9.0;
	event.particle_weight = 1.0;
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
    TallyEvent event;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyFactoryTest
//---------------------------------------------------------------------------//
// Tests Tally constructor for default number of energy bins
TEST_F(TallyFactoryTest, DefaultNumEnergyBins)
{
  input.tally_type = "cell_coll";
  tally = Tally::create_tally(input);
  const TallyData& data = tally->getTallyData();
  EXPECT_EQ(1, data.get_num_energy_bins());
}
//---------------------------------------------------------------------------//
// Tests Tally constructor for non-default number of energy bins
TEST_F(TallyFactoryTest, NumEnergyBins)
{
  input.tally_type = "cell_coll";

  // increase number of energy bins from default
  input.energy_bin_bounds.push_back(12.0);
  input.energy_bin_bounds.push_back(15.0);
  input.energy_bin_bounds.push_back(21.0);

  tally = Tally::create_tally(input);
  const TallyData& data = tally->getTallyData();

  // should have 4 energy bins plus one total bin
  EXPECT_EQ(5, data.get_num_energy_bins());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateTrackLengthMeshTally)
{
  input.tally_type = "unstr_track";
  input.options.insert(std::make_pair("inp", "../unstructured_mesh.h5m"));
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("unstr_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//

TEST_F(TallyFactoryTest, CreateKDETrackMeshTally)
{
  input.tally_type = "kde_track";
  input.options.insert(std::make_pair("inp", "../unstructured_mesh.h5m"));
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateKDESubtrackMeshTally)
{
  input.tally_type = "kde_subtrack";
  input.options.insert(std::make_pair("inp", "../unstructured_mesh.h5m"));
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_subtrack", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateKDECollisionMeshTally)
{
  input.tally_type = "kde_coll";
  input.options.insert(std::make_pair("inp", "../unstructured_mesh.h5m"));
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("kde_coll", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateCellTrackTally)
{
  input.tally_type = "cell_track";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("cell_track", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, CreateCellCollisionTally)
{
  input.tally_type = "cell_coll";
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally != NULL);
  EXPECT_EQ("cell_coll", tally->get_tally_type());
}
//---------------------------------------------------------------------------//
TEST_F(TallyFactoryTest, TallyTypeNotSet)
{
  tally = Tally::create_tally(input);
  EXPECT_TRUE(tally == NULL);
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyEnergyBinTest
//---------------------------------------------------------------------------//
TEST_F(TallyEnergyBinTest, EnergyNotInBounds)
{
  tally = Tally::create_tally(input);
  event.particle_energy = -10.0;
  tally->compute_score(event);
  tally->end_history();

  const TallyData& data = tally->getTallyData();
  for (int i=0; i<5; ++i)
  {
     std::pair<double, double> result = data.get_data(0,i);
     EXPECT_DOUBLE_EQ(0.0, result.first);
     EXPECT_DOUBLE_EQ(0.0, result.second);
  }

  event.particle_energy = 69.9;
  tally->compute_score(event);
  tally->end_history();
  for (int i=0; i<5; ++i)
  {
     std::pair<double, double> result = data.get_data(0,i);
     EXPECT_DOUBLE_EQ(0.0, result.first);
     EXPECT_DOUBLE_EQ(0.0, result.second);
  }
}
//---------------------------------------------------------------------------//
//  Bin boundaries for 4 bins:
//	0.0 10.0  20.3  34.5  67.9
//	bin 0 includes 0
//	bin 1 includes 10.0
//	bin 2 includes 20.3
//	bin 3 includes 34.5 AND 67.9
//---------------------------------------------------------------------------//
TEST_F(TallyEnergyBinTest, EnergyInBounds)
{
  tally = Tally::create_tally(input);
  event.particle_energy = 0.0;
  tally->compute_score(event);
  tally->end_history();

  const TallyData& data = tally->getTallyData();
  std::pair<double, double> result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(9.0, result.first);
  EXPECT_DOUBLE_EQ(81.0, result.second);

  event.particle_energy = 10.0;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,1);
  EXPECT_DOUBLE_EQ(9.0, result.first);
  EXPECT_DOUBLE_EQ(81.0, result.second);

  event.particle_energy = 20.3;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,2);
  EXPECT_DOUBLE_EQ(9.0, result.first);
  EXPECT_DOUBLE_EQ(81.0, result.second);

  event.particle_energy = 25.87;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,2);
  EXPECT_DOUBLE_EQ(18.0, result.first);
  EXPECT_DOUBLE_EQ(162.0, result.second);

  event.particle_energy = 34.5;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,3);
  EXPECT_DOUBLE_EQ(9.0, result.first);
  EXPECT_DOUBLE_EQ(81.0, result.second);

  event.particle_energy = 67.9;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,3);
  EXPECT_DOUBLE_EQ(18.0, result.first);
  EXPECT_DOUBLE_EQ(162.0, result.second);

  // check total energy bin (i.e. bin 4)
  result = data.get_data(0,4);
  EXPECT_DOUBLE_EQ(54.0, result.first);
  EXPECT_DOUBLE_EQ(486.0, result.second);
}
//---------------------------------------------------------------------------//
// Tests single energy bin with (min_energy, max_energy) and min_energy != 0
TEST_F(TallyEnergyBinTest, NonZeroEnergyBounds)
{
  // reset energy_bin_bounds to include one energy bin
  input.energy_bin_bounds.clear();
  input.energy_bin_bounds.push_back(5.0);
  input.energy_bin_bounds.push_back(17.0);

  tally = Tally::create_tally(input);

  event.particle_energy = 3.3;
  tally->compute_score(event);
  tally->end_history();
  const TallyData& data = tally->getTallyData();
  std::pair<double, double> result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(0.0, result.first);
  EXPECT_DOUBLE_EQ(0.0, result.second);

  event.particle_energy = 5.0;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(9.0, result.first);
  EXPECT_DOUBLE_EQ(81.0, result.second);

  event.particle_energy = 8.2;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(18.0, result.first);
  EXPECT_DOUBLE_EQ(162.0, result.second);

  event.particle_energy = 17.0;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(27.0, result.first);
  EXPECT_DOUBLE_EQ(243.0, result.second);

  event.particle_energy = 18.0;
  tally->compute_score(event);
  tally->end_history();
  result = data.get_data(0,0);
  EXPECT_DOUBLE_EQ(27.0, result.first);
  EXPECT_DOUBLE_EQ(243.0, result.second);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_Tally.cpp
