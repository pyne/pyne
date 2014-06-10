// MCNP5/dagmc/test/test_CellTally.cpp

#include "gtest/gtest.h"

#include "moab/CartVect.hpp"

#include "../CellTally.hpp"
#include "../TallyEvent.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
// Base test fixture that defines common input parameters
class CellTallyTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // Create default CellTally object
        input.tally_id = 1;
        input.tally_type = "cell_coll";
        input.energy_bin_bounds.push_back(0.0);
        input.energy_bin_bounds.push_back(10.0);
        input.multiplier_id = -1;  // indicates not using multipliers
        cell_tally1 = new CellTally(input, TallyEvent::NONE);

        // Create a collision-based CellTally object
        input.tally_id = 2;
        std::multimap<std::string, std::string> options;
        options.insert(std::make_pair("cell", "10"));
        options.insert(std::make_pair("volume", "0.01"));
        input.options = options;
        cell_tally2 = new CellTally(input, TallyEvent::COLLISION); 

        // Create a tracklength-based CellTally object
        input.tally_id = 0;
        options.insert(std::make_pair("cell", "45"));
        options.insert(std::make_pair("volume", "2.36"));
        input.options = options;
	input.multiplier_id = 2;  // Use a multiplier for this tally
        cell_tally3 = new CellTally(input, TallyEvent::TRACK); 
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete cell_tally1;
        delete cell_tally2;
        delete cell_tally3;
    }

  protected:
    // data needed for each test
    TallyInput input;
    CellTally* cell_tally1;
    CellTally* cell_tally2;
    CellTally* cell_tally3;
};
//---------------------------------------------------------------------------//
// SIMPLE TESTS
//---------------------------------------------------------------------------//
// Test parsing of an invalid cell id
TEST(CellTallyInputTest, InvalidCellID)
{
    TallyInput input;

    // add general input data
    input.tally_id = 1;
    input.energy_bin_bounds.push_back(0.0);
    input.energy_bin_bounds.push_back(10.0);

    std::multimap<std::string, std::string> options;
    options.insert(std::make_pair("cell", "hello10"));
    input.options = options;

    CellTally* cell_tally;
    EXPECT_NO_THROW(cell_tally = new CellTally(input, TallyEvent::NONE));
    EXPECT_EQ(1, cell_tally->get_cell_id());
}
//---------------------------------------------------------------------------//
// Test parsing of an invalid volume
TEST(CellTallyInputTest, InvalidCellVolume)
{
    TallyInput input;

    // add general input data
    input.tally_id = 1;
    input.energy_bin_bounds.push_back(0.0);
    input.energy_bin_bounds.push_back(10.0);

    std::multimap<std::string, std::string> options;
    options.insert(std::make_pair("volume", "Omaha! Omaha!"));
    input.options = options;

    CellTally* cell_tally;
    EXPECT_NO_THROW(cell_tally = new CellTally(input, TallyEvent::NONE));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: CellTallyTest
//---------------------------------------------------------------------------//
// Test that cell_id mismatch results in no score being computed
TEST_F(CellTallyTest, NotInCell)
{
    TallyEvent event;
    event.current_cell = -1;

    // cell_tally1 has a default cell id of 1, since it was not set
    EXPECT_NO_THROW(cell_tally1->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally1->end_history());

    const TallyData& data1 = cell_tally1->getTallyData();
    std::pair<double, double> result = data1.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
    EXPECT_EQ(1, cell_tally1->get_cell_id());
    
    // Collision-based tally
    EXPECT_NO_THROW(cell_tally2->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally2->end_history());

    const TallyData& data2 = cell_tally2->getTallyData();
    result = data2.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
    EXPECT_EQ(10, cell_tally2->get_cell_id());
    
    // Track-based tally
    EXPECT_NO_THROW(cell_tally3->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
    EXPECT_EQ(45, cell_tally3->get_cell_id());
}
//---------------------------------------------------------------------------//
// Test the case of a match between the current_cell and the cell_id
TEST_F(CellTallyTest, InCell)
{
    TallyEvent event;
    event.current_cell = 1;  // default
    event.type         = TallyEvent::NONE;
    event.particle_energy     = 5.3;
    event.particle_weight     = 1.0;

    // No tally type defined
    EXPECT_NO_THROW(cell_tally1->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally1->end_history());

    const TallyData& data1 = cell_tally1->getTallyData();
    std::pair<double, double> result = data1.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
    EXPECT_EQ(1, cell_tally1->get_cell_id());

    // Collision-based tally
    event.type         = TallyEvent::COLLISION;
    event.current_cell = 10;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.total_cross_section = 8.1;

    EXPECT_NO_THROW(cell_tally2->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally2->end_history());

    const TallyData& data2 = cell_tally2->getTallyData();
    result = data2.get_data(0,0);
   
    EXPECT_NEAR(0.123457, result.first, 1e-6);
    EXPECT_NEAR(0.015242, result.second, 1e-6); 
    EXPECT_EQ(10, cell_tally2->get_cell_id());
    
    // Track-based tally
    event.type         = TallyEvent::TRACK;
    event.current_cell = 45;
    event.track_length = 7.9;
    event.direction = moab::CartVect(0.0, 0.0, 1.0);
    EXPECT_NO_THROW(cell_tally3->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(7.9, result.first);
    EXPECT_DOUBLE_EQ(62.41, result.second);
    EXPECT_EQ(45, cell_tally3->get_cell_id());
}
//---------------------------------------------------------------------------//
// Test the response to different event types
//---------------------------------------------------------------------------//
// Test TallyEvent::NONE
TEST_F(CellTallyTest,NullEventScore) 
{
    TallyEvent event;
    event.type             = TallyEvent::NONE;
    event.particle_weight  = 1.0;

    // No tally type defined
    event.current_cell = 1;  
    EXPECT_NO_THROW(cell_tally1->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally1->end_history());

    const TallyData& data1 = cell_tally1->getTallyData();
    std::pair<double, double> result = data1.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);

    // Collision-based tally
    event.current_cell = 10;
    EXPECT_NO_THROW(cell_tally2->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally2->end_history());

    const TallyData& data2 = cell_tally2->getTallyData();
    result = data2.get_data(0,0);
   
    EXPECT_NEAR(0.0, result.first,  1e-6);
    EXPECT_NEAR(0.0, result.second, 1e-6); 

    // Track-based tally
    event.current_cell = 45;
    EXPECT_NO_THROW(cell_tally3->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0,  result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
}
//---------------------------------------------------------------------------//
// Test TallyEvent::COLLISION
TEST_F(CellTallyTest,CollisionEventScore) 
{

    TallyEvent event;
    event.type             = TallyEvent::COLLISION;
    event.particle_weight  = 1.0;
    event.particle_energy  = 5.3;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.total_cross_section = 6.5;

    // No tally type defined
    event.current_cell = 1;  
    EXPECT_NO_THROW(cell_tally1->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally1->end_history());

    const TallyData& data1 = cell_tally1->getTallyData();
    std::pair<double, double> result = data1.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);

    // Collision-based tally
    event.current_cell = 10;
    EXPECT_NO_THROW(cell_tally2->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally2->end_history());

    const TallyData& data2 = cell_tally2->getTallyData();
    result = data2.get_data(0,0);
   
    EXPECT_NEAR(0.153846, result.first,  1e-6);
    EXPECT_NEAR(0.023669, result.second, 1e-6); 

    // Track-based tally
    event.current_cell = 45;
    EXPECT_NO_THROW(cell_tally3->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0,  result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);
}
//---------------------------------------------------------------------------//
// Test TallyEvent::TRACK
TEST_F(CellTallyTest,TrackEventScore) 
{
    TallyEvent event;
    event.type             = TallyEvent::TRACK;
    event.particle_weight  = 1.0;
    event.particle_energy  = 5.3;
    event.position  = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.0, 1.0, 0.0);
    event.track_length     = 2.8;

    // No tally type defined
    event.current_cell = 1;  
    EXPECT_NO_THROW(cell_tally1->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally1->end_history());

    const TallyData& data1 = cell_tally1->getTallyData();
    std::pair<double, double> result = data1.get_data(0,0);
    EXPECT_DOUBLE_EQ(0.0, result.first);
    EXPECT_DOUBLE_EQ(0.0, result.second);

    // Collision-based tally
    event.current_cell = 10;
    EXPECT_NO_THROW(cell_tally2->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally2->end_history());

    const TallyData& data2 = cell_tally2->getTallyData();
    result = data2.get_data(0,0);
   
    EXPECT_NEAR(0.0, result.first,  1e-6);
    EXPECT_NEAR(0.0, result.second, 1e-6); 

    // Track-based tally
    event.current_cell = 45;
    EXPECT_NO_THROW(cell_tally3->compute_score(event)); 
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(2.8,  result.first);
    EXPECT_DOUBLE_EQ(7.84, result.second);
} 
//---------------------------------------------------------------------------//
// Tests tally multiplier works
TEST_F(CellTallyTest, TallyMultiplier)
{
    TallyEvent event;
    event.type             = TallyEvent::TRACK;
    event.particle_weight  = 1.1;
    event.particle_energy  = 5.3;
    event.position  = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.0, 1.0, 0.0);
    event.track_length     = 2.8;
    event.current_cell = 45;

    // add some tally multipliers
    event.multipliers.push_back(2.4);
    event.multipliers.push_back(10.0);
    event.multipliers.push_back(12.9);

    // Using cell_tally3, which has tally multiplier id = 2
    EXPECT_NO_THROW(cell_tally3->compute_score(event));
    EXPECT_NO_THROW(cell_tally3->end_history());

    const TallyData& data3 = cell_tally3->getTallyData();
    std::pair<double, double> result = data3.get_data(0,0);
    EXPECT_DOUBLE_EQ(39.732, result.first);
    EXPECT_DOUBLE_EQ(1578.631824, result.second);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_CellTally.cpp
