// MCNP5/dagmc/test/test_TallyEvent.cpp

#include <vector>

#include "gtest/gtest.h"

#include "moab/CartVect.hpp"

#include "../TallyEvent.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class TallyEventTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        score_multiplier = 0.0;

        // set up event 1 with no multipliers
        event1 = new TallyEvent();
        event1->particle_weight = 1.0;

        // set up event 2 with only one multiplier
        std::vector<double> multipliers;
        multipliers.push_back(10);

        event2 = new TallyEvent();
        event2->particle_weight = 4.68;
        event2->multipliers = multipliers;

        // set up event 3 with five multipliers
        multipliers.push_back(0.0);
        multipliers.push_back(1.0);
        multipliers.push_back(-3.1);
        multipliers.push_back(7.9);

        event3 = new TallyEvent();
        event3->particle_weight = 0.53;
        event3->multipliers = multipliers;
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete event1;
        delete event2;
        delete event3;
    }

  protected:
    // data needed for each test
    double score_multiplier;
    TallyEvent* event1;
    TallyEvent* event2;
    TallyEvent* event3;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: TallyEventTest
//---------------------------------------------------------------------------//
// Tests only particle weight is returned for invalid multiplier index
TEST_F(TallyEventTest, InvalidMultiplierIndex)
{
    // test case with no multipliers
    EXPECT_NO_THROW(score_multiplier = event1->get_score_multiplier(-3));
    EXPECT_DOUBLE_EQ(1.0, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event1->get_score_multiplier(0));
    EXPECT_DOUBLE_EQ(1.0, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event1->get_score_multiplier(4));
    EXPECT_DOUBLE_EQ(1.0, score_multiplier);

    // test case with one multiplier
    EXPECT_NO_THROW(score_multiplier = event2->get_score_multiplier(-3));
    EXPECT_DOUBLE_EQ(4.68, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event2->get_score_multiplier(1));
    EXPECT_DOUBLE_EQ(4.68, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event2->get_score_multiplier(4));
    EXPECT_DOUBLE_EQ(4.68, score_multiplier);

    // test case with five multipliers
    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(-3));
    EXPECT_DOUBLE_EQ(0.53, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(5));
    EXPECT_DOUBLE_EQ(0.53, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(8));
    EXPECT_DOUBLE_EQ(0.53, score_multiplier);
}
//---------------------------------------------------------------------------//
// Tests correct multiplier returned for valid multiplier index
TEST_F(TallyEventTest, ValidMultiplierIndex)
{
    // test case with no multipliers
    EXPECT_NO_THROW(score_multiplier = event1->get_score_multiplier(-1));
    EXPECT_DOUBLE_EQ(1.0, score_multiplier);

    // test case with one multiplier
    EXPECT_NO_THROW(score_multiplier = event2->get_score_multiplier(-1));
    EXPECT_DOUBLE_EQ(4.68, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event2->get_score_multiplier(0));
    EXPECT_DOUBLE_EQ(46.8, score_multiplier);

    // test case with five multipliers
    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(-1));
    EXPECT_DOUBLE_EQ(0.53, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(0));
    EXPECT_DOUBLE_EQ(5.3, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(1));
    EXPECT_DOUBLE_EQ(0.0, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(2));
    EXPECT_DOUBLE_EQ(0.53, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(3));
    EXPECT_DOUBLE_EQ(-1.643, score_multiplier);

    EXPECT_NO_THROW(score_multiplier = event3->get_score_multiplier(4));
    EXPECT_DOUBLE_EQ(4.187, score_multiplier);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_TallyEvent.cpp
