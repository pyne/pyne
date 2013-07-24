// MCNP5/dagmc/test/test_KDENeighborhood.cpp

#include <cmath>

#include "gtest/gtest.h"

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Types.hpp"

#include "../KDENeighborhood.hpp"
#include "../TallyEvent.hpp"

const double PI = 3.14159265359;

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class PointWithinMaxRadiusTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // define bandwidth and maximum radius
        moab::CartVect bandwidth(0.1, 0.1, 0.1);
        max_radius = bandwidth.length();

        // define fake kd-tree variables for constructor
        moab::AdaptiveKDTree* tree = NULL;
        moab::EntityHandle tree_root = 0;

        // create first neighborhood region
        event1.type = TallyEvent::TRACK;
        event1.position = moab::CartVect(0.0, 0.0, 0.0);
        event1.direction = moab::CartVect(1.0, 0.0, 0.0);
        event1.track_length = 1.0;
        region1 = new KDENeighborhood(event1, bandwidth, *tree, tree_root);

        // create second neighborhood region
        event2.type = TallyEvent::TRACK;
        event2.position = moab::CartVect(-0.1, 0.1, -0.3);
        event2.direction = moab::CartVect(-0.5, -0.3, sqrt(0.66));
        event2.track_length = 2.5;
        region2 = new KDENeighborhood(event2, bandwidth, *tree, tree_root);

        // create third neighborhood region
        event3.type = TallyEvent::COLLISION;
        event3.position = moab::CartVect(0.0, 0.0, 0.0);
        event3.total_cross_section = 0.1;
        region3 = new KDENeighborhood(event3, bandwidth, *tree, tree_root);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete region1;
        delete region2;
        delete region3;
    }

  protected:
    // data needed for each test
    double max_radius;
    TallyEvent event1;
    TallyEvent event2;
    TallyEvent event3;
    KDENeighborhood* region1;
    KDENeighborhood* region2;
    KDENeighborhood* region3;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: PointWithinMaxRadiusTest
//---------------------------------------------------------------------------//
TEST_F(PointWithinMaxRadiusTest, PointOutsideMaxRadius)
{
    // test point outside max radius in region1
    moab::CartVect point1(0.7, -0.02, 2 * max_radius);
    EXPECT_FALSE(region1->point_within_max_radius(point1));

    // test point outside max radius in region2
    moab::CartVect point2(-0.7, 2 * max_radius, 0.67588);
    EXPECT_FALSE(region2->point_within_max_radius(point2));

    // test point outside max radius in region3 (max radius should be zero)
    moab::CartVect point3(0.06, -0.03, 0.02);
    EXPECT_FALSE(region3->point_within_max_radius(point3));
}
//---------------------------------------------------------------------------//
TEST_F(PointWithinMaxRadiusTest, PointAtMaxRadius)
{
    // test point at max radius in region1
    moab::CartVect point1(-0.1, max_radius * cos(PI/4), max_radius * sin(PI/4));
    EXPECT_FALSE(region1->point_within_max_radius(point1));

    // test point at max radius in region2
    moab::CartVect point2(-0.1 + 0.514495755428 * max_radius,
                           0.1 - 0.857492925713 * max_radius,
                          -0.3);
    EXPECT_FALSE(region2->point_within_max_radius(point2));

    // test point at max radius in region3 (max radius should be zero)
    moab::CartVect point3(0.0, 0.0, 0.0);
    EXPECT_FALSE(region3->point_within_max_radius(point3));
}
//---------------------------------------------------------------------------//
TEST_F(PointWithinMaxRadiusTest, PointOnTrack)
{
    // test point on track segment in region1
    moab::CartVect point1(0.0, 0.0, 0.0);
    EXPECT_TRUE(region1->point_within_max_radius(point1));

    // test point on track segment in region2
    moab::CartVect point2(-1.35, -0.65, 1.73100960116);
    EXPECT_TRUE(region2->point_within_max_radius(point2));

    // no test point in region 3 can be on a track as it is a collision
}
//---------------------------------------------------------------------------//
TEST_F(PointWithinMaxRadiusTest, PointInsideMaxRadius)
{
    // test point inside max radius in region1
    moab::CartVect point1(1.05, 0.09, -0.04);
    EXPECT_TRUE(region1->point_within_max_radius(point1));

    // test point inside max radius in region2
    moab::CartVect point2(-0.65, -0.1, 0.59);
    EXPECT_TRUE(region2->point_within_max_radius(point2));

    // no test point in region 3 can be inside max radius
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDENeighborhood.cpp
