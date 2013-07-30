// MCNP5/dagmc/test/test_KDENeighborhood.cpp

#include <cassert>
#include <cmath>
#include <set>

#include "gtest/gtest.h"

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Types.hpp"

#include "../KDENeighborhood.hpp"
#include "../TallyEvent.hpp"

const double PI = 3.14159265359;

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class GetPointsTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        mbi = new moab::Core();

        // load the mesh data
        moab::ErrorCode rval = mbi->load_mesh("../structured_mesh.h5m");
        assert(rval == moab::MB_SUCCESS);

        // get all of the mesh nodes from the MOAB root set
        moab::EntityHandle root_set = 0;
        moab::Range mesh_nodes;
        rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
        assert(rval == moab::MB_SUCCESS);

        // build a kd-tree from all of the mesh nodes
        kd_tree = new moab::AdaptiveKDTree(mbi);
        rval = kd_tree->build_tree(mesh_nodes, kd_tree_root);
        assert(rval == moab::MB_SUCCESS);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kd_tree;
        delete mbi;
    }

  protected:
    // data needed for each test
    moab::Interface* mbi;
    moab::AdaptiveKDTree* kd_tree;
    moab::EntityHandle kd_tree_root;

    // helper function to check all points are valid
    bool check_all_points(const std::set<moab::EntityHandle>& points,
                          const KDENeighborhood& region)
    {
        std::set<moab::EntityHandle>::iterator it;

        for (it = points.begin(); it != points.end(); ++it)
        {
            moab::CartVect coords(0.0, 0.0, 0.0);
            moab::EntityHandle point = *it;
            moab::ErrorCode rval = mbi->get_coords(&point, 1, coords.array());
            assert(rval == moab::MB_SUCCESS);

            if (!region.point_in_region(coords)) return false;
        }

        return true;
    }
};
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
// SIMPLE TESTS
//---------------------------------------------------------------------------//
TEST(KDENeighborhoodDeathTest, InvalidTallyEvent)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // define a null tally event and fake kd-tree
    TallyEvent event;
    event.type = TallyEvent::NONE;
    moab::CartVect bandwidth(0.1, 0.1, 0.1);
    moab::AdaptiveKDTree* kd_tree = NULL;
    moab::EntityHandle kd_tree_root = 0;

    // make sure KDENeighborhood returns error for NULL tally event
    EXPECT_EXIT(KDENeighborhood(event, bandwidth, *kd_tree, kd_tree_root),
                ::testing::ExitedWithCode(EXIT_FAILURE),
                "\nError: Could not define neighborhood for tally event");
}
//---------------------------------------------------------------------------//
TEST(KDENeighborhoodTest, ValidPointsForCollision)
{
    // define fake kd-tree variables
    moab::AdaptiveKDTree* kd_tree = NULL;
    moab::EntityHandle kd_tree_root = 0;

    // define a neighborhood region using a collision event
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    moab::CartVect bandwidth(0.1, 0.2, 0.3);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // check all corner points are valid
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, -0.2, -0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, 0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, -0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, 0.2, -0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, -0.2, -0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, 0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, -0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, 0.2, -0.3)));
    
    // check some edge points are valid
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, -0.2, 0.1667)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.1, -0.05, -0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.0, -0.2, -0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, 0.2, -0.01)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, 0.13, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.09, 0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.07, -0.2, 0.3)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.1, 0.0, -0.3)));

    // check some interior points are valid
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.0, 0.0, 0.0)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(0.001, 0.19, 0.2)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.06, -0.1, -0.1)));
    EXPECT_TRUE(region1.point_in_region(moab::CartVect(-0.0, 0.05, -0.23)));
}
//---------------------------------------------------------------------------//
TEST(KDENeighborhoodTest, InvalidPointsForCollision)
{
    // define fake kd-tree variables
    moab::AdaptiveKDTree* kd_tree = NULL;
    moab::EntityHandle kd_tree_root = 0;

    // define a neighborhood region using a collision event
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    moab::CartVect bandwidth(0.1, 0.2, 0.3);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // check some points with two valid coordinates
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.1, -0.2, 0.5)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.1, -0.3, -0.3)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(0.2, -0.2, -0.3)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(0.1, 0.2, -0.7)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(0.1, -1.8, 0.3)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.1, 0.9, 0.3)));

    // check some points with one valid coordinate
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.1, 5.1, -2.9)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(1.7, 0.0, 10.4)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.5, -0.8, 0.3)));

    // check some points with no valid coordinates
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(0.6, -0.9, 1.5)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.7, 0.4, -0.5)));
    EXPECT_FALSE(region1.point_in_region(moab::CartVect(-0.1, -0.3, -0.3)));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: GetPointsTest
//---------------------------------------------------------------------------//
// Tests all points are returned if neighborhood is conformal to mesh
TEST_F(GetPointsTest, GetAllPointsInBox)
{
    // define neighborhood conformal to mesh based on collision event
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(2.5, 0.0, 0.0);
    moab::CartVect bandwidth(2.5, 0.5, 0.5);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // test if all points are returned and check all are valid
    std::set<moab::EntityHandle> points1 = region1.get_points();
    EXPECT_EQ(2025, points1.size());
    EXPECT_TRUE(check_all_points(points1, region1));

    // change to a conformal neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.position[0] = 2.0;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 1.0;
    bandwidth[0] = 2.0;
    KDENeighborhood region2(event, bandwidth, *kd_tree, kd_tree_root);

    // test if all points are returned and check all are valid
    std::set<moab::EntityHandle> points2 = region2.get_points();
    EXPECT_EQ(2025, points2.size());
    EXPECT_TRUE(check_all_points(points2, region2));
}
//---------------------------------------------------------------------------//
// Tests no points are returned if neighborhood exists outside mesh
TEST_F(GetPointsTest, GetNoPointsOutsideBox)
{
    // define neighborhood outside the mesh based on collision event
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(-5.0, 0.0, 0.0);
    moab::CartVect bandwidth(1.0, 0.5, 0.5);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // test if no points are returned
    std::set<moab::EntityHandle> points1 = region1.get_points();
    EXPECT_EQ(0, points1.size());

    // change to neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 1.0;
    KDENeighborhood region2(event, bandwidth, *kd_tree, kd_tree_root);

    // test if no points are returned
    std::set<moab::EntityHandle> points2 = region2.get_points();
    EXPECT_EQ(0, points2.size());
}
//---------------------------------------------------------------------------//
// Tests no points are returned if mesh cell is bigger than neighborhood
TEST_F(GetPointsTest, GetNoPointsInBox)
{
    // define neighborhood within mesh cell based on collision event
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(2.6, -0.06, 0.06);
    moab::CartVect bandwidth(0.05, 0.05, 0.05);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // test if no points are returned
    std::set<moab::EntityHandle> points1 = region1.get_points();
    EXPECT_EQ(0, points1.size());

    // change to neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 0.1;
    KDENeighborhood region2(event, bandwidth, *kd_tree, kd_tree_root);

    // test if no points are returned
    std::set<moab::EntityHandle> points2 = region2.get_points();
    EXPECT_EQ(0, points2.size());
}
//---------------------------------------------------------------------------//
TEST_F(GetPointsTest, GetPointsInBox)
{
    // define neighborhood using a collision event (region inside mesh)
    TallyEvent event;
    event.type = TallyEvent::COLLISION;
    event.position = moab::CartVect(0.2, -0.2, 0.2);
    moab::CartVect bandwidth(0.2, 0.2, 0.2);
    KDENeighborhood region1(event, bandwidth, *kd_tree, kd_tree_root);

    // test number of points returned and check all are valid
    std::set<moab::EntityHandle> points1 = region1.get_points();
    EXPECT_EQ(32, points1.size());
    EXPECT_TRUE(check_all_points(points1, region1));

    // change to neighborhood based on track event (region overlaps z-mesh)
    double uvw_val = 1.0/sqrt(2.0);
    event.type = TallyEvent::TRACK;
    event.direction = moab::CartVect(uvw_val, 0.0, -1.0 * uvw_val);
    event.track_length = 2.3;
    KDENeighborhood region2(event, bandwidth, *kd_tree, kd_tree_root);

    // test number of points returned and check all are valid
    std::set<moab::EntityHandle> points2 = region2.get_points();
    EXPECT_EQ(320, points2.size());
    EXPECT_TRUE(check_all_points(points2, region2));
}
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
