// MCNP5/dagmc/test/test_KDENeighborhood.cpp

#include <cassert>
#include <cmath>
#include <set>

#include "gtest/gtest.h"

#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Types.hpp"

#include "../KDENeighborhood.hpp"
#include "../TallyEvent.hpp"

const double PI = 3.14159265359;

//---------------------------------------------------------------------------//
// HELPER METHODS
//---------------------------------------------------------------------------//
// loads default mesh into the mbi instance and gets all mesh nodes
void load_default_mesh(moab::Interface* mbi, moab::Range& mesh_nodes)
{
    moab::ErrorCode rval = mbi->load_mesh("../structured_mesh.h5m");
    assert(rval == moab::MB_SUCCESS);

    // copy all mesh nodes into Range
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);
}
//---------------------------------------------------------------------------//
// check all points are valid for given neighborhood region
bool check_all_points(const KDENeighborhood& region,
                      const std::set<moab::EntityHandle>& points)
{
    std::set<moab::EntityHandle>::iterator it;

    for (it = points.begin(); it != points.end(); ++it)
    {
        if (!region.is_calculation_point(*it)) return false;
    }

    return true;
}
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

        // load the default mesh and get all mesh nodes
        moab::Range mesh_nodes;
        load_default_mesh(mbi, mesh_nodes);

        // create neighborhood regions with and without kd-trees
        region1 = new KDENeighborhood(mbi, mesh_nodes, false);
        region2 = new KDENeighborhood(mbi, mesh_nodes, true);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete mbi;
        delete region1;
        delete region2;
    }

  protected:
    // data needed for each test
    moab::Interface* mbi;
    KDENeighborhood* region1;
    KDENeighborhood* region2;
};
//---------------------------------------------------------------------------//
class IsCalculationPointTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        mbi = new moab::Core();

        // load the default mesh
        rval = mbi->load_mesh("../structured_mesh.h5m");
        assert(rval == moab::MB_SUCCESS);

        // set up the TallyEvent and bandwidth
        event.type = TallyEvent::COLLISION;
        event.position = moab::CartVect(0.0, 0.0, 0.0);
        bandwidth = moab::CartVect(0.1, 0.2, 0.3);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete mbi;
    }

  protected:
    // data needed for each test
    moab::Interface* mbi;
    moab::CartVect bandwidth;
    moab::ErrorCode rval;
    TallyEvent event;
};
//---------------------------------------------------------------------------//
// SIMPLE TESTS
//---------------------------------------------------------------------------//
// Tests construction with no MOAB instance
TEST(KDENeighborhoodDeathTest, NullMOABInstance)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // define a NULL MOAB instance
    moab::Interface* mbi = NULL;
    moab::Range mesh_nodes;

    // check no errors occur if the kd-tree is set to off
    EXPECT_NO_THROW(KDENeighborhood(mbi, mesh_nodes, false));

    // check program exits if kd-tree is set on with NULL MOAB instance
    EXPECT_EXIT(KDENeighborhood(mbi, mesh_nodes, true),
                ::testing::ExitedWithCode(EXIT_FAILURE),
                "\nError: invalid moab::Interface for building KD-tree");
}
//---------------------------------------------------------------------------//
// Tests NULL tally event passed to the update_neighborhood method
TEST(KDENeighborhoodDeathTest, InvalidTallyEvent)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    moab::Core mb_core;
    moab::Interface* mbi = &mb_core;

    // load the default mesh and get all mesh nodes
    moab::Range mesh_nodes;
    load_default_mesh(mbi, mesh_nodes);

    // set up a NULL tally event
    moab::CartVect bandwidth(0.1, 0.1, 0.1);
    TallyEvent event;
    event.type = TallyEvent::NONE;
    
    // create neighborhood regions with and without kd-trees
    KDENeighborhood region1(mbi, mesh_nodes, false);
    KDENeighborhood region2(mbi, mesh_nodes, true);

    // check default number of points
    EXPECT_EQ(2025, region1.get_points().size());
    EXPECT_EQ(0, region2.get_points().size());

    // check NULL tally event behavior
    EXPECT_NO_THROW(region1.update_neighborhood(event, bandwidth));
    EXPECT_EXIT(region2.update_neighborhood(event, bandwidth),
                ::testing::ExitedWithCode(EXIT_FAILURE),
                 "\nError: Could not define neighborhood for tally event");
}
//---------------------------------------------------------------------------//
// Tests construction with an empty moab::Range containing no mesh nodes
TEST(KDENeighborhoodTest, NoMeshNodes)
{
    moab::Core mb_core;
    moab::Interface* mbi = &mb_core;

    // load default mesh
    moab::ErrorCode rval = mbi->load_mesh("../structured_mesh.h5m");
    assert(rval == moab::MB_SUCCESS);

    // define a Range containing no mesh nodes
    moab::Range mesh_nodes;

    // make sure no errors occur on construction with or without kd-tree
    EXPECT_NO_THROW(KDENeighborhood(mbi, mesh_nodes, false));
    EXPECT_NO_THROW(KDENeighborhood(mbi, mesh_nodes, true));

    // check default number of points in region1
    KDENeighborhood region1(mbi, mesh_nodes, false);
    std::set<moab::EntityHandle> points1 = region1.get_points();
    EXPECT_EQ(0, points1.size());

    // check default number of points in region2
    KDENeighborhood region2(mbi, mesh_nodes, true);
    std::set<moab::EntityHandle> points2 = region2.get_points();
    EXPECT_EQ(0, points2.size());
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
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 and region2 return all points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(2025, region2->get_points().size());
    
    // change to a conformal neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.position[0] = 2.0;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 1.0;
    bandwidth[0] = 2.0;
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 and region2 return all points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(2025, region2->get_points().size());
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
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points and region2 returns no points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(0, region2->get_points().size());

    // change to neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 1.0;
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points and region2 returns no points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(0, region2->get_points().size());
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
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points and region2 returns no points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(0, region2->get_points().size());

    // change to neighborhood based on track event
    event.type = TallyEvent::TRACK;
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 0.1;
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points and region2 returns no points
    EXPECT_EQ(2025, region1->get_points().size());
    EXPECT_EQ(0, region2->get_points().size());
}
//---------------------------------------------------------------------------//
// Tests correct points are returned for normal cases
TEST_F(GetPointsTest, GetPointsInBox)
{
    // define neighborhood using a track event (region overlaps z-mesh)
    TallyEvent event;
    double uvw_val = 1.0/sqrt(2.0);
    event.type = TallyEvent::TRACK;
    event.position = moab::CartVect(0.2, -0.2, 0.2);
    event.direction = moab::CartVect(uvw_val, 0.0, -1.0 * uvw_val);
    event.track_length = 2.3;
    moab::CartVect bandwidth(0.2, 0.2, 0.2);
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points
    EXPECT_EQ(2025, region1->get_points().size());

    // test number of points returned by region2 and check all are valid
    std::set<moab::EntityHandle> points1 = region2->get_points();
    EXPECT_EQ(320, points1.size());
    EXPECT_TRUE(check_all_points(*region2, points1));

    // change to neighborhood based on collision event (region inside mesh)
    event.type = TallyEvent::COLLISION;
    region1->update_neighborhood(event, bandwidth);
    region2->update_neighborhood(event, bandwidth);

    // test region1 still returns all points
    EXPECT_EQ(2025, region1->get_points().size());

    // test number of points returned by region2 and check all are valid
    std::set<moab::EntityHandle> points2 = region2->get_points();
    EXPECT_EQ(32, points2.size());
    EXPECT_TRUE(check_all_points(*region2, points2));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: IsCalculationPointTest
//---------------------------------------------------------------------------//
// Tests all points are calculation points when no kd-tree is used
TEST_F(IsCalculationPointTest, NoKDTree)
{
    // copy all mesh nodes into Range and convert into set
    moab::Range mesh_nodes;
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);

    std::set<moab::EntityHandle> mesh_set(mesh_nodes.begin(), mesh_nodes.end());
    EXPECT_EQ(2025, mesh_set.size());

    // create neighborhood with no kd-tree and check default calculation points
    KDENeighborhood region(mbi, mesh_nodes, false);

    EXPECT_EQ(2025, region.get_points().size());
    EXPECT_TRUE(check_all_points(region, mesh_set));

    // update neighborhood region and check it still includes all points
    region.update_neighborhood(event, bandwidth);

    EXPECT_EQ(2025, region.get_points().size());
    EXPECT_TRUE(check_all_points(region, mesh_set));
}
//---------------------------------------------------------------------------//
// Tests corners of neighborhood are valid calculation points
TEST_F(IsCalculationPointTest, ValidCornerPoints)
{
    // define additional valid calculation points at corners of neighborhood
    moab::Range corner_points;
    double corner_coords[] =
        { -0.1, -0.2, -0.3,
           0.1,  0.2,  0.3,
          -0.1, -0.2,  0.3,
          -0.1,  0.2, -0.3,
           0.1, -0.2,  0.3,
          -0.1,  0.2,  0.3,
           0.1, -0.2,  0.3,
           0.1,  0.2, -0.3 };

    // add new corner points to the MOAB instance and convert Range to set
    rval = mbi->create_vertices(corner_coords, 8, corner_points);
    assert(rval == moab::MB_SUCCESS);

    std::set<moab::EntityHandle> corner_set(corner_points.begin(),
                                            corner_points.end());

    EXPECT_EQ(8, corner_set.size());

    // copy all mesh nodes into Range
    moab::Range mesh_nodes;
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);

    EXPECT_EQ(2033, mesh_nodes.size());

    // create neighborhood with kd-tree and check default calculation points
    KDENeighborhood region(mbi, mesh_nodes, true);
    EXPECT_EQ(0, region.get_points().size());

    // update neighborhood region and check it includes corner set
    region.update_neighborhood(event, bandwidth);

    EXPECT_TRUE(region.get_points().size() > 0);
    EXPECT_TRUE(check_all_points(region, corner_set));
}
//---------------------------------------------------------------------------//
// Tests points along edges of neighborhood are valid calculation points
TEST_F(IsCalculationPointTest, ValidEdgePoints)
{
    // define additional valid calculation points along edges of neighborhood
    moab::Range edge_points;
    double edge_coords[] =
        { -0.1,  -0.2,   0.1667,
          -0.1,  -0.05, -0.3,
           0.0,  -0.2,  -0.3,
           0.1,   0.2,  -0.01,
           0.1,   0.13,  0.3,
          -0.09,  0.2,   0.3,
           0.07, -0.2,   0.3,
           0.1,  0.0,   -0.3 };

    // add new edge points to the MOAB instance and convert Range to set
    rval = mbi->create_vertices(edge_coords, 8, edge_points);
    assert(rval == moab::MB_SUCCESS);

    std::set<moab::EntityHandle> edge_set(edge_points.begin(),
                                          edge_points.end());

    EXPECT_EQ(8, edge_set.size());

    // copy all mesh nodes into Range
    moab::Range mesh_nodes;
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);

    EXPECT_EQ(2033, mesh_nodes.size());

    // create neighborhood with kd-tree and check default calculation points
    KDENeighborhood region(mbi, mesh_nodes, true);
    EXPECT_EQ(0, region.get_points().size());

    // update neighborhood region and check it includes edge set
    region.update_neighborhood(event, bandwidth);

    EXPECT_TRUE(region.get_points().size() > 0);
    EXPECT_TRUE(check_all_points(region, edge_set));
}
//---------------------------------------------------------------------------//
// Tests points inside neighborhood are valid calculation points
TEST_F(IsCalculationPointTest, ValidInteriorPoints)
{
    // define additional valid points within interior of neighborhood
    moab::Range interior_points;
    double interior_coords[] =
        {  0.0,    0.0,   0.0,
           0.001,  0.19,  0.2,
          -0.06,  -0.1,  -0.1,
          -0.0,    0.05, -0.23 };

    // add new interior points to the MOAB instance and convert Range to set
    rval = mbi->create_vertices(interior_coords, 4, interior_points);
    assert(rval == moab::MB_SUCCESS);

    std::set<moab::EntityHandle> interior_set(interior_points.begin(),
                                              interior_points.end());

    EXPECT_EQ(4, interior_set.size());

    // copy all mesh nodes into Range
    moab::Range mesh_nodes;
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);

    EXPECT_EQ(2029, mesh_nodes.size());

    // create neighborhood with kd-tree and check default calculation points
    KDENeighborhood region(mbi, mesh_nodes, true);
    EXPECT_EQ(0, region.get_points().size());

    // update neighborhood region and check it includes edge set
    region.update_neighborhood(event, bandwidth);

    EXPECT_TRUE(region.get_points().size() > 0);
    EXPECT_TRUE(check_all_points(region, interior_set));
}
//---------------------------------------------------------------------------//
// Tests points that are NOT in neighborhood are NOT valid calculation points
TEST_F(IsCalculationPointTest, InvalidPoints)
{
    // define additional points that are invalid calculation points
    moab::Range invalid_points;
    double invalid_coords[] =
        { -0.1, -0.2,  0.5,
          -0.1, -0.3, -0.3,
           0.2, -0.2, -0.3,
           0.1,  0.2, -0.7,
           0.1,  0.9,  0.3, 
          -0.1,  5.1, -2.9,
           1.7,  0.0, 10.4,
          -0.5, -0.8,  0.3,
           0.6, -0.9,  1.5,
          -0.7,  0.4, -0.5,
          -0.1, -0.3, -0.3 };
 
    // add new invalid points to the MOAB instance and convert Range to set
    rval = mbi->create_vertices(invalid_coords, 11, invalid_points);
    assert(rval == moab::MB_SUCCESS);

    std::set<moab::EntityHandle> invalid_set(invalid_points.begin(),
                                             invalid_points.end());

    EXPECT_EQ(11, invalid_set.size());

    // copy all mesh nodes into Range
    moab::Range mesh_nodes;
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);

    EXPECT_EQ(2036, mesh_nodes.size());

    // create neighborhood with kd-tree and check default calculation points
    KDENeighborhood region(mbi, mesh_nodes, true);
    EXPECT_EQ(0, region.get_points().size());

    // update neighborhood region and check invalid points are all still invalid
    region.update_neighborhood(event, bandwidth);

    std::set<moab::EntityHandle>::iterator it;
    for(it = invalid_set.begin(); it != invalid_set.end(); ++it)
    {
        EXPECT_FALSE(region.is_calculation_point(*it));
    }
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDENeighborhood.cpp
