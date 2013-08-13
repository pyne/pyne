// MCNP5/dagmc/test/test_KDEMeshTally.cpp

#include <cmath>

#include "gtest/gtest.h"

#include "../KDEMeshTally.hpp"
#include "../TallyEvent.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
// Base test fixture that defines common input parameters
class KDEMeshTallyTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // define general tally options
        std::multimap<std::string, std::string> options;
        options.insert(std::make_pair("hx", "0.1"));
        options.insert(std::make_pair("hy", "0.1"));
        options.insert(std::make_pair("hz", "0.1"));

        // add general input data
        input.tally_id = 1;
        input.input_filename = "../structured_mesh.h5m";
        input.energy_bin_bounds.push_back(0.0);
        input.energy_bin_bounds.push_back(10.0);
        input.total_energy_bin = false;
        input.options = options;
    }

  protected:
    // data needed for each test
    MeshTallyInput input;
    KDEMeshTally* kde_tally;

    // wrapper for the KDEMeshTally::integral_track_score method
    double test_integral_track_score(const moab::CartVect& coords,
                                     const TallyEvent& event)
    {
        KDEMeshTally::CalculationPoint X;
        X.coords[0] = coords[0];
        X.coords[1] = coords[1];
        X.coords[2] = coords[2];
        return kde_tally->integral_track_score(X, event);
    }

    // wrapper for the KDEMeshTally::subtrack_score method
    double test_subtrack_score(const moab::CartVect& coords,
                               const std::vector<moab::CartVect>& points)
    {
        KDEMeshTally::CalculationPoint X;
        X.coords[0] = coords[0];
        X.coords[1] = coords[1];
        X.coords[2] = coords[2];
        return kde_tally->subtrack_score(X, points);
    }

    // accessor method to change the bandwidth value
    void change_bandwidth(const moab::CartVect& new_bandwidth)
    {
        kde_tally->bandwidth = new_bandwidth;
    }
};
//---------------------------------------------------------------------------//
// Tests the private integral_track_score method in KDEMeshTally
class KDEIntegralTrackTest : public KDEMeshTallyTest
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // set up default input parameters for kde mesh tally
        KDEMeshTallyTest::SetUp();

        // create kde mesh tally
        kde_tally = new KDEMeshTally(input, KDEMeshTally::INTEGRAL_TRACK);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kde_tally;
    }
};
//---------------------------------------------------------------------------//
// Tests the private subtrack_score method in KDEMeshTally
class KDESubtrackTest : public KDEMeshTallyTest
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // set up default input parameters for kde mesh tally
        KDEMeshTallyTest::SetUp();

        // create kde mesh tally
        kde_tally = new KDEMeshTally(input, KDEMeshTally::SUB_TRACK);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kde_tally;
    }

  protected:
    // data needed for each test
    std::vector<moab::CartVect> points;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: KDEIntegralTrackTest
//---------------------------------------------------------------------------//
// Tests cases that have a valid [Smin, Smax] interval with Smin != Smax
TEST_F(KDEIntegralTrackTest, ValidLimits)
{
    // set up tally event
    TallyEvent event;
    event.type = TallyEvent::TRACK;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    double U_val = 1.0/sqrt(3.0);
    event.direction = moab::CartVect(U_val, U_val, U_val);
    event.track_length = 1.0;

    // Case 1: Smin, Smax < 0
    moab::CartVect coords1(-0.215, -0.157, -0.1);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords1, event));

    // Case 2: Smin < 0, Smax = 0
    event.direction = moab::CartVect(0.1, -sqrt(0.18), 0.9);
    event.track_length = 0.123;
    moab::CartVect coords2(0.06, 0.1, -0.08);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords2, event));

    // Case 3: Smin < 0, 0 < Smax < dic
    event.direction = moab::CartVect(0.9, -sqrt(0.18), 0.1);
    event.track_length = 0.2;
    moab::CartVect coords3(0.01, -0.015147186258, 0.085);
    EXPECT_NEAR(12.445817, test_integral_track_score(coords3, event), 1e-6);

    // Case 4: Smin = 0, 0 < Smax < dic
    event.position = moab::CartVect(1.0, 2.15, -1.45);
    event.direction = moab::CartVect(1.0, 0.0, 0.0);
    event.track_length = 1.0;
    moab::CartVect coords4(1.0, 2.2, -1.5);
    EXPECT_NEAR(15.820313, test_integral_track_score(coords4, event), 1e-6);

    // Case 5: 0 < Smin, Smax < dic
    event.position = moab::CartVect(0.14, 0.0, -0.28787753827);
    event.direction = moab::CartVect(-0.2, 0.0, sqrt(0.96));
    change_bandwidth(moab::CartVect(0.1, 0.2, 0.3));
    moab::CartVect coords5(0.0, 0.0, 0.0);
    EXPECT_NEAR(10.300734, test_integral_track_score(coords5, event), 1e-6);

    // Case 6: 0 < Smin < dic, Smax = dic
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.0, sqrt(0.5), -sqrt(0.5));
    change_bandwidth(moab::CartVect(0.1, 0.1, 0.1));
    moab::CartVect coords6(0.0, 0.67781745931, -0.67781745931);
    EXPECT_NEAR(48.320769, test_integral_track_score(coords6, event), 1e-6);

    // Case 7: 0 < Smin < dic, Smax > dic
    event.position = moab::CartVect(-0.2, -0.4, -0.1);
    event.direction = moab::CartVect(-0.7, 0.2, -sqrt(0.47));
    moab::CartVect coords7(-0.93, -0.14, -0.817008914036);
    EXPECT_NEAR(8.444942, test_integral_track_score(coords7, event), 1e-6);

    // Case 8: Smin = dic, Smax > dic
    event.position = moab::CartVect(-0.5, 0.5, -0.5);
    event.direction = moab::CartVect(0.5, 0.1, sqrt(0.74));
    event.track_length = 0.5;
    change_bandwidth(moab::CartVect(0.5, 0.5, 0.5));
    moab::CartVect coords8(-0.5, 1.05, 0.344093010682);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords8, event));

    // Case 9: Smin, Smax > dic
    event.position = moab::CartVect(-1.0, -2.0, -3.0);
    event.direction = moab::CartVect(-0.5, 0.5, sqrt(0.5));
    event.track_length = 0.725;
    change_bandwidth(moab::CartVect(0.1, 0.1, 0.1));
    moab::CartVect coords9(-1.85, -1.2, -1.83933982822);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords9, event));

    // Case 10: Smin = 0, Smax = dic
    event.position = moab::CartVect(0.5, 0.4, 0.2);
    event.direction = moab::CartVect(0.0, 0.0, -1.0);
    event.track_length = 1.0;
    change_bandwidth(moab::CartVect(0.2, 0.4, 0.5));
    moab::CartVect coords10(0.6, 0.33, -0.3);
    EXPECT_NEAR(5.111938, test_integral_track_score(coords10, event), 1e-6); 

    // Case 11: Smin < 0, Smax > dic
    event.position = moab::CartVect(-1.0, -1.0, -1.0);
    event.direction = moab::CartVect(0.8, sqrt(0.2), 0.4);
    change_bandwidth(moab::CartVect(1.0, 1.0, 1.0));
    moab::CartVect coords11(-0.4, -0.2683281573, -0.4);
    EXPECT_NEAR(0.240948, test_integral_track_score(coords11, event), 1e-6);

    // Case 12: Smin < 0, Smax = dic
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(-sqrt(0.5), 0.5, 0.5);
    moab::CartVect coords12(-0.292893218813, -0.5, 0.5);
    EXPECT_NEAR(0.150014, test_integral_track_score(coords12, event), 1e-6);

    // Case 13: Smin = 0, Smax > dic
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(-sqrt(0.63), -0.1, -0.6);
    event.track_length = 1.3;
    moab::CartVect coords13(-1.0, 0.0, -1.0);
    EXPECT_NEAR(0.257392, test_integral_track_score(coords13, event), 1e-6);
}
//---------------------------------------------------------------------------//
// Tests cases that do not have a valid [Smin, Smax] interval
TEST_F(KDEIntegralTrackTest, InvalidLimits)
{
    // set up tally event
    TallyEvent event;
    event.type = TallyEvent::TRACK;
    event.position = moab::CartVect(0.0, 0.0, 0.0);
    event.direction = moab::CartVect(0.5, sqrt(0.5), 0.5);
    event.track_length = 1.0;

    // Case 14: No overlaps for Sx, Sy, Sz
    moab::CartVect coords1(0.1, 0.453553390593, -0.15);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords1, event));

    // Case 15: Only Sx-Sy overlaps
    moab::CartVect coords2(0.0, 0.1, 0.3);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords2, event));

    // Case 16: Only Sx-Sy and Sy-Sz overlap
    moab::CartVect coords3(0.0, 0.1, 0.225);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords3, event));

    // Case 17: Sx & Sz are identical, but Sy only overlaps at boundary
    moab::CartVect coords4(0.1, -0.1, 0.1);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords4, event));

    // Case 18: Sx, Sy, Sz only overlap at boundaries
    moab::CartVect coords5(0.0, 0.241421356237, -0.2);
    EXPECT_DOUBLE_EQ(0.0, test_integral_track_score(coords5, event));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: KDESubtrackTest
//---------------------------------------------------------------------------//
TEST_F(KDESubtrackTest, NoSubtracks)
{
    // verify no score is returned for a few calculation points
    moab::CartVect coords1(0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords1, points));

    moab::CartVect coords2(-1.5, -0.7, -0.18);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords2, points));

    moab::CartVect coords3(5.13, -9.27, 0.0);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords3, points));
}
//---------------------------------------------------------------------------//
TEST_F(KDESubtrackTest, OneSubtrack)
{
    // add only one subtrack point
    points.push_back(moab::CartVect(-1.0, 0.0, 2.0));

    // verify no score is returned for calculation points outside neighborhood
    moab::CartVect coords1(-2.0, 5.7, -4.2);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords1, points));

    moab::CartVect coords2(-2.0, 0.0, -4.2);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords2, points));

    moab::CartVect coords3(-0.9, 0.0, 3.0);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords3, points));

    // verify score returned for calculation points inside neighborhood
    moab::CartVect coords4(-1.0, 0.0, 2.0);
    EXPECT_DOUBLE_EQ(421.875, test_subtrack_score(coords4, points));

    moab::CartVect coords5(-0.9, 0.1, 1.9);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords5, points));

    moab::CartVect coords6(-1.033, -0.01, 1.99);
    EXPECT_NEAR(368.451750, test_subtrack_score(coords6, points), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(KDESubtrackTest, MultipleSubtracks)
{
    // add multiple subtrack points
    points.push_back(moab::CartVect(-0.05, -0.05, -0.05));
    points.push_back(moab::CartVect(0.0, 0.0, 0.0));
    points.push_back(moab::CartVect(0.1, 0.1, 0.1));

    // verify no score is returned for calculation points outside neighborhood
    moab::CartVect coords1(-1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords1, points));

    moab::CartVect coords2(0.8, 0.0, 4.9);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords2, points));

    moab::CartVect coords3(0.0, -5.2, 0.1);
    EXPECT_DOUBLE_EQ(0.0, test_subtrack_score(coords3, points));

    // verify score returned for calculation points inside neighborhood
    moab::CartVect coords4(0.0, 0.0, 0.0);
    EXPECT_NEAR(199.951172, test_subtrack_score(coords4, points), 1e-6);

    moab::CartVect coords5(-0.02, 0.01, 0.03);
    EXPECT_NEAR(151.105500, test_subtrack_score(coords5, points), 1e-6);

    moab::CartVect coords6(0.01, 0.02, 0.03);
    EXPECT_NEAR(143.051063, test_subtrack_score(coords6, points), 1e-6);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDEMeshTally.cpp
