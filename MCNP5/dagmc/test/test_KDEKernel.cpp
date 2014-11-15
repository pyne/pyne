// MCNP5/dagmc/test/test_KDEKernel.cpp

#include "gtest/gtest.h"
#include "../KDEKernel.hpp"

//---------------------------------------------------------------------------//
// MOCK OBJECTS
//---------------------------------------------------------------------------//
class MockEpanechnikovKernel : public KDEKernel
{
  public:
    // Constructor
    MockEpanechnikovKernel(){}

    // >>> DERIVED PUBLIC INTEFACE from KDEKernel.hpp

    // evaluates the second-order Epanechnikov kernel K(u)
    double evaluate(double u) const
    {
        // test if outside domain u = [-1, 1]
        if (u < -1.0 || u > 1.0) return 0.0;

        // compute K(u)
        return 0.75 * (1 - u * u);
    }

    // not implemented
    std::string get_kernel_name() const {}
    int get_order() const {}
    int get_min_quadrature(unsigned int i) const {}

    // integrates the ith moment function
    double integrate_moment(double a, double b, unsigned int i) const
    {
        // create the ith moment function
        MomentFunction moment(i, *this);
        
        // define the quadrature set for integrating the ith moment function
        unsigned int n = 2 + (i/2);
        Quadrature quadrature(n);

        // evaluate the integral
        return quadrature.integrate(a, b, moment);
    }
};
//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
// Tests KDEKernel factory method
class KDEKernelTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        kernel = NULL;
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kernel;
    }

  protected:
    // data needed for each test
    KDEKernel* kernel;
};
//---------------------------------------------------------------------------//
// Tests 1D boundary kernel implementation
class BoundaryKernel1DTest : public KDEKernelTest
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        kernel = new MockEpanechnikovKernel();
        side = 0;
        p_ratio = 1;
    }

  protected:
    // data needed for each test
    unsigned int side;
    double p_ratio;
};
//---------------------------------------------------------------------------//
// Tests both 2D and 3D boundary kernel implementation.
class BoundaryKernel3DTest : public KDEKernelTest
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        kernel = new MockEpanechnikovKernel();

        // set up parameters that will produce valid moments
        u.push_back(0.1);
        u.push_back(-0.5);
        u.push_back(1.0);

        p.push_back(0.3);
        p.push_back(0.0);
        p.push_back(0.6);

        sides.push_back(1);
        sides.push_back(0);
        sides.push_back(1);
    }

  protected:
    // data needed for each test
    std::vector<double> u;
    std::vector<double> p;
    std::vector<unsigned int> sides;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: KDEKernelTest
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateInvalidKernelType)
{
    kernel = KDEKernel::createKernel("invalid_kernel");
    EXPECT_TRUE(kernel == NULL);
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateInvalidHigherOrderKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov", 5);
    EXPECT_TRUE(kernel == NULL);
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateInvalidZeroOrderKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov", 0);
    EXPECT_TRUE(kernel == NULL);
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateInvalidTypeAndOrderKernel)
{
    kernel = KDEKernel::createKernel("invalid_kernel", 5);
    EXPECT_TRUE(kernel == NULL);
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateUniformKernel)
{
    kernel = KDEKernel::createKernel("uniform");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order uniform", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());
    EXPECT_EQ(1, kernel->get_min_quadrature(0));
    EXPECT_EQ(1, kernel->get_min_quadrature(1));
    EXPECT_EQ(2, kernel->get_min_quadrature(2));
    EXPECT_EQ(2, kernel->get_min_quadrature(3));
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateEpanechnikovKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order epanechnikov", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());
    EXPECT_EQ(2, kernel->get_min_quadrature(0));
    EXPECT_EQ(2, kernel->get_min_quadrature(1));
    EXPECT_EQ(3, kernel->get_min_quadrature(2));
    EXPECT_EQ(3, kernel->get_min_quadrature(3));
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateBiweightKernel)
{
    kernel = KDEKernel::createKernel("biweight");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order biweight", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());
    EXPECT_EQ(3, kernel->get_min_quadrature(0));
    EXPECT_EQ(3, kernel->get_min_quadrature(1));
    EXPECT_EQ(4, kernel->get_min_quadrature(2));
    EXPECT_EQ(4, kernel->get_min_quadrature(3));
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateTriweightKernel)
{
    kernel = KDEKernel::createKernel("triweight");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order triweight", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());
    EXPECT_EQ(4, kernel->get_min_quadrature(0));
    EXPECT_EQ(4, kernel->get_min_quadrature(1));
    EXPECT_EQ(5, kernel->get_min_quadrature(2));
    EXPECT_EQ(5, kernel->get_min_quadrature(3));
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateHigherOrderKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov", 6);
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("6th-order epanechnikov", kernel->get_kernel_name());
    EXPECT_EQ(6, kernel->get_order());
    EXPECT_EQ(4, kernel->get_min_quadrature(0));
    EXPECT_EQ(4, kernel->get_min_quadrature(1));
    EXPECT_EQ(5, kernel->get_min_quadrature(2));
    EXPECT_EQ(5, kernel->get_min_quadrature(3));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: BoundaryKernel1DTest
//---------------------------------------------------------------------------//
// Tests point located at the max distance from the LOWER boundary
// i.e. distance = bandwidth, side = 0 (LOWER)
TEST_F(BoundaryKernel1DTest, EvaluatePointAtLowerMax)
{
    // test correction factor is one over the domain u = [-1, 1]
    double u1[] = {-1.0, -0.8, -0.6, -0.4, -0.2,
                    0.0,  0.2,  0.4,  0.6,  0.8, 1.0};

    for (int i = 0; i < 11; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(1.0, value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [-1, 1]
    double u2[] = {-2.0, 2.0};

    double value2 = kernel->boundary_correction(&u2[0], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);

    value2 = kernel->boundary_correction(&u2[1], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);
}
//---------------------------------------------------------------------------//
// Tests point located at the max distance from the UPPER boundary
// i.e. distance = bandwidth, side = 1 (UPPER)
TEST_F(BoundaryKernel1DTest, EvaluatePointAtUpperMax)
{
    side = 1;

    // test correction factor is one over the domain u = [-1, 1]
    double u1[] = {-1.0, -0.8, -0.6, -0.4, -0.2,
                    0.0,  0.2,  0.4,  0.6,  0.8, 1.0};

    for (int i = 0; i < 11; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(1.0, value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [-1, 1]
    double u2[] = {-2.0, 2.0};

    double value2 = kernel->boundary_correction(&u2[0], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);

    value2 = kernel->boundary_correction(&u2[1], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);
}
//---------------------------------------------------------------------------//
// Tests point located at half the max distance from the LOWER boundary
// i.e. distance = 0.5 * bandwidth, side = 0 (LOWER)
TEST_F(BoundaryKernel1DTest, EvaluatePointAtHalfLowerMax)
{
    p_ratio = 0.5;

    // define valid u and reference values
    double u1[] = {-1.0, -0.8, -0.6, -0.4,
                   -0.2,  0.0,  0.2,  0.5};

    double ref[] = {0.220500, 0.440999, 0.661499, 0.881998,
                    1.102498, 1.322997, 1.543497, 1.874246};

    // test correction factor over the domain u = [-1, 0.5]
    for (int i = 0; i < 8; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(ref[i], value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [-1, 0.5]
    double u2[] = {-2.0, 0.8, 1.0, 2.0};

    for (int i = 0; i < 4; ++i)
    {
        double value2 = kernel->boundary_correction(&u2[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(0.0, value2);
    }
}
//---------------------------------------------------------------------------//
// Tests point located at half the max distance from the UPPER boundary
// i.e. distance = 0.5 * bandwidth, side = 1 (UPPER)
TEST_F(BoundaryKernel1DTest, EvaluatePointAtHalfUpperMax)
{
    p_ratio = 0.5;
    side = 1;

    // define valid u and reference values
    double u1[] = {-0.5, -0.2, 0.0, 0.2,
                    0.4,  0.6, 0.8, 1.0};

    double ref[] = {1.874246, 1.543497, 1.322997, 1.102498,
                    0.881998, 0.661499, 0.440999, 0.220500};

    // test correction factor over the domain u = [-0.5, 1]
    for (int i = 0; i < 8; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(ref[i], value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [-0.5, 1]
    double u2[] = {-2.0, -1.0, -0.8, 2.0};

    for (int i = 0; i < 4; ++i)
    {
        double value2 = kernel->boundary_correction(&u2[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(0.0, value2);
    }
}
//---------------------------------------------------------------------------//
// Tests point located on the LOWER boundary
// i.e. distance = 0, side = 0 (LOWER)
TEST_F(BoundaryKernel1DTest, EvaluatePointOnLowerBoundary)
{
    p_ratio = 0.0;

    // define valid u and reference values
    double u1[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0};

    double ref[] = {-5.894737, -3.368421, -0.842105,
                     1.684211,  4.210526,  6.736842};

    // test correction factor over the domain u = [-1, 0]
    for (int i = 0; i < 6; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(ref[i], value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [-1, 0]
    double u2[] = {-2.0, 0.5, 1.0, 2.0};

    for (int i = 0; i < 4; ++i)
    {
        double value2 = kernel->boundary_correction(&u2[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(0.0, value2);
    }
}
//---------------------------------------------------------------------------//
// Tests point located on the UPPER boundary
// i.e. distance = 0, side = 1 (UPPER)
TEST_F(BoundaryKernel1DTest, EvaluatePointOnUpperBoundary)
{
    p_ratio = 0.0;
    side = 1;

    // define valid u and reference values
    double u1[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    double ref[] = { 6.736842,  4.210526,  1.684211,
                    -0.842105, -3.368421, -5.894737};

    // test correction factor over the domain u = [0, 1]
    for (int i = 0; i < 6; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_NEAR(ref[i], value1, 1e-6);
    }

    // test correction factor is zero outside the domain u = [0, 1]
    double u2[] = {-2.0, -0.5, -1.0, 2.0};

    for (int i = 0; i < 4; ++i)
    {
        double value2 = kernel->boundary_correction(&u2[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(0.0, value2);
    }
}
//---------------------------------------------------------------------------//
// Tests point that is outside the max distance from a boundary
TEST_F(BoundaryKernel1DTest, EvaluatePointOutsideMaxDistance)
{
    p_ratio = 1.5;

    // test correction factor is always one over the domain u = [-1, 1]
    double u1[] = {-1.0, -0.8, -0.6, -0.4, -0.2,
                    0.0,  0.2,  0.4,  0.6,  0.8, 1.0};

    for (int i = 0; i < 11; ++i)
    {
        double value1 = kernel->boundary_correction(&u1[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(1.0, value1);
    }

    // test correction factor is zero outside the domain u = [-1, 1]
    double u2[] = {-2.0, 2.0};

    double value2 = kernel->boundary_correction(&u2[0], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);

    value2 = kernel->boundary_correction(&u2[1], &p_ratio, &side, 1);
    EXPECT_DOUBLE_EQ(0.0, value2);
}
//---------------------------------------------------------------------------//
// Tests point that is at a negative distance from a boundary
TEST_F(BoundaryKernel1DTest, EvaluateNegativeDistanceRatio)
{
    p_ratio = -1.5;

    // test correction factor is always zero if p ratio is negative
    double u2[] = {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0};

    for (int i = 0; i < 7; ++i)
    {
        double value2 = kernel->boundary_correction(&u2[i], &p_ratio, &side, 1);
        EXPECT_DOUBLE_EQ(0.0, value2);
    }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: BoundaryKernel3DTest
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, AllValidKernelMoments)
{
    // 1D correction
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_NEAR(1.275992, value, 1e-6);

    // 3D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_NEAR(-0.183359, value, 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, ReduceCorrectionMatrixTo2D)
{
    // change 3D correction to include one point outside max distance
    p[2] = 2.3;

    // 1D correction
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_NEAR(1.275992, value, 1e-6);

    // 3D correction reduces to the 2D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_NEAR(1.275992, value, 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, ReduceCorrectionMatrixTo1D)
{
    // change 2D and 3D corrections to include points outside max distance
    p[1] = 1.1;
    p[2] = 7.9;

    // 1D correction
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D correction reduces to the 1D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 3D correction reduces to the 1D correction
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_NEAR(1.737159, value, 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, AllInvalidKernelMoments)
{
    // change u so that all moments will be invalid
    u[0] = -0.6;
    u[1] = 0.5;
    u[2] = -1.0;

    // test correction factor is always zero
    for (int i = 1; i <= 3; ++i)
    {
        double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], i);
        EXPECT_DOUBLE_EQ(0.0, value);
    }
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForU)
{
    // change p so that only the first set of moments will be invalid
    p[0] = -1.5;

    // test correction factor is always zero
    for (int i = 1; i <= 3; ++i)
    {
        double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], i);
        EXPECT_DOUBLE_EQ(0.0, value);
    }
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForV)
{
    // change u so that only the second set of moments will be invalid
    u[1] = 0.5;

    // 1D correction should still work
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D and 3D correction factors should be zero
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForW)
{
    // change p so that only the third set of moments will be invalid
    p[2] = -5;

    // 1D correction should still work
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D correction should still work
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_NEAR(1.275992, value, 1e-6);

    // 3D correction factor should be zero
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForUV)
{
    // change p so that the first and second set of moments will be invalid
    p[0] = -1.0;
    p[1] = 2.0;

    // test correction factor is always zero
    for (int i = 1; i <= 3; ++i)
    {
        double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], i);
        EXPECT_DOUBLE_EQ(0.0, value);
    }
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForUW)
{
    // change u so that the first and third set of moments will be invalid
    u[0] = 1.2;
    u[2] = -0.7;

    // test correction factor is always zero
    for (int i = 1; i <= 3; ++i)
    {
        double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], i);
        EXPECT_DOUBLE_EQ(0.0, value);
    }
}
//---------------------------------------------------------------------------//
TEST_F(BoundaryKernel3DTest, InvalidKernelMomentsForVW)
{
    // change u/p so that the second and third set of moments will be invalid
    p[1] = -1.2;
    u[2] = 10.5;

    // 1D correction should still work
    double value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 1);
    EXPECT_NEAR(1.737159, value, 1e-6);

    // 2D and 3D correction factors should be zero
    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 2);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->boundary_correction(&u[0], &p[0], &sides[0], 3);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDEKernel.cpp
