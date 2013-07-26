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
class BoundaryKernelTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        kernel = new MockEpanechnikovKernel();
        side = 0;
        distance = 0.5;
        max_distance = 0.5;
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kernel;
    }

  protected:
    // data needed for each test
    KDEKernel* kernel;
    unsigned int side;
    double distance;
    double max_distance;
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
TEST_F(KDEKernelTest, CreateUniformKernel)
{
    kernel = KDEKernel::createKernel("uniform");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order uniform", kernel->get_kernel_name());
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateEpanechnikovKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order epanechnikov", kernel->get_kernel_name());
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateBiweightKernel)
{
    kernel = KDEKernel::createKernel("biweight");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order biweight", kernel->get_kernel_name());
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateTriweightKernel)
{
    kernel = KDEKernel::createKernel("triweight");
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("2nd-order triweight", kernel->get_kernel_name());
}
//---------------------------------------------------------------------------//
TEST_F(KDEKernelTest, CreateHigherOrderKernel)
{
    kernel = KDEKernel::createKernel("epanechnikov", 6);
    EXPECT_TRUE(kernel != NULL);
    EXPECT_EQ("6th-order epanechnikov", kernel->get_kernel_name());
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: BoundaryKernelTest
//---------------------------------------------------------------------------//
// Tests point located at the max distance from the LOWER boundary
// i.e. distance = max_distance, side = 0 (LOWER)
TEST_F(BoundaryKernelTest, EvaluatePointAtLowerMax)
{
    // test evaluation over the domain u = [-1, 1]
    double value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    value = kernel->evaluate(-0.8, max_distance, distance, side);
    EXPECT_NEAR(0.27, value, 1e-6);

    value = kernel->evaluate(-0.6, max_distance, distance, side);
    EXPECT_NEAR(0.48, value, 1e-6);

    value = kernel->evaluate(-0.4, max_distance, distance, side);
    EXPECT_NEAR(0.63, value, 1e-6);

    value = kernel->evaluate(-0.2, max_distance, distance, side);
    EXPECT_NEAR(0.72, value, 1e-6);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(0.75, value, 1e-6);

    value = kernel->evaluate(0.2, max_distance, distance, side);
    EXPECT_NEAR(0.72, value, 1e-6);

    value = kernel->evaluate(0.4, max_distance, distance, side);
    EXPECT_NEAR(0.63, value, 1e-6);

    value = kernel->evaluate(0.6, max_distance, distance, side);
    EXPECT_NEAR(0.48, value, 1e-6);

    value = kernel->evaluate(0.8, max_distance, distance, side);
    EXPECT_NEAR(0.27, value, 1e-6);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    // test evaluation outside the domain u = [-1, 1]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point located at the max distance from the UPPER boundary
// i.e. distance = max_distance, side = 1 (UPPER)
TEST_F(BoundaryKernelTest, EvaluatePointAtUpperMax)
{
    side = 1;

    // test evaluation over the domain u = [-1, 1]
    double value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    value = kernel->evaluate(-0.8, max_distance, distance, side);
    EXPECT_NEAR(0.27, value, 1e-6);

    value = kernel->evaluate(-0.6, max_distance, distance, side);
    EXPECT_NEAR(0.48, value, 1e-6);

    value = kernel->evaluate(-0.4, max_distance, distance, side);
    EXPECT_NEAR(0.63, value, 1e-6);

    value = kernel->evaluate(-0.2, max_distance, distance, side);
    EXPECT_NEAR(0.72, value, 1e-6);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(0.75, value, 1e-6);

    value = kernel->evaluate(0.2, max_distance, distance, side);
    EXPECT_NEAR(0.72, value, 1e-6);

    value = kernel->evaluate(0.4, max_distance, distance, side);
    EXPECT_NEAR(0.63, value, 1e-6);

    value = kernel->evaluate(0.6, max_distance, distance, side);
    EXPECT_NEAR(0.48, value, 1e-6);

    value = kernel->evaluate(0.8, max_distance, distance, side);
    EXPECT_NEAR(0.27, value, 1e-6);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    // test evaluation outside the domain u = [-1, 1]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point located at half the max distance from the LOWER boundary
// i.e. distance = 0.5 * max_distance, side = 0 (LOWER)
TEST_F(BoundaryKernelTest, EvaluatePointAtHalfLowerMax)
{
    distance = 0.25;

    // test evaluation over the domain u = [-1, 0.5]
    double value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    value = kernel->evaluate(-0.8, max_distance, distance, side);
    EXPECT_NEAR(0.119070, value, 1e-6);

    value = kernel->evaluate(-0.6, max_distance, distance, side);
    EXPECT_NEAR(0.317519, value, 1e-6);

    value = kernel->evaluate(-0.4, max_distance, distance, side);
    EXPECT_NEAR(0.555659, value, 1e-6);

    value = kernel->evaluate(-0.2, max_distance, distance, side);
    EXPECT_NEAR(0.793798, value, 1e-6);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(0.992248, value, 1e-6);

    value = kernel->evaluate(0.2, max_distance, distance, side);
    EXPECT_NEAR(1.111318, value, 1e-6);

    value = kernel->evaluate(0.5, max_distance, distance, side);
    EXPECT_NEAR(1.054264, value, 1e-6);

    // test evaluation outside the domain u = [-1, 0.5]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.8, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point located at half the max distance from the UPPER boundary
// i.e. distance = 0.5 * max_distance, side = 1 (UPPER)
TEST_F(BoundaryKernelTest, EvaluatePointAtHalfUpperMax)
{
    distance = 0.25;
    side = 1;

    // test evaluation over the domain u = [-0.5, 1]
    double value = kernel->evaluate(-0.5, max_distance, distance, side);
    EXPECT_NEAR(1.054264, value, 1e-6);

    value = kernel->evaluate(-0.2, max_distance, distance, side);
    EXPECT_NEAR(1.111318, value, 1e-6);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(0.992248, value, 1e-6);

    value = kernel->evaluate(0.2, max_distance, distance, side);
    EXPECT_NEAR(0.793798, value, 1e-6);

    value = kernel->evaluate(0.4, max_distance, distance, side);
    EXPECT_NEAR(0.555659, value, 1e-6);

    value = kernel->evaluate(0.6, max_distance, distance, side);
    EXPECT_NEAR(0.317519, value, 1e-6);

    value = kernel->evaluate(0.8, max_distance, distance, side);
    EXPECT_NEAR(0.119070, value, 1e-6);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    // test evaluation outside the domain u = [-0.5, 1]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-0.8, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point located on the LOWER boundary
// i.e. distance = 0, side = 0 (LOWER)
TEST_F(BoundaryKernelTest, EvaluatePointOnLowerBoundary)
{
    distance = 0.0;

    // test evaluation over the domain u = [-1, 0]
    double value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    value = kernel->evaluate(-0.8, max_distance, distance, side);
    EXPECT_NEAR(-0.909474, value, 1e-6);

    value = kernel->evaluate(-0.6, max_distance, distance, side);
    EXPECT_NEAR(-0.404211, value, 1e-6);

    value = kernel->evaluate(-0.4, max_distance, distance, side);
    EXPECT_NEAR(1.061053, value, 1e-6);

    value = kernel->evaluate(-0.2, max_distance, distance, side);
    EXPECT_NEAR(3.031579, value, 1e-6);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(5.052632, value, 1e-6);

    // test evaluation outside the domain u = [-1, 0]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point located on the UPPER boundary
// i.e. distance = 0, side = 1 (UPPER)
TEST_F(BoundaryKernelTest, EvaluatePointOnUpperBoundary)
{
    distance = 0.0;
    side = 1;

    // test evaluation over the domain u = [0, 1]
    double value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_NEAR(5.052632, value, 1e-6);

    value = kernel->evaluate(0.2, max_distance, distance, side);
    EXPECT_NEAR(3.031579, value, 1e-6);

    value = kernel->evaluate(0.4, max_distance, distance, side);
    EXPECT_NEAR(1.061053, value, 1e-6);

    value = kernel->evaluate(0.6, max_distance, distance, side);
    EXPECT_NEAR(-0.404211, value, 1e-6);

    value = kernel->evaluate(0.8, max_distance, distance, side);
    EXPECT_NEAR(-0.909474, value, 1e-6);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_NEAR(0.0, value, 1e-6);

    // test evaluation outside the domain u = [0, 1]
    value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point that is outside the max distance from a boundary
TEST_F(BoundaryKernelTest, EvaluatePointOutsideMaxDistance)
{
    distance = 1.5;

    // test evaluation over various values to verify always outside domain
    double value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//
// Tests point that is at a negative distance from a boundary
TEST_F(BoundaryKernelTest, EvaluateNegativeDistanceRatio)
{
    distance = -1.5;

    // test evaluation over various values to verify always outside domain
    double value = kernel->evaluate(-2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(-0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(0.5, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(1.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);

    value = kernel->evaluate(2.0, max_distance, distance, side);
    EXPECT_DOUBLE_EQ(0.0, value);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDEKernel.cpp
