// MCNP5/dagmc/test/test_PolynomialKernel.cpp

#include "gtest/gtest.h"
#include "../PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class PolynomialKernelTest : public ::testing::Test
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
    PolynomialKernel* kernel;
};
//---------------------------------------------------------------------------//
class IntegrateMomentTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        i = 0;
        a = -1.0;
        b = 1.0;
        kernel1 = new PolynomialKernel(0, 1);
        kernel2 = new PolynomialKernel(1, 1);
        kernel3 = new PolynomialKernel(2, 1);
        kernel4 = new PolynomialKernel(1, 2);
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete kernel1;
        delete kernel2;
        delete kernel3;
        delete kernel4;
    }

  protected:
    // data needed for each test
    int i;
    double a, b;
    PolynomialKernel* kernel1;
    PolynomialKernel* kernel2;
    PolynomialKernel* kernel3;
    PolynomialKernel* kernel4;
};
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: PolynomialKernelTest
//---------------------------------------------------------------------------//
// Tests the basic evaluate method of a 2nd-order uniform kernel
TEST_F(PolynomialKernelTest, EvaluateUniformKernel)
{
    kernel = new PolynomialKernel(0, 1);
    EXPECT_EQ("2nd-order uniform", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());

    // test evaluation over the domain
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(-1.0));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(-0.8));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(-0.6));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(-0.4));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(-0.2));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(0.0));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(0.2));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(0.4));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(0.6));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(0.8));
    EXPECT_DOUBLE_EQ(0.5, kernel->evaluate(1.0));

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate method of a 2nd-order epanechnikov kernel
TEST_F(PolynomialKernelTest, EvaluateEpanechnikovKernel)
{
    kernel = new PolynomialKernel(1, 1);
    EXPECT_EQ("2nd-order epanechnikov", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());

    // test evaluation over the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-1.0));
    EXPECT_DOUBLE_EQ(0.27, kernel->evaluate(-0.8));
    EXPECT_DOUBLE_EQ(0.48, kernel->evaluate(-0.6));
    EXPECT_DOUBLE_EQ(0.63, kernel->evaluate(-0.4));
    EXPECT_DOUBLE_EQ(0.72, kernel->evaluate(-0.2));
    EXPECT_DOUBLE_EQ(0.75, kernel->evaluate(0.0));
    EXPECT_DOUBLE_EQ(0.72, kernel->evaluate(0.2));
    EXPECT_DOUBLE_EQ(0.63, kernel->evaluate(0.4));
    EXPECT_DOUBLE_EQ(0.48, kernel->evaluate(0.6));
    EXPECT_DOUBLE_EQ(0.27, kernel->evaluate(0.8));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(1.0));

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate method of a 2nd-order biweight kernel
TEST_F(PolynomialKernelTest, EvaluateBiweightKernel)
{
    kernel = new PolynomialKernel(2, 1);
    EXPECT_EQ("2nd-order biweight", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.0000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(0.1215, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(0.3840, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.6615, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(0.8640, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(0.9375, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(0.8640, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.6615, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(0.3840, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(0.1215, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.0000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate method of a 2nd-order triweight kernel
TEST_F(PolynomialKernelTest, EvaluateTriweightKernel)
{
    kernel = new PolynomialKernel(3, 1);
    EXPECT_EQ("2nd-order triweight", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.00000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(0.05103, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(0.28672, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.64827, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(0.96768, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(1.09375, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(0.96768, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.64827, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(0.28672, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(0.05103, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.00000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate method of a 2nd-order general polynomial kernel
TEST_F(PolynomialKernelTest, EvaluateGeneralPolynomialKernel)
{
    kernel = new PolynomialKernel(4, 1);
    EXPECT_EQ("2nd-order polynomial (s = 4)", kernel->get_kernel_name());
    EXPECT_EQ(2, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(0.020667, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(0.206438, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.612615, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.045094, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(1.230469, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.045094, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.612615, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(0.206438, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(0.020667, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 4th-order epanechnikov kernel
TEST_F(PolynomialKernelTest, Evaluate4thOrderEpanechnikov)
{
    kernel = new PolynomialKernel(1, 2);
    EXPECT_EQ("4th-order epanechnikov", kernel->get_kernel_name());
    EXPECT_EQ(4, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.00000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(-0.24975, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(0.14400, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.74025, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.22400, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(1.40625, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.22400, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.74025, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(0.14400, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(-0.24975, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.00000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 6th-order epanechnikov kernel
TEST_F(PolynomialKernelTest, Evaluate6thOrderEpanechnikov)
{
    kernel = new PolynomialKernel(1, 3);
    EXPECT_EQ("6th-order epanechnikov", kernel->get_kernel_name());
    EXPECT_EQ(6, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(-0.100879, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(-0.399840, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.359966, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.517040, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(2.050781, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.517040, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.359966, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(-0.399840, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(-0.100879, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 4th-order biweight kernel
TEST_F(PolynomialKernelTest, Evaluate4thOrderBiweight)
{
    kernel = new PolynomialKernel(2, 2);
    EXPECT_EQ("4th-order biweight", kernel->get_kernel_name());
    EXPECT_EQ(4, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(-0.195615, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(-0.053760, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.601965, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.330560, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(1.640625, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.330560, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.601965, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(-0.053760, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(-0.195615, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 6th-order biweight kernel
TEST_F(PolynomialKernelTest, Evaluate6thOrderBiweight)
{
    kernel = new PolynomialKernel(2, 3);
    EXPECT_EQ("6th-order biweight", kernel->get_kernel_name());
    EXPECT_EQ(6, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(0.063245, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(-0.382234, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.115126, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.534982, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(2.307129, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.534982, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.115126, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(-0.382234, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(0.063245, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 4th-order triweight kernel
TEST_F(PolynomialKernelTest, Evaluate4thOrderTriweight)
{
    kernel = new PolynomialKernel(3, 2);
    EXPECT_EQ("4th-order triweight", kernel->get_kernel_name());
    EXPECT_EQ(4, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(-0.115966, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(-0.154829, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(0.452168, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.393459, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(1.845703, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.393459, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(0.452168, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(-0.154829, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(-0.115966, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// Tests the basic evaluate order of a 6th-order triweight kernel
TEST_F(PolynomialKernelTest, Evaluate6thOrderTriweight)
{
    kernel = new PolynomialKernel(3, 3);
    EXPECT_EQ("6th-order triweight", kernel->get_kernel_name());
    EXPECT_EQ(6, kernel->get_order());

    // test evaluation over the domain
    EXPECT_NEAR(0.000000, kernel->evaluate(-1.0), 1e-6);
    EXPECT_NEAR(0.092135, kernel->evaluate(-0.8), 1e-6);
    EXPECT_NEAR(-0.289530, kernel->evaluate(-0.6), 1e-6);
    EXPECT_NEAR(-0.081026, kernel->evaluate(-0.4), 1e-6);
    EXPECT_NEAR(1.513645, kernel->evaluate(-0.2), 1e-6);
    EXPECT_NEAR(2.537842, kernel->evaluate(0.0), 1e-6);
    EXPECT_NEAR(1.513645, kernel->evaluate(0.2), 1e-6);
    EXPECT_NEAR(-0.081026, kernel->evaluate(0.4), 1e-6);
    EXPECT_NEAR(-0.289530, kernel->evaluate(0.6), 1e-6);
    EXPECT_NEAR(0.092135, kernel->evaluate(0.8), 1e-6);
    EXPECT_NEAR(0.000000, kernel->evaluate(1.0), 1e-6);

    // test evaluation outside the domain
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(-2.0));
    EXPECT_DOUBLE_EQ(0.0, kernel->evaluate(2.0));
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: IntegrateMomentTest
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, Integrate0thMoment)
{
    EXPECT_NEAR(1.000000, kernel1->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(1.000000, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(1.000000, kernel3->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(1.000000, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, Integrate1stMoment)
{
    i = 1;
    EXPECT_NEAR(0.000000, kernel1->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel3->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, Integrate2ndMoment)
{
    i = 2;
    EXPECT_NEAR(0.333333, kernel1->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.200000, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.142857, kernel3->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, Integrate3rdMoment)
{
    i = 3;
    EXPECT_NEAR(0.000000, kernel1->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel3->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.000000, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, Integrate4thMoment)
{
    i = 4;
    EXPECT_NEAR(0.200000, kernel1->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.085714, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(0.047619, kernel3->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.047619, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, IntegrateWithInvalidLimits)
{
    i = 4;

    // case 1 (a < -1.0, b < -1.0)
    a = -5.7;
    b = -2.2;
    EXPECT_DOUBLE_EQ(0.0, kernel2->integrate_moment(a, b, i));
    EXPECT_DOUBLE_EQ(0.0, kernel4->integrate_moment(a, b, i));

    // case 2 (a < -1.0, b = -1.0)
    b = -1.0;
    EXPECT_DOUBLE_EQ(0.0, kernel2->integrate_moment(a, b, i));
    EXPECT_DOUBLE_EQ(0.0, kernel4->integrate_moment(a, b, i));

    // case 3 (a = 1.0, b > 1.0)
    a = 1.0;
    b = 3.9;
    EXPECT_DOUBLE_EQ(0.0, kernel2->integrate_moment(a, b, i));
    EXPECT_DOUBLE_EQ(0.0, kernel4->integrate_moment(a, b, i));

    // case 4 (a > 1.0, b > 1.0)
    a = 1.8;
    EXPECT_DOUBLE_EQ(0.0, kernel2->integrate_moment(a, b, i));
    EXPECT_DOUBLE_EQ(0.0, kernel4->integrate_moment(a, b, i));
}
//---------------------------------------------------------------------------//
TEST_F(IntegrateMomentTest, IntegrateWithValidLimits)
{
    i = 4;

    // case 1 (a = -1.0, b < 1.0)
    b = 0.6;
    EXPECT_NEAR(0.051522, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.017011, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 2 (a = -1.0, b > 1.0)
    b = 4.4;
    EXPECT_NEAR(0.085714, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.047619, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 3 (a > -1.0, b > 1.0)
    a = 0.0;
    EXPECT_NEAR(0.042857, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.023810, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 4 (a > -1.0, b < 1.0)
    b = 0.9;
    EXPECT_NEAR(0.037327, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.012966, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 5 (a > -1.0, b = 1.0)
    b = 1.0;
    EXPECT_NEAR(0.042857, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.023810, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 6 (a < -1.0, b = 1.0)
    a = -7.5;
    EXPECT_NEAR(0.085714, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.047619, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 7 (a < -1.0, b < 1.0)
    b = 0.0;
    EXPECT_NEAR(0.042857, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.023810, kernel4->integrate_moment(a, b, i), 1e-6);

    // case 8 ( a < -1.0, b > 1.0)
    b = 5.3;
    EXPECT_NEAR(0.085714, kernel2->integrate_moment(a, b, i), 1e-6);
    EXPECT_NEAR(-0.047619, kernel4->integrate_moment(a, b, i), 1e-6);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_PolynomialKernel.cpp
