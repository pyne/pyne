// MCNP5/dagmc/test/test_Quadrature.cpp

#include <cassert>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "../Quadrature.hpp"

//---------------------------------------------------------------------------//
// MOCK OBJECTS
//---------------------------------------------------------------------------//
class PolynomialFunction : public Function
{
  public:
    /**
     * \brief Defines an nth order polynomial function
     * \param coefficients the ith element represents coefficient for x^(i)
     * \param n the order of the polynomial
     */
    PolynomialFunction(std::vector<double> coefficients, unsigned int n)
        : coefficients(coefficients), order(n)
    {
        assert(coefficients.size() == n + 1);
    }

    /**
     * \brief evalute the polynomial p(x)
     * \param x the value at which this polynomial will be evaluated
     * \return the polynomial evaluation p(x)
     */
    double evaluate(double x) const
    {
        double value = 0.0;

        for (unsigned int i = 0 ; i <= order; ++i)
        {
            value += pow(x, i) * coefficients[i];
        }

        return value;
    }

  private:
    std::vector<double> coefficients;
    unsigned int order;
};
//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class QuadratureTest : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
        // define integration limits
        a = -2.0;
        b = 3.6;

        // define default polynomial order and coefficients
        order = 0;
        coefficients.push_back(1.0);

        // set Quadrature and Function objects to NULL
        quadrature = NULL;
        function = NULL;
    }

    // deallocate memory resources
    virtual void TearDown()
    {
        delete quadrature;
        delete function;
    }

  protected:
    // data needed for each test
    double a, b;
    unsigned int order;
    std::vector<double> coefficients;
    Quadrature* quadrature;
    Function* function;
};
//---------------------------------------------------------------------------//
// SIMPLE TESTS
//---------------------------------------------------------------------------//
TEST(InvalidQuadratureTest, InvalidQuadPoints)
{
    Quadrature quadrature(50);
    EXPECT_EQ(10, quadrature.get_num_quad_points());
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: QuadratureTest
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, IntegrateZeroFunction)
{
    // set up polynomial function
    coefficients[0] = 0.0;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(1);
    EXPECT_EQ(1, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(0.0, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, IntegrateConstant)
{
    // set up polynomial function
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(1);
    EXPECT_EQ(1, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(5.6, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate1stOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    order = 1;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(1);
    EXPECT_EQ(1, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(-3.36, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate2ndOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    order = 2;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(2);
    EXPECT_EQ(2, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(51.296, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate3rdOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    order = 3;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(2);
    EXPECT_EQ(2, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(-100.6656, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate4thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    order = 4;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(3);
    EXPECT_EQ(3, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(535.99616, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate5thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    coefficients.push_back(-6.0);
    order = 5;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(3);
    EXPECT_EQ(3, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(-1576.786176, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate6thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    coefficients.push_back(-6.0);
    coefficients.push_back(7.0);
    order = 6;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(4);
    EXPECT_EQ(4, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(6387.630234, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate7thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    coefficients.push_back(-6.0);
    coefficients.push_back(7.0);
    coefficients.push_back(-8.0);
    order = 7;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(4);
    EXPECT_EQ(4, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(-21567.468841, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate8thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    coefficients.push_back(-6.0);
    coefficients.push_back(7.0);
    coefficients.push_back(-8.0);
    coefficients.push_back(9.0);
    order = 8;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(5);
    EXPECT_EQ(5, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(80504.4878275, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, Integrate9thOrderPolynomial)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    coefficients.push_back(-4.0);
    coefficients.push_back(5.0);
    coefficients.push_back(-6.0);
    coefficients.push_back(7.0);
    coefficients.push_back(-8.0);
    coefficients.push_back(9.0);
    coefficients.push_back(-10.0);
    order = 9;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(5);
    EXPECT_EQ(5, quadrature->get_num_quad_points());

    // evaluate the integral
    EXPECT_NEAR(-284087.356179,quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//
TEST_F(QuadratureTest, ChangeQuadratureSet)
{
    // set up polynomial function
    coefficients.push_back(-2.0);
    coefficients.push_back(3.0);
    order = 2;
    function = new PolynomialFunction(coefficients, order);

    // set up quadrature
    quadrature = new Quadrature(2);
    EXPECT_EQ(2, quadrature->get_num_quad_points());

    // change the quadrature set
    quadrature->change_quadrature_set(1);
    EXPECT_EQ(1, quadrature->get_num_quad_points());

    // reevaluate the integral using the lower order quadrature set
    EXPECT_NEAR(7.392, quadrature->integrate(a, b, *function), 1e-6);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/test/test_KDEKernel.cpp
