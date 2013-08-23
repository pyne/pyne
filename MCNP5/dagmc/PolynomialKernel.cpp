// MCNP5/dagmc/PolynomialKernel.cpp

#include <cassert>
#include <cmath>
#include <sstream>

#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
PolynomialKernel::PolynomialKernel(unsigned int s, unsigned int r)
    : s(s), r(r), multiplier(compute_multiplier())
{
    assert(r > 0);

    // populate coefficients vector for higher order kernel functions
    if (r > 1)
    {
        for (unsigned int k = 0; k < r; ++k)
        {
            double value = pow(-1, k) * pochhammer(0.5 + s + r, k);
            value /= factorial(k) * factorial(r - 1 - k) * pochhammer(1.5, k);
            coefficients.push_back(value);
        }

        assert(coefficients.size() == r);
    }

    // set quadrature for integrating the 0th moment function
    quadrature = new Quadrature(s + r);
}
//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
PolynomialKernel::~PolynomialKernel()
{
    delete quadrature;
}
//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from KDEKernel.hpp
//---------------------------------------------------------------------------//
double PolynomialKernel::evaluate(double u) const
{
    // test if outside kernel function domain [-1.0, 1.0]
    if (u < -1.0 || u > 1.0) return 0.0;

    // evaluate a 2nd-order kernel function
    double value = multiplier;

    if (s == 1) // epanechnikov kernel
    {
        value *= (1 - u * u);
    }
    else if (s > 1)
    {
        value *= pow(1 - u * u, s);
    }

    // multiply value by second polynomial for kernels of higher order
    if (r > 1)
    {
        double sum = 0.0;

        for (unsigned int k = 0; k < r; ++k)
        {
            sum += coefficients[k] * pow(u, 2 * k);
        }

        value *= sum;
    }

    return value;
}
//---------------------------------------------------------------------------//
std::string PolynomialKernel::get_kernel_name() const
{
    // determine the order of this kernel and add to kernel name
    std::stringstream kernel_name;
    kernel_name << 2 * r;

    if (r == 1) kernel_name << "nd-order ";
    else kernel_name << "th-order ";

    // determine type of this kernel and return its full name
    switch (s)
    {
      case 0:
        kernel_name << "uniform";
        break;

      case 1:
        kernel_name << "epanechnikov";
        break;

      case 2:
        kernel_name << "biweight";
        break;

      case 3:
        kernel_name << "triweight";
        break;

      default:
        kernel_name << "polynomial (s = " << s << ")";
    }

    return kernel_name.str();
}
//---------------------------------------------------------------------------//
int PolynomialKernel::get_order() const
{
    return 2 * r;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::integrate_moment(double a,
                                          double b,
                                          unsigned int i) const
{
    assert(quadrature != NULL);

    double value = 0.0;

    // check if integral limits are within the domain u = [-1, 1]
    if (a < 1.0 && b > -1.0)
    {
        // create the ith moment function
        MomentFunction moment(i, *this);

        // define the quadrature set for integrating the ith moment function
        unsigned int n = s + r + (i/2);

        if (quadrature->get_num_quad_points() != n)
        {
            quadrature->change_quadrature_set(n);
        }

        // modify integration limits if needed
        if (a < -1.0) a = -1.0;
        if (b > 1.0) b = 1.0;

        // evaluate the integral
        value = quadrature->integrate(a, b, moment);
    }

    return value;
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
double PolynomialKernel::compute_multiplier()
{
    // compute common multiplier term for 2nd-order polynomial kernel
    double value = pochhammer(0.5, s + 1) / factorial(s);

    // add extra factor to multiplier for kernels of higher order
    if (r > 1)
    {
        unsigned int n = r - 1;
        value *= pochhammer(1.5, n) * pochhammer(1.5 + s, n);
        value /= pochhammer(s + 1, n);
    }

    return value;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::pochhammer(double x, unsigned int n) const
{
    // set default result for (x)_0 = 1
    double value = 1.0;

    // compute (x)_n
    if (n > 0)
    {
        while (n != 1)
        {
            --n;
            value *= (x + n);
        }

        value *= x;
    }

    return value;
}
//---------------------------------------------------------------------------//
long double PolynomialKernel::factorial(unsigned int n) const
{
    // set default result for 0! = 1! = 1
    long double value = 1.0;

    // compute n!
    if (n > 1)
    {
        while (n != 1)
        {
            value *= n;
            --n;
        }
    }

    return value;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/PolynomialKernel.cpp
