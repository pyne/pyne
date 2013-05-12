// MCNP5/dagmc/PolynomialKernel.cpp

#include <cmath>
#include <sstream>

#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
PolynomialKernel::PolynomialKernel(unsigned int s, unsigned int r)
    : s(s), r(r), multiplier(set_multiplier())
{
    // TODO populate coefficients vector for higher order kernel functions
}
//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from KDEKernel.hpp
//---------------------------------------------------------------------------//
double PolynomialKernel::evaluate(double u)
{
    // test if outside kernel function domain [-1.0, 1.0]
    if (u < -1.0 || u > 1.0)
    {
        return 0.0;
    }

    // evaluate kernel function
    double value = multiplier;

    if (s == 1) // epanechnikov kernel
    {
        value *= (1 - u * u);
    }
    else if (s > 1)
    {
        value *= pow(1 - u * u, s);
    }

    return value;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::evaluate(double u, KDEKernel::Boundary side)
{
    return 0.0;
}
//---------------------------------------------------------------------------//
std::string PolynomialKernel::get_kernel_name()
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

      default:
        return "Not a valid polynomial kernel";
    }

    return kernel_name.str();
}
//---------------------------------------------------------------------------//
// PROTECTED FUNCTIONS
//---------------------------------------------------------------------------//
double PolynomialKernel::set_multiplier()
{
    double value = 1.0;

    if (s == 0) value = 0.5;
    if (s == 1) value = 0.75;

    return value;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::pochhammer(double x, unsigned int n)
{
    return 0.0;
}
//---------------------------------------------------------------------------//
long double PolynomialKernel::factorial(unsigned int n)
{
    return 0.0;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/PolynomialKernel.cpp
