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
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
double PolynomialKernel::set_multiplier()
{
    // set multiplier for 2nd-order kernel function
    double value = pochhammer(0.5, s + 1) / factorial(s);

    // TODO add extra multiplier for kernels of higher order

    return value;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::pochhammer(double x, unsigned int n)
{
    // set default result for (x)0 = 1
    double value = 1.0;

    // compute (x)n
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
long double PolynomialKernel::factorial(unsigned int n)
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
