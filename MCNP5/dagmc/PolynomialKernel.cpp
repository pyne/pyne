// MCNP5/dagmc/PolynomialKernel.cpp

#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
PolynomialKernel::PolynomialKernel(unsigned int s, unsigned int r)
    : KDEKernel(2 * r), multiplier(1.0)
{
    this->s = s;
    this->r = r;
}
//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE from KDEKernel.hpp
//---------------------------------------------------------------------------//
double PolynomialKernel::evaluate(double u)
{
    return 0.0;
}
//---------------------------------------------------------------------------//
double PolynomialKernel::evaluate_boundary(double u, KDEKernel::Boundary side)
{
    return 0.0;
}
//---------------------------------------------------------------------------//
std::string PolynomialKernel::get_kernel_name()
{
    return "Not Implemented Yet";
}
//---------------------------------------------------------------------------//
// PROTECTED FUNCTIONS
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
