// MCNP5/dagmc/KDEKernel.cpp

#include <cassert>
#include <cmath>

#include "KDEKernel.hpp"
#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// FACTORY METHOD
//---------------------------------------------------------------------------//
KDEKernel* KDEKernel::createKernel(const std::string& type, unsigned int order)
{
    KDEKernel* kernel = NULL;

    // order must be a multiple of 2 as only symmetric kernels are supported
    if (order % 2 == 0 && order > 0)
    {
        int r = order / 2;

        if (type == "uniform")
        {
            kernel = new PolynomialKernel(0, r);
        }
        else if (type == "epanechnikov")
        {
            kernel = new PolynomialKernel(1, r);
        }
        else if (type == "biweight")
        {
            kernel = new PolynomialKernel(2, r);
        }
        else if (type == "triweight")
        {
            kernel = new PolynomialKernel(3, r);
        }
    }     

    return kernel;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
double KDEKernel::evaluate(double u,
                           double bandwidth,
                           double distance,
                           unsigned int side) const
{
    assert(side <= 1);

    // compute the scaled distance from the boundary
    double p = distance / bandwidth;

    // test if outside domain p = [0, 1]
    if (p < 0.0 || p > 1.0) return 0.0;

    // determine the integration limits
    double lower_limit = -1.0;
    double upper_limit = 1.0;

    if (p < 1.0)
    {
        if (side == 0) // side == LOWER
        {
            upper_limit = p;
        }
        else // side == UPPER
        {
            lower_limit = -1.0 * p;
        }
    }

    // test if outside domain u = [lower_limit, upper_limit]
    if (u < lower_limit || u > upper_limit) return 0.0;

    // evaluate the moment functions ai(p)
    double a0 = this->integrate_moment(lower_limit, upper_limit, 0);
    double a1 = this->integrate_moment(lower_limit, upper_limit, 1);
    double a2 = this->integrate_moment(lower_limit, upper_limit, 2);

    // compute the value of the boundary kernel
    double value = (a2 - a1 * u) * this->evaluate(u);
    value /= a0 * a2 - a1 * a1;

    return value;
}
//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
double KDEKernel::MomentFunction::evaluate(double x) const
{
    if (moment_index == 0)
    {
        return kernel.evaluate(x);
    }
    else
    { 
        return pow(x, moment_index) * kernel.evaluate(x);
    }
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDEKernel.cpp
