// MCNP5/dagmc/KDEKernel.cpp

#include <cassert>

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
// TODO implement this method for 2D, and 3D cases
double KDEKernel::boundary_correction(const double* u,
                                      const double* p,
                                      const unsigned int* side,
                                      unsigned int num_corrections) const
{
    assert(num_corrections <= 3);

    // use 1D boundary kernel factor if correction is only needed in 1D
    double correction_factor = 0.0;

    if (num_corrections == 1)
    {
        // compute partial moments ai(p)
        std::vector<double> ai;
        compute_moments(u[0], p[0], side[0], ai);

        // compute boundary correction only if all three moments were computed
        if (ai.size() == 3)
        {
            correction_factor = (ai[2] - ai[1] * u[0]);
            correction_factor /= ai[0] * ai[2] - ai[1] * ai[1];
        }
        // else outside boundary kernel domain
    }
    else
    {
        std::cout << "Not implemented yet" << std::endl;
    }

    return correction_factor;
}
//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
void KDEKernel::compute_moments(double u,
                                double p,
                                unsigned int side,
                                std::vector<double>& moments) const
{
    assert(side <= 1);
    assert(moments.empty());

    // test if outside domain p = [0, 1]
    if (p < 0.0 || p > 1.0) return;

    // determine the integration limits
    double u_min = -1.0;
    double u_max = 1.0;

    if (p < 1.0)
    {
        if (side == 0) // side == LOWER
        {
            u_max = p;
        }
        else // side == UPPER
        {
            u_min = -1.0 * p;
        }
    }

    // test if outside domain u = [u_min, u_max]
    if (u < u_min || u > u_max) return;

    // evaluate the partial moment functions ai(p) and add to moments vector
    moments.push_back(this->integrate_moment(u_min, u_max, 0));
    moments.push_back(this->integrate_moment(u_min, u_max, 1));
    moments.push_back(this->integrate_moment(u_min, u_max, 2));
}
//---------------------------------------------------------------------------//
double KDEKernel::MomentFunction::evaluate(double x) const
{
    double value = kernel.evaluate(x);

    if (moment_index > 0)
    {
        for (unsigned int i = 0; i < moment_index; ++i)
        {
            value *= x;
        } 

    }

    return value;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDEKernel.cpp
