// MCNP5/dagmc/KDEKernel.cpp

#include "KDEKernel.hpp"
#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// FACTORY METHOD
//---------------------------------------------------------------------------//
KDEKernel* KDEKernel::createKernel(std::string type, unsigned int order)
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

// end of MCNP5/dagmc/KDEKernel.cpp
