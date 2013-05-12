// MCNP5/dagmc/KDEKernel.cpp

#include "KDEKernel.hpp"
#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// FACTORY METHOD
//---------------------------------------------------------------------------//
KDEKernel* KDEKernel::createKernel(std::string type, unsigned int order)
{
    KDEKernel* kernel = NULL;
    int s = 0;
    int r = order / 2;

    if (type == "uniform" && r == 1)
    {
        kernel = new PolynomialKernel(s, r);
    }
    else if (type == "epanechnikov" && r == 1)
    {
        s = 1;
        kernel = new PolynomialKernel(s, r);
    }

    return kernel;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDEKernel.cpp
