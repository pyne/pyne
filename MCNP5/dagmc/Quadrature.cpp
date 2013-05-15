// MCNP5/dagmc/Quadrature.cpp

#include <cassert>
#include <iostream>

#include "Quadrature.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
Quadrature::Quadrature(unsigned int n) : num_quad_points(n)
{
    // set quadrature points and weights
    set_up_quadrature();

    assert(quad_points.size() == quad_weights.size());
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
double Quadrature::integrate(double a, double b, const Function& f) const
{
    // define scaling constants
    double c1 = 0.5 * (b - a);
    double c2 = 0.5 * (b + a);

    // sum contributions for all quadrature points
    double sum = 0;

    for (unsigned int i = 0; i < quad_points.size(); ++i)
    {
        // defined scaled quadrature point
        double x = c1 * quad_points[i] + c2;

        // add contribution for this quadrature point to the sum
        sum += quad_weights[i] * f.evaluate(x);
    } 

    // return the value of the definite integral of f(x) from a to b
    return c1 * sum;
}
//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
void Quadrature::set_up_quadrature()
{
    switch (num_quad_points)
    {
      case 1:
        quad_points.push_back(0.000000000000);
        quad_weights.push_back(2.000000000000);
        break;

      case 2:
        quad_points.push_back(0.577350269190);
        quad_points.push_back(-0.577350269190);

        quad_weights.push_back(1.000000000000);
        quad_weights.push_back(1.000000000000);
        break;

      case 3:
        quad_points.push_back(0.000000000000);
        quad_points.push_back(0.774596669241);
        quad_points.push_back(-0.774596669241);

        quad_weights.push_back(0.888888888889);
        quad_weights.push_back(0.555555555556);
        quad_weights.push_back(0.555555555556);
        break;

      case 4:
        quad_points.push_back(0.339981043585);
        quad_points.push_back(-0.339981043585);
        quad_points.push_back(0.861136311594);
        quad_points.push_back(-0.861136311594);

        quad_weights.push_back(0.652145154863);
        quad_weights.push_back(0.652145154863);
        quad_weights.push_back(0.347854845137);
        quad_weights.push_back(0.347854845137);
        break;

      default:
        std::cerr << "Warning: " << num_quad_points << " quadrature points "
                  << "is not supported" << std::endl;
        std::cerr << "    using default value of n = 4" << std::endl;

        num_quad_points = 4;
        set_up_quadrature();
    }    
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/Quadrature.cpp
