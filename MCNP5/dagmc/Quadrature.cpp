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
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void Quadrature::change_quadrature_set(unsigned int new_n)
{
    // reset the original quadrature points and weights
    quad_points.clear();
    quad_weights.clear();

    // set the new quadrature points and weights
    num_quad_points = new_n;
    set_up_quadrature();
}
//---------------------------------------------------------------------------//
double Quadrature::integrate(double a, double b, const Function& f) const
{
    assert(b > a);

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
unsigned int Quadrature::get_num_quad_points() const
{
    return num_quad_points;
}
//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
void Quadrature::set_up_quadrature()
{
    switch (num_quad_points)
    {
      case 1:
        quad_points.push_back(0.0000000000000000);
        quad_weights.push_back(2.0000000000000000);
        break;

      case 2:
        quad_points.push_back(-0.5773502691896257);
        quad_points.push_back(0.5773502691896257);

        quad_weights.push_back(1.0000000000000000);
        quad_weights.push_back(1.0000000000000000);
        break;

      case 3:
        quad_points.push_back(0.0000000000000000);
        quad_points.push_back(-0.7745966692414834);
        quad_points.push_back(0.7745966692414834);

        quad_weights.push_back(0.8888888888888888);
        quad_weights.push_back(0.5555555555555556);
        quad_weights.push_back(0.5555555555555556);
        break;

      case 4:
        quad_points.push_back(-0.3399810435848563);
        quad_points.push_back(0.3399810435848563);
        quad_points.push_back(-0.8611363115940526);
        quad_points.push_back(0.8611363115940526);

        quad_weights.push_back(0.6521451548625461);
        quad_weights.push_back(0.6521451548625461);
        quad_weights.push_back(0.3478548451374538);
        quad_weights.push_back(0.3478548451374538);
        break;

      case 5:
        quad_points.push_back(0.0000000000000000);
        quad_points.push_back(-0.5384693101056831);
        quad_points.push_back(0.5384693101056831);
        quad_points.push_back(-0.9061798459386640);
        quad_points.push_back(0.9061798459386640);

        quad_weights.push_back(0.5688888888888889);
        quad_weights.push_back(0.4786286704993665);
        quad_weights.push_back(0.4786286704993665);
        quad_weights.push_back(0.2369268850561891);
        quad_weights.push_back(0.2369268850561891);
        break;

      case 6:
        quad_points.push_back(0.6612093864662645);
        quad_points.push_back(-0.6612093864662645);
        quad_points.push_back(-0.2386191860831969);
        quad_points.push_back(0.2386191860831969);
        quad_points.push_back(-0.9324695142031521);
        quad_points.push_back(0.9324695142031521);

        quad_weights.push_back(0.3607615730481386);
        quad_weights.push_back(0.3607615730481386);
        quad_weights.push_back(0.4679139345726910);
        quad_weights.push_back(0.4679139345726910);
        quad_weights.push_back(0.1713244923791704);
        quad_weights.push_back(0.1713244923791704);
        break;

      case 7:
        quad_points.push_back(0.0000000000000000);
        quad_points.push_back(0.4058451513773972);
        quad_points.push_back(-0.4058451513773972);
        quad_points.push_back(-0.7415311855993945);
        quad_points.push_back(0.7415311855993945);
        quad_points.push_back(-0.9491079123427585);
        quad_points.push_back(0.9491079123427585);

        quad_weights.push_back(0.4179591836734694);
        quad_weights.push_back(0.3818300505051189);
        quad_weights.push_back(0.3818300505051189);
        quad_weights.push_back(0.2797053914892766);
        quad_weights.push_back(0.2797053914892766);
        quad_weights.push_back(0.1294849661688697);
        quad_weights.push_back(0.1294849661688697);
        break;

      case 8:
        quad_points.push_back(-0.1834346424956498);
        quad_points.push_back(0.1834346424956498);
        quad_points.push_back(-0.5255324099163290);
        quad_points.push_back(0.5255324099163290);
        quad_points.push_back(-0.7966664774136267);
        quad_points.push_back(0.7966664774136267);
        quad_points.push_back(-0.9602898564975363);
        quad_points.push_back(0.9602898564975363);

        quad_weights.push_back(0.3626837833783620);
        quad_weights.push_back(0.3626837833783620);
        quad_weights.push_back(0.3137066458778873);
        quad_weights.push_back(0.3137066458778873);
        quad_weights.push_back(0.2223810344533745);
        quad_weights.push_back(0.2223810344533745);
        quad_weights.push_back(0.1012285362903763);
        quad_weights.push_back(0.1012285362903763);
        break;

      case 9:
        quad_points.push_back(0.0000000000000000);
        quad_points.push_back(-0.8360311073266358);
        quad_points.push_back(0.8360311073266358);
        quad_points.push_back(-0.9681602395076261);
        quad_points.push_back(0.9681602395076261);
        quad_points.push_back(-0.3242534234038089);
        quad_points.push_back(0.3242534234038089);
        quad_points.push_back(-0.6133714327005904);
        quad_points.push_back(0.6133714327005904);

        quad_weights.push_back(0.3302393550012598);
        quad_weights.push_back(0.1806481606948574);
        quad_weights.push_back(0.1806481606948574);
        quad_weights.push_back(0.0812743883615744);
        quad_weights.push_back(0.0812743883615744);
        quad_weights.push_back(0.3123470770400029);
        quad_weights.push_back(0.3123470770400029);
        quad_weights.push_back(0.2606106964029354);
        quad_weights.push_back(0.2606106964029354);
        break;

      case 10:
        quad_points.push_back(-0.1488743389816312);
        quad_points.push_back(0.1488743389816312);
        quad_points.push_back(-0.4333953941292472);
        quad_points.push_back(0.4333953941292472);
        quad_points.push_back(-0.6794095682990244);
        quad_points.push_back(0.6794095682990244);
        quad_points.push_back(-0.8650633666889845);
        quad_points.push_back(0.8650633666889845);
        quad_points.push_back(-0.9739065285171717);
        quad_points.push_back(0.9739065285171717);

        quad_weights.push_back(0.2955242247147529);
        quad_weights.push_back(0.2955242247147529);
        quad_weights.push_back(0.2692667193099963);
        quad_weights.push_back(0.2692667193099963);
        quad_weights.push_back(0.2190863625159820);
        quad_weights.push_back(0.2190863625159820);
        quad_weights.push_back(0.1494513491505806);
        quad_weights.push_back(0.1494513491505806);
        quad_weights.push_back(0.0666713443086881);
        quad_weights.push_back(0.0666713443086881);
        break;

      default:
        std::cerr << "Warning: " << num_quad_points << " quadrature points "
                  << "is not supported" << std::endl;
        std::cerr << "    using default value of n = 10" << std::endl;

        num_quad_points = 10;
        set_up_quadrature();
    } 

    assert(quad_points.size() == num_quad_points);
    assert(quad_points.size() == quad_weights.size());   
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/Quadrature.cpp
