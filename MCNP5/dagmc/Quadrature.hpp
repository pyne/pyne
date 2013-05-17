// MCNP5/dagmc/Quadrature.hpp

#ifndef DAGMC_QUADRATURE_H
#define DAGMC_QUADRATURE_H

#include <vector>

/**
 * \class Function
 * \brief Defines an abstract function interface
 *
 * Function defines an abstract interface for creating functions f(x) that can
 * be integrated by the integrate() method in the Quadrature class.
 */
class Function
{
  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~Function(){}

    /**
     * \brief evaluate the function f(x)
     * \param x the value at which this function will be evaluated
     * \return the function evaluation f(x)
     */
    virtual double evaluate(double x) const = 0;
};

/**
 * \class Quadrature
 * \brief Defines a Quadrature scheme for computing definite integrals of f(x)
 *
 * Quadrature is a class that represents a Gaussian Quadrature scheme based on
 * n quadrature points.  This scheme can be used to compute definite integrals
 * for functions that implement the Function interface.
 *
 * NOTE: Polynomials of order 2n - 1 or less are integrated exactly by an
 * n-point Quadrature.  However, the quadrature points and weights are only
 * exact to 12 significant figures.  This may limit the final accuracy of
 * the results due to floating point arithmetic.  Use a higher order
 * Quadrature if you need more accuracy.
 */
class Quadrature
{
  public:
    /**
     * \brief Constructor
     * \param n the number of points to use with this Quadrature
     */
    Quadrature(unsigned int n);

    // >>> PUBLIC INTERFACE

    /**
     * \brief changes the set of quadrature points for this Quadrature
     * \param new_n the new number of quadrature points
     */
    void change_quadrature_set(unsigned int new_n);

    /**
     * \brief computes the definite integral of the function f(x) from a to b
     * \param a the lower integration limit
     * \param b the upper integration limit
     * \param f the function to be integrated
     * \return the value of the definite integral of f(x) from a to b
     */
    double integrate(double a, double b, const Function& f) const;

    /**
     * \brief get_num_quad_points()
     * \return the number of quadrature points used by this Quadrature
     */
    unsigned int get_num_quad_points() const;

  private:
    unsigned int num_quad_points;
    std::vector<double> quad_points;
    std::vector<double> quad_weights;

    /**
     * \brief set up the quadrature points and weights for this Quadrature
     */
    void set_up_quadrature();
};

#endif // DAGMC_QUADRATURE_H

// end of MCNP5/dagmc/Quadrature.hpp
