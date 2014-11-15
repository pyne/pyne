// MCNP5/dagmc/Quadrature.hpp

#ifndef DAGMC_QUADRATURE_HPP
#define DAGMC_QUADRATURE_HPP

#include <vector>

//===========================================================================//
/**
 * \class Function
 * \brief Defines an abstract function interface
 *
 * Function defines an abstract interface for creating functions f(x) that can
 * be integrated by the integrate() method in the Quadrature class.
 */
//===========================================================================//
class Function
{
  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~Function(){}

    /**
     * \brief Evaluate this Function f
     * \param[in] x the value at which f will be evaluated
     * \return f(x)
     */
    virtual double evaluate(double x) const = 0;
};

//===========================================================================//
/**
 * \class Quadrature
 * \brief Defines a Quadrature scheme for computing definite integrals of f(x)
 *
 * Quadrature is a class that represents a Gaussian Quadrature scheme based on
 * n quadrature points.  This scheme can be used to compute definite integrals
 * for classes that implement the Function interface.
 *
 * NOTE: Polynomials of order 2n - 1 or less are integrated exactly by an
 * n-point Quadrature.  However, the quadrature points and weights are only
 * exact to 16 significant figures.  This may limit the final accuracy of
 * the results due to floating point arithmetic.  Use a higher order
 * Quadrature if you need more accuracy.
 */
//===========================================================================//
class Quadrature
{
  public:
    /**
     * \brief Constructor
     * \param[in] n the number of points to use with this Quadrature
     */
    explicit Quadrature(unsigned int n);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Changes the set of quadrature points for this Quadrature
     * \param[in] new_n the new number of quadrature points
     */
    void change_quadrature_set(unsigned int new_n);

    /**
     * \brief Computes the definite integral of a Function f
     * \param[in] a, b the lower and upper integration limits
     * \param[in] f the Function to be integrated
     * \return definite integral of f(x) for [a, b]
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

    // >>> PRIVATE METHODS

    /**
     * \brief Set up the quadrature points and weights for this Quadrature
     */
    void set_up_quadrature();
};

#endif // DAGMC_QUADRATURE_HPP

// end of MCNP5/dagmc/Quadrature.hpp
