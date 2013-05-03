// MCNP5/dagmc/PolynomialKernel.hpp

#ifndef DAGMC_POLYNOMIAL_KERNEL_H
#define DAGMC_POLYNOMIAL_KERNEL_H

#include <vector>
#include <string>

#include "KDEKernel.hpp"

/**
 * \class PolynomialKernel
 * \brief Defines a polynomial-based kernel function
 * 
 * TODO add detailed description of class
 */
class PolynomialKernel : public KDEKernel
{
  public:
    /**
     * \brief Constructor
     * \param s the level of smoothness of the kernel
     * \param r defines a kernel of 2rth-order
     */
    PolynomialKernel(unsigned int s, unsigned int r);

    // >>> DERIVED PUBLIC INTERFACE from KDEKernel.hpp

    /**
     * \brief evaluate the kernel function K(u)
     * \param u the value at which the kernel will be evaluated
     * \return the kernel function evaluation K(u)
     */
    virtual double evaluate(double u);

    /**
     * \brief evaluate the boundary kernel function Kb(u)
     * \param u the value at which the boundary kernel will be evaluated
     * \param side the location of the boundary
     * \return the boundary kernel function evaluation Kb(u)
     */
    virtual double evaluate_boundary(double u, KDEKernel::Boundary side);

    /**
     * \brief get_kernel_name()
     * \return string representing kernel name
     */
    virtual std::string get_kernel_name();
  
  private:
    /// Smoothness factor for this kernel
    unsigned int s;

    /// Related to the order of this kernel (order = 2r)
    unsigned int r;

    /// Represents constant multiplier of kernel function
    const double multiplier;

    /// Coefficients of the polynomial generated for kernels of order > 2
    std::vector<double> coefficients;

    // >>> PROTECTED FUNCTIONS

    /**
     * \brief evaluates the pochhammer symbol (x)n
     * \param x the value for which to evaluate the pochhammer symbol
     * \param n any non-negative integer (n >= 0)
     * \return the pochhammer symbol evaluation (x)n
     *
     * Computes (x)n = x(x + 1)(x + 2)...(x + n - 1) where (x)0 = 1.
     */
    double pochhammer(double x, unsigned int n);

    /**
     * \brief evaluates n!
     * \param n any non-negative integer (n >=0)
     * \return the n! evaluation
     *
     * Computes n! = n(n-1)(n-2)...(2)(1) where 0! = 1.
     */
    long double factorial(unsigned int n);
};

#endif // DAGMC_POLYNOMIAL_KERNEL_H

// end of MCNP5/dagmc/PolynomialKernel.hpp
