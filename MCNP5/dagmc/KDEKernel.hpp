// MCNP5/dagmc/KDEKernel.hpp

#ifndef DAGMC_KDE_KERNEL_H
#define DAGMC_KDE_KERNEL_H

#include <string>

#include "Quadrature.hpp"

/**
 * \class KDEKernel
 * \brief Defines a general 1D kernel function interface
 *
 * KDEKernel is an abstract Base class that defines the methods that are
 * typically needed to implement a 1D kernel function K(u) for use with a
 * Kernel Density Estimator.
 *
 * Some of the most commonly used KDEKernel objects can be created via the
 * createKernel() factory method.  This method allocates memory for storing
 * the object, so it will need to be deleted once it is no longer needed to
 * prevent memory leaks.
 *
 * =========================
 * Derived Class Information
 * =========================
 *
 * The following methods must be implemented by all Derived classes
 *
 *     1) evaluate(double)
 *     2) get_kernel_name()
 */
class KDEKernel
{
  protected:
    /**
     * \brief Constructor
     */
    KDEKernel(){}

  public:
    /**
     * \brief Virtual destructor
     */
    virtual ~KDEKernel(){}

    // >>> FACTORY METHOD

    /**
     * \brief creates a KDEKernel of the specified type and order
     * \param type the name of the kernel type
     * \param order the order of the kernel
     * \return pointer to the new KDEKernel that was created
     *
     * Valid kernel types include
     *     "uniform": Uniform kernel (s = 0)
     *     "epanechnikov": Epanechnikov kernel (s = 1)
     *     "biweight": Biweight kernel (s = 2)
     *     "triweight": Triweight kernel (s = 3)
     *
     * Only symmetric kernel functions of these types are supported, which
     * means that the order must be a multiple of 2.
     *
     * NOTE: if an invalid kernel is requested, a NULL pointer is returned.
     */
    static KDEKernel* createKernel(std::string type, unsigned int order = 2);

    // >>> PUBLIC INTERFACE

    /**
     * \brief evaluate the kernel function
     * \param u the value at which the kernel will be evaluated
     * \return K(u)
     */
    virtual double evaluate(double u) const = 0;

    /**
     * \brief get_kernel_name()
     * \return string representing kernel name
     */
    virtual std::string get_kernel_name() const = 0;

  protected:
    /**
     * \class MomentFunction
     * \brief Defines the ith moment function for a general kernel object
     */
    class MomentFunction : public Function
    {
      public:
        /**
         * \brief Constructor
         * \param i the index representing the ith moment function
         * \param kernel the kernel for which the moment function is desired
         */
        MomentFunction(unsigned int i, const KDEKernel& kernel)
            : moment_index(i), kernel(kernel) {}

        /**
         * \brief evaluates the ith moment function
         * \param x the value at which this moment function will be evaluated
         * \return x^i * K(x) where K(x) is the kernel function
         */
        double evaluate(double x) const;

      private:
        unsigned int moment_index;
        const KDEKernel& kernel;
    };
};

#endif // DAGMC_KDE_KERNEL_H

// end of MCNP5/dagmc/KDEKernel.hpp
