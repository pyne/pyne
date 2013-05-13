// MCNP5/dagmc/KDEKernel.hpp

#ifndef DAGMC_KDE_KERNEL_H
#define DAGMC_KDE_KERNEL_H

#include <string>

/**
 * \class KDEKernel
 * \brief Defines a general 1D kernel function interface
 *
 * KDEKernel is an abstract Base class that defines the methods that are
 * typically needed to implement a 1D kernel function for use with a Kernel
 * Density Estimator.
 *
 * KDEKernel objects should be created via the createKernel() factory method.
 * This method assigns memory for storing the object, so it will need to be
 * deleted once it is no longer needed to prevent memory leaks.
 *
 * The following functions must be implemented in all Derived classes
 *
 *     1) evaluate(double)
 *     2) evaluate(double, KDEKernel::Boundary)
 *     3) get_kernel_name()
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
     * \brief Defines the location of the boundary for 1D geometries
     */
    enum Boundary {LEFT = 0, RIGHT = 1};

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
     * \brief evaluate the kernel function K(u)
     * \param u the value at which the kernel will be evaluated
     * \return the kernel function evaluation K(u)
     */
    virtual double evaluate(double u) = 0;

    /**
     * \brief evaluate the boundary kernel function Kb(u)
     * \param u the value at which the boundary kernel will be evaluated
     * \param side the location of the boundary
     * \return the boundary kernel function evaluation Kb(u)
     */
    virtual double evaluate(double u, KDEKernel::Boundary side) = 0;

    /**
     * \brief get_kernel_name()
     * \return string representing kernel name
     */
    virtual std::string get_kernel_name() = 0;
};

#endif // DAGMC_KDE_KERNEL_H

// end of MCNP5/dagmc/KDEKernel.hpp
