// MCNP5/dagmc/KDEKernel.hpp

#ifndef DAGMC_KDE_KERNEL_H
#define DAGMC_KDE_KERNEL_H

#include "string"

/**
 * \class KDEKernel
 * \brief Defines a general kernel function interface
 *
 * KDEKernel is a Base class that defines the variables and methods that are
 * typically needed to implement a kernel function for use with a Kernel
 * Density Estimator.
 *
 * KDEKernel objects can be created via the createKernel() factory method.
 * This method assigns memory for storing the object, so it will need to be
 * deleted once it is no longer needed to prevent memory leaks.
 *
 * Currently, no kernel functions are supported.
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
     * \param order the order of the kernel
     */
    explicit KDEKernel(unsigned int order) : kernel_order(order) {}

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
  
  protected:
    /// order of this kernel
    unsigned int kernel_order;
};

#endif // DAGMC_KDE_KERNEL_H

// end of MCNP5/dagmc/KDEKernel.hpp
