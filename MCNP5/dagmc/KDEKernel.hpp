// KDEKernel.hpp

#ifndef KDEKERNEL_H
#define KDEKERNEL_H

// forward declarations
class string;

/**  
 * A class that represents a one-dimensional kernel function k(u), which is a
 * function whose integral over the entire domain is equal to one.
 */
class KDEKernel {

  public:

    /**
     * An enumerative type that specifies the kernel function being used by a
     * Kernel object.
     */
    enum KernelType { EPANECHNIKOV = 0, UNIFORM = 1 };
    static const char* const kernel_names[];
    
     /**
      * Constructs a default kernel based on the Epanechnikov kernel function.
      */
    KDEKernel() : type( EPANECHNIKOV ) {}
    
     /**
      * Constructs a kernel based on the kernel function k.
      */
    KDEKernel( KernelType k ) : type( k ) {}

    /**
     * Returns the type of kernel being used by this Kernel object.
     *
     * @return the kernel type
     */
    KernelType get_type();

    /** 
     * \brief Changes the type of this kernel function
     * \param new_type name of the new kernel type in lower case
     *
     * Available options include "epanechnikov" and "uniform".
     */
    void change_type(const std::string& new_type);

    /**
     * \brief Changes the type of this kernel function
     * \param new_type enum value representing new kernel type
     */
    void change_type(KernelType new_type);

   /**
     * Evaluates the kernel function at the parameter u.
     *
     * @param u the value at which the kernel is evaluated
     * @return the kernel function evaluation k(u)
     */
    double evaluate( double u );
  
  private:

    KernelType type;  // the kernel function being used

    /**
     * Evaluates the kernel function at the parameter u using the Epanechnikov
     * kernel. 
     *
     * @param u the value at which the kernel is evaluated
     * @return the epanechnikov kernel function evaluation k(u)
     */
    double epanechnikov( double u );

    /**
     * Evaluates the kernel function at the parameter u using the Uniform
     * kernel. 
     *
     * @param u the value at which the kernel is evaluated
     * @return the uniform kernel function evaluation k(u)
     */
    double uniform( double u );

};

#endif
