// KDEKernel.hpp

#ifndef KDEKERNEL_H
#define KDEKERNEL_H

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
    enum KernelType { EPANECHNIKOV = 0, UNIFORM };

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
     * Changes the type of kernel function being used by this Kernel object to
     * the type specified by the parameter k.
     *
     * @param k the kernel type
     */
    void change_type( KernelType k );

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
