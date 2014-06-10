// MCNP5/dagmc/KDEKernel.hpp

#ifndef DAGMC_KDE_KERNEL_HPP
#define DAGMC_KDE_KERNEL_HPP

#include <string>
#include <vector>

#include "Quadrature.hpp"

// LAPACK routine for solving Ax = b using symmetric matrix with packed storage
extern "C" {
    void dspsv_(char* UPLO, int* N, int* NRHS, double* AP, int* IPIV, double* B,
                int* LDB, int* INFO); 
}

//===========================================================================//
/**
 * \class KDEKernel
 * \brief Defines a general 1D kernel function interface
 *
 * KDEKernel is an abstract Base class that defines the methods that are
 * typically needed to implement a 1D kernel function K(u) for use with a
 * Kernel Density Estimator.
 *
 * Some of the most commonly used kernel functions can be created via the
 * createKernel() factory method.  This method allocates memory for storing
 * the object, so it will need to be deleted once it is no longer needed to
 * prevent memory leaks.
 *
 * Once a kernel K(u) has been created, it can then be evaluated using the
 * evaluate(double u) method.
 *
 * If a calculation point lies within one bandwidth of an external boundary,
 * then K(u) should be multiplied by the boundary correction factor computed
 * by the boundary_correction method.  This fixes the boundary bias issue that
 * would otherwise occur, but is currently only valid for 2nd-order kernels.
 *
 * =======================
 * Derived Class Interface
 * =======================
 *
 * The following methods must be implemented by all Derived classes
 *
 *     1) evaluate(double u)
 *     2) get_kernel_name()
 *     3) get_order()
 *     4) get_min_quadrature(unsigned int i)
 *     5) integrate_moment(double a, double b, unsigned int i)
 *
 * To assist Derived classes in implementing method 5, there is a protected
 * MomentFunction class defined within KDEKernel that implements the Function
 * interface.  This class can be used to create general moment functions that
 * can be integrated using the integrate method in the Quadrature class.  See
 * the PolynomialKernel implementation for an example.
 */
//===========================================================================//
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
     * \brief Creates a KDEKernel of the specified type and order
     * \param[in] type the name of the kernel type
     * \param[in] order the order of the kernel
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
    static KDEKernel* createKernel(const std::string& type,
                                   unsigned int order = 2);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Evaluate this kernel function K
     * \param[in] u the value at which K will be evaluated
     * \return K(u)
     */
    virtual double evaluate(double u) const = 0;

    /**
     * \brief get_kernel_name()
     * \return string representing kernel name
     */
    virtual std::string get_kernel_name() const = 0;

    /**
     * \brief get_order()
     * \return the order of the kernel
     */
    virtual int get_order() const = 0;

    /**
     * \brief Gets minimum number of quadrature points based on the ith moment
     * \param[in] i the index representing the ith moment function
     * \return number of quadrature points needed to integrate the ith moment
     */
    virtual int get_min_quadrature(unsigned int i) const = 0;

    /**
     * \brief Integrates the ith moment function for this kernel
     * \param[in] a, b the lower and upper integration limits
     * \param[in] i the index representing the ith moment function
     * \return definite integral of the ith moment function for [a, b]
     */
    virtual double integrate_moment(double a, double b, unsigned int i) const = 0;

    /**
     * \brief Evaluate the boundary correction factor for this kernel function K
     * \param[in] u value(s) at which the kernel is to be evaluated
     * \param[in] p ratio(s) of distance from the boundary divided by bandwidth
     * \param[in] side the location(s) of the boundary (0 = LOWER, 1 = UPPER)
     * \param[in] num_corrections number of dimensions requiring correction
     * \return the boundary correction factor for K(u,v,w)
     *
     * Corrects for the boundary bias using the boundary kernel method.  Note
     * that this method works for 1D, 2D, and 3D boundary corrections, as long
     * as the 1D kernel function K is assumed to be the same in all dimenions.
     *
     * The size of u, p, and side are assumed to be equal to num_corrections.
     * Though the 3D case will work for all boundary and interior calculation
     * points, for efficiency purposes values should only be added for the
     * dimensions that actually need correcting.  For example, if only u and
     * w need correcting for K(u, v, w), then
     *
     *                   u = (u, w)
     *                   p = (p(u), p(w))
     *                side = (side(u), side(w))
     *     num_corrections = 2
     *
     * The form of the boundary correction factor is a0 + a1*u + a2*v + a3*w.
     * It will still need to be multiplied by K(u,v,w) to form a valid kernel
     * contribution.
     *
     * If a correction factor of 0.0 is returned, this means that u lies
     * outside the valid domain for the boundary kernel and there is no valid
     * kernel contribution.
     */
    virtual double boundary_correction(const double* u,
                                       const double* p,
                                       const unsigned int* side,
                                       unsigned int num_corrections) const;

  protected:
    /**
     * \brief Computes partial moments ai(p) for this kernel up to i = 2
     * \param[in] u the value at which the kernel is to be evaluated
     * \param[in] p ratio of the distance from the boundary divided by bandwidth
     * \param[in] side the location of the boundary (0 = LOWER, 1 = UPPER)
     * \param[out] moments an empty vector that will store the new ai(p) values
     * \return true if moments are defined for boundary kernel; false otherwise
     *
     * The partial moments ai(p) are integrals of the ith moment function of a
     * kernel K(u).  They are used to determine the correction factor for the
     * boundary kernel method.  Currently, this method only supports second-
     * order kernel functions.
     *
     * Note that the integration limits of a_i(p) depend on whether a LOWER or
     * UPPER boundary is used.  If LOWER, then the integration is performed on
     * [-1, p].  If UPPER, then the integration is performed on [-p, 1]. If
     * p >= 1 moments will be always be defined on the domain [-1, 1].
     */
    bool compute_moments(double u,
                         double p,
                         unsigned int side,
                         std::vector<double>& moments) const;

    /**
     * \brief Sets up the 3x3 matrix needed to solve for the 2D boundary kernel
     * \param[in] ai_u the set of moments for the u-dimension
     * \param[in] ai_v the set of moments for the v-dimension
     * \param[out] matrix the 3x3 correction matrix
     */
    void get_correction_matrix2D(const std::vector<double>& ai_u,
                                 const std::vector<double>& ai_v,
                                 std::vector<double>& matrix) const;

    /**
     * \brief Sets up the 4x4 matrix needed to solve for the 3D boundary kernel
     * \param[in] ai_u the set of moments for the u-dimension
     * \param[in] ai_v the set of moments for the v-dimension
     * \param[in] ai_w the set of moments for the w-dimension
     * \param[out] matrix the 4x4 correction matrix
     */
    void get_correction_matrix3D(const std::vector<double>& ai_u,
                                 const std::vector<double>& ai_v,
                                 const std::vector<double>& ai_w,
                                 std::vector<double>& matrix) const;

    /**
     * \brief Solve a symmetric matrix system Ax = b
     * \param[in/out] A an NxN symmetric matrix
     * \param[in/out] b the right-hand side vector
     * \return true if matrix system was solved; false otherwise
     *
     * The matrix A should be in lower triangular format, stored by columns.
     *
     * On exit, A will be overwritten by the diagonal matrix obtained through
     * the factorization method that was used to solve the matrix system.  The
     * vector b will also be overwritten with the solution x.
     */
    bool solve_symmetric_matrix(std::vector<double>& A,
                                std::vector<double>& b) const; 

    /**
     * \class MomentFunction
     * \brief Defines the ith moment function for a general kernel object
     */
    class MomentFunction : public Function
    {
      public:
        /**
         * \brief Constructor
         * \param[in] i the index representing the ith moment function
         * \param[in] kernel the kernel for which the moment function is desired
         */
        MomentFunction(unsigned int i, const KDEKernel& kernel)
            : moment_index(i), kernel(kernel) {}

        /**
         * \brief Evaluates the ith moment function
         * \param[in] x the value at which this moment function will be evaluated
         * \return x^i * K(x) where K is the kernel function
         */
        double evaluate(double x) const;

      private:
        unsigned int moment_index;
        const KDEKernel& kernel;
    };
};

#endif // DAGMC_KDE_KERNEL_HPP

// end of MCNP5/dagmc/KDEKernel.hpp
