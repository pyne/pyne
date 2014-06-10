// MCNP5/dagmc/KDEKernel.cpp

#include <cassert>
#include <iostream>

#include "KDEKernel.hpp"
#include "PolynomialKernel.hpp"

//---------------------------------------------------------------------------//
// FACTORY METHOD
//---------------------------------------------------------------------------//
KDEKernel* KDEKernel::createKernel(const std::string& type, unsigned int order)
{
    KDEKernel* kernel = NULL;

    // order must be a multiple of 2 as only symmetric kernels are supported
    if (order % 2 == 0 && order > 0)
    {
        int r = order / 2;

        if (type == "uniform")
        {
            kernel = new PolynomialKernel(0, r);
        }
        else if (type == "epanechnikov")
        {
            kernel = new PolynomialKernel(1, r);
        }
        else if (type == "biweight")
        {
            kernel = new PolynomialKernel(2, r);
        }
        else if (type == "triweight")
        {
            kernel = new PolynomialKernel(3, r);
        }
    }     

    return kernel;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
double KDEKernel::boundary_correction(const double* u,
                                      const double* p,
                                      const unsigned int* side,
                                      unsigned int num_corrections) const
{
    assert(num_corrections <= 3);
    assert(num_corrections > 0);

    // compute partial moments ai(p) for first dimension
    std::vector <double> ai_u;
    bool valid_moments = compute_moments(u[0], p[0], side[0], ai_u);

    // check within boundary kernel domain
    if (!valid_moments) return 0.0;

    // solve for the boundary correction factor
    if (num_corrections == 1)
    {
        double correction_factor = (ai_u[2] - ai_u[1] * u[0]);
        correction_factor /= ai_u[0] * ai_u[2] - ai_u[1] * ai_u[1];
        return correction_factor;
    }
    else  // correction needed in more than one dimension
    {
        // compute partial moments ai(p) for second dimension
        std::vector<double> ai_v;
        valid_moments = compute_moments(u[1], p[1], side[1], ai_v);

        // check still within boundary kernel domain
        if (!valid_moments) return 0.0;

        // create coefficients vector initially with right-hand side values
        std::vector<double> coefficients;
        coefficients.push_back(1.0);
        coefficients.push_back(0.0);
        coefficients.push_back(0.0);

        // solve for the coefficients of the boundary correction factor
        bool solved = false;
        std::vector<double> correction_matrix;

        if (num_corrections == 2)
        {
            // solve 3x3 matrix system to get coefficients
            get_correction_matrix2D(ai_u, ai_v, correction_matrix);
            solved = solve_symmetric_matrix(correction_matrix, coefficients);

            if (!solved) return 0.0;
        }
        else  // correction needed in all three dimensions
        {
            coefficients.push_back(0.0);

            // compute partial moments ai(p) for third dimension
            std::vector<double> ai_w;
            valid_moments = compute_moments(u[2], p[2], side[2], ai_w); 

            // check still within boundary kernel domain
            if (!valid_moments) return 0.0;

            // solve 4x4 matrix system to get coefficients
            get_correction_matrix3D(ai_u, ai_v, ai_w, correction_matrix);
            solved = solve_symmetric_matrix(correction_matrix, coefficients);

            if (!solved) return 0.0;
        }

        // compute the boundary correction factor from coefficients
        double correction_factor = coefficients[0];

        for (unsigned int i = 1; i <= num_corrections; ++i)
        {
            correction_factor += u[i-1] * coefficients[i];
        }

        return correction_factor;
    }
}
//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
bool KDEKernel::compute_moments(double u,
                                double p,
                                unsigned int side,
                                std::vector<double>& moments) const
{
    assert(side <= 1);
    assert(moments.empty());

    // make sure p is not negative
    if (p < 0.0) return false;

    // determine the integration limits
    double u_min = -1.0;
    double u_max = 1.0;

    if (p < 1.0)
    {
        if (side == 0) // side == LOWER
        {
            u_max = p;
        }
        else // side == UPPER
        {
            u_min = -1.0 * p;
        }
    }

    // test if outside domain u = [u_min, u_max]
    if (u < u_min || u > u_max) return false;

    // evaluate the partial moment functions ai(p) and add to moments vector
    moments.push_back(this->integrate_moment(u_min, u_max, 0));
    moments.push_back(this->integrate_moment(u_min, u_max, 1));
    moments.push_back(this->integrate_moment(u_min, u_max, 2));

    return true;
}
//---------------------------------------------------------------------------//
void KDEKernel::get_correction_matrix2D(const std::vector<double>& ai_u,
                                        const std::vector<double>& ai_v,
                                        std::vector<double>& matrix) const
{
    assert(matrix.empty());

    // populate matrix elements in lower triangular format using moments
    matrix.push_back(ai_u[0] * ai_v[0]);
    matrix.push_back(ai_u[1] * ai_v[0]);
    matrix.push_back(ai_u[0] * ai_v[1]);
    matrix.push_back(ai_u[2] * ai_v[0]);
    matrix.push_back(ai_u[1] * ai_v[1]);
    matrix.push_back(ai_u[0] * ai_v[2]);
}
//---------------------------------------------------------------------------//
void KDEKernel::get_correction_matrix3D(const std::vector<double>& ai_u,
                                        const std::vector<double>& ai_v,
                                        const std::vector<double>& ai_w,
                                        std::vector<double>& matrix) const
{
    assert(matrix.empty());

    // populate matrix elements in lower triangular format using moments
    matrix.push_back(ai_u[0] * ai_v[0] * ai_w[0]);
    matrix.push_back(ai_u[1] * ai_v[0] * ai_w[0]);
    matrix.push_back(ai_u[0] * ai_v[1] * ai_w[0]);
    matrix.push_back(ai_u[0] * ai_v[0] * ai_w[1]);
    matrix.push_back(ai_u[2] * ai_v[0] * ai_w[0]);
    matrix.push_back(ai_u[1] * ai_v[1] * ai_w[0]);
    matrix.push_back(ai_u[1] * ai_v[0] * ai_w[1]);
    matrix.push_back(ai_u[0] * ai_v[2] * ai_w[0]);
    matrix.push_back(ai_u[0] * ai_v[1] * ai_w[1]);
    matrix.push_back(ai_u[0] * ai_v[0] * ai_w[2]);
}
//---------------------------------------------------------------------------//
bool KDEKernel::solve_symmetric_matrix(std::vector<double>& A,
                                       std::vector<double>& b) const
{
    int n = b.size();

    // LAPACK routine DSPSV input variables
    char uplo = 'L';  // store symmetric matrix in lower triangular format    
    int nrhs = 1;
    int ldb = n;
    int info = 0;
    std::vector<int> ipiv(n, 0);

    dspsv_(&uplo, &n, &nrhs, &A[0], &ipiv[0], &b[0], &ldb, &info);

    if (info == 0)
    {
        return true;
    }
    else  // matrix system could not be solved
    {
        std::cerr << "Warning: LAPACK could not solve symmetric matrix system"
                  << std::endl;
        return false;
    }
}
//---------------------------------------------------------------------------//
double KDEKernel::MomentFunction::evaluate(double x) const
{
    double value = kernel.evaluate(x);

    if (moment_index > 0)
    {
        for (unsigned int i = 0; i < moment_index; ++i)
        {
            value *= x;
        } 

    }

    return value;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/KDEKernel.cpp
