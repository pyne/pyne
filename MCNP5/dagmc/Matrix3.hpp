/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**\file Matrix3.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#ifndef MB_MATRIX_3_HPP
#define MB_MATRIX_3_HPP

#include "moab/CartVect.hpp"
#include "moab/Types.hpp"
#include <iosfwd>
#include <limits>

#ifdef _MSC_VER
# define finite _finite
#endif

namespace moab {

class Matrix3 
{
  double d[9];
public:

  inline Matrix3() 
    {}
  
  inline Matrix3( double diagonal )
    { 
      d[0] = d[4] = d[8] = diagonal;
      d[1] = d[2] = d[3] = 0;
      d[5] = d[6] = d[7] = 0;
    }
    
  inline Matrix3( const CartVect& diagonal )
    { 
      d[0] = diagonal[0];
      d[4] = diagonal[1],
      d[8] = diagonal[2];
      d[1] = d[2] = d[3] = 0;
      d[5] = d[6] = d[7] = 0;
    }
  
  
  inline Matrix3( const CartVect& row0,
                    const CartVect& row1,
                    const CartVect& row2 )
    {
      row0.get( d );
      row1.get( d+3 );
      row2.get( d+6 );
    }
  
  inline Matrix3( const double* v )
    { 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
    }
  
  inline Matrix3( double v00, double v01, double v02,
                    double v10, double v11, double v12,
                    double v20, double v21, double v22 )
    {
      d[0] = v00; d[1] = v01; d[2] = v02;
      d[3] = v10; d[4] = v11; d[5] = v12;
      d[6] = v20; d[7] = v21; d[8] = v22;
    }
  
  inline Matrix3( const Matrix3& m )
    {
      d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2];
      d[3] = m.d[3]; d[4] = m.d[4]; d[5] = m.d[5];
      d[6] = m.d[6]; d[7] = m.d[7]; d[8] = m.d[8];
    }
   
  inline Matrix3& operator=( const Matrix3& m )
    {
      d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2];
      d[3] = m.d[3]; d[4] = m.d[4]; d[5] = m.d[5];
      d[6] = m.d[6]; d[7] = m.d[7]; d[8] = m.d[8];
      return *this;
    }
  
  inline Matrix3& operator=( const double* v )
    { 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
      return *this;
    }

  inline double* operator[]( unsigned i )
    { return d + 3*i; }
  inline const double* operator[]( unsigned i ) const
    { return d + 3*i; }
  inline double& operator()(unsigned r, unsigned c)
    { return d[3*r+c]; }
  inline double operator()(unsigned r, unsigned c) const
    { return d[3*r+c]; }
  
  inline Matrix3& operator+=( const Matrix3& m )
    {
      d[0] += m.d[0]; d[1] += m.d[1]; d[2] += m.d[2];
      d[3] += m.d[3]; d[4] += m.d[4]; d[5] += m.d[5];
      d[6] += m.d[6]; d[7] += m.d[7]; d[8] += m.d[8];
      return *this;
    }
  
  inline Matrix3& operator-=( const Matrix3& m )
    {
      d[0] -= m.d[0]; d[1] -= m.d[1]; d[2] -= m.d[2];
      d[3] -= m.d[3]; d[4] -= m.d[4]; d[5] -= m.d[5];
      d[6] -= m.d[6]; d[7] -= m.d[7]; d[8] -= m.d[8];
      return *this;
    }
  
  inline Matrix3& operator*=( double s )
    {
      d[0] *= s; d[1] *= s; d[2] *= s;
      d[3] *= s; d[4] *= s; d[5] *= s;
      d[6] *= s; d[7] *= s; d[8] *= s;
      return *this;
    }
  
  inline Matrix3& operator/=( double s )
    {
      d[0] /= s; d[1] /= s; d[2] /= s;
      d[3] /= s; d[4] /= s; d[5] /= s;
      d[6] /= s; d[7] /= s; d[8] /= s;
      return *this;
    }
  
  inline Matrix3& operator*=( const Matrix3& m );
  
  inline double determinant() const;
  
  inline Matrix3 inverse() const;
  
    // invert matrix without re-calculating the
    // reciprocal of the determinant.
  inline Matrix3 inverse( double inverse_det ) const;
  
  inline bool positive_definite() const;
  inline bool positive_definite( double& determinant_out ) const;
  
  inline Matrix3 transpose() const;
  
  inline bool invert();
  
    // Calculate determinant of 2x2 submatrix composed of the
    // elements not in the passed row or column.
  inline double subdet( int r, int c ) const;
};

inline Matrix3 operator+( const Matrix3& a, const Matrix3& b )
  { return Matrix3(a) += b; }

inline Matrix3 operator-( const Matrix3& a, const Matrix3& b )
  { return Matrix3(a) -= b; }

inline Matrix3 operator*( const Matrix3& a, const Matrix3& b )
{
  return Matrix3( a(0,0) * b(0,0) + a(0,1) * b(1,0) + a(0,2) * b(2,0),
                    a(0,0) * b(0,1) + a(0,1) * b(1,1) + a(0,2) * b(2,1),
                    a(0,0) * b(0,2) + a(0,1) * b(1,2) + a(0,2) * b(2,2),
                    a(1,0) * b(0,0) + a(1,1) * b(1,0) + a(1,2) * b(2,0),
                    a(1,0) * b(0,1) + a(1,1) * b(1,1) + a(1,2) * b(2,1),
                    a(1,0) * b(0,2) + a(1,1) * b(1,2) + a(1,2) * b(2,2),
                    a(2,0) * b(0,0) + a(2,1) * b(1,0) + a(2,2) * b(2,0),
                    a(2,0) * b(0,1) + a(2,1) * b(1,1) + a(2,2) * b(2,1),
                    a(2,0) * b(0,2) + a(2,1) * b(1,2) + a(2,2) * b(2,2) );
}

inline Matrix3& Matrix3::operator*=( const Matrix3& m )
  { return *this = Matrix3(*this) * m; }

inline Matrix3 outer_product( const CartVect& u,
                                const CartVect& v )
{
  return Matrix3( u[0] * v[0], u[0] * v[1], u[0] * v[2],
                    u[1] * v[0], u[1] * v[1], u[1] * v[2],
                    u[2] * v[0], u[2] * v[1], u[2] * v[2] );
}

inline CartVect operator*( const CartVect& v, const Matrix3& m )
{
  return CartVect( v[0] * m(0,0) + v[1] * m(1,0) + v[2] * m(2,0),
                     v[0] * m(0,1) + v[1] * m(1,1) + v[2] * m(2,1),
                     v[0] * m(0,2) + v[1] * m(1,2) + v[2] * m(2,2) );
}

inline CartVect operator*( const Matrix3& m, const CartVect& v )
{
  return CartVect( v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2),
                     v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2),
                     v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2) );
} 

inline double Matrix3::determinant() const
{
  return d[0] * d[4] * d[8] 
       + d[1] * d[5] * d[6]
       + d[2] * d[3] * d[7]
       - d[0] * d[5] * d[7]
       - d[1] * d[3] * d[8]
       - d[2] * d[4] * d[6];
}

inline bool Matrix3::positive_definite( double& det ) const
{
  double subdet6 = d[1]*d[5]-d[2]*d[4];
  double subdet7 = d[2]*d[3]-d[0]*d[5];
  double subdet8 = d[0]*d[4]-d[1]*d[3];
  det = d[6]*subdet6 + d[7]*subdet7 + d[8]*subdet8;
  return d[0] > 0 && subdet8 > 0 && det > 0;
}

inline bool Matrix3::positive_definite() const
{
  double d;
  return positive_definite(d);
}

inline Matrix3 Matrix3::inverse( double i ) const
{
  return Matrix3( i * (d[4] * d[8] - d[5] * d[7]),
                    i * (d[2] * d[7] - d[8] * d[1]),
                    i * (d[1] * d[5] - d[4] * d[2]),
                    i * (d[5] * d[6] - d[8] * d[3]),
                    i * (d[0] * d[8] - d[6] * d[2]),
                    i * (d[2] * d[3] - d[5] * d[0]),
                    i * (d[3] * d[7] - d[6] * d[4]),
                    i * (d[1] * d[6] - d[7] * d[0]),
                    i * (d[0] * d[4] - d[3] * d[1]) );
}  

inline Matrix3 Matrix3::inverse() const
{
  return inverse( 1.0 / determinant() );
}

inline bool Matrix3::invert()
{
  double i = 1.0 / determinant();
  if (!finite(i) || fabs(i) < std::numeric_limits<double>::epsilon())
    return false;
  *this = inverse( i );
  return true;
}

inline Matrix3 Matrix3::transpose() const
{
  return Matrix3( d[0], d[3], d[6],
                    d[1], d[4], d[7],
                    d[2], d[5], d[8] );
}

inline double Matrix3::subdet( int r, int c ) const
{
  const int r1 = (r+1)%3, r2 = (r+2)%3;
  const int c1 = (c+1)%3, c2 = (c+2)%3;
  return d[3*r1+c1]*d[3*r2+c2] - d[3*r1+c2]*d[3*r2+c1];
}
                         
ErrorCode EigenDecomp( const Matrix3& a, 
                         double Eigenvalues[3],
                         CartVect Eigenvectors[3] );

std::ostream& operator<<( std::ostream&, const Matrix3& );
  
} // namespace moab

#endif
