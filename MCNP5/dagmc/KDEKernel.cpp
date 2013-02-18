// KDEKernel.cpp

#include <cmath>

#include "KDEKernel.hpp"

// initializing static variable
const char* const KDEKernel::kernel_names[] = { "epanechnikov", "uniform" };

//-----------------------------------------------------------------------------
KDEKernel::KernelType KDEKernel::get_type()
{

  return type;

}
//-----------------------------------------------------------------------------
void KDEKernel::change_type( KernelType k )
{

  type = k;

}
//-----------------------------------------------------------------------------
double KDEKernel::evaluate( double u )
{

  double value = 0;

  switch ( type ) {
  
    case EPANECHNIKOV :

      value = epanechnikov( u );
      break;

    case UNIFORM :

      value = uniform( u );
      break;

  }

  return value;

}
//-----------------------------------------------------------------------------
double KDEKernel::epanechnikov( double u )
{

  double value = 0;
  
  if ( fabs( u ) <= 1 )
    value = 0.75 * ( 1 - pow( u, 2 ) );


  return value;

}
//-----------------------------------------------------------------------------
double KDEKernel::uniform( double u )
{

  double value = 0;
  
  if ( fabs( u ) <= 1 )
    value = 0.5;

  return value;

}
//-----------------------------------------------------------------------------
