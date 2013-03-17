// KDECollision.cpp

#include <cmath>
#include <iostream>

#include "moab/CartVect.hpp"

#include "KDECollision.hpp"
#include "KDEKernel.hpp"

//-----------------------------------------------------------------------------
KDECollision::KDECollision( const moab::CartVect & collision_loc,
                            const moab::CartVect & bandwidth,
                            KDEKernel* k )
: Xic( collision_loc ), H( bandwidth ), kernel( k )
{

  // check bandwidth values are valid
  for ( int i = 0 ; i < 3 ; ++i ) {

    if ( H[i] <= 0 ) {

      if ( H[i] == 0 )
        H[i] = 1;
      
      H[i] = fabs( H[i] );
      
      std::cerr << "\nError! KDE bandwidth value h(" << i << ") is invalid:\n";
      std::cerr << " using h(" << i << ") = " << H[i] << " in the tally. \n\n";
      
    }

  }
  
  set_neighborhood();

}
//-----------------------------------------------------------------------------
moab::CartVect KDECollision::get_collision()
{
  
  return Xic;

}
//-----------------------------------------------------------------------------
moab::CartVect KDECollision::get_bandwidth()
{
  
  return H;

}
//-----------------------------------------------------------------------------
void KDECollision::get_neighborhood( double box_min[3], double box_max[3] )
{

  for ( int i = 0 ; i < 3 ; ++i ) {

    box_min[i] = min[i];
    box_max[i] = max[i];

  }

}
//-----------------------------------------------------------------------------
void KDECollision::change_collision( const moab::CartVect & newXic )
{

  this->Xic = newXic;
   set_neighborhood();

}
//-----------------------------------------------------------------------------
void KDECollision::change_bandwidth( const moab::CartVect & newH )
{

  this->H = newH;

  // check bandwidth values are valid
  for ( int i = 0 ; i < 3 ; ++i ) {

    if ( H[i] <= 0 ) {

      if ( H[i] == 0 )
        H[i] = 1;

      H[i] = fabs( H[i] );
      
      std::cerr << "\nError! KDE bandwidth value h(" << i << ") is invalid:\n";
      std::cerr << " using h(" << i << ") = " << H[i] << " in the tally. \n\n";
      
    }

  }

  set_neighborhood();

}
//-----------------------------------------------------------------------------
double KDECollision::compute_contribution( const moab::CartVect & X )
{
  
  double u = 0;
  double h = 1;
  double k_xyz = 1;
  
  // computes the 3D kernel contribution using product of 1D kernels
  for ( int i = 0 ; i < 3 ; ++i ) {
    
    u = ( X[i] - Xic[i] ) / H[i] ;
    k_xyz *= kernel->evaluate( u );
    h *= H[i];

  }      
  
  k_xyz = k_xyz / h;  

  return k_xyz;

}
//-----------------------------------------------------------------------------
double KDECollision::compute_contribution( const double coords[3] )
{
  moab::CartVect X( coords );
  
  return compute_contribution( X );

}
//-----------------------------------------------------------------------------
void KDECollision::set_neighborhood()
{

  for ( int i = 0 ; i < 3 ; ++i ) { 

    min[i] = Xic[i] - H[i];
    max[i] = Xic[i] + H[i];

  }

}
//-----------------------------------------------------------------------------
