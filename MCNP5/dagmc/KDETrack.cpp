// KDETrack.cpp

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "moab/CartVect.hpp"

#include "KDETrack.hpp"
#include "KDEKernel.hpp"

// quadrature points and weights for the integrate_path_kernel function
const double quad_points[4] = { 0.339981043585, -0.339981043585,
                                0.861136311594, -0.861136311594 };

const double quad_weights[4] = { 0.652145154863, 0.652145154863,
                                 0.347854845137, 0.347854845137 };

// initializing static variable
bool KDETrack::seed_is_set = false;
 
//-----------------------------------------------------------------------------
KDETrack::KDETrack( KDEKernel* k )
{
  
  // set up the track segment
  moab::CartVect Xo( 0, 0, 0 );
  track.Xo = Xo;

  moab::CartVect Uo( 1, 0, 0 );
  track.Uo = Uo;

  track.length = 1.0;
  
  // set the bandwidth, kernel function and neighborhood
  H = Xo + 0.1;

  kernel = k;

  set_neighborhood();

}
//-----------------------------------------------------------------------------
KDETrack::KDETrack( const moab::CartVect & start_point,
                    const moab::CartVect & direction,
                    const moab::CartVect & bandwidth,
                    const double track_length,
                    KDEKernel* k,
                    unsigned int numSubtracks )
: H( bandwidth ), kernel( k )
{
  
  // set up the track segment
  track.Xo = start_point;
  track.Uo = direction;
  track.length = track_length;

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

  // store random points along the track if subtrack estimator is requested
  if ( numSubtracks > 0 )
    subtrack_points = choose_points( numSubtracks );

}
//-----------------------------------------------------------------------------
moab::CartVect KDETrack::get_start_point()
{
  
  return track.Xo;

}
//-----------------------------------------------------------------------------
moab::CartVect KDETrack::get_direction()
{
  
  return track.Uo;

}
//-----------------------------------------------------------------------------
moab::CartVect KDETrack::get_bandwidth()
{
  
  return H;

}
//-----------------------------------------------------------------------------
double KDETrack::get_length()
{
  
  return track.length;

}
//-----------------------------------------------------------------------------
void KDETrack::get_neighborhood( double box_min[3],
                                 double box_max[3],
                                 double & max_radius )
{

  for ( int i = 0 ; i < 3 ; ++i ) {

    box_min[i] = min[i];
    box_max[i] = max[i];

  }

  max_radius = radius;
  
}
//-----------------------------------------------------------------------------
void KDETrack::change_track_segment( const moab::CartVect & newXo,
                                     const moab::CartVect & newUo,
                                     const double newLength )
{

  track_segment newTrack;
  newTrack.Xo = newXo;
  newTrack.Uo = newUo;
  newTrack.length = newLength;

  this->track = newTrack;

  // update neighborhood to reflect new track segment
  set_neighborhood();

  // update subtrack_points to reflect new track segment (if non-empty)
  if ( !subtrack_points.empty() )
    subtrack_points = choose_points( subtrack_points.size() );

}
//-----------------------------------------------------------------------------
void KDETrack::change_bandwidth( const moab::CartVect & newH )
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

  // update neighborhood to reflect new bandwidth
  set_neighborhood();

}
//-----------------------------------------------------------------------------
bool KDETrack::point_within_max_radius( const moab::CartVect & X )
{

  // create a vector from the point Xo to the point X
  moab::CartVect B;

  for ( int i = 0 ; i < 3 ; ++i )
    B[i] = X[i] - track.Xo[i];

  // compute perpendicular distance from X to the line defined by the track
  double r = ( track.Uo * B ).length();
  
  // return true if r is less than max_radius
  if ( r < radius )
    return true;
  else
    return false;

}
//-----------------------------------------------------------------------------
double KDETrack::compute_contribution( const moab::CartVect & X )
{

  if ( subtrack_points.empty() ) {

    // use integral estimator if subtrack_points is empty
    return integrate_path_kernel( X );

  }
  else {

    // otherwise use the subtrack estimator
    return sum_subtracks( X );

  }

}
//-----------------------------------------------------------------------------
double KDETrack::compute_contribution( const double coords[3] )
{

  moab::CartVect X( coords );
  
  return compute_contribution( X );

}
//-----------------------------------------------------------------------------
std::vector<moab::CartVect> KDETrack::choose_points( int p )
{

  // check number of subtrack points is valid
  if ( p <= 0 ) {

    p = 1;
    std::cerr << "\nError! Number of subtrack points is invalid:\n";
    std::cerr << " using p = 1. \n\n";

  }

  std::vector<moab::CartVect> random_points;

  // set random number generator seed if it has not yet been set
  if ( !seed_is_set ) {
    srand( time(NULL) );
    seed_is_set = true;
  }

  // compute subtrack length, assumed to be equal for all subtracks
  double subtrack_length = track.length / p;

  // set the starting point to the beginning of the track
  moab::CartVect start_point = track.Xo;
 
  // choose a random position along each subtrack
  for ( int i = 0 ; i < p ; ++i ) {
    
    double path_length = rand() * subtrack_length / RAND_MAX;
    
    // add the coordinates of the corresponding point
    random_points.push_back( start_point + path_length * track.Uo );
    
    // shift starting point to the next subtrack
    start_point += subtrack_length * track.Uo;
    
  }
  
  return random_points;

}
//-----------------------------------------------------------------------------
void KDETrack::set_neighborhood()
{

  for ( int i = 0 ; i < 3 ; ++i ) {

    // default case where Uo_coord is zero
    min[i] = -H[i] + track.Xo[i];
    max[i] = H[i] + track.Xo[i];
  
    // adjust for direction being positive or negative
    double Uo_coord = track.Uo[i];
  
    if ( Uo_coord > 0 )
      max[i] += track.length * Uo_coord;
    else if ( Uo_coord < 0 )
      min[i] += track.length * Uo_coord;

  }

  // set maximum radius around the track
  radius = H.length();

}
//-----------------------------------------------------------------------------
bool KDETrack::set_integral_limits( const moab::CartVect & X,
                                    double & lower, double & upper )
{
  
  // set initial integral limits to the full track length (default values)
  lower = 0;
  upper = track.length;
  
  // check limits against the valid path length interval for each dimension
  for ( int i = 0 ; i < 3 ; ++i ) {

    double Uo_coord = track.Uo[i];
    double path_min = lower;
    double path_max = upper;
  
    // compute the valid path length interval S = [ path_min, path_max ]
    if ( Uo_coord > 0 ) {
    
      path_min = ( -H[i] + X[i] - track.Xo[i] ) / Uo_coord;
      path_max = ( H[i] + X[i] - track.Xo[i] ) / Uo_coord;

    }
    else if ( Uo_coord < 0 ) {
  
      path_min = ( H[i] + X[i] - track.Xo[i] ) / Uo_coord;
      path_max = ( -H[i] + X[i] - track.Xo[i] ) / Uo_coord;
  
    }

    // set lower limit to highest minimum
    if ( path_min > lower )
      lower = path_min;

    // set upper limit to lowest maximum
    if ( path_max < upper )
      upper = path_max;

  }
  
  // only add limits if the upper limit is greater than the lower limit
  if ( lower < upper )
    return true;
  else
    return false;
    
}
//-----------------------------------------------------------------------------
// NOTE: integrate_path_kernel uses the 4-point gaussian quadrature method
double KDETrack::integrate_path_kernel( const moab::CartVect & X )
{

  // determine the limits of integration
  double lower = 0;
  double upper = 0;
  
  bool limits_set = set_integral_limits( X, lower, upper );

  // compute value of the integral only if valid limits exist
  if ( limits_set ) {

    // define scaling constants
    double c1 = 0.5 * ( upper - lower );
    double c2 = 0.5 * ( upper + lower );
    
    // sum contributions for all quadrature points
    double sum = 0;
    
    for ( int i = 0 ; i < 4 ; ++i ) {
  
      // define scaled quadrature point
      double s = c1 * quad_points[i] + c2;
      
      // compute the value of the kernel function k(X,s)
      double h = 1;
      double k_xyz = 1;
    
      for ( int j = 0 ; j < 3 ; ++j ) {

        double u = ( X[j] - track.Xo[j] - s * track.Uo[j] ) / H[j];
        
        k_xyz *= kernel->evaluate( u );
        h *= H[j];
       
      }
      
      k_xyz = k_xyz / h;
      
      // multiply by quadrature weight and add to sum
      sum += quad_weights[i] * k_xyz;

    }

    // return value of the integral
    return c1 * sum;

  }
  else
    return 0;

}
//-----------------------------------------------------------------------------
double KDETrack::sum_subtracks( const moab::CartVect & X )
{

  double sum = 0;

  std::vector<moab::CartVect>::iterator i;

  for ( i = subtrack_points.begin() ; i != subtrack_points.end() ; ++i ) {

    // compute the value of the kernel function k(X)
    double h = 1; 
    double k_xyz = 1;

    for ( int j = 0 ; j < 3 ; ++j ) {
      
      double u = ( X[j] - (*i)[j] ) / H[j];
 
      k_xyz *= kernel->evaluate( u );
      h *= H[j];

    }

    k_xyz = k_xyz / h;

    // add kernel contribution for subtrack point to sum
    sum += k_xyz;

  }
  
  // normalize by the total number of subtrack points
  sum = sum / subtrack_points.size();

  return sum;

}
//-----------------------------------------------------------------------------
