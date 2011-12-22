// KDETrack.hpp

#ifndef KDETRACK_H
#define KDETRACK_H

#include <vector>

#include "moab/CartVect.hpp"

#include "KDEKernel.hpp"

/** 
 * A class that represents the kernel density estimator contribution for one 
 * track segment from a particle's history in a Monte Carlo particle transport
 * simulation. 
 */
class KDETrack {

  public:

    /**
     * Default constructor.  Constructs a kernel density estimator object for a
     * track segment of length 1.0 beginning at (0,0,0) with direction (1,0,0)
     * and bandwidth (0.1,0.1,0.1).  All contributions are computed using the
     * Epanechnikov kernel and the integral estimator.
     */
    KDETrack();

    /**
     * Constructs a kernel density estimator object for a track segment based
     * on the bandwidth H.  The default mode for this function is to compute
     * all contributions using the Epanechnikov kernel and the integral
     * estimator.
     *
     * NOTE: To use the subtrack estimator instead of the integral estimator,
     * simply add a non-zero value for the "numSubtracks" parameter.
     *
     * @param start_point the starting location of the track (xo, yo, zo)
     * @param direction the direction the particle is traveling (uo, vo, wo)
     * @param bandwidth the set of bandwidth values (hx, hy, hz)
     * @param track_length the length of the track segment
     * @param numSubtracks (optional) the number of subtracks to be used
     * @param k (optional) the KernelType function to be used in the computation
     */
    KDETrack( const moab::CartVect & start_point,
              const moab::CartVect & direction,
              const moab::CartVect & bandwidth,
              const double track_length,
              unsigned int numSubtracks = 0,
              KDEKernel::KernelType k = KDEKernel::EPANECHNIKOV );

    /**
     * Copy constructor.
     */
    KDETrack( const KDETrack & obj );

    /**
     * Destructor.
     */
    ~KDETrack();

    /**
     * Overload assignment =.
     */
    KDETrack & operator=( const KDETrack & obj );

    /**
     * @return the starting point Xo = (xo, yo, zo)
     */
    moab::CartVect get_start_point();
    
    /**
     * @return the direction Uo = (uo, vo, wo)
     */
    moab::CartVect get_direction();

    /**
     * @return the set of bandwidth values H = (hx, hy, hz)
     */
    moab::CartVect get_bandwidth();

    /**
     * @return the length of the track segment
     */
    double get_length();

    /**
     * Gets the parameters of the neighborhood for which the contribution
     * may be non-zero.  This includes the coordinates of the minimum and
     * maximum corners, in addition to the maximum radius around the track.
     *
     * @param box_min returns the minimum corner of the neighborhood
     * @param box_max returns the maximum corner of the neighborhood
     * @param max_radius returns the maximum radius around the track
     */
    void get_neighborhood( double box_min[3],
                           double box_max[3],
                           double & max_radius );
    
    /**
     * Changes the starting point, direction and track length of the track
     * segment.  Also updates the neighborhood and stored subtrack points
     * (if non-empty) to reflect the new track values.
     *
     * @param newXo the new starting point (xo, yo, zo)
     * @param newUo the new direction (uo, vo, wo)
     * @param newLength the new track lenth
     */
    void change_track_segment( const moab::CartVect & newXo,
                               const moab::CartVect & newUo,
                               const double newLength );

    /**
     * Changes the set of bandwidth values and updates the neighborhood to
     * reflect the new values.
     *
     * @param newH the new set of bandwidth values (hx, hy, hz)
     */
    void change_bandwidth( const moab::CartVect & newH );

    /**
     * Tests if the perpendicular distance from the calculation point X to
     * the line defined by the track segment is less than the maximum radius
     * of the neighborhood around the track.  Any point that fails this test
     * is guaranteed to have a zero contribution.
     *
     * @param X the calculation point to be tested
     * @return true if the calculation point lies within the maximum radius
     */
    bool point_within_max_radius( const moab::CartVect & X );

    /**
     * Computes the kernel contribution K_xyz for the track segment using
     * the given calculation point X.
     *
     * @param X the calculation point to be evaluated
     * @return the kernel density estimator contribution at the given point
     */
    double compute_contribution( const moab::CartVect & X );

    /**
     * Computes the kernel contribution K_xyz for the track segment using
     * the given calculation point coordinates.
     *
     * @param coords the calculation point to be evaluated
     * @return the kernel density estimator contribution at the given point
     */
    double compute_contribution( const double coords[3] );

    /**
     * Subdivides the track segment into p subtracks of equal length and
     * randomly choses the coordinates of a point from each subtrack.  This
     * function is called in the constructor if the subtrack estimator is
     * requested (i.e. numSubtracks > 0 ) and the results are stored in the
     * vector subtrack_points.  However, it may also be used to obtain a
     * set of p random locations along the track for computing an optimal
     * bandwidth value for either subtrack or track tallies.
     *
     * @param p the number of subtracks to be used
     * @return the vector containing p randomly chosen points along the track
     */
    std::vector<moab::CartVect> choose_points( int p );

  private:

    /**
     * A struct that represents a track segment.
     */
    struct track_segment {

      moab::CartVect Xo;  
      moab::CartVect Uo;    
      double length;

    };
 
    track_segment track;  // stores information about the track segment
    moab::CartVect H;     // the set of bandwidth values (hx, hy, hz)
    KDEKernel* kernel;    // the kernel to be used in the computation

    // neighborhood parameters
    double min[3];        // minimum corner
    double max[3];        // maximum corner
    double radius;    // maximum radius around the track

    /**
     * Sets the neighborhood around the track segment for which a non-zero
     * contribution may occur.
     */
    void set_neighborhood();

    // subtrack estimator parameters
    static bool seed_is_set;
    std::vector<moab::CartVect> subtrack_points;  // random points along track

    /**
     * The set_integral_limits function determines the limits of integration
     * for the integral of the 3D path-length dependent kernel function k(X,s)
     * with respect to path length s for the given calculation point X.
     *
     * @param X the calculation point to be evaluated
     * @param lower returns the lower limit of integration
     * @param upper returns the upper limit of integration
     * @return true if valid limits of integration exist
     *
     * NOTE: if this function returns false, that means there are no valid
     * limits of integration within the range from 0 to total track length.
     * This essentially means that the resulting integrand over this range
     * would have been zero.
     */
    bool set_integral_limits( const moab::CartVect & X,
                              double & lower, double & upper );

    /**
     * The integrate_path_kernel function computes the integral of the 3D
     * path-length dependent kernel function k(X,s) with respect to s for
     * the given calculation point X, using the limits of integration as
     * determined by the set_integral_limits function.
     *
     * NOTE: this function is only called by the compute_contribution function
     * if numSubtracks is set to zero in the KDETrack constructor.
     *
     * @param X the calculation point to be evaluated
     * @return the value of the path-length dependent kernel integral
     */
    double integrate_path_kernel( const moab::CartVect & X );

    /**
     * The sum_subtracks function computes the average 3D kernel contribution
     * for the number of subtracks requested for the given calculation point X,
     * using the randomly chosen points along the track that were computed by
     * the choose_points function.
     *
     * NOTE: this function is only called by the compute_contribution function
     * if numSubtracks is set to a non-zero value in the KDETrack constructor.
     *
     * @param X the calculation point to be evaluated
     * @return the average kernel contribution along a subtrack
     */
    double sum_subtracks( const moab::CartVect & X );

};

#endif
