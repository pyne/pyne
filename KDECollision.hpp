// KDECollision.hpp

#ifndef KDECOLLISION_H
#define KDECOLLISION_H

#include "moab/CartVect.hpp"

#include "KDEKernel.hpp"

/** 
 * A class that represents the kernel density estimator contribution for one 
 * collision from a particle's history in a Monte Carlo particle transport
 * simulation. 
 */
class KDECollision {

  public:

    /**
     * Constructs a kernel density estimator object for the collision point
     * Xic based on the bandwidth H.
     *
     * @param collision the collision point (xic, yic, zic)
     * @param bandwidth the set of bandwidth values (hx, hy, hz)
     * @param k the kernel function to be used in the computation
     */
    KDECollision( const moab::CartVect & collision,
                  const moab::CartVect & bandwidth,
                  KDEKernel* k );

    /**
     * @return the collision point Xic = (xic, yic, zic)
     */
    moab::CartVect get_collision();
    
    /**
     * @return the set of bandwidth values H = (hx, hy, hz)
     */
    moab::CartVect get_bandwidth();

    /**
     * Gets the coordinates of the minimum and maximum corners of the
     * collision neighborhood for which the contribution is non-zero.
     *
     * @param box_min gets the minimum corner of the neighborhood
     * @param box_max gets the maximum corner of the neighborhood
     */
    void get_neighborhood( double box_min[3], double box_max[3] );
    
    /**
     * Changes the collision point and updates the neighborhood to reflect
     * the new values.
     *
     * @param newXic the new collision point (xic, yic, zic)
     */
    void change_collision( const moab::CartVect & newXic );

    /**
     * Changes the set of bandwidth values and updates the neighborhood to
     * reflect the new values.
     *
     * @param newH the new set of bandwidth values (hx, hy, hz)
     */
    void change_bandwidth( const moab::CartVect & newH );

    /**
     * Computes the kernel contribution K_xyz for the collision location Xic
     * using the given calculation point X.
     *
     * @param X the calculation point to be evaluated
     * @return the kernel density estimator contribution at the given point
     */
    double compute_contribution( const moab::CartVect & X );

    /**
     * Computes the kernel contribution K_xyz for the collision location Xic
     * using the given calculation point coordinates.
     *
     * @param coords the calculation point to be evaluated
     * @return the kernel density estimator contribution at the given point
     */
    double compute_contribution( const double coords[3] );

  private:
    
    moab::CartVect Xic;   // the collision point (xic, yic, zic)
    moab::CartVect H;     // the set of bandwidth values (hx, hy, hz)
    KDEKernel* kernel;    // the kernel to be used in the computation
    double min[3];        // minimum corner of collision neighborhood
    double max[3];        // maximum corner of collision neighborhood

    /**
     * Sets the neighborhood around the collision point for which a non-zero
     * contribution occurs.
     */
    void set_neighborhood();

};

#endif
