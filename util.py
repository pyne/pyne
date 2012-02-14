import numpy

from dagmc import *

def find_graveyard_inner_box():
    """Estimate the dimension of the inner wall of the graveyard, assuming box shape.

    Return the the (low, high) xyz coordinates of the inner wall of the graveyard volume.
    This assumes the graveyard volume is a box-shaped, axis-aligned shell.  This function
    is slow and should not be called many times where performance is desirable.
    """
    volumes = get_volume_list()
    graveyard = 0
    for v in volumes:
        if volume_is_graveyard( v ): 
            graveyard = v
            break
    if graveyard == 0:
        raise DagmcError( 'Could not find a graveyard volume' )

    xyz_lo, xyz_hi = volume_boundary( graveyard )
    xyz_mid = numpy.array( [ (hi+lo)/2.0 for (hi,lo) in zip( xyz_hi, xyz_lo) ], dtype=numpy.float64 )

    result_lo = numpy.array( [0]*3, dtype=numpy.float64 )
    result_hi = numpy.array( [0]*3, dtype=numpy.float64 )

    for i in range(0,3):
        uvw = [0,0,0]
        uvw[i] = 1
        lo_mid = xyz_mid.copy()
        lo_mid[i] = xyz_lo[i]
        _, dist = fire_one_ray( graveyard, lo_mid, uvw )
        result_lo[i] = lo_mid[i] + dist
        uvw[i] = -1
        hi_mid = xyz_mid.copy()
        hi_mid[i] = xyz_hi[i]
        _, dist = fire_one_ray( graveyard, hi_mid, uvw )
        result_hi[i] = hi_mid[i] - dist
   
    return result_lo, result_hi
