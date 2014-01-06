import numpy

from dagmc import *
from bridge import DagmcError

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

def get_material_set(**kw):
    """Return all material IDs used in the geometry as a set of integers
    
    If the keyword argument 'with_rho' is True, the set will contain (int,float) tuples
    containing material ID and density
    """
    mat_ids = set()
    volumes = get_volume_list()
    for v in volumes:
        d = volume_metadata( v )
        if( kw.get('with_rho') is True ):
            # rho is undefined for the void material and dagmc may return anything.
            if d['material'] == 0:
                mat_ids.add( (d['material'], 0.0) )
            else:
                mat_ids.add( (d['material'], d['rho']) )
        else:
            mat_ids.add( d['material'] )
    return mat_ids
