from __future__ import print_function

#### start __init__

import ctypes
import ctypes.util
import sys
import os.path

# Attempt to find the dagmc_bridge shared library
# ctypes.util.find_library is (sometimes) useful, but does not work the
# same way on all platforms.
_dirname = os.path.dirname(__file__)
if sys.platform == 'darwin':
    _bridgepath = ctypes.util.find_library(
                    os.path.join(_dirname, 'dagmc_bridge', 'dagmc_bridge'))
elif sys.platform.startswith('linux'):
    _bridgepath = os.path.join(_dirname, 'dagmc_bridge', 'dagmc_bridge.so')
else:
    print('Warning: unknown platform, dagmc_bridge lib import may fail',
          file=sys.stderr)
    _bridgepath = ctypes.util.find_library('dagmc_bridge')


daglib = ctypes.CDLL(_bridgepath)


# Set the entity handle type; I don't know a cleaner way than to ask the daglib
# to return the byte width of the type
def get_EH_type():
    EH_size = daglib.dag_ent_handle_size()
    if EH_size == 8:
        EH_t = ctypes.c_uint64
    elif EH_size == 4:
        EH_t = ctypes.c_uint32
    else:
        raise TypeError("Unrecognized entity handle size in dagmc library: "
                         + str(EH_size))
    return type("EntityHandle", (EH_t,), {})


EntityHandle = get_EH_type()
_ErrorCode = type("ErrorCode", (ctypes.c_int,), {})


class DagmcError(Exception):
    pass


#### end __init__


#### start bridge

import ctypes
from contextlib import contextmanager
from numpy import float64
from numpy.ctypeslib import ndpointer

from pydagmc import daglib, EntityHandle, _ErrorCode, DagmcError

lib = daglib

lib.dag_version.restype = ctypes.c_float
lib.dag_version.argtypes = []

lib.dag_rev_version.restype = ctypes.c_uint
lib.dag_rev_version.argtypes = []


def _geom_dim_check(result, func, arguments):
    if arguments[0] not in (2, 3):
        raise DagmcError('Incorrect geometric dimension: ' + str(arguments[0]))
    return result

lib.geom_id_list.restype = ctypes.POINTER(ctypes.c_int)
lib.geom_id_list.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
lib.geom_id_list.errcheck = _geom_dim_check

lib.handle_from_id.restype = EntityHandle
lib.handle_from_id.argtypes = [ctypes.c_int, ctypes.c_int]
lib.handle_from_id.errcheck = _geom_dim_check

lib.id_from_handle.restype = ctypes.c_int
lib.id_from_handle.argtypes = [EntityHandle]


def _moab_error_check(result, func, arguments):
    if result.value != 0:
        raise DagmcError("Error code " + str(result.value) +
                         " returned from " + func.__name__)
    return result


def _returns_moab_errors(function):
    function.restype = _ErrorCode
    function.errcheck = _moab_error_check

_returns_moab_errors(lib.dag_load)
lib.dag_load.argtypes = [ctypes.c_char_p]

_vec3 = ndpointer(dtype=float64, shape=(3,), flags='CONTIGUOUS')

_returns_moab_errors(lib.dag_pt_in_vol)
lib.dag_pt_in_vol.argtypes = [EntityHandle, _vec3,
                              ctypes.POINTER(ctypes.c_int),
                              _vec3, ctypes.c_void_p]

lib.dag_alloc_ray_history.restype = ctypes.c_void_p
lib.dag_alloc_ray_history.argtypes = []

lib.dag_dealloc_ray_history.restype = None
lib.dag_dealloc_ray_history.argtypes = [ctypes.c_void_p]

lib.dag_dealloc_ray_buffer.restype = None
lib.dag_dealloc_ray_buffer.argtypes = [ctypes.c_void_p]

@contextmanager
def _ray_history():
    history = lib.dag_alloc_ray_history()
    yield history
    lib.dag_dealloc_ray_history(history)

_returns_moab_errors(lib.dag_ray_fire)
lib.dag_ray_fire.argtypes = [EntityHandle, _vec3, _vec3,
                             ctypes.POINTER(EntityHandle),
                             ctypes.POINTER(ctypes.c_double),
                             ctypes.c_void_p, ctypes.c_double]

_returns_moab_errors(lib.dag_ray_follow)
lib.dag_ray_follow.argtypes = [EntityHandle, _vec3, _vec3, ctypes.c_double,
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.POINTER(EntityHandle)), 
                               ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
                               ctypes.POINTER(ctypes.POINTER(EntityHandle)),
                               ctypes.c_void_p]
                               

_returns_moab_errors(lib.dag_next_vol)
lib.dag_next_vol.argtypes = [EntityHandle, EntityHandle,
                             ctypes.POINTER(EntityHandle)]

lib.vol_is_graveyard.restype = ctypes.c_int
lib.vol_is_graveyard.argtypes = [EntityHandle]

lib.vol_is_implicit_complement.restype = ctypes.c_int
lib.vol_is_implicit_complement.argtypes = [EntityHandle]

_returns_moab_errors(lib.get_volume_metadata)
lib.get_volume_metadata.argtypes = [EntityHandle, ctypes.POINTER(ctypes.c_int),
                                    ctypes.POINTER(ctypes.c_double),
                                    ctypes.POINTER(ctypes.c_double)]

_returns_moab_errors(lib.get_volume_boundary)
lib.get_volume_boundary.argtypes = [EntityHandle, _vec3, _vec3]

### end bridge

import sys
import ctypes
import numpy
from numpy.linalg import norm

surf_id_to_handle = {}
surf_handle_to_id = {}
vol_id_to_handle = {}
vol_handle_to_id = {}


def versions():
    """Return a (str,int) tuple: the version and SVN revision of the 
    active DagMC C++ library.
    """
    return ('{0:.4}'.format(bridge.lib.dag_version()), int(bridge.lib.dag_rev_version()))


def load(filename):
    """Load a given filename into DagMC"""
    global surf_id_to_handle, surf_handle_to_id, vol_id_to_handle, vol_handle_to_id
    bridge.lib.dag_load(filename)

    def get_geom_list(dim):
        
        count = ctypes.c_int(0)
        
        list_p = bridge.lib.geom_id_list(dim, ctypes.byref(count))

        r_dict_forward = {}
        for i in range(0,count.value):
            eh = bridge.lib.handle_from_id(dim, list_p[i]).value 
            if eh == 0:
                raise bridge.DagmcError(
                        '{0} ID {1} has no entity handle'.format(
                         {2:'surf',3:'vol'}[dim], list_p[i]))
            else:
                r_dict_forward[ list_p[i] ] = eh

        r_dict_backward = dict((v,k) for k,v in r_dict_forward.iteritems())
        
        return r_dict_forward, r_dict_backward

    surf_id_to_handle, surf_handle_to_id = get_geom_list(2)
    vol_id_to_handle, vol_handle_to_id  = get_geom_list(3)


def get_surface_list():
    """return a list of valid surface IDs"""
    return surf_id_to_handle.keys()


def get_volume_list():
    """return a list of valid volume IDs"""
    return vol_id_to_handle.keys()


def volume_is_graveyard(vol_id):
    """True if the given volume id is a graveyard volume"""
    eh = vol_id_to_handle[ vol_id ]
    result = bridge.lib.vol_is_graveyard(eh)
    return (result != 0)


def volume_is_implicit_complement(vol_id):
    """True if the given volume id is the implicit complement volume"""
    eh = vol_id_to_handle[ vol_id ]
    result = bridge.lib.vol_is_implicit_complement(eh)
    return (result != 0)


def volume_metadata(vol_id):
    """Get the metadata of the given volume id

    returns a dictionary containing keys 'material', 'rho', and 'imp', corresponding
    to the DagmcVolData struct in DagMC.hpp

    """
    eh = vol_id_to_handle[ vol_id ]
    mat = ctypes.c_int()
    rho = ctypes.c_double()
    imp = ctypes.c_double()

    bridge.lib.get_volume_metadata(eh, mat, rho, imp)

    return {'material':mat.value, 'rho':rho.value, 'imp':imp.value}


def volume_boundary(vol_id):
    """Get the lower and upper boundary of a volume in (x,y,z) coordinates.

    Return the lower and upper coordinates of an axis-aligned bounding box for the given
    volume.  The returned box may or may not be the minimal bounding box for the volume.
    Return (xyz low) and (xyz high) as numpy arrays.
    """
    eh = vol_id_to_handle[ vol_id ]
    low = numpy.array([0,0,0], dtype=numpy.float64)
    high = numpy.array([0,0,0], dtype=numpy.float64)

    bridge.lib.get_volume_boundary(eh, low, high)
    return low, high


def point_in_volume(vol_id, xyz, uvw=[1,0,0]):
    """Determine whether the given point, xyz, is in the given volume.
    
    If provided, uvw is used to determine the ray fire direction for the underlying 
    query.  Otherwise, a random direction will be chosen. 
    
    """
    xyz = numpy.array(xyz, dtype=numpy.float64)
    uvw = numpy.array(uvw, dtype=numpy.float64)

    eh = vol_id_to_handle[ vol_id ]
    result = ctypes.c_int(-2)

    bridge.lib.dag_pt_in_vol(eh, xyz, ctypes.byref(result), uvw, None)

    return (result.value == 1)


def find_volume(xyz, uvw=[1,0,0]):
    """Determine which volume the given point is in.

    Return a volume id.  If no volume contains the point, a DagmcError may be raised,
    or the point may be reported to be part of the implicit complement.

    This function may be slow if many volumes exist.

    """
    xyz = numpy.array(xyz, dtype=numpy.float64)
    uvw = numpy.array(uvw, dtype=numpy.float64)

    for eh, vol_id in vol_handle_to_id.iteritems():
        result = ctypes.c_int(-2)
        bridge.lib.dag_pt_in_vol(eh, xyz, ctypes.byref(result), uvw, None)
        if result.value == 1:
            return vol_id
    
    raise bridge.DagmcError("The point {0} does not appear to be in any volume".format(xyz))


def fire_one_ray(vol_id, xyz, uvw):
    """Fire a ray from xyz, in the direction uvw, at the specified volume

    uvw must represent a unit vector.

    Only intersections that *exit* the volume will be detected.  Entrance intersections
    are not detected.  In most cases, you should only 
    call this function with arguments for which point_in_volume would return True.

    Returns a (surface id, distance) tuple, or None if no intersection detected.

    If a ray in a given direction will traverse several volumes in a row, ray_iterator should
    be used instead.
    """
    xyz = numpy.array(xyz, dtype=numpy.float64)
    uvw = numpy.array(uvw, dtype=numpy.float64)

    eh = vol_id_to_handle[ vol_id ]
    
    surf_result = bridge.EntityHandle(0)
    dist_result = ctypes.c_double(0.0)

    bridge.lib.dag_ray_fire(eh, xyz, uvw, 
                             ctypes.byref(surf_result), ctypes.byref(dist_result),
                             None, 0.0)

    if(surf_result.value != 0):
        return (surf_handle_to_id[ surf_result.value ], dist_result.value)
    else:
        return None


def ray_iterator_slow(init_vol_id, startpoint, direction, **kw):
    """Return an iterator for a ray in a single direction.

    The iterator will yield a series of tuples (vol,dist,surf), indicating the next
    volume intersected, the distance to the next intersection (from the last intersection),
    and the surface intersected.  Stops iterating when no further intersections are 
    detected along the ray.  This is the only way to traverse volumes along a given ray.

    Keyword arguments:
    yield_xyz: results will contain a fourth tuple element, being the xyz position of the 
               intersection
    dist_limit: distance at which to consider the ray ended
    """

    eh = bridge.EntityHandle(vol_id_to_handle[ init_vol_id ])
    xyz = numpy.array(startpoint, dtype=numpy.float64)
    uvw = numpy.array(direction, dtype=numpy.float64)

    use_dist_limit = ('dist_limit' in kw)
    dist_limit = kw.get('dist_limit',0.0)

    with bridge._ray_history() as history:
        
        surf = bridge.EntityHandle(0)
        dist_result = ctypes.c_double(0.0)
        while eh != 0:

            if use_dist_limit and dist_limit <= 0:
                break

            bridge.lib.dag_ray_fire(eh, xyz, uvw, 
                                     ctypes.byref(surf), ctypes.byref(dist_result),
                                     history, dist_limit)

            if surf.value == 0:
                break

            # eh = the new volume
            bridge.lib.dag_next_vol(surf, eh, ctypes.byref(eh))
            xyz += uvw * dist_result.value
            if use_dist_limit:
                dist_limit -= dist_result.value

            newvol = vol_handle_to_id[eh.value]
            dist = dist_result.value
            newsurf = surf_handle_to_id[ surf.value ]
            
            if kw.get('yield_xyz', False) :
                yield (newvol, dist, newsurf, xyz)
            else: 
                yield (newvol, dist, newsurf)

def ray_iterator(init_vol_id, startpoint, direction, **kw):

    eh = bridge.EntityHandle(vol_id_to_handle[ init_vol_id ])
    xyz = numpy.array(startpoint, dtype=numpy.float64)
    uvw = numpy.array(direction, dtype=numpy.float64)

    dist_limit = kw.get('dist_limit',0.0)

    surfs = ctypes.POINTER(bridge.EntityHandle)()
    vols = ctypes.POINTER(bridge.EntityHandle)()
    dists = ctypes.POINTER(ctypes.c_double)()
    buf = ctypes.c_void_p()
    x = ctypes.c_int(0)

    bridge.lib.dag_ray_follow(eh, xyz, uvw, dist_limit, 
                               ctypes.byref(x), ctypes.byref(surfs), 
                               ctypes.byref(dists), ctypes.byref(vols), buf)
    for i in range(x.value):
        vol_id = vol_handle_to_id[vols[i].value]
        surf_id = surf_handle_to_id[surfs[i].value]
        if kw.get('yield_xyz', False):
            xyz += uvw * dists[i]
            yield (vol_id, dists[i], surf_id, xyz)
        else:
            yield (vol_id, dists[i], surf_id)
    
    bridge.lib.dag_dealloc_ray_buffer(buf)

def tell_ray_story(startpoint, direction, output=sys.stdout, **kw):
    """Write a human-readable history of a ray in a given direction.

    The history of the ray from startpoint in direction is written to the given output file.
    The initial volume in which startpoint resides will be determined, and 
    the direction argument will be normalized to a unit vector.

    kw args are passed on to underlying call to ray_iterator

    """
    xyz = numpy.array(startpoint, dtype=numpy.float64)
    uvw = numpy.array(direction, dtype=numpy.float64) 
    uvw /= norm(uvw)

    def pr(*args): 
        print(*args, file=output)

    def vol_notes(v):
        notes = []
        md = volume_metadata(v)
        if md['material'] == 0:
            notes.append('void')
        else:
            notes.append('mat='+str(md['material']))
            notes.append('rho='+str(md['rho']))
        if volume_is_graveyard(v): 
            notes.append('graveyard')
        if volume_is_implicit_complement(v):
            notes.append('implicit complement')
        return '({0})'.format(', '.join(notes))

    pr('Starting a ray at',xyz,'in the direction',uvw)
    
    if 'dist_limit' in kw:
        pr('with a dist_limit of', kw['dist_limit'])

    first_volume = find_volume(xyz, uvw)
    
    pr('The ray starts in volume', first_volume, vol_notes(first_volume))

    kwargs = kw.copy()
    kwargs['yield_xyz'] = True

    for (vol, dist, surf, xyz) in ray_iterator(first_volume, xyz, uvw, **kwargs):

        pr('  next intersection at distance',dist,'on surface',surf)
        pr('  new xyz =', xyz)
        pr('proceeding into volume', vol, vol_notes(vol))

    pr('No more intersections')

#### start util

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

#### start util
