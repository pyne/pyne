from __future__ import print_function

#from pyne cimport cpp_dagmc_bridge
from pyne cimport cpp_dagmc_bridge

cimport numpy as np
import numpy as np

# Python imports
import sys
from contextlib import contextmanager
from copy import copy
from random import uniform

from numpy.linalg import norm

np.import_array()

from pyne.mesh import Mesh

# Set the entity handle type; I don't know a cleaner way than to ask the daglib
# to return the byte width of the type
def get_entity_handle_type():
    eh_size = cpp_dagmc_bridge.dag_ent_handle_size()
    if eh_size == 8:
        eh_t = np.uint64
    elif eh_size == 4:
        eh_t = np.uint32
    else:
        raise TypeError("Unrecognized entity handle size in dagmc library: "
                         + str(eh_size))
    return type("EntityHandle", (eh_t,), {})

EntityHandle = get_entity_handle_type()
_ErrorCode = type("ErrorCode", (np.int,), {})

class DagmcError(Exception):
    pass


def dag_version():
    """Returns the DagMC version."""
    return cpp_dagmc_bridge.dag_version()


def dag_rev_version():
    """Returns the DagMC version number."""
    return int(cpp_dagmc_bridge.dag_rev_version())


def geom_id_list(int dimension):
    """Generates a list of geometry ids.
    """
    cdef int number_of_items, i
    cdef const int * crtn
    if dimension != 2 and dimension != 3:
        raise DagmcError('Incorrect geometric dimension: ' + str(dimension))
    crtn = cpp_dagmc_bridge.geom_id_list(dimension, &number_of_items)
    rtn = [int(crtn[i]) for i in range(number_of_items)]
    return number_of_items, rtn

def handle_from_id(int dimension, int id):
    """Get entity from id number."""
    cdef cpp_dagmc_bridge.EntityHandle crtn
    if dimension != 2 and dimension != 3:
        raise DagmcError('Incorrect geometric dimension: ' + str(dimension))
    crtn = cpp_dagmc_bridge.handle_from_id(dimension, id)
    rtn = EntityHandle(crtn)
    return rtn


def id_from_handle(eh):
    """Get id from entity handle."""
    cdef int eh_id
    if not isinstance(eh, EntityHandle):
        eh = EntityHandle(eh)
    eh_id = cpp_dagmc_bridge.id_from_handle(<cpp_dagmc_bridge.EntityHandle> eh)
    return eh_id

def dag_load(str filename):
    """Loads a file."""
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef const char* cfilename
    bytes_filename = filename.encode('ascii')
    cfilename = bytes_filename
    crtn = cpp_dagmc_bridge.dag_load(cfilename)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    rtn = _ErrorCode(crtn)
    return rtn


cdef class RayHistory(object):
    """A dumb holder for history pointers."""
    cdef void * ptr


def dag_pt_in_vol(vol, np.ndarray[np.float64_t, ndim=1] pt, 
                  np.ndarray[np.float64_t, ndim=1] dir, RayHistory history=None):
    cdef int result
    cdef cpp_dagmc_bridge.ErrorCode crtn
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    if pt.ndim != 1 and pt.shape[0] != 3:
        raise ValueError("pt must have shape=(3,)")
    if dir.ndim != 1 and dir.shape[0] != 3:
        raise ValueError("dir must have shape=(3,)")
    if history is None:
        history = RayHistory()
    crtn = cpp_dagmc_bridge.dag_pt_in_vol(<cpp_dagmc_bridge.EntityHandle> vol, 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(pt), &result, 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(dir), history.ptr)
    return result


def dag_alloc_ray_history():
    """Allocates a new ray history object."""
    cdef RayHistory history = RayHistory()
    history.ptr = cpp_dagmc_bridge.dag_alloc_ray_history()
    return history


def dag_dealloc_ray_history(RayHistory history):
    """Frees an existing ray history object."""
    cpp_dagmc_bridge.dag_dealloc_ray_history(history.ptr)


cdef class RayBuffer(object):
    """A dumb holder for data buffer pointers."""
    cdef void * ptr


def dag_dealloc_ray_buffer(RayBuffer data_buffers):
    """Frees an existing ray buffers object."""
    cpp_dagmc_bridge.dag_dealloc_ray_buffer(data_buffers.ptr)



@contextmanager
def _ray_history():
    history = dag_alloc_ray_history()
    yield history
    dag_dealloc_ray_history(history)


def dag_ray_fire(vol, np.ndarray[np.float64_t, ndim=1] ray_start, 
                 np.ndarray[np.float64_t, ndim=1] ray_dir,
                 RayHistory history=None, double distance_limit=0.0):
    cdef cpp_dagmc_bridge.EntityHandle next_surf 
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef double next_surf_dist = 0.0 
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    if ray_start.ndim != 1 and ray_start.shape[0] != 3:
        raise ValueError("ray_start must have shape=(3,)")
    if ray_dir.ndim != 1 and ray_dir.shape[0] != 3:
        raise ValueError("ray_dir must have shape=(3,)")
    if history is None:
        history = RayHistory()
    crtn = cpp_dagmc_bridge.dag_ray_fire(<cpp_dagmc_bridge.EntityHandle> vol, 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(ray_start), 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(ray_dir),
                &next_surf, &next_surf_dist, history.ptr, distance_limit)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    return EntityHandle(next_surf), next_surf_dist


def dag_ray_follow(firstvol, np.ndarray[np.float64_t, ndim=1] ray_start, 
                   np.ndarray[np.float64_t, ndim=1] ray_dir, double distance_limit, 
                   RayBuffer data_buffers):
    cdef int i
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef int num_intersections = 0
    cdef cpp_dagmc_bridge.EntityHandle * surfs
    cdef double * distances
    cdef cpp_dagmc_bridge.EntityHandle * volumes
    cdef np.npy_intp shape[1]
    shape[0] = 3
    if not isinstance(firstvol, EntityHandle):
        firstvol = EntityHandle(firstvol)
    if ray_start.ndim != 1 and ray_start.shape[0] != 3:
        raise ValueError("ray_start must have shape=(3,)")
    if ray_dir.ndim != 1 and ray_dir.shape[0] != 3:
        raise ValueError("ray_dir must have shape=(3,)")
    crtn = cpp_dagmc_bridge.dag_ray_follow(<cpp_dagmc_bridge.EntityHandle> firstvol, 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(ray_start), 
                <cpp_dagmc_bridge.vec3> np.PyArray_DATA(ray_dir), distance_limit, 
                &num_intersections, &surfs, &distances, &volumes, data_buffers.ptr)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    pysurfs = []
    pydistances = []
    pyvolumes = []
    for i in range(num_intersections):
        pysurfs.append(EntityHandle(surfs[i]))
        pydistances.append(float(distances[i]))
        pyvolumes.append(EntityHandle(volumes[i]))
    return num_intersections, pysurfs, pydistances, pyvolumes


def  dag_next_vol(surface, volume): 
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef cpp_dagmc_bridge.EntityHandle next_vol
    if not isinstance(surface, EntityHandle):
        surface = EntityHandle(surface)
    if not isinstance(volume, EntityHandle):
        volume = EntityHandle(volume)
    crtn = cpp_dagmc_bridge.dag_next_vol(<cpp_dagmc_bridge.EntityHandle> surface, 
                                         <cpp_dagmc_bridge.EntityHandle> volume,
                                         &next_vol)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    return EntityHandle(next_vol)


def vol_is_graveyard(vol):
    """True if the given volume id is a graveyard volume"""
    cdef int crtn 
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    crtn = cpp_dagmc_bridge.vol_is_graveyard(<cpp_dagmc_bridge.EntityHandle> vol)
    print(crtn)
    return bool(crtn)
    

def vol_is_implicit_complement(vol):
    """True if the given volume id is the implicit complement volume"""
    cdef int crtn 
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    crtn = cpp_dagmc_bridge.vol_is_implicit_complement(
                                <cpp_dagmc_bridge.EntityHandle> vol)
    return bool(crtn)


def get_volume_metadata(vol):
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef int material
    cdef double density
    cdef double importance
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    crtn = cpp_dagmc_bridge.get_volume_metadata(<cpp_dagmc_bridge.EntityHandle> vol,
                                                &material, &density, &importance)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    return material, density, importance


def get_volume_boundary(vol):
    cdef int i
    cdef cpp_dagmc_bridge.ErrorCode crtn
    cdef cpp_dagmc_bridge.vec3 minpt  
    cdef cpp_dagmc_bridge.vec3 maxpt 
    cdef np.ndarray pyminpt = np.empty(3, dtype=np.float64)
    cdef np.ndarray pymaxpt = np.empty(3, dtype=np.float64)
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> 3
    if not isinstance(vol, EntityHandle):
        vol = EntityHandle(vol)
    crtn = cpp_dagmc_bridge.get_volume_boundary(<cpp_dagmc_bridge.EntityHandle> vol,
                                                minpt, maxpt)
    if crtn != 0:
        raise DagmcError("Error code " + str(crtn))
    for i in range(3):
        pyminpt[i] = minpt[i]
        pymaxpt[i] = maxpt[i]
    return pyminpt, pymaxpt
    


### end bridge

surf_id_to_handle = {}
surf_handle_to_id = {}
vol_id_to_handle = {}
vol_handle_to_id = {}

def versions():
    """Return a (str, int) tuple: the version and SVN revision of the 
    active DagMC C++ library.
    """
    return ('{0:.4}'.format(dag_version()), int(dag_rev_version()))


def load(filename):
    """Load a given filename into DagMC"""
    global surf_id_to_handle, surf_handle_to_id, vol_id_to_handle, vol_handle_to_id
    dag_load(filename)

    def get_geom_list(dim):
        cdef int count
        count, list_p = geom_id_list(dim)

        r_dict_forward = {}
        for i in range(0, count):
            eh = handle_from_id(dim, list_p[i])
            if eh == 0:
                raise DagmcError('{0} ID {1} has no entity handle'.format({2:'surf',
                                 3:'vol'}[dim], list_p[i]))
            else:
                r_dict_forward[list_p[i]] = eh
        r_dict_backward = dict((v,k) for k,v in r_dict_forward.items())
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
    eh = vol_id_to_handle[vol_id]
    return vol_is_graveyard(eh)


def volume_is_implicit_complement(vol_id):
    """True if the given volume id is the implicit complement volume"""
    eh = vol_id_to_handle[vol_id]
    return vol_is_implicit_complement(eh)


def volume_metadata(vol_id):
    """Get the metadata of the given volume id

    returns a dictionary containing keys 'material', 'rho', and 'imp', corresponding
    to the DagmcVolData struct in DagMC.hpp

    """
    eh = vol_id_to_handle[vol_id]
    mat, rho, imp = get_volume_metadata(eh)
    return {'material': mat, 'rho': rho, 'imp': imp}


def volume_boundary(vol_id):
    """Get the lower and upper boundary of a volume in (x,y,z) coordinates.

    Return the lower and upper coordinates of an axis-aligned bounding box for 
    the given volume.  The returned box may or may not be the minimal bounding 
    box for the volume. Return (xyz low) and (xyz high) as np arrays.
    """
    eh = vol_id_to_handle[vol_id]
    low, high = get_volume_boundary(eh)
    return low, high


def point_in_volume(vol_id, xyz, uvw=[1,0,0]):
    """Determine whether the given point, xyz, is in the given volume.
    
    If provided, uvw is used to determine the ray fire direction for the underlying 
    query.  Otherwise, a random direction will be chosen. 
    
    """
    xyz = np.array(xyz, dtype=np.float64)
    uvw = np.array(uvw, dtype=np.float64)
    eh = vol_id_to_handle[vol_id]
    result = dag_pt_in_vol(eh, xyz, uvw)
    return (result == 1)


def find_volume(xyz, uvw=[1,0,0]):
    """Determine which volume the given point is in.

    Return a volume id.  If no volume contains the point, a DagmcError may be raised,
    or the point may be reported to be part of the implicit complement.

    This function may be slow if many volumes exist.

    """
    xyz = np.array(xyz, dtype=np.float64)
    uvw = np.array(uvw, dtype=np.float64)
    for eh, vol_id in vol_handle_to_id.iteritems():
        result = dag_pt_in_vol(eh, xyz, uvw)
        if result == 1:
            return vol_id    
    raise DagmcError("The point {0} does not appear to be in any volume".format(xyz))


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
    xyz = np.array(xyz, dtype=np.float64)
    uvw = np.array(uvw, dtype=np.float64)
    eh = vol_id_to_handle[vol_id]
    surf_result, dist_result = dag_ray_fire(eh, xyz, uvw)
    if(surf_result != 0):
        return (surf_handle_to_id[surf_result], dist_result)
    else:
        return None


def ray_iterator_slow(init_vol_id, startpoint, direction, **kw):
    """Return an iterator for a ray in a single direction.

    The iterator will yield a series of tuples (vol,dist,surf), indicating the next
    volume intersected, the distance to the next intersection (from the last 
    intersection), and the surface intersected.  Stops iterating when no further 
    intersections are detected along the ray.  This is the only way to traverse 
    volumes along a given ray.

    Keyword arguments:
    yield_xyz: results will contain a fourth tuple element, being the xyz 
               position of the intersection
    dist_limit: distance at which to consider the ray ended
    """
    eh = EntityHandle(vol_id_to_handle[init_vol_id])
    xyz = np.array(startpoint, dtype=np.float64)
    uvw = np.array(direction, dtype=np.float64)

    use_dist_limit = ('dist_limit' in kw)
    dist_limit = kw.get('dist_limit', 0.0)

    with _ray_history() as history:
        while eh != 0:
            if use_dist_limit and dist_limit <= 0:
                break
            surf, dist_result = dag_ray_fire(eh, xyz, uvw, history, dist_limit)
            if surf == 0:
                break

            # eh = the new volume
            eh = dag_next_vol(surf, eh)
            xyz += uvw * dist_result
            if use_dist_limit:
                dist_limit -= dist_result

            newvol = vol_handle_to_id[eh]
            dist = dist_result
            newsurf = surf_handle_to_id[surf]
            
            if kw.get('yield_xyz', False) :
                yield (newvol, dist, newsurf, xyz)
            else: 
                yield (newvol, dist, newsurf)


def ray_iterator(init_vol_id, startpoint, direction, **kw):
    cdef int i, x
    eh = EntityHandle(vol_id_to_handle[init_vol_id])
    xyz = np.array(startpoint, dtype=np.float64)
    uvw = np.array(direction, dtype=np.float64)
    dist_limit = kw.get('dist_limit', 0.0)
    buf = RayBuffer()

    x, surfs, dists, vols = dag_ray_follow(eh, xyz, uvw, dist_limit, buf)
    for i in range(x):
        vol_id = vol_handle_to_id[vols[i]]
        surf_id = surf_handle_to_id[surfs[i]]
        if kw.get('yield_xyz', False):
            xyz += uvw * dists[i]
            yield (vol_id, dists[i], surf_id, xyz)
        else:
            yield (vol_id, dists[i], surf_id)
    
    dag_dealloc_ray_buffer(buf)


def tell_ray_story(startpoint, direction, output=sys.stdout, **kw):
    """Write a human-readable history of a ray in a given direction.

    The history of the ray from startpoint in direction is written 
    to the given output file.
    The initial volume in which startpoint resides will be determined, and 
    the direction argument will be normalized to a unit vector.

    kw args are passed on to underlying call to ray_iterator

    """
    xyz = np.array(startpoint, dtype=np.float64)
    uvw = np.array(direction, dtype=np.float64) 
    uvw /= norm(uvw)

    def pr(*args): 
        print(*args, file=output)

    def vol_notes(v):
        notes = []
        md = volume_metadata(v)
        if md['material'] == 0:
            notes.append('void')
        else:
            notes.append('mat=' + str(md['material']))
            notes.append('rho=' + str(md['rho']))
        if volume_is_graveyard(v): 
            notes.append('graveyard')
        if volume_is_implicit_complement(v):
            notes.append('implicit complement')
        return '({0})'.format(', '.join(notes))

    pr('Starting a ray at', xyz, 'in the direction', uvw)
    
    if 'dist_limit' in kw:
        pr('with a dist_limit of', kw['dist_limit'])

    first_volume = find_volume(xyz, uvw)
    
    pr('The ray starts in volume', first_volume, vol_notes(first_volume))

    kwargs = kw.copy()
    kwargs['yield_xyz'] = True

    for (vol, dist, surf, xyz) in ray_iterator(first_volume, xyz, uvw, **kwargs):

        pr('  next intersection at distance', dist, 'on surface', surf)
        pr('  new xyz =', xyz)
        pr('proceeding into volume', vol, vol_notes(vol))

    pr('No more intersections')

#### start util

def find_graveyard_inner_box():
    """Estimate the dimension of the inner wall of the graveyard, assuming box shape.

    Return the the (low, high) xyz coordinates of the inner wall of the 
    graveyard volume.
    This assumes the graveyard volume is a box-shaped, axis-aligned shell.  
    This function is slow and should not be called many times where performance 
    is desirable.
    """
    cdef int i
    volumes = get_volume_list()
    graveyard = 0
    for v in volumes:
        if volume_is_graveyard(v): 
            graveyard = v
            break
    if graveyard == 0:
        raise DagmcError('Could not find a graveyard volume')

    xyz_lo, xyz_hi = volume_boundary(graveyard)
    xyz_mid = np.array([(hi + lo) / 2.0 for (hi, lo) in zip(xyz_hi, xyz_lo)], 
                       dtype=np.float64)

    result_lo = np.array([0] * 3, dtype=np.float64 )
    result_hi = np.array([0] * 3, dtype=np.float64 )

    for i in range(0, 3):
        uvw = [0, 0, 0]
        uvw[i] = 1
        lo_mid = xyz_mid.copy()
        lo_mid[i] = xyz_lo[i]
        _, dist = fire_one_ray(graveyard, lo_mid, uvw)
        result_lo[i] = lo_mid[i] + dist
        uvw[i] = -1
        hi_mid = xyz_mid.copy()
        hi_mid[i] = xyz_hi[i]
        _, dist = fire_one_ray(graveyard, hi_mid, uvw)
        result_hi[i] = hi_mid[i] - dist
   
    return result_lo, result_hi

def get_material_set(**kw):
    """Return all material IDs used in the geometry as a set of integers
    
    If the keyword argument 'with_rho' is True, the set will contain (int, float) 
    tuples containing material ID and density
    """
    mat_ids = set()
    volumes = get_volume_list()
    for v in volumes:
        d = volume_metadata(v)
        if(kw.get('with_rho') is True):
            # rho is undefined for the void material and dagmc may return anything.
            if d['material'] == 0:
                mat_ids.add((d['material'], 0.0))
            else:
                mat_ids.add((d['material'], d['rho']))
        else:
            mat_ids.add(d['material'])
    return mat_ids

#### start util
def discretize_geom(mesh, num_rays, grid=False):
    """This function reads in a structured, axis-aligned, PyNE mesh object,
    then uses Monte Carlo ray tracing to determing the volume fraction of each
    geometry volume within each mesh volume element of the mesh. Note that a 
    DAGMC geometry must already be loaded into memory.
 
    Parameters
    ----------
    mesh : PyNE Mesh object
        The Mesh object to tag with cell volume fractions. The mesh must be a 
        Cartesean structured mesh that overlays the geometry.
    num_rays : int
        The number of rays to fire in each mesh row for each direction.
    grid : boolean
        If false, rays starting points are chosen randomly (on the boundary) 
        for each mesh row. If true, a linear spaced grid of starting points is 
        chosen, with dimension sqrt(num_rays) x sqrt(num_rays)

    Returns
    -------
    results : list
        List with an entry for each mesh volume element in the order determined
        by the structured_ordering attribute of the PyNE mesh object. Each entry
        in the list is a dictionary maps geometry volume numbers to its volume 
        fraction within the mesh volume element.
    uncs : list
        Same structure as "results" but contains the standard error for the 
        volume fractions.
     """
    # Ensure input is valid
    mesh._structured_check()
    #load(filename)
    divs = [mesh.structured_get_divisions(a) for a in 'xyz']
    num_ves = (len(divs[0]) - 1) * (len(divs[1]) - 1 ) * (len(divs[2]) - 1)
    results = [{} for x in range(0, num_ves)]
    uncs = [{} for x in range(0, num_ves)]

    # direction indicies: x = 0, y = 1, z = 2
    dis = [0, 1, 2]
    # Iterate over all directions indicies
    for di in dis:
        # For each direction, the remaining two directions define the sampling 
        # surface. These two directions are the values in s_dis (surface
        # direction indices)
        s_dis = copy(dis)
        s_dis.remove(di)

        # iterate through all all the sampling planes perpendicular to di,
        # creating a _MeshRow in each, and subsequently evaluating that row.
        for a in range(0, len(divs[s_dis[0]]) - 1):
            for b in range(0, len(divs[s_dis[1]]) - 1):
                s_min_0 = divs[s_dis[0]][a]
                s_max_0 = divs[s_dis[0]][a + 1]
                s_min_1 = divs[s_dis[1]][b]
                s_max_1 = divs[s_dis[1]][b + 1]

                #print("direction {0}".format(di))
                #print("bounds {0}, {1}, {2}, {3}".format(s_min_0, s_max_0, s_min_1, s_max_1))

                # create a lines of starting points to fire rays for this
                # particular mesh row
                if not grid:
                    start_points = _rand_start(num_rays, di, divs[di][0],
                                               s_dis[0], s_min_0, s_max_0, 
                                               s_dis[1], s_min_1, s_max_1)
                else:
                    start_points = _grid_start(num_rays, di, divs[di][0],
                                               s_dis[0], s_min_0, s_max_0, 
                                               s_dis[1], s_min_1, s_max_1)
                #print(start_points)
                
                row_samples = _evaluate_row(di, divs[di], start_points)
                #for res in row_results:
                #    print(res)

                if di == 0:
                    ves = mesh.structured_iterate_hex('x', y=a, z=b)
                elif di == 1:
                    ves = mesh.structured_iterate_hex('y', x=a, z=b)
                elif di == 2:
                    ves = mesh.structured_iterate_hex('z', x=a, y=b)

                idx_tag = mesh.mesh.getTagHandle("idx")
                idx = []
                for ve in ves:
                    #print(idx_tag[ve])
                    idx.append(idx_tag[ve])

                # Calculate means. Simotaneous populate uncs with a list of all
                # samples so the standard error can later be calculated.
                for i, samples in zip(idx, row_samples):
                    for vol, val in samples.iteritems():
                        if vol not in results[i].keys():
                           results[i][vol] = 0
                           uncs[i][vol] = []

                        results[i][vol] += np.sum(samples[vol])/(3 * num_rays)
                        uncs[i][vol] += samples[vol]

    # calculate standard errors
    for i in range(0, num_ves):
        for vol in uncs[i].keys():
            # calculate sample standard deviation
            var = (np.sum([(results[i][vol] - x)**2 for x in uncs[i][vol]])
                  /(len(uncs[i][vol]) - 1))
            stdev = np.sqrt(var)
            uncs[i][vol] = stdev/np.sqrt(len(uncs[i][vol]))

    return results, uncs
                        

def _evaluate_row(di, divs, start_points):
    """Private function that fires rays down a single mesh row and returns the
    results.

    Parameters
    ----------
    di : int
        The direction index of the current firing direction.
    divs : list
        The mesh boundaries perpendicular to the firing direction.
    start_points : list
        The xyz points describing where rays should start from.

    Returns
    -------
    samples : list
        One entry for each mesh volume element in the firing direction. Each
        entry is a dictionary that maps geometry volume numbers to a list of
        normalized track length samples (ratio of track length to mesh volume 
        element width).
    """

    # Total number of rays fired: multiply by 3 to account for 3 directions
    num_rays = len(start_points) * 3
    # number of volume elements in this mesh row
    num_ve = len(divs) - 1
    direction = [0, 0, 0]
    direction[di] = 1
    samples = [{} for x in range(0, len(divs) - 1)]
    width = [divs[x] - divs[x - 1] for x in range(1, len(divs))]
    # fire ray for each starting point
    for point in start_points:
        vol = find_volume(point, direction)
        mesh_dist = width[0]
        ve_count = 0
        complete = False
        #print("\nStarting at point {0}, in vol {1}".format(point, vol))
        # track a single ray down the mesh row and tally accordingly
        for next_vol, distance, _ in ray_iterator(vol, point, direction):
            if complete:
                break
            #print("next_vol {0} distance {1}".format(next_vol, distance))
            # volume extends past mesh boundary
            while distance >= mesh_dist:
                # check to see if current volume has already by tallied
                if vol not in samples[ve_count].keys():
                    #print("hello {0}".format(distance))
                    samples[ve_count][vol] = []

                #print("hello Travel {0} in ve_count {1}".format(mesh_dist, ve_count))
                samples[ve_count][vol].append(mesh_dist/width[ve_count])
                distance -= mesh_dist

                # if not on the last volume element, continue into the next
                # volume element
                if ve_count == num_ve - 1:
                    complete = True
                    break
                else:
                    ve_count += 1
                    mesh_dist = width[ve_count]

            # volume does not extend past mesh volume
            if distance < mesh_dist and distance > 1E-10 and not complete:
                # check to see if current volume has already by tallied
                if vol not in samples[ve_count].keys():
                    #print("goob {0}".format(distance))
                    samples[ve_count][vol] = []

                #print("goob Travel {0} in ve_count {1}".format(distance, ve_count))
                samples[ve_count][vol].append(distance/width[ve_count])
                mesh_dist -= distance
            
            vol = next_vol

    return samples

def _rand_start(num_rays, di_fire, min_fire, di_1, min_1, max_1, 
                                             di_2, min_2, max_2):
    """Private function for randomly generating ray starting points for a single
    mesh row.
    
    Parameters
    ----------
    num_rays : int
        Number of starting points to generate.
    di_fire : int
        The direction index of ray fire.
    min_fire : float
        The location of fire plane perpedicular to the firing direction, i.e.
        the plane for which all starting locations are on.
    di_1 : int
        The direction index of one of the directions that makes up the rectangle
        for which starting points are sampled from.
    min_1 : float
        The location that bounds the sampling rectangle from the left in di_1.
    max_1 : float
        The location that bounds the sampling rectangle from the right in di_1.
    di_2 : int
        The direction index of one the directions (that is not di_1) that makes
        up the rectangle for which starting points are sampled from.
    min_2 : float
        The location that bounds the sampling rectangle from the left in di_2.
    max_2 : float
        The location that bounds the sampling rectangle from the right in di_2.
  
    Returns
    -------
    start_points : list
        The xyz values of the randomly sampled starting points.
    """
    start_points = []
    ray_count = 0
    while ray_count < num_rays:
        start_point = [0]*3
        start_point[di_fire] = min_fire
        start_point[di_1] = uniform(min_1, max_1)
        start_point[di_2] = uniform(min_2, max_2)
        start_points.append(start_point)
        ray_count += 1

    return start_points

def _grid_start(num_rays, di_fire, min_fire, di_1, min_1, max_1,
                                             di_2, min_2, max_2):
    """Private function for generating a uniform grid of start points for a
    sinlge mesh row. The grid is a rectangle, each side divided equally into 
    sqrt(num_rays) rays. The boundaries (values specified by min_1, min_2,
    max_1, max_2) are not included in the sampling grid.
    
    Parameters
    ----------
    num_rays : int
        Number of starting points to generate. Must be a perfect square.
    di_fire : int
        The direction index of ray fire.
    min_fire : float
        The location of fire plane perpedicular to the firing direction, i.e.
        the plane for which all starting locations are on.
    di_1 : int
        The direction index of one of the directions that makes up the rectangle
        for which starting points are sampled from.
    min_1 : float
        The location that bounds the sampling rectangle from the left in di_1.
    max_1 : float
        The location that bounds the sampling rectangle from the right in di_1.
    di_2 : int
        The direction index of one the directions (that is not di_1) that makes
        up the rectangle for which starting points are sampled from.
    min_2 : float
        The location that bounds the sampling rectangle from the left in di_2.
    max_2 : float
        The location that bounds the sampling rectangle from the right in di_2.

    Returns
    -------
    start_points : list
        The xyz values of the starting points.
    """
    # test to see if num_rays is a perfect square
    if int(np.sqrt(num_rays))**2 != num_rays:
        raise ValueError("For rays fired in a grid, "
                         "num_rays must be a perfect square.")
    else:
       square_dim = int(np.sqrt(num_rays))
 
    step_1 = (max_1 - min_1)/(float(square_dim) + 1)
    step_2 = (max_2 - min_2)/(float(square_dim) + 1)
    range_1 = np.linspace(min_1 + step_1, max_1, square_dim, endpoint=False)
    range_2 = np.linspace(min_2 + step_2, max_2, square_dim, endpoint=False)
 
    start_points = []
    for point_1 in range_1:
        for point_2 in range_2:
            start_point = [0]*3
            start_point[di_fire] = min_fire
            start_point[di_1] = point_1
            start_point[di_2] = point_2
            start_points.append(start_point)
 
    return start_points
