from __future__ import print_function, division, unicode_literals

# Python imports
import sys
from contextlib import contextmanager
from warnings import warn
from pyne.utils import QA_warn

cimport numpy as np
import numpy as np

from pyne cimport cpp_dagmc_bridge
from pyne.mesh import Mesh
from numpy.linalg import norm
from pyne.material import Material
from pyne.material_library import MaterialLibrary

np.import_array()

QA_warn(__name__)

if sys.version_info[0] >= 3:
    unichr = chr

# Mesh specific imports
from pyne.mesh import HAVE_PYMOAB
from pymoab import core as mb_core, types

if not HAVE_PYMOAB:
    warn("The PyMOAB optional dependency could not be imported. "
         "Some aspects of dagmc module may be incomplete",
         ImportWarning)

# Globals
VOL_FRAC_TOLERANCE = 1E-10 # The maximum volume fraction to be considered valid

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
    return type(str("EntityHandle"), (eh_t,), {})

EntityHandle = get_entity_handle_type()
_ErrorCode = type(str("ErrorCode"), (np.int,), {})

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
    return list(surf_id_to_handle.keys())

def get_volume_list():
    """return a list of valid volume IDs"""
    return list(vol_id_to_handle.keys())


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


def cell_material_assignments(hdf5):
    """Get dictionary of cell to material assignments

    Parameters:
    -----------
    hdf5 : string
        Path to hdf5 material-laden geometry

    Returns:
    --------
    mat_assigns : dict
        Dictionary of the cell to material assignments. Keys are cell
        numbers and values are material names
    """
    # Load the geometry as a pymoab instance
    dag_geom = mb_core.Core()
    dag_geom.load_file(hdf5)

    # Get tag handle
    cat_tag = dag_geom.tag_get_handle(types.CATEGORY_TAG_NAME)
    id_tag = dag_geom.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
    name_tag = dag_geom.tag_get_handle(types.NAME_TAG_NAME)

    # Get list of materials and list of cells
    mat_assigns={}

    # Assign the implicit complement to vacuum
    # NOTE: This is a temporary work-around and it is just assumed that there
    # is no material already assigned to the implicit complement volume.
    implicit_vol = find_implicit_complement()
    mat_assigns[implicit_vol] = "mat:Vacuum"

   # Get all meshsets that topologically represent a "Group"
    group_meshsets = dag_geom.get_entities_by_type_and_tag(0, types.MBENTITYSET,
                                                           [cat_tag,], ["Group",])

    # loop over all group_sets in model
    for group_set in group_meshsets:
        group_members = dag_geom.get_entities_by_handle(group_set)
        name = dag_geom.tag_get_data(name_tag, group_set, flat = True)[0]
        # if group_set is a group with a material name_tag, loop over group
        # members and assign name to cell
        # in general `name` may be a non-string or string-like type,
        # so convert to str for safety
        if 'mat:' in str(name):
            for group_member in group_members:
                try:
                    cell = dag_geom.tag_get_data(id_tag, group_member, flat = True)[0]
                    mat_assigns[cell] = name
                except:
                    pass

    return mat_assigns

def cell_materials(hdf5, **kwargs):
    """Obtain a material object for each cell in a DAGMC material-laden
    geometry, tagged in UWUW format [1], i.e. "mat:<name>/rho:<density>" or
    "mat:<name>".

    Parameters:
    -----------
    hdf5 : string
        Path to hdf5 material-laden geometry
    datapath: str, optional, default ='/materials',
        The path in the heirarchy to the material data table in the HDF5 file.

    Returns:
    --------
    cell_mats : dict
        Dictionary that maps cells numbers to PyNE Material objects.

    [1] http://svalinn.github.io/DAGMC/usersguide/uw2.html
    """
    datapath = kwargs.get('datapath', '/materials')

    # void material
    void_mat = Material({}, density = 0.0, metadata={'name': 'void',
                                                      'mat_number': 0})
    # strings that specify that a region is void
    void_names = ['vacuum', 'graveyard', 'void']

    ml = MaterialLibrary()
    ml.from_hdf5(hdf5, datapath=datapath)
    mat_assigns = cell_material_assignments(hdf5)
    cell_mats = {}
    for cell_num, mat_name in mat_assigns.items():
        if cell_num is None:
            continue
        elif np.any([x in mat_name.lower() for x in void_names]):
            cell_mats[cell_num] = void_mat
        else:
            cell_mats[cell_num] = ml[mat_name]

    return cell_mats

def find_implicit_complement():
    """Find the implicit complement and return the volume id.
    Note that a DAGMC geometry must already be loaded into memory.
    """
    volumes = get_volume_list()
    for vol in volumes:
        if volume_is_implicit_complement(vol):
            return vol


def _tag_to_string(tag):
    """Convert ascii to string
    """
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii and join to string
            a.append(str(unichr(part)))
            string = ''.join(a)
    return string


#### start util
def discretize_geom(mesh, **kwargs):
    """discretize_geom(mesh, **kwargs)
    This function discretizes a geometry (by geometry cell) onto a superimposed
    mesh. If the mesh is structured, ray_discretize() is called and Monte Carlo
    ray tracing is used to determine the volume fractions of each geometry cell
    within each mesh volume element of the mesh. If the mesh is not structured,
    cell_at_ve_centers is called and mesh volume elements are assigned a
    geometry cell based off of what geometry cell occupies the center of the
    mesh volume element. The output of cell_at_ve_centers is then put in the
    same structured array format used by ray_discretize(). Note that a DAGMC
    geometry must already be loaded into memory.

    Parameters
    ----------
    mesh : PyNE Mesh
        A Cartesian, structured, axis-aligned Mesh that superimposed the
        geometry.
    num_rays : int, optional, default = 10
        Structured mesh only. The number of rays to fire in each mesh row for
        each direction.
    grid : boolean, optional, default = False
        Structured mesh only. If false, rays starting points are chosen randomly
        (on the boundary) for each mesh row. If true, a linearly spaced grid of
        starting points is used, with dimension sqrt(num_rays) x sqrt(num_rays).
        In this case, "num_rays" must be a perfect square.

    Returns
    -------
    results : structured array
        Stores in a one dimensional array, each entry containing the following
        fields:
        :idx: int
            The volume element index.
        :cell: int
            The geometry cell number.
        :vol_frac: float
            The volume fraction of the cell withing the mesh ve.
        :rel_error: float
            The relative error associated with the volume fraction.
        This array is returned in sorted order with respect to idx and cell, with
        cell changing fastest.
    """
    if mesh.structured:
       num_rays = kwargs['num_rays'] if 'num_rays' in kwargs else 10
       grid = kwargs['grid'] if 'grid' in kwargs else False
       results = ray_discretize(mesh, num_rays, grid)
    else:
       if kwargs:
           raise ValueError("No valid key word arguments for unstructed mesh.")
       cells = cells_at_ve_centers(mesh)
       # Use str for python2/3 compatibility
       results = np.zeros(len(mesh), dtype=[(str('idx'), np.int64),
                                            (str('cell'), np.int64),
                                            (str('vol_frac'), np.float64),
                                            (str('rel_error'), np.float64)])
       for i, cell in enumerate(cells):
           results[i] = (i, cells[i], 1.0, 1.0)

    return results

def cells_at_ve_centers(mesh):
    """cells_at_ve_centers(mesh)
    This function reads in any PyNE Mesh object and finds the geometry cell
    at the point in the center of each mesh volume element. A DAGMC geometry
    must be loaded prior to using this function.

    Parameters
    ----------
    mesh : PyNE Mesh
        Any Mesh that is superimposed over the geometry.

    Returns
    -------
    cells : list
        The cell numbers of the geometry cells that occupy the center of the
        mesh volume element, in the order of the mesh idx.
    """
    cells = []
    for i, mat, ve in mesh:
        center = mesh.ve_center(ve)
        cell = find_volume(center)
        cells.append(cell)

    return cells

def ray_discretize(mesh, num_rays=10, grid=False):
    """ray_discretize(mesh, num_rays=10, grid=False)
    This function discretizes a geometry (by geometry cell) onto a
    superimposed, structured, axis-aligned mesh using the method described in
    [1]. Ray tracing is used to sample track lengths in geometry cells in mesh
    volume elements, and volume fractions are determined statiscally. Rays are
    fired down entire mesh rows, in three directions: x, y, and z. Rays starting
    points can be chosen randomly, or on a uniform grid. Note that a DAGMC
    geometry must already be loaded into memory.

    [1] Moule, D. and Wilson, P., Mesh Generation methods for Deterministic
    Radiation Transport Codes, Transacztions of the American Nuclear Society,
    104, 407--408, (2009).

    Parameters
    ----------
    mesh : PyNE Mesh
        A Cartesian, structured, axis-aligned Mesh that superimposed the
        geometry.
    num_rays : int, optional, default = 10
        The number of rays to fire in each mesh row for each direction.
    grid : boolean, optional, default = False
        If false, rays starting points are chosen randomly (on the boundary)
        for each mesh row. If true, a linearly spaced grid of starting points is
        used, with dimension sqrt(num_rays) x sqrt(num_rays). In this case,
        "num_rays" must be a perfect square.

    Returns
    -------
    results : structured array
        Stores in a one dimensional array, each entry containing the following
        fields:
        :idx: int
            The mesh volume element index.
        :cell: int
            The geometry cell number.
        :vol_frac: float
            The volume fraction of the cell withing the mesh volume element.
        :rel_error: float
            The relative error associated with the volume fraction.
        This array is returned in sorted order with respect to idx and cell, with
        cell changing fastest.
    """
    mesh._structured_check()
    # add the str here to prevent the 'xyz' be transfered to ascii
    divs = [mesh.structured_get_divisions(x) for x in str('xyz')]
    num_ves = (len(divs[0])-1)*(len(divs[1])-1)*(len(divs[2])-1)
    #  Stores a running tally of sums of x and sums of x^2 for each ve
    mesh_sums = [{} for x in range(0, num_ves)]
    #  The length of the output array (will vary because different ve contain
    #  different numbers of geometry cells.
    len_count = 0

    #  Direction indicies: x = 0, y = 1, z = 2
    dis = [0, 1, 2]
    #  Iterate over all directions indicies
    for di in dis:
        #  For each direction, the remaining two directions define the sampling
        #  surface. These two directions are the values in s_dis (surface
        #  direction indices)
        s_dis = [0, 1, 2]
        s_dis.remove(di)
        #  Iterate through all all the sampling planes perpendicular to di,
        #  creating a _MeshRow in each, and subsequently evaluating that row.
        for a in range(0, len(divs[s_dis[0]]) - 1):
            for b in range(0, len(divs[s_dis[1]]) - 1):
                mesh_row = _MeshRow()
                mesh_row.di = di
                mesh_row.divs = divs[di]
                mesh_row.num_rays = num_rays
                mesh_row.s_dis_0 = s_dis[0]
                mesh_row.s_min_0 = divs[s_dis[0]][a]
                mesh_row.s_max_0 = divs[s_dis[0]][a + 1]
                mesh_row.s_dis_1 = s_dis[1]
                mesh_row.s_min_1 = divs[s_dis[1]][b]
                mesh_row.s_max_1 = divs[s_dis[1]][b + 1]

                #  Create a lines of starting points to fire rays for this
                #  particular mesh row.
                if not grid:
                    mesh_row._rand_start()
                else:
                    mesh_row._grid_start()

                #  Create a list of mesh idx corresponding to this mesh row.
                if di == 0:
                    ves = mesh.structured_iterate_hex('x', y=a, z=b)
                elif di == 1:
                    ves = mesh.structured_iterate_hex('y', x=a, z=b)
                elif di == 2:
                    ves = mesh.structured_iterate_hex('z', x=a, y=b)

                idx_tag = mesh.idx
                idx = []
                for ve in ves:
                    idx.append(idx_tag[ve][0])

                #  Fire rays.
                row_sums = mesh_row._evaluate_row()

                #  Add row results to the full mesh sum matrix.
                for j, ve_sums in enumerate(row_sums):
                   for cell in ve_sums.keys():
                       if ve_sums[cell][0] < VOL_FRAC_TOLERANCE:
                           continue
                       if cell not in mesh_sums[idx[j]].keys():
                           mesh_sums[idx[j]][cell] = [0, 0]
                           len_count += 1

                       mesh_sums[idx[j]][cell][0] += ve_sums[cell][0]
                       mesh_sums[idx[j]][cell][1] += ve_sums[cell][1]


    #  Create structured array
    total_rays = num_rays*3 # three directions
    # Use str for python2/3 compatibility
    results = np.zeros(len_count, dtype=[(str('idx'), np.int64),
                                         (str('cell'), np.int64),
                                         (str('vol_frac'), np.float64),
                                         (str('rel_error'), np.float64)])

    row_count = 0
    total_rays = num_rays*3
    for i, ve_sums in enumerate(mesh_sums):
       for vol in ve_sums.keys():
           vol_frac = ve_sums[vol][0]/total_rays
           rel_error = np.sqrt((ve_sums[vol][1])/(ve_sums[vol][0])**2
                                - 1.0/total_rays)
           results[row_count] = (i, vol, vol_frac, rel_error)
           row_count += 1

    results.sort()

    return results

class _MeshRow():
    """A class to store data and fire rays down a single mesh row.

    Attributes
    ----------
    di : int
        The direction index of the current firing direction.
    divs : list
        The mesh boundaries perpendicular to the firing direction.
    start_points : list
        The xyz points describing where rays should start from.
    num_rays : int
        The number of rays to fire down the mesh row.
    s_dis_0 : int
        The first of two directions that define the surface for which rays
        are fired.
    s_min_0 : float
        The location of the plane the bounds the firing surface from the left
        in direction s_dis_0.
    s_max_0 : float
        The location of the plane the bounds the firing surface from the right
        in direction s_dis_0.
    s_dis_1 : int
        The second of two directions that define the surface for which rays
        are fired.
    s_min_1 : float
        The location of the plane the bounds the firing surface from the left
        in direction s_dis_1.
    s_max_1 : float
        The location of the plane the bounds the firing surface from the right
        in direction s_dis_1.
       """
    def __init__(self):
        pass

    def _rand_start(self):
        """Private function for randomly generating ray starting points to
        populate self.starting_points
        """
        self.start_points = []
        ray_count = 0
        while ray_count < self.num_rays:
            start_point = [0]*3
            start_point[self.di] = self.divs[0]
            start_point[self.s_dis_0] = np.random.uniform(self.s_min_0,
                                                          self.s_max_0)
            start_point[self.s_dis_1] = np.random.uniform(self.s_min_1,
                                                          self.s_max_1)
            self.start_points.append(start_point)
            ray_count += 1

    def _grid_start(self):
        """Private function for generating a uniform grid of ray starting
        points to populate self.starting_points.
        """
        #  Test to see if num_rays is a perfect square.
        if int(np.sqrt(self.num_rays))**2 != self.num_rays:
            raise ValueError("For rays fired in a grid, "
                             "num_rays must be a perfect square.")
        else:
           square_dim = int(np.sqrt(self.num_rays))

        step_1 = (self.s_max_0-self.s_min_0)/(float(square_dim) + 1)
        step_2 = (self.s_max_1 - self.s_min_1)/(float(square_dim) + 1)
        range_1 = np.linspace(self.s_min_0 + step_1, self.s_max_0, square_dim,
                              endpoint=False)
        range_2 = np.linspace(self.s_min_1 + step_2, self.s_max_1, square_dim,
                              endpoint=False)

        self.start_points = []
        for point_1 in range_1:
            for point_2 in range_2:
                start_point = [0]*3
                start_point[self.di] = self.divs[0]
                start_point[self.s_dis_0] = point_1
                start_point[self.s_dis_1]= point_2
                self.start_points.append(start_point)

    def _evaluate_row(self):
        """Private function that fires rays down a single mesh row and returns the
        results.

        Returns
        -------
        row_sums : list
            A list with one entry per mesh volume element in the mesh row. Each
            entry is a dictionary that maps geometry cell number to a list
            containing two statistical results. The first result is the sum of
            all the samples and the second result is the sum of all the squares
            of all of the samples.
        """

        # Number of volume elements in this mesh row.
        num_ve = len(self.divs) - 1
        direction = [0, 0, 0]
        direction[self.di] = 1
        row_sums = [{} for x in range(0, len(self.divs) - 1)]
        width = [self.divs[x] - self.divs[x - 1] for x in range(1, len(self.divs))]

        #  Find the first volume the first point is located in.
        vol = find_volume(self.start_points[0], direction)
        #  Fire ray for each starting point.
        for point in self.start_points:
            #  Check to see if the staring point is in the same volume as the
            #  last staring point to avoid expensive find_volume calls.
            if not point_in_volume(vol, point, direction):
                vol = find_volume(point, direction)

            mesh_dist = width[0]
            ve_count = 0
            complete = False
            #  Track a single ray down the mesh row and tally accordingly.
            for next_vol, distance, _ in ray_iterator(vol, point, direction):
                if complete:
                    break

                #  Volume extends past mesh boundary
                while distance >= mesh_dist:
                    #  Check to see if current volume has already by tallied
                    if vol not in row_sums[ve_count].keys():
                        row_sums[ve_count][vol] = [0, 0]

                    sample = mesh_dist/width[ve_count]
                    row_sums[ve_count][vol][0] += sample
                    row_sums[ve_count][vol][1] += sample**2
                    distance -= mesh_dist

                    #  If not on the last volume element, continue into the next
                    #  volume element.
                    if ve_count == num_ve - 1:
                        complete = True
                        break
                    else:
                        ve_count += 1
                        mesh_dist = width[ve_count]

                #  Volume does not extend past mesh volume.
                if distance < mesh_dist and distance > VOL_FRAC_TOLERANCE \
                    and not complete:
                    #  Check to see if current volume has already by tallied
                    if vol not in row_sums[ve_count].keys():
                        row_sums[ve_count][vol] = [0, 0]

                    sample = distance/width[ve_count]
                    row_sums[ve_count][vol][0] += sample
                    row_sums[ve_count][vol][1] += sample**2
                    mesh_dist -= distance

                vol = next_vol

        return row_sums
