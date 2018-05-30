from __future__ import division
import os
import math

from distutils.dir_util import remove_tree

from pyne._utils import fromstring_split, fromstring_token, endftod,\
                        use_fast_endftod, fromendf_tok, toggle_warnings,\
                        use_warnings, fromendl_tok


class QAWarning(UserWarning):
    pass

time_conv_dict = {'as': 1e-18,
                  'attosec': 1e-18,
                  'attosecond': 1e-18,
                  'attoseconds': 1e-18,
                  'fs': 1e-15,
                  'femtosec': 1e-15,
                  'femtosecond': 1e-15,
                  'femtoseconds': 1e-15,
                  'ps': 1e-12,
                  'picosec': 1e-12,
                  'picosecond': 1e-12,
                  'picoseconds': 1e-12,
                  'ns': 1e-9,
                  'nanosec': 1e-9,
                  'nanosecond': 1e-9,
                  'nanoseconds': 1e-9,
                  'us': 1e-6,
                  'microsec': 1e-6,
                  'microsecond': 1e-6,
                  'microseconds': 1e-6,
                  'ms': 1e-3,
                  'millisec': 1e-3,
                  'millisecond': 1e-3,
                  'milliseconds': 1e-3,
                  's': 1.0,
                  'sec': 1.0,
                  'second': 1.0,
                  'seconds': 1.0,
                  'm': 60.0,
                  'min': 60.0,
                  'minute': 60.0,
                  'minutes': 60.0,
                  'h': 3600.0,
                  'hour': 3600.0,
                  'hours': 3600.0,
                  'd': 86400.0,
                  'day': 86400.0,
                  'days': 86400.0,
                  'w': 86400.0*7.0,
                  'week': 86400.0*7.0,
                  'weeks': 86400.0*7.0,
                  'y': 86400.0*365.25,
                  'year': 86400.0*365.25,
                  'years': 86400.0*365.25,
                  'c': 86400.0*365.25*100,
                  'century': 86400.0*365.25*100,
                  'centuries': 86400.0*365.25*100,
                  }


def to_sec(input_time, units):
    """Converts a time with units to seconds.

    Parameters
    ----------
    input_time : number
        Time value in [units].
    units : str
        Units flag, eg 'min', 'ms', 'days'

    Returns
    -------
    sec_time : float
        Time value in [sec].

    """
    conv = time_conv_dict.get(units.lower(), None)
    if conv:
        sec_time = input_time * conv
        return sec_time
    else:
        raise ValueError('Invalid units: {0}'.format(units))

barn_conv_dict = {
    'mb': 1E-3,
    'ub': 1E-6,
    'microbarn': 1E-6,
    'b': 1.0,
    'barn': 1.0,
    'barns': 1.0,
    'kb': 1E+3,
    'kilobarn': 1E+3,
    'cm2': 1E+24,
    'cm^2': 1E+24,
    }


def to_barns(xs, units):
    """Converts a cross section with units to barns.

    Parameters
    ----------
    xs :
        Cross section value in [units].
    units : str
        Units flag, eg 'b', 'microbarn'.

    Returns
    -------
    barn_xs :
        Cross section value in [barns].

    """
    return xs * barn_conv_dict[units.lower()]


def from_barns(xs, units):
    """Converts a cross section from barns to units.

    Parameters
    ----------
    xs :
        Cross section value in [barns].
    units : str
        Units flag, eg 'b', 'microbarn'.

    Returns
    -------
    unit_xs :
        Cross section value in [units].

    """
    return xs / barn_conv_dict[units.lower()]


#########################
### message functions ###
#########################

USE_COLOR = (os.name is 'posix')


def message(s):
    """Formats a message for printing.  If on a posix system the message will
    be in color.

    """
    head = "\033[1;32m" if USE_COLOR else "*** MESSAGE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def failure(s):
    """Formats a fail message for printing.  If on a posix system the message
    will be in color.

    """
    head = "\033[1;31m" if USE_COLOR else "*** FAILURE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def warning(s):
    """Formats a warning message for printing. If on a posix system the message
    will be in color.
    """
    head = "\033[1;33m" if USE_COLOR else "*** WARNING ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg

##################################
### Path manipulation routines ###
##################################


def remove(path):
    """Removes a path, or recursively a directory, or does nothing
    if path is neither a file nor a directory.

    """
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        remove_tree(path, verbose=False)
    else:
        pass

##################################
### Facet voxel traverse tools ###
##################################


class Point(object):
    """
    Point class represents x, y and z coordinates.
    Coordinate of a dimension could be any float, including 'inf' and '-inf'.
    """
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """
        Return a string with readable information of the coordinate.
        """
        return ''.join(['(',str(self.x), ', ',
                        str(self.y), ', ',
                        str(self.z), ')'])

def _distance(start, end):
    """
    Calculate the distance between two points.

    Parameters:
    -----------
    start : Point
        The start point of the line segment.
    end : Point
        The end point of the line segment.

    Return:
    -------
    dist : float
        The distance between the start point and the end point.
    """
    dist = math.sqrt(math.pow(end.x - start.x, 2) +
                      math.pow(end.y - start.y, 2) +
                      math.pow(end.z - start.z, 2))
    return dist

def _is_point_in_mesh(bounds, point):
    """
    Check whether a point is in the mesh.

    Parameters:
    -----------
    bounds : list
        Boundaries. List of the boundaries along x, y, and z directions.
    point : Point
        The point to be checked.

    Return:
    -------
    True : when point in the mesh.
    False : when point outside the mesh.
    """
    if bounds[0][0] <= point.x <= bounds[0][-1] and \
       bounds[1][0] <= point.y <= bounds[1][-1] and \
       bounds[2][0] <= point.z <= bounds[2][-1]:
        return True
    else:
        return False

def _find_voxel_idx_1d(bound, cor, v_1d):
    """
    Find the voxel id (in a given dimension).

    Parameters:
    -----------
    bound : list
        Boundary of the mesh in a given direction.
    cor : float
        Coordinate of the point in the given direction.
    v_1d : float
        Vector attribute in the given direction.

    Return:
    -------
    idx : int
        Mesh index of the given direction. This value is -1 when the point not
        in the mesh.
    ----------
    """
    idx = -1
    if v_1d > 0:
        for i in range(len(bound)-1):
            if bound[i] <= cor < bound[i+1]:
                idx = i
                break
    else:
        for i in range(len(bound)-1, -1, -1):
            if bound[i] < cor <= bound[i+1]:
                idx = i
                break
    # idx could be -1, that means point outside the mesh
    return idx

def _cal_dir_v(start, end):
    """
    Calculate the direction from start to end (start -> end).
    
    Parameters:
    -----------
    start : Point
        Start point.
    end :
        End point.

    Return:
    -------
    v : list
        direction from start to end, an unit vector.
    """
    v_x, v_y, v_z = 0.0, 0.0, 0.0
    d_x = end.x - start.x
    d_y = end.y - start.y
    d_z = end.z - start.z
    dist = math.sqrt(d_x * d_x + d_y * d_y + d_z * d_z)
    v_x = d_x / dist
    v_y = d_y / dist
    v_z = d_z / dist
    return [v_x, v_y, v_z]

def _find_next_grid_1d(cor, v_1d, bound):
    """
    Find the next grid coordinate of a direction.

    Parameters:
    -----------
    cor : float
        Coordinate of current point in the direction.
    v_1d : float
        Direction value.
    bound : list of float
        Boundary values of current direction.

    Return:
    -------
    n_grid : float
        The coordinate value of next mesh grid.
    """
    n_grid = 0.0
    flag = False
    if v_1d > 0:
        for i in range(len(bound)):
            if cor < bound[i]:
                n_grid = bound[i]
                flag = True
                break
    else:
        for i in range(len(bound)-1, -1, -1):
            if cor > bound[i]:
                n_grid = bound[i]
                flag = True
                break
    if flag:
        return n_grid
    else:
        # point goes out of the mesh
        return v_1d * float('inf')


def _cal_max_travel_length_in_current_voxel(tp, end, v, bounds):
    """
    Calculate the maximum travel length in current voxel before reach next
    boundary.
    
    Parameters:
    -----------
    tp : Point
        Temporary point position.
    end : Point
        End point
    v : list
        Direction vector from start to end: [v_x, v_y, v_z].
    bounds : list
        Boundaries of the mesh grid.

    Return:
    -------
    t_maxs : list
        The maximum travel lenght in current voxel for 3 directions. Format: 
        [t_max_x, t_max_y, t_max_z]. These value could be both positive and
        negtive. And could be 'inf' or '-inf'.
    """
    n_x_grid = _find_next_grid_1d(tp.x, v[0], bounds[0])
    n_y_grid = _find_next_grid_1d(tp.y, v[1], bounds[1])
    if v[0] == 0.0:
        t_max_x = float('inf')
    else:
        t_max_x = (n_x_grid - tp.x) / v[0]
    if v[1] == 0.0:
        t_max_y = float('inf')
    else:
        t_max_y = (n_y_grid - tp.y) / v[1]
    return [t_max_x, t_max_y]

def _move_to_next_voxel(x_idx, y_idx, z_idx, tp, t_temp, t_maxs, v):
    """
    Move the point to next voxel.

    Parameters:
    -----------
    x_idx : int
        Mesh idx in x direction.
    y_idx : int
        Mesh idx in y direction.
    z_idx : int
        Mesh idx in z direction.
    tp : Point
        Temporary point position.
    t_temp: float
        Temporary travel lenght from the start.
    t_maxs: list
        The maximum travel length in current voxel for 3 directions.
    v : list
        Direction vector from start to end.

    Returns:
    --------
    x_idx : int
        Modified x_idx.
    y_idx : int
        Modified y_idx.
    z_idx : int
        Modified z_idx.
    tp : Point
        Point after moving along the direction of v.
    t_temp : float
        Modifed travel length from the start.
    """
    t_min = min(t_maxs)
    if t_maxs[0] == t_min:
        x_idx += 1
    elif t_maxs[1] == t_min:
        y_idx += 1
    elif t_maxs[2] == t_min:
        z_idx += 1
    tp.x += t_min * v[0]
    tp.y += t_min * v[1]
    tp.z += t_min * v[2]
    t_temp += t_min
    return x_idx, y_idx, z_idx, tp, t_temp

def _move_to_boundary(tp, t_min, v):
    """
    Move the point to the nearest boundary when Initialize the problem.

    Parameters:
    -----------
    tp : Point
        Temporary point. Actually, it is also the start point.
    t_min : float
        Distance to the nearest (positive) boundary. For a 3D mesh, there are 6
        outside boundaries. Therefore there are 6 value of the distance from the
        start to the boundaries. This t_min is the smallest one among the
        positive distance. (The distance is negtive when the direction is
        negtive.)
    v : list
        Direction vector.

    Return:
    -------
    tp : Point
        Modified point. This point lies on the nearest boundary. Therefore, we
        moved point from outside the mesh to the mesh boundary (or the extension
        of the boundary, which still outside the mesh).
    """

    # calculate the new coordinate
    tp.x += t_min * v[0]
    tp.y += t_min * v[1]
    return tp

def _ray_voxel_traverse(bounds, start, end):
    """
    Voxel traversal algorithm.
    Return the voxels that intersect with the line segment.

    Parameters:
    -----------
    bounds : list
        Boundaries of the mesh grid.
    start : Point
        Start point.
    end : Point
        End point.

    Return:
    -------
    voxles : list
        List of voxel indices. These voxels intersect with the line segment.
    """
    # Initialization
    # Find the voxel that start point in
    tp = Point(start.x, start.y, start.z)
    x_idx, y_idx, z_idx = -1, -1, -1
    voxel_list = []
    v = _cal_dir_v(start, end)
    t_end = _distance(start, end)
    t_temp = 0.0
    while t_temp < t_end:
        if _is_point_in_mesh(bounds, tp):
            # the point is in the mesh
            # calculate the voxel id
            x_idx = _find_voxel_idx_1d(bounds[0], tp.x, v[0])
            y_idx = _find_voxel_idx_1d(bounds[1], tp.y, v[1])
            z_idx = _find_voxel_idx_1d(bounds[2], tp.z, v[2])
            if -1 in (x_idx, y_idx, z_idx):
                # point outside the mesh
                return voxel_list
            # add current voxel
            voxel_list.append([x_idx, y_idx, z_idx])
            t_maxs = _cal_max_travel_length_in_current_voxel(tp, end, v, bounds)
            x_idx, y_idx, z_idx, tp, t_temp = _move_to_next_voxel(
                x_idx, y_idx, z_idx, tp, t_temp, t_maxs, v)
        else:
            # Check whether the ray intersecs the mesh. Calculate the maxinum travel
            # length of the ray to the outside boundaries of the mesh. If all the
            # length are negtive, then there is no intersection, return. Ohterwise,
            # choose the minimum length and move the point to that position.
            # left, right, back, front, down, up
            t_maxs = [0.0] * 4
            if v[0] == 0.0:
                t_maxs[0] = float('inf')
                t_maxs[1] = float('inf')
            else:
                t_maxs[0] = (bounds[0][0] - tp.x) / v[0]
                t_maxs[1] = (bounds[0][-1] - tp.x) / v[0]

            if v[1] == 0.0:
                t_maxs[2] = float('inf')
                t_maxs[3] = float('inf')
            else:
                t_maxs[2] = (bounds[1][0] - tp.y) / v[1]
                t_maxs[3] = (bounds[1][-1] - tp.y) / v[1]

            if v[2] == 0.0:
                t_maxs[4] = float('inf')
                t_maxs[5] = float('inf')
            else:
                t_maxs[4] = (bounds[2][0] - tp.z) / v[2]
                t_maxs[5] = (bounds[2][-1] - tp.z) / v[2]

            if all(t_maxs) < 0:
                # current point outside the mesh, not any intersection
                return voxel_list
            else:
                t_min = min(i for i in t_maxs if i > 0)
                t_temp += t_min
                # Move the point to the boundary
                # calculate the new coordinate
                tp.x += t_min * v[0]
                tp.y += t_min * v[1]
                tp.z += t_min * v[2]
    return voxel_list
 
def _cal_min_grid_step(bounds):
    """
    Calculate the minimum grid step. Used for scan.

    Parameters:
    -----------
    bounds : list
        Boundaries of mesh grid.

    Return:
    -------
    min_step: float
        The minimum grid step.
    """
    min_step = float('inf')
    for i in range(len(bounds)):
        for j in range(len(bounds[i])-1):
            step = bounds[i][j+1] - bounds[i][j]
            if min_step > step:
                min_step = step
    return min_step

def _divide_edge(start, end, step):
    """
    Divide the edge into several parts. Return the divide points.

    Parameters:
    -----------
    start : Point
        Start point.
    end : Point
        End point.
    step : float
        The minimum grid step. Divider used to divide the line segment.

    Return:
    -------
    intert_points : list of points
        List of points. These points are located on the line segment from start
        to end.
    """
    dist = _distance(start, end)
    v = _cal_dir_v(start, end)
    t_temp = 0.0
    tp = Point(start.x, start.y, start.z)
    insert_points = [start]
    while t_temp < dist:
        tp.x += step * v[0]
        tp.y += step * v[1]
        tp.z += step * v[2]
        insert_points.append(Point(tp.x, tp.y, tp.z))
        t_temp += step
    if end not in insert_points:
        insert_points.append(end)
    return insert_points

def _create_rays_from_points(start, insert_points):
    """
    Create rays from points. In which, 'start' is the start point of all the
    rays, and the point in the `insert_points' are the end points.

    Parameters:
    -----------
    start : Point
        Start point for all the rays (line segments).
    insert_points : list of points
        List of points.

    Return:
    -------
    rays: list
        List of rays (line segments).
    """
    rays = []
    for p in insert_points:
        rays.append([start, p])
    return rays

def _create_rays_from_triangle_facet(A, B, C, min_step):
    """
    Create rays form triangle for scanning.
    Divide the shortest edge of the triangle with the min_step.

    Parameters:
    -----------
    A : Point
        Vertex of the triangle
    B : Point
        Vertex of the triangle
    C : Point
        Vertex of the triangle
    min_step: float
        The minimum step of the mesh grid. Used to divide the shortest edge
        of the facet triangle. Use this value could make sure that every voxel
        that intersect with the facet has intersection with some of the rays.

    Return:
    -------
    rays : list of rays, [start, end] pairs
        List of the created rays. Used to represent the facet.
    """
    a = _distance(B, C)
    b = _distance(A, C)
    c = _distance(A, B)
    min_edge = min([a, b, c])
    if a == min_edge:
        # divide BC
        insert_points = _divide_edge(B, C, min_step)
        rays = _create_rays_from_points(A, insert_points)
    elif b == min_edge:
        # divide AC
        insert_points = _divide_edge(A, C, min_step)
        rays = _create_rays_from_points(B, insert_points)
    else:
        # divide AB
        insert_points = _divide_edge(A, B, min_step)
        rays = _create_rays_from_points(C, insert_points)
    return rays

def _merge_voxel_list(old_list, new_list):
    """
    Merge two voxel lists.

    Parameters:
    -----------
    old_list : list
        Old list of voxels.
    new_list : list
        New list of the voxels.

    Return:
    -------
    old_list : list
        Merged list.
    """
    for v in new_list:
        if v not in old_list:
            old_list.append(v)
    return old_list

def _facet_voxel_traverse(A, B, C, bounds):
    """
    Facet voxel traverse. Calculate all the voxles that intersect with the
    facet. Facet is represented with 3 triangle verteice.

    Parameters:
    -----------
    A : Point
        Vertex of the triangle facet
    B : Point
        Vertex of the triangle facet
    C : Point
        Vertex of the triangle facet
    bounds : list
        Boundaries of the mesh grid.

    Return:
    -------
    voxel_list: list
        List of voxels that intersect with the facet.
    """
    min_step = _cal_min_grid_step(bounds)
    rays = _create_rays_from_triangle_facet(A, B, C, min_step)
    voxel_list = []
    for ray in rays:
        temp_voxel_list = _ray_voxel_traverse(bounds, ray[0], ray[1])
        # merge two lists
        voxel_list = _merge_voxel_list(voxel_list, temp_voxel_list)
    return voxel_list

