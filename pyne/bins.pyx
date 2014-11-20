"""Tools to generate and handle various binning structures."""
from warnings import warn

cimport cython

cimport numpy as np
import numpy as np
from numpy import logspace

from pyne.utils import QAWarning

warn(__name__ + " is not yet QA compliant.", QAWarning)

def ninespace(start, stop, num=50, endpoint=True):
    """Splits the range into one-minus-log-uniform bins defined by num points.
    In the vernacular, the space is 'split in the nines'.  Note that this assumes
    base = 10.0.

    Parameters
    ----------
    start : number
        The starting value of the sequence.
    stop : number
        The final value of the sequence, unless endpoint
        is False.  In that case, num + 1 values are spaced over the
        interval in nines-space, of which all but the last (a sequence of
        length num) are returned.
    num : integer, optional
        Number of samples to generate. See endpoint.
    endpoint : boolean, optional
        If true, stop is the last sample. Otherwise, it is not included.

    Returns
    -------
    samples : ndarray
        num samples, equally spaced in the nines.

    Examples
    --------
    >>> ninespace(0.9, 0.9999, 4)
        array([0.9, 0.99, 0.999, 0.9999])

    """
    log_start = np.log10(1.0 - start)
    log_stop  = np.log10(1.0 - stop)
    samples = 1.0 - logspace(log_start, log_stop, num, endpoint)
    return samples


def stair_step(x, y):
    """Makes a 1d data set of boundaries (x) and cell-centered values (y) 
    into stair step arrays of the same length.  The returned arrays are 
    suitable for graphing.  This is especially useful in energy vs. spectrum
    data where there are G+1 boundaries and G data points.

    Parameters
    ----------
    x : sequence
        Data of length G+1.
    y : sequence
        Data of length G.

    Returns
    -------
    xss : ndarray 
        Stair-step version of x data, length 2G.
    yss : ndarray 
        Stair-step version of y data, length 2G.

    Examples
    --------
    >>> x = [0.1, 1.0, 10.0, 100.0]
    >>> y = [2.0, 3.0, 4.0]
    >>> bins.stair_step(x, y)
    (array([   0.1,    1. ,    1. ,   10. ,   10. ,  100. ]),
    array([ 2.,  2.,  3.,  3.,  4.,  4.]))

    """
    # Grab the number of points, G
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    G = len(y)
    assert G + 1 == len(x)

    # Initilaize data
    xss = np.empty(2*G, dtype=x.dtype)
    yss = np.empty(2*G, dtype=y.dtype)

    # Set values
    xss[:-1:2] = x[:-1]
    xss[1::2] = x[1:]
 
    yss[::2] = y
    yss[1::2] = y

    return xss, yss

@cython.boundscheck(False)
def pointwise_linear_collapse(np.ndarray[np.float64_t, ndim=1] x_g, 
                              np.ndarray[np.float64_t, ndim=1] x, 
                              np.ndarray[np.float64_t, ndim=1] y):
    """Collapses pointwise data to G groups based on a linear interpolation
    between the points. This is useful for collapsing cross section data.

    Parameters
    ----------
    x_g : array-like
        Group boundaries, length G+1 for G groups, must be monotonic in the 
        same direction as x.
    x : array-like
        Pointwise abscissa to be collapsed, must be monotonic in the same direction
        as x_g and have the same length as y.
    y : array-like
        Pointwise data to be interpolated, must have the same length as x.

    Returns
    -------
    y_g : np.ndarray
        The group collapsed data, length G. 
    """
    cdef int G = x_g.shape[0] - 1
    cdef int N = x.shape[0]
    cdef int g0, g1  # current group index
    cdef int n0, n1  # current point index
    cdef double val, ylower, yupper
    cdef np.ndarray[np.float64_t, ndim=1] y_g = np.empty(G, dtype='float64')
    reversed = False
    if x_g[0] > x_g[-1]:
        # monotonically decreaing, make increasing so logic is simpler
        x_g = x_g[::-1]
        x = x[::-1]
        y = y[::-1]
        reversed = True
    n0 = 0
    n1 = 1
    for g0 in range(G):
        g1 = g0 + 1
        val = 0.0
        while x[n1] <= x_g[g1] and n1 < N:
            if x_g[g0] <= x[n0]:
                # interpolation fully in group, can take midpoint
                val += 0.5 * (y[n1] + y[n0]) * (x[n1] - x[n0])
            else: 
                # lower bound intersection
                ylower = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g0] - x[n0]) + y[n0]
                val += 0.5 * (y[n1] + ylower) * (x[n1] - x_g[g0])
            n0 += 1
            n1 += 1
        # upper bound intersection
        if x_g[g1] < x[n1]:
            if x_g[g0] <= x[n0]:
                yupper = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g1] - x[n0]) + y[n0]
                val += 0.5 * (yupper + y[n0]) * (x_g[g1] - x[n0])
            else: 
                yupper = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g1] - x[n0]) + y[n0]
                ylower = ((y[n1] - y[n0])/(x[n1] - x[n0]))*(x_g[g0] - x[n0]) + y[n0]
                val += 0.5 * (yupper + ylower) * (x_g[g1] - x_g[g0])
        y_g[g0] = val / (x_g[g1] - x_g[g0])
    if reversed:
        y_g = y_g[::-1]
    return y_g
    
