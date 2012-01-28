"""Module for tools to generate and handle various binning structures."""

import numpy as np
from numpy import linspace, logspace


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
