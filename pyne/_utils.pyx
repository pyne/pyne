"""Cython-based utils to be imported into utils."""
from libc.stdlib cimport malloc, free

cimport numpy as np
import numpy as np

from std cimport atof, strtok, strcpy


def fromstring_split(char * s, sep=None, dtype=float):
    """A replacement for numpy.fromstring() using the Python str.split() 
    and np.array().

    Parameters
    ----------
    s : str
        String of data.
    sep : str or None
        String of separator characters, has the same meaning as in 
        str.split().
    dtype : np.dtype
        Numpy dtype to cast elements enough.

    Returns
    -------
    data : ndarry, 1d
        Will always return a 1d array of dtype.  You must reshape to the 
        appropriate shape.

    See Also
    --------
    fromstring_token : May faster depending on the data. 

    """
    cdef list rawdata
    rawdata = s.split(sep)
    return np.array(rawdata, dtype=dtype)


def fromstring_token(char * s, char * sep=" ", bint inplace=False, int maxsize=-1):
    """A replacement for numpy.fromstring() using the C standard
    library atof() and strtok() functions.

    Parameters
    ----------
    s : str
        String of data.
    sep : str
        String of separator characters.  Unlike numpy.fromstring(), 
        all characters are separated on independently.
    inplace : bool
        Whether s should tokenized in-place or whether a copy should 
        be made.  If done in-place, the first instance of sep between 
        any tokens will replaced with the NULL character.
    maxsize : int
        Specifies the size of the array to pre-allocate.  If negative,
        this will be set to the maximum possible number of elements, 
        ie len(s)/2 + 1.

    Returns
    -------
    data : ndarry, 1d, float64
        Will always return a 1d float64 array.  You must cast and reshape 
        to the appropriate type and shape.

    See Also
    --------
    fromstring_split : May faster depending on the data.

    """
    cdef char* cstring
    cdef char* cs
    cdef char* csep
    cdef int i, I
    cdef np.ndarray[np.float64_t, ndim=1] cdata

    I = len(s)
    csep = sep

    if inplace:
        cs = s
    else:
        cs = <char *> malloc(I * sizeof(char))
        strcpy(cs, s)

    if maxsize < 0:
        maxsize = (I / 2) + 1

    data = np.empty(maxsize, dtype=np.float64)
    cdata = data

    i = 0
    cstring = strtok(cs, csep)
    while cstring != NULL:
        cdata[i] = atof(cstring)
        cstring = strtok(NULL, csep)
        i += 1

    if not inplace:
        free(cs)

    data = data[:i].copy()
    return data
