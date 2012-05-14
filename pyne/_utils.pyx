"""Cython-based utils to be imported into utils."""
from libc.stdlib cimport malloc, free

cimport numpy as np
import numpy as np

from std cimport atof, strtok, strcpy


def fastfromstring(char * s, char * sep=" ", bint inplace=False):
    """A quick replacement for numpy.fromstring() using the C standard
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

    Notes
    -----
    Will always return a 1d float64 array.  You must cast and reshape 
    to the appropriate type and shape.
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

    I = (I / 2) + 1
    data = np.empty(I, dtype=np.float64)
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


