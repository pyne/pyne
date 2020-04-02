from pyne cimport stlcontainers
from libc.stdlib cimport free
from libcpp.string cimport string as std_string
import numpy as np
cimport numpy as np

import collections
from pyne import stlcontainers

from warnings import warn
from pyne.utils import QAWarning

warn(__name__ + " is not yet QA compliant.", QAWarning)

cdef class PointSource:
    """

    Attributes
    ----------


    Methods
    -------

    Notes
    -----
    This class was defined in source.h

    The class is found in the "pyne" namespace"""

    # constuctors
    def __cinit__(self, x=0.0, y=0.0, z=0.0, i=0.0, j=0.0, k=0.0, E=14., p="n", weight=1):
        p_bytes = p.encode()
        self._inst = new cpp_source.PointSource(x, y, z, i, j, k, E, std_string(<char *> p_bytes), weight)
        self._free_inst = True

    def mcnp(self, version=5):
        cdef std_string card
        card = self._inst.mcnp(version)
        return card.decode()
    
    
    def __dealloc__(self):
        """PointSource C++ destructor."""
        if self._free_inst:
            del self._inst
    
    # attributes
    property x:
        """no docstring for x, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).x

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).x = <double> value
    
    property y:
        """no docstring for y, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).y

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).y = <double> value
    
    property z:
        """no docstring for z, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).z

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).z = <double> value
    
    property i:
        """no docstring for i, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).i

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).i = <double> value
    
    property j:
        """no docstring for j, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).j

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).j = <double> value
    
    property k:
        """no docstring for k, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).k

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).k = <double> value
    
    property E:
        """no docstring for E, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).E

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).E = <double> value
    
    property weight:
        """no docstring for weight, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).weight

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).weight = <double> value
    
    property particle:
        """no docstring for particle, please file a bug report!"""
        def __get__(self):
            return bytes(<char *> (<cpp_source.PointSource *> self._inst).particle.c_str()).decode()

        def __set__(self, value):
            cdef char * value_proxy
            value_bytes = value.encode()
            (<cpp_source.PointSource *> self._inst).particle = std_string(<char *> value_bytes)
