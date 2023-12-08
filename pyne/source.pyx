from pyne cimport stlcontainers
from libc.stdlib cimport free
from libcpp.string cimport string as std_string
import numpy as np
cimport numpy as np

from pyne import stlcontainers

from pyne.utils import QA_warn

QA_warn(__name__)

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
    def __cinit__(self, x=0.0, y=0.0, z=0.0, u=0.0, v=0.0, w=0.0, E=14., p="Neutron", weight=1):
        p_bytes = p.encode()
        self._inst = new cpp_source.PointSource(x, y, z, u, v, w, E, std_string(<char *> p_bytes), weight)
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
    
    property u:
        """no docstring for u, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).u

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).u = <double> value
    
    property v:
        """no docstring for j, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).v

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).v = <double> value
    
    property w:
        """no docstring for w, please file a bug report!"""
        def __get__(self):
            return (<cpp_source.PointSource *> self._inst).w

        def __set__(self, value):
            (<cpp_source.PointSource *> self._inst).w = <double> value
    
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
