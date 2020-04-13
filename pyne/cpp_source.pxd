from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as cstr
from libcpp.vector cimport vector

cdef extern from "source.h" namespace "pyne":

    cdef cppclass PointSource:
        # constructors
        PointSource() except +
        PointSource(double) except +
        PointSource(double, double) except +
        PointSource(double, double, double) except +
        PointSource(double, double, double, double) except +
        PointSource(double, double, double, double, double) except +
        PointSource(double, double, double, double, double, double) except +
        PointSource(double, double, double, double, double, double, double) except +
        PointSource(double, double, double, double, double, double, double, cstr) except +
        PointSource(double, double, double, double, double, double, double, cstr, double) except +

        # attributes
        double x
        double y
        double z
        double u
        double v
        double w
        double E
        cstr particle
        double weight

        # methods
        cstr mcnp() except +
        cstr mcnp(int) except +
        pass
