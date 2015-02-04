"""Cython header for enrichment library."""
from libcpp.string cimport string as std_string

from pyne cimport cpp_material

cdef extern from "enrichment_cascade.h" namespace "pyne::enrichment":

    cdef cppclass Cascade:
        # Constructors
        EnrichmentParameters() except +

        # Attributes
        double alpha
        double Mstar

        int j
        int k

        double N
        double M

        double x_feed_j
        double x_prod_j
        double x_tail_j

        cpp_material.Material mat_feed
        cpp_material.Material mat_prod
        cpp_material.Material mat_tail

        double l_t_per_feed
        double swu_per_feed
        double swu_per_prod

        void _reset_xjs() except +

cdef extern from "enrichment_symbolic.h" namespace "pyne::enrichment":

    Cascade solve_symbolic(Cascade &) except +

cdef extern from "enrichment.h" namespace "pyne::enrichment":

    Cascade _fill_default_uranium_cascade() except +
    extern Cascade default_uranium_cascade

    double feed_per_prod(double, double, double) except +
    double feed_per_tail(double, double, double) except +
    double prod_per_feed(double, double, double) except +
    double prod_per_tail(double, double, double) except +
    double tail_per_feed(double, double, double) except +
    double tail_per_prod(double, double, double) except +
    double value_func(double) except +
    double swu_per_feed(double, double, double) except +
    double swu_per_prod(double, double, double) except +
    double swu_per_tail(double, double, double) except +
    
    double alphastar_i(double, double, double) except +

    void _recompute_nm(Cascade &, double) except +
    void _recompute_prod_tail_mats(Cascade &) except +
    Cascade _norm_comp_secant(Cascade &, double, int) except +
    double _deltaU_i_OverG(Cascade &, int) except +

    Cascade solve_numeric(Cascade &) except +
    Cascade solve_numeric(Cascade &, double) except +
    Cascade solve_numeric(Cascade &, double, int) except +

    Cascade multicomponent(Cascade &, char *) except +
    Cascade multicomponent(Cascade &, char *, double) except +
    Cascade multicomponent(Cascade &, char *, double, int) except +
    Cascade multicomponent(Cascade &) except +
    Cascade multicomponent(Cascade &, std_string) except +
    Cascade multicomponent(Cascade &, std_string, double) except +
    Cascade multicomponent(Cascade &, std_string, double, int) except +
