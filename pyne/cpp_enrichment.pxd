"""Cython header for enrichment library."""
from pyne cimport std
from pyne cimport cpp_material

cdef extern from "enrichment.h" namespace "pyne::enrichment":

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

    Cascade _fill_default_uranium_cascade() except +
    extern Cascade default_uranium_cascade

    double prod_per_feed(double, double, double) except +
    double tail_per_feed(double, double, double) except +
    double tail_per_prod(double, double, double) except +

    double alphastar_i(double, double, double) except +

    void _recompute_nm(Cascade &, double) except +
    void _recompute_prod_tail_mats(Cascade &) except +
    Cascade _norm_comp_secant(Cascade &, double, int) except +
    double _deltaU_i_OverG(Cascade &, int) except +

    Cascade ltot_per_feed(Cascade &, double, int) except +
    Cascade multicomponent(Cascade &, double, int) except +
