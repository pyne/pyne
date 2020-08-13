"""The enrichment module contains tools for defining and manipulating 
enrichment cascades.  The Cascade class is a simple container for storing 
parameters that define an enrichment setup.  These include feed, product, 
and tail materials, target enrichments, and separation factors.  The main 
functions in this module compute the total flow rate and separation factors
from an initial cascade.  Other helper functions compute relative flow rates 
and nuclide-specific separation factors.
"""
from __future__ import unicode_literals

# Cython imports
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free
from libcpp.string cimport string as std_string

from warnings import warn
from pyne.utils import QA_warn

from pyne cimport nucname
from pyne import nucname
from pyne cimport stlcontainers as conv
cimport pyne.cpp_material
cimport pyne.material
import pyne.material
from pyne cimport cpp_enrichment


QA_warn(__name__)


#####################
### Cascade Class ###
#####################


cdef class Cascade:
    """This class is a container for enrichment cascade parameters that 
    define the perfomance of a separations plant. Instances of this class 
    are passed into and out of many enrichment functions.  
    """

    def __cinit__(self, **kwargs):
        self._inst = new cpp_enrichment.Cascade()

    def __init__(self, **kwargs):
        """__init__(self, **kwargs)

        Parameters
        ----------
        kwargs : optional
            Any keyword argument which is supplied is applied as an attribute
            to this instance.
        """
        self._mat_feed = None
        self._mat_prod = None
        self._mat_tail = None
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __dealloc__(self):
        del self._inst

    def __repr__(self):
        attrs = [a for a in dir(self) if not a.startswith('_')]
        attr_reps = [a + '=' + repr(getattr(self, a)) for a in attrs]
        r = self.__class__.__name__ + '(' + ', '.join(attr_reps) + ')'
        return r

    #
    # Class Attributes
    #
    property alpha:
        """The :math:`\\alpha` attribute specifies the overall stage separation factor
        for the cascade.  This should be set on initialization.  Values should be
        greater than one.  Values less than one represent de-enrichment."""
        def __get__(self):
            return self._inst.alpha

        def __set__(self, value):
            self._inst.alpha = <double> value

    property Mstar:
        """This is the mass separation factor :math:`M^*`.  On initialization, this 
        should be in the ballpark of the optimized result of the Mstar value.  However, 
        this must always have a value between the weights of the j and k key components.
        """
        def __get__(self):
            return self._inst.Mstar

        def __set__(self, value):
            self._inst.Mstar = <double> value

    property j:
        """This is an integer in id-form that represents the jth key component.
        This nuclide is preferentially enriched in the product stream. For standard 
        uranium cascades j is 922350 (ie U-235).
        """
        def __get__(self):
            return self._inst.j

        def __set__(self, value):
            self._inst.j = nucname.id(value)

    property k:
        """This is an integer in id-form that represents the kth key component.
        This nuclide is preferentially enriched in the tails stream. For standard 
        uranium cascades k is 922380 (ie U-238).
        """
        def __get__(self):
            return self._inst.k

        def __set__(self, value):
            self._inst.k = nucname.id(value)

    property N:
        """The number of enriching stages."""
        def __get__(self):
            return self._inst.N

        def __set__(self, value):
            self._inst.N = <double> value

    property M:
        """The number of stripping stages."""
        def __get__(self):
            return self._inst.M

        def __set__(self, value):
            self._inst.M = <double> value

    property x_feed_j:
        """This is the target enrichment of the jth isotope in the
        feed stream mat_feed.  The :math:`x^F_j` value should be 
        set prior to solving for the remainder of the cascade.  For 
        typical uranium vectors, this value is about U-235 = 0.00711.
        """
        def __get__(self):
            return self._inst.x_feed_j

        def __set__(self, value):
            self._inst.x_feed_j = <double> value

    property x_prod_j:
        """This is the target enrichment of the jth isotope in the
        product stream mat_prod.  The :math:`x^P_j` value should be 
        set prior to solving for the remainder of the cascade.  For 
        typical uranium vectors, this value is about U-235 = 0.05.
        """
        def __get__(self):
            return self._inst.x_prod_j

        def __set__(self, value):
            self._inst.x_prod_j = <double> value

    property x_tail_j:
        """This is the target enrichment of the jth isotope in the
        Tails stream mat_tail.  The :math:`x^T_j` value should be 
        set prior to solving for the remainder of the cascade. For 
        typical uranium vectors, this value is about U-235 = 0.0025.
        """
        def __get__(self):
            return self._inst.x_tail_j

        def __set__(self, value):
            self._inst.x_tail_j = <double> value

    property mat_feed:
        """Feed material to be enriched.  Often set at initialization.
        """
        def __get__(self):
            cdef pyne.material._Material mat_feed_proxy
            if self._mat_feed is None:
                mat_feed_proxy = pyne.material.Material(free_mat=False)
                mat_feed_proxy.mat_pointer = &(<cpp_enrichment.Cascade *> self._inst).mat_feed
                self._mat_feed = mat_feed_proxy
            return self._mat_feed
    
        def __set__(self, value):
            cdef pyne.material._Material value_proxy
            value_proxy = pyne.material.Material(value, free_mat=not isinstance(value, pyne.material._Material))
            (<cpp_enrichment.Cascade *> self._inst).mat_feed = value_proxy.mat_pointer[0]
            self._mat_feed = None

    property mat_prod:
        """Product (enriched) material.
        """
        def __get__(self):
            cdef pyne.material._Material mat_prod_proxy
            if self._mat_prod is None:
                mat_prod_proxy = pyne.material.Material(free_mat=False)
                mat_prod_proxy.mat_pointer = &(<cpp_enrichment.Cascade *> self._inst).mat_prod
                self._mat_prod = mat_prod_proxy
            return self._mat_prod
    
        def __set__(self, value):
            cdef pyne.material._Material value_proxy
            value_proxy = pyne.material.Material(value, free_mat=not isinstance(value, pyne.material._Material))
            (<cpp_enrichment.Cascade *> self._inst).mat_prod = value_proxy.mat_pointer[0]
            self._mat_prod = None

    property mat_tail:
        """Tails (de-enriched) material.
        """
        def __get__(self):
            cdef pyne.material._Material mat_tail_proxy
            if self._mat_tail is None:
                mat_tail_proxy = pyne.material.Material(free_mat=False)
                mat_tail_proxy.mat_pointer = &(<cpp_enrichment.Cascade *> self._inst).mat_tail
                self._mat_tail = mat_tail_proxy
            return self._mat_tail
    
        def __set__(self, value):
            cdef pyne.material._Material value_proxy
            value_proxy = pyne.material.Material(value, free_mat=not isinstance(value, pyne.material._Material))
            (<cpp_enrichment.Cascade *> self._inst).mat_tail = value_proxy.mat_pointer[0]
            self._mat_tail = None

    property l_t_per_feed:
        """Total flow rate (:math:`L_t`) per feed flow rate.  This is a 
        characteristic of the cascade as a whole.  As such it is this 
        quatity which is minimized in any real cascade.
        """
        def __get__(self):
            return self._inst.l_t_per_feed

        def __set__(self, value):
            self._inst.l_t_per_feed = <double> value

    property swu_per_feed:
        """The seperative work units (SWU) per unit mass of feed material. 
        """
        def __get__(self):
            return self._inst.swu_per_feed

        def __set__(self, value):
            self._inst.swu_per_feed = <double> value

    property swu_per_prod:
        """The seperative work units (SWU) per unit mass of prod material. 
        """
        def __get__(self):
            return self._inst.swu_per_prod

        def __set__(self, value):
            self._inst.swu_per_prod = <double> value

    # Class methods
    def _reset_xjs(self):
        """Sets the x_feedeed_j, x_prod_j:, and x_tail_j attributes to their
        values in the mat_feed, mat_prod, and mat_tail materials.
        """
        self._inst._reset_xjs()

def default_uranium_cascade():
    """Returns a copy of a default uranium enrichment cascade, which has 
    sensible initial values for this very common case.

    The values of this instance of Cascade are as follows:

    .. code-block:: python

        duc = pyne.enrichment.Cascade(N=30.0, M=10.0, alpha=1.05, Mstar=236.5, 
                j=922350, k=922380, x_feed_j=0.0072, x_prod_j=0.05, x_tail_j=0.0025,
                l_t_per_feed=0.0, swu_per_feed=0.0, swu_per_prod=0.0, 
                mat_feed=pyne.material.Material({922340: 5.5e-05, 922350: 0.0072, 
                                                 922380: 0.992745}, 1.0, 
                                                'Natural Uranium', 1.0), 
                mat_prod=pyne.material.Material({}, -1.0, '', -1.0), 
                mat_tail=pyne.material.Material({}, -1.0, '', -1.0))

    Returns
    -------
    duc : Cascade
        As defined above.

    """
    cdef cpp_enrichment.Cascade cpp_duc = cpp_enrichment._fill_default_uranium_cascade()
    cdef Cascade duc = Cascade()
    duc._inst[0] = cpp_duc
    return duc

def feed(double x_feed, double x_prod, double x_tail, double product=0, 
         double tails=0):
    """feed(x_feed, x_prod, x_tail, product=0, tails=0)
    Calculates the feed quantity in kg from either the product or tails.

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Product enrichment.
    x_tail : float
        Feed enrichment.
    product : float, optional
        Quantity of product in kg
    tails : float, optional
        Quantity of tails in kg

    Returns
    -------
    feed : float
        Feed quantity
    """
    if product > 0:
        return product * cpp_enrichment.feed_per_prod(x_feed, x_prod, x_tail)
    else:
        return tails * cpp_enrichment.feed_per_tail(x_feed, x_prod, x_tail)
    
def product(double x_feed, double x_prod, double x_tail, double feed=0, 
            double tails=0):
    """product(x_feed, x_prod, x_tail, feed=0, tails=0)
    Calculates the product quantity in kg from either the feed or tails.

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Product enrichment.
    x_tail : float
        Product enrichment.
    feed : float, optional
        Quantity of feed in kg
    tails : float, optional
        Quantity of tails in kg

    Returns
    -------
    product : float
        Product quantity
    """
    if feed > 0:
        return feed * cpp_enrichment.prod_per_feed(x_feed, x_prod, x_tail)
    else:
        return tails * cpp_enrichment.prod_per_tail(x_feed, x_prod, x_tail)

def tails(double x_feed, double x_prod, double x_tail, double feed=0, 
          double product=0):
    """tails(x_feed, x_prod, x_tail, feed=0, product=0)
    Calculates the tails quantity in kg from either the feed or product.

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Tails enrichment.
    x_tail : float
        Tails enrichment.
    feed : float, optional
        Quantity of feed in kg
    product : float, optional
        Quantity of product in kg

    Returns
    -------
    tails : float
        Tails quantity
    """
    if feed > 0:
        return feed * cpp_enrichment.tail_per_feed(x_feed, x_prod, x_tail)
    else:
        return product * cpp_enrichment.tail_per_prod(x_feed, x_prod, x_tail)

def value_func(double x):
    """value_func(x)
    Calculates the value or separation potential of an assay.

    .. math::

        V(x) = (2x - 1) \\log{\\frac{x}{x - 1}}

    Parameters
    ----------
    x : float
        assay enrichment.
    
    Returns
    -------
    val : float
        As calculated above.
    """
    return cpp_enrichment.value_func(x)

def swu(double x_feed, double x_prod, double x_tail, double feed=0, 
        double product=0, double tails=0):
    """swu(x_feed, x_prod, x_tail, feed=0, product=0, tails=0)
    Calculates the SWU required to reach a given quantity of an enrichment
    level. One of feed, product, or tails must be provided.

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Product enrichment.
    x_tail : float
        Feed enrichment.
    feed : float, optional
        Quantity of feed in kg
    product : float, optional
        Quantity of product in kg
    tails : float, optional
        Quantity of tails in kg

    Returns
    -------
    SWU : float
        SWU required
    """
    if feed > 0:
        return feed * cpp_enrichment.swu_per_feed(x_feed, x_prod, x_tail)
    elif product > 0:
        return product * cpp_enrichment.swu_per_prod(x_feed, x_prod, x_tail)
    else:
        return tails * cpp_enrichment.swu_per_tail(x_feed, x_prod, x_tail)
        
def prod_per_feed(double x_feed, double x_prod, double x_tail):
    """prod_per_feed(x_feed, x_prod, x_tail)
    Calculates the product over feed enrichment ratio.

    .. math::

        \\frac{p}{f} = \\frac{(x_f - x_t)}{(x_p - x_t)}

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
       Product enrichment.
    x_tail : float
        Tails enrichment.

    Returns
    -------
    pfratio : float
        As calculated above.

    """
    return cpp_enrichment.prod_per_feed(x_feed, x_prod, x_tail)


def tail_per_feed(double x_feed, double x_prod, double x_tail):
    """tail_per_feed(x_feed, x_prod, x_tail)
    Calculates the tails over feed enrichment ratio.

    .. math::

        \\frac{t}{f} = \\frac{(x_f - x_p)}{(x_t - x_p)}

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Product enrichment.
    x_tail : float
        Tails enrichment.

    Returns
    -------
    tfratio : float
        As calculated above.

    """
    return cpp_enrichment.tail_per_feed(x_feed, x_prod, x_tail)


def tail_per_prod(double x_feed, double x_prod, double x_tail):
    """tail_per_prod(x_feed, x_prod, x_tail)
    Calculates the tails over product enrichment ratio.

    .. math::

        \\frac{t}{p} = \\frac{(x_f - x_p)}{(x_t - x_f)}

    Parameters
    ----------
    x_feed : float
        Feed enrichment.
    x_prod : float
        Product enrichment.
    x_tail : float
        Tails enrichment.

    Returns
    -------
    tpratio : float
        As calculated above.

    """
    return cpp_enrichment.tail_per_prod(x_feed, x_prod, x_tail)


def alphastar_i(double alpha, double Mstar, double M_i):
    """alphastar_i(alpha, Mstar, M_i)
    Calculates the stage separation factor for a nuclide i of atomic mass :math:`M_i`.

    .. math::

        \\alpha^*_i = \\alpha^{(M^* - M_i)}

    Parameters
    ----------
    alpha : float
        Stage separation factor.
    Mstar : float
        Mass separation factor.      
    M_i : float
        Atomic mass of the ith nuclide.

    Returns
    -------
    astar_i : float
        As calculated above.

    """
    return cpp_enrichment.alphastar_i(alpha, Mstar, M_i)


def solve_symbolic(Cascade orig_casc):
    """solve_symbolic(orig_casc)
    Computes the cascade parameters based on a given initial state.

    Parameters
    ----------
    orig_casc : Cascade
        A cascade to compute the l_t_per_feed, swu_per_feed, swu_per_prod,
        mat_prod, and mat_tail attributes for.

    Returns
    -------
    casc : Cascade
        A new cascade object, copied from the original, with the appropriate
        attributes computed.

    """
    cdef Cascade casc = Cascade()
    cdef cpp_enrichment.Cascade ccasc = cpp_enrichment.solve_symbolic(orig_casc._inst[0])
    casc._inst[0] = ccasc
    return casc


def solve_numeric(Cascade orig_casc, double tolerance=1.0E-7, int max_iter=100):
    """solve_numeric(orig_casc, tolerance=1.0E-7, max_iter=100)
    Calculates the total flow rate (:math:`L_t`) over the feed flow 
    rate (:math:`F`).

    Parameters
    ----------
    orig_casc : Cascade
        A cascade to compute the l_t_per_feed, swu_per_feed, swu_per_prod,
        mat_prod, and mat_tail attributes for.
    tolerance : float, optional
        Numerical tolerance for solvers, default=1E-7.
    max_iter : int, optional
        Maximum number of iterations for underlying solvers, default=100.

    Returns
    -------
    casc : Cascade
        A new cascade object, copied from the original, with the appropriate
        attributes computed.

    """
    cdef Cascade casc = Cascade()
    cdef cpp_enrichment.Cascade ccasc = cpp_enrichment.solve_numeric(orig_casc._inst[0], tolerance, max_iter)
    casc._inst[0] = ccasc
    return casc


def multicomponent(Cascade orig_casc, solver="symbolic", 
                   double tolerance=1.0E-7, int max_iter=100):
    """multicomponent(orig_casc, solver="symbolic", tolerance=1.0E-7, max_iter=100)
    Calculates the optimal value of Mstar by minimzing the seperative power.
    The minimizing the seperative power is equivelent to minimizing :math:`L_t/F`,
    or the total flow rate for the cascade divided by the feed flow rate. 
    Note that orig_casc.Mstar represents an intial guess at what Mstar might be.
    This function is appropriate for feed materials with more than 2 nuclides 
    (i.e. multicomponent).

    Parameters
    ----------
    orig_casc : Cascade
        A cascade to optimize.
    solver : str, optional
        Flag for underlying cascade solver function to use. Current options 
        are either "symbolic" or "numeric".
    tolerance : float, optional
        Numerical tolerance for underlying solvers, default=1E-7.
    max_iter : int, optional
        Maximum number of iterations for underlying solvers, default=100.

    Returns
    -------
    casc : Cascade
        A new cascade object, copied from the original, which has been optimized
        to minimize flow rates.  Correct values of product and tails materials
        are also computed on this instance.

    """
    cdef char * csolver
    s_bytes = solver.encode('UTF-8')
    csolver = s_bytes
    cdef Cascade casc = Cascade()
    cdef std_string strsolver = std_string(csolver)
    cdef cpp_enrichment.Cascade ccasc = cpp_enrichment.multicomponent(\
                                    orig_casc._inst[0], strsolver, tolerance, max_iter)
    casc._inst[0] = ccasc
    return casc
