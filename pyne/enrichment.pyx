"""Python wrapper for enrichment."""
# Cython imports
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

from pyne cimport std

from pyne cimport nucname
from pyne import nucname

from pyne cimport stlconverters as conv

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_enrichment



#####################
### Cascade Class ###
#####################


cdef class Cascade:
    """This class is a container for enrichment cascade parameters which 
    defines the perfomance of a separations plant. Instances of this class 
    are passed into and out of many enrichment functions.  
    """

    def __cinit__(self):
        self._inst = new cpp_enrichment.Cascade()

    def __dealloc__(self):
        del self._inst

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
        """This is an integer in zzaaam-form that represents the jth key component.
        This nuclide is preferentially enriched in the product stream. For standard 
        uranium cascades j is 922350 (ie U-235).
        """
        def __get__(self):
            return self._inst.j

        def __set__(self, value):
            self._inst.j = nucname.zzaaam(value)

    property k:
        """This is an integer in zzaaam-form that represents the kth key component.
        This nuclide is preferentially enriched in the waste stream. For standard 
        uranium cascades k is 922380 (ie U-238).
        """
        def __get__(self):
            return self._inst.k

        def __set__(self, value):
            self._inst.k = nucname.zzaaam(value)

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
        waste stream mat_tail.  The :math:`x^T_j` value should be 
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
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self._inst.mat_feed
            return pymat

        def __set__(self, mat):
            cdef pyne.material._Material pymat
            if isinstance(mat, pyne.material._Material):
                pymat = mat
            else:
                pymat = pyne.material.Material(mat)
            self._inst.mat_feed = <pyne.cpp_material.Material> pymat.mat_pointer[0]

    property mat_prod:
        """Product (enriched) material.
        """
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self._inst.mat_prod
            return pymat

        def __set__(self, mat):
            cdef pyne.material._Material pymat
            if isinstance(mat, pyne.material._Material):
                pymat = mat
            else:
                pymat = pyne.material.Material(mat)
            self._inst.mat_prod = <pyne.cpp_material.Material> pymat.mat_pointer[0]

    property mat_tail:
        """Tails (de-enriched) material.
        """
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self._inst.mat_tail
            return pymat

        def __set__(self, mat):
            cdef pyne.material._Material pymat
            if isinstance(mat, pyne.material._Material):
                pymat = mat
            else:
                pymat = pyne.material.Material(mat)
            self._inst.mat_tail = <pyne.cpp_material.Material> pymat.mat_pointer[0]

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
        """Sets the x_feed_j, x_prod_j, and x_tail_j attributes to their
        values in the mat_feed, mat_prod, and mat_tail materials.
        """
        self._inst._reset_xjs()

def default_uranium_cascade():
    """Returns a copy of a default uranium enrichment cascade, which has 
    sensible initial values for this very common case.

    The values of this instance of Cascade are as follows::

        duc = pyne.enrichment.Cascade()

        duc.alpha = 1.05
        duc.Mstar = 236.5

        duc.j = 922350
        duc.k = 922380

        duc.N = 30.0
        duc.M = 10.0

        duc.x_feed_j = 0.0072
        duc.x_prod_j = 0.05
        duc.x_tail_j = 0.0025

        duc.mat_feed = pyne.material.Material({922340: 0.000055, 922350: 0.00720, 
                                    922380: 0.992745}, 1.0, "Natural Uranium", 1.0)
        duc.mat_prod = pyne.material.Material()
        duc.mat_tail = pyne.material.Material()

        duc.l_t_per_feed = 0.0
        duc.swu_per_feed = 0.0
        duc.swu_per_prod = 0.0

    Returns
    -------
    duc : Cascade
        As defined above.

    """
    cdef cpp_enrichment.Cascade cpp_duc = cpp_enrichment._fill_default_uranium_cascade()
    cdef Cascade duc = Cascade()
    duc._inst[0] = cpp_duc
    return duc



cdef class Enrichment(fccomp.FCComp):
    """Enrichment Fuel Cycle Component Class.  Daughter of FCComp.

    Parameters
    ----------
    enrich_params : EnrichmentParameters, optional 
        This specifies how the enrichment cascade should be set up.  It is a EnrichmentParameters
        instance.  If enrich_params is not specified, then the cascade is initialized with values 
        from uranium_enrichment_defaults().
    name : str 
        The name of the enrichment fuel cycle component instance.

    """

    def __cinit__(self, *args, **kwargs):
        pass

    def __init__(self, enrich_params=None, char * name="", *args, **kwargs):
        cdef EnrichmentParameters enr_par

        if enrich_params is None:
            self._inst = new cpp_enrichment.Enrichment(std.string(name))
        elif isinstance(enrich_params, EnrichmentParameters):
            enr_par = enrich_params
            self._inst = new cpp_enrichment.Enrichment(<cpp_enrichment.EnrichmentParameters> enr_par.ptr[0], std.string(name))


    #
    # Class Attributes
    #

    # Enrichment Attributes

    property alpha_0:
        """The :math:`\\alpha_0` attribute specifies the overall stage separation factor
        for the cascade.  This should be set on initialization.  Values should be
        greater than one.  Values less than one represent de-enrichment."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).alpha_0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).alpha_0 = <double> value


    property Mstar_0:
        """The :math:`M^*_0` represents a first guess at what the `Mstar` should be.
        The value of Mstar_0 on initialization should be in the ballpark
        of the optimized result of the Mstar attribute.  However, :math:`M^*_0` must
        always have a value between the weights of the j and k key components."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).Mstar_0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).Mstar_0 = <double> value


    property Mstar:
        """The :math:`M^*` attribute represents the mass for which the adjusted
        stage separation factor, :math:`\\alpha^*_i`, is equal to one.  It is this
        value that is varied to achieve an optimized enrichment cascade."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).Mstar

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).Mstar = <double> value


    property mat_tail:
        """In addition to the mat_feed and mat_prod materials, Enrichment
        also has a tails or waste stream that is represented by this attribute.
        The mass of this material and the ms_prod product material should always 
        add up to the mass of the mat_feed feed stock."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).mat_tail
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_enrichment.Enrichment *> self._inst).mat_tail = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property j:
        """This is an integer in zzaaam-form that represents the jth key component.
        This nuclide is preferentially enriched in the product stream.
        For standard uranium cascades j is 922350 (ie U-235)."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).j = nucname.zzaaam(value)


    property k:
        """This is an integer in zzaaam-form that represents the kth key component.
        This nuclide is preferentially enriched in the waste stream.
        For standard uranium cascades k is 922380 (ie U-238)."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).k

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).k = nucname.zzaaam(value)


    property x_prod_j:
        """This is the target enrichment of the jth isotope in the
        product stream mat_prod.  The :math:`x^P_j` value is set by 
        the user at initialization or run-time.  For typical uranium 
        vectors, this value is about U-235 = 0.05."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).x_prod_j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).x_prod_j = <double> value


    property x_tail_j:
        """This is the target enrichment of the jth isotope in the
        waste stream ms_tail.  The :math:`x^W_j` value is set by the 
        user at initialization or runtime.  For typical uranium vectors,
        this value is about U-235 = 0.0025."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).x_tail_j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).x_tail_j = <double> value


    property N:
        """This is the number of enriching stages present in an ideal cascade.
        Along with Mstar and M, this number is optimized to ensure that a product 
        enrichment of x_prod_j is attained."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).N

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).N = <double> value


    property M:
        """This is the number of stripping stages present in an ideal cascade.
        Along with Mstar and N, this number is optimized to ensure that a waste 
        enrichment of x_tail_j is attained."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).M

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).M = <double> value


    property N0:
        """This is the number of enriching stages initially guessed by the user."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).N0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).N0 = <double> value


    property M0:
        """This is the number of stripping stages initially guessed by the user."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).M0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).M0 = <double> value


    property TotalPerFeed:
        """This represents the total flow rate of the cascade divided by the 
        feed flow rate.  As such, it shows the mass of material needed in the
        cascade to enrich an additional kilogram of feed.  Symbolically,
        the total flow rate is given as :math:`L` while the feed rate is
        :math:`F`.  Therefore, this quantity is sometimes seen as 'L-over-F'
        or as 'L/F'.  TotalPerFeed is the value that is minimized to form an 
        optimized cascade."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).TotalPerFeed

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).TotalPerFeed = <double> value


    property SWUperFeed:
        """This value denotes the number of separative work units (SWU) required
        per kg of feed for the specified cascade."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).SWUperFeed

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).SWUperFeed = <double> value


    property SWUperProduct:
        """This value is the number of separative work units (SWU) required
        to produce 1 [kg] of product in the specified cascade."""
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).SWUperProduct

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).SWUperProduct = <double> value



    #
    # Class Methods
    # 

    def initialize(self, EnrichmentParameters enrich_params):
        """The initialize method takes an enrichment parameter object and sets
        the corresponding Enrichment attributes to the same value.

        Parameters
        ----------
        enrich_params : EnrichmentParameters
            A class containing the values to (re-)initialize an Enrichment cascade with.

        """
        cdef EnrichmentParameters enr_par = enrich_params
        (<cpp_enrichment.Enrichment *> self._inst).initialize(<cpp_enrichment.EnrichmentParameters> enr_par.ptr[0])


    def calc_params(self):
        """This sets the Enrichment parameters to the following 
        values::

                self.params_prior_calc["MassFeed"] = self.mat_feed.mass
                self.params_after_calc["MassFeed"] = 0.0

                self.params_prior_calc["MassProduct"] = 0.0
                self.params_after_calc["MassProduct"] = self.mat_prod.mass

                self.params_prior_calc["MassTails"] = 0.0
                self.params_after_calc["MassTails"] = self.mat_tail.mass

                self.params_prior_calc["N"] = self.N
                self.params_after_calc["N"] = self.N

                self.params_prior_calc["M"] = self.M
                self.params_after_calc["M"] = self.M

                self.params_prior_calc["Mstar"] = self.Mstar
                self.params_after_calc["Mstar"] = self.Mstar

                self.params_prior_calc["TotalPerFeed"] = self.TotalPerFeed
                self.params_after_calc["TotalPerFeed"] = self.TotalPerFeed

                self.params_prior_calc["SWUperFeed"] = self.SWUperFeed
                self.params_after_calc["SWUperFeed"] = 0.0

                self.params_prior_calc["SWUperProduct"] = 0.0
                self.params_after_calc["SWUperProduct"] = self.SWUperProduct

        """
        (<cpp_enrichment.FCComp *> self._inst).calc_params()


    def calc(self, input=None):
        """This method performs an optimization calculation on M* and solves for 
        appropriate values for all Enrichment attributes.  This includes the 
        product and waste streams flowing out of the the cascade as well.

        Parameters
        ----------
        input : dict or Material or None, optional
            If input is present, it is set as the component's mat_feed.  If input is 
            a nuclide mapping (zzaaam keys, float values), it is first converted into a 
            Material before being set as mat_feed.

        Returns
        -------
        output : Material
            mat_prod

        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_enrichment.FCComp *> self._inst).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output


    def PoverF(self, double x_F, double x_P, double x_W):
        """Solves for the product over feed enrichment ratio.

        .. math::

            \\frac{p}{f} = \\frac{(x_F - x_W)}{(x_P - x_W)}

        Parameters
        ----------
        x_F : float
            Feed enrichment.
        x_P : float
            Product enrichment.
        x_W : float
            Waste enrichment.

        Returns
        -------
        pfratio : float
            As calculated above.

        """
        return (<cpp_enrichment.Enrichment *> self._inst).PoverF(x_F, x_P, x_W)


    def WoverF(self, double x_F, double x_P, double x_W):
        """Solves for the waste over feed enrichment ratio.

        .. math::

            \\frac{p}{f} = \\frac{(x_F - x_P)}{(x_W - x_P)}

        Parameters
        ----------
        x_F : float
            Feed enrichment.
        x_P : float
            Product enrichment.
        x_W : float
            Waste enrichment.

        Returns
        -------
        wfratio : float
            As calculated above.

        """
        return (<cpp_enrichment.Enrichment *> self._inst).WoverF(x_F, x_P, x_W)

