#!/usr/bin/env python

"""The :py:mod:`cards` module can be imported as such::

    from pyne.simplesim import cards

The user does not interact with this module. Rather, the user interacts with
the :py:mod:`definition` module. Instances of classes in the :py:mod:`definition` module
contains instances of the classes in this module.  Below is the reference for
the module.

"""

# TODO in the Sphinx documentation, create a developers and a user's version of
# this documentation.
# TODO in the Sphinx documentation, provide a list of all classes without any
# docstrings. autosummary
# autosummary, :nosignatures: pyne.simplesim.cards
# TODO use sphinx domains where possible instead of double single quotes.
# TODO write a development guide next to the usersguide.
# TODO make error messages valuable: give back to the user their input.
# TODO sphinx inheritance diagrams.
# TODO if i make a material card, make sure to revert back to the Material
# object that Anthony has once they implement a comment.
# Maybe I place an underscore for abstract base classes so that the user
# doesn't see them, but I want them to see them...

import abc

import numpy as np

from pyne import material

class ICard(object):
    """This class is not used by the user. Abstract base class for all cards.
    All cards have a name property and a comment() method.

    """
    # This line makes this class an abstract base class (ABC).
    __metaclass__ = abc.ABCMeta
    # Useful for cell cards and the TMP card in MCNP.
    # k is Boltzmann's constant in eV/K, and T is temperature in K
    kelvin2kT = 8.6173423e-11
    # These are provided in the order in which they are listed in the MCNP
    # manual.
    mcnp_particle = {'neutron': 'N',
                     'anti-neutron': '-N',
                     'photon': 'P',
                     'electron': 'E',
                     'positron': '-E',
                     'muon': '|',
                     'anti-muon': '-|',
                     'tau': '*',
                     'electron neutrino': 'U',
                     'anti-electron neutrino': '-U',
                     'muon neutrino': 'V',
                     'tau neutrino': 'W',
                     'proton': 'H',
                     'anti-proton': '-H',
                     'lambda': 'L',
                     'sigma+': '+',
                     'sigma-': '-',
                     'cascade0': 'X',
                     'cascade-': 'Y',
                     'omega': 'O',
                     'lambdac+': 'C',
                     'cascadec+': '!',
                     'cascadec0': '?',
                     'lambdab0': '<',
                     'pion+': '/',
                     'pion-': '-/',
                     'neutral pion': 'Z',
                     'kaon+': 'K',
                     'kaon-': '-K',
                     'K0 short': '%',
                     'K0 long': '^',
                     'D+': 'G',
                     'D0': '@',
                     'Ds+': 'f',
                     'B+': '>',
                     'B0': 'B',
                     'Bs0': 'Q',
                     'deuteron': 'D',
                     'triton': 'T',
                     'helium-3': 'S',
                     'helium-4': 'A',
                     'heavy ions': '#'
                     }
    
    def __init__(self, name, unique=False):
        """
        Parameters
        ----------
        name : str
            Name of this instance of this card, used for referencing this card
            from others and modifying this card after it has been added to a
            definition (see :py:mod:`pyne.simplesim.definition`). There is no
            known restriction on the characters that can be used in the name,
            but the name must be unique.

        """
        self._unique = unique
        if self._unique:
            self._name = name
        else:
            self.name = name

    # All subclasses must define a comment() method.
    @abc.abstractmethod
    def comment(self):
        """All cards define a comment describing the content of the card."""
        raise NotImplementedError

    #@abc.abstractmethod
    def mcnp(self, float_format, sim):
        # sim is an instance of
        # :py:class:`pyne.simplesim.definition.SimulationDefinition`.
        raise NotImplementedError

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if self._unique:
            raise StandardError("This is a unique card, meaning only one card"
                    " of this type can be found in a ``definition``. "
                    "Accordingly, the name is read-only.")
        if value == '':
            raise ValueError("The ``name`` property of the cell cannot "
                    "be empty.")
        self._name = value

   
class Cell(ICard):
    """A cell is a region of space filled with a material. If requesting a void
    cell, the ``material``, ``density``, and ``density_units`` attributes are
    all set to None (as by default).

    """
    def __init__(self, name, region, material=None, density=None,
                 density_units=None, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region` subclass
            Defines the region of space that this cell occupies (see
            :py:class:`Region`).
        material : :py:class:`pyne.material.Material`, None for void
            A material definition using the :py:mod:`pyne.material` module.
            For use here, the material's :py:attr:`name` property must be set
            to something other than '' and must be unique. See
            :py:class:`pyne.material.Material`.
        density : float, None for void
            Density for the material, in units of density_units.
        density_units : str, None for void
            Either 'g/cm^3', or 'atoms/b/cm'.

        Examples
        --------
        Suppose we have surface cards ``surfA`` and ``surfB``, and material
        cards ``matA`` and ``matB``. The following creates a void cell on the
        negative side of ``surfA`` and the positive side of ``surfB`` (see
        :py:class:`Region` to learn how to create regions)::

            cellA = Cell('A', surfA.neg & surfB.pos)

        The following cell is filled with ``matA``::

            cellB = Cell('B', surfA.neg & surfB.pos, matA, 10.0, 'g/cm^3')
            cellB = Cell('B', surfA.neg & surfB.pos, matA, 0.5, 'atoms/b/cm')

        Note that if a material is specified, a density and density units must
        also be provided The following is not allowed::

            cellB = Cell('B', surfA.neg & surfB.pos, matA)
            cellB = Cell('B', surfA.neg & surfB.pos, matA, density=1)
            cellB = Cell('B', surfA.neg & surfB.pos, matA, 
                    density_units='g/cm^3')
            cellB = Cell('B', surfA.neg & surfB.pos, density=1)
            cellB = Cell('B', surfA.neg & surfB.pos, 
                    density_units='g/cm^3')
        
        """
        # TODO decide how I will do cross-referencing.

        super(Cell, self).__init__(name, *args, **kwargs)
        self.region = region
        self.material = material
        self.density = density
        self.density_units = density_units
        if ((self.material and not (self.density and self.density_units)) or
                (self.density and not (self.density_units and self.material)) or
                (self.density_units and not (self.density and self.material))):
            raise ValueError("If specifying a material, ``material``, "
                    "``density``, and ``density_units`` must all be "
                    "specified.")

    def comment(self):
        string = "Cell '%s':" % self.name
        string += " region %s, " % self.region.comment()
        if self.material and self.density and self.density_units:
            string += "material '%s'" % self.material.name
            string += " density %.5e %s" % (self.density, self.density_units)
        else:
            string += "void"
        return string

    def mcnp(self, float_format, sim):
        # Card number.
        string = "%i" % sim.sys.cell_num(self.name)
        if self.material and self.density and self.density_units:
            # Material number.
            string += " %i " % sim.sys.material_num(self.material.name)
            # Density, with units prefix.
            string += self._mcnp_density_prefix(self.density_units)
            string += formatstr % self.density
        else:
            # Void.
            string += " 0"
        # Print surfaces.
        string += " %s" % self.region.mcnp(sim)
        return string

    def _mcnp_density_prefix(self, density_units):
        """In MCNP, mass densities are represented by prepending a minus sign.

        Parameters
        ----------
        density_units : str
            Must be either 'g/cm^3' or 'atoms/b/cm'.

        Returns
        -------
        prefix : str
            Returns a minus sign for mass densities and an empty string for
            number densities.

        """
        if density_units == 'g/cm^3':
            return '-'
        elif density_units == 'atoms/b/cm':
            return ''
        else:
            raise Exception("Invalid string to specify density units: %s."
                    % density_units)

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, obj):
        self._region = obj

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, obj):
        # This check is redundant.
        if obj and type(obj) is not material.Material:
            raise ValueError("The property ``material`` must be instance "
                    "of ``pyne.material.Material``. User provided {0}.".format(
                    obj))
        if obj and obj.name == '':
            raise ValueError("The ``name`` property of the material cannot "
                    "be empty.")
        self._material = obj

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value

    @property
    def density_units(self):
        return self._density_units

    @density_units.setter
    def density_units(self, value):
        if (value and value != 'g/cm^3' and value != 'atoms/b/cm'):
            raise ValueError("The property ``density_units`` must be either "
                    "'g/cm^3' or 'atoms/b/cm'. User provided "
                    "'{0}'".format(value))
        self._density_units = value


class CellMCNP(Cell):
    """A cell card with keyword options that are available in MCNP. Thus, it
    only makes sense to use this card if writing an input for MCNP. A number of
    the keyword arguments are for a particular particle. The particles
    available are given in :py:attr:`mcnp_particle`. The user provides the full
    name of the particle, as given as keys in :py:attr:`mcnp_particle`. The
    card will then use the appropriate particle designator when writing the
    card.
    
    The U, LAT, and FILL keywords are not available; as this functionality
        should be obtained by using Universe and Lattice cards.

    Note this card was written with MCNPX version 2.7 in mind.

    .. inheritance-diagram:: pyne.simplesim.cards.CellMCNP

    """
    # TODO Sphinx documentation should not list all keyword arguments.
    # TODO let the user add their own through keyword arguments...
    # user_string

    def __init__(self, name, region, material=None,
                 density=None, density_units=None,
                 temperature=None, volume=None,
                 importance=None,
                 exp_transform=None,
                 forced_coll=None,
                 weight_win_bound=None,
                 dxtran_contrib=None,
                 photon_weight=None,
                 fission_turnoff=None,
                 det_contrib=None,
                 transformation=None,
                 user_custom=None
                 ):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        region : :py:class:`Region`
            See :py:class:`Cell`.
        material : :py:class:`pyne.material.Material`, None for void
            See :py:class:`Cell`.
        density : float, None for void
            See :py:class:`Cell`.
        density_units : str, None for void
            See :py:class:`Cell`.
        temperature : float, otional [Kelvin]
            Temperature of the cell for thermal treatment, **TMP**.
        volume : float, optional [cm^3]
            Volume of the cell, **VOL**.
        importance : 2-element tuple of str and int, optional
            Particle importance, **IMP**, for variance reduction. The 1st
            element is a particle name (see :py:attr:`mcnp_particle`). The
            2nd (int) is the cell importance. 
            To specify this input for more
            than one particle, provide a list of these tuples.
        exp_transform : 2-element tuple of str and float, optional
            An exponential transform, **EXT**. The 1st element is a particle
            name (see :py:attr:`mcnp_particle`). The 2nd element (float) is the
            parameter 'a' as described in the MCNP manual. This is a form of
            variance reduction. 
            To specify this input for more
            than one particle, provide a list of these tuples.
        forced_coll : 2-element tuple of str and float, optional
            Forced collisions, **FCL**. The 1st element is a particle name (see
            :py:attr:`mcnp_particle`). The 2nd element (float) is forced
            collision parameter, between -1 and 1. See MCNP manual.
            To specify this input for more
            than one particle, provide a list of these tuples.
        weight_win_bound : 3-element tuple of str, int and str/float, optional
            Weight window lower bound, **WWN**. The 1st element is a particle
            name (see :py:attr:`mcnp_paticle`). The 2nd element (int) of the tuple
            is the energy/time index, and the 3rd element is 'kill' to kill
            particles entering the cell, 'nogame', or a lower bound as a
            non-negative float. See MCNP manual.
            To specify this input for more
            than one particle, or energy/time index, provide a list of these tuples.
        dxtran_contrib : 3-element tuple of str, int/None and float, optional
            DXTRAN Contribution, **DXC**. The 1st element is a particle name
            (see :py:attr:`mcnp_particle`). The 2nd element (int) is an index of a
            DXTRAN sphere (on the DXTRANSphere card), or None if this is to
            apply to all DXTRAN spheres. The 3rd element (float) is the probability
            of contribution to the DXTRAN sphere(s). Only for neutrons and
            photons.
            To specify this input for more than one particle or DXTRAN sphere,
            provide a list of these tuples.
        photon_weight : 0, '-inf', or float; optional
            Photon weight, **PWT**. Relative threshold weight of photons that
            are produced at neutron collisions. With a value of 0, one photon
            is generated per collision; with a value of '-inf', no photons are
            generated. In general, this input can be a positive or negative
            float.
        fission_turnoff : str, optional
            Fission turnoff, **NONU**. The allowed values are: 'capture-gamma',
            'real-gamma', and 'capture-nogamma'. See MCNP manual.
        det_contrib : tuple of str and float, optional
            Detector contribution, **PD**. The 1st element (str) is the name of the
            tally obtaining contribution from this cell, and the 2nd element
            (float) is the probability of contribution.
            To specify this input for more than one tally,
            provide a list of these tuples.
        transformation : str or 4-element tuple, optional
            Cell transformation, **TRCL**. If str, it is the name of a
            :py:class:`Transformation` card (**TR**). If tuple, the 1st element
            is a displacement vector as a 3-element list or
            :py:class:`np.array`; the 2nd element is a rotation matrix as a 3 x
            3 list, :py:class:`np.array`, or :py:class:`np.matrix`; the 3rd
            element is the string 'aux-in-main' or 'main-in-aux'; the 4th
            element is a bool, 'cosines' if rotation matrix is in cosines and
            'degrees' if rotation matrix is in degrees. See
            :py:class:`Transformation`.
        user_custom : str, optional
            This string is appended to the end of the mcnp card. It is possible
            (likely) that the user will want to append a string to the end of
            the MCNP card, given the limited support of keyword
            arguments. This is perhaps most useful if the user wants to specify
            a keyword for a particle that is not supported by any of the
            kwargs.

        Examples
        --------
        TODO

        """
        # TODO allow user to use defalt arguments on the transformation card
        # (e.g. not require the specification of cosines/degrees).
        # TODO allow use of U, LAT, and FILL keywords?
        super(CellMCNP, self).__init__(name, region, material, density,
                                       density_units)
        # Assign keyword arguments.
        self.temperature = temperature
        self.volume = volume
        self.importance = importance
        self.exp_transform = exp_transform
        self.forced_coll = forced_coll
        self.weight_win_bound = weight_win_bound
        self.dxtran_contrib = dxtran_contrib
        self.photon_weight = photon_weight
        self.fission_turnoff = fission_turnoff
        self.det_contrib = det_contrib
        self.transformation = transformation

    def comment(self):
        string = super(CellMCNP, self).comment()
        float_format = "%.5e"
        # temperature
        if self.temperature:
            string += "TMP="
            string += float_format % self.temperature
            string += " K"
        # volume
        if self.volume:
            string += "VOL="
            string += float_format % self.volume
            string += " cm^3"
        # importance
        if self.importance:
            importance = self._make_list(self.importance)
            for entry in importance:
                string += (" IMP:%s=%i" % 
                        (self.mcnp_particle(entry[0]), entry[1]))
        # exp_transform
        if self.exp_transform:
            exp_transform = self._make_list(self.exp_transform)
            for entry in exp_transform:
                string += " EXT:%s=" % self.mcnp_particle(entry[0])
                string += float_format % entry[1]
        # forced_coll
        if self.forced_coll:
            forced_coll = self._make_list(self.forced_coll)
            for entry in forced_coll:
                string += " FCL:%s=" % self.mcnp_particle(entry[0])
                string += float_format % entry[1]
        # weight_win_bound
        if self.neutron_weight_win_bound:
            string += " WWN%i:N=" % self.neutron_weight_win_bound[0]
            string += str(self.neutron_weight_win_bound[1])
        if self.weight_win_bound:
            weight_win_bound = self._make_list(self.weight_win_bound)
            for entry in weight_win_bound:
                string += (" WWN%i:%s=" % 
                        (entry[1], self.mcnp_particle(entry[0])))
                if type(entry[2]) is str: string += entry[2]
                else: string += float_format % entry[2]
        # dxtran_contrib
        if self.dxtran_contrib:
            dxtran_contrib = self._make_list(self.dxtran_contrib)
            for entry in dxtran_contrib:
                string += (" DXC%i:%s=" % 
                    (entry[1], self.mcnp_particle(entry[0])))
                string += float_format % entry[2]
        # photon_weight
        if self.photon_weight: string += " PWT={0}".format(self.photon_weight)
        # fission_turnoff
        if self.fission_turnoff: string += " NONU=%s" % self.fission_turnoff
        # det_contrib
        if self.det_contrib:
            det_contrib = self._make_list(self.det_contrib)
            for entry in det_contrib:
                string += " PD for tally '%s'=" % entry[0]
                string += float_format % entry[1]
        # transform
        if self.transformation:
            # For brevity.
            transform = self.transformation
            string += " TRCL "
            if type(transform) is str:
                string += "'%s'" % transform
            else:
                tempcard = Transformation('temp', transform[0],
                    transform[1], transform[3], transform[4])
                string += tempcard._comment_data()
        # user_custom
        if self.user_custom: string += " and user's custom input."
        string += "."
        return string

    def mcnp(self, float_format, sim):
        string = super(CellMCNP, self).mcnp(float_format, sim)
        # temperature
        if self.temperature:
            string += "TMP="
            string += float_format % (self.temperature * self.kelvin2kT)
        # volume
        if self.volume:
            string += "VOL="
            string += float_format % self.volume
        # importance
        if self.importance:
            importance = self._make_list(self.importance)
            for entry in importance:
                string += (" IMP:%s=%i" % 
                        (self.mcnp_particle(entry[0]), entry[1]))
        # exp_transform
        if self.exp_transform:
            exp_transform = self._make_list(self.exp_transform)
            for entry in exp_transform:
                string += " EXT:%s=" % self.mcnp_particle(entry[0])
                string += float_format % entry[1]
        # forced_coll
        if self.forced_coll:
            forced_coll = self._make_list(self.forced_coll)
            for entry in forced_coll:
                string += " FCL:%s=" % self.mcnp_particle(entry[0])
                string += float_format % entry[1]
        # weight_win_bound
        if self.weight_win_bound:
            weight_win_bound = self._make_list(self.weight_win_bound)
            for entry in weight_win_bound:
                string += (" WWN%i:%s=" % 
                        (entry[1], self.mcnp_particle(entry[0])))
                if entry[2] == 'kill': string += '-1'
                elif entry[2] == 'nogame': string =+ '0'
                else: string += float_format % entry[2]
        # dxtran_contrib
        if self.dxtran_contrib:
            dxtran_contrib = self._make_list(self.dxtran_contrib)
            for entry in dxtran_contrib:
                string += (" DXC%i:%s=" % 
                    (entry[1], self.mcnp_particle(entry[0])))
                string += float_format % entry[2]
        # photon_weight
        if self.photon_weight:
            string += " PWT="
            if self.photon_weight == '-inf': string += "-1.0E6"
            else: string += float_format % self.photon_weight
        # fission_turnoff TODO move exception to setter?
        if self.fission_turnoff:
            string += " NONU="
            if self.fission_turnoff == 'capture-gamma': string += "0"
            elif self.fission_turnoff == 'real-gamma': string += "1"
            elif self.fission_turnoff == 'capture-nogamma': string += "2"
            else: 
                raise ValueError("Expected 'capture-gamma', 'real-gamma', "
                        "'capture-nogamma'. User provided '%s'." %
                        self.fission_turnoff)
        # det_contrib
        if self.det_contrib:
            det_contrib = self._make_list(self.det_contrib)
            for entry in det_contrib:
                string += " PD%i=" % sim.tally_num(entry[0])
                string += float_format % entry[1]
        # transform
        if self.transformation:
            # For brevity.
            transform = self.transformation
            string += " "
            if type(transform) is str:
                string += "TRCL=%i" % sim.transformation_num(transform)
            else:
                if transform[4] == 'degrees': string += "*"
                tempcard = Transformation('temp', transform[0],
                    transform[1], transform[3], transform[4])
                string += "TRCL (%s)" % tempcard._mcnp_lineup(float_format)
                
        # user_custom
        if self.user_custom: string += " %s" % self.user_custom
                
    def _make_list(self, arg):
        if arg is list: return arg
        else: return [arg]

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if value is not None:
            if value < 200:
                raise UserWarning("Temperature set as less than 200 K. "
                        "Are you trying to specify temperature in degrees "
                        "Celcius, etc.? User provided %.4f." % value)
            if value < 1:
                raise UserWarning("Temperature set as less than 1 K. "
                        "Are you trying to specify temperature as 'kT'? "
                        "User provided %.4f." % value)
        self._temperature = value

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, value):
        if value is not None and value < 0:
            raise ValueError("The ``volume`` property cannot be negative. "
                "User provided %.4f." % value)

    @property
    def importance(self):
        return self._importance

    @importance.setter
    def importance(self, value):
        self._importance = value

    @property
    def exp_transform(self):
        return self._exp_transform

    @exp_transform.setter
    def exp_transform(self, value):
        self._exp_transform = value

    @property
    def forced_coll(self):
        return self._forced_coll

    @forced_coll.setter
    def forced_coll(self, value):
        self._forced_coll = value

    @property
    def weight_win_bound(self):
        return self._weight_win_bound

    @weight_win_bound.setter
    def weight_win_bound(self, value):
        self._weight_win_bound = value

    @property
    def dxtran_contrib(self):
        return self._dxtran_contrib

    @dxtran_contrib.setter
    def dxtran_contrib(self, value):
        self._dxtran_contrib = value
    
    @property
    def photon_weight(self):
        return self._photon_weight
    
    @photon_weight.setter
    def photon_weight(self, value):
        self._photon_weight = value
    
    @property
    def fission_turnoff(self):
        return self._fission_turnoff
    
    @fission_turnoff.setter
    def fission_turnoff(self, value):
        self._fission_turnoff = value

    @property
    def det_contrib(self):
        return self._det_contrib
    
    @det_contrib.setter
    def det_contrib(self, value):
        self._det_contrib = value

    @property
    def transformation(self):
        return self._transformation
    
    @transformation.setter
    def transformation(self, value):
        self._transformation = value

    @property
    def user_custom(self):
        return self._user_custom
    
    @user_custom.setter
    def user_custom(self, value):
        self._user_custom = value
    

class IUniverse(ICard):
    """This class is not used by the user. Abstract base class for all
    universe cards.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        pass


class UniverseByRegion(IUniverse):
    """

    """
    # TODO
    def __init__(self, name, region):
        pass


class UniverseByLattice(IUniverse):
    """

    """
    # TODO
    def __init__(self, name, lattice):
        pass


class Lattice(ICard):
    """

    """
    # TODO 
    def __init__(self, name, geom, universe):
        pass


class LatticeByArray(ICard):
    """

    """
    # TODO support of 3D arrays.
    def __init(self, name, geom, xindices, yindices, zindices,
               universe_array):
        pass
        

class ISurface(ICard):
    """This class is not used by the user. Abstract base class for all
    surface cards.

    The Surface superclass contains properties to set the surface as reflecting
    or white. For codes other than MCNPX, reflecting or white surfaces may be
    specified on a separate boundary condition card (i.e. in Serpent) or may
    not even be available. For other codes, then, the appropriate :py:mod:`inputfile`
    class needs to pick up this information and print the appropriate string to
    the code's input file, or in the latter case return an exception.

    .. inheritance-diagram:: pyne.simplesim.cards.ISurface

    """
    # TODO support rotation.
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, reflecting, white, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        reflecting : bool, optional
            The surface has a reflective boundary condition.
        white : bool, optional
            The surface has a white boundary condition (reflection with a
            cosine distribution with respect to the surface normal).

        """
        super(ISurface, self).__init__(name, *args, **kwargs)
        self.reflecting = reflecting
        self.white = white 
        if self.reflecting and self.white:
            raise ValueError("The user set the surface to be reflecting AND "
                    "white, but can only be neither or one of the two.")

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError

    @abc.abstractmethod
    def shift(self, vector):
        """Translates the surface. This is an abstract method, and must be
        defined by each surface card. Shifts may not be permitted for all
        surfaces in all directions, and in such cases an exception is raised.
        
        Parameters
        ----------
        vector : 3-element list or :py:class:`np.array`, float [centimeters]
            The elements specify translation in the x, y, and z directions, 
            in this order.

        Examples
        --------
        Both of the following lines shifts the surface along the x axis by 3 cm::

            surf.shift([3, 0, 0])
            surf.shift(np.array([3, 0, 0]))
        
        """
        raise NotImplementedError

    @abc.abstractmethod
    def stretch(self, vector):
        """Stretches (scales) the surface from the origin. This is an abstract
        method, and must be defined by each surface card. Stretches may not be
        permitted for all surfaces in all directions, and in such cases an
        exception is raised.

        Parameters
        ----------
        vector : 3-element list or :py:class`np.array`, float [unitless]
            The elements specify a stretch in the x, y, and z directions, in
            this order. A zero in any of the directions indicates that no
            stretch is done in that direction. Negative values are allowed, and
            represent reflections.

        Examples
        --------
        Both of the following lines stretch the surface along the y axis by a
        factor of 2. The x and z directions are unaffected::

            surf.stretch([0, 2, 0])
            surf.stretch(np.array([0, 2, 0]))
        
        """
        raise NotImplementedError

    @property
    def neg(self):
        """A property that creates and returns a
        :py:class:`RegionLeaf` that can then be used in
        boolean arithmetic with subclasses of
        :py:class:`Region`. The region is define as the
        space on the side of the surface that has a negative sense.

        In the expected typical usage of the :py:mod:`pyne.simplesim` package,
        regions are constructed using these properties.

        For more information, see :py:class:`Region` and
        :ref:`usersguide_simplesim`.

        Examples
        --------
        The following shows a simple case of how a more complex region can be
        constructed from regions returned by this property::

            reg1 = surf1.neg
            reg2 = surf2.neg
            reg3 = reg1 & reg2
            reg4 = reg1 | reg2

        """
        return RegionLeaf(self, False)

    @property
    def pos(self):
        """Similar to :py:attr:`neg`, except the resulting
        :py:class:`RegionLeaf` is on the side of the surface with a positive
        sense.
        
        """
        return RegionLeaf(self, True)
    
    @property
    def reflecting(self):
        return self._reflecting

    @reflecting.setter
    def reflecting(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("The property ``reflecting`` must be "
                    "None or of boolean type. User provided "
                    "{0}.".format(value))
        self._reflecting = value

    @property
    def white(self):
        return self._white

    @white.setter
    def white(self, value):
        if value is not None and type(value) is not bool:
            raise TypeError("The property ``white`` must be "
                    "None or of boolean type. User provided "
                    "{0}.".format(value))
        self._white = value


class IAxisSurface(ISurface):
    """This class is not used by the user. Abstract base class for all simple
    axis-aligned surfaces. Accordingly, such classes share the cartesian_axis
    property.

    .. inheritance-diagram:: pyne.simplesim.cards.IAxisSurface
    
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, cartesian_axis, reflecting, white, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            Either 'x', 'X', 'y', 'Y', 'z', or 'Z'. Regardless of input, it is
            stored as lower-case. The meaning depends on the surface.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        """
        super(IAxisSurface, self).__init__(name, reflecting, white, *args,
                                           **kwargs)
        self.cartesian_axis = cartesian_axis

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError
    
    @abc.abstractmethod
    def shift(self, vector):
        """See :py:meth:`ISurface.shift`."""
        raise NotImplementedError

    @abc.abstractmethod
    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`."""
        raise NotImplementedError

    @property
    def cartesian_axis(self):
        return self._cartesian_axis

    @cartesian_axis.setter
    def cartesian_axis(self, value):
        if type(value) is not str:
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "a string. User provided {0}.".format(value))
        if (value.lower() != 'x' and value.lower() != 'y' and 
                value.lower() != 'z'):
            raise ValueError("AxisCylinder's cartesian_axis property must be "
                    "'x', 'X', 'y', 'Y', or 'z', 'Z'. "
                    "User provided '{0}'.".format(value))
        self._cartesian_axis = value.lower()


class AxisCylinder(IAxisSurface):
    """Cylinder aligned with and centered on one of the Cartesian axes.
    
    .. inheritance-diagram:: pyne.simplesim.cards.AxisCylinder
    
    """

    # TODO if this is shifted, then it becomes not an axis-cylinder.
    def __init__(self, name, cartesian_axis, radius,
                 reflecting=None, white=None):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            The axis with which the cylinder is aligned and centered.
            See :py:class`IAxisSurface`.
        radius : float [centimeters]
            Radius of the cylinder.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        Examples
        --------
        The following creates a cylinder aligned with the z axis, going through
        centered on the z axis, with a radius of 0.4 cm::

            cyl = AxisCylinder('mycyl', 'z', 0.4)

        """
        super(AxisCylinder, self).__init__(name, cartesian_axis, 
                                           reflecting, white)
        self.radius = radius

    def comment(self):
        return ("Axis cylinder '%s': aligned and centered on %s axis, "
                "with radius %.4f cm (diameter %.4f cm)." %
                        (self.name, self.cartesian_axis,
                        self.radius, 2 * self.radius))
    
    def shift(self, vector):
        """See :py:meth:`ISurface.shift`. Axis cylinders can only be shifted along
        their axis, and even in such cases the shift has no effect. However,
        such a shift must be permitted in case this surface is part of a region
        that is being shifted.

        Examples
        --------
        The following is okay (where we have imported :py:mod:`numpy` as ``np``)::

            cyl = AxisCylinder('mycyl', 'z', 0.4)
            cyl.shift([0, 0, 3])
            cyl.shift(np.array([0, 0, 3]))

        The following do not work:

            cyl.shift([3, 0, 0])
            cyl.shift([0, 3, 3])
            
        
        """
        # Flag for exception.
        iserror = False
        if self.cartesian_axis == 'x' and (vector[1] != 0 or vector[2] != 0):
            iserror = True
            dirs = ('x', 'y', 'z')
        if self.cartesian_axis == 'y' and (vector[0] != 0 or vector[2] != 0):
            iserror = True
            dirs = ('y', 'x', 'z')
        if self.cartesian_axis == 'z' and (vector[0] != 0 or vector[1] != 0):
            iserror = True
            dirs = ('z', 'x', 'y')
        if iserror:
            raise ValueError("A cylinder aligned with the %s axis cannot "
                    "be shifted in the %s or %s directions." % dirs)

    def stretch(self, vector):
        """See :py:meth:`ISurface:stretch`. Axis cylinders can be stretched in
        the direction of their axis, which has no effect (permitted in case
        this surface is part of a region that is being stretched), or can be
        stretched `uniformly` in the plane perpendicular to its axis.
        
        Examples
        --------
        The following stretches are okay for a cylinder aligned with the x axis
        (where we have imported :py:mod`numpy` as ``np``)::
            
            cyl = AxisCylinder('mycyl', 'z', 0.4)
            cyl.stretch([0, 0, 2])
            cyl.stretch([3, 3, 0])
            cyl.stretch(np.array([3, 3, 2]))

        However, the following would cause the cylinder to lose its
        circular cross section, which cannot be accommodated::

            cyl.stretch([0, 3, 0])
            cyl.stretch([2, 3, 1])
        
        """
        # TODO allow some slop between the same two values for a uniform
        # perpendicular stretch.
        # Flag for exception.
        iserror = False
        # 'out' is used in the exception below.
        if self.cartesian_axis == 'x':
            if vector[1] != vector[2]:
                iserror = True
                out = ('y', vector[1], 'z', vector[2], 'x')
            elif vector[1] != 0:
                self.radius *= vector[1]
        if self.cartesian_axis == 'y':
            if vector[0] != vector[2]:
                iserror = True
                out = ('x', vector[0], 'z', vector[2], 'y')
            elif vector[0] != 0:
                self.radius *= vector[0]
        if self.cartesian_axis == 'z':
            if vector[0] != vector[1]:
                iserror = True
                out = ('x', vector[0], 'y', vector[1], 'z')
            elif vector[0] != 0:
                self.radius *= vector[0]
        if iserror:
            raise ValueError("Stretches perpendicular to the axis must be "
                    "uniform in the two perpendicular directions. User "
                    "provided %s stretch %.4f and %s stretch %.4f for a "
                    "%s-aligned cylinder." % out)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value <= 0:
            raise ValueError("The ``radius`` property must be "
                    "positive. User provided %.4f." % value)
        self._radius = value


class AxisPlane(IAxisSurface):
    """Plane perpendicular to one of the Cartesian axes.

    .. inheritance-diagram:: pyne.simplesim.cards.AxisCylinder
    
    """

    def __init__(self, name, cartesian_axis, position,
                 reflecting=None, white=None):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        cartesian_axis : str
            The axis to which the plane is perpendicular.
            See :py:class:`IAxisSurface`.
        position : float [centimeters]
            Position of the plane along :py:attr:`cartesian_axis`.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        Examples
        --------
        The following creates a plane perpendicular to the x axis, 3 cm along
        the positive x axis with a reflecting boundary condition::

           plane = AxisPlane('myplane', 'x', 3, reflecting=True) 

        """
        super(AxisPlane, self).__init__(name, cartesian_axis,
                                    reflecting, white)
        self.position = position
    
    def comment(self):
        return "Axis plane '%s': %s = %.4f cm." % (self.name, self.cartesian_axis,
                self.position)

    def shift(self, vector):
        """See :py:meth:`ISurface.shift`. Axis planes can be shifted in any
        direction, but only shifts along their axis have an effect.

        Examples
        --------
        The following has the effect of shifting the plane's position to x = 6
        cm::

           plane = AxisPlane('myplane', 'x', 3)
           plane.shift([3, 0, 0])

        The following has no effect, but is allowed::

           plane.shift([0, 3, 2])

        """
        if self.cartesian_axis == 'x':
            self.position += vector[0]
        elif self.cartesian_axis == 'y':
            self.position += vector[1]
        elif self.cartesian_axis == 'z':
            self.position += vector[2]

    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`. Axis planes can be stretched in any
        direction, but only stretches along their axis have an effect. The
        position of the plane is scaled by the stretch factor.

        Examples
        --------
        The following has the effect of moving the plane's position to x = 9 cm::

            plane = AxisPlane('myplane', 'x', 3)
            plane.stretch([3, 0, 0])
        
        The following has no effect, but is allowed::

            plane.stretch([0, 3, 2])

        """
        if self.cartesian_axis == 'x' and vector[0] != 0:
            self.position *= vector[0]
        elif self.cartesian_axis == 'y' and vector[1] != 0:
            self.position *= vector[1]
        elif self.cartesian_axis == 'z' and vector[2] != 0:
            self.position *= vector[2]

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value


class IMacrobody(ISurface):
    """This class is not used by the user. Abstract base class for all
    macrobody cards. Macrobodies are an MCNP concept.

    .. inheritance-diagram:: pyne.simplesim.cards.IMacrobody

    """
    __metaclass__ = abc.ABCMeta

    # TODO abstract method for obtaining "sub"-surfaces.
    def __init__(self, name, reflecting, white, *args, **kwargs):
        """

        """
        super(IMacrobody, self).__init__(name, reflecting, white, *args,
                                         **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class Parallelepiped(IMacrobody):
    """Rectangular parallelepiped in which all surfaces are parallel to the
    cartesian axes.

    .. inheritance-diagram::pyne.simplesim.cards.Parallelepiped

    """
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
                 reflecting=False, white=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        xmin, xmax, ymin, ymax, zmin, zmax : float [centimeters]
            Bounds of the parallelepiped in the given direction. The *min value
            must be less than the *max value. Setting both min and max in a
            given direction to 0 indicates the parallelepiped is infinite in
            that direction.
        reflecting : bool, optional
            See :py:class:`ISurface`
        white : bool, optional
            See :py:class:`ISurface`

        Examples
        --------
        The following creates a cube centered at the origin with 4 cm sides::

            pp = Parallelepiped('mypp', -2, 2, -2, 2, -2, 2)

        """
        super(Parallelepiped, self).__init__(name, reflecting, white)
        self.xlims = np.array([xmin, xmax])
        self.ylims = np.array([ymin, ymax])
        self.zlims = np.array([zmin, zmax])

    def comment(self):
        return ("Parallelepiped '%s': [%.4f, %.4f] x [%.4f, %.4f] x " \
               "[%.4f, %.4f] cm." % (self.name, self.xlims[0], self.xlims[1],
                   self.ylims[0], self.ylims[1], self.zlims[0], self.zlims[1]))

    def shift(self, vector):
        """See :py:meth:`ISurface.shift`.
        
        Examples
        --------
        The following::

            pp = Parallelepiped('mypp', -2, 2, -2, 2, -2, 2)
            pp.shift([2, 0, 0])

        creates a parallelepiped bounded by [0, 4] x [-2, 2] x [-2, 2].

        """
        self.xlims += vector[0]
        self.ylims += vector[1]
        self.zlims += vector[2]

    def stretch(self, vector):
        """See :py:meth:`ISurface.stretch`. Handling reflections (negative
        stretch factors) requires additional consideration for this surface,
        but is implemented.
        
        Examples
        --------
        The following::

            pp = Parallelepiped('mypp', 0, 4, -2, 2, -2, 2)
            pp.stretch([2, 0, 3])

        creates a parallelepiped bounded by [0, 8] x [-2, 2] x [-6, 6].
        Consider the reflection of the following parallelepiped
        about the z axis::
        
            pp = Parallelepiped('mypp', 0, 4, -2, 2, -3, 6)
            pp.stretch([0, 0, -1])

        This results in bounds of [0, 4] x [-2, 2] x [-6, 3].
        
        """
        if vector[0] != 0:
            if vector[0] > 0:
                self.xlims *= vector[0]
            else:
                # Stretch factor is negative, swap limits.
                self.xlims = vector[0] * self.xlims[::-1]
        if vector[1] != 0:
            if vector[1] > 0:
                self.ylims *= vector[1]
            else:
                self.ylims = vector[1] * self.ylims[::-1]
        if vector[2] != 0:
            if vector[2] > 0: 
                self.zlims *= vector[2]
            else:
                self.zlims = vector[2] * self.zlims[::-1]

    @property
    def xlims(self):
        return self._xlims

    @xlims.setter
    def xlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of xmin, %.4f, is greater than "
                    "that of xmax, %.4f." % (value[0], value[1]))
        self._xlims = value

    @property
    def ylims(self):
        return self._ylims

    @ylims.setter
    def ylims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of ymin, %.4f, is greater than "
                    "that of ymax, %.4f." % (value[0], value[1]))
        self._ylims = value

    @property
    def zlims(self):
        return self._zlims

    @zlims.setter
    def zlims(self, value):
        if value[0] > value[1]:
            raise ValueError("The value of zmin, %.4f, is greater than "
                    "that of zmax, %.4f." % (value[0], value[1]))
        self._zlims = value


class Cuboid(Parallelepiped):
    """Same exact thing as a :py:class:`Parallelepiped`. This class is provided
    because the name is shorter, and thus may be preferred by those who fancy
    brevity.

    """
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
                 reflecting=False, white=False):
        super(Cuboid, self).__init__(name, xmin, xmax, ymin, ymax, zmin, zmax,
                                     reflecting, white)


class IRegion(ICard):
    """This class is not used by the user. Abstract base class for
    all regions.

    Represents a volume (space) confined by unions and intersections of
    surfaces."""
    # TODO transformation functions
    # Cell cards are then formed by a region and a material.

    # TODO Complement functionality can be added by overloading the
    # __not__ operator and defining a complement boolean property that is set
    # by the __not__ operator.
    # TODO add transformation methods.
    # TODO describe how parent works. it is not actually needed...
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(IRegion, self).__init__(name, *args, **kwargs)
        self.parent = None

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError

    def set(self, region):
        # TODO copy?
        if issubclass(region, RegionLeaf):
            self = RegionLeaf(region.surface, region.pos_sense, self.name)
        elif issubclass(region, RegionOr):
            self = RegionOr(region.left_child, region.right_child, self.name)
        elif issubclass(region, RegionAnd):
            self = RegionAnd(region.left_child, region.right_child, self.name)
        self.parent = region.parent

    def __and__(self, arg):
        return self.intersect(arg)

    def __or__(self, arg):
        return self.union(arg)

    def intersect(self, arg):
        return RegionAnd(self, arg)

    def union(self, arg):
        return RegionOr(self, arg)

    def shift(self, vector):
        """The surfaces themselves are modified; copies are not made.

        """
        # TODO walk.

    def stretch(self, vector):
        """

        """
        # TODO walk.

    def walk(self, leaf_func, and_func=None, or_func=None):
        """

        """
        # TODO
        if isinstance(self, RegionLeaf):
            leaf_func.im_func(leaf_func.im_self, self)
        else:
            self.left_child.walk(leaf_func)
            if and_func and isinstance(self, and_func):
                and_func.im_func(and_func.im_self, self)
            if or_func and isinstance(self, or_func):
                or_func.im_func(or_func.im_func, self)
            self.right_child.walk(leaf_func)
        
    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value


class IRegionBool(IRegion):
    """This class is not used by the user. Abstract base class for
    :py:class:`RegionAnd` and :py:class:`RegionOr`.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, left_child, right_child, name='<Empty>', *args, **kwargs):
        super(IRegionBool, self).__init__(name, *args, **kwargs)
        self.left_child = left_child
        self.right_child = right_child
        self.left_child.parent = self
        self.right_child.parent = self

    @abc.abstractmethod
    def comment(self, midchar):
        return ("(" + self.left_child.comment() + " " + midchar + " " +
                self.right_child.comment() + ")")

    @abc.abstractmethod
    def mcnp(self, sim, midchar):
        return ("(" + self.left_child.mcnp(sim) + midchar + 
                self.right_child.mcnp(sim) + ")")

    @property
    def left_child(self):
        return self._left_child

    @left_child.setter
    def left_child(self, value):
        self._left_child = value

    @property
    def right_child(self):
        return self._right_child

    @right_child.setter
    def right_child(self, value):
        self._right_child = value


class RegionAnd(IRegionBool):
    """

    """
    def comment(self):
        return super(RegionAnd, self).comment('&')

    def mcnp(self, sim):
        return super(RegionAnd, self).mcnp(sim, ' ')


class RegionOr(IRegionBool):
    """

    """
    def comment(self):
        return super(RegionOr, self).comment('|')

    def mcnp(self, sim):
        return super(RegionOr, self).mcnp(sim, ':')


class RegionLeaf(IRegion):
    """
    """

    def __init__(self, surface, pos_sense, name='<Empty>'):
        # TODO Default name is an empty string.
        super(RegionLeaf, self).__init__(name)
        self.surface = surface
        self.pos_sense = pos_sense

    def comment(self):
        if self.pos_sense:
            prefix = '+'
        else:
            prefix = '-'
        return prefix + self.surface.name

    def mcnp(self, sim):
        if self.pos_sense:
            prefix = ''
        else:
            prefix = '-'
        return "%s%i" % (prefix, sim.sys.surface_num(self.surface.name))

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, value):
        self._surface = value

    @property
    def pos_sense(self):
        return self._pos_sense

    @pos_sense.setter
    def pos_sense(self, value):
        # TODO this is probably not okay by the proponents of duck typing.
        if type(value) is not bool:
            raise TypeError("User provided a value for pos_sense that is "
                    "not of boolean type.")
        self._pos_sense = value


class IMisc(ICard):
    """ """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(IMisc, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class ISource(ICard):
    """ """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, *args, **kwargs):
        super(ISource, self).__init__(name, *args, **kwargs)

    @abc.abstractmethod
    def comment(self):
        raise NotImplementedError


class Criticality(ISource):
    """A criticality source of neutrons. Unique card with name `criticality`.
    In MCNP, this is the **KCODE** card.
    
    .. inheritance-diagram:: pyne.simplesim.cards.Criticality

    """

    # TODO in the example, include the resulting MCNP output?
    def __init__(self, n_histories=1000, keff_guess=1.0,
            n_skip_cycles=30, n_cycles=130):
        """
        Parameters
        ----------
        n_histories : int, optional
            Number of particle histories to run in each cycle.
        keff_guess : float, optional
            Initial guess for the effective multiplication constant of the
            system.
        n_skip_cycles : int, optional
            The number of cycles to skip.
        n_cycles : int, optional
            The total number of cycles to simulate (skipped + active).

        Examples
        --------
        The following::

            critsrc = Criticality(2000, 1.5, 30, 300)

        creates a criticality source with 2000 histories per cycle, an initial
        k_eff guess of 1.5, 30 skipped cyles, and 300 total cycles.
        
        """
        super(Criticality, self).__init__('criticality', unique=True) 
        self.n_histories = n_histories
        self.keff_guess = keff_guess
        self.n_skip_cycles = n_skip_cycles
        self.n_cycles = n_cycles

    def comment(self):
        return ("Criticality source '%s': n_histories: %i, keff_guess: %.4f"
                ", n_skip_cycles: %i, n_cycles: %i." % (self.name,
                    self.n_histories, self.keff_guess, self.n_skip_cycles,
                    self.n_cycles))

    @property
    def n_histories(self):
        return self._n_histories

    @n_histories.setter
    def n_histories(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_histories`` must be an "
                    "integer. User provided %.4f." % value)
        if value <= 0:
            raise ValueError("The property ``n_histories`` must be positive. "
                "User provided %i." % value)
        self._n_histories = value

    @property
    def keff_guess(self):
        return self._keff_guess
 
    @keff_guess.setter
    def keff_guess(self, value):
        if value < 0:
            raise ValueError("The property ``keff_guess`` must be "
                    "non-negative. User provided %.4f." % value)
        self._keff_guess = value

    @property
    def n_skip_cycles(self):
        return self._n_skip_cycles

    @n_skip_cycles.setter
    def n_skip_cycles(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_skip_cycles`` must be an "
                    "integer. User provided %.4f." % value)
        if value <= 0:
            raise ValueError("The property ``n_skip_cycles`` must be positive. "
                "User provided %i." % value)
        self._n_skip_cycles = value

    @property
    def n_cycles(self):
        return self._n_cycles

    @n_cycles.setter
    def n_cycles(self, value):
        if type(value) is not int:
            raise ValueError("The property ``n_cycles`` must be an "
                    "integer. User provided %.4f." % value)
        if value < self.n_skip_cycles:
            raise ValueError("The property ``n_cycles`` must be equal to or "
                    "greater than ``n_skip_cycles``. User provided %i." %
                    value)
        # If n_cycles is greater or equal to n_skip_cycles, it is positive.
        self._n_cycles = value


class CriticalityPoints(ISource):
    """Initial source points for neutrons generated by a criticality source.
    Unique card with name 'criticalitypoints'. In MCNP, this is the **KSRC**
    card.

    .. inheritance-diagram:: pyne.simplesim.cards.CriticalityPoints
    
    """
    def __init__(self, points=[[0, 0, 0]]):
        """
        Parameters
        ----------
        points : list of 3-element lists, optional [centimeters]
            A list of 3-element lists (or numpy arrays) specifying
            initial neutron source points in 3-D space.

        Examples
        --------
        The following specifies two initial source points at (1, 2, 3) cm and
        at (3.141..., 2.718..., 0) cm, where we have imported :py:mod:`numpy`
        as ``np``::

            critpts = CriticalityPoints([[1, 2, 3],
                                        np.array([np.pi, np.e, 0])])

        """
        super(CriticalityPoints, self).__init__('criticalitypoints',
                unique=True)
        self.points = points

    def comment(self):
        string = "Criticality points 'criticalitypoints': "
        counter = 0
        for point in self.points:
            counter += 1
            string += "(%.4f, %.4f, %.4f)" % tuple(point)
            if counter < len(self.points):
                string += ", "
            else:
                string += "."
        return string

    @property
    def points(self):
        return self._points

    @points.setter
    def points(self, value):
        for point in value:
            if len(point) is not 3:
                raise ValueError("Length of all point lists/arrays must be 3."
                        " User provided a point %s." % str(point))
        self._points = value


class ITally(ICard):
    """This class is not used by the user. Abstract base class for
    tallies (MCNP) or detectors (Serpent).

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`. Used for, e.g., tally multiplier cards.
        particle : str
            Either 'neutron', 'photon', electron', or 'proton'.

        """
        super(ITally, self).__init__(name, *args, **kwargs)
        self.particle = particle

    @abc.abstractmethod
    def comment(self, title):
        string = "%s tally '%s' of " % (title, self.name)
        if type(self.particle) is not list:
            string += self.particle
            if self.particle != 'all':
                string += "s"
        else:
            pcounter = 0
            for part in self.particle:
                pcounter += 1
                string += "%ss" % part
                if pcounter < len(self.particle):
                    string += ", "
        string += ": "
        return string

    @property
    def particle(self):
        return self._particle

    @particle.setter
    def particle(self, value):
        if (value != 'neutron' and value != 'photon' and 
                value != 'electron' and value != 'proton'):
            raise ValueError("The property ``particle`` must be 'neutron', "
                    "'photon', 'electron', or 'proton'. "
                    "User provided '%s'." % value)
        self._particle = value


class ICellSurfTally(ITally):
    """This class is not used by the user. Abstract base class for
    tallies over cells and surfaces, as opposed to detector tallies. In MCNP,
    these are the **F1**, **F2**, **F4**, **F6**, **F7** and **F8** tallies.
    Some of these are for cells, and some are for surfaces.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, cards, alt_units=False, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : :py:class:`Cell` or :py:class:`ISurface`, list, list of lists
            If tallying 1 cell/surface, the input is that cell/surface card. If
            tallying multiple cells/surfaces, the individual cell/surface cards
            are provided in a list. To obtain the average tally across multiple
            cells/surfaces, these cell/surface cards are provided in their own
            list, within the outer list. To avoid ambiguity, if only one set of
            averages is desired, then this set must be nested in two lists. See
            the examples.
        alt_units : bool, optional
            If set to True and the tally can use alternative units, alternative
            units are used for the tally. See subclasses.

        Examples
        --------
        The following gives the tally in cell A::

            tally = CellFlux('fuel', 'neutron', cellA)

        The following two cards give the tally in surface A and B, and
        the average across surfaces B and C::

            tally = SurfaceFlux('fuel', 'photon', [surfA, surfB, [surfB,
                    surfC]], average=True)

        To obtain the average across surface A and B only, a nested list is
        required to distinguish the case of tallying on A and B individually::

            tally = SurfaceFlux('fuel', 'neutron', [[surfA, surfB]])

        """
        super(ICellSurfTally, self).__init__(name, particle, *args, **kwargs)
        self.cards = cards
        self.alt_units = alt_units

    @abc.abstractmethod
    def comment(self, title, union_type, card_type):
        # card_type is either 'cell', or 'surface'
        # Assuming the user has provided objects of the appropriate type; the
        # issubclass check was not working with little effort. TODO
        string = super(ICellSurfTally, self).comment(title)

        if card_type == 'cell':
            classcheck = Cell
        elif card_type == 'surface':
            classcheck = ISurface
        if type(self.cards) is not list: # issubclass(self.cards, classcheck):
            string += "%s '%s'" % (card_type, self.cards.name)
        elif type(self.cards) is list:
            if type(self.cards[0]) is not list:
                string += "%ss " % card_type
            outcounter = 0
            for obj in self.cards:
                outcounter += 1
                if type(obj) is not list: # issubclass(obj, classcheck):
                    string += "'%s'" % obj.name
                elif type(obj) is list:
                    string += "%s in " % union_type
                    incounter = 0
                    for avgobj in obj:
                        incounter += 1
                        # Must be a cell/surface card.
                        string += "'%s'" % avgobj.name
                        if incounter < len(obj):
                            string += ", "
                # TODO an anti-duck-typing exception:
                #else:
                #    raise ValueError("Expected {0} or list, got {1}.".format(
                #            card_type, type(obj)))
                if outcounter < len(self.cards):
                    string += "; "
        # TODO an anti-duck-typing exception:
        #else:
        #    raise ValueError("Expected {0} or list, got {1}.".format(
        #            card_type, type(self.cards)))
        return string

    def _unique_card_list(self):
        # Returns a unique list of all the cards provided in self.cards.
        # This method is called by
        # :py:class:`pyne.simplesim.SimulationDefinition` for error-checking.
        # TODO this is ideally recursive, and maybe can be implemented in a
        # cleaner way.
        # Assuming the user has provided objects of the appropriate type; the
        # issubclass check was not working with little effort. TODO. See
        # comment().
        if type(self.cards) is not list: # issubclass(self.cards, ?):
            return [self.cards]
        elif type(self.cards) is list:
            card_list = list()
            for obj in self.cards:
                if type(obj) is list:
                    for avgobj in obj:
                        # issubclass(avgobj, IAverageTally)
                        if (type(avgobj) is not list and
                                avgobj not in card_list):
                            card_list += [avgobj]
                elif obj not in card_list:
                    # Not a list, must be a cell/surface.
                    card_list += [obj]
            return card_list
        else:
            raise ValueError("Expected cell, surface, or list,"
                " got {0}.".format(type(self.cards)))

    @property
    def cards(self):
        return self._cards

    @cards.setter
    def cards(self, value):
        self._cards = value

    @property
    def alt_units(self):
        return self._alt_units

    @alt_units.setter
    def alt_units(self, value):
        self._alt_units = value


class SurfaceCurrent(ICellSurfTally):
    """Surface current tally. In MCNP, this is the **F1** card.

    .. inheritance-diagram:: pyne.simplesim.cards.SurfaceCurrent

    """
    def __init__(self, name, particle, cards, total=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : :py:class:`ISurface`, list, list of lists
            See :py:class:`ICellSurfTally`.
        total : bool, optional
            Include a tally for the total current across all surfaces
            specified on this card (note: NOT across all surfaces in the
            `problem`).
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.

        Examples
        --------
        The following requests the tally in surface A, surface B, as well as the
        total across A and B::

            tally = SurfaceCurrent('fuel', 'electron', [surfA, surfB],
                    total=True)

        In the following, the tally is also weighted by particle energy::

            tally = SurfaceCurrent('fuel', 'photon', [[surfA, surfB]],
                    alt_units=True)

        """
        super(SurfaceCurrent, self).__init__(name, particle, cards, alt_units)
        self.total = total

    def comment(self):
        string = super(SurfaceCurrent, self).comment("Surface current", 'total',
                'surface')
        if self.total:
            string += "; and total of all provided."
        else:
            string += "."
        return string


class IAverageTally(ICellSurfTally):
    """This class is not used by the user. Abstract base class for
    tallies of averaged quantities. In MCNP, these are the **F2**, **F4**,
    **F6**, **F7** and **F8** tallies. Some of these are for cells, and some
    are for surfaces.

    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, particle, cards, average=False, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : :py:class:`Cell` or :py:class:`ISurface`, list, list of lists
            See :py:class:`ICellSurfTally`.
        average : bool, optional
            Include a tally for the average flux across all cells/surfaces
            specified on this card (note: NOT across all cells in the
            `problem`).
        alt_units : bool, optional
            If set to True and the tally can use alternative units, alternative
            units are used for the tally. See subclasses.

        Examples
        --------
        The following requests the tally in cell A, cell B, as well as the
        average across A and B::

            tally = CellEnergyDeposition('fuel', 'neutron', [cellA, cellB],
                    average=True)

        """
        super(IAverageTally, self).__init__(name, particle, cards, *args,
                                            **kwargs)
        self.average = average

    @abc.abstractmethod
    def comment(self, title, card_type):
        avgstr = 'avg.'
        string = super(IAverageTally, self).comment(title, avgstr, card_type)
        if self.average:
            string += "; and %s of all provided." % avgstr
        else:
            string += "."
        return string

    @property
    def average(self):
        return self._average

    @average.setter
    def average(self, value):
        self._average = value


class SurfaceFlux(IAverageTally):
    """Surface flux tally. In MCNP, this is the **F2** card.

    .. inheritance-diagram:: pyne.simplesim.cards.SurfaceFlux

    """
    def __init__(self, name, particle, cards, average=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : :py:class:`ISurface`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.

        Examples
        --------
        The following requests the tally in surface A, surface B, as well as
        the average across A and B::

            tally = SurfaceFlux('fuel', 'electron', [surfA, surfB],
                    average=True)

        In the following, the tally is also weighted by particle energy::

            tally = SurfaceFlux('fuel', 'proton', [[surfA, surfB]],
                    alt_units=True)
        
        See base classes for more examples.

        """
        super(SurfaceFlux, self).__init__(name, particle, cards, average,
                                          alt_units)

    def comment(self):
        return super(SurfaceFlux, self).comment("Surface flux", 'surface')

    @property
    def total(self):
        return self._total

    @total.setter
    def total(self, value):
        self._total = value


class CellFlux(IAverageTally):
    """Cell flux tally. In MCNP, this is the **F4** card.

    .. inheritance-diagram:: pyne.simplesim.cards.CellFlux

    """
    def __init__(self, name, particle, cards, average=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particle : str
            See :py:class:`ITally`.
        cards : :py:class:`Cell`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, Tally is additionally weighted by particle energy.

        Examples
        --------
        The following requests the tally in cell A, cell B, as well as the
        average across A and B::

            tally = CellFlux('fuel', 'electron', [cellA, cellB],
                    average=True)

        In the following, the tally is also weighted by particle energy::

            tally = CellFlux('fuel', 'proton', [[cellA, cellB]],
                    alt_units=True)
        
        See base classes for more examples.

        """
        super(CellFlux, self).__init__(name, particle, cards, average,
                                       alt_units)

    def comment(self):
        return super(CellFlux, self).comment("Cell flux", 'cell')


class CellEnergyDeposition(IAverageTally):
    """Energy deposition tally. In MCNP, this is the **F6** card. In MCNP, it
    is not permitted to use a particle 'all' and also to use alternative units.

    .. inheritance-diagram:: pyne.simplesim.cards.CellEnergyDeposition

    """
    # TODO in mcnp input, prevent particle all and alt_units
    def __init__(self, name, particles, cards, average=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. For this tally, the user can specify the
            particle type as a list of strs to tally more than one type of
            particle. Also, the additional value of 'all' is allowed, and
            specifies collision heating. As is expected, 'all' cannot be
            provided as part of a list.
        cards : :py:class:`Cell`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are MeV/g and the alternative units are jerks/g.

        Examples
        --------
        The following requests the energy deposited by neutrons in cell A::

            tally = CellEnergyDeposition('energy', 'neutron', cellA)


        The following requests the energy deposited by neutrons and protons in
        cell A::

            tally = CellEnergyDeposition('energy', ['neutron', 'proton'],
                    cellA)
        
        The following requests the energy deposited by all particles in cell
        A::

            tally = CellEnergyDeposition('energy', 'all', cellA)

        The following are not allowed in MCNP::
            
            tally = CellEnergyDeposition('energy', ['neutron', 'all'], cellA)
            tally = CellEnergyDeposition('energy', 'all', cellA, alt_units=True)

        See base classes for more examples.

        """
        super(CellEnergyDeposition, self).__init__(name, particles, cards,
                average, alt_units)
        # TODO move this error check to the MCNP method.
        if self.particle == 'all' and self.alt_units:
            raise ValueError("The particle cannot be 'all' if alt_units is "
                    "True.")

    def comment(self):
        return super(CellEnergyDeposition, self).comment("Energy deposition",
                'cell')

    @property
    def particle(self):
        return self._particle

    @particle.setter
    def particle(self, value):
        if type(value) is list:
            for string in value:
                if (string != 'neutron' and string != 'photon' and 
                        string != 'electron' and string != 'proton'):
                    raise ValueError("The ``particle`` list must "
                            "contain only 'neutron', 'photon', 'electron',"
                            " or 'proton'. User provided '%s'." % string)
        else:
            # A single string is provided.
            if (value != 'neutron' and value != 'photon' and 
                    value != 'electron' and value != 'proton' and 
                    value != 'all'):
                raise ValueError("The property ``particle`` must be "
                        "'neutron', 'photon', 'electron', 'proton', or 'all'."
                        "User provided '%s'." % value)
        self._particle = value


class CellFissionEnergyDeposition(IAverageTally):        
    """Fission energy deposition tally. In MCNP, this is the **F7** card. The
    particle is necessarily neutron.

    .. inheritance-diagram:: pyne.simplesim.cards.CellFissionEnergyDeposition

    """
    # TODO prevent user from specifying a different particle.
    def __init__(self, name, cards, average=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        cards : :py:class:`Cell`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are MeV/g and the alternative units are jerks/g.

        Examples
        --------
        The following requests the tally in cell A, cell B, as well as the
        average across A and B::

            tally = CellFissionEnergyDeposition('fuel', [cellA,
                    cellB], average=True)

        In the following, the alternate units are used::

            tally = CellFissionEnergyDeposition('fuel', [[cellA,
                    cellB]], alt_units=True)
        
        See base classes for more examples.

        """
        super(CellFissionEnergyDeposition, self).__init__(name, 'neutron',
                cards, average, alt_units)

    def comment(self):
        return super(CellFissionEnergyDeposition, self).comment(
                "Fission energy deposition", 'cell')


class CellPulseHeight(IAverageTally):
    """Pulse height tally in cells. In MCNP, this is the **F8** card. For a
    charge deposition tally, see :py:class:`CellChargeDeposition`.

    .. inheritance-diagram:: pyne.simplesim.cards.CellPulseHeight

    """
    def __init__(self, name, particles, cards, average=False, alt_units=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. Multiple particles can be provided in a
            list of str. In MCNP, if only 'proton', or 'electron' is
            specified, both are automatically included.
        cards : :py:class:`Cell`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.
        alt_units : bool, optional
            If True, alternative units are used for the tally. In MCNP, the
            default units are pulses and the alternative units are MeV.

        Examples
        --------
        The following requests the tally in cell A and cell B for both protons
        and electrons, and requests units of MeV::

            tally = CellPulseHeight('fuel', ['proton', 'electron'], [cellA,
                    cellB], alt_units=True)

        See base classes for more examples.

        """
        super(CellPulseHeight, self).__init__(name, particles, cards, average,
                                              alt_units)

    def comment(self):
        return super(CellPulseHeight, self).comment("Pulse height", 'cell')

    @property
    def particle(self):
        return self._particle

    @particle.setter
    def particle(self, value):
        # Copied from CellEnergyDeposition.particle
        if type(value) is list:
            for string in value:
                self._assert_particle(string)
        else:
            # A single string is provided.
            self._assert_particle(string)
        self._particle = value

    def _assert_particle(self, string):
        if (string != 'neutron' and string != 'photon' and 
                string != 'electron' and string != 'proton'):
            raise ValueError("The ``particle`` list must "
                    "contain only 'neutron', 'photon', 'electron',"
                    " or 'proton'. User provided '%s'." % string)

class CellChargeDeposition(CellPulseHeight):
    """Charge deposition tally in cells. In MCNP, this is the **+F8** card. No
    alternative units are available.

    .. inheritance-diagram:: pyne.simplesim.cards.CellChargeDeposition

    """
    # TODO it doesn't make sense that the user can provide the particle type
    # here, at least for MCNP.
    def __init__(self, name, particles, cards, average=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. Multiple particles can be provided in a
            list of str. In MCNP, if only 'proton', or 'electron' is
            specified, both are automatically included.
        cards : :py:class:`Cell`, list, list of lists
            See :py:class:`IAverageTally`.
        average : bool, optional
            See :py:class:`IAverageTally`.

        Examples
        --------
        The following requests the tally in cell A and cell B for both protons
        and electrons::

            tally = CellChargeDeposition('fuel', ['proton', 'electron'], [cellA,
                    cellB])

        See base classes for more examples.

        """
        super(CellPulseHeight, self).__init__(name, particles, cards, average)

    def comment(self):
        return super(CellPulseHeight, self).comment("Charge deposition", 'cell')


class RepeatedStructure(IAverageTally):
    # TODO
    pass


class IDetector(ITally):
    def __init__(self, name, particle, spec, sep_direct=True, *args, **kwargs):
        super(IDetector, self).__init__(name, particle, *args, **kwargs)
        self.spec = spec
        self.sep_direct = sep_direct

    @abc.abstractmethod
    def comment(self, name):
        string = super(IDetector, self).comment(name)
        if type(self.spec) is tuple:
            string += self._tuple_tostring(self.spec)
        else:
            counter = 0
            for point in self.spec:
                counter += 1
                string += self._tuple_tostring(point)
                if counter < len(self.spec):
                    string += "; "
        if not self.sep_direct: 
            dircontr = 'not '
        else:
            dircontr = ''
        return string + "; direct contrib is %sseparate." % dircontr
    
    @abc.abstractmethod
    def _tuple_tostring(self):
        raise NotImplementedError

    @property
    def spec(self):
        return self._spec

    @spec.setter
    def spec(self, value):
        self._spec = value

    @property
    def sep_direct(self):
        return self._sep_direct

    @sep_direct.setter
    def sep_direct(self, value):
        self._sep_direct = value


class PointDetector(IDetector):
    """A point detector tally. In MCNP, this is the **F5** card. This is not to
    be confused with the more general use of the term `Detector` in Serpent.

    .. inheritance-diagram:: pyne.simplesim.cards.PointDetector

    """
    # TODO ideally *args would be used to let the user specify any number of
    # points.
    # TODO I wish we could avoid the use of negative numbers to signal
    # somethign semantic other than a negative number, but other alternatives
    # here seem to not be as clean or easy or general.
    def __init__(self, name, particle, spec, sep_direct=True):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. In MCNP for this tally, only neutrons and
            photons are allowed.
        spec : tuple, list of tuples [centimeters/mean free paths]
            The tuple has 2 elements: a 3-element list of floats and a float.
            The 3-element list provides the location of the point detector, and
            the float is the radius of a sphere of exclusion. The list can also
            be a :py:mod:`numpy` array. By default, the units for the radius is
            also centimeters, but can be changed to mean free paths by
            providing a negative radius. If requesting multiple point
            detectors, a list of point-radius tuples can be provided.
        sep_direct : bool, optional
            In MCNP, the direct contribution to the tally is printed
            separately. Set to False to disable the separate printing. This is
            a property of the undocumented :py:class:`IDetector`.

        Examples
        --------
        The following creates a single point detector at the origin, without a
        sphere of exclusion::

            det = PointDetector('point', 'neutron', ([0, 0, 0], 0))

        The following creates a detector at (1, 1, 1) cm with a sphere of
        exclusion with a radius of 1 cm, where we have imported :py:mod:`numpy`
        as ``np``)::

            det = PointDetector('point', 'neutron', (np.array([1, 1, 1]), 1))

        The radius for the sphere of exclusion here is 3 mfp::

            det = PointDetector('point', 'neutron', ([1, 0, 0], -3))

        This is an example of requesting two point detectors::
        
            det = PointDetector('point', 'photon', [([0, 0, 0],  0),
                                                     ([1, 0, 0], -3)])

        Here, it is requested that the direct contribution is not tallied
        separately::

            det = PointDetector('point', 'photon', ([0, 0, 0], 0),
                    sep_direct=False)

        """
        super(PointDetector, self).__init__(name, particle, spec, sep_direct)

    def comment(self):
        return super(PointDetector, self).comment("Point detector")

    def _tuple_tostring(self, apoint):
        numbertuple = tuple(apoint[0]) + (abs(apoint[1]),)
        string = "point (%.4f, %.4f, %.4f) cm, radius %.4f " % numbertuple
        if apoint[1] < 0:
            string += 'mfp'
        else:
            string += 'cm'
        return string


class RingDetector(IDetector):
    """A ring detector tally. In MCNP, this is the **F5a** card. This is not to
    be confused with the more general use of the term `Detector` in Serpent.

    .. inheritance-diagram:: pyne.simplesim.cards.RingDetector

    """
    # TODO use *args instead of these silly lists.
    def __init__(self, name, particle, spec, sep_direct=True):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ITally`.
        particles : str, list of str
            See :py:class:`ITally`. In MCNP for this tally, only neutrons and
            photons are allowed.
        spec : tuple, list of tuples [centimeters/mean free paths]
            The tuple has 4 elements: a Cartesian axis string ('x', 'y', 'z'),
            a position (float) along that axis, the radius (float) of the ring, and the
            radius (float) of the sphere of exclusion. A negative radius for the sphere
            changes the units to mean free paths. To request multiple ring
            detectors, a list of these tuples can be provided. The Cartesian
            axis strings can be upper or lower case ('x', 'X', 'y', 'Y', 'z',
            'Z').
        sep_direct : bool, optional
            In MCNP, the direct contribution to the tally is printed
            separately. Set to False to disable the separate printing. This is
            a property of the undocumented :py:class:`IDetector`.

        Examples
        --------
        The following creates a single ring detector at x = 10.0 cm, with a
        2.0 cm radius, and a 1.0 cm radius sphere of exclusion::

            det = RingDetector('ring', 'neutron', ('x', 10.0, 2.0,  1.0))

        In the following, the sphere of exclusion has a radius of 1.0 mfp::

            det = RingDetector('ring', 'neutron', ('x', 10.0, 2.0, -1.0))

        This is an example of requesting two ring detectors::

            det = RingDetector('ring', 'neutron', [('x', 10.0, 2.0, -1.0),
                                                   ('y', 20.0, 3.0, 1.0)])

        Here it is requested that the direct contribution is not tallied
        separately::
        
            det = RingDetector('ring', 'neutron', ('x', 10.0, 2.0, -1.0), 
                    sep_direct=False)

        """
        super(RingDetector, self).__init__(name, particle, spec, sep_direct)

    def comment(self):
        return super(RingDetector, self).comment("Ring detector")

    def _tuple_tostring(self, aring):
        string = ("ring %s = %.4f cm, radius %.4f cm, s.o.e. "
                "radius %.4f " %
                (aring[0], aring[1], aring[2], abs(aring[3])))
        if aring[3] < 0:
            string += 'mfp'
        else:
            string += 'cm'
        return string


class EnergyGrid(IMisc):
    """Energy grid for tallies. In MCNP, this is the **E** card.

    .. inheritance-diagram:: pyne.simplesim.cards.EnergyGrid

    """
    def __init__(self, name, tally, energies):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        tally : :py:class:`ITally`, None
            The tally for which this energy grid should apply. If requesting
            for this grid to apply to all tallies, then this is None.
        energies : list, :py:class:`np.array`
            The upper bounds of the energy groups.

        Examples
        --------
        """
        super(EnergyGrid, self).__init__(name)
        self.tally = tally
        self.energies = energies

    def comment(self):
        string = "Energy grid '%s' for " % self.name
        if self.tally is None:
            string += "all tallies"
        else:
            string += "tally %s" % self.tally.name
        return string + ": %i groups." % len(self.energies)

    @property
    def energies(self):
        return self._energies

    @energies.setter
    def energies(self, value):
        self._energies = value


class Comment(ITally):
    pass


class Mode(IMisc):
    pass


class NeutronPhysics(IMisc):
    pass


class PhotonPhysics(IMisc):
    pass


class ElectronPhysics(IMisc):
    pass


class ProtonPhysics(IMisc):
    pass


class Transformation(IMisc):
    """A coordinate transformation. In MCNP, this is the **TR** card. MCNP
    allows the specification of less than 9 values in the rotation matrix; this
    is not supported here.

    .. inheritance-diagram:: pyne.simplesim.cards.Transformation

    """
    # TODO support for MCNP's less-than-9-element transformation matrices.
    def __init__(self, name, displacement, rotation, aux_in_main=True,
            degrees=False):
        """
        Parameters
        ----------
        name : str
            See :py:class:`ICard`.
        displacement : 3-element list, :py:class:`np.array`, float [cm]
            Displacement vector.
        rotation : 3 x 3 list, :py:class:`np.array`, :py:class:`np.matrix`
            Rotation matrix.
        aux_in_main : bool, optional
            The displacement vector is in the main coordinate system and points
            to the origin of the transformed/auxiliary coordinate system. If
            False, then the displacement vector is in the
            transformed/auxiliary coordiante system and points to the origin of
            the main coordinate system.
        degrees : bool, optional
            If True, then the elements of the rotation matrix are in degrees.

        """
        super(Transformation, self).__init__(name)
        self.displacement = displacement
        self.rotation = rotation
        self.aux_in_main = aux_in_main

    def comment(self):
        string = "Transformation '%s': pos. of " % self.name
        string += self._comment_data()
        return string

    def _comment_data(self):
        if self.aux_in_main:
            string += "aux origin in main"
        else:
            string += "main origin in aux"
        string += " (%.5e, %.5e, %.5e) cm" % tuple(self.displacement)
        dirs = ['x', 'y', 'z']
        for idx in range(3):
            string += " %s' <%.5e, %.5e, %.5e>" % tuple(dirs[idx],
                    self.rotation[0][0],
                    self.rotation[1][0], 
                    self.rotation[2][0])
        if self.degrees: string += " deg"
        string += "."
        return string

    def mcnp(self, float_format, sim):
        string = ""
        if self.degrees: string += "*"        
        string += "TR%i" % sim.transformation_num(self.name)
        string += self._mcnp_lineup(float_format)
        return string

    def _mcnp_lineup(self, float_format): 
        # Needed by CellMCNP.
        formatstr = " %s %s %s" % 3 * (float_format,)
        string += formatstr % tuple(self.displacement)
        for i_row in range(3):
            for i_col in range(3):
                string += " "
                string += float_format % self.rotation[i_row][i_col]
        if self.aux_in_main: string += "1"
        else: string += "-1"


        return string

    @property
    def displacement(self):
        return self._displacement
    
    @displacement.setter
    def displacement(self, value):
        self._displacement = value

    @property
    def rotation(self):
        return self._rotation
    
    @rotation.setter
    def rotation(self, value):
        self._rotation = value

    @property
    def aux_in_main(self):
        return self._aux_in_main
    
    @aux_in_main.setter
    def aux_in_main(self, value):
        self._aux_in_main = value
