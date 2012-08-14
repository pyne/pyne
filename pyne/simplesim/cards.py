#!/usr/bin/env python

"""The ``cards`` module can be imported as such::

    from pyne.simulation import cards

The user does not interact with this module. Rather, the user interacts with
the ``definition`` module. Instances of classes in the ``definition`` module
contains instances of the classes in this module.  Below is the reference for
the module.

"""

import abc

from pyne import material

class Card(object):
    """Abstract superclass for all cards. All cards have a name and number
    property, as well as a comment method.
    
    """

    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def __init__(self, name, number):
        name(name)
        number(number)

    @abc.abstractmethod
    def comment(self):
        return

    @property
    def number(self):
        return self._number

    @number.setter
    def number(self, value):
        self._number = value
       
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

class CellVoid(Card):
    """
    """
    def __init__(self, name, number, region):
        """
        """
        self.region(region)

   def comment(self):
       return

   @property
   def region(self):
       return self._region

   
class Cell(CellVoid):
    """
    """

    def __init__(self, name, number, region, material, density, density_units):
        """
        """
        super(Cell, self).__init__(name, number, region)
        self.material(material)
        self.density(density)
        self.density_units(density_units)

    def comment(self):
        return

    @property
    def material(self):
        return self._material

    @property
    def density(self):
        """Density for the material, in units of self.density_units."""
        return self._density

    @property
    def density_units(self):
        """Either 'g/cm^3', or 'atoms/b/cm'."""
        return self._density_units


class CellSimpleVoid(CellVoid):
    """
    """

    def __init__(self, name, number, neg_surfs, pos_surfs):
        """There must be thought given to how default arguments are treated. It
        seems to make more sense that there are no kwargs in this module, that
        they are all in ``definition``, since ``definition`` contains the
        methods that the user calls.

        """
        # TODO this could be an overloaded constructor to CellVoid().
        # Create a region using the negative and positive surfaces.
        super(CellSimpleVoid, self).__init__(name, number, region)
        return

    def comment(self):
        return

#    @property
#    def neg_surface_cards(self):
#        """Surfaces with a negative sense for this cell."""
#        return self._neg_surface_cards
#
#    @property
#    def pos_surface_cards(self):
#        """Surfaces with a positive sense for this cell."""
#        return self._pos_surface_cards

class CellSimple(Cell):
    """
    """
    def __init__(self, name, number, neg_surfs, pos_surfs, material, density,
            density_units):
        # Create region from neg and pos surfaces.
        super(CellSimple, self).__init__(name, number, region,
                                         material, density, density_units)
        return


class CellSimpleVoidMCNP(CellSimpleVoid):
    """
    """

    def __init__(self, name, number, neg_surfs, pos_surfs,
                 temperature=None, volume=None,
                 neutron_imp=None, photon_imp=None, electron_imp=None,
                 proton_imp=None):
        super(CellSimpleVoidMCNP, self).__init__(name, number, neg_surfs,
                pos_surfs)
        # Assign keyword arguments.
        self.temperature(temperature)
        self.volume(volume)
        self.neutron_imp(neutron_imp)
        self.photon_imp(photon_imp)
        self.electron_imp(electron_imp)
        self.proton_imp(proton_imp)

    @property
    def temperature(self):
        """Temperature of the cell, in Kelvin."""
        return self._temperature

    @property
    def volume(self):
        """Volume of the cell, in cm^3. Optional."""
        return self._volume

    @property
    def neutron_imp(self):
        return self._neutron_imp

    @property
    def photon_imp(self):
        return self._photon_imp

    @property
    def electron_imp(self):
        return self._electron_imp

    @property
    def proton_imp(self):
        return self._proton_imp


class CellSimpleMCNP(CellSimpleVoidMCNP, CellSimple):
    """
    """

    def __init__(self, name, number, neg_surfs, pos_surfs, material,
                 temperature=None, volume=None,
                 neutron_imp=None, photon_imp=None, electron_imp=None,
                 proton_imp=None):
        """
        """
        # Based on Python's Method Resolution Order (MRO), the constructor for
        # CellSimpleVoidMCNP is called because it is listed first above.
        super(CellSimpleMCNP, self).__init__(name, number, neg_surfs,
                                             pos_surfs,
                                             temperature=temperature,
                                             volume=volume,
                                             neutron_imp=neutron_imp,
                                             photon_imp=photon_imp,
                                             electron_imp=electron_imp,
                                             proton_imp=proton_imp)
        # The following fields are not initialized via the constructor above.
        self.material(material)
        self.density(density)
        self.density_units(density_units)
        return

    def comment(self):
        return
        

class Surface(Card):
    """
    """

    def __init__(self, name, number, reflecting, white):
        """TODO The Surface superclass contains fields for reflecting and white
        surfaces. For codes other than MCNPX, reflecting or white surfaces may
        be specified on a separate boundary condition card (Serpent) or may not
        even be available. For other codes, then, the appropriate ``inputfile``
        class needs to pick up this information and print the appropriate
        string to the code's input file, or in the latter case return an
        exception.

        """
        return

    def comment(self):
        """This should be abstract."""
        return

    @property
    def reflecting(self):
        return self._reflecting

    @property
    def white(self):
        return self._white


class AxisCylinder(Surface):
    """
    """

    def __init__(self, name, number, cartesian_axis, radius,
                 reflecting, white):
        """
        """
        return

    def comment(self):
        return

    @property
    def radius(self):
        return self._radius

    @property
    def cartesian_axis(self):
        return self._cartesian_axis


class Macrobody(Surface):
    def __init__(self):
        pass

class Parallelepiped(Macrobody):
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
            reflecting=False, white=False):


class Cuboid(Parallelepiped):
    """
    """
    def __init__(self, name, xmin, xmax, ymin, ymax, zmin, zmax,
            reflecting=False, white=False):
        """
        """
        super(Cuboid, self).__init__(name, xmin, xmax, ymin, ymax, zmin, zmax,
                                     reflecting, white)

class Region(Card):
    """Represents a volume (space) confined by unions and intersections of
    surfaces."""
    # Cell cards are then formed by a region and a material.
    def __init__(self, )
        """
        """


class Option(Card):
    """All cards classified as data cards in MCNP, except for materials."""

    def __init__(self, name, number):
        """
        """
        return

    def comment(self):
        """This should be abstract as well."""
        return






