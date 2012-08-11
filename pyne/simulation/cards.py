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

    def __init__(self, name, number, neg_surfs, pos_surfs, neutron_imp=None):
        """There must be thought given to how default arguments are treated. It
        seems to make more sense that there are no kwargs in this module, that
        they are all in ``definition``, since ``definition`` contains the
        methods that the user calls.

        """
        return

    def comment(self):
        return
    
    @property
    def neg_surface_cards(self):
        """Surfaces with a negative sense for this cell."""
        return self._neg_surface_cards

    @property
    def pos_surface_cards(self):
        """Surfaces with a positive sense for this cell."""
        return self._pos_surface_cards

    @property
    def neutron_importance(self):
        """MCNPX only."""
        return self._neutron_importance


class Cell(CellVoid):
    """
    """

    def __init__(self, name, number, material, density, density_units,
                 neg_surf_cards, pos_surf_cards, neutron_imp, temp, vol):
        """

        """
        return

    def comment(self):
        return

    @property
    def material(self):
        """A pyne.material object specifying the material filling the cell."""
        return self._materialcard

    @property
    def density(self):
        """Density for the material, in units of self.density_units."""
        return self._density

    @property
    def density_units(self):
        """Either 'g/cm^3', or 'atoms/b/cm'."""
        return self._density_units

    @property
    def temperature(self):
        """Temperature of the cell, in Kelvin."""
        return self._temperature

    @property
    def volume(self):
        """Volume of the cell, in cm^3. Optional."""
        return self._volume
        

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


class Data(Card):
    """All cards classified as data cards in MCNP, except for materials."""

    def __init__(self, name, number):
        """
        """
        return

    def comment(self):
        """This should be abstract as well."""
        return






