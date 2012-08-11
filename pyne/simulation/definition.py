#!/usr/bin/env python

"""Another possible name is latticesim or latticesimulation.

"""

import abc
import collections

from pyne import material

class SimulationDefinition(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def _create_new(self):
        """Definition started from scratch. Initialize all fields. """
        return

    @abc.abstractmethod
    def _open(self, fname):
        """Open object data from a JSON file."""
        return

    @abc.abstractmethod
    def save(self, fname):
        """Save object data to a JSON file."""
        return


class SystemDefinition(SimulationDefinition):
    """This class creates a reactor definition as is done in MCNPX: homogeneous
    regions in space in the reactor, called cells, are defined through the
    intersection, union, etc of surfaces and are filled by materials. The
    definition of materials is done using the `material` module of PyNE.
    
    """

    def __init__(self, fname=None):
        """Creates a new reactor definition or loads one from a JSON file."""

        if fname is not None:
            self._open(fname)
        else:
            self._create_new()

    def _create_new(self):
        self.surfaces = collections.OrderedDict()
        self.materials = collections.OrderedDict()
        self.cells = collections.OrderedDict()

    def save(self, fname):
        """Saves reactor definition to a JSON file."""

    def _open(self, fname):


class OptionsDefinition(SimulationDefinition):
    """This may need to be subclassed for different codes."""

    def __init__(self):
        return

    def add_criticality_source(self):
        return

    def add_criticality_points(self):
        return


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

    def __init__(self, name, number):
        """
        """

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

    def __init__(self, name, number):
        """

        """

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

    def __init__(self, name, number):
        return











