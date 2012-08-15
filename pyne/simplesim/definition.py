#!/usr/bin/env python

"""The ``definition`` module can be imported as such::

    from pyne.simplesim import definition

Below is the reference for this module.


Another possible name is latticesim or latticesimulation.

"""

import abc
import collections

from pyne import material
from pyne.simulation import cards

class IDefinition(object):
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

    def _unique(self, card_type, name):
        """Checks that the name on a card has not already been used for another
        card.

        """
        if card_type == "cell":
            dict_to_check = self.cells
        elif card_type == "surface":
            dict_to_check = self.surfaces
        elif card_type == "material":
            dict_to_check = self.materials
        elif card_type == "tally":
            dict_to_check = self.tallies
        else:
            raise ValueError("The input ``card_type`` must be either "
                    "'cell', 'surface', or 'material'.")
        if name in dict_to_check:
            raise Exception("The %s name %s has already been used for "
                    "another %s" % (card_type, name, card_type))


class SystemDefinition(IDefinition):
    """This class creates a system definition as is done in MCNPX: homogeneous
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

    def add_cell(self, cell):
        """

        """
        self._unique("cell", name)
        # Check for the material.
        if cell.material.name not in self.materials:
            raise ValueError("Material %s is not found in this system
            definition." % cell.material.name)
        # Check for all required cells.
        # Okay, all checks passed.
        self.cell[cell.name] = cell

    def add_material(self, material):
        if material.name == None or material.name == '':
            raise ValueError("The ``name`` property of the material cannot "
                    "be empty.")
        self._unique("material", name)
        self.materials[material.name] = material

    def save(self, fname):
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        return

    def _open(self, fname):
        return


class SimulationDefinition(IDefinition):
    """This is basically where all the data cards are stored. The easy name for
    this class is either OptionsDefinition (Serpent) or DataDefinition (MCNP),
    but I'm not too happy with either. I'd like any ideas for this. This may
    need to be subclassed for different codes, because different codes do not
    provide the same options.
    
    """

    def __init__(self, systemdef, fname=None):
        """Creates a new options definition or loads one from a JSON file."""
        self.sys = systemdef
        if fname is not None:
            self._open(fname)
        else:
            self._create_new()

    def _open(self, fname):
        return

    def add_source(self, source):
        return

    def add_card(self, card):
        return

    def add_tally(self, card):
        return

    def _create_new(self):
        """Initialize any attributes/properties."""

    def save(self, fname):
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        return

    @property
    def sys(self):
        return self._sys

    @sys.setter
    def sys(self, value):
        self._sys = value

class MCNPSimulation(SimulationDefinition):
    pass
