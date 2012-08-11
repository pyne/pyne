#!/usr/bin/env python

"""The ``definition`` module can be imported as such::

    from pyne.simulation import definition

Below is the reference for this module.


Another possible name is latticesim or latticesimulation.

"""

import abc
import collections

from pyne import material
from pyne.simulation import cards

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
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        return

    def _open(self, fname):
        return

    def add_cell(self):
        return

    def add_cylinder(self):
        return


class OptionsDefinition(SimulationDefinition):
    """This is basically where all the data cards are stored. The easy name for
    this class is either OptionsDefinition (Serpent) or DataDefinition (MCNP),
    but I'm not too happy with either. I'd like any ideas for this. This may
    need to be subclassed for different codes, because different codes do not
    provide the same options."""

    def __init__(self, fname=None):
        """Creates a new options definition or loads one from a JSON file."""

        if fname is not None:
            self._open(fname)
        else:
            self._create_new()

    def _create_new(self):
        """Initialize any attributes/properties."""

    def save(self, fname):
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        return

    def _open(self, fname):
        return

    def add_criticality_source(self):
        return

    def add_criticality_points(self):
        return

    def add_tally_cellflux(self):
        return

