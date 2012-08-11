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
    def save(self, fname):
        """Save object data to a JSON file."""
        return

    @abc.abstractmethod
    def _open(self, fname):
        """Open object data from a JSON file."""
        return


class SystemDefinition(SimulationDefinition):
    """This class creates a reactor definition as is done in MCNPX: homogeneous
    regions in space in the reactor, called cells, are defined through the
    intersection, union, etc of surfaces and are filled by materials. The
    definition of materials is done using the `material` module of PyNE.
    
    """

    def __init__(self, fname=None):
        """Creates a new reactor definition or loads one from a JSON file."""

        if fname not None:
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
        pass

    def add_criticality_source(self):

    def add_criticality_points(self):


class Card(object):
    # each card should still have a comment string.

