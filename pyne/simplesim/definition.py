#!/usr/bin/env python

"""The ``definition`` module can be imported as such::

    from pyne.simplesim import definition

Below is the reference for this module.


Another possible name is latticesim or latticesimulation.

"""

import abc
import collections
import pickle
import json

import numpy as np

from pyne import material
from pyne.simulation import cards

class IDefinition(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, fname=None, verbose=True):
        """

        """
        self.verbose = verbose
        if fname is not None:
            if self.verbose:
                print "Opening definition stored in %s." % fname
            self._open(fname)
        else:
            if self.verbose:
                print "Creating a new definition."
            self._create_new()

    @abc.abstractmethod
    def _create_new(self):
        """Definition started from scratch. Initialize all fields. """
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, fname):
        """Save object data to a JSON file."""
        fid = open(fname, 'w')
        #pickle.dump(self, fid)
        fid.write(DefinitionEncoder().encode(self))
        fid.close()

    @abc.abstractmethod
    def _open(self, fname):
        """Open object data from a JSON file."""
        # TODO Ensure the file exists.
        fid = open(fname)
        self = pickle.load(fname)
        fid.close()

    def _assert_unique(self, card_type, card):
        """Checks that the name on a card has not already been used for another
        card.

        """
        # TODO Don't require the user to pass the card_type string.
        if card_type == 'cell':
            dict_to_check = self.cells
        elif card_type == 'surface':
            dict_to_check = self.surfaces
        elif card_type == 'material':
            dict_to_check = self.materials
        elif card_type == 'source':
            dict_to_check = self.source
        elif card_type == 'tally':
            dict_to_check = self.tally
        elif card_type == 'misc':
            dict_to_check = self.misc
        else:
            raise ValueError("The input ``card_type`` must be either "
                    "'cell', 'surface', 'material', 'source', "
                    "'tally', or 'misc'.")
        if card.name in dict_to_check:
            raise Exception("The %s name %s has already been used for "
                    "another %s" % (card_type, card.name, card_type))

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = verbose

    @verbose.setter
    def verbose(self, value):
        # TODO boolean, check or rely on ducktyping?
        if type(value) is not bool:
            raise TypeError("The ``verbose`` property must be of type bool.")
        self._verbose = value


class SystemDefinition(IDefinition):
    """This class creates a system definition as is done in MCNPX: homogeneous
    regions in space in the reactor, called cells, are defined through the
    intersection, union, etc of surfaces and are filled by materials. The
    definition of materials is done using the `material` module of PyNE.

    """
    def __init__(self, fname=None, verbose=True):
        """Creates a new reactor definition or loads one from a JSON file."""
        super(SystemDefinition, self).__init__(fname, verbose)

    def _create_new(self):
        self.surfaces = collections.OrderedDict()
        self.materials = collections.OrderedDict()
        self.cells = collections.OrderedDict()

    def add_cell(self, cell):
        """

        """
        if self.verbose:
            print "Adding cell %s." % cell.name
        # Only add the cell if a cell with the same name hasn't been added
        # already.
        self._assert_unique('cell', cell)
        # Add all surfaces that aren't already added. Do this by walking the
        # region tree and calling _add_unique_surfaces() at the leaves.
        cell.region.walk(self._add_unique_surfaces)
        # Only add the material if this is not a void cell and if it doesn't
        # already exist.
        if (hasattr(cell, 'material') and
                cell.material.name not in self.materials):
            self.add_material(cell.material)
        # Constituent surfaces and material have been added, so we can added
        # the cell.
        self.cells[cell.name] = cell

    def add_surface(self, surface):
        """This method is only used by the user for surfaces that are not on a
        cell card. Surfaces on a cell card are added automatically. Duplication
        is checked by card name, not by using Python's 'is' operator to compare
        the objects.

        """
        self._assert_unique('surface', surface)
        self.surfaces[surface.name] = surface

    def _add_unique_surfaces(self, regionleaf):
        name = regionleaf.surface.name
        if self.verbose:
            print "Trying to add surface %s..." % name
        if name not in self.surfaces:
            self.add_surface(regionleaf.surface)
            if self.verbose:
                print "Surface %s added successfully." % name
        else:
            if self.verbose:
                print "Surface %s already exists in the definition." % name

    def add_material(self, material):
        """This method is only used by the user for materials that are not on a
        cell card.  Materials on cell cards are added automatically.

        """
        if material.name == None or material.name == '':
            raise ValueError("The ``name`` property of the material cannot "
                    "be empty.")
        self._assert_unique('material', material)
        self.materials[material.name] = material

    def save(self, fname):
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        super(SystemDefinition, self).save(fname)

    def _open(self, fname):
        super(SystemDefinition, self)._open(fname)


class SimulationDefinition(IDefinition):
    """This is basically where all the data cards are stored. The easy name for
    this class is either OptionsDefinition (Serpent) or DataDefinition (MCNP),
    but I'm not too happy with either. I'd like any ideas for this. This may
    need to be subclassed for different codes, because different codes do not
    provide the same options.

    """
    # TODO when adding a criticality points card, check that a criticality card
    # has been added as well.

    def __init__(self, systemdef, fname=None, verbose=True):
        """Creates a new options definition or loads one from a JSON file."""
        super(SimulationDefinition, self).__init__(fname, verbose)
        self.sys = systemdef

    def _create_new(self):
        """Initialize any attributes/properties."""
        self.source = collections.OrderedDict()
        self.misc = collections.OrderedDict()
        self.tally = collections.OrderedDict()

    def add_card(self, card):
        if issubclass(card, cards.ISource):
            self._add_source(card)
        elif issubclass(card, cards.ITally):
            self._add_tally(card)
        elif issubclass(card, cards.IMisc):
            self._add_misc(card)

    def _add_source(self, card):
        self._assert_unique('source', card)
        self.source[card.name] = card

    def _add_tally(self, card):
        self._assert_unique('tally', card)
        self.tally[card.name] = card

    def _add_misc(self, card):
        self._assert_unique('misc', card)
        self.misc[card.name] = card

    def save(self, fname):
        """Saves definition to a JSON file. It is unlikely that the class will
        be amenable to json.dump()."""
        return

    def _open(self, fname):
        pass

    @property
    def sys(self):
        return self._sys

    @sys.setter
    def sys(self, value):
        self._sys = value


class MCNPSimulation(SimulationDefinition):
    """ """




class DefinitionEncoder(json.JSONEncoder):
    # TODO circular reference issue.
    def default(self, obj):
        try:
            if hasattr(obj, '__dict__'):
                print "a"
                print type(obj)
                print obj.__dict__
                return obj.__dict__
            if isinstance(obj, np.ndarray):
                print "b"
                print type(obj)
                print obj.tolist()
                return obj.tolist()
            #if isinstance(obj, cards.Cell) or issubclass(obj, cards.Cell):
            #    print "cell"
            #    return ''
            else:
                print "c"
                print type(obj)
                return json.JSONEncoder.default(self, obj)
        except:
            print "exception: "
            print type(obj)
            print obj
            raise

